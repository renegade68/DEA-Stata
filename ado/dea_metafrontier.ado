*! version 1.4.0  20jul2026
*! Nonparametric DEA metafrontier: group and pooled frontiers

capture program drop dea_metafrontier
program define dea_metafrontier, rclass
    version 16.0

    local raw0 `"`0'"'
    if strpos(`"`raw0'"', "=") > 0 {
        local _inputs ""
        gettoken _w 0 : 0, parse("=,")
        while "`_w'" != "=" & "`_w'" != "" {
            if "`_w'" == "," {
                di as error "syntax error: comma before '='"
                exit 198
            }
            local _inputs `_inputs' `_w'
            gettoken _w 0 : 0, parse("=,")
        }
        local _cpos = strpos(`"`0'"', ",")
        if `_cpos' > 0 {
            local _pre  = substr(`"`0'"', 1, `_cpos')
            local _post = substr(`"`0'"', `_cpos' + 1, .)
            local 0 `"`_pre' inputs(`_inputs') `_post'"'
        }
        else local 0 `"`0', inputs(`_inputs')"'
    }

    syntax [varlist(numeric)] [if] [in], GROUP(varname) ///
        [Inputs(varlist numeric) Outputs(varlist numeric) DMU(varname) ///
         METHOD(string) RTS(string) ORT(string) ORIentation(string) ///
         GENerate(name) ENGINE(string) TOLerance(real 1e-7) ///
         MINOBS(integer 0) POOLTEST(integer 0) SEED(integer 0) ///
         BREPS(integer 0) Level(real 95) ///
         SAVING(string asis) REPLACE DETail]

    if "`outputs'" == "" local outputs `varlist'
    if "`inputs'" == "" | "`outputs'" == "" {
        di as error "specify inputs and outputs, or use: dea_metafrontier inputs = outputs, group()"
        exit 198
    }
    if "`ort'" != "" & "`orientation'" != "" {
        di as error "specify only one of ort() or orientation()"
        exit 198
    }
    if "`ort'" != "" local orientation "`ort'"
    if "`method'" == ""      local method "dea"
    if "`rts'" == ""         local rts "vrs"
    if "`orientation'" == "" local orientation "input"
    if "`generate'" == ""    local generate "_mf"
    if "`engine'" == ""      local engine "auto"

    local method      = lower("`method'")
    local rts         = lower("`rts'")
    local orientation = lower("`orientation'")
    if inlist("`rts'", "c", "ccr") local rts "crs"
    if inlist("`rts'", "v", "bcc") local rts "vrs"
    if inlist("`rts'", "d", "nirs") local rts "drs"
    if inlist("`rts'", "i", "ndrs") local rts "irs"
    if inlist("`orientation'", "i", "in")  local orientation "input"
    if inlist("`orientation'", "o", "out") local orientation "output"

    if "`method'" != "dea" {
        di as error "method() currently supports only dea"
        exit 198
    }
    if !inlist("`rts'", "crs", "vrs", "drs", "irs") {
        di as error "rts() must be crs, vrs, drs, or irs"
        exit 198
    }
    if !inlist("`orientation'", "input", "output") {
        di as error "orientation() must be input or output"
        exit 198
    }
    if !inlist("`engine'", "auto", "simplex", "ipm") {
        di as error "engine() must be auto, simplex, or ipm"
        exit 198
    }
    if `tolerance' <= 0 {
        di as error "tolerance() must be positive"
        exit 198
    }
    if `minobs' < 0 {
        di as error "minobs() must be nonnegative"
        exit 198
    }
    if `pooltest' < 0 {
        di as error "pooltest() must be nonnegative"
        exit 198
    }
    if `breps' != 0 & `breps' < 50 {
        di as error "breps() must be zero or at least 50"
        exit 198
    }
    if `level' <= 0 | `level' >= 100 {
        di as error "level() must be between 0 and 100"
        exit 198
    }
    if `seed' < 0 {
        di as error "seed() must be nonnegative"
        exit 198
    }
    if "`replace'" != "" & `"`saving'"' == "" {
        di as error "replace requires saving()"
        exit 198
    }

    local v_group "`generate'_group"
    local v_meta  "`generate'_meta"
    local v_tgr   "`generate'_tgr"
    local v_gap   "`generate'_gap"
    local bootstrap_vars ""
    if `breps' > 0 {
        local v_tgr_bc "`generate'_tgr_bc"
        local v_tgr_se "`generate'_tgr_se"
        local v_tgr_lo "`generate'_tgr_lo"
        local v_tgr_hi "`generate'_tgr_hi"
        local bootstrap_vars `v_tgr_bc' `v_tgr_se' `v_tgr_lo' `v_tgr_hi'
    }
    foreach v in `v_group' `v_meta' `v_tgr' `v_gap' `bootstrap_vars' {
        confirm name `v'
        capture confirm new variable `v'
        if _rc {
            di as error "variable `v' already exists"
            exit 110
        }
    }

    marksample touse, novarlist
    markout `touse' `inputs' `outputs'
    markout `touse' `group', strok
    if "`dmu'" != "" markout `touse' `dmu', strok
    quietly count if `touse'
    if r(N) == 0 error 2000
    local N = r(N)

    foreach v of varlist `inputs' `outputs' {
        quietly count if `touse' & `v' <= 0
        if r(N) {
            di as error "inputs and outputs must be strictly positive; check `v'"
            exit 459
        }
    }

    tempvar gid gscore mscore
    quietly egen long `gid' = group(`group') if `touse', label
    quietly summarize `gid' if `touse', meanonly
    local G = r(max)
    local m : word count `inputs'
    local s : word count `outputs'
    local recommended_N = 3 * (`m' + `s')
    local min_group_N = .

    forvalues g = 1/`G' {
        quietly count if `touse' & `gid' == `g'
        if missing(`min_group_N') | r(N) < `min_group_N' local min_group_N = r(N)
        if `minobs' > 0 & r(N) < `minobs' {
            di as error "group `g' has " r(N) " observations; minobs(`minobs') required"
            exit 459
        }
    }

    quietly generate double `v_group' = .
    forvalues g = 1/`G' {
        quietly dea `outputs' if `touse' & inlist(`gid', `g'), inputs(`inputs') ///
            model(te) rts(`rts') ort(`orientation') stage(1) ///
            generate(`gscore') engine(`engine') tolerance(`tolerance')
        quietly replace `v_group' = `gscore' if `touse' & `gid' == `g'
        capture drop `gscore'
    }

    quietly dea `outputs' if `touse', inputs(`inputs') model(te) rts(`rts') ///
        ort(`orientation') stage(1) generate(`mscore') ///
        engine(`engine') tolerance(`tolerance')
    quietly generate double `v_meta' = `mscore' if `touse'
    quietly generate double `v_tgr'  = `v_meta' / `v_group' if `touse'
    quietly generate double `v_gap'  = 1 - `v_tgr' if `touse'

    label variable `v_group' "Group-frontier DEA efficiency"
    label variable `v_meta'  "Pooled metafrontier DEA efficiency"
    label variable `v_tgr'   "Technology gap ratio"
    label variable `v_gap'   "Technology gap (1-TGR)"

    quietly count if `touse' & (`v_group' <= 0 | missing(`v_group'))
    if r(N) {
        di as error "one or more group-frontier problems were not solved"
        exit 430
    }
    quietly count if `touse' & (`v_meta' <= 0 | missing(`v_meta'))
    if r(N) {
        di as error "one or more pooled-frontier problems were not solved"
        exit 430
    }
    quietly count if `touse' & `v_meta' > `v_group' + max(`tolerance', 1e-6)
    if r(N) {
        di as error "pooled efficiency exceeds group efficiency beyond tolerance"
        exit 459
    }

    quietly summarize `v_group' if `touse', meanonly
    local mean_group = r(mean)
    quietly summarize `v_meta' if `touse', meanonly
    local mean_meta = r(mean)
    quietly summarize `v_tgr' if `touse', meanonly
    local mean_tgr = r(mean)
    local pool_stat = max(0, 1 - `mean_tgr')

    local mean_tgr_bc = .
    local mean_tgr_se = .
    if `breps' > 0 {
        quietly generate double `v_tgr_bc' = . if `touse'
        quietly generate double `v_tgr_se' = . if `touse'
        quietly generate double `v_tgr_lo' = . if `touse'
        quietly generate double `v_tgr_hi' = . if `touse'

        local old_bs_rngstate `"`c(rngstate)'"'
        if `seed' > 0 quietly set seed `seed'
        capture findfile dea_bootstrap.ado
        if _rc {
            di as error "dea_bootstrap.ado is required for breps()"
            exit 601
        }
        quietly run `"`r(fn)'"'
        mata: _dea_mf_bootstrap("`inputs'", "`outputs'", "`gid'", ///
            "`touse'", "`v_group'", "`v_meta'", "`v_tgr'", ///
            "`v_tgr_bc'", "`v_tgr_se'", "`v_tgr_lo'", "`v_tgr_hi'", ///
            "`orientation'", "`rts'", `breps', `level', `tolerance')
        local mean_tgr_bc = scalar(__mf_mean_tgr_bc)
        local mean_tgr_se = scalar(__mf_mean_tgr_se)
        scalar drop __mf_mean_tgr_bc __mf_mean_tgr_se
        quietly set rngstate `old_bs_rngstate'
        label variable `v_tgr_bc' "Bias-corrected technology gap ratio"
        label variable `v_tgr_se' "Bootstrap standard error of TGR"
        label variable `v_tgr_lo' "Bootstrap TGR CI lower bound"
        label variable `v_tgr_hi' "Bootstrap TGR CI upper bound"
    }

    matrix _mf_summary = J(`G', 4, .)
    forvalues g = 1/`G' {
        quietly count if `touse' & `gid' == `g'
        matrix _mf_summary[`g',4] = r(N)
        local j 0
        foreach v in `v_group' `v_meta' `v_tgr' {
            local ++j
            quietly summarize `v' if `touse' & `gid' == `g', meanonly
            matrix _mf_summary[`g',`j'] = r(mean)
        }
    }
    matrix colnames _mf_summary = mean_group mean_meta mean_tgr N

    local p_pool = .
    if `pooltest' > 0 {
        local old_rngstate `"`c(rngstate)'"'
        if `seed' > 0 quietly set seed `seed'

        tempvar original_order random_order perm_rank perm_gid perm_group perm_gap
        quietly generate long `original_order' = _n
        local exceed = 0
        forvalues b = 1/`pooltest' {
            quietly generate double `random_order' = runiform() if `touse'
            sort `touse' `random_order'
            quietly generate long `perm_rank' = sum(`touse')
            quietly generate long `perm_gid' = .
            local lo = 1
            forvalues g = 1/`G' {
                local ng = _mf_summary[`g',4]
                local hi = `lo' + `ng' - 1
                quietly replace `perm_gid' = `g' if `touse' & ///
                    inrange(`perm_rank', `lo', `hi')
                local lo = `hi' + 1
            }
            sort `original_order'

            quietly generate double `perm_group' = .
            forvalues g = 1/`G' {
                quietly dea `outputs' if `touse' & inlist(`perm_gid', `g'), ///
                    inputs(`inputs') model(te) rts(`rts') ///
                    ort(`orientation') stage(1) ///
                    generate(`gscore') engine(`engine') tolerance(`tolerance')
                quietly replace `perm_group' = `gscore' if `touse' & `perm_gid' == `g'
                capture drop `gscore'
            }
            quietly generate double `perm_gap' = max(0, 1 - `v_meta'/`perm_group') if `touse'
            quietly summarize `perm_gap' if `touse', meanonly
            if r(mean) >= `pool_stat' - `tolerance' local ++exceed
            drop `random_order' `perm_rank' `perm_gid' `perm_group' `perm_gap'
        }
        sort `original_order'
        local p_pool = (1 + `exceed') / (1 + `pooltest')
        quietly set rngstate `old_rngstate'
    }

    di as text "{hline 72}"
    di as text "DEA-based metafrontier"
    di as text "Groups: " as result `G' as text "   DMUs: " as result `N' ///
        as text "   RTS: " as result upper("`rts'") ///
        as text "   Orientation: " as result "`orientation'"
    di as text "{hline 72}"
    di as text "Mean group efficiency: " as result %9.6f `mean_group'
    di as text "Mean meta efficiency:  " as result %9.6f `mean_meta'
    di as text "Mean TGR:              " as result %9.6f `mean_tgr'
    if `breps' > 0 {
        di as text "Mean bias-corrected TGR: " as result %9.6f `mean_tgr_bc'
        di as text "Mean TGR bootstrap SE:   " as result %9.6f `mean_tgr_se' ///
            as text " (" `breps' " replications, " %4.1f `level' "% CI)"
    }
    if `pooltest' > 0 {
        di as text "Poolability permutation test: statistic = " ///
            as result %9.6f `pool_stat' as text ", p = " ///
            as result %7.4f `p_pool' as text " (" `pooltest' " permutations)"
    }
    if `min_group_N' < `recommended_N' {
        di as text "Sample-size note: smallest group has " as result `min_group_N' ///
            as text "; the descriptive 3*(inputs+outputs) rule is " ///
            as result `recommended_N' as text "."
    }
    if "`detail'" != "" {
        tabstat `v_group' `v_meta' `v_tgr' if `touse', by(`group') statistics(mean)
        if "`dmu'" != "" {
            list `dmu' `group' `v_group' `v_meta' `v_tgr' if `touse', ///
                noobs separator(0)
        }
    }

    if `"`saving'"' != "" {
        preserve
        quietly keep if `touse'
        local savevars `group' `inputs' `outputs' ///
            `v_group' `v_meta' `v_tgr' `v_gap' `bootstrap_vars'
        if "`dmu'" != "" local savevars `dmu' `savevars'
        keep `savevars'
        if "`replace'" != "" save `saving', replace
        else save `saving'
        restore
    }

    return clear
    return scalar N = `N'
    return scalar N_groups = `G'
    return scalar mean_group = `mean_group'
    return scalar mean_meta = `mean_meta'
    return scalar mean_tgr = `mean_tgr'
    return scalar min_group_N = `min_group_N'
    return scalar recommended_N = `recommended_N'
    return scalar pool_stat = `pool_stat'
    return scalar p_pool = `p_pool'
    return scalar pool_reps = `pooltest'
    return scalar bootstrap_reps = `breps'
    return scalar level = `level'
    return scalar mean_tgr_bc = `mean_tgr_bc'
    return scalar mean_tgr_se = `mean_tgr_se'
    return matrix group_summary = _mf_summary
    return local inputs "`inputs'"
    return local outputs "`outputs'"
    return local group "`group'"
    return local method "dea"
    return local rts "`rts'"
    return local orientation "`orientation'"
    return local generate "`generate'"
    return local saving `"`saving'"'
end

mata:

void _dea_mf_bootstrap(string scalar inputs,
                       string scalar outputs,
                       string scalar gid_var,
                       string scalar touse,
                       string scalar group_var,
                       string scalar meta_var,
                       string scalar tgr_var,
                       string scalar bc_var,
                       string scalar se_var,
                       string scalar lo_var,
                       string scalar hi_var,
                       string scalar orientation,
                       string scalar rts,
                       real scalar B,
                       real scalar level,
                       real scalar tol)
{
    real matrix X, Y, Xb, Yb, tg_boot
    real colvector gid, obs, eg, em, tg, beta, ratio, eg_b, em_b, idx
    real colvector bc, se, lo, hi
    real rowvector xs, ys, z, dev
    real scalar n, m, s, G, b, g, i, k, h, alpha, ilo, ihi, mb, a, c

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    gid = st_data(., gid_var, touse)
    eg = st_data(., group_var, touse)
    em = st_data(., meta_var, touse)
    tg = st_data(., tgr_var, touse)
    obs = selectindex(st_data(., touse))
    n = rows(X); m = cols(X); s = cols(Y); G = max(gid)

    xs = mean(X); ys = mean(Y)
    for (k=1; k<=m; k++) if (xs[k] < 1e-12) xs[k] = 1
    for (k=1; k<=s; k++) if (ys[k] < 1e-12) ys[k] = 1
    X = X :/ xs
    Y = Y :/ ys
    tg_boot = J(n, B, .)

    for (b=1; b<=B; b++) {
        Xb = X; Yb = Y
        for (g=1; g<=G; g++) {
            idx = selectindex(gid :== g)
            h = _sw_bandwidth(eg[idx])
            beta = _sw_sample(eg[idx], h)
            if (orientation == "input") {
                ratio = eg[idx] :/ beta
                Xb[idx,.] = ratio :* X[idx,.]
            }
            else {
                ratio = beta :/ eg[idx]
                Yb[idx,.] = ratio :* Y[idx,.]
            }
        }
        eg_b = J(n,1,.); em_b = J(n,1,.)
        for (i=1; i<=n; i++) {
            idx = selectindex(gid :== gid[i])
            eg_b[i] = _dea_bs_solve(X[i,.], Y[i,.], Xb[idx,.], Yb[idx,.],
                                     orientation, rts, tol)
            em_b[i] = _dea_bs_solve(X[i,.], Y[i,.], Xb, Yb,
                                     orientation, rts, tol)
            if (eg_b[i] != . & em_b[i] != . & eg_b[i] > tol)
                tg_boot[i,b] = em_b[i] / eg_b[i]
        }
    }

    bc=J(n,1,.); se=J(n,1,.); lo=J(n,1,.); hi=J(n,1,.)
    alpha=(100-level)/100
    for (i=1; i<=n; i++) {
        z = select(tg_boot[i,.], tg_boot[i,.] :!= .)
        if (cols(z)==0) continue
        mb = mean(z')
        bc[i] = 2*tg[i] - mb
        se[i] = sqrt(variance(z'))
        dev = sort((z :- tg[i])',1)'
        ilo = max((1,ceil(cols(dev)*alpha/2)))
        ihi = min((cols(dev),floor(cols(dev)*(1-alpha/2))))
        a = tg[i]-dev[ihi]
        c = tg[i]-dev[ilo]
        lo[i] = min((1,max((0,min((a,c))))))
        hi[i] = max((0,min((1,max((a,c))))))
        bc[i] = min((1,max((0,bc[i]))))
    }
    for (i=1; i<=n; i++) {
        st_store(obs[i],st_varindex(bc_var),bc[i])
        st_store(obs[i],st_varindex(se_var),se[i])
        st_store(obs[i],st_varindex(lo_var),lo[i])
        st_store(obs[i],st_varindex(hi_var),hi[i])
    }
    st_numscalar("__mf_mean_tgr_bc",mean(select(bc,bc:!=.)))
    st_numscalar("__mf_mean_tgr_se",mean(select(se,se:!=.)))
}

end
