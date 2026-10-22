*! version 1.0.0  20jul2026
*! Global metafrontier Malmquist productivity index

capture program drop dea_mf_malmquist
program define dea_mf_malmquist, rclass
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
            local _pre = substr(`"`0'"',1,`_cpos')
            local _post = substr(`"`0'"',`_cpos'+1,.)
            local 0 `"`_pre' inputs(`_inputs') `_post'"'
        }
        else local 0 `"`0', inputs(`_inputs')"'
    }

    syntax [varlist(numeric)] [if] [in], GROUP(varname) DMU(varname) ///
        TIME(varname numeric) [Inputs(varlist numeric) Outputs(varlist numeric) ///
        RTS(string) ORT(string) ORIentation(string) ENGINE(string) ///
        GENerate(name) TOLerance(real 1e-7) SAVing(string asis) REPLACE DETail]

    if "`outputs'" == "" local outputs `varlist'
    if "`inputs'" == "" | "`outputs'" == "" {
        di as error "specify inputs and outputs, or use: dea_mf_malmquist inputs = outputs"
        exit 198
    }
    if "`ort'" != "" & "`orientation'" != "" {
        di as error "specify only one of ort() or orientation()"
        exit 198
    }
    if "`ort'" != "" local orientation "`ort'"
    if "`rts'" == "" local rts "crs"
    if "`orientation'" == "" local orientation "input"
    if "`engine'" == "" local engine "auto"
    if "`generate'" == "" local generate "_mfm"
    local rts = lower("`rts'")
    local orientation = lower("`orientation'")
    if inlist("`rts'","c","ccr") local rts "crs"
    if inlist("`rts'","v","bcc") local rts "vrs"
    if inlist("`orientation'","i","in") local orientation "input"
    if inlist("`orientation'","o","out") local orientation "output"
    if !inlist("`rts'","crs","vrs") {
        di as error "rts() must be crs or vrs"
        exit 198
    }
    if !inlist("`orientation'","input","output") {
        di as error "orientation() must be input or output"
        exit 198
    }
    if !inlist("`engine'","auto","simplex","ipm") {
        di as error "engine() must be auto, simplex, or ipm"
        exit 198
    }
    if `tolerance' <= 0 {
        di as error "tolerance() must be positive"
        exit 198
    }
    if "`replace'" != "" & `"`saving'"' == "" {
        di as error "replace requires saving()"
        exit 198
    }

    local v_ge "`generate'_group_eff"
    local v_me "`generate'_meta_eff"
    local v_tg "`generate'_tgr"
    local v_gp "`generate'_gmpi"
    local v_tc "`generate'_tgch"
    local v_mp "`generate'_mmpi"
    local v_id "`generate'_identity_gap"
    foreach v in `v_ge' `v_me' `v_tg' `v_gp' `v_tc' `v_mp' `v_id' {
        confirm name `v'
        capture confirm new variable `v'
        if _rc {
            di as error "variable `v' already exists"
            exit 110
        }
    }

    marksample touse, novarlist
    markout `touse' `inputs' `outputs' `time'
    markout `touse' `dmu' `group', strok
    quietly count if `touse'
    if r(N)==0 error 2000
    local N = r(N)
    foreach v of varlist `inputs' `outputs' {
        quietly count if `touse' & `v'<=0
        if r(N) {
            di as error "inputs and outputs must be strictly positive; check `v'"
            exit 459
        }
    }

    tempvar gid order dup nper gmin gmax ge_tmp me_tmp
    quietly generate long `order' = _n
    quietly egen long `gid' = group(`group') if `touse', label
    quietly summarize `gid' if `touse', meanonly
    local G = r(max)
    quietly levelsof `time' if `touse', local(periods)
    local T : word count `periods'
    if `T'<2 {
        di as error "at least two periods are required"
        exit 198
    }
    quietly duplicates tag `dmu' `time' if `touse', generate(`dup')
    quietly count if `touse' & `dup'>0
    if r(N) {
        di as error "dmu() and time() must uniquely identify observations"
        exit 459
    }
    quietly bysort `dmu': egen long `nper' = total(`touse')
    quietly count if `touse' & `nper'!=`T'
    if r(N) {
        di as error "a balanced panel over all selected periods is required"
        exit 459
    }
    quietly bysort `dmu': egen long `gmin' = min(cond(`touse',`gid',.))
    quietly bysort `dmu': egen long `gmax' = max(cond(`touse',`gid',.))
    quietly count if `touse' & `gmin'!=`gmax'
    if r(N) {
        di as error "group membership must remain constant within DMU"
        exit 459
    }

    quietly generate double `v_ge' = .
    forvalues g=1/`G' {
        quietly dea `outputs' if `touse' & `gid'==`g', inputs(`inputs') ///
            model(te) rts(`rts') orientation(`orientation') stage(1) ///
            engine(`engine') tolerance(`tolerance') generate(`ge_tmp')
        quietly replace `v_ge' = `ge_tmp' if `touse' & `gid'==`g'
        capture drop `ge_tmp'
    }
    quietly dea `outputs' if `touse', inputs(`inputs') model(te) rts(`rts') ///
        orientation(`orientation') stage(1) engine(`engine') ///
        tolerance(`tolerance') generate(`me_tmp')
    quietly generate double `v_me' = `me_tmp' if `touse'
    quietly generate double `v_tg' = `v_me'/`v_ge' if `touse'

    sort `dmu' `time'
    quietly by `dmu': generate double `v_gp' = `v_ge'/`v_ge'[_n-1] ///
        if `touse' & _n>1
    quietly by `dmu': generate double `v_tc' = `v_tg'/`v_tg'[_n-1] ///
        if `touse' & _n>1
    quietly by `dmu': generate double `v_mp' = `v_me'/`v_me'[_n-1] ///
        if `touse' & _n>1
    quietly generate double `v_id' = abs(`v_mp'-`v_gp'*`v_tc') if `touse' & `v_mp'<.
    sort `order'

    label variable `v_ge' "Group global-frontier efficiency"
    label variable `v_me' "Meta global-frontier efficiency"
    label variable `v_tg' "Technology gap ratio"
    label variable `v_gp' "Group global Malmquist productivity change"
    label variable `v_tc' "Technology-gap change"
    label variable `v_mp' "Metafrontier global Malmquist index"
    label variable `v_id' "Absolute MFM identity residual"

    quietly count if `touse' & (`v_tg'<=0 | `v_tg'>1+max(`tolerance',1e-6))
    if r(N) {
        di as error "one or more technology gap ratios are outside (0,1]"
        exit 459
    }
    quietly summarize `v_id' if `touse', meanonly
    local max_gap = r(max)
    if `max_gap'>max(`tolerance',1e-8) {
        di as error "metafrontier Malmquist decomposition identity failed"
        exit 459
    }
    quietly count if `touse' & `v_mp'<.
    local N_trans = r(N)
    quietly summarize `v_gp' if `touse', meanonly
    local mean_gp = r(mean)
    quietly summarize `v_tc' if `touse', meanonly
    local mean_tc = r(mean)
    quietly summarize `v_mp' if `touse', meanonly
    local mean_mp = r(mean)

    di as text "{hline 72}"
    di as text "Global metafrontier Malmquist index"
    di as text "DMUs: " as result `=`N'/`T'' as text "   Groups: " as result `G' ///
        as text "   Periods: " as result `T' as text "   Transitions: " as result `N_trans'
    di as text "Mean group productivity change: " as result %9.6f `mean_gp'
    di as text "Mean technology-gap change:     " as result %9.6f `mean_tc'
    di as text "Mean meta productivity change:  " as result %9.6f `mean_mp'
    di as text "Identity: MMPI = GMPI x TGCH; max residual = " as result %9.2e `max_gap'
    di as text "{hline 72}"
    if "`detail'" != "" {
        list `dmu' `time' `group' `v_gp' `v_tc' `v_mp' if `touse' & `v_mp'<., ///
            noobs separator(0)
    }

    if `"`saving'"' != "" {
        preserve
        quietly keep if `touse'
        keep `dmu' `time' `group' `inputs' `outputs' `v_ge' `v_me' `v_tg' ///
            `v_gp' `v_tc' `v_mp' `v_id'
        if "`replace'" != "" save `saving', replace
        else save `saving'
        restore
    }

    return clear
    return scalar N = `N'
    return scalar N_groups = `G'
    return scalar N_periods = `T'
    return scalar N_transitions = `N_trans'
    return scalar mean_group_mpi = `mean_gp'
    return scalar mean_tgch = `mean_tc'
    return scalar mean_meta_mpi = `mean_mp'
    return scalar max_identity_gap = `max_gap'
    return local inputs "`inputs'"
    return local outputs "`outputs'"
    return local group "`group'"
    return local dmu "`dmu'"
    return local time "`time'"
    return local rts "`rts'"
    return local orientation "`orientation'"
    return local index "global_metafrontier"
    return local generate "`generate'"
    return local saving `"`saving'"'
end
