*! version 1.0.1  02nov2026
*! Retry failed short regression fits on the SAME pseudo-data before redrawing.
*! Sequential Simar-Wilson Algorithm 2. Called after the original DEA and MLE.
program define _dea_sw2007_seq, eclass
    version 16.0
    syntax, INputs(varlist numeric) OUtputs(varlist numeric) ENVars(varlist numeric) ///
        TOUSE(varname) SCORES(varname) INITIAL(name) VINITIAL(name) SIGMA(real) ///
        REPS(integer) BCREPS(integer) LEVEL(real) TAU(real) ORIentation(string) ///
        RTS(string) TOLerance(real) GENerate(name) [TIME(varname) NOBCorrected TRACE]

    tempvar delta delta_bc xb draw boot total esample
    qui gen double `delta' = `tau'/`scores' if `touse'
    qui predict double `xb' if `touse', xb
    qui gen double `draw' = .
    qui gen double `boot' = .
    qui gen double `total' = 0 if `touse'
    qui count if `touse'
    local n = r(N)
    qui count if `touse' & `delta' > `tau'
    local n_initial = r(N)
    local p = colsof(`initial')
    local k = `p'-1
    local names : colnames `initial'
    local eqs : coleq `initial'
    local bc_ok = 0
    local bc_try = 0
    di as text "Step 3: `bcreps' successful pseudo-frontier replications..."
    while `bc_ok' < `bcreps' & `bc_try' < 2*`bcreps' {
        local ++bc_try
        mata: _sw07_seq_draw("`xb'", `sigma', `tau', "`touse'", "`draw'")
        qui count if `touse' & (missing(`draw') | `draw' <= `tau')
        if r(N) continue
        qui replace `boot' = .
        if "`time'" == "" {
            mata: _sw07_seq_lp("`inputs'", "`outputs'", "`scores'", "`draw'", ///
                `tau', "`touse'", "`boot'", "`orientation'", "`rts'", `tolerance')
        }
        else {
            qui levelsof `time' if `touse', local(periods)
            foreach t of local periods {
                tempvar period_sample
                qui gen byte `period_sample' = `touse' & `time' == `t'
                mata: _sw07_seq_lp("`inputs'", "`outputs'", "`scores'", "`draw'", ///
                    `tau', "`period_sample'", "`boot'", "`orientation'", "`rts'", `tolerance')
                drop `period_sample'
            }
        }
        // Accept or reject a WHOLE replication; every DMU has the same B1.
        qui count if `touse' & (missing(`boot') | `boot' <= 0)
        if r(N) continue
        qui replace `total' = `total' + `boot' if `touse'
        local ++bc_ok
        if "`trace'" != "" di as text "  DEA replication `bc_ok' (attempt `bc_try')"
    }
    if `bc_ok' != `bcreps' {
        di as error "Could not obtain `bcreps' valid DEA replications in " 2*`bcreps' " attempts"
        exit 498
    }
    // Correction must be made on the regression's inefficiency scale.
    qui gen double `delta_bc' = 2*`delta' - `total'/`bcreps' if `touse'
    qui count if `touse' & (missing(`delta_bc') | `delta_bc' <= `tau')
    if r(N) {
        di as error "Corrected inverse scores must exceed tau(); no scores were clipped"
        exit 498
    }
    di as text "Steps 4-5: correct inverse scores and refit truncated regression..."
    qui truncreg `delta_bc' `envars' if `touse', ll(`tau') iterate(200)
    if !e(converged) exit 430
    if e(N) != `n' | colsof(e(b)) != `p' exit 498
    tempname point modelV draws cov mean ci pivot bc blo bhi sigmas betas
    matrix `point' = e(b)
    matrix `modelV' = e(V)
    local sigma_refit = e(sigma)
    local iterations = e(ic)
    local ll = e(ll)
    qui gen byte `esample' = e(sample)
    drop `xb'
    qui predict double `xb' if `touse', xb
    matrix `draws' = J(`reps', `p', .)
    local reg_ok = 0
    local reg_try = 0
    local reg_fits = 0
    local reg_retried = 0
    local reg_recovered = 0
    local maxiter = max(25,3*`iterations')
    di as text "Step 6: `reps' successful parametric regression replications..."
    while `reg_ok' < `reps' & `reg_try' < 2*`reps' {
        local ++reg_try
        mata: _sw07_seq_draw("`xb'", `sigma_refit', `tau', "`touse'", "`draw'")
        qui count if `touse' & (missing(`draw') | `draw' <= `tau')
        if r(N) continue
        local ++reg_fits
        capture qui truncreg `draw' `envars' if `touse', ll(`tau') iterate(`maxiter')
        local fit_ok = (_rc == 0)
        if `fit_ok' {
            local fit_ok = (e(converged) == 1 & e(N) == `n' & ///
                colsof(e(b)) == `p' & !matmissing(e(b)))
        }
        // A short iteration cap is a speed safeguard, not evidence that this
        // pseudo-sample has no valid MLE. Refit it once with the same 200-step
        // limit used for the original/refitted model. Do not draw new errors.
        if !`fit_ok' & `maxiter' < 200 {
            local ++reg_retried
            local ++reg_fits
            capture qui truncreg `draw' `envars' if `touse', ll(`tau') iterate(200)
            local fit_ok = (_rc == 0)
            if `fit_ok' {
                local fit_ok = (e(converged) == 1 & e(N) == `n' & ///
                    colsof(e(b)) == `p' & !matmissing(e(b)))
            }
            if `fit_ok' local ++reg_recovered
        }
        if !`fit_ok' continue
        local ++reg_ok
        matrix `draws'[`reg_ok',1] = e(b)
        if "`trace'" != "" di as text "  Regression replication `reg_ok' (attempt `reg_try')"
    }
    if `reg_ok' != `reps' {
        di as error "Could not obtain `reps' converged regression replications in " 2*`reps' " attempts"
        exit 498
    }
    // Empirical inverse-CDF percentiles; pivot intervals are also returned.
    mata: _sw07_seq_summary("`draws'", "`point'", `level', "`cov'", "`mean'", "`ci'", "`pivot'")
    foreach m in draws cov mean ci pivot {
        matrix colnames ``m'' = `names'
        matrix coleq ``m'' = `eqs'
    }
    matrix rownames `cov' = `names'
    matrix roweq `cov' = `eqs'
    matrix rownames `ci' = lower upper
    matrix rownames `pivot' = lower upper
    matrix `bc' = (2*`point'[1,1..`k'] - `mean'[1,1..`k'])'
    matrix `blo' = `ci'[1,1..`k']'
    matrix `bhi' = `ci'[2,1..`k']'
    matrix `betas' = `draws'[1..`reps',1..`k']
    matrix `sigmas' = `draws'[1..`reps',`p']
    // Write outputs only after the complete procedure has succeeded.
    capture confirm variable `generate'_eff, exact
    if !_rc drop `generate'_eff
    qui gen double `generate'_eff = `scores' if `touse'
    label var `generate'_eff "DEA efficiency (original)"
    if "`nobcorrected'" == "" {
        capture confirm variable `generate'_eff_bc, exact
        if !_rc drop `generate'_eff_bc
        qui gen double `generate'_eff_bc = `tau'/`delta_bc' if `touse'
        label var `generate'_eff_bc "Efficiency from SW2007 corrected inverse score"
    }
    ereturn post `point' `cov', esample(`esample') obs(`n') depname(inefficiency)
    ereturn matrix b_initial = `initial'
    ereturn matrix V_initial = `vinitial'
    ereturn matrix V_modelbased = `modelV'
    ereturn matrix b_reps = `draws'
    ereturn matrix b_boot = `betas'
    ereturn matrix sigma_boot = `sigmas'
    ereturn matrix b_mean = `mean'
    ereturn matrix ci_percentile = `ci'
    ereturn matrix ci_pivot = `pivot'
    ereturn matrix b_bc = `bc'
    matrix `bc' = e(b_bc)
    ereturn matrix beta_bc = `bc'
    ereturn matrix ci_lo = `blo'
    ereturn matrix ci_hi = `bhi'
    ereturn scalar sigma0 = `sigma'
    ereturn scalar sigma = `sigma_refit'
    ereturn scalar N_initial = `n_initial'
    ereturn scalar N_trunc = `n'
    ereturn scalar reps = `reps'
    ereturn scalar bcreps = `bcreps'
    ereturn scalar boot_success = `reg_ok'
    ereturn scalar bc_success = `bc_ok'
    ereturn scalar boot_failed = `reg_try'-`reg_ok'
    ereturn scalar boot_fit_attempts = `reg_fits'
    ereturn scalar boot_retried = `reg_retried'
    ereturn scalar boot_recovered = `reg_recovered'
    ereturn scalar bc_failed = `bc_try'-`bc_ok'
    ereturn scalar dea_reps = `bc_try'
    ereturn scalar level = `level'
    ereturn scalar tau = `tau'
    ereturn scalar ll = `ll'
    ereturn scalar converged = 1
    ereturn local cmd "dea_sw2007"
    ereturn local backend "internal"
    ereturn local procedure "sequential"
    ereturn local envars "`envars'"
    ereturn local orientation "`orientation'"
    ereturn local rts "`rts'"
    ereturn local vce "bootstrap"
    ereturn local vcetype "Bootstrap"
    di as text "Algorithm 2 completed: e(b) is the corrected-score fit; e(V) is bootstrap covariance."
    matrix list e(b), noheader format(%10.5f)
    di as text "Percentile intervals: e(ci_percentile); pivot intervals: e(ci_pivot)."
end

mata:
// Inverse survival CDF avoids cancellation in 1-Phi(a). For extreme lower
// truncation use an exponential proposal with exact rejection (no clipping).
void _sw07_seq_draw(string scalar xbvar, real scalar sigma, real scalar tau,
    string scalar touse, string scalar dest)
{
    real colvector ix, mu, out, u
    real scalar i, a, rate, z, j, p
    ix=selectindex(st_data(.,touse)); mu=st_data(ix,xbvar)
    out=J(rows(ix),1,.); u=runiform(rows(ix),1)
    for (i=1;i<=rows(ix);i++) {
        a=(tau-mu[i])/sigma
        if (a<=8) {
            p=(1-u[i])*normal(-a)
            if (p>0 & p<1) out[i]=mu[i]-sigma*invnormal(p)
        }
        else if (a<.) {
            rate=(a+sqrt(a*a+4))/2
            for(j=1;j<=10000;j++) {
                z=a-ln(1-runiform(1,1))/rate
                if (runiform(1,1)<=exp(-.5*(z-rate)^2)) {
                    out[i]=tau+sigma*(z-a)
                    break
                }
            }
        }
    }
    st_store(ix,dest,out)
}

// Return delta**=tau/theta** on the SAME scale used by truncated regression.
void _sw07_seq_lp(string scalar ins, string scalar outs, string scalar scores,
    string scalar draw, real scalar tau, string scalar touse, string scalar dest,
    string scalar orientation, string scalar rts, real scalar tol)
{
    real colvector ix, theta, d, ratio, result
    real matrix X,Y,XX,YY
    real scalar i,th
    struct RSM_Prob scalar prob
    struct DEA_Info scalar info
    struct RSM_Result scalar res
    ix=selectindex(st_data(.,touse))
    X=st_data(ix,tokens(ins)); Y=st_data(ix,tokens(outs))
    X=X:/mean(X); Y=Y:/mean(Y)
    theta=st_data(ix,scores); d=st_data(ix,draw)
    ratio=tau:/(theta:*d)
    XX=X; YY=Y
    if (orientation=="output") YY=Y:*ratio
    else XX=X:/ratio
    result=J(rows(ix),1,.)
    for(i=1;i<=rows(ix);i++) {
        _dea_lp_stage1(prob,info,XX,YY,X[i,]',Y[i,]',orientation,rts)
        res=_rsm_two_phase(prob,tol,"")
        if(res.status!=0) continue
        th=_dea_extract_theta(res,info)
        if(th<=0 | th>=.) continue
        result[i]=(orientation=="output" ? tau*th : tau/th)
    }
    st_store(ix,dest,result)
}

void _sw07_seq_summary(string scalar bm,string scalar pm,real scalar level,
    string scalar vm,string scalar mm,string scalar cm,string scalar pivotm)
{
    real matrix B,C,P
    real rowvector point
    real colvector sorted
    real scalar j,n,lo,hi
    B=st_matrix(bm); point=st_matrix(pm); n=rows(B)
    lo=max((1,ceil(n*(100-level)/200)))
    hi=min((n,ceil(n*(1-(100-level)/200))))
    C=J(2,cols(B),.)
    for(j=1;j<=cols(B);j++) {
        sorted=sort(B[.,j],1)
        C[.,j]=(sorted[lo]\sorted[hi])
    }
    P=(2*point-C[2,.])\(2*point-C[1,.])
    st_matrix(vm,variance(B)); st_matrix(mm,mean(B))
    st_matrix(cm,C); st_matrix(pivotm,P)
}
end
