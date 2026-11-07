*! version 1.0.4  05nov2026
*! Public worker for explicit convex joint weak-disposal DDF.
program define _dea_weakddf, rclass
    version 16.0
    syntax [if] [in], Inputs(varlist numeric) Outputs(varlist numeric) BAD(varlist numeric) ///
        GENerate(name) TOLerance(real) ENGine(string) DIRection(string) RTS(string) ///
        [DMU(varname) WEAKMATrix SAVing(string) REPLACE REPLACEVARS DEBUG]
    local names `generate' `generate'_status
    if "`dmu'"=="" {
        confirm new variable _dmu
        local reserved _dmu
        if `: list reserved in names' {
            di as error "_dmu is reserved for the automatic identifier"
            exit 110
        }
    }
    _dea_output_check, names(`names') protect(`inputs' `outputs' `bad' `dmu') `replacevars'
    local existing `r(existing)'
    marksample touse, novarlist
    markout `touse' `inputs' `outputs' `bad'
    if "`dmu'"!="" markout `touse' `dmu', strok
    quietly count if `touse'
    local n=r(N)
    if `n'==0 error 2000
    foreach v of varlist `inputs' `outputs' `bad' {
        quietly count if `touse' & `v'<0
        if r(N) {
            di as error "disposability(weak) requires nonnegative inputs, good and bad outputs"
            exit 459
        }
        quietly summarize `v' if `touse', meanonly
        local scale=r(max)
        quietly count if `touse' & `v'>0 & `v'/`scale'==0
        if r(N) {
            di as error "column scaling underflows for `v'; respecify the numerical units/range"
            exit 459
        }
    }
    foreach side in inputs outputs {
        tempvar positive
        quietly generate byte `positive'=0 if `touse'
        foreach v of local `side' {
            quietly replace `positive'=1 if `touse' & `v'>0
        }
        quietly count if `touse' & !`positive'
        if r(N) {
            di as error "disposability(weak) requires a positive input and good output in each selected DMU"
            exit 459
        }
    }
    local m : word count `inputs'
    local s : word count `outputs'
    local q : word count `bad'
    if "`weakmatrix'"!="" & max(`n',`m',`s',`q',10)>c(matsize) {
        di as error "weakmatrix requires matsize >= max(N, input/good/bad count, 10)"
        exit 908
    }
    if `"`saving'"'!="" {
        if !regexm(lower(`"`saving'"'), "\.dta$") local saving `"`saving'.dta"'
        if "`replace'"=="" confirm new file `"`saving'"'
    }
    capture mata: assert(_dea_weak_ready()==6 & _dea_simplex_audit_ready()==16 & findexternal("_rsm_duplicate_ratio()")!=NULL & findexternal("_rsm_duplicate_map()")!=NULL & findexternal("_rsm_duplicate_problem()")!=NULL & findexternal("_rsm_duplicate_accept()")!=NULL & findexternal("_rsm_duplicate_retry()")!=NULL & findexternal("_rsm_two_phase_core()")!=NULL & findexternal("_rsm_solve_lp()")!=NULL & findexternal("_rsm_quality_fast()")!=NULL & findexternal("_rsm_affine_quality()")!=NULL & findexternal("_rsm_affine_refine()")!=NULL & findexternal("_rsm_success_fastread()")!=NULL & findexternal("_rsm_success_reset()")!=NULL & findexternal("_rsm_success_read()")!=NULL & findexternal("_rsm_success_counts_read()")!=NULL & findexternal("_rsm_success_check()")!=NULL & findexternal("_rsm_success_apply()")!=NULL & findexternal("_rsm_xgrow()")!=NULL & findexternal("_rsm_xadd()")!=NULL & findexternal("_rsm_xmul()")!=NULL & findexternal("_rsm_xnormalize()")!=NULL & findexternal("_rsm_rowmap()")!=NULL & findexternal("_rsm_row_empty()")!=NULL & findexternal("_rsm_row_accept()")!=NULL & findexternal("_rsm_row_retry()")!=NULL & findexternal("_dea_link_reset()")!=NULL & findexternal("_dea_link_read()")!=NULL & findexternal("_dea_link_store()")!=NULL & findexternal("_dea_link_small()")!=NULL & findexternal("_dea_link_residual()")!=NULL & findexternal("_dea_link_basis()")!=NULL & findexternal("_dea_link_apply()")!=NULL & findexternal("_dea_link_stage2()")!=NULL)
    if _rc {
        quietly findfile _dea_weak.ado
        quietly run "`r(fn)'"
    }
    capture mata: assert(_dea_simplex_audit_ready()==16 & findexternal("_rsm_duplicate_ratio()")!=NULL & findexternal("_rsm_duplicate_map()")!=NULL & findexternal("_rsm_duplicate_problem()")!=NULL & findexternal("_rsm_duplicate_accept()")!=NULL & findexternal("_rsm_duplicate_retry()")!=NULL & findexternal("_rsm_two_phase_core()")!=NULL & findexternal("_rsm_solve_lp()")!=NULL & findexternal("_rsm_quality_fast()")!=NULL & findexternal("_rsm_affine_quality()")!=NULL & findexternal("_rsm_affine_refine()")!=NULL & findexternal("_rsm_success_fastread()")!=NULL & findexternal("_rsm_success_reset()")!=NULL & findexternal("_rsm_success_read()")!=NULL & findexternal("_rsm_success_counts_read()")!=NULL & findexternal("_rsm_success_check()")!=NULL & findexternal("_rsm_success_apply()")!=NULL & findexternal("_rsm_xgrow()")!=NULL & findexternal("_rsm_xadd()")!=NULL & findexternal("_rsm_xmul()")!=NULL & findexternal("_rsm_xnormalize()")!=NULL & findexternal("_rsm_rowmap()")!=NULL & findexternal("_rsm_row_empty()")!=NULL & findexternal("_rsm_row_accept()")!=NULL & findexternal("_rsm_row_retry()")!=NULL & findexternal("_dea_link_reset()")!=NULL & findexternal("_dea_link_read()")!=NULL & findexternal("_dea_link_store()")!=NULL & findexternal("_dea_link_small()")!=NULL & findexternal("_dea_link_residual()")!=NULL & findexternal("_dea_link_basis()")!=NULL & findexternal("_dea_link_apply()")!=NULL & findexternal("_dea_link_stage2()")!=NULL)
    if _rc {
        di as error "Incompatible simplex library; reinstall the complete matched DEA package."
        exit 498
    }

    tempvar score status
    quietly generate double `score'=.
    quietly generate double `status'=.
    local destinations ""
    forvalues j=1/12 {
        tempname mat`j'
        local destinations `destinations' `mat`j''
    }
    mata: _dea_weak_compute("`inputs'","`outputs'","`bad'","`touse'", ///
        "`score' `status'","`direction'","`rts'","`engine'",`tolerance',"`weakmatrix'"!="", ///
        "`debug'"!="","`destinations'")
    label variable `score' "Joint weak-disposal directional distance beta (0 = efficient)"
    label variable `status' "0 accepted weak-disposal LP; 1 numerical result unavailable"
    // Save before changing any existing output column in the user's dataset.
    if `"`saving'"'!="" {
        preserve
        foreach v of local existing {
            quietly drop `v'
        }
        rename `score' `generate'
        rename `status' `generate'_status
        local fileid `dmu'
        if "`dmu'"=="" {
            quietly generate long _dmu=_n if `touse'
            local fileid _dmu
        }
        quietly keep if `touse'
        keep `fileid' `inputs' `outputs' `bad' `names'
        capture noisily save `"`saving'"', `replace'
        local saverc=_rc
        restore
        if `saverc' exit `saverc'
    }
    foreach v of local existing {
        quietly drop `v'
    }
    rename `score' `generate'
    rename `status' `generate'_status
    if "`dmu'"=="" quietly generate long _dmu=_n if `touse'
    local failed=`mat1'[1,2]
    di as text "Joint weak-disposal DDF (`rts'; direction `direction'; beta = 0 is efficient)"
    di as text "  Observations: `n'; available scores: " as result `n'-`failed'
    local lp_backend "package_simplex"
    if "`engine'"=="ipm" local lp_backend "stata_linearprogram"
    if "`engine'"=="auto" local lp_backend "simplex_then_native"
    di as text "  Backend: `lp_backend'; accepted simplex/native: " as result `mat1'[1,6] "/" `mat1'[1,7]
    if `failed' di as text "  `failed' LP results unavailable; see `generate'_status (not an infeasibility certificate)."
    return scalar N=`n'
    return scalar N_complete=`n'-`failed'
    return scalar n_failed=`failed'
    return scalar N_eff=`mat1'[1,3]
    return scalar n_native_retry=`mat1'[1,4]
    return scalar n_native_recovered=`mat1'[1,5]
    return scalar n_simplex=`mat1'[1,6]
    return scalar n_native=`mat1'[1,7]
    return scalar n_simplex_retry=`mat1'[1,8]
    return scalar n_simplex_recovered=`mat1'[1,9]
    return scalar n_native_attempted=`mat1'[1,10]
    return scalar n_dual_attempted=`mat1'[1,11]
    return scalar n_dual_recovered=`mat1'[1,12]
    return scalar tolerance=`tolerance'
    return scalar acceptance_tolerance=min(`tolerance',1e-9)
    quietly summarize `generate' if `touse', meanonly
    return scalar mean_beta=r(mean)
    return scalar mean_eff=r(mean)
    return local model "undesirable"
    return local model_requested "undesirable"
    return local formulation "convex_joint_weak_disposal"
    return local disposability "weak"
    return local lp_backend "`lp_backend'"
    return local direction "`direction'"
    return local orientation "directional"
    return local rts "`rts'"
    return local engine "`engine'"
    return local inputs "`inputs'"
    return local outputs "`outputs'"
    return local bad "`bad'"
    return local generate "`generate'"
    if "`weakmatrix'"!="" {
        matrix colnames `mat10'=status backend simplex_retry native_retry eq_abs ineq_abs bound_abs eq_scaled ineq_scaled bound_scaled
        matrix colnames `mat11'=observation
        matrix colnames `mat12'=state upper_bound dual_solver_status dual_margin_lower normalization_lower corrections dual_attempts
        return matrix dual_recovery=`mat12'
        local matnames active idle intensity target_inputs target_outputs target_bad input_slack output_slack diagnostics obs_index
        forvalues j=2/11 {
            local k=`j'-1
            local name : word `k' of `matnames'
            if inlist(`j',5,8) matrix colnames `mat`j''=`inputs'
            if inlist(`j',6,9) matrix colnames `mat`j''=`outputs'
            if `j'==7 matrix colnames `mat`j''=`bad'
            return matrix `name'=`mat`j''
        }
    }
end
