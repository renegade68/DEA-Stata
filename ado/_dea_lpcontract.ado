*! version 0.4.9  05nov2026
*! Loader for the internal canonical LP contract
version 16.0

// Load the shared acceptance criterion before compiling solver functions.
capture mata: assert(_dea_lpcheck_ok(J(1,6,0),1e-8))
if _rc {
    quietly findfile _dea_lpcheck.ado
    quietly run "`r(fn)'"
}

capture mata: assert(_dea_native_retry_ready())
if _rc {
    quietly findfile _dea_native_retry.ado
    quietly run "`r(fn)'"
}

capture mata: assert(_dea_simplex_audit_ready()==16 & findexternal("_rsm_duplicate_ratio()")!=NULL & findexternal("_rsm_duplicate_map()")!=NULL & findexternal("_rsm_duplicate_problem()")!=NULL & findexternal("_rsm_duplicate_accept()")!=NULL & findexternal("_rsm_duplicate_retry()")!=NULL & findexternal("_rsm_two_phase_core()")!=NULL & findexternal("_rsm_solve_lp()")!=NULL & findexternal("_rsm_quality_fast()")!=NULL & findexternal("_rsm_affine_quality()")!=NULL & findexternal("_rsm_affine_refine()")!=NULL & findexternal("_rsm_success_fastread()")!=NULL & findexternal("_rsm_success_reset()")!=NULL & findexternal("_rsm_success_read()")!=NULL & findexternal("_rsm_success_counts_read()")!=NULL & findexternal("_rsm_success_check()")!=NULL & findexternal("_rsm_success_apply()")!=NULL & findexternal("_rsm_xgrow()")!=NULL & findexternal("_rsm_xadd()")!=NULL & findexternal("_rsm_xmul()")!=NULL & findexternal("_rsm_xnormalize()")!=NULL & findexternal("_rsm_rowmap()")!=NULL & findexternal("_rsm_row_empty()")!=NULL & findexternal("_rsm_row_accept()")!=NULL & findexternal("_rsm_row_retry()")!=NULL & findexternal("_dea_link_reset()")!=NULL & findexternal("_dea_link_read()")!=NULL & findexternal("_dea_link_store()")!=NULL & findexternal("_dea_link_small()")!=NULL & findexternal("_dea_link_residual()")!=NULL & findexternal("_dea_link_basis()")!=NULL & findexternal("_dea_link_apply()")!=NULL & findexternal("_dea_link_stage2()")!=NULL)
if _rc {
    capture findfile _dea_simplex.ado
    if _rc {
        di as error "_dea_lpcontract: _dea_simplex.ado not found"
        exit 601
    }
    quietly run "`r(fn)'"
}

capture mata: assert(_dea_simplex_audit_ready()==16 & findexternal("_rsm_duplicate_ratio()")!=NULL & findexternal("_rsm_duplicate_map()")!=NULL & findexternal("_rsm_duplicate_problem()")!=NULL & findexternal("_rsm_duplicate_accept()")!=NULL & findexternal("_rsm_duplicate_retry()")!=NULL & findexternal("_rsm_two_phase_core()")!=NULL & findexternal("_rsm_solve_lp()")!=NULL & findexternal("_rsm_quality_fast()")!=NULL & findexternal("_rsm_affine_quality()")!=NULL & findexternal("_rsm_affine_refine()")!=NULL & findexternal("_rsm_success_fastread()")!=NULL & findexternal("_rsm_success_reset()")!=NULL & findexternal("_rsm_success_read()")!=NULL & findexternal("_rsm_success_counts_read()")!=NULL & findexternal("_rsm_success_check()")!=NULL & findexternal("_rsm_success_apply()")!=NULL & findexternal("_rsm_xgrow()")!=NULL & findexternal("_rsm_xadd()")!=NULL & findexternal("_rsm_xmul()")!=NULL & findexternal("_rsm_xnormalize()")!=NULL & findexternal("_rsm_rowmap()")!=NULL & findexternal("_rsm_row_empty()")!=NULL & findexternal("_rsm_row_accept()")!=NULL & findexternal("_rsm_row_retry()")!=NULL & findexternal("_dea_link_reset()")!=NULL & findexternal("_dea_link_read()")!=NULL & findexternal("_dea_link_store()")!=NULL & findexternal("_dea_link_small()")!=NULL & findexternal("_dea_link_residual()")!=NULL & findexternal("_dea_link_basis()")!=NULL & findexternal("_dea_link_apply()")!=NULL & findexternal("_dea_link_stage2()")!=NULL)
if _rc {
    di as error "Incompatible simplex library; reinstall the complete matched DEA package."
    exit 498
}

capture findfile _dea_lpcontract_impl.do
if _rc {
    di as error "_dea_lpcontract: implementation file not found"
    exit 601
}
quietly run "`r(fn)'"

capture mata: assert(_dea_lpcontract_ready()==9 & findexternal("_rsm_public_success()")!=NULL & findexternal("_dea_lpcontract_simplex()")!=NULL & findexternal("_dea_lpcontract_dispatch()")!=NULL)
if _rc {
    di as error "Incompatible LP adapter; reinstall the complete matched DEA package."
    exit 498
}
