*! version 1.0.0 08nov2026
*! Load before replacing outputs; restore shared solver state on every load exit.
program define _dea_crossgame_load
    version 16.0
    capture mata: assert(_c57_context_ready()==57 & findexternal("_c57_context_save()")!=NULL & findexternal("_c57_context_restore()")!=NULL & findexternal("_c57_context_defaults()")!=NULL)
    if _rc {
        quietly findfile _dea_crosseff_context.ado
        quietly run "`r(fn)'"
    }
    tempname snapshot
    mata: `snapshot'=_c57_context_save()
    local reload=0
    capture noisily {
        capture mata: assert(_dea_nonte_lp_ready()==1 & findexternal("_dea_nonte_lp_ready()")!=NULL & findexternal("_dea_nonte_optimize()")!=NULL)
        if _rc {
            quietly findfile _dea_nonte_lp.ado
            quietly run "`r(fn)'"
            local reload=1
        }
        mata: assert(_dea_nonte_lp_ready()==1 & findexternal("_dea_nonte_lp_ready()")!=NULL & findexternal("_dea_nonte_optimize()")!=NULL)
        capture mata: assert(_dea_simplex_audit_ready()==16 & findexternal("_rsm_solve_lp()")!=NULL & findexternal("_dea_rec_abort()")!=NULL & findexternal("_rsm_audit_reset()")!=NULL & findexternal("_rsm_success_reset()")!=NULL & findexternal("_dea_link_reset()")!=NULL & findexternal("_dea_nonte_simplex()")!=NULL)
        if _rc {
            quietly findfile _dea_simplex.ado
            quietly run "`r(fn)'"
            local reload=1
        }
        mata: assert(_dea_simplex_audit_ready()==16 & findexternal("_rsm_solve_lp()")!=NULL & findexternal("_dea_rec_abort()")!=NULL & findexternal("_rsm_audit_reset()")!=NULL & findexternal("_rsm_success_reset()")!=NULL & findexternal("_dea_link_reset()")!=NULL & findexternal("_dea_nonte_simplex()")!=NULL)
        mata: _dea_rec_abort(); _rsm_audit_reset(0); _rsm_success_reset(1,0); _dea_link_reset(2,0)
        mata: _c57_context_defaults(`snapshot')
        capture mata: assert(_c70_objective_ready()==7001 & findexternal("_c70_objective_ready()")!=NULL & findexternal("_c70_scale()")!=NULL & findexternal("_c70_numeric()")!=NULL & findexternal("_c70_counts()")!=NULL & findexternal("_c70_detail()")!=NULL)
        if _rc {
            quietly findfile _dea_crosseff_objective.ado
            quietly run "`r(fn)'"
            local reload=1
        }
        mata: assert(_c70_objective_ready()==7001 & findexternal("_c70_objective_ready()")!=NULL & findexternal("_c70_scale()")!=NULL & findexternal("_c70_numeric()")!=NULL & findexternal("_c70_counts()")!=NULL & findexternal("_c70_detail()")!=NULL)
        capture mata: assert(_c58_cert_ready()==6204 & findexternal("_c58f_residual_sign()")!=NULL & findexternal("_c58f_quotient()")!=NULL & findexternal("_c58f_update_quotient()")!=NULL & findexternal("_c58f_restore()")!=NULL & findexternal("_c58f_rowcol()")!=NULL & findexternal("_c58f_split()")!=NULL & findexternal("_c58f_wide_product()")!=NULL & findexternal("_c58f_ratio()")!=NULL & findexternal("_c58f_identity_needs()")!=NULL & findexternal("_c58f_identity_cell()")!=NULL & findexternal("_c58_value()")!=NULL & findexternal("_c58_point()")!=NULL & findexternal("_c58_fit()")!=NULL & findexternal("_c58_map()")!=NULL & findexternal("_c60_certdiag()")!=NULL & findexternal("_c58_cert_ready()")!=NULL & findexternal("_c58d_double()")!=NULL & findexternal("_c58d_bits()")!=NULL & findexternal("_c58d_valid()")!=NULL & findexternal("_c58r_repair()")!=NULL & findexternal("_c58r_plan()")!=NULL & findexternal("_c58r_ratio_compare()")!=NULL & findexternal("_c58r_rowcol()")!=NULL & findexternal("_c58u_update()")!=NULL & findexternal("_c58b_solve()")!=NULL & findexternal("_c58b_certify()")!=NULL & findexternal("_c58b_pack()")!=NULL & findexternal("_c58b_pick()")!=NULL & findexternal("_c58i_rows()")!=NULL & findexternal("_c58b_shift()")!=NULL & findexternal("_c58b_decode()")!=NULL & findexternal("_c58b_div()")!=NULL & findexternal("_c58b_mul()")!=NULL & findexternal("_c58b_add()")!=NULL & findexternal("_c58b_neg()")!=NULL & findexternal("_c58b_cmp()")!=NULL & findexternal("_c58b_int()")!=NULL & findexternal("_c58b_tick()")!=NULL & findexternal("_c58b_norm()")!=NULL & findexternal("_C58_Cert()")!=NULL & findexternal("_C58_Work()")!=NULL & findexternal("_c61_pick()")!=NULL & findexternal("_c61_basis()")!=NULL & findexternal("_c61_active()")!=NULL & findexternal("_c61_recover()")!=NULL & findexternal("_c62_empty_trace()")!=NULL & findexternal("_c62_initial()")!=NULL & findexternal("_c62_ray()")!=NULL)
        if _rc {
            quietly findfile _dea_crosseff_cert.ado
            quietly run "`r(fn)'"
            local reload=1
        }
        mata: assert(_c58_cert_ready()==6204 & findexternal("_c58f_residual_sign()")!=NULL & findexternal("_c58f_quotient()")!=NULL & findexternal("_c58f_update_quotient()")!=NULL & findexternal("_c58f_restore()")!=NULL & findexternal("_c58f_rowcol()")!=NULL & findexternal("_c58f_split()")!=NULL & findexternal("_c58f_wide_product()")!=NULL & findexternal("_c58f_ratio()")!=NULL & findexternal("_c58f_identity_needs()")!=NULL & findexternal("_c58f_identity_cell()")!=NULL & findexternal("_c58_value()")!=NULL & findexternal("_c58_point()")!=NULL & findexternal("_c58_fit()")!=NULL & findexternal("_c58_map()")!=NULL & findexternal("_c60_certdiag()")!=NULL & findexternal("_c58_cert_ready()")!=NULL & findexternal("_c58d_double()")!=NULL & findexternal("_c58d_bits()")!=NULL & findexternal("_c58d_valid()")!=NULL & findexternal("_c58r_repair()")!=NULL & findexternal("_c58r_plan()")!=NULL & findexternal("_c58r_ratio_compare()")!=NULL & findexternal("_c58r_rowcol()")!=NULL & findexternal("_c58u_update()")!=NULL & findexternal("_c58b_solve()")!=NULL & findexternal("_c58b_certify()")!=NULL & findexternal("_c58b_pack()")!=NULL & findexternal("_c58b_pick()")!=NULL & findexternal("_c58i_rows()")!=NULL & findexternal("_c58b_shift()")!=NULL & findexternal("_c58b_decode()")!=NULL & findexternal("_c58b_div()")!=NULL & findexternal("_c58b_mul()")!=NULL & findexternal("_c58b_add()")!=NULL & findexternal("_c58b_neg()")!=NULL & findexternal("_c58b_cmp()")!=NULL & findexternal("_c58b_int()")!=NULL & findexternal("_c58b_tick()")!=NULL & findexternal("_c58b_norm()")!=NULL & findexternal("_C58_Cert()")!=NULL & findexternal("_C58_Work()")!=NULL & findexternal("_c61_pick()")!=NULL & findexternal("_c61_basis()")!=NULL & findexternal("_c61_active()")!=NULL & findexternal("_c61_recover()")!=NULL & findexternal("_c62_empty_trace()")!=NULL & findexternal("_c62_initial()")!=NULL & findexternal("_c62_ray()")!=NULL)
        capture mata: assert(_c70_exact_ready()==7003 & findexternal("_c70_exact_ready()")!=NULL & findexternal("_c70_decode()")!=NULL & findexternal("_c70_build()")!=NULL & findexternal("_c70_value()")!=NULL & findexternal("_c70_fit()")!=NULL & findexternal("_c70b_certify()")!=NULL & findexternal("_c70b_solve()")!=NULL & findexternal("_c70u_update()")!=NULL & findexternal("_c70r_repair()")!=NULL & findexternal("_c70_recover()")!=NULL & findexternal("_c70_ray()")!=NULL)
        if _rc {
            quietly findfile _dea_crosseff_exactobjective.ado
            quietly run "`r(fn)'"
            local reload=1
        }
        mata: assert(_c70_exact_ready()==7003 & findexternal("_c70_exact_ready()")!=NULL & findexternal("_c70_decode()")!=NULL & findexternal("_c70_build()")!=NULL & findexternal("_c70_value()")!=NULL & findexternal("_c70_fit()")!=NULL & findexternal("_c70b_certify()")!=NULL & findexternal("_c70b_solve()")!=NULL & findexternal("_c70u_update()")!=NULL & findexternal("_c70r_repair()")!=NULL & findexternal("_c70_recover()")!=NULL & findexternal("_c70_ray()")!=NULL)
        capture mata: assert(_c107_ready()==10702 & findexternal("_C107_Result()")!=NULL & findexternal("_C107_Cost()")!=NULL & findexternal("_c108_needed()")!=NULL & findexternal("_c108_join()")!=NULL & findexternal("_c107_ready()")!=NULL & findexternal("_c107_valid()")!=NULL & findexternal("_c107_gcd()")!=NULL & findexternal("_c107_div()")!=NULL & findexternal("_c107_cost()")!=NULL & findexternal("_c107_attempt()")!=NULL)
        if _rc {
            quietly findfile _dea_crosseff_costunits.ado
            quietly run "`r(fn)'"
            local reload=1
        }
        mata: assert(_c107_ready()==10702 & findexternal("_C107_Result()")!=NULL & findexternal("_C107_Cost()")!=NULL & findexternal("_c108_needed()")!=NULL & findexternal("_c108_join()")!=NULL & findexternal("_c107_ready()")!=NULL & findexternal("_c107_valid()")!=NULL & findexternal("_c107_gcd()")!=NULL & findexternal("_c107_div()")!=NULL & findexternal("_c107_cost()")!=NULL & findexternal("_c107_attempt()")!=NULL)
        capture mata: assert(_c84_ready()==8404 & findexternal("_c84_work()")!=NULL & findexternal("_c84_log()")!=NULL & findexternal("_c84_scale()")!=NULL & findexternal("_c84_coefficient()")!=NULL & findexternal("_c84_build()")!=NULL & findexternal("_c84_dyad()")!=NULL & findexternal("_c84_number()")!=NULL & findexternal("_c84_abs()")!=NULL & findexternal("_c84_sub()")!=NULL & findexternal("_c84_div()")!=NULL & findexternal("_c84_gt()")!=NULL & findexternal("_c84_value()")!=NULL & findexternal("_c84_dot()")!=NULL & findexternal("_c84_objective()")!=NULL & findexternal("_c84_point()")!=NULL & findexternal("_c84_guard()")!=NULL & findexternal("_c84_ready()")!=NULL & findexternal("_c84_mul()")!=NULL & findexternal("_c84_peer()")!=NULL & findexternal("_c84_multiplier()")!=NULL & findexternal("_C84_Model()")!=NULL & findexternal("_C84_Q()")!=NULL)
        if _rc {
            quietly findfile _dea_crosseff_raw.ado
            quietly run "`r(fn)'"
            local reload=1
        }
        mata: assert(_c84_ready()==8404 & findexternal("_c84_work()")!=NULL & findexternal("_c84_log()")!=NULL & findexternal("_c84_scale()")!=NULL & findexternal("_c84_coefficient()")!=NULL & findexternal("_c84_build()")!=NULL & findexternal("_c84_dyad()")!=NULL & findexternal("_c84_number()")!=NULL & findexternal("_c84_abs()")!=NULL & findexternal("_c84_sub()")!=NULL & findexternal("_c84_div()")!=NULL & findexternal("_c84_gt()")!=NULL & findexternal("_c84_value()")!=NULL & findexternal("_c84_dot()")!=NULL & findexternal("_c84_objective()")!=NULL & findexternal("_c84_point()")!=NULL & findexternal("_c84_guard()")!=NULL & findexternal("_c84_ready()")!=NULL & findexternal("_c84_mul()")!=NULL & findexternal("_c84_peer()")!=NULL & findexternal("_c84_multiplier()")!=NULL & findexternal("_C84_Model()")!=NULL & findexternal("_C84_Q()")!=NULL)
        capture mata: assert(_g114_core_ready()==11401 & findexternal("_dea_cg_lp()")!=NULL & findexternal("_dea_crossgame_compute()")!=NULL & findexternal("_dea_cg_peers()")!=NULL & findexternal("_dea_cg_average()")!=NULL & findexternal("_dea_cg_show()")!=NULL)
        if _rc {
            quietly findfile _dea_crossgame_core.ado
            quietly run "`r(fn)'"
            local reload=1
        }
        mata: assert(_g114_core_ready()==11401)
        capture mata: assert(_g114_ready()==11401 & findexternal("_g114_opt()")!=NULL & findexternal("_g114_optc()")!=NULL & findexternal("_g114_guard()")!=NULL & findexternal("_g114_recover()")!=NULL & findexternal("_g114_dispatch()")!=NULL & findexternal("_g114_report()")!=NULL)
        if _rc | `reload' {
            quietly findfile _dea_crossgame_recovery.ado
            quietly run "`r(fn)'"
        }
        mata: assert(_g114_ready()==11401 & findexternal("_g114_opt()")!=NULL & findexternal("_g114_optc()")!=NULL & findexternal("_g114_guard()")!=NULL & findexternal("_g114_recover()")!=NULL & findexternal("_g114_dispatch()")!=NULL & findexternal("_g114_report()")!=NULL)
    }
    local rc=_rc
    mata: _c57_context_restore(`snapshot')
    mata: mata drop `snapshot'
    if `rc' exit `rc'
end
