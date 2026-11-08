*! version 1.10.3 08nov2026
*! Peer-only objectives, exact intended-sum certification and build diagnostics.
*! Cross-efficiency: native IPM default and bounded explicit simplex certificates.
*! 2026-11-02: explicit input VRS Lim-Zhu peer formula; preserve legacy default.
*! 2026-11-02: scale-free nonnegative weights, validated LPs, and stage diagnostics.
*! 2026-11-01: accept the minimum option abbreviations advertised in help.
*! 2026-11-01: use shared Form A/B parser; preserve if/in expressions.
*! 2026-11-01: preflight output names; explicit replacevars; protect model data.
*! dea_crosseff.ado — Cross-Efficiency DEA
*! Requires: Stata 16+ (LinearProgram)
*!
*! AUTHORS:
*! --------
*!   Choonjoo Lee (corresponding author)
*!     Korea National Defense University
*!     Seoul, Korea
*!     bloom.rampike@gmail.com
*!
*! REFERENCES:
*! -----------
*!   Doyle, J., & Green, R. (1994). Efficiency and cross-efficiency in
*!     DEA: Derivations, meanings and uses. Journal of the Operational
*!     Research Society, 45(5), 567-578. doi:10.1057/jors.1994.84
*!
*!   Lim, S., & Zhu, J. (2015). DEA cross-efficiency evaluation under
*!     variable returns to scale. Journal of the Operational Research
*!     Society, 66(3), 476-487. doi:10.1057/jors.2014.13
*!
*!   Sexton, T. R., Silkman, R. H., & Hogan, A. J. (1986). Data
*!     envelopment analysis: Critique and extensions. New Directions
*!     for Program Evaluation, 1986(32), 73-105. doi:10.1002/ev.1441
*!
*! Changelog:
*!   1.4.0  25apr2026  Implement method(benevolent/aggressive): 2-stage LP
*!   1.3.1  25apr2026  VRS CE peer-eval: use u*y/v*x only (no w), avoids negative CE
*!   1.3.0  25apr2026  Add VRS and output orientation support
*!   1.2.2  25apr2026  Fix string DMU: remove [varlist] from syntax (Stata defaults to
*!                     all dataset vars causing r(109)); Form A transform puts outputs
*!                     into outputs() option; markout excludes dmu variable
*!   1.1.0  11feb2026  Unified syntax: supports both
*!                        dea_crosseff x = y, [options]  (Form A)
*!                        dea_crosseff, inputs(x) outputs(y) [options] (Form B)
*!   1.10.2 08nov2026 Exact dyadic factoring within existing certificate arithmetic limits.
*!   1.10.1 08nov2026 Scoped finite decoder for failed raw primaries; matched helper ABI.
*!   1.0.0  06feb2026  Initial packaged version; based on earlier teaching/seminar versions

/**
 * CROSS-EFFICIENCY DEA
 * ====================
 *
 * Cross-efficiency evaluates each DMU using the optimal weights
 * of all other DMUs, providing a peer evaluation mechanism.
 *
 * Key concepts:
 *   - Self-appraisal: DMU evaluated with its own optimal weights
 *   - Peer-appraisal: DMU evaluated with other DMUs' weights
 *   - Cross-efficiency score: Average of peer appraisals
 *
 * Advantages:
 *   - Provides peer rankings; ties and alternative optimal weights can remain
 *   - Reduces unrealistic weight flexibility
 *   - Identifies "maverick" DMUs with atypical weights
 *
 * Methods:
 *   - Aggressive: Secondary goal minimizes cross-efficiency of others
 *   - Benevolent: Secondary goal maximizes cross-efficiency of others
 *   - Arbitrary: Uses primary optimal weights only
 */

capture program drop dea_crosseff
program define dea_crosseff, rclass
    version 16.0

    // =========================================================================
    // Unified Syntax Parser
    //   Form A: dea_crosseff x1 x2 = y1 y2 [if] [in], [options]
    //     → transforms to: , inputs(x1 x2) outputs(y1 y2) [options]
    //   Form B: dea_crosseff [if] [in], inputs() outputs() [options]
    //
    // No [varlist] in syntax: avoids r(109) from Stata defaulting to all
    // dataset vars (including string dmu) when no varlist is specified.
    // =========================================================================
    _dea_parse_io `0'
    local 0 `"`r(args)'"'

    syntax [if] [in], ///
        [DMU(varname)           ///
         Inputs(varlist numeric) ///
         Outputs(varlist numeric) ///
         Rts(string) ///
         ORT(string) ORientation(string) ///
         METHOD(string) ///
         ENGine(string) ///
         VRSCross(string) ///
         GENerate(string) ///
         REPLACEVARS      ///
         TOLerance(real 1e-7) ///
         MATrix ///
         DETail ///
         DEBUG OPTIMAdiag CERange]

    if "`cerange'"!="" local optimadiag "optimadiag"

    if missing(`tolerance') | `tolerance'<=0 {
        di as error "tolerance() must be positive and finite"
        exit 198
    }
    capture mata: assert(_dea_nonte_lp_ready()==1)
    if _rc {
        quietly findfile _dea_nonte_lp.ado
        quietly run "`r(fn)'"
    }

    // Validate required options (Form A already placed both into Inputs/Outputs options)
    if "`outputs'" == "" {
        di as error "outputs() required, or use: dea_crosseff inputs = outputs, [options]"
        exit 198
    }
    if "`inputs'" == "" {
        di as error "inputs() required, or use: dea_crosseff inputs = outputs, [options]"
        exit 198
    }

    // Validate dmu is a real variable (if provided)
    if "`dmu'" != "" {
        capture confirm variable `dmu'
        if _rc != 0 {
            di as error "dmu(`dmu') is not a variable in the dataset"
            exit 111
        }
    }

    // Reconcile ort() / orientation()
    if "`ort'" != "" & "`orientation'" != "" {
        di as error "Specify only one of: ort() or orientation()"
        exit 198
    }
    if "`ort'" != "" local orientation "`ort'"

    // Defaults
    if "`rts'" == "" local rts "crs"
    if "`orientation'" == "" local orientation "input"
    if "`method'" == "" local method "arbitrary"
    if "`vrscross'" == "" local vrscross "legacy"
    if "`generate'" == "" local generate "_ce_eff"

    // Abbreviations
    if inlist("`orientation'", "i", "in")  local orientation "input"
    if inlist("`orientation'", "o", "out") local orientation "output"
    if "`rts'" == "c" local rts "crs"
    if "`rts'" == "v" local rts "vrs"

    // Validate
    if !inlist("`rts'", "crs", "vrs") {
        di as error "rts() must be crs or vrs"
        exit 198
    }
    if !inlist("`orientation'", "input", "output") {
        di as error "orientation() must be input or output"
        exit 198
    }
    if !inlist("`method'", "arbitrary", "benevolent", "aggressive") {
        di as error "method() must be arbitrary, benevolent, or aggressive"
        exit 198
    }
    if !inlist("`vrscross'", "legacy", "limzhu") {
        di as error "vrscross() must be legacy or limzhu"
        exit 198
    }
    if "`vrscross'"=="limzhu" & ("`rts'"!="vrs" | "`orientation'"!="input") {
        di as error "vrscross(limzhu) requires rts(vrs) and input orientation"
        exit 198
    }

    if "`engine'"=="" local engine "ipm"
    local engine=lower("`engine'")
    if !inlist("`engine'","ipm","simplex") {
        di as error "engine() must be ipm or simplex"
        exit 198
    }

    // Validate every output before any existing variable can be replaced.
    _dea_output_check, names(`generate' `generate'_self `generate'_maverick) ///
        protect(`inputs' `outputs' `dmu') `replacevars'
    local _replace_outputs "`r(existing)'"

    // Mark sample — do NOT include dmu in markout (string DMU causes r(109))
    marksample touse
    markout `touse' `inputs' `outputs'
    if "`dmu'" != "" {
        qui replace `touse' = 0 if missing(`dmu')   // handles both string("") and numeric(.)
    }

    // Count
    qui count if `touse'
    local n = r(N)
    if `n'==0 error 2000
    foreach _side in inputs outputs {
        tempvar _positive_side
        quietly generate byte `_positive_side'=0 if `touse'
        foreach _v of local `_side' {
            quietly count if `touse' & `_v'<0
            if r(N) {
                di as error "cross-efficiency quantities must be nonnegative"
                exit 459
            }
            quietly replace `_positive_side'=1 if `touse' & `_v'>0
        }
        quietly count if `touse' & `_positive_side'==0
        if r(N) {
            di as error "each selected observation requires a positive input and output"
            exit 459
        }
    }

    capture mata: assert(_c70_objective_ready()==7001 & findexternal("_c70_scale()")!=NULL & findexternal("_c70_numeric()")!=NULL & findexternal("_c70_counts()")!=NULL & findexternal("_c70_detail()")!=NULL)
    if _rc {
        quietly findfile _dea_crosseff_objective.ado
        quietly run "`r(fn)'"
    }
    mata: assert(_c70_objective_ready()==7001 & findexternal("_c70_scale()")!=NULL & findexternal("_c70_numeric()")!=NULL & findexternal("_c70_counts()")!=NULL & findexternal("_c70_detail()")!=NULL)

    // Read-only global-scaling diagnostics; do not pass these copies to any solver.
    capture mata: assert(_c60_diag_ready()==7001 & findexternal("_c60_prepare()")!=NULL & findexternal("_c60_loss()")!=NULL & findexternal("_c60_computation()")!=NULL & findexternal("_c60_empty_cert()")!=NULL & findexternal("_c60_loss_detail()")!=NULL & findexternal("_c60_cert_detail()")!=NULL & findexternal("_c60_stop_label()")!=NULL & findexternal("_c62_empty_recovery()")!=NULL & findexternal("_c62_recovery_detail()")!=NULL & findexternal("_c66_native_detail()")!=NULL)
    if _rc {
        quietly findfile _dea_crosseff_diagnostics.ado
        quietly run "`r(fn)'"
    }
    mata: assert(_c60_diag_ready()==7001 & findexternal("_c60_prepare()")!=NULL & findexternal("_c60_loss()")!=NULL & findexternal("_c60_computation()")!=NULL & findexternal("_c60_empty_cert()")!=NULL & findexternal("_c60_loss_detail()")!=NULL & findexternal("_c60_cert_detail()")!=NULL & findexternal("_c60_stop_label()")!=NULL & findexternal("_c62_empty_recovery()")!=NULL & findexternal("_c62_recovery_detail()")!=NULL & findexternal("_c66_native_detail()")!=NULL)
    tempname norm_counts norm_inputs norm_outputs norm_obs computation_diag cert_p_diag cert_s_diag cert_map_diag recovery_p_diag recovery_s_diag
    mata: _c60_prepare("`inputs'","`outputs'","`touse'","`matrix'"!="","`norm_inputs'","`norm_outputs'","`norm_obs'","`norm_counts'")
    local norm_lost=`norm_counts'[1,1]
    local norm_in=`norm_counts'[1,2]
    local norm_out=`norm_counts'[1,3]
    local norm_obs_count=`norm_counts'[1,4]

    local _optima_supported=0
    if "`optimadiag'"!="" {
        local _optima_m : word count `inputs'
        local _optima_s : word count `outputs'
        local _optima_supported=(`n'<=32 & `_optima_m'+`_optima_s'+2*("`rts'"=="vrs")<=16 & "`vrscross'"=="legacy")
    }

    // Load before replacing any outputs. Preserve shared solver diagnostics,
    // including on a failed helper load. A cold native call never loads RSM.
    tempname solver_snapshot
    if "`engine'"=="simplex" | "`optimadiag'"!="" {
        capture mata: assert(_c57_context_ready()==57)
        if _rc {
            quietly findfile _dea_crosseff_context.ado
            quietly run "`r(fn)'"
        }
        mata: `solver_snapshot'=_c57_context_save()
        capture noisily {
            capture mata: assert(_dea_simplex_audit_ready()==16 & findexternal("_rsm_solve_lp()")!=NULL)
            if _rc {
                quietly findfile _dea_simplex.ado
                quietly run "`r(fn)'"
            }
            mata: assert(_dea_simplex_audit_ready()==16)
            mata: _dea_rec_abort(); _rsm_audit_reset(0); _rsm_success_reset(1,0); _dea_link_reset(2,0)
            mata: _c57_context_defaults(`solver_snapshot')
            local reload_exact=0
            capture mata: assert(_c58_cert_ready()==6203 & findexternal("_c58_fit()")!=NULL & findexternal("_c58b_decode()")!=NULL & findexternal("_c60_certdiag()")!=NULL & findexternal("_c61_recover()")!=NULL & findexternal("_c61_pick()")!=NULL & findexternal("_c61_basis()")!=NULL & findexternal("_c61_active()")!=NULL & findexternal("_c62_ray()")!=NULL & findexternal("_c62_initial()")!=NULL & findexternal("_c62_empty_trace()")!=NULL & findexternal("_c58f_split()")!=NULL & findexternal("_c58f_wide_product()")!=NULL & findexternal("_c58f_ratio()")!=NULL & findexternal("_c58f_identity_needs()")!=NULL & findexternal("_c58f_identity_cell()")!=NULL & findexternal("_c58f_restore()")!=NULL & findexternal("_c58f_rowcol()")!=NULL & findexternal("_c58f_quotient()")!=NULL & findexternal("_c58f_update_quotient()")!=NULL & findexternal("_c58f_residual_sign()")!=NULL)
            if _rc {
                capture mata: mata drop _c58_multiplier()
                capture mata: mata drop _c57_compute()
                quietly findfile _dea_crosseff_cert.ado
                quietly run "`r(fn)'"
                local reload_exact=1
            }
            mata: assert(_c58_cert_ready()==6203 & findexternal("_c58_fit()")!=NULL & findexternal("_c58b_decode()")!=NULL & findexternal("_c60_certdiag()")!=NULL & findexternal("_c61_recover()")!=NULL & findexternal("_c61_pick()")!=NULL & findexternal("_c61_basis()")!=NULL & findexternal("_c61_active()")!=NULL & findexternal("_c62_ray()")!=NULL & findexternal("_c62_initial()")!=NULL & findexternal("_c62_empty_trace()")!=NULL & findexternal("_c58f_split()")!=NULL & findexternal("_c58f_wide_product()")!=NULL & findexternal("_c58f_ratio()")!=NULL & findexternal("_c58f_identity_needs()")!=NULL & findexternal("_c58f_identity_cell()")!=NULL & findexternal("_c58f_restore()")!=NULL & findexternal("_c58f_rowcol()")!=NULL & findexternal("_c58f_quotient()")!=NULL & findexternal("_c58f_update_quotient()")!=NULL & findexternal("_c58f_residual_sign()")!=NULL)
            capture mata: assert(_c70_exact_ready()==7002 & findexternal("_c70_build()")!=NULL & findexternal("_c70_fit()")!=NULL & findexternal("_c70_value()")!=NULL & findexternal("_c70_decode()")!=NULL & findexternal("_c70b_certify()")!=NULL & findexternal("_c70b_solve()")!=NULL & findexternal("_c70u_update()")!=NULL & findexternal("_c70r_repair()")!=NULL & findexternal("_c70_recover()")!=NULL & findexternal("_c70_ray()")!=NULL)
            if _rc | `reload_exact' {
                quietly findfile _dea_crosseff_exactobjective.ado
                quietly run "`r(fn)'"
            }
            mata: assert(_c70_exact_ready()==7002 & findexternal("_c70_build()")!=NULL & findexternal("_c70_fit()")!=NULL & findexternal("_c70_value()")!=NULL & findexternal("_c70_decode()")!=NULL & findexternal("_c70b_certify()")!=NULL & findexternal("_c70b_solve()")!=NULL & findexternal("_c70u_update()")!=NULL & findexternal("_c70r_repair()")!=NULL & findexternal("_c70_recover()")!=NULL & findexternal("_c70_ray()")!=NULL)
            capture mata: assert(_c84_ready()==8403 & findexternal("_c84_work()")!=NULL & findexternal("_c84_log()")!=NULL & findexternal("_c84_scale()")!=NULL & findexternal("_c84_coefficient()")!=NULL & findexternal("_c84_build()")!=NULL & findexternal("_c84_dyad()")!=NULL & findexternal("_c84_number()")!=NULL & findexternal("_c84_abs()")!=NULL & findexternal("_c84_sub()")!=NULL & findexternal("_c84_div()")!=NULL & findexternal("_c84_gt()")!=NULL & findexternal("_c84_value()")!=NULL & findexternal("_c84_dot()")!=NULL & findexternal("_c84_objective()")!=NULL & findexternal("_c84_point()")!=NULL & findexternal("_c84_guard()")!=NULL & findexternal("_c84_ready()")!=NULL & findexternal("_c84_mul()")!=NULL & findexternal("_c84_peer()")!=NULL & findexternal("_c84_multiplier()")!=NULL & findexternal("_C84_Model()")!=NULL & findexternal("_C84_Q()")!=NULL)
            if _rc | `reload_exact' {
                quietly findfile _dea_crosseff_raw.ado
                quietly run "`r(fn)'"
            }
            mata: assert(_c84_ready()==8403 & findexternal("_c84_work()")!=NULL & findexternal("_c84_log()")!=NULL & findexternal("_c84_scale()")!=NULL & findexternal("_c84_coefficient()")!=NULL & findexternal("_c84_build()")!=NULL & findexternal("_c84_dyad()")!=NULL & findexternal("_c84_number()")!=NULL & findexternal("_c84_abs()")!=NULL & findexternal("_c84_sub()")!=NULL & findexternal("_c84_div()")!=NULL & findexternal("_c84_gt()")!=NULL & findexternal("_c84_value()")!=NULL & findexternal("_c84_dot()")!=NULL & findexternal("_c84_objective()")!=NULL & findexternal("_c84_point()")!=NULL & findexternal("_c84_guard()")!=NULL & findexternal("_c84_ready()")!=NULL & findexternal("_c84_mul()")!=NULL & findexternal("_c84_peer()")!=NULL & findexternal("_c84_multiplier()")!=NULL & findexternal("_C84_Model()")!=NULL & findexternal("_C84_Q()")!=NULL)
            capture mata: assert(_c57_simplex_ready()==8402 & findexternal("_c57_multiplier()")!=NULL & findexternal("_c58_multiplier()")!=NULL & findexternal("_c57_compute()")!=NULL)
            if _rc {
                quietly findfile _dea_crosseff_simplex.ado
                quietly run "`r(fn)'"
            }
            mata: assert(_c57_simplex_ready()==8402 & findexternal("_c57_multiplier()")!=NULL & findexternal("_c58_multiplier()")!=NULL & findexternal("_c57_compute()")!=NULL)
            if "`cerange'"!="" {
                capture mata: assert(_c76_ready()==7601 & _c77_ready()==7701 & findexternal("_c76_ready()")!=NULL & findexternal("_c76_budget()")!=NULL & findexternal("_c76_begin()")!=NULL & findexternal("_c76_allow()")!=NULL & findexternal("_c76_charge()")!=NULL & findexternal("_c76_work()")!=NULL & findexternal("_c76_empty()")!=NULL & findexternal("_c76_endempty()")!=NULL & findexternal("_c76_result()")!=NULL & findexternal("_c76_certpoint()")!=NULL & findexternal("_c76_pack()")!=NULL & findexternal("_c76_dot()")!=NULL & findexternal("_c76_cmp()")!=NULL & findexternal("_c76_check()")!=NULL & findexternal("_c76_positive()")!=NULL & findexternal("_c76_homogenize()")!=NULL & findexternal("_c76_dehomogenize()")!=NULL & findexternal("_c76_display()")!=NULL & findexternal("_c76_solve()")!=NULL & findexternal("_c76_endpoint()")!=NULL & findexternal("_c76_cell()")!=NULL & findexternal("_c76_coeff()")!=NULL & findexternal("_c76_original()")!=NULL & findexternal("_c76_peer()")!=NULL & findexternal("_C76_Budget()")!=NULL & findexternal("_C76_Result()")!=NULL & findexternal("_C76_End()")!=NULL & findexternal("_C76_Point()")!=NULL & findexternal("_c77_ready()")!=NULL & findexternal("_c77_start()")!=NULL & findexternal("_c77_begin()")!=NULL & findexternal("_c77_face()")!=NULL & findexternal("_c77_finish()")!=NULL)
                if _rc | `reload_exact' {
                    quietly findfile _dea_crosseff_cerange.ado
                    quietly run "`r(fn)'"
                }
                mata: assert(_c76_ready()==7601 & _c77_ready()==7701 & findexternal("_c76_ready()")!=NULL & findexternal("_c76_budget()")!=NULL & findexternal("_c76_begin()")!=NULL & findexternal("_c76_allow()")!=NULL & findexternal("_c76_charge()")!=NULL & findexternal("_c76_work()")!=NULL & findexternal("_c76_empty()")!=NULL & findexternal("_c76_endempty()")!=NULL & findexternal("_c76_result()")!=NULL & findexternal("_c76_certpoint()")!=NULL & findexternal("_c76_pack()")!=NULL & findexternal("_c76_dot()")!=NULL & findexternal("_c76_cmp()")!=NULL & findexternal("_c76_check()")!=NULL & findexternal("_c76_positive()")!=NULL & findexternal("_c76_homogenize()")!=NULL & findexternal("_c76_dehomogenize()")!=NULL & findexternal("_c76_display()")!=NULL & findexternal("_c76_solve()")!=NULL & findexternal("_c76_endpoint()")!=NULL & findexternal("_c76_cell()")!=NULL & findexternal("_c76_coeff()")!=NULL & findexternal("_c76_original()")!=NULL & findexternal("_c76_peer()")!=NULL & findexternal("_C76_Budget()")!=NULL & findexternal("_C76_Result()")!=NULL & findexternal("_C76_End()")!=NULL & findexternal("_C76_Point()")!=NULL & findexternal("_c77_ready()")!=NULL & findexternal("_c77_start()")!=NULL & findexternal("_c77_begin()")!=NULL & findexternal("_c77_face()")!=NULL & findexternal("_c77_finish()")!=NULL)
            }
            if "`optimadiag'"!="" {
                capture mata: assert(_c73_ready()==7301 & _c77_optima_ready()==7701 & findexternal("_c73_limits()")!=NULL & findexternal("_c73_allow()")!=NULL & findexternal("_c73_charge()")!=NULL & findexternal("_c73_solve()")!=NULL & findexternal("_c73_second()")!=NULL & findexternal("_c73_fraction()")!=NULL & findexternal("_c73_range()")!=NULL & findexternal("_c73_link()")!=NULL & findexternal("_c73_face()")!=NULL & findexternal("_c73_compute()")!=NULL)
                if _rc | `reload_exact' {
                    quietly findfile _dea_crosseff_optima.ado
                    quietly run "`r(fn)'"
                }
                mata: assert(_c73_ready()==7301 & _c77_optima_ready()==7701 & findexternal("_c73_limits()")!=NULL & findexternal("_c73_allow()")!=NULL & findexternal("_c73_charge()")!=NULL & findexternal("_c73_solve()")!=NULL & findexternal("_c73_second()")!=NULL & findexternal("_c73_fraction()")!=NULL & findexternal("_c73_range()")!=NULL & findexternal("_c73_link()")!=NULL & findexternal("_c73_face()")!=NULL & findexternal("_c73_compute()")!=NULL)
            }
        }
        local load_rc=_rc
        if `load_rc' {
            mata: _c57_context_restore(`solver_snapshot')
            mata: mata drop `solver_snapshot'
            exit `load_rc'
        }
    }

    if "`engine'"=="simplex" | "`optimadiag'"!="" {
        mata: _c57_context_restore(`solver_snapshot')
        mata: mata drop `solver_snapshot'
    }

    // Display header
    di as text ""
    di as text "{hline 70}"
    di as text "Cross-Efficiency DEA"
    di as text "{hline 70}"
    di as text "Method:        " as result "`method'"
    di as text "RTS:           " as result "`rts'"
    di as text "Orientation:   " as result "`orientation'"
    if "`rts'"=="vrs" di as text "VRS peer rule: " as result "`vrscross'"
    di as text "DMUs:          " as result "`n'"
    di as text "Inputs:        " as result "`inputs'"
    di as text "Outputs:       " as result "`outputs'"
    di as text "{hline 70}"

    // The sample may refer to an existing output; mark it before replacement.
    foreach _ov of local _replace_outputs {
        quietly drop `_ov'
    }

    // Create output variables
    qui gen double `generate' = . if `touse'
    qui gen double `generate'_self = . if `touse'
    qui gen double `generate'_maverick = . if `touse'

    // Create DMU identifier if not provided
    if "`dmu'" == "" {
        tempvar dmu_temp
        qui gen `dmu_temp' = _n if `touse'
        local dmu "`dmu_temp'"
    }

    // Call Mata function.  The tempname keeps the optional matrix alive long
    // enough for the public r(ce_matrix) contract to be restored below.
    tempname ce_matrix_return ce_v_return ce_u_return ce_w_return ce_status_return ce_peers_return
    tempname ce_diag_return ce_link_return native_numeric_return objective_return
    tempname raw_recovery_return raw_scaling_return raw_peer_return
    capture noisily {
    if "`engine'"=="simplex" {
        mata: `solver_snapshot'=_c57_context_save()
        mata: _dea_rec_abort(); _rsm_audit_reset(0); _rsm_success_reset(1,0); _dea_link_reset(2,0)
    }
    if "`engine'"=="simplex" {
    mata: _c57_compute("`inputs'", "`outputs'", ///
        "`dmu'", "`touse'", ///
        "`generate'", "`generate'_self", "`generate'_maverick", ///
        "`method'", "`rts'", "`orientation'", "`vrscross'", `tolerance', ///
        "`matrix'" != "", "`detail'" != "", "`debug'" != "", ///
        "`ce_matrix_return'", "`ce_v_return'", "`ce_u_return'", ///
        "`ce_w_return'", "`ce_status_return'", "`ce_peers_return'", "`ce_diag_return'", "`ce_link_return'", "`cert_p_diag'", "`cert_s_diag'", "`cert_map_diag'", "`recovery_p_diag'", "`recovery_s_diag'", "`objective_return'", "`raw_recovery_return'", "`raw_scaling_return'", "`raw_peer_return'", `_optima_supported')
    }
    else {
    mata: _dea_crosseff_compute("`inputs'", "`outputs'", ///
        "`dmu'", "`touse'", ///
        "`generate'", "`generate'_self", "`generate'_maverick", ///
        "`method'", "`rts'", "`orientation'", "`vrscross'", `tolerance', ///
        "`matrix'" != "", "`detail'" != "", "`debug'" != "", ///
        "`ce_matrix_return'", "`ce_v_return'", "`ce_u_return'", ///
        "`ce_w_return'", "`ce_status_return'", "`ce_peers_return'", "`native_numeric_return'", "`objective_return'", `_optima_supported')
    }
    }
    local compute_rc=_rc
    if "`engine'"=="simplex" {
        mata: _c57_context_restore(`solver_snapshot')
        mata: mata drop `solver_snapshot'
    }
    if `compute_rc' exit `compute_rc'
    foreach _k in attempted build_failed numeric_rounded numeric_lost exact_verified warning {
        local objective_`_k'=r(n_objective_`_k')
    }
    foreach _k in attempted recovered unavailable primary_unknown export_obstructed secondary_unknown secondary_unbounded peers_evaluated peers_unavailable peers_range peers_arithmetic {
        local raw_`_k'=0
        if "`engine'"=="simplex" local raw_`_k'=r(n_raw_`_k')
    }
    local _failed_primary=r(n_failed_primary)
    local _failed_secondary=r(n_failed_secondary)
    local _failed_ce=r(n_failed_ce)
    local _missing_pairs=r(n_missing_peer_pairs)
    local _missing_diagonal=r(n_missing_diagonal)
    local native_bad_primary=0
    local native_bad_secondary=0
    if "`engine'"=="ipm" {
        local native_bad_primary=r(n_native_nonfinite_primary)
        local native_bad_secondary=r(n_native_nonfinite_secondary)
    }

    local n_cert_primary_attempted=0
    local n_cert_primary_verified=0
    local n_cert_secondary_attempted=0
    local n_cert_secondary_verified=0
    local n_cert_repair_eq_stops=0
    local n_cert_repair_both_stops=0
    local n_cert_mapping_failed=0
    local n_recovery_requested=0
    local n_recovery_finite=0
    local n_recovery_unknown=0
    local n_recovery_skipped=0
    local n_recovery_equality=0
    local n_recovery_auxiliary=0
    local n_exact_unbounded=0
    local recovery_integer_ops=0
    local recovery_steps=0
    local certificate_integer_ops=0
    if "`engine'"=="simplex" {
        local n_cert_primary_attempted=r(n_cert_primary_attempted)
        local n_cert_primary_verified=r(n_cert_primary_verified)
        local n_cert_secondary_attempted=r(n_cert_secondary_attempted)
        local n_cert_secondary_verified=r(n_cert_secondary_verified)
        local n_cert_repair_eq_stops=r(n_cert_repair_eq_stops)
        local n_cert_repair_both_stops=r(n_cert_repair_both_stops)
        local n_cert_mapping_failed=r(n_cert_mapping_failed)
        local n_recovery_requested=r(n_recovery_requested)
        local n_recovery_finite=r(n_recovery_finite)
        local n_recovery_unknown=r(n_recovery_unknown)
        local n_recovery_skipped=r(n_recovery_skipped)
        local n_recovery_equality=r(n_recovery_equality)
        local n_recovery_auxiliary=r(n_recovery_auxiliary)
        local n_exact_unbounded=r(n_exact_unbounded)
        local recovery_integer_ops=r(recovery_integer_ops)
        local recovery_steps=r(recovery_steps)
        local certificate_integer_ops=r(certificate_integer_ops)
    }

    // Diagnostics run after the original result is secured; they never write weights.
    tempname optima_faces optima_coords optima_denoms optima_counts ce_range_table ce_range_counts
    if "`optimadiag'"!="" {
        mata: `solver_snapshot'=_c57_context_save()
        capture noisily {
            mata: _dea_rec_abort(); _rsm_audit_reset(0); _rsm_success_reset(1,0); _dea_link_reset(2,0)
            mata: _c73_compute("`inputs'","`outputs'","`touse'","`rts'","`orientation'","`method'","`vrscross'",`tolerance',`_optima_supported',"`ce_v_return'","`ce_u_return'","`ce_w_return'","`ce_status_return'","`optima_faces'","`optima_coords'","`optima_denoms'","`optima_counts'","`cerange'"!="","`ce_range_table'","`ce_range_counts'")
        }
        local optima_rc=_rc
        mata: _c57_context_restore(`solver_snapshot')
        mata: mata drop `solver_snapshot'
        if `optima_rc' {
            di as error "optimadiag internal evaluation failed; existing DEA output variables were retained."
            exit `optima_rc'
        }
    }

    // Summary statistics
    qui summ `generate' if `touse'
    local mean_ce = r(mean)
    local min_ce = r(min)
    local max_ce = r(max)

    qui summ `generate'_self if `touse'
    local mean_self = r(mean)

    di as text ""
    di as text "{hline 70}"
    di as text "Summary"
    di as text "{hline 70}"
    di as text "Mean self-appraisal:      " as result %9.4f `mean_self'
    di as text "Mean cross-efficiency:    " as result %9.4f `mean_ce'
    di as text "Min cross-efficiency:     " as result %9.4f `min_ce'
    di as text "Max cross-efficiency:     " as result %9.4f `max_ce'
    di as text "{hline 70}"

    if `raw_recovered'>0 di as text "Note: " `raw_recovered' " primary solve(s) recovered with certified original-unit weights."
    if `raw_unavailable'>0 di as text "Note: " `raw_unavailable' " raw recovery attempt(s) remain unavailable; this is not an original-model infeasibility claim."
    if `raw_export_obstructed'>0 di as text "Note: " `raw_export_obstructed' " raw recovery attempt(s) encountered weights that could not be represented in Stata; this alone does not establish model infeasibility."
    if `raw_secondary_unknown'>0 di as error "Warning: " `raw_secondary_unknown' " raw secondary solve(s) remain unverified or unexportable; certified primary weights retained."
    if `raw_peers_range'>0 di as error "Warning: " `raw_peers_range' " raw peer ratio(s), including diagonal, cannot be represented without range loss."
    if `raw_peers_arithmetic'>0 di as error "Warning: " `raw_peers_arithmetic' " raw peer ratio(s), including diagonal, exceeded the arithmetic budget."
    if `_failed_primary'>0 di as error "Warning: " `_failed_primary' " primary weight solve(s) unavailable."
    if `_failed_secondary'>0 di as error "Warning: " `_failed_secondary' " secondary solve(s) unavailable; primary weights retained."
    if `_missing_pairs'>0 di as error "Warning: " `_missing_pairs' " peer ratio(s) unavailable; averages use available peers."
    if "`rts'"=="vrs" & "`vrscross'"=="legacy" di as text "Note: VRS peer ratios omit the intercept; this is not the Lim-Zhu correction."
    if "`vrscross'"=="limzhu" {
        di as text "Note: Lim-Zhu input peer formula; diagonal is translated appraisal, not the self score."
        di as text "      Peer means exclude the diagonal; undefined or unstable ratios remain missing."
        if "`method'"!="arbitrary" di as text "      The existing secondary surrogate is retained; it does not optimize corrected peer ratios."
    }

    if `norm_lost'>0 {
        di as error "Warning: " `norm_lost' " positive quantity value(s) became zero during global column scaling in " `norm_obs_count' " selected observation(s)."
        di as text "      This is a numerical representation issue; it does not establish infeasibility of the original data."
        if "`detail'"!="" mata: _c60_loss_detail("`inputs'","`outputs'","`touse'","`dmu'")
    }
    if `_failed_primary'>0 di as text "Note: Missing numerical primary weights do not establish infeasibility of the original-data model."
    if "`engine'"=="simplex" {
        if `n_cert_mapping_failed'>0 di as error "Warning: Original-data certificate mapping unavailable for " `n_cert_mapping_failed' " observation(s); see certification_map_reason with matrix."
        if (`n_cert_primary_attempted'-`n_cert_primary_verified'+`n_cert_secondary_attempted'-`n_cert_secondary_verified')>0 di as error "Warning: Some original-data certificates are unverified; matrix returns stage-specific reasons and repair stops."
    }

    if `n_exact_unbounded'>0 di as text "Note: " as result `n_exact_unbounded' as text " original-data LP(s) have an exact feasible-point and improving-ray certificate; primary fallback retained."

    if `native_bad_primary'>0 di as error "Warning: " `native_bad_primary' " primary LP(s) skipped because prepared coefficients exceeded the finite numeric range."
    if `native_bad_secondary'>0 di as error "Warning: " `native_bad_secondary' " secondary LP(s) skipped because prepared coefficients exceeded the finite numeric range; primary weights retained."
    if `_failed_ce'>0 di as error "Warning: " `_failed_ce' " observation(s) have no available peer mean."
    if `objective_build_failed'>0 di as error "Warning: " `objective_build_failed' " peer objective build(s) unavailable; primary weights retained."
    if `objective_warning'>0 & `objective_build_failed'==0 di as error "Warning: Some numerical peer objectives were rounded or lost positive terms; see objective diagnostics."
    local numerical_warning=(`raw_peers_unavailable'>0 | (`objective_warning'>0 | `_failed_primary'>0 | `_failed_secondary'>0 | `_missing_pairs'>0 | `_failed_ce'>0 | `norm_lost'>0 | `n_cert_mapping_failed'>0 | (`n_cert_primary_attempted'-`n_cert_primary_verified'+`n_cert_secondary_attempted'-`n_cert_secondary_verified')>0))
    if `numerical_warning' di as text "      Use detail and matrix for affected observations and stage diagnostics."

    if "`optimadiag'"!="" {
        if `optima_counts'[1,4]>0 & "`cerange'"=="" di as text "Note: Multiple optimal weights verified on " `optima_counts'[1,4] " original-data face(s); CE/rank variability is not determined by this diagnostic."
        if `optima_counts'[1,5]>0 di as error "Warning: Uniqueness remains unverified on " `optima_counts'[1,5] " face(s); no unique-weight claim is made."
        if `optima_counts'[1,6]>0 di as error "Warning: optimadiag is outside its supported dimensions or peer rule; DEA results retained."
        if `optima_counts'[1,7]>0 di as error "Warning: Some diagnostic objectives are exactly unbounded; no finite optimal-face claim is made."
        if `optima_counts'[1,12]>0 di as error "Warning: Some peer denominators are identically zero on an exact optimal face; those ratios are undefined."
        if `optima_counts'[1,13]>0 | `optima_counts'[1,15]>0 di as error "Warning: Peer denominator availability can depend on the optimal weights."
        if `optima_counts'[1,14]>0 di as error "Warning: Some denominator states remain unverified within the diagnostic budget."
        if `optima_counts'[1,9]>0 di as error "Warning: Some stored numerical weights are off the exact diagnostic face; existing DEA results retained."
        if `optima_counts'[1,18]>0 di as error "Warning: Some stored-weight links to diagnostic faces remain unverified."
    }

    if "`cerange'"!="" {
        local cerange_warning=0
        if `ce_range_counts'[1,4]>0 {
            local cerange_warning=1
            di as error "Warning: Some peer CE ratios vary across exact optimal weights; stored DEA results are unchanged."
        }
        if `ce_range_counts'[1,5]>0 {
            local cerange_warning=1
            di as error "Warning: Some exact-face CE ratios have an identically zero denominator and are undefined."
        }
        if `ce_range_counts'[1,6]>0 {
            local cerange_warning=1
            di as error "Warning: Some CE ranges remain unverified within the diagnostic limits."
        }
        if `ce_range_counts'[1,7]>0 {
            local cerange_warning=1
            di as error "Warning: Some CE ranges are unsupported by the current dimension, peer-rule, mapping or denominator restrictions; DEA results retained."
        }
        if `ce_range_counts'[1,8]>0 {
            local cerange_warning=1
            di as error "Warning: Some CE ranges have no verified finite optimal face; no substitute face is used."
        }
        if `ce_range_counts'[1,9]>0 {
            local cerange_warning=1
            di as error "Warning: Some CE bounds are limits and are not attained by finite positive-denominator weights."
        }
        if `ce_range_counts'[1,10]>0 {
            local cerange_warning=1
            di as error "Warning: Some CE ranges are exactly unbounded."
        }
        if `ce_range_counts'[1,11]>0 {
            local cerange_warning=1
            di as error "Warning: Some exact-face CE denominators can be zero; ranges use only positive denominators."
        }
        if `ce_range_counts'[1,12]>0 {
            local cerange_warning=1
            di as error "Warning: The numerical peer cutoff can include or exclude a cell depending on optimal weights."
        }
        if `ce_range_counts'[1,13]>0 {
            local cerange_warning=1
            di as error "Warning: Some exact CE bounds cannot be represented faithfully as Stata numeric values."
        }
        if `ce_range_counts'[1,14]>0 {
            local cerange_warning=1
            di as error "Warning: Some CE bound displays remain unverified."
        }
        if `ce_range_counts'[1,15]>0 {
            local cerange_warning=1
            di as error "Warning: Some stored numerical weights are off the exact CE diagnostic face."
        }
        if `ce_range_counts'[1,16]>0 {
            local cerange_warning=1
            di as error "Warning: Some exact-face CE ranges have no stored numerical weights to link."
        }
        if `ce_range_counts'[1,17]>0 {
            local cerange_warning=1
            di as error "Warning: Some stored-weight links to CE diagnostic faces remain unverified."
        }
        if `ce_range_counts'[1,18]>0 {
            local cerange_warning=1
            di as error "Warning: Some reconstructed CE endpoint links remain unverified."
        }
    }

    // Return values
    if "`cerange'"!="" {
        return scalar cerange_version=7701
        return local cerange_scope "bounded_original_exact_faces"
        return local cerange_engine "independent_simplex_exact"
        return scalar cerange_mean_ranges_available=0
        return scalar cerange_rank_ranges_available=0
        return scalar has_cerange_warning=`cerange_warning'
        return scalar n_cerange_requested=`ce_range_counts'[1,1]
        return scalar n_cerange_verified=`ce_range_counts'[1,2]
        return scalar n_cerange_fixed=`ce_range_counts'[1,3]
        return scalar n_cerange_variable=`ce_range_counts'[1,4]
        return scalar n_cerange_undefined=`ce_range_counts'[1,5]
        return scalar n_cerange_unknown=`ce_range_counts'[1,6]
        return scalar n_cerange_unsupported=`ce_range_counts'[1,7]
        return scalar n_cerange_face_unavailable=`ce_range_counts'[1,8]
        return scalar n_cerange_limit_only=`ce_range_counts'[1,9]
        return scalar n_cerange_unbounded=`ce_range_counts'[1,10]
        return scalar n_cerange_den_zero_possible=`ce_range_counts'[1,11]
        return scalar n_cerange_cutoff_variable=`ce_range_counts'[1,12]
        return scalar n_cerange_display_loss=`ce_range_counts'[1,13]
        return scalar n_cerange_display_unknown=`ce_range_counts'[1,14]
        return scalar n_cerange_offface=`ce_range_counts'[1,15]
        return scalar n_cerange_no_weights=`ce_range_counts'[1,16]
        return scalar n_cerange_link_unknown=`ce_range_counts'[1,17]
        return scalar n_cerange_raw_link_unknown=`ce_range_counts'[1,18]
        return scalar n_cerange_lp_calls=`ce_range_counts'[1,19]
        return scalar n_cerange_integer_ops=`ce_range_counts'[1,20]
        if "`matrix'"!="" {
            matrix colnames `ce_range_table'=obs face peer_obs face_state range_state denominator_state cutoff_state lower upper lower_state upper_state lower_display upper_display raw_lower_link raw_upper_link stored_weight_link public_stage reason row_lp_calls row_integer_ops
            return matrix ce_ranges=`ce_range_table'
        }
    }
    if "`optimadiag'"!="" {
        return scalar optima_diagnostic_version=7301
        return local optima_scope "bounded_original_exact_faces"
        return local optima_engine "independent_simplex_exact"
        return local optima_link_scope "stored_weights_only"
        return scalar optima_ce_ranges_available=("`cerange'"!="")
        return scalar optima_rank_ranges_available=0
        return scalar n_optima_requested=`optima_counts'[1,1]
        return scalar n_optima_certified=`optima_counts'[1,2]
        return scalar n_optima_unique=`optima_counts'[1,3]
        return scalar n_optima_multiple=`optima_counts'[1,4]
        return scalar n_optima_unknown=`optima_counts'[1,5]
        return scalar n_optima_unsupported=`optima_counts'[1,6]
        return scalar n_optima_unbounded=`optima_counts'[1,7]
        return scalar n_optima_linked=`optima_counts'[1,8]
        return scalar n_optima_offface=`optima_counts'[1,9]
        return scalar n_optima_no_weights=`optima_counts'[1,10]
        return scalar n_optima_den_positive=`optima_counts'[1,11]
        return scalar n_optima_den_zero=`optima_counts'[1,12]
        return scalar n_optima_den_possible_zero=`optima_counts'[1,13]
        return scalar n_optima_den_unknown=`optima_counts'[1,14]
        return scalar n_optima_cutoff_variable=`optima_counts'[1,15]
        return scalar n_optima_lp_calls=`optima_counts'[1,16]
        return scalar n_optima_integer_ops=`optima_counts'[1,17]
        return scalar n_optima_link_unknown=`optima_counts'[1,18]

        if "`matrix'"!="" {
            matrix colnames `optima_faces'=obs face face_state weight_state coords_done den_done weight_link public_stage reason row_lp_calls row_integer_ops
            matrix colnames `optima_coords'=obs face coordinate min_state max_state range_state lower upper
            matrix colnames `optima_denoms'=obs face peer_obs min_state max_state denominator_state cutoff_state lower upper
            return matrix optima_faces=`optima_faces'
            return matrix optima_coordinates=`optima_coords'
            return matrix optima_denominators=`optima_denoms'
        }
    }
    return scalar objective_diagnostic_version=7001
    foreach _k in attempted build_failed numeric_rounded numeric_lost exact_verified warning {
        return scalar n_objective_`_k'= `objective_`_k''
    }
    return scalar numerical_diagnostic_version=6601
    return scalar has_numerical_warning=`numerical_warning'
    return scalar n_native_nonfinite_primary=`native_bad_primary'
    return scalar n_native_nonfinite_secondary=`native_bad_secondary'
    if "`engine'"=="ipm" return local native_numerical_scope "prepared_array_preflight"
    else return local native_numerical_scope "not_evaluated_simplex"
    if "`engine'"=="simplex" {
        return scalar raw_recovery_version=8401
        return local raw_recovery_scope "failed_primary_only_certified_raw_weights"
        return local raw_peer_scope "recovered_rows_only_exported_weights"
        foreach _k in attempted recovered unavailable primary_unknown export_obstructed secondary_unknown secondary_unbounded peers_evaluated peers_unavailable peers_range peers_arithmetic {
            return scalar n_raw_`_k'=`raw_`_k''
        }
        if "`matrix'"!="" {
            matrix colnames `raw_recovery_return'=obs attempted attempts representation primary_status primary_reason primary_export secondary_status secondary_reason secondary_export final_stage stop eq_error iq_error obj_error scale_ops export_ops export_obstructed
            return matrix raw_recovery_diagnostics=`raw_recovery_return'
            return matrix raw_scaling_diagnostics=`raw_scaling_return'
            return matrix raw_peer_status=`raw_peer_return'
        }
    }
    return scalar diagnostic_version=6201
    return scalar n_recovery_requested=`n_recovery_requested'
    return scalar n_recovery_finite=`n_recovery_finite'
    return scalar n_recovery_unknown=`n_recovery_unknown'
    return scalar n_recovery_skipped=`n_recovery_skipped'
    return scalar n_recovery_equality=`n_recovery_equality'
    return scalar n_recovery_auxiliary=`n_recovery_auxiliary'
    return scalar n_exact_unbounded=`n_exact_unbounded'
    return scalar recovery_integer_ops=`recovery_integer_ops'
    return scalar recovery_steps=`recovery_steps'
    return scalar certificate_integer_ops=`certificate_integer_ops'
    if "`engine'"=="simplex" return local recovery_scope "bounded_exact_failure_paths"
    else return local recovery_scope "not_evaluated_native"
    return scalar n_normalization_lost_positive=`norm_lost'
    return scalar n_normalization_lost_inputs=`norm_in'
    return scalar n_normalization_lost_outputs=`norm_out'
    return scalar n_normalization_lost_obs=`norm_obs_count'
    return scalar n_numeric_primary_unavailable=`_failed_primary'
    return local normalization_diagnostic_scope "global_column_max"
    if "`engine'"=="simplex" return local certificate_diagnostic_scope "bounded_original_data_optimal_face"
    else return local certificate_diagnostic_scope "not_evaluated_native"
    return scalar n_cert_primary_attempted=`n_cert_primary_attempted'
    return scalar n_cert_primary_verified=`n_cert_primary_verified'
    return scalar n_cert_secondary_attempted=`n_cert_secondary_attempted'
    return scalar n_cert_secondary_verified=`n_cert_secondary_verified'
    return scalar n_cert_repair_eq_stops=`n_cert_repair_eq_stops'
    return scalar n_cert_repair_both_stops=`n_cert_repair_both_stops'
    return scalar n_cert_mapping_failed=`n_cert_mapping_failed'
    return local engine "`engine'"
    return local engine_requested "`engine'"
    return local actual_engine "`engine'"
    return local method "`method'"
    return local rts "`rts'"
    return local orientation "`orientation'"
    return scalar n = `n'
    return scalar n_failed_primary=`_failed_primary'
    return scalar n_failed_secondary=`_failed_secondary'
    return scalar n_failed_ce=`_failed_ce'
    return scalar n_missing_peer_pairs=`_missing_pairs'
    return scalar n_missing_diagonal=`_missing_diagonal'
    return scalar tolerance=`tolerance'
    return local vrscross "`vrscross'"
    if "`vrscross'"=="limzhu" return local peer_rule "limzhu_input"
    else return local peer_rule "ratio_no_intercept"
    return scalar mean_ce = `mean_ce'
    return scalar mean_self = `mean_self'
    if "`matrix'" != "" {
        matrix colnames `ce_v_return'=`inputs'
        matrix colnames `ce_u_return'=`outputs'
        matrix colnames `ce_w_return'=intercept
        matrix colnames `ce_status_return'=status
        matrix colnames `ce_peers_return'=n_peers
        if "`engine'"=="simplex" {
            matrix colnames `ce_diag_return'=p_status p_solver p_obj p_iters p_scaled p_absolute s_status s_solver s_obj s_iters s_scaled s_absolute
            matrix colnames `ce_link_return'=p_cert s_cert applied reason p_ops s_ops p_solver s_solver p_reason face_rows
            return matrix link_diagnostics=`ce_link_return'
            return local link_rule "original_data_certified_face"
            return matrix solver_diagnostics=`ce_diag_return'
        }
        mata: _c60_computation("`ce_status_return'","`norm_inputs'","`norm_outputs'","`norm_obs'","`computation_diag'")
        if "`engine'"!="simplex" {
            mata: st_matrix("`cert_p_diag'",_c60_empty_cert(`n'));st_matrix("`cert_s_diag'",_c60_empty_cert(`n'));st_matrix("`cert_map_diag'",J(`n',1,.))
        }
        if "`engine'"!="simplex" {
            mata: st_matrix("`recovery_p_diag'",_c62_empty_recovery(`n'));st_matrix("`recovery_s_diag'",_c62_empty_recovery(`n'))
        }
        matrix colnames `recovery_p_diag'=requested outcome path stop initial_status initial_reason original_steps eq_attempted eq_steps aux_attempted aux_steps final_steps total_steps ray_verified initial_ops ray_ops eq_ops aux_ops final_ops total_ops maxlimbs returned_stop eq_status eq_reason eq_stop aux_status aux_reason aux_stop final_status final_reason final_stop
        matrix colnames `recovery_s_diag'=requested outcome path stop initial_status initial_reason original_steps eq_attempted eq_steps aux_attempted aux_steps final_steps total_steps ray_verified initial_ops ray_ops eq_ops aux_ops final_ops total_ops maxlimbs returned_stop eq_status eq_reason eq_stop aux_status aux_reason aux_stop final_status final_reason final_stop
        return matrix primary_recovery_diagnostics=`recovery_p_diag'
        return matrix secondary_recovery_diagnostics=`recovery_s_diag'
        if "`engine'"!="ipm" mata: st_matrix("`native_numeric_return'",J(`n',4,.))
        mata: st_matrix("`native_numeric_return'",(st_matrix("`norm_obs'"),st_matrix("`native_numeric_return'")))
        matrix colnames `native_numeric_return'=obs primary_state primary_rc secondary_state secondary_rc
        return matrix native_numerical_diagnostics=`native_numeric_return'
        mata: st_matrix("`objective_return'",(st_matrix("`norm_obs'"),st_matrix("`objective_return'")))
        matrix colnames `objective_return'=obs attempted build_state numeric_state scope build_ops max_limbs
        return matrix secondary_objective_diagnostics=`objective_return'

        matrix colnames `computation_diag'=obs primary_available weight_stage lost_positive_cells
        matrix colnames `norm_inputs'=`inputs'
        matrix colnames `norm_outputs'=`outputs'
        matrix colnames `cert_p_diag'=attempted status reason repair_stop repair_steps maxlimbs eq_violation iq_violation dual_violation
        matrix colnames `cert_s_diag'=attempted status reason repair_stop repair_steps maxlimbs eq_violation iq_violation dual_violation
        matrix colnames `cert_map_diag'=map_reason
        return matrix computation_diagnostics=`computation_diag'
        return matrix normalization_input_loss=`norm_inputs'
        return matrix normalization_output_loss=`norm_outputs'
        return matrix primary_cert_diagnostics=`cert_p_diag'
        return matrix secondary_cert_diagnostics=`cert_s_diag'
        return matrix certification_map_reason=`cert_map_diag'
        return matrix ce_matrix = `ce_matrix_return'
        return matrix weights_input=`ce_v_return'
        return matrix weights_output=`ce_u_return'
        return matrix intercept=`ce_w_return'
        return matrix stage_status=`ce_status_return'
        return matrix peer_count=`ce_peers_return'
    }
end


// ============================================================================
// MATA FUNCTIONS
// ============================================================================

capture mata: mata drop _dea_dmu_read()
capture mata: mata drop _dea_crosseff_compute()
capture mata: mata drop _dea_crosseff_multiplier()
capture mata: mata drop _dea_crosseff_lzratio()
capture mata: mata drop _dea_rank()

mata:
mata set matastrict on

// ─────────────────────────────────────────────────────────────────────────────
// _dea_dmu_read()
//
// Reads a DMU variable that may be string or numeric.
// Returns: dmu_ids (numeric, 1-indexed for strings; original values for numeric)
//          dmu_labels (string representations for display)
// ─────────────────────────────────────────────────────────────────────────────
void _dea_dmu_read(string scalar  dmu_var,
                   string scalar  touse,
                   real colvector dmu_ids,
                   string colvector dmu_labels)
{
    real scalar n, i, j
    string colvector dmu_str, uniq

    if (substr(st_vartype(dmu_var), 1, 3) == "str") {
        dmu_str    = st_sdata(., dmu_var, touse)
        n          = rows(dmu_str)
        uniq       = sort(uniqrows(dmu_str), 1)
        dmu_ids    = J(n, 1, .)
        dmu_labels = dmu_str
        for (i = 1; i <= n; i++) {
            for (j = 1; j <= rows(uniq); j++) {
                if (dmu_str[i] == uniq[j]) {
                    dmu_ids[i] = j
                    break
                }
            }
        }
    }
    else {
        dmu_ids    = st_data(., dmu_var, touse)
        dmu_labels = strtrim(strofreal(dmu_ids,"%21.17g"))
    }
}


/**
 * _dea_crosseff_compute()
 *
 * Computes Cross-Efficiency scores. Supports string and numeric DMU identifiers.
 */
void _dea_crosseff_compute(string scalar inputs,
                            string scalar outputs,
                            string scalar dmu_var,
                            string scalar touse,
                            string scalar eff_ce_var,
                            string scalar eff_self_var,
                            string scalar eff_maverick_var,
                            string scalar method,
                            string scalar rts,
                            string scalar orientation,
                            string scalar vrscross,
                            real scalar tol,
                            real scalar show_matrix,
                            real scalar show_detail,
                            real scalar show_debug,
                            string scalar ce_matrix_name,
                            string scalar v_name, string scalar u_name,
                            string scalar w_name, string scalar status_name,
                            string scalar peers_name, string scalar native_numeric_name,string scalar objective_name,|real scalar capture_optima)
{
    // Load data
    real matrix X, Y
    real rowvector sx, sy
    real scalar jj
    real colvector dmu_ids, obs_indices
    string colvector dmu_labels
    string rowvector input_names, output_names

    input_names  = tokens(inputs)
    output_names = tokens(outputs)

    X           = st_data(., input_names,  touse)
    Y           = st_data(., output_names, touse)
    obs_indices = selectindex(st_data(., touse))

    // Read DMU — handles string and numeric
    _dea_dmu_read(dmu_var, touse, dmu_ids, dmu_labels)

    real scalar n, m, s
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    // A change of measurement units leaves the normalized LP unchanged.
    sx=colmax(X); sy=colmax(Y)
    for (jj=1; jj<=m; jj++) if (sx[jj]==0) sx[jj]=1
    for (jj=1; jj<=s; jj++) if (sy[jj]==0) sy[jj]=1
    X=X:/sx; Y=Y:/sy

    if (show_debug) {
        printf("  n=%g, m=%g (inputs), s=%g (outputs)\n", n, m, s)
    }

    // Step 1: Compute optimal weights for each DMU using multiplier model
    // For input-oriented CRS: max sum(u*y0) s.t. sum(v*x0)=1, sum(u*yj)-sum(v*xj)<=0

    real matrix V_opt, U_opt   // Optimal weights: V (inputs), U (outputs)
    real colvector W_opt       // VRS free intercept (0 for CRS)
    V_opt = J(n, m, .)
    U_opt = J(n, s, .)
    W_opt = J(n, 1, 0)

    real colvector self_eff, stage_status, peer_count
    self_eff = J(n, 1, .)
    stage_status=J(n,1,0); peer_count=J(n,1,0)

    real scalar i, j
    real rowvector x0, y0, v_i, u_i
    real scalar eff_i, w_i, stage_i
    real rowvector numeric_i,objective_i
    real matrix objective_diag
    objective_diag=J(n,6,0)
    real matrix numeric_diag
    numeric_diag=J(n,4,.)

    for (i = 1; i <= n; i++) {
        x0 = X[i, .]
        y0 = Y[i, .]

        // Solve multiplier model (Stage 1 + optional Stage 2 for benevolent/aggressive)
        _dea_crosseff_multiplier(x0, y0, X, Y, rts, orientation, method, tol, v_i, u_i, w_i, eff_i, stage_i, numeric_i, i, objective_i)

        V_opt[i, .] = v_i
        U_opt[i, .] = u_i
        W_opt[i]    = w_i
        self_eff[i] = eff_i
        stage_status[i]=stage_i
        numeric_diag[i,.]=numeric_i
        objective_diag[i,.]=objective_i

        if (show_debug) {
            printf("  DMU %g: self-eff = %9.4f\n", dmu_ids[i], eff_i)
        }
    }

    // Step 2: Compute cross-efficiency matrix
    // CE[i,j] = efficiency of DMU j using weights of DMU i
    real matrix CE
    CE = J(n, n, .)

    real rowvector v_i2, u_i2
    real scalar num, denom

    // Legacy: u*y/(v*x). Lim-Zhu (2015), equation (4), input VRS:
    // u*y/(v*x-w), where our input multiplier w equals -xi in the paper.
    // Apply the chosen formula to the diagonal too; do not replace it by self.
    for (i = 1; i <= n; i++) {
        v_i2 = V_opt[i, .]
        u_i2 = U_opt[i, .]

        for (j = 1; j <= n; j++) {
            num   = u_i2 * Y[j, .]'
            denom = v_i2 * X[j, .]'
            // CRS permits genuinely small positive peer denominators. Under
            // the legacy VRS ratio convention retain its tolerance cutoff:
            // the intercept can carry the entire output-oriented self value.
            if (vrscross=="limzhu") CE[i,j]=_dea_crosseff_lzratio(denom,num,W_opt[i],tol)
            else CE[i, j] = (denom>(rts=="vrs" ? tol : 0) & denom<. & num>=0 & num<. ? num/denom : .)
        }
    }

    // Step 3: Compute cross-efficiency scores (column means excluding self)
    real colvector ce_scores
    ce_scores = J(n, 1, .)

    real scalar sum_ce, count_ce

    for (j = 1; j <= n; j++) {
        sum_ce = 0
        count_ce = 0

        for (i = 1; i <= n; i++) {
            if (i != j && CE[i, j] < .) {
                sum_ce = sum_ce + CE[i, j]
                count_ce++
            }
        }

        peer_count[j]=count_ce
        if (count_ce > 0) {
            ce_scores[j] = sum_ce / count_ce
        }
    }

    // Step 4: Compute maverick index (self-appraisal vs cross-efficiency)
    real colvector maverick
    maverick = J(n,1,.)
    for (i=1; i<=n; i++) {
        if (self_eff[i]<. & ce_scores[i]>0 & ce_scores[i]<.) maverick[i]=self_eff[i]/ce_scores[i]
    }

    // Compute max label width for aligned display
    real scalar lw
    string scalar fmt_dmu
    lw = 10
    for (i = 1; i <= n; i++) {
        if (strlen(dmu_labels[i]) > lw) lw = strlen(dmu_labels[i])
    }
    fmt_dmu = "%" + strofreal(lw) + "s"

    // Display results — use string labels
    printf("\n")
    printf("Cross-Efficiency Results\n")
    printf("{hline 70}\n")
    printf(fmt_dmu + " | %10s %10s %10s %10s\n",
           "DMU", "Self", "Cross-Eff", "Maverick", "Rank")
    printf(fmt_dmu + " | %10s %10s %10s %10s\n",
           "-":*lw, "----------", "----------", "----------", "----------")

    // Rank by cross-efficiency
    real colvector ranks
    ranks = _dea_rank(ce_scores)

    for (i = 1; i <= n; i++) {
        printf(fmt_dmu + " | %10.4f %10.4f %10.4f %10.0f\n",
               dmu_labels[i], self_eff[i], ce_scores[i], maverick[i], ranks[i])
    }
    printf("{hline 70}\n")

    // Show cross-efficiency matrix if requested
    if (show_matrix) {
        printf("\nCross-Efficiency Matrix (row i evaluates column j)\n")
        printf("{hline 70}\n")

        // Header
        printf(fmt_dmu + " |", "Evaluator")
        for (j = 1; j <= n; j++) {
            printf(" %7s", dmu_labels[j])
        }
        printf("\n")

        printf(fmt_dmu + " |", "-":*lw)
        for (j = 1; j <= n; j++) {
            printf(" -------")
        }
        printf("\n")

        // Data
        for (i = 1; i <= n; i++) {
            printf(fmt_dmu + " |", dmu_labels[i])
            for (j = 1; j <= n; j++) {
                printf(" %7.4f", CE[i, j])
            }
            printf("\n")
        }
        printf("{hline 70}\n")
        st_matrix(objective_name,objective_diag)
        st_matrix(ce_matrix_name, CE)
        st_matrix(v_name,V_opt:/sx)
        st_matrix(u_name,U_opt:/sy)
        st_matrix(w_name,W_opt)
        st_matrix(status_name,stage_status)
        st_matrix(native_numeric_name,numeric_diag)
        st_matrix(peers_name,peer_count)
    }

    if (!show_matrix & capture_optima==1) {
        st_matrix(v_name,V_opt:/sx);st_matrix(u_name,U_opt:/sy)
        st_matrix(w_name,W_opt);st_matrix(status_name,stage_status)
    }

    if (show_detail) {
        printf("\nWeight stages: 0=primary unavailable, 1=primary only, 2=secondary accepted, 3=primary fallback\n")
        for (i=1; i<=n; i++) printf(fmt_dmu+" | stage %g | available peers %g\n",dmu_labels[i],stage_status[i],peer_count[i])
    }
    _c70_counts(objective_diag)
    if(show_detail) _c70_detail(objective_diag,obs_indices,dmu_labels)
    if(show_detail) _c66_native_detail(numeric_diag,obs_indices,dmu_labels)
    st_numscalar("r(n_native_nonfinite_primary)",sum(numeric_diag[.,1]:==2))
    st_numscalar("r(n_native_nonfinite_secondary)",sum(numeric_diag[.,3]:==2))
    st_numscalar("r(n_failed_primary)",sum(stage_status:==0))
    st_numscalar("r(n_failed_secondary)",sum(stage_status:==3))
    st_numscalar("r(n_failed_ce)",sum(ce_scores:>=.))
    st_numscalar("r(n_missing_peer_pairs)",n*(n-1)-sum(peer_count))
    st_numscalar("r(n_missing_diagonal)",sum(diagonal(CE):>=.))

    // Store results in Stata
    for (i = 1; i <= n; i++) {
        st_store(obs_indices[i], st_varindex(eff_ce_var), ce_scores[i])
        st_store(obs_indices[i], st_varindex(eff_self_var), self_eff[i])
        st_store(obs_indices[i], st_varindex(eff_maverick_var), maverick[i])
    }
}

// Safely evaluate the translated input VRS ratio with zero weight bounds.
// Positive data do not exclude zero weighted output or a 0/0 denominator.
// A relative cancellation check is homogeneous in (v,u,w). No clipping,
// epsilon imputation, or artificial positive weight floor is applied.
real scalar _dea_crosseff_lzratio(real scalar vx, real scalar uy,
                                 real scalar w, real scalar tol)
{
    real scalar den, scale, score
    if (vx>=. | uy>=. | w>=. | vx<0 | uy<0) return(.)
    den=vx-w
    scale=max((abs(vx),abs(w),uy))
    if (den>=. | den<=0 | den<=tol*scale) return(.)
    score=uy/den
    if (score>=. | score<0 | score>1+tol) return(.)
    return(score)
}


/**
 * _dea_crosseff_multiplier()
 *
 * Solves the multiplier (dual) DEA model to get optimal weights.
 *
 * For input-oriented CRS:
 *   max  sum(u * y0)
 *   s.t. sum(v * x0) = 1
 *        sum(u * yj) - sum(v * xj) <= 0  for all j
 *        u, v >= 0
 */
void _dea_crosseff_multiplier(real rowvector x_original,
                               real rowvector y_original,
                               real matrix X_original,
                               real matrix Y_original,
                               string scalar rts,
                               string scalar orientation,
                               string scalar method,
                               real scalar tol,
                               real rowvector v_opt,
                               real rowvector u_opt,
                               real scalar w_opt,
                               real scalar eff, real scalar stage_status, real rowvector numeric_diag,real scalar index,real rowvector objective_diag)
{
    real scalar n, m, s, j
    real matrix X, Y
    real rowvector x0, y0, sx, sy
    real colvector row_scale
    // Reparameterize weights around the evaluated row, then balance each
    // peer constraint. This preserves the LP and keeps normalization usable
    // when DMUs have very different magnitudes within the same column.
    sx=x_original; sy=y_original
    for (j=1; j<=cols(sx); j++) if (sx[j]<=0) sx[j]=1
    for (j=1; j<=cols(sy); j++) if (sy[j]<=0) sy[j]=1
    X=X_original:/sx; Y=Y_original:/sy
    x0=x_original:/sx; y0=y_original:/sy
    n = rows(X)
    m = cols(X)
    s = cols(Y)

    real scalar is_vrs, is_output, nvars, rc
    real rowvector c, lb, ub
    real matrix    A_ineq, A_eq
    real colvector b_ineq, b_eq
    class LinearProgram scalar S
    real rowvector sol

    stage_status=0
    numeric_diag=(0,.,0,.)
    objective_diag=J(1,6,0)
    is_vrs    = (rts == "vrs")
    is_output = (orientation == "output")
    nvars     = m + s + is_vrs        // add 1 free variable w if VRS

    // ── Objective (Stage 1: self-efficiency) ─────────────────────────
    if (!is_output) {
        // Input: min −u·y₀ [− w]
        c = (J(1, m, 0), -y0)
        if (is_vrs) c = (c, -1)
    }
    else {
        // Output: min v·x₀ [+ w]
        c = (x0, J(1, s, 0))
        if (is_vrs) c = (c, 1)
    }

    // ── Equality constraint ──────────────────────────────────────────
    if (!is_output) {
        // Input: v·x₀ = 1
        A_eq = !is_vrs ? (x0, J(1, s, 0)) : (x0, J(1, s, 0), 0)
    }
    else {
        // Output: u·y₀ = 1
        A_eq = !is_vrs ? (J(1, m, 0), y0) : (J(1, m, 0), y0, 0)
    }
    b_eq = 1

    // ── Inequality constraints: [−X, Y {, ±1}]·z ≤ 0 ────────────────
    if (!is_vrs) {
        A_ineq = (-X, Y)                          // CRS: −v·xⱼ + u·yⱼ ≤ 0
    }
    else if (!is_output) {
        A_ineq = (-X, Y, J(n, 1,  1))            // Input VRS:  +w term (u·y+w ≤ v·x)
    }
    else {
        A_ineq = (-X, Y, J(n, 1, -1))            // Output VRS: −w term (u·y−w ≤ v·x)
    }
    row_scale=rowmax(abs(A_ineq))
    for (j=1; j<=n; j++) if (row_scale[j]==0) row_scale[j]=1
    A_ineq=A_ineq:/row_scale
    b_ineq = J(n, 1, 0)

    // ── Bounds: nonnegative v, u; w free ────────────────────────────────────
    lb = !is_vrs ? J(1, nvars, 0) : (J(1, m+s, 0), .)
    ub = J(1, nvars, .)

    // ── Stage 1: Solve self-efficiency ───────────────────────────────
    // Missing bounds below represent unbounded variables by design. All
    // coefficient/RHS arrays must be finite before any LinearProgram setter.
    if(hasmissing(c) | hasmissing(A_eq) | hasmissing(b_eq) | hasmissing(A_ineq) | hasmissing(b_ineq)) {
        numeric_diag[1]=2
        v_opt=J(1,m,.);u_opt=J(1,s,.);w_opt=.;eff=.
        return
    }
    S = LinearProgram()
    S.setMaxOrMin("min")
    S.setCoefficients(c)
    S.setEquality(A_eq, b_eq)
    S.setInequality(A_ineq, b_ineq)
    S.setBounds(lb, ub)
    S.setTrace("off")

    rc=_dea_nonte_optimize(S,tol,1e-12)
    numeric_diag[2]=rc;numeric_diag[1]=3

    if (rc==0) {
        sol   = S.parameters()
        v_opt = sol[1..m]:/sx
        u_opt = sol[(m+1)..(m+s)]:/sy
        w_opt = (is_vrs ? sol[nvars] : 0)
        if (!is_output) {
            eff = -S.value()           // Input:  −(−u·y₀ [−w]) = u·y₀ [+w]
        }
        else {
            eff = (S.value()>0 & S.value()<. ? 1/S.value() : .)   // Output: 1/(v·x₀ [+w])
        }
    }
    else {
        v_opt = J(1, m, .)
        u_opt = J(1, s, .)
        w_opt = .
        eff   = .
        return
    }

    if (eff<=0 | eff>=.) {
        numeric_diag[1]=4
        v_opt=J(1,m,.); u_opt=J(1,s,.); w_opt=.; eff=.
        return
    }
    stage_status=1
    numeric_diag[1]=1

    // ── Stage 2: Benevolent / Aggressive secondary LP ────────────────
    // Preserve the existing surrogate objectives (not a sum of peer ratios).
    // Input uses aggregate peer slack without the intercept; output uses peer inputs.
    if (method != "arbitrary" && eff < .) {
        stage_status=3
        real scalar eff_prime
        real rowvector c2, A_eq2_row,safe_c2
        real matrix    A_eq2
        real colvector b_eq2
        class LinearProgram scalar S2
        real rowvector sol2

        eff_prime = eff

        // Second equality: lock self-efficiency at eff_prime
        if (!is_output) {
            // u·y₀ [+w] = eff_prime
            A_eq2_row = !is_vrs ? (J(1,m,0), y0) : (J(1,m,0), y0, 1)
            b_eq2 = (1 \ eff_prime)
        }
        else {
            // v·x₀ [+w] = 1/eff_prime
            A_eq2_row = !is_vrs ? (x0, J(1,s,0)) : (x0, J(1,s,0), 1)
            b_eq2 = (1 \ (1/eff_prime))
        }
        A_eq2 = (A_eq \ A_eq2_row)

        // Secondary objective
        real rowvector sum_y, sum_x
        sum_y = colsum(Y) - y0    // Σ_{j≠k} y_j  (1×s)
        sum_x = colsum(X) - x0    // Σ_{j≠k} x_j  (1×m)

        if (!is_output) {
            // The input secondary goal uses aggregate peer constraint slack
            // d_j = v*x_j - u*y_j.  Benevolent minimizes sum d_j; aggressive
            // maximizes it, while preserving the Stage-1 self-efficiency.
            c2 = !is_vrs ? (sum_x, -sum_y) : (sum_x, -sum_y, 0)
            if (method == "aggressive") c2 = -c2
        }
        else {
            // Benevolent: min v·sum_x  (minimise peer input → higher CE for peers)
            // Aggressive: max v·sum_x  → min −v·sum_x
            c2 = !is_vrs ? (sum_x, J(1,s,0)) : (sum_x, J(1,s,0), 0)
            if (method == "aggressive") c2 = -c2
        }

        _c70_numeric(X,Y,index,is_output,is_vrs,method=="aggressive",c2,safe_c2,objective_diag)
        if(objective_diag[2]!=1 | objective_diag[3]==3) return
        c2=safe_c2
        if(hasmissing(c2) | hasmissing(A_eq2) | hasmissing(b_eq2) | hasmissing(A_ineq) | hasmissing(b_ineq)) {
            numeric_diag[3]=2
            return
        }
        S2 = LinearProgram()
        S2.setMaxOrMin("min")
        S2.setCoefficients(c2)
        S2.setEquality(A_eq2, b_eq2)
        S2.setInequality(A_ineq, b_ineq)
        S2.setBounds(lb, ub)
        S2.setTrace("off")

        rc=_dea_nonte_optimize(S2,tol,1e-12)
        numeric_diag[4]=rc;numeric_diag[3]=(rc==0 ? 1 : 3)

        if (rc==0) {
            sol2  = S2.parameters()
            v_opt = sol2[1..m]:/sx
            u_opt = sol2[(m+1)..(m+s)]:/sy
            w_opt = (is_vrs ? sol2[nvars] : 0)
            stage_status=2
            // eff unchanged — keeps Stage 1 self-efficiency
        }
        // If Stage 2 fails: keep Stage 1 solution (already in v_opt, u_opt, w_opt)
    }
}


/**
 * _dea_rank()
 *
 * Returns ranks (1 = highest value)
 */
real colvector _dea_rank(real colvector x)
{
    real scalar n, i, j, rank
    real colvector ranks

    n = rows(x)
    ranks = J(n, 1, .)

    for (i = 1; i <= n; i++) {
        if (x[i] == .) {
            ranks[i] = .
            continue
        }

        rank = 1
        for (j = 1; j <= n; j++) {
            if (j != i && x[j] != . && x[j] > x[i]) {
                rank++
            }
        }
        ranks[i] = rank
    }

    return(ranks)
}

end
