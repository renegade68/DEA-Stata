*! version 2.4.0  08nov2026
*! Validated scale-free LPs; complete Nash iterations; matrix/detail and diagnostics.
*! 2026-11-01: honor if/in for evaluation and reference samples; reject empty samples.
*! 2026-11-01: preflight output names; explicit replacevars; protect model data.
*! dea_crossgame.ado — Cross-Efficiency Game DEA
*! Includes: Aggressive, Benevolent, Nash Equilibrium strategies
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
*!   Liang, L., Wu, J., Cook, W. D., & Zhu, J. (2008). The DEA game
*!     cross-efficiency model and its Nash equilibrium. Operations
*!     Research, 56(5), 1278-1288. doi:10.1287/opre.1070.0487
*!
*!   Wu, J., Liang, L., & Yang, F. (2009). Achievement and benchmarking
*!     of countries at the Summer Olympics using cross efficiency
*!     evaluation method. European Journal of Operational Research,
*!     197(2), 722-730. doi:10.1016/j.ejor.2008.06.030
*!
*! Changelog:
*!   2.2.0  03may2026  Implement Liang et al. (2008) model (5) for strategy(nash)
*!   2.1.0  28apr2026  Add Form A syntax (x1 x2 = y1 y2)

/**
 * CROSS-EFFICIENCY GAME DEA
 * =========================
 *
 * Models DMU competition/cooperation through cross-efficiency evaluation.
 *
 * Strategies:
 *   - Aggressive: Minimize the weighted-output secondary surrogate
 *   - Benevolent: Maximize the weighted-output secondary surrogate
 *   - Nash: Iterate model (5), with explicit numerical stopping status
 *
 * References:
 *   - Doyle & Green (1994): Cross-efficiency derivations, meanings, and uses
 *   - Liang et al. (2008): DEA game cross-efficiency model and Nash equilibrium
 *   - Wu et al. (2009): Cross-efficiency benchmarking application
 */

capture program drop dea_crossgame
program define dea_crossgame, rclass
    version 16.0

    // =========================================================================
    // Unified Syntax Parser (Form A / Form B)
    //   Form A: dea_crossgame x1 x2 = y1 y2 [if] [in], [options]
    //   Form B: dea_crossgame [if] [in], inputs(x1 x2) outputs(y1 y2) [options]
    // =========================================================================
    _dea_parse_io `0'
    local 0 `"`r(args)'"'

    syntax [if] [in], ///
        Inputs(varlist numeric) ///
        Outputs(varlist numeric) ///
        [DMU(varname) ///
         STRAtegy(string) ///
         MAXiter(integer 100) ///
         GENerate(string) ///
         REPLACEVARS      ///
         TOLerance(real 1e-6) ///
         MATrix ///
         DETail ///
         DEBUG]

    if missing(`tolerance') | `tolerance'<=0 {
        di as error "tolerance() must be positive and finite"
        exit 198
    }
    if missing(`maxiter') | `maxiter'<1 {
        di as error "maxiter() must be a positive integer"
        exit 198
    }
    capture mata: assert(_dea_nonte_lp_ready()==1)
    if _rc {
        quietly findfile _dea_nonte_lp.ado
        quietly run "`r(fn)'"
    }

    // Defaults
    if "`strategy'" == "" local strategy "all"
    if "`generate'" == "" local generate "_cg"

    // Validate strategy
    if !inlist("`strategy'", "aggressive", "benevolent", "nash", "all") {
        di as error "strategy() must be aggressive, benevolent, nash, or all"
        exit 198
    }

    local _planned ""
    if inlist("`strategy'", "all", "aggressive") local _planned `_planned' `generate'_agg
    if inlist("`strategy'", "all", "benevolent") local _planned `_planned' `generate'_ben
    if inlist("`strategy'", "all", "nash") local _planned `_planned' `generate'_nash
    // Validate every output before any existing variable can be replaced.
    _dea_output_check, names(`_planned') ///
        protect(`inputs' `outputs' `dmu') `replacevars'
    local _replace_outputs "`r(existing)'"

    // Mark sample — do NOT include dmu in markout (string DMU causes r(109))
    marksample touse, novarlist
    markout `touse' `inputs' `outputs'
    if "`dmu'"!="" quietly replace `touse'=0 if missing(`dmu')

    // Count
    qui count if `touse'
    local n = r(N)
    if `n' == 0 {
        di as error "No observations in estimation sample"
        exit 2000
    }
    // Aggressive/benevolent averages exclude self; all requests both of them.
    // Check the marked sample before IDs, output replacement, or any LP solve.
    if `n' < 2 & inlist("`strategy'", "aggressive", "benevolent", "all") {
        di as error "strategy(`strategy') requires at least 2 DMUs after if/in and missing-value exclusions"
        exit 2001
    }

    foreach _side in inputs outputs {
        tempvar _positive
        quietly generate byte `_positive'=0 if `touse'
        foreach _v of local `_side' {
            quietly count if `touse' & `_v'<0
            if r(N) {
                di as error "crossgame quantities must be nonnegative"
                exit 459
            }
            quietly replace `_positive'=1 if `touse' & `_v'>0
        }
        quietly count if `touse' & `_positive'==0
        if r(N) {
            di as error "each selected observation requires a positive input and output"
            exit 459
        }
    }

    _dea_crossgame_load

    // Display header
    di as text ""
    di as text "{hline 70}"
    di as text "Cross-Efficiency Game DEA"
    di as text "{hline 70}"
    di as text "Strategy:      " as result "`strategy'"
    di as text "DMUs:          " as result "`n'"
    di as text "Max iterations:" as result "`maxiter'" as text " (for Nash)"
    di as text "{hline 70}"

    // Create DMU identifier if not provided
    if "`dmu'" == "" {
        tempvar dmu_temp
        qui gen `dmu_temp' = _n if `touse'
        local dmu "`dmu_temp'"
    }

    // Compute in tempvars. Errors or user Break preserve existing output columns.
    tempvar workagg workben worknash
    if inlist("`strategy'","all","aggressive") quietly generate double `workagg'=.
    if inlist("`strategy'","all","benevolent") quietly generate double `workben'=.
    if inlist("`strategy'","all","nash") quietly generate double `worknash'=.
    tempname cg snapshot
    mata: `snapshot'=_c57_context_save()
    capture noisily {
        mata: _dea_rec_abort(); _rsm_audit_reset(0); _rsm_success_reset(1,0); _dea_link_reset(2,0)
        mata: _dea_crossgame_compute("`inputs'", "`outputs'", ///
            "`dmu'", "`touse'", "`workagg'", "`workben'", "`worknash'", ///
            "`strategy'", `maxiter', `tolerance', ///
            "`matrix'" != "", "`detail'" != "", "`debug'" != "", "`cg'")
    }
    local compute_rc=_rc
    mata: _c57_context_restore(`snapshot')
    mata: mata drop `snapshot'
    if `compute_rc' {
        foreach _m in meta self v u initial agg va ua sa pa ben vb ub sb pb nash history trace sourcea sourceb recovery recovery_counts {
            capture matrix drop `cg'_`_m'
        }
        exit `compute_rc'
    }
    // Publish only after successful computation, including diagnostics/matrices.
    foreach _ov of local _replace_outputs {
        quietly drop `_ov'
    }
    if inlist("`strategy'","all","aggressive") rename `workagg' `generate'_agg
    if inlist("`strategy'","all","benevolent") rename `workben' `generate'_ben
    if inlist("`strategy'","all","nash") rename `worknash' `generate'_nash

    di as text ""
    di as text "{hline 70}"

    // Counts are available without matrix/detail. Diagnostics never create files.
    return scalar recovery_version=11401
    return scalar n_recovery_attempted=`cg'_recovery_counts[1,1]
    return scalar n_recovered=`cg'_recovery_counts[1,2]
    return scalar n_recovery_unresolved=`cg'_recovery_counts[1,3]
    return scalar n_recovery_limit=`cg'_recovery_counts[1,4]
    return scalar n_native_failed_agg=`cg'_recovery_counts[1,5]
    return scalar n_native_failed_ben=`cg'_recovery_counts[1,6]
    return scalar n_native_failed_game=`cg'_recovery_counts[1,7]
    return scalar n_recovered_agg=`cg'_recovery_counts[1,8]
    return scalar n_recovered_ben=`cg'_recovery_counts[1,9]
    return scalar n_recovered_game=`cg'_recovery_counts[1,10]
    return scalar n_primary_target_changed=`cg'_recovery_counts[1,11]
    return scalar max_primary_target_change=`cg'_recovery_counts[1,12]
    return scalar max_recovery_residual=`cg'_recovery_counts[1,13]
    return scalar n_recovery_rows=`cg'_recovery_counts[1,14]
    return scalar n_recovery_rows_omitted=`cg'_recovery_counts[1,15]
    if `cg'_recovery_counts[1,14]>0 {
        matrix colnames `cg'_recovery=iteration evaluator observation mode peer peer_observation native_converged native_error native_value native_rc self_initial self_certified self_delta theta recovered reason primary_cert primary_reason secondary_cert secondary_reason primary_ops secondary_ops max_limbs primary_pivots secondary_pivots export_eq export_ineq export_gap export_bounds export_error objective_value iteration_accepted
        return matrix recovery_diagnostics=`cg'_recovery
    }

    // Return values
    return local strategy "`strategy'"
    return scalar n = `n'
    return scalar tolerance=`tolerance'
    return scalar maxiter=`maxiter'
    return scalar n_failed_primary=`cg'_meta[1,1]
    if inlist("`strategy'","all","aggressive") {
        return scalar n_failed_agg=`cg'_meta[1,2]
        return scalar n_missing_agg=`cg'_meta[1,10]
    }
    if inlist("`strategy'","all","benevolent") {
        return scalar n_failed_ben=`cg'_meta[1,3]
        return scalar n_missing_ben=`cg'_meta[1,11]
    }
    if inlist("`strategy'","all","nash") {
        return scalar nash_status=`cg'_meta[1,4]
        return scalar nash_converged=(`cg'_meta[1,4]==0)
        return scalar nash_iterations=`cg'_meta[1,5]
        return scalar nash_attempts=`cg'_meta[1,6]
        return scalar nash_failed_pairs=`cg'_meta[1,7]
        return scalar nash_last_change=`cg'_meta[1,8]
        return scalar nash_initial_missing=`cg'_meta[1,9]
    }
    if "`matrix'"!="" {
        matrix colnames `cg'_v=`inputs'
        matrix colnames `cg'_u=`outputs'
        return matrix self_eff=`cg'_self
        return matrix weights_input=`cg'_v
        return matrix weights_output=`cg'_u
        return matrix ce_initial=`cg'_initial
        if inlist("`strategy'","all","aggressive") {
            matrix colnames `cg'_va=`inputs'
            matrix colnames `cg'_ua=`outputs'
            return matrix ce_agg=`cg'_agg
            return matrix weights_agg_input=`cg'_va
            return matrix weights_agg_output=`cg'_ua
            return matrix status_agg=`cg'_sa
            return matrix source_agg=`cg'_sourcea
            return matrix count_agg=`cg'_pa
        }
        if inlist("`strategy'","all","benevolent") {
            matrix colnames `cg'_vb=`inputs'
            matrix colnames `cg'_ub=`outputs'
            return matrix ce_ben=`cg'_ben
            return matrix weights_ben_input=`cg'_vb
            return matrix weights_ben_output=`cg'_ub
            return matrix status_ben=`cg'_sb
            return matrix source_ben=`cg'_sourceb
            return matrix count_ben=`cg'_pb
        }
        if inlist("`strategy'","all","nash") {
            matrix colnames `cg'_trace=iteration max_change missing_pairs complete
            return matrix ce_nash=`cg'_nash
            return matrix nash_history=`cg'_history
            return matrix nash_trace=`cg'_trace
        }
    }
    // Mata-created suffixes are not individually registered Stata tempnames.
    // Return copies above, then remove every private workspace matrix.
    foreach _m in meta self v u initial agg va ua sa pa ben vb ub sb pb nash history trace sourcea sourceb recovery recovery_counts {
        capture matrix drop `cg'_`_m'
    }
end
