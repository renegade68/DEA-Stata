*! version 1.3.0  09nov2026
*! 2026-11-01: use shared Form A/B parser; preserve if/in expressions.
*! 2026-11-01: preflight output names; explicit replacevars; protect model data.
*! dea_nash.ado — Nash Bargaining Game DEA
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
*!   Du, J., Liang, L., Chen, Y., Cook, W. D., & Zhu, J. (2011).
*!     A bargaining game model for measuring performance of two-stage
*!     network structures. European Journal of Operational Research,
*!     210(2), 390-397. doi:10.1016/j.ejor.2010.08.025
*!
*!   Nash, J. F. (1950). The bargaining problem. Econometrica, 18(2),
*!     155-162. doi:10.2307/1907266
*!
*!   Zhu, J. (2014). Quantitative Models for Performance Evaluation and
*!     Benchmarking: Data Envelopment Analysis with Spreadsheets (3rd ed.).
*!     Cham: Springer. doi:10.1007/978-3-319-06647-9
*!
*! Changelog:
*!   1.1.1  02may2026  Return Nash-product summaries and tighten sample handling
*!   1.1.0  28apr2026  Add Form A syntax (x1 x2 = y1 y2)

*! 2026-11-02: explicit conventions, validated shared LPs, complete-sample diagnostics.

capture program drop dea_nash
program define dea_nash, rclass
    version 16.0
    _dea_parse_io `0'
    local 0 `"`r(args)'"'
    syntax [if] [in], Inputs(varlist numeric) Outputs(varlist numeric) ///
        [DMU(varname) RTS(string) GENerate(string) REPLACEVARS ///
         TOLerance(real 1e-7) MATrix DETail DEBUG]
    if "`rts'"=="" local rts "crs"
    if "`generate'"=="" local generate "_nash"
    if !inlist("`rts'","crs","vrs") {
        di as error "rts() must be crs or vrs"
        exit 198
    }
    if missing(`tolerance') | `tolerance'<=0 {
        di as error "tolerance() must be positive and finite"
        exit 198
    }
    _dea_output_check, names(`generate'_disagree `generate'_coop `generate'_bargain `generate'_gain) protect(`inputs' `outputs' `dmu') `replacevars'
    local existing "`r(existing)'"
    marksample touse, novarlist
    markout `touse' `inputs' `outputs'
    if "`dmu'"!="" quietly replace `touse'=0 if missing(`dmu')
    quietly count if `touse'
    local n=r(N)
    if `n'==0 {
        di as error "no observations"
        exit 2000
    }
    foreach side in inputs outputs {
        tempvar positive
        quietly generate byte `positive'=0 if `touse'
        foreach v of local `side' {
            quietly count if `touse' & `v'<0
            if r(N) {
                di as error "DEA quantities must be nonnegative"
                exit 459
            }
            quietly replace `positive'=1 if `touse' & `v'>0
        }
        quietly count if `touse' & `positive'==0
        if r(N) {
            di as error "each selected observation requires a positive input and a positive output"
            exit 459
        }
    }
    // Reload the ado loader itself: an earlier installation may be cached.
    quietly findfile _dea_games_load.ado
    local gamesloader "`r(fn)'"
    capture program drop _dea_games_load
    quietly run "`gamesloader'"
    quietly _dea_games_load
    if "`dmu'"=="" {
        tempvar autoid
        quietly generate long `autoid'=_n if `touse'
        local dmu "`autoid'"
    }
    foreach v of local existing {
        quietly drop `v'
    }
    foreach suffix in disagree coop bargain gain {
        quietly generate double `generate'_`suffix'=.
    }
    di as text "DEA historical bargaining-style allocation (input orientation; `rts')"
    di as text "Heuristic scenario and capped proportional rule; not a Nash-product optimizer."
    tempname gm
    tempname gmcontext
    mata: `gmcontext'=_c57_context_save()
    mata: _dea_rec_abort(); _rsm_audit_reset(0); _rsm_success_reset(1,0); _dea_link_reset(2,0)
    capture noisily mata: _dea_games_compute("`inputs'","`outputs'","`dmu'","`touse'","`generate'", ///
        "nash","`rts'",`tolerance',"`matrix'"!="","`detail'"!="","`debug'"!="","`gm'")
    local gmrc=_rc
    mata: _c57_context_restore(`gmcontext')
    mata: mata drop `gmcontext'
    if `gmrc' exit `gmrc'
    return scalar n_certified=`gm'_summary[1,1]
    return scalar n_unverified=`gm'_summary[1,2]
    return scalar n_unavailable=`gm'_summary[1,3]
    return scalar n_native_failed=`gm'_summary[1,4]
    return scalar n_native_boundary=`gm'_summary[1,5]
    return scalar rank_quality=`gm'_summary[1,6]
    return scalar rank_reason=`gm'_summary[1,7]
    return scalar n_withheld_results=`gm'_summary[1,8]
    return scalar n_withheld_aggregates=`gm'_summary[1,9]
    return scalar n_unavailable_coalitions=`gm'_summary[1,10]
    return scalar n_numerical_coalitions=`gm'_summary[1,11]
    return scalar n_certified_coalitions=`gm'_summary[1,12]
    return scalar n_withheld_coalitions=`gm'_summary[1,13]
    return scalar allocation_quality=`gm'_summary[1,14]
    return local rts "`rts'"
    return local orientation "input"
    return local formulation "historical_capped_proportional"
    return scalar n=`n'
    return scalar n_lp=`n'
    return scalar tolerance=`tolerance'
    return scalar n_failed=`gm'_meta[1,1]
    return scalar n_boundary=`gm'_meta[1,2]
    return scalar complete=`gm'_meta[1,3]
    return scalar total_gain=`gm'_alloc[1,1]
    return scalar allocated_gain=`gm'_alloc[1,2]
    return scalar unallocated_gain=`gm'_alloc[1,3]
    return scalar mean_gain=`gm'_alloc[1,4]
    return scalar nash_product=`gm'_alloc[1,5]
    return scalar log_product=`gm'_alloc[1,6]
    return scalar n_zero_gain=`gm'_alloc[1,7]
    return scalar n_capped=`gm'_alloc[1,8]
    return scalar product_underflow=`gm'_alloc[1,9]
    return scalar allocation_complete=`gm'_alloc[1,10]
    return scalar total_gain_quality=`gm'_aggregate[1,1]
    return scalar total_gain_reason=`gm'_aggregate[1,2]
    return scalar allocated_gain_quality=`gm'_aggregate[2,1]
    return scalar allocated_gain_reason=`gm'_aggregate[2,2]
    return scalar unallocated_gain_quality=`gm'_aggregate[3,1]
    return scalar unallocated_gain_reason=`gm'_aggregate[3,2]
    return scalar mean_gain_quality=`gm'_aggregate[4,1]
    return scalar mean_gain_reason=`gm'_aggregate[4,2]
    return scalar nash_product_quality=`gm'_aggregate[5,1]
    return scalar nash_product_reason=`gm'_aggregate[5,2]
    return scalar log_product_quality=`gm'_aggregate[6,1]
    return scalar log_product_reason=`gm'_aggregate[6,2]
    if "`matrix'"!="" {
        matrix colnames `gm'_results=disagree coop bargain gain
        matrix colnames `gm'_lp=raw_score score accepted endpoint_adjusted
        return matrix results=`gm'_results
        return matrix lp_diagnostics=`gm'_lp
        matrix colnames `gm'_quality=disagree coop bargain gain
        matrix colnames `gm'_reason=disagree coop bargain gain
        matrix colnames `gm'_aggregate=quality reason value relative_error export_ops max_limbs
        matrix rownames `gm'_aggregate=total_gain allocated_gain unallocated_gain mean_gain nash_product log_product
        return matrix score_diagnostics=`gm'_scorediag
        return matrix native_lp_diagnostics=`gm'_native
        return matrix result_quality=`gm'_quality
        return matrix result_reason=`gm'_reason
        return matrix aggregate_diagnostics=`gm'_aggregate
        return matrix numerical_results=`gm'_numerical
        return matrix numerical_aggregates=`gm'_numaggregate
        if `gm'_summary[1,1]>0 {
            return matrix lambda=`gm'_lambda
            return matrix weights=`gm'_weights
        }
    }
    foreach suffix in meta alloc results lp summary aggregate scorediag native quality reason numerical numaggregate lambda weights {
        capture matrix drop `gm'_`suffix'
    }
end
