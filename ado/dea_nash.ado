*! version 1.2.0  02nov2026
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
    capture mata: assert(_dea_games_ready()==1)
    if _rc {
        quietly findfile _dea_games.ado
        quietly run "`r(fn)'"
    }
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
    mata: _dea_games_compute("`inputs'","`outputs'","`dmu'","`touse'","`generate'", ///
        "nash","`rts'",`tolerance',"`matrix'"!="","`detail'"!="","`debug'"!="","`gm'")
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
    if "`matrix'"!="" {
        matrix colnames `gm'_results=disagree coop bargain gain
        matrix colnames `gm'_lp=raw_score score accepted endpoint_adjusted
        return matrix results=`gm'_results
        return matrix lp_diagnostics=`gm'_lp
    }
    foreach suffix in meta alloc results lp {
        capture matrix drop `gm'_`suffix'
    }
end
