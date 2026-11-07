*! version 1.2.0  02nov2026
*! 2026-11-01: honor if/in for evaluation and reference samples; reject empty samples.
*! 2026-11-01: preflight output names; explicit replacevars; protect model data.
*! dea_cooperative.ado — Cooperative Game DEA (Shapley Value)
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
*!   Nakabayashi, K., & Tone, K. (2006). Egoist's dilemma: A DEA game.
*!     Omega, 34(2), 135-148. doi:10.1016/j.omega.2004.08.003
*!
*!   Lozano, S. (2012). Information sharing in DEA: A cooperative game
*!     theory approach. European Journal of Operational Research,
*!     222(3), 558-565. doi:10.1016/j.ejor.2012.05.014
*!
*!   Li, Y., & Liang, L. (2010). A Shapley value index on the importance
*!     of variables in DEA models. Expert Systems with Applications,
*!     37(9), 6287-6292. doi:10.1016/j.eswa.2010.02.093
*!
*! Changelog:
*!   1.1.1  02may2026  Return r(v_grand) and display coalition values with detail
*!   1.1.0  28apr2026  Add Form A syntax (x1 x2 = y1 y2)

*! 2026-11-02: explicit conventions, validated shared LPs, complete-sample diagnostics.

capture program drop dea_cooperative
program define dea_cooperative, rclass
    version 16.0
    _dea_parse_io `0'
    local 0 `"`r(args)'"'
    syntax [if] [in], Inputs(varlist numeric) Outputs(varlist numeric) ///
        [DMU(varname) RTS(string) GENerate(string) REPLACEVARS ///
         TOLerance(real 1e-7) MATrix DETail DEBUG]
    if "`rts'"=="" local rts "crs"
    if "`generate'"=="" local generate "_coop"
    if !inlist("`rts'","crs","vrs") {
        di as error "rts() must be crs or vrs"
        exit 198
    }
    if missing(`tolerance') | `tolerance'<=0 {
        di as error "tolerance() must be positive and finite"
        exit 198
    }
    _dea_output_check, names(`generate'_indiv `generate'_shapley `generate'_contrib `generate'_rank) protect(`inputs' `outputs' `dmu') `replacevars'
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
    foreach suffix in indiv shapley contrib rank {
        quietly generate double `generate'_`suffix'=.
    }
    di as text "Additive DEA game (input orientation; `rts')"
    di as text "Fixed full-sample frontier: v(S)=sum of member efficiencies."
    tempname gm
    mata: _dea_games_compute("`inputs'","`outputs'","`dmu'","`touse'","`generate'", ///
        "cooperative","`rts'",`tolerance',"`matrix'"!="","`detail'"!="","`debug'"!="","`gm'")
    return local rts "`rts'"
    return local orientation "input"
    return local formulation "additive_fixed_frontier"
    return scalar n=`n'
    return scalar n_lp=`n'
    return scalar tolerance=`tolerance'
    return scalar n_failed=`gm'_meta[1,1]
    return scalar n_boundary=`gm'_meta[1,2]
    return scalar complete=`gm'_meta[1,3]
    return scalar v_grand=`gm'_meta[1,4]
    return scalar detail_shown=`gm'_meta[1,5]
    if "`matrix'"!="" {
        matrix colnames `gm'_results=indiv shapley contrib rank
        matrix colnames `gm'_lp=raw_score score accepted endpoint_adjusted
        return matrix results=`gm'_results
        return matrix lp_diagnostics=`gm'_lp
    }
    foreach suffix in meta alloc results lp {
        capture matrix drop `gm'_`suffix'
    }
end
