*! version 1.2.0  02nov2026
*! 2026-11-01: preflight output names; explicit replacevars; protect model data.
*! dea_stackelberg.ado — Stackelberg Game DEA
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
*!   Liang, L., Cook, W. D., & Zhu, J. (2008). DEA models for two-stage
*!     processes: Game approach and efficiency decomposition. Naval Research
*!     Logistics, 55(7), 643-653. doi:10.1002/nav.20308
*! -----------
*!   Liang, L., Yang, F., Cook, W. D., & Zhu, J. (2006). DEA models for
*!     supply chain efficiency evaluation. Annals of Operations Research,
*!     145(1), 35-49. doi:10.1007/s10479-006-0026-7
*!
*!   Li, Y., Chen, Y., Liang, L., & Xie, J. (2012). DEA models for
*!     extended two-stage network structures. Omega, 40(5), 611-618.
*!     doi:10.1016/j.omega.2011.11.007
*!
*! Changelog:
*!   1.1.1  02may2026  Tighten help alignment and support string DMU identifiers
*!   1.1.0  28apr2026  Remove spurious [varlist] from syntax

*! 2026-11-02: validated/scaled LPs, explicit formulations, diagnostics and matrices.

capture program drop dea_stackelberg
program define dea_stackelberg, rclass
    version 16.0
    syntax [if] [in], ///
        Inputs(varlist numeric) ///
        INTERmediates(varlist numeric) ///
        Outputs(varlist numeric) ///
        [DMU(varname) ///
         LEADer(string) ///
         RTS(string) ///
         METHOD(string) ///
         MATrix ///
         GENerate(string) ///
         REPLACEVARS      ///
         TOLerance(real 1e-7) ///
         DETail ///
         DEBUG]
    
    // Defaults
    if "`leader'" == "" local leader "stage1"
    local leader = lower("`leader'")
    if "`rts'" == "" local rts "crs"
    local rts = lower("`rts'")
    if "`generate'" == "" local generate "_stack"
    
    // Validate
    if !inlist("`rts'", "crs", "vrs") {
        di as error "rts() must be crs or vrs"
        exit 198
    }
    if !inlist("`leader'", "stage1", "stage2") {
        di as error "leader() must be stage1 or stage2"
        exit 198
    }
    
    if "`method'"=="" local method "independent"
    if !inlist("`method'","independent","noncooperative") {
        di as error "method() must be independent or noncooperative"
        exit 198
    }
    if "`method'"=="noncooperative" & "`rts'"!="crs" {
        di as error "method(noncooperative) currently requires rts(crs)"
        exit 198
    }
    local orientation "input"
    local mode "`method'"
    local formulation "independent_ordering"
    if "`method'"=="noncooperative" local formulation "noncooperative_crs"
    if missing(`tolerance') | `tolerance'<=0 {
        di as error "tolerance() must be positive and finite"
        exit 198
    }
    capture mata: assert(_dea_ts_ready()==1)
    if _rc {
        quietly findfile _dea_twostage.ado
        quietly run "`r(fn)'"
    }
    _dea_output_check, names(`generate'_leader `generate'_follower `generate'_overall `generate'_indep1 `generate'_indep2) ///
        protect(`inputs' `intermediates' `outputs' `dmu') `replacevars'
    local _replace_outputs "`r(existing)'"
    marksample touse, novarlist
    markout `touse' `inputs' `intermediates' `outputs'
    if "`dmu'"!="" quietly replace `touse'=0 if missing(`dmu')
    quietly count if `touse'
    local n=r(N)
    if `n'==0 {
        di as error "no observations"
        exit 2000
    }
    foreach _side in inputs intermediates outputs {
        tempvar positive
        quietly generate byte `positive'=0 if `touse'
        foreach v of local `_side' {
            quietly count if `touse' & `v'<0
            if r(N) {
                di as error "two-stage quantities must be nonnegative"
                exit 459
            }
            quietly replace `positive'=1 if `touse' & `v'>0
        }
        quietly count if `touse' & `positive'==0
        if r(N) {
            di as error "each selected observation requires a positive input, intermediate, and output"
            exit 459
        }
    }
    if "`dmu'"=="" {
        tempvar dmu_temp
        quietly generate long `dmu_temp'=_n if `touse'
        local dmu "`dmu_temp'"
    }
    foreach v of local _replace_outputs {
        quietly drop `v'
    }
    quietly generate double `generate'_leader=.
    quietly generate double `generate'_follower=.
    quietly generate double `generate'_overall=.
    quietly generate double `generate'_indep1=.
    quietly generate double `generate'_indep2=.
    di as text "Two-stage DEA: " as result "`formulation'" as text " (RTS: `rts'; orientation: `orientation')"
    di as text "Leader: `leader'"
    if "`method'"=="independent" di as text "Note: independent stage ordering; leader() changes labels only."
    tempname ts
    mata: _dea_ts_compute("`inputs'","`intermediates'","`outputs'","`dmu'","`touse'", ///
        "`generate'","stack","`mode'","`leader'","`rts'","`orientation'",`tolerance', ///
        "`matrix'"!="","`detail'"!="","`debug'"!="","`ts'")
    return local method "`method'"
    return local rts "`rts'"
    return local orientation "`orientation'"
    return local formulation "`formulation'"
    return scalar n=`n'
    return scalar tolerance=`tolerance'
    return scalar n_failed_overall=`ts'_meta[1,1]
    return scalar n_failed_stage1=`ts'_meta[1,2]
    return scalar n_failed_stage2=`ts'_meta[1,3]
    return local leader "`leader'"
    return scalar n_failed_follower_lp=`ts'_meta[1,5]
    return scalar n_failed_indep1=`ts'_meta[1,6]
    return scalar n_failed_indep2=`ts'_meta[1,7]
    return scalar n_single_link=`ts'_meta[1,8]
    if "`matrix'"!="" {
        matrix colnames `ts'_scores=overall stage1 stage2
        matrix colnames `ts'_status=overall stage1 stage2
        return matrix efficiencies=`ts'_scores
        return matrix status=`ts'_status
        if inlist("`mode'","relational","noncooperative") {
            matrix colnames `ts'_v=`inputs'
            matrix colnames `ts'_w=`intermediates'
            matrix colnames `ts'_u=`outputs'
            return matrix weights_input=`ts'_v
            return matrix weights_intermediate=`ts'_w
            return matrix weights_output=`ts'_u
            if "`mode'"=="relational" return matrix decomposition_status=`ts'_joint
            else {
                return matrix follower_status=`ts'_joint
                return matrix follower_route=`ts'_route
            }
        }
        if "`mode'"!="relational" {
            return matrix lambda_stage1=`ts'_l1
            return matrix lambda_stage2=`ts'_l2
        }
    }
    foreach suffix in meta scores status v w u joint l1 l2 route {
        capture matrix drop `ts'_`suffix'
    }
end
