*! version 1.2.0  02nov2026
*! 2026-11-01: accept the minimum option abbreviations advertised in help.
*! 2026-11-01: use shared Form A/B parser; preserve if/in expressions.
*! 2026-11-01: preflight output names; explicit replacevars; protect model data.
*! dea_network.ado — Network DEA (Two-Stage)
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
*!   Fare, R., & Grosskopf, S. (2000). Network DEA. Socio-Economic
*!     Planning Sciences, 34(1), 35-49. doi:10.1016/S0038-0121(99)00012-9
*!
*!   Kao, C., & Hwang, S.-N. (2008). Efficiency decomposition in
*!     two-stage data envelopment analysis: An application to non-life
*!     insurance companies in Taiwan. European Journal of Operational
*!     Research, 185(1), 418-429. doi:10.1016/j.ejor.2006.11.041
*!
*!   Cooper, W. W., Seiford, L. M., & Tone, K. (2007). Data Envelopment
*!     Analysis: A Comprehensive Text with Models, Applications, References
*!     and DEA-Solver Software (2nd ed.). New York: Springer.
*!     doi:10.1007/978-0-387-45283-8
*!
*! Changelog:
*!   1.1.2  03may2026  Implement Kao-Hwang CRS/input shared intermediate-weight
*!                     relational multiplier model for method(multiplicative)
*!   1.1.1  02may2026  Fix string DMU handling and pass orientation to stage LPs
*!   1.1.0  11feb2026  Unified syntax: dea_network x = y | inputs()/outputs()
*!   1.0.0  06feb2026  Initial packaged version; based on earlier teaching/seminar versions

*! 2026-11-02: validated/scaled LPs, explicit formulations, diagnostics and matrices.

capture program drop dea_network
program define dea_network, rclass
    version 16.0
    _dea_parse_io `0'
    local 0 `"`r(args)'"'

    syntax [varlist(numeric default=none)] [if] [in], ///
        INTermediates(varlist numeric) ///
        [Inputs(varlist numeric) ///
         Outputs(varlist numeric) ///
         DMU(varname) ///
         Rts(string) ///
         ORT(string) ORientation(string) ///
         METHOD(string) ///
         GENerate(string) ///
         REPLACEVARS      ///
         TOLerance(real 1e-7) ///
         MATrix ///
         DETail ///
         DEBUG]
    
    // Resolve outputs
    if "`outputs'" == "" {
        if "`varlist'" != "" {
            local outputs `varlist'
        }
        else {
            di as error "outputs() required, or use: dea_network inputs = outputs, intermediates() [options]"
            exit 198
        }
    }
    if "`inputs'" == "" {
        di as error "inputs() required, or use: dea_network inputs = outputs, intermediates() [options]"
        exit 198
    }
    if "`dmu'" != "" {
        confirm variable `dmu'
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
    if "`method'" == "" local method "multiplicative"
    if "`generate'" == "" local generate "_net_eff"
    
    // Abbreviations for orientation
    if inlist("`orientation'", "i", "in")  local orientation "input"
    if inlist("`orientation'", "o", "out") local orientation "output"
    
    // Validate
    if !inlist("`rts'", "crs", "vrs") {
        di as error "rts() must be crs or vrs"
        exit 198
    }
    if !inlist("`orientation'", "input", "output") {
        di as error "orientation() must be input or output"
        exit 198
    }
    if !inlist("`method'", "multiplicative", "additive", "independent") {
        di as error "method() must be multiplicative, additive, or independent"
        exit 198
    }
    
    local formulation "independent_product"
    local legacy_fallback=0
    if "`method'"=="additive" local formulation "independent_mean"
    if "`method'"=="multiplicative" {
        if "`rts'"=="crs" & "`orientation'"=="input" local formulation "relational_crs"
        else local legacy_fallback=1
    }
    local mode "product"
    if "`formulation'"=="relational_crs" local mode "relational"
    if "`formulation'"=="independent_mean" local mode "mean"
    local leader ""
    if missing(`tolerance') | `tolerance'<=0 {
        di as error "tolerance() must be positive and finite"
        exit 198
    }
    capture mata: assert(_dea_ts_ready()==1)
    if _rc {
        quietly findfile _dea_twostage.ado
        quietly run "`r(fn)'"
    }
    _dea_output_check, names(`generate'_overall `generate'_stage1 `generate'_stage2) ///
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
    quietly generate double `generate'_overall=.
    quietly generate double `generate'_stage1=.
    quietly generate double `generate'_stage2=.
    di as text "Two-stage DEA: " as result "`formulation'" as text " (RTS: `rts'; orientation: `orientation')"
    if `legacy_fallback' di as text "Note: this multiplicative option combination uses independent stages; no shared-weight network constraint."
    tempname ts
    mata: _dea_ts_compute("`inputs'","`intermediates'","`outputs'","`dmu'","`touse'", ///
        "`generate'","network","`mode'","`leader'","`rts'","`orientation'",`tolerance', ///
        "`matrix'"!="","`detail'"!="","`debug'"!="","`ts'")
    quietly summarize `generate'_overall if `touse'
    local mean_overall=r(mean)
    quietly summarize `generate'_stage1 if `touse'
    local mean_stage1=r(mean)
    quietly summarize `generate'_stage2 if `touse'
    local mean_stage2=r(mean)
    return local method "`method'"
    return local rts "`rts'"
    return local orientation "`orientation'"
    return local formulation "`formulation'"
    return scalar n=`n'
    return scalar tolerance=`tolerance'
    return scalar n_failed_overall=`ts'_meta[1,1]
    return scalar n_failed_stage1=`ts'_meta[1,2]
    return scalar n_failed_stage2=`ts'_meta[1,3]
    return scalar legacy_fallback=`legacy_fallback'
    return scalar n_decomp_fallback=`ts'_meta[1,4]
    return scalar mean_overall=`mean_overall'
    return scalar mean_stage1=`mean_stage1'
    return scalar mean_stage2=`mean_stage2'
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
            else return matrix follower_status=`ts'_joint
        }
        if "`mode'"!="relational" {
            return matrix lambda_stage1=`ts'_l1
            return matrix lambda_stage2=`ts'_l2
        }
    }
    foreach suffix in meta scores status v w u joint l1 l2 {
        capture matrix drop `ts'_`suffix'
    }
end
