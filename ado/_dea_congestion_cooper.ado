*! version 1.1.0 06nov2026
*! Public worker; all result names are validated before computation/publication.
program define _dea_congestion_cooper, rclass
    version 16.0
    syntax [if] [in], Inputs(varlist numeric) Outputs(varlist numeric) ///
        GENerate(name) TOLerance(real) ENGine(string) CONGWeight(string) ///
        [CONGPolicy(string) CONGAnchor(string) CONGREFerence(string) CONGSelectlimit(integer 16) DMU(varname) CONGMATrix CONGRanges SAVing(string) REPLACE REPLACEVARS DEBUG]
    local _c53_dmu ""
    if "`dmu'"!="" local _c53_dmu "dmu(`dmu')"
    if "`engine'"=="simplex" {
        _dea_congestion_exact `if' `in', inputs(`inputs') outputs(`outputs') generate(`generate') tolerance(`tolerance') engine(simplex) congweight(`congweight') congpolicy(`congpolicy') conganchor(`conganchor') congreference(`congreference') congselectlimit(`congselectlimit') `congmatrix' `congranges' `_c53_dmu' saving(`saving') `replace' `replacevars' `debug'
        return add
        exit
    }
    if "`engine'"!="ipm" | "`congpolicy'`conganchor'`congreference'"!="" | `congselectlimit'!=16 exit 198
    if "`congranges'"!="" local congmatrix "congmatrix"
    local names `generate' `generate'_phi `generate'_eff `generate'_total `generate'_technical `generate'_status
    if "`dmu'"=="" {
        confirm new variable _dmu
        local reserved _dmu
        if `: list reserved in names' {
            di as error "_dmu is reserved for the automatic identifier"
            exit 110
        }
    }
    _dea_output_check, names(`names') protect(`inputs' `outputs' `dmu') `replacevars'
    local existing `r(existing)'
    marksample touse, novarlist
    markout `touse' `inputs' `outputs'
    if "`dmu'"!="" markout `touse' `dmu', strok
    quietly count if `touse'
    local n=r(N)
    if `n'==0 error 2000
    foreach v of varlist `inputs' `outputs' {
        quietly count if `touse' & `v'<=0
        if r(N) {
            di as error "congmethod(cooper) requires strictly positive inputs and outputs"
            exit 459
        }
    }
    foreach v of varlist `inputs' `outputs' {
        quietly summarize `v' if `touse', meanonly
        if r(min)/r(max)==0 {
            di as error "Cooper column scaling underflows for `v'; respecify the numerical units/range"
            exit 459
        }
    }
    capture mata: assert(_dea_cooper_weights_ready()==5401 & findexternal("_dea_cooper_weights_safe()")!=NULL)
    if _rc {
        quietly findfile _dea_cooper_weights.ado
        quietly run "`r(fn)'"
    }
    tempname weights_ok
    mata: st_numscalar("`weights_ok'",_dea_cooper_weights_safe(st_data(.,tokens("`inputs'"),"`touse'"),st_data(.,tokens("`outputs'"),"`touse'"),"`congweight'"))
    if `weights_ok'!=1 {
        di as error "Cooper objective weights overflow or underflow when scaled; respecify numerical units or congweight()"
        exit 459
    }
    local m : word count `inputs'
    local s : word count `outputs'
    if "`congmatrix'"!="" & max(`n',`m'+`s',7)>c(matsize) {
        di as error "congmatrix/congranges require matsize >= max(N, input+output count, 7)"
        exit 908
    }
    if `"`saving'"'!="" {
        if !regexm(lower(`"`saving'"'), "\.dta$") local saving `"`saving'.dta"'
        if "`replace'"=="" confirm new file `"`saving'"'
    }
    capture mata: assert(_dea_cooper_ready()==1)
    if _rc {
        quietly findfile _dea_cooper.ado
        quietly run "`r(fn)'"
    }
    local resultvars ""
    forvalues j=1/6 {
        tempvar v`j'
        quietly generate double `v`j''=.
        local resultvars `resultvars' `v`j''
    }
    local destinations ""
    forvalues j=1/16 {
        tempname mat`j'
        local destinations `destinations' `mat`j''
    }
    mata: _dea_cooper_compute("`inputs'","`outputs'","`touse'","`resultvars'", ///
        "`congweight'",`tolerance',"`congmatrix'"!="","`congranges'"!="", ///
        "`debug'"!="","`destinations'")

    label variable `v1' "Cooper mean relative input congestion (0 = none)"
    label variable `v2' "BCC output expansion phi"
    label variable `v3' "Reciprocal BCC output expansion (1/phi)"
    label variable `v4' "Mean relative total input slack"
    label variable `v5' "Mean relative technical input slack"
    label variable `v6' "0 complete; 1 radial/2 slack/3 decomposition LP unavailable"

    // Save before publishing in memory. A failed save restores existing data.
    if `"`saving'"'!="" {
        preserve
        foreach v of local existing {
            quietly drop `v'
        }
        forvalues j=1/6 {
            local name : word `j' of `names'
            rename `v`j'' `name'
        }
        local fileid `dmu'
        if "`dmu'"=="" {
            quietly generate long _dmu=_n if `touse'
            local fileid _dmu
        }
        quietly keep if `touse'
        keep `fileid' `inputs' `outputs' `names'
        capture noisily save `"`saving'"', `replace'
        local saverc=_rc
        restore
        if `saverc' exit `saverc'
    }
    foreach v of local existing {
        quietly drop `v'
    }
    forvalues j=1/6 {
        local name : word `j' of `names'
        rename `v`j'' `name'
    }
    if "`dmu'"=="" quietly generate long _dmu=_n if `touse'
    local failed=`mat1'[1,2]
    di as text "Cooper hybrid congestion (VRS; output expansion; `congweight' slack weights)"
    di as text "  Observations: `n'; complete decompositions: " as result `n'-`failed'
    if `failed' di as text "  Unavailable LP stages: `failed' DMUs; see `generate'_status."
    return scalar N=`n'
    return scalar N_complete=`n'-`failed'
    return scalar n_failed=`failed'
    return scalar n_efficient=`mat1'[1,3]
    return scalar N_eff=`mat1'[1,3]
    return scalar n_congested=`mat1'[1,4]
    return scalar n_endpoint_adjusted=`mat1'[1,5]
    return scalar n_native_retry=`mat1'[2,1]
    return scalar n_native_recovered=`mat1'[2,2]
    return scalar tolerance=`tolerance'
    return scalar acceptance_tolerance=min(`tolerance',1e-9)
    foreach pair in congestion phi {
        local var `generate'
        if "`pair'"=="phi" local var `generate'_phi
        quietly summarize `var' if `touse', meanonly
        return scalar mean_`pair'=r(mean)
    }
    return local formulation "cooper_two_stage_vrs"
    return local lp_backend "stata_linearprogram"
    return local actual_engine "ipm"
    return local actual_backend "stata_linearprogram_with_native_retry"
    return local congmethod "cooper"
    return local congweight "`congweight'"
    return local orientation "output"
    return local rts "vrs"
    return local engine "`engine'"
    return local projection_policy "selected_slack_optimum"
    if "`congmatrix'"!="" {
        matrix colnames `mat10'=radial slack decomposition adjusted_coordinates projection_adjustment equivalent_retry recovered_retry
        matrix colnames `mat11'=observation
        local matnames input_total input_technical input_congestion output_slack target_inputs target_outputs lambda_stage1 lambda_stage2 status obs_index
        forvalues j=2/11 {
            local k=`j'-1
            local name : word `k' of `matnames'
            if inlist(`j',2,3,4,6) matrix colnames `mat`j''=`inputs'
            if inlist(`j',5,7) matrix colnames `mat`j''=`outputs'
            return matrix `name'=`mat`j''
        }
    }
    if "`congranges'"!="" {
        return scalar n_projection_multiple=`mat1'[1,6]
        return scalar n_decomposition_multiple=`mat1'[1,7]
        return scalar n_range_failed=`mat1'[1,8]
        return scalar n_range_complete=`mat1'[1,9]
        return local range_scope "phase1_optimal_face; phase2_conditional_on_selected_projection"
        return matrix projection_min=`mat12'
        return matrix projection_max=`mat13'
        matrix colnames `mat14'=`inputs'
        matrix colnames `mat15'=`inputs'
        return matrix congestion_min=`mat14'
        return matrix congestion_max=`mat15'
        return matrix range_status=`mat16'
    }
end
