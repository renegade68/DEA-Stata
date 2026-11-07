*! version 1.0.0 06nov2026
capture program drop _dea_congestion_exact
capture program drop _dea_congestion_exact_run
program define _dea_congestion_exact, rclass
    version 16.0
    capture mata: assert(_c53_context_ready()==53)
    if _rc {
        quietly findfile _dea_cooper_context.ado
        quietly run "`r(fn)'"
    }
    tempname snapshot
    mata: `snapshot'=_c53_context_save()
    capture noisily _dea_congestion_exact_run `0' snapshot(`snapshot')
    local result_rc=_rc
    mata: _c53_context_restore(`snapshot')
    mata: mata drop `snapshot'
    if `result_rc' exit `result_rc'
    return add
end
program define _dea_congestion_exact_run, rclass
    version 16.0
    syntax [if] [in], Inputs(varlist numeric) Outputs(varlist numeric) ///
        GENerate(name) TOLerance(real) ENGine(string) CONGWeight(string) ///
        [SNAPshot(string) CONGPolicy(string) CONGAnchor(string) CONGREFerence(string) CONGSelectlimit(integer 16) AUDit(string) DMU(varname) CONGMATrix CONGRanges SAVing(string) REPLACE REPLACEVARS DEBUG]
    if "`engine'"!="simplex" exit 198
    if !inlist("`congweight'","raw","normalized") exit 198
    if "`congpolicy'"=="" local congpolicy "stage1"
    if !inlist("`congpolicy'","stage1","stage2") | !inlist(`congselectlimit',16,128) exit 198
    if ("`congpolicy'"=="stage2" & "`congreference'"=="") | ("`congpolicy'"!="stage2" & "`congreference'"!="") exit 198
    foreach ref in conganchor congreference {
        if "``ref''"!="" {
            capture confirm name ``ref''
            if _rc exit 198
            confirm matrix ``ref''
        }
    }
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
    local m : word count `inputs'
    local s : word count `outputs'
    if 2*`n'+`m'+`s'+3>64 {
        di as error "bounded Cooper certificate requires 2*N+inputs+outputs+3<=64; use engine(ipm) for larger models"
        exit 198
    }
    foreach ref in conganchor congreference {
        if "``ref''"!="" {
            if rowsof(``ref'')!=`n' | colsof(``ref'')!=`n' exit 198
        }
    }
    if "`congmatrix'"!="" & max(`n',`m'+`s',18)>c(matsize) {
        di as error "congmatrix/congranges require matsize >= max(N, input+output count, 18)"
        exit 908
    }
    if `"`saving'"'!="" {
        if !regexm(lower(`"`saving'"'), "\.dta$") local saving `"`saving'.dta"'
        if "`replace'"=="" confirm new file `"`saving'"'
    }
    capture mata: assert(_c53_data_ready()==5401)
    if _rc {
        quietly findfile _dea_cooper_data.ado
        quietly run "`r(fn)'"
    }
    mata: _dea_rec_abort(); _rsm_audit_reset(0); _rsm_success_reset(1,0); _dea_link_reset(2,0)
    mata: _c53_context_defaults(`snapshot')
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
    if `"`audit'"'!="" {
        forvalues j=1/`n' {
            quietly mkdir "`audit'/dmu_`j'"
        }
    }
    local resultvars ""
    forvalues j=1/6 {
        tempvar v`j'
        quietly generate double `v`j''=.
        local resultvars `resultvars' `v`j''
    }
    local destinations ""
    forvalues j=1/17 {
        tempname mat`j'
        local destinations `destinations' `mat`j''
    }
    mata: _c53_compute("`inputs'","`outputs'","`touse'","`resultvars'", ///
        "`congweight'",`tolerance',"`congmatrix'"!="","`congranges'"!="", ///
        "`debug'"!="","`destinations'","`congpolicy'","`conganchor'","`congreference'",`congselectlimit',"`audit'")


    label variable `v1' "Cooper mean relative input congestion (0 = none)"
    label variable `v2' "BCC output expansion phi"
    label variable `v3' "Reciprocal BCC output expansion (1/phi)"
    label variable `v4' "Mean relative total input slack"
    label variable `v5' "Mean relative technical input slack"
    label variable `v6' "0 complete; 1 radial/2 slack-anchor/3 decomposition/4 numeric export unavailable"

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
    return scalar acceptance_tolerance=0
    return local tolerance_role "summary_classification_only"
    foreach pair in congestion phi {
        local var `generate'
        if "`pair'"=="phi" local var `generate'_phi
        quietly summarize `var' if `touse', meanonly
        return scalar mean_`pair'=r(mean)
    }
    return local formulation "cooper_two_stage_vrs"
    return local lp_backend "bounded_integer_primal_dual_after_rsm"
    return local actual_engine "simplex"
    return local actual_backend "rsm_proposal_and_bounded_integer_basis_repair"
    return local certificate_scope "internal_exact_point"
    return local certificate_engine "bounded_integer_primal_dual"
    return local reference_backend "explicit_matrix_or_unanchored"
    return local congpolicy "`congpolicy'"
    return scalar selection_maxpivots=`congselectlimit'
    return scalar n_selection_not_requested=`mat1'[3,1]
    return scalar n_selection_certified=`mat1'[3,2]
    return scalar n_selection_missing_reference=`mat1'[3,3]
    return scalar n_selection_uncertified=`mat1'[3,4]
    return scalar n_selection_mismatch=`mat1'[3,5]
    return scalar n_conversion_failed=`mat1'[3,6]
    return scalar n_unanchored=`mat1'[3,7]
    return local integration_status "integrated"

    return local congmethod "cooper"
    return local congweight "`congweight'"
    return local orientation "output"
    return local rts "vrs"
    return local engine "`engine'"
    return local projection_policy "exact_slack_vertex"
    if "`conganchor'"!="" return local projection_policy "nearest_reference_stage1_on_exact_slack_face"
    return local decomposition_policy "`congpolicy'"

    if "`congmatrix'"!="" {
        matrix colnames `mat17'=usable base_certificate selection_status selection_applied selection_certificate conversion_status rounded_values selection_maxpivots reference_distance ranges_complete range_conversion_status range_rounded_values radial slack anchor decomposition projection congestion
        return matrix certificate_status=`mat17'
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
        return local range_scope "projection_on_slack_face_congestion_on_original_technical_face"
        return matrix projection_min=`mat12'
        return matrix projection_max=`mat13'
        matrix colnames `mat14'=`inputs'
        matrix colnames `mat15'=`inputs'
        return matrix congestion_min=`mat14'
        return matrix congestion_max=`mat15'
        matrix colnames `mat16'=projection congestion
        return matrix range_status=`mat16'
    }
end
