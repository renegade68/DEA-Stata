*! version 1.0.0  20jul2026
*! Tolerance-grid stability diagnostics for DEA reference classifications
program define dea_reference_sensitivity, rclass
    version 16.0
    syntax [if] [in], INPUTS(varlist numeric min=1) OUTPUTS(varlist numeric min=1) ///
        [OPTGRID(numlist) PEERGRID(numlist) ///
        RTS(string) ORIentation(string) ENGINE(string) SCALING(string) ///
        DMU(varname) EVALuate(varname numeric) CANDidates(varname numeric) ///
        TOLerance(real 1e-8) SAVing(string) REPLACE TRACE]

    if "`optgrid'"=="" local optgrid "1e-6 1e-5 1e-4"
    if "`peergrid'"=="" local peergrid "1e-8 1e-7 1e-6"
    foreach z of numlist `optgrid' `peergrid' {
        if `z'<=0 {
            di as error "dea_reference_sensitivity: tolerance grids must be positive"
            exit 198
        }
    }
    if "`replace'"!="" & "`saving'"=="" {
        di as error "dea_reference_sensitivity: replace requires saving()"
        exit 198
    }
    local forward ""
    foreach o in rts orientation engine scaling dmu evaluate candidates {
        if "``o''"!="" local forward `"`forward' `o'(``o'')"'
    }
    if "`trace'"!="" local forward `"`forward' trace"'

    tempfile one all
    local scenario=0
    local total_aux=0
    foreach ot of numlist `optgrid' {
        foreach pt of numlist `peergrid' {
            local ++scenario
            quietly dea_reference `if' `in', inputs(`inputs') outputs(`outputs') ///
                `forward' tolerance(`tolerance') opttolerance(`ot') ///
                peertolerance(`pt') saving("`one'") replace
            if r(aux_failures)>0 {
                di as error "dea_reference_sensitivity: auxiliary failure in scenario `scenario'"
                exit 498
            }
            local total_aux=`total_aux'+r(aux_lp_solves)
            preserve
            quietly use "`one'", clear
            quietly keep dmu_index peer_index lambda_min lambda_max reference_status
            quietly generate int scenario=`scenario'
            quietly generate double opt_tolerance=`ot'
            quietly generate double peer_tolerance=`pt'
            quietly generate byte possible=reference_status>0 if !missing(reference_status)
            quietly generate byte necessary=reference_status==2 if !missing(reference_status)
            quietly generate byte optional=reference_status==1 if !missing(reference_status)
            if `scenario'>1 quietly append using "`all'"
            quietly save "`all'", replace
            restore
        }
    }

    preserve
    quietly use "`all'", clear
    quietly collapse (mean) possible_frequency=possible ///
        necessary_frequency=necessary optional_frequency=optional ///
        (min) status_min=reference_status lambda_min_grid=lambda_min ///
        (max) status_max=reference_status lambda_max_grid=lambda_max ///
        (count) scenarios=reference_status, by(dmu_index peer_index)
    quietly generate byte classification_changed=status_min!=status_max
    quietly generate str12 stability="stable"
    quietly replace stability="unstable" if classification_changed
    quietly count if classification_changed
    local unstable=r(N)
    quietly count if necessary_frequency==1
    local robust_necessary=r(N)
    quietly count if possible_frequency==1 & necessary_frequency<1
    local robust_possible=r(N)
    quietly count if possible_frequency==0
    local robust_impossible=r(N)
    tempname STAB
    mkmat possible_frequency necessary_frequency optional_frequency ///
        classification_changed, matrix(`STAB')
    matrix colnames `STAB'=possible_frequency necessary_frequency optional_frequency changed
    if "`saving'"!="" {
        if "`replace'"!="" quietly save `"`saving'"', replace
        else quietly save `"`saving'"'
    }
    restore

    return scalar scenarios=`scenario'
    return scalar unstable_pairs=`unstable'
    return scalar robust_necessary=`robust_necessary'
    return scalar robust_possible=`robust_possible'
    return scalar robust_impossible=`robust_impossible'
    return scalar aux_lp_solves=`total_aux'
    return matrix stability=`STAB'
    return local optgrid "`optgrid'"
    return local peergrid "`peergrid'"
    return local saving `"`saving'"'
    return local cmd "dea_reference_sensitivity"
    di as text "Reference tolerance sensitivity: " as result `scenario' ///
        as text " scenarios; unstable pairs: " as result `unstable'
end
