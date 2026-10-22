*! version 1.14.0  20jul2026
*! Optimal-set diagnostics for radial two-stage DEA reference weights
program define dea_reference, rclass
    version 16.0
    syntax [if] [in], INPUTS(varlist numeric min=1) ///
        OUTPUTS(varlist numeric min=1) ///
        [RTS(string) ORIentation(string) ENGINE(string) SCALING(string) ///
         DMU(varname) EVALuate(varname numeric) CANDidates(varname numeric) ///
         LEXicographic(string) LEXOrder(varname numeric) ///
         POLICYWeight(varname numeric) POLICYMatrix(name) POLICYSense(string) ///
         SPARSE BROAD RELativeinterior BASIC AUDIT RTSONLY COMPARE INTERCEPTBounds ///
         MAXNodes(integer 10000) MAXAUXLP(integer -1) INTTol(real 1e-7) ///
         TOLerance(real 1e-8) PEERTOLerance(real 1e-7) ///
         SLACKNormalization(string) RADIALTolerance(real -1) ///
         OPTTOLerance(real -1) GENerate(name) SAVing(string) REPLACE TRACE]

    if "`rts'" == "" local rts "vrs"
    local rts = lower("`rts'")
    if "`rts'" == "nirs" local rts "drs"
    if "`rts'" == "ndrs" local rts "irs"
    if !inlist("`rts'", "crs", "vrs", "drs", "irs") {
        di as error "dea_reference: rts() must be crs, vrs, drs, or irs"
        exit 198
    }
    if "`orientation'" == "" local orientation "input"
    local orientation = lower("`orientation'")
    if inlist("`orientation'", "i", "in", "inp", "inpu") local orientation "input"
    if inlist("`orientation'", "o", "out", "outp", "outpu") local orientation "output"
    if !inlist("`orientation'", "input", "output") {
        di as error "dea_reference: orientation() must be input or output"
        exit 198
    }
    if "`engine'" == "" local engine "auto"
    local engine = lower("`engine'")
    if !inlist("`engine'", "simplex", "ipm", "auto") {
        di as error "dea_reference: engine() must be simplex, ipm, or auto"
        exit 198
    }
    if "`scaling'" == "" local scaling "row"
    local scaling = lower("`scaling'")
    if !inlist("`scaling'", "none", "row") {
        di as error "dea_reference: scaling() must be none or row"
        exit 198
    }
    if `tolerance' <= 0 | `peertolerance' <= 0 {
        di as error "dea_reference: tolerances must be positive"
        exit 198
    }
    if `opttolerance' == 0 | `opttolerance' < -1 {
        di as error "dea_reference: opttolerance() must be positive or omitted"
        exit 198
    }
    if "`slacknormalization'" == "" local slacknormalization "range"
    local slacknormalization = lower("`slacknormalization'")
    if !inlist("`slacknormalization'", "raw", "range") {
        di as error "dea_reference: slacknormalization() must be raw or range"
        exit 198
    }
    if `radialtolerance' < -1 {
        di as error "dea_reference: radialtolerance() must be nonnegative or omitted"
        exit 198
    }
    if `radialtolerance' < 0 local radialtolerance = 100*`tolerance'
    if "`replace'" != "" & "`saving'" == "" {
        di as error "dea_reference: replace requires saving()"
        exit 198
    }
    if "`rtsonly'" != "" {
        if "`rts'" != "crs" {
            di as error "dea_reference: rtsonly requires rts(crs)"
            exit 198
        }
        if "`candidates'`lexicographic'`lexorder'`policyweight'`policymatrix'`sparse'`broad'`relativeinterior'`basic'`audit'`saving'" != "" {
            di as error "dea_reference: rtsonly may not be combined with peer-selection or saving options"
            exit 198
        }
    }
    if "`generate'" != "" & "`rtsonly'" == "" {
        di as error "dea_reference: generate() currently requires rtsonly"
        exit 198
    }
    if "`compare'" != "" & "`rtsonly'" == "" {
        di as error "dea_reference: compare requires rtsonly"
        exit 198
    }
    if "`interceptbounds'" != "" & "`rtsonly'" == "" {
        di as error "dea_reference: interceptbounds requires rtsonly"
        exit 198
    }
    local lexicographic = lower("`lexicographic'")
    if "`lexicographic'" != "" & !inlist("`lexicographic'", "min", "max") {
        di as error "dea_reference: lexicographic() must be min or max"
        exit 198
    }
    if "`lexorder'" != "" & "`lexicographic'" == "" {
        di as error "dea_reference: lexorder() requires lexicographic()"
        exit 198
    }
    if "`policyweight'" != "" & "`policymatrix'" != "" {
        di as error "dea_reference: specify only one of policyweight() or policymatrix()"
        exit 198
    }
    if "`policysense'" != "" & "`policyweight'`policymatrix'" == "" {
        di as error "dea_reference: policysense() requires policyweight() or policymatrix()"
        exit 198
    }
    if "`policysense'" == "" local policysense "min"
    local policysense = lower("`policysense'")
    if !inlist("`policysense'", "min", "max") {
        di as error "dea_reference: policysense() must be min or max"
        exit 198
    }
    if "`audit'" != "" {
        local sparse "sparse"
        local broad "broad"
        local relativeinterior "relativeinterior"
        local basic "basic"
    }
    if "`sparse'" != "" & "`candidates'" != "" {
        di as error "dea_reference: sparse may not be combined with candidates()"
        exit 198
    }
    if "`broad'" != "" & "`candidates'" != "" {
        di as error "dea_reference: broad may not be combined with candidates()"
        exit 198
    }
    if "`relativeinterior'" != "" & "`candidates'" != "" {
        di as error "dea_reference: relativeinterior may not be combined with candidates()"
        exit 198
    }
    if `maxnodes' < 1 | `inttol' <= 0 {
        di as error "dea_reference: maxnodes() and inttol() must be positive"
        exit 198
    }
    if `maxauxlp'==0 | `maxauxlp'<-1 {
        di as error "dea_reference: maxauxlp() must be positive or -1 for unlimited"
        exit 198
    }

    marksample touse, novarlist
    markout `touse' `inputs' `outputs'
    if "`lexorder'" != "" markout `touse' `lexorder'
    if "`policyweight'" != "" {
        quietly count if `touse' & missing(`policyweight')
        if r(N)>0 {
            di as error "dea_reference: policyweight() must be complete on the marked technology sample"
            exit 459
        }
    }
    if "`dmu'" != "" {
        capture confirm string variable `dmu'
        if !_rc markout `touse' `dmu', strok
        else markout `touse' `dmu'
    }
    quietly count if `touse'
    if r(N)==0 {
        di as error "dea_reference: no complete observations"
        exit 2000
    }
    tempvar evaluse canduse
    quietly generate byte `evaluse' = `touse'
    quietly generate byte `canduse' = `touse'
    if "`evaluate'" != "" quietly replace `evaluse' = ///
        `touse' & !missing(`evaluate') & `evaluate'!=0
    if "`candidates'" != "" quietly replace `canduse' = ///
        `touse' & !missing(`candidates') & `candidates'!=0
    quietly count if `evaluse'
    local N_eval = r(N)
    quietly count if `canduse'
    local N_cand = r(N)
    if `N_eval'==0 | `N_cand'==0 {
        di as error "dea_reference: evaluate() and candidates() must each select at least one observation"
        exit 2000
    }

    capture findfile _dea_tecontract.ado
    if _rc {
        di as error "dea_reference: _dea_tecontract.ado not found"
        exit 601
    }
    tempvar efficiency refstatus
    quietly _dea_tecontract if `touse', inputs(`inputs') outputs(`outputs') ///
        rts(`rts') orientation(`orientation') stage(2) engine(`engine') ///
        scaling(`scaling') tolerance(`tolerance') ///
        slacknormalization(`slacknormalization') radialtolerance(`radialtolerance') ///
        generate(`efficiency') status(`refstatus') `trace'

    local N = r(N)
    local compare_lp_planned = cond("`compare'"!="",2*`N',0)
    if r(n_failed)>0 | r(stage2_failures)>0 {
        di as error "dea_reference: baseline two-stage DEA failed for one or more DMUs"
        exit 498
    }
    tempname RAD TS SEL SW RADADJ PX PY BCCPX BCCPY
    matrix `RAD' = r(radial_stage2)
    matrix `RADADJ' = r(radial_adjustment)
    matrix `TS' = r(total_slack)
    matrix `SEL' = r(lambda)
    matrix `SW' = r(slack_weights)
    matrix `PX' = r(projected_inputs)
    matrix `PY' = r(projected_outputs)
    if "`interceptbounds'" != "" {
        tempvar __bcc_eff __bcc_status
        quietly _dea_tecontract if `touse', inputs(`inputs') outputs(`outputs') ///
            rts(vrs) orientation(`orientation') stage(2) engine(`engine') ///
            scaling(`scaling') tolerance(`tolerance') ///
            slacknormalization(`slacknormalization') radialtolerance(`radialtolerance') ///
            generate(`__bcc_eff') status(`__bcc_status') `trace'
        if r(n_failed)>0 | r(stage2_failures)>0 {
            di as error "dea_reference: BCC projection failed for intercept bounds"
            exit 498
        }
        matrix `BCCPX' = r(projected_inputs)
        matrix `BCCPY' = r(projected_outputs)
    }
    else {
        matrix `BCCPX' = `PX'
        matrix `BCCPY' = `PY'
    }
    tempname VRSRAD NIRSRAD
    if "`compare'" != "" {
        tempvar __rts_vrs_eff __rts_vrs_status __rts_nirs_eff __rts_nirs_status
        quietly _dea_tecontract if `touse', inputs(`inputs') outputs(`outputs') ///
            rts(vrs) orientation(`orientation') stage(1) engine(`engine') ///
            scaling(`scaling') tolerance(`tolerance') ///
            generate(`__rts_vrs_eff') status(`__rts_vrs_status') `trace'
        if r(n_failed)>0 {
            di as error "dea_reference: standard VRS comparison failed"
            exit 498
        }
        matrix `VRSRAD' = r(radial)
        quietly _dea_tecontract if `touse', inputs(`inputs') outputs(`outputs') ///
            rts(drs) orientation(`orientation') stage(1) engine(`engine') ///
            scaling(`scaling') tolerance(`tolerance') ///
            generate(`__rts_nirs_eff') status(`__rts_nirs_status') `trace'
        if r(n_failed)>0 {
            di as error "dea_reference: standard NIRS comparison failed"
            exit 498
        }
        matrix `NIRSRAD' = r(radial)
    }
    local m : word count `inputs'
    local s : word count `outputs'
    local nv = `N' + `m' + `s'
    local nineq = 2 + inlist("`rts'", "drs", "irs")
    local planned = cond("`rtsonly'"!="", 0, 2*`N_eval'*`N_cand')
    if "`rts'" == "crs" local planned = `planned' + 2*`N_eval'
    if "`interceptbounds'" != "" local planned = `planned' + 2*`N_eval'
    if "`lexicographic'" != "" local planned = `planned' + `N_eval'*`N'
    if "`relativeinterior'" != "" local planned = `planned' + `N_eval'*(1+`nv'+`nineq')
    if "`basic'" != "" local planned = `planned' + `N_eval'
    if "`policyweight'`policymatrix'" != "" local planned = `planned' + `N_eval'
    if `maxauxlp'>0 & `planned'>`maxauxlp' {
        di as error "dea_reference: planned auxiliary LP count `planned' exceeds maxauxlp(`maxauxlp')"
        exit 498
    }
    if "`policymatrix'" != "" {
        capture confirm matrix `policymatrix'
        if _rc {
            di as error "dea_reference: policymatrix() matrix not found"
            exit 111
        }
        if !inlist(rowsof(`policymatrix'), 1, `N') | colsof(`policymatrix')!=`N' {
            di as error "dea_reference: policymatrix() must be 1 x N or N x N for the marked technology sample"
            exit 503
        }
    }

    mata: _dea_reference_solve("`inputs'", "`outputs'", "`touse'", ///
        "`evaluse'", "`canduse'", ///
        "`RAD'", "`TS'", "`SW'", "`BCCPX'", "`BCCPY'", ///
        "`rts'", "`orientation'", "`engine'", ///
        "`scaling'", `tolerance', `peertolerance', `opttolerance', ///
        "`SEL'", "`lexicographic'", "`lexorder'", "`sparse'", "`broad'", ///
        "`relativeinterior'", "`basic'", "`policyweight'", "`policymatrix'", ///
        "`policysense'", "`audit'", "`rtsonly'", "`interceptbounds'", ///
        `maxnodes', `inttol', "`trace'")

    capture matrix drop __dea_ref_milpc __dea_ref_milpaeq ///
        __dea_ref_milpbeq __dea_ref_milpaub __dea_ref_milpbub ///
        __dea_ref_milplb __dea_ref_milpub __dea_ref_milpints

    local failures = r(aux_failures)
    local possible_count = r(n_possible)
    local necessary_count = r(n_necessary)
    local optional_count = r(n_optional)
    local impossible_count = r(n_impossible)
    local opttolret = r(opt_tolerance)
    local lexfailures = r(lex_failures)
    local sparsefailures = r(sparse_failures)
    local sparsenodes = r(sparse_nodes)
    local broadfailures = r(broad_failures)
    local relativefailures = r(relative_failures)
    local basicfailures = r(basic_failures)
    local policyfailures = r(policy_failures)
    local auxsolves = r(aux_lp_solves)
    local sparselps = r(sparse_lp_solves)
    local totalsolves = `auxsolves' + `sparselps'
    tempname LMIN LMAX POSS NEC CLASS ASTAT EIDX CIDX SSEL LEX LSP SSP SSTAT LBR BSP BSTAT LRI SRI RSTAT RWIT LBA SBA BASTAT LPO SPO POSTAT POBJ ASUP AL1 AMET DASUP DAL1 DAMET RTSB RTSC RTSSTAT RTSM RTSS RTSCMP DRTSB DRTSC DRTSCMP NDIS U0B U0C U0S U0A DU0B DU0A
    matrix `LMIN' = r(lambda_min)
    matrix `LMAX' = r(lambda_max)
    matrix `RTSB' = r(rts_sum_bounds)
    matrix `RTSC' = r(rts_robust_code)
    matrix `RTSSTAT' = r(rts_bound_status)
    matrix `RTSM' = r(rts_margin)
    matrix `RTSS' = r(rts_stable)
    matrix `U0B' = r(u0_bounds)
    matrix `U0C' = r(u0_rts_code)
    matrix `U0S' = r(u0_bound_status)
    matrix `U0A' = r(u0_rts_agree)
    matrix colnames `U0B' = u0_min u0_max
    matrix colnames `U0C' = u0_rts
    matrix colnames `U0S' = status
    matrix colnames `U0A' = agree
    matrix colnames `RTSB' = sum_lambda_min sum_lambda_max
    matrix colnames `RTSC' = robust_rts
    matrix colnames `RTSSTAT' = status
    matrix colnames `RTSM' = margin
    matrix colnames `RTSS' = stable
    matrix `POSS' = r(possible)
    matrix `NEC' = r(necessary)
    matrix `CLASS' = r(reference_status)
    matrix `ASTAT' = r(aux_status)
    matrix `EIDX' = r(evaluated_index)
    matrix `CIDX' = r(candidate_index)
    matrix `RTSCMP' = J(`N_eval',3,.)
    if "`compare'" != "" {
        mata: st_matrix("`RTSCMP'", _dea_rts_compare( ///
            st_matrix("`RAD'"),st_matrix("`VRSRAD'"),st_matrix("`NIRSRAD'"), ///
            st_matrix("`EIDX'"),st_matrix("`RTSC'"),`tolerance'))
    }
    matrix colnames `RTSCMP' = standard robust agree
    mata: st_numscalar("`NDIS'",sum(st_matrix("`RTSCMP'")[,3]:==0))
    matrix `SSEL' = r(lambda_selected_subset)
    matrix `LEX' = r(lambda_lexicographic)
    matrix `LSP' = r(lambda_sparse)
    matrix `SSP' = r(sparse_support)
    matrix `SSTAT' = r(sparse_status)
    matrix `LBR' = r(lambda_broad)
    matrix `BSP' = r(broad_support)
    matrix `BSTAT' = r(broad_status)
    matrix `LRI' = r(lambda_relative)
    matrix `SRI' = r(relative_solution)
    matrix `RSTAT' = r(relative_status)
    matrix `RWIT' = r(relative_witnesses)
    matrix `LBA' = r(lambda_basic)
    matrix `SBA' = r(basic_solution)
    matrix `BASTAT' = r(basic_status)
    matrix `LPO' = r(lambda_policy)
    matrix `SPO' = r(policy_solution)
    matrix `POSTAT' = r(policy_status)
    matrix `POBJ' = r(policy_objective)
    matrix `ASUP' = r(audit_support)
    matrix `AL1' = r(audit_l1_selected)
    matrix `AMET' = r(audit_identification)
    matrix colnames `ASUP' = selected basic sparse broad relative policy
    matrix colnames `AL1' = basic sparse broad relative policy
    matrix colnames `AMET' = possible missing_selected mean_range max_range
    matrix `DASUP' = `ASUP'
    matrix `DAL1' = `AL1'
    matrix `DAMET' = `AMET'
    matrix `DRTSB' = `RTSB'
    matrix `DRTSC' = `RTSC'
    matrix `DRTSCMP' = `RTSCMP'
    matrix `DU0B' = `U0B'
    matrix `DU0A' = `U0A'
    return clear

    if "`generate'" != "" {
        local gen_suffixes "lb ub code status margin stable"
        if "`compare'" != "" local gen_suffixes "`gen_suffixes' standard robust agree"
        foreach suffix of local gen_suffixes {
            capture confirm new variable `generate'_`suffix'
            if _rc {
                di as error "dea_reference: variable `generate'_`suffix' already exists"
                exit 110
            }
        }
        tempvar __rts_eval_order
        quietly generate long `__rts_eval_order' = sum(`evaluse')
        quietly generate double `generate'_lb = .
        quietly generate double `generate'_ub = .
        quietly generate byte `generate'_code = .
        quietly generate int `generate'_status = .
        quietly generate double `generate'_margin = .
        quietly generate byte `generate'_stable = .
        if "`compare'" != "" {
            quietly generate byte `generate'_standard = .
            quietly generate byte `generate'_robust = .
            quietly generate byte `generate'_agree = .
        }
        quietly replace `generate'_lb = el(`RTSB',`__rts_eval_order',1) if `evaluse'
        quietly replace `generate'_ub = el(`RTSB',`__rts_eval_order',2) if `evaluse'
        quietly replace `generate'_code = el(`RTSC',`__rts_eval_order',1) if `evaluse'
        quietly replace `generate'_status = el(`RTSSTAT',`__rts_eval_order',1) if `evaluse'
        quietly replace `generate'_margin = el(`RTSM',`__rts_eval_order',1) if `evaluse'
        quietly replace `generate'_stable = el(`RTSS',`__rts_eval_order',1) if `evaluse'
        if "`compare'" != "" {
            quietly replace `generate'_standard = el(`RTSCMP',`__rts_eval_order',1) if `evaluse'
            quietly replace `generate'_robust = el(`RTSCMP',`__rts_eval_order',2) if `evaluse'
            quietly replace `generate'_agree = el(`RTSCMP',`__rts_eval_order',3) if `evaluse'
        }
        label variable `generate'_lb "Minimum sum(lambda) on CCR optimal set"
        label variable `generate'_ub "Maximum sum(lambda) on CCR optimal set"
        label variable `generate'_code "Robust RTS code (-1 DRS, 0 CRS, 1 IRS)"
        label variable `generate'_status "Robust RTS bound status (0 success)"
        label variable `generate'_margin "RTS interval separation from classification boundary"
        label variable `generate'_stable "RTS margin exceeds local optimality tolerance"
        if "`compare'" != "" {
            label variable `generate'_standard "Standard CRS-VRS-NIRS RTS code"
            label variable `generate'_robust "Optimal-set robust RTS code"
            label variable `generate'_agree "Standard and robust RTS codes agree"
        }
    }

    if "`saving'" != "" {
        if "`dmu'" != "" {
            tempfile dmu_map peer_map
            preserve
            quietly keep if `touse'
            quietly keep `dmu'
            quietly generate long dmu_index = _n
            quietly rename `dmu' dmu_id
            quietly save "`dmu_map'"
            quietly rename dmu_index peer_index
            quietly rename dmu_id peer_id
            quietly save "`peer_map'"
            restore
        }
        preserve
        quietly clear
        quietly set obs `=`N_eval'*`N_cand''
        quietly generate long dmu_index = .
        quietly generate long peer_index = .
        quietly generate double lambda_min = .
        quietly generate double lambda_max = .
        quietly generate double lambda_selected = .
        quietly generate byte reference_status = .
        quietly generate double lambda_lexicographic = .
        quietly generate double lambda_sparse = .
        quietly generate byte sparse_selected = .
        quietly generate double lambda_broad = .
        quietly generate byte broad_selected = .
        quietly generate double lambda_relative = .
        quietly generate double lambda_basic = .
        quietly generate double lambda_policy = .
        quietly generate str10 status = ""
        quietly replace dmu_index = el(`EIDX',ceil(_n/`N_cand'),1)
        quietly replace peer_index = el(`CIDX',mod(_n-1,`N_cand')+1,1)
        mata: st_store(., "lambda_min", vec(st_matrix("`LMIN'")') )
        mata: st_store(., "lambda_max", vec(st_matrix("`LMAX'")') )
        mata: st_store(., "lambda_selected", ///
            vec(st_matrix("`SEL'")[st_matrix("`EIDX'"),st_matrix("`CIDX'")]') )
        mata: st_store(., "reference_status", vec(st_matrix("`CLASS'")') )
        if "`lexicographic'" != "" ///
            mata: st_store(., "lambda_lexicographic", ///
                vec(st_matrix("`LEX'")[.,st_matrix("`CIDX'")]') )
        if "`sparse'" != "" {
            mata: st_store(., "lambda_sparse", ///
                vec(st_matrix("`LSP'")[.,st_matrix("`CIDX'")]') )
            quietly replace sparse_selected = lambda_sparse>`peertolerance' ///
                if !missing(lambda_sparse)
        }
        if "`broad'" != "" {
            mata: st_store(., "lambda_broad", ///
                vec(st_matrix("`LBR'")[.,st_matrix("`CIDX'")]') )
            quietly replace broad_selected = reference_status>0 ///
                if !missing(lambda_broad) & !missing(reference_status)
        }
        if "`relativeinterior'" != "" ///
            mata: st_store(., "lambda_relative", ///
                vec(st_matrix("`LRI'")[.,st_matrix("`CIDX'")]') )
        if "`basic'" != "" ///
            mata: st_store(., "lambda_basic", ///
                vec(st_matrix("`LBA'")[.,st_matrix("`CIDX'")]') )
        if "`policyweight'`policymatrix'" != "" ///
            mata: st_store(., "lambda_policy", ///
                vec(st_matrix("`LPO'")[.,st_matrix("`CIDX'")]') )
        quietly replace status = "impossible" if reference_status==0
        quietly replace status = "optional" if reference_status==1
        quietly replace status = "necessary" if reference_status==2
        quietly replace status = "failed" if missing(reference_status)
        capture label drop dea_ref_status
        label define dea_ref_status 0 "impossible" 1 "optional" 2 "necessary"
        label values reference_status dea_ref_status
        if "`dmu'" != "" {
            quietly merge m:1 dmu_index using "`dmu_map'", keep(master match) nogen
            assert !missing(dmu_id)
            quietly merge m:1 peer_index using "`peer_map'", keep(master match) nogen
            assert !missing(peer_id)
            order dmu_index dmu_id peer_index peer_id
        }
        if "`replace'" != "" save "`saving'", replace
        else save "`saving'"
        restore
    }

    return scalar N = `N'
    return scalar N_evaluated = `N_eval'
    return scalar N_candidates = `N_cand'
    return scalar aux_failures = `failures'
    return scalar lex_failures = `lexfailures'
    return scalar sparse_failures = `sparsefailures'
    return scalar sparse_nodes = `sparsenodes'
    return scalar broad_failures = `broadfailures'
    return scalar relative_failures = `relativefailures'
    return scalar basic_failures = `basicfailures'
    return scalar policy_failures = `policyfailures'
    return scalar lp_solves_planned = `planned'
    return scalar aux_lp_solves = `auxsolves'
    return scalar sparse_lp_solves = `sparselps'
    return scalar total_lp_solves = `totalsolves'
    return scalar rts_compare_lp_planned = `compare_lp_planned'
    return scalar maxauxlp = `maxauxlp'
    return scalar n_possible = `possible_count'
    return scalar n_necessary = `necessary_count'
    return scalar n_optional = `optional_count'
    return scalar n_impossible = `impossible_count'
    return scalar tolerance = `tolerance'
    return scalar peer_tolerance = `peertolerance'
    return scalar opt_tolerance = `opttolret'
    return scalar radial_tolerance = `radialtolerance'
    return matrix lambda_selected = `SEL'
    return matrix lambda_selected_subset = `SSEL'
    return matrix lambda_lexicographic = `LEX'
    return matrix lambda_sparse = `LSP'
    return matrix sparse_support = `SSP'
    return matrix sparse_status = `SSTAT'
    return matrix lambda_broad = `LBR'
    return matrix broad_support = `BSP'
    return matrix broad_status = `BSTAT'
    return matrix lambda_relative = `LRI'
    return matrix relative_solution = `SRI'
    return matrix relative_status = `RSTAT'
    return matrix relative_witnesses = `RWIT'
    return matrix lambda_basic = `LBA'
    return matrix basic_solution = `SBA'
    return matrix basic_status = `BASTAT'
    return matrix lambda_policy = `LPO'
    return matrix policy_solution = `SPO'
    return matrix policy_status = `POSTAT'
    return matrix policy_objective = `POBJ'
    return matrix audit_support = `ASUP'
    return matrix audit_l1_selected = `AL1'
    return matrix audit_identification = `AMET'
    return matrix lambda_min = `LMIN'
    return matrix lambda_max = `LMAX'
    return matrix rts_sum_bounds = `RTSB'
    return matrix rts_robust_code = `RTSC'
    return matrix rts_bound_status = `RTSSTAT'
    return matrix rts_margin = `RTSM'
    return matrix rts_stable = `RTSS'
    return matrix rts_compare = `RTSCMP'
    return matrix u0_bounds = `U0B'
    return matrix u0_rts_code = `U0C'
    return matrix u0_bound_status = `U0S'
    return matrix u0_rts_agree = `U0A'
    return matrix radial_adjustment = `RADADJ'
    if "`compare'" != "" {
        return scalar n_rts_disagree = `NDIS'
    }
    return matrix possible = `POSS'
    return matrix necessary = `NEC'
    return matrix reference_status = `CLASS'
    return matrix aux_status = `ASTAT'
    return matrix evaluated_index = `EIDX'
    return matrix candidate_index = `CIDX'
    return local cmd "dea_reference"
    return local model "te"
    return local stage "2"
    return local inputs "`inputs'"
    return local outputs "`outputs'"
    return local rts "`rts'"
    return local orientation "`orientation'"
    return local engine_requested "`engine'"
    return local scaling "`scaling'"
    return local slack_normalization "`slacknormalization'"
    return local dmu "`dmu'"
    return local evaluate "`evaluate'"
    return local candidates "`candidates'"
    return local lexicographic "`lexicographic'"
    return local lexorder "`lexorder'"
    return local sparse "`sparse'"
    return local broad "`broad'"
    return local relativeinterior "`relativeinterior'"
    return local basic "`basic'"
    return local policyweight "`policyweight'"
    return local policymatrix "`policymatrix'"
    return local policysense "`policysense'"
    return local audit "`audit'"
    return local rtsonly "`rtsonly'"
    return local generate "`generate'"
    return local compare "`compare'"
    return local interceptbounds "`interceptbounds'"
    return local saving "`saving'"

    di as text "Two-stage optimal-set reference diagnostics"
    di as text "  technology DMUs: " as result %9.0f `N' ///
        as text "   evaluated: " as result %9.0f `N_eval' ///
        as text "   candidates: " as result %9.0f `N_cand'
    di as text "  " ///
        as text "   auxiliary failures: " as result %9.0f `failures'
    if "`rtsonly'" == "" di as text "  possible: " as result %9.0f `possible_count' ///
        as text "   necessary: " as result %9.0f `necessary_count' ///
        as text "   optional: " as result %9.0f `optional_count'
    if "`lexicographic'" != "" di as text "  lexicographic(`lexicographic') failures: " ///
        as result %9.0f `lexfailures'
    if "`sparse'" != "" di as text "  sparse-selection failures: " ///
        as result %9.0f `sparsefailures' as text "   B&B nodes: " ///
        as result %9.0f `sparsenodes'
    if "`broad'" != "" di as text "  broad-support failures: " ///
        as result %9.0f `broadfailures'
    if "`relativeinterior'" != "" di as text "  relative-interior failures: " ///
        as result %9.0f `relativefailures'
    if "`basic'" != "" di as text "  simplex-basic failures: " ///
        as result %9.0f `basicfailures'
    if "`policyweight'`policymatrix'" != "" di as text "  policy-selection failures: " ///
        as result %9.0f `policyfailures'
    if "`audit'" != "" {
        di as text "  benchmark audit support counts:"
        matrix list `DASUP', format(%9.3g)
        di as text "  L1 distance from selected decomposition:"
        matrix list `DAL1', format(%9.3g)
        di as text "  optimal-set identification diagnostics:"
        matrix list `DAMET', format(%9.3g)
    }
    if "`rts'" == "crs" {
        di as text "  robust RTS bounds (min/max sum(lambda)):"
        matrix list `DRTSB', format(%10.6f)
        di as text "  RTS codes: -1 = DRS, 0 = CRS, 1 = IRS"
        matrix list `DRTSC', format(%4.0f)
    }
    if "`compare'" != "" {
        di as text "  standard/robust RTS comparison:"
        matrix list `DRTSCMP', format(%4.0f)
        di as text "  disagreements: " as result %9.0f `NDIS'
        di as text "  standard VRS/NIRS comparison LPs planned: " ///
            as result %9.0f `compare_lp_planned'
    }
    if "`interceptbounds'" != "" {
        di as text "  BCC intercept bounds and agreement with CCR interval:"
        matrix list `DU0B', format(%10.6f)
        matrix list `DU0A', format(%4.0f)
    }
    di as text "  LP calls planned/continuous/sparse/total: " ///
        as result %9.0f `planned' " / " %9.0f `auxsolves' " / " ///
        %9.0f `sparselps' " / " %9.0f `totalsolves'
end

capture mata: mata drop _dea_reference_solve()
capture mata: mata drop _dea_rts_compare()
mata:
mata set matastrict on

real matrix _dea_rts_compare(real matrix crsin, real matrix vrsin,
    real matrix nirsin, real matrix eiin, real matrix robustin,
    real scalar tol)
{
    real colvector crs, vrs, nirs, ei, robust, standard, agree
    real scalar a, k
    crs=vec(crsin); vrs=vec(vrsin); nirs=vec(nirsin)
    ei=vec(eiin); robust=vec(robustin)
    standard=J(rows(ei),1,.)
    for (a=1; a<=rows(ei); a++) {
        k=ei[a]
        if (crs[k]<. & vrs[k]<. & abs(crs[k]-vrs[k])<tol) standard[a]=0
        else if (nirs[k]<. & vrs[k]<. & abs(nirs[k]-vrs[k])<tol) standard[a]=-1
        else if (vrs[k]<.) standard[a]=1
    }
    agree=(standard:==robust)
    return((standard,robust,agree))
}

void _dea_reference_solve(string scalar inputs, string scalar outputs,
    string scalar touse, string scalar evaluse, string scalar canduse,
    string scalar radialname, string scalar totalslackname, string scalar slackweightname,
    string scalar projectedxname, string scalar projectedyname,
    string scalar rts, string scalar orientation, string scalar engine,
    string scalar scaling, real scalar tol, real scalar peertol,
    real scalar useropttol, string scalar selectedname,
    string scalar lexsense, string scalar lexorder, string scalar sparse,
    string scalar broad, string scalar relativeinterior, string scalar basic,
    string scalar policyweight, string scalar policymatrix, string scalar policysense,
    string scalar audit, string scalar rtsonly, string scalar interceptbounds,
    real scalar maxnodes, real scalar inttol, string scalar trace)
{
    real matrix X, Y, Aeq, Aub, Lmin, Lmax, possible, necessary, classes, auxstat, Llex
    real matrix RtsBounds, ProjectedX, ProjectedY, Aeq0, Aub0, U0Bounds
    real matrix Lsparse, Lbroad, Lrelative, Srelative, Lbasic, Sbasic, Lpolicy, Spolicy
    real matrix selected, milpAeq, milpAub, milpx, AuditSupport, AuditL1, AuditIdentification
    real colvector SparseStatus, BroadStatus, RelativeStatus, RelativeWitnesses
    real colvector BasicStatus, PolicyStatus, PolicyObjective, RtsCode, RtsStatus, RtsMargin, RtsStable
    real colvector U0Code, U0Status, U0Agree, beq0, bub0
    real matrix policyweights
    real colvector radial, totalslack, x0, y0, beq, bub, milpbub, ei, ci, lexpriority, lexidx
    real rowvector c, lb, ub, bigm, milpc, milplb, milpub, ints, broadacc, riacc, riub
    real rowvector slackweights, c0, lb0, ub0, objrow
    real scalar n, ne, nc, m, s, nv, nv0, a, b, p, k, j, fallback, objscale, rowscalemin, rowscalemax
    real scalar opttol, localopttol, fixtol, failures, lexfailures, np, nn, no, ni
    real scalar sparsefailures, sparsenodes, sparselps, sparsestatus, auxsolves
    real scalar broadfailures, broadstatus, broadcount
    real scalar relativefailures, relativestatus, relativecount, relativeslack
    real scalar basicfailures
    real scalar policyfailures
    struct _DEA_LPResult scalar outmin, outmax, outlex, outbase, outri, outbasic, outpolicy
    struct _DEA_LPResult scalar outu0min, outu0max

    X=st_data(.,tokens(inputs),touse); Y=st_data(.,tokens(outputs),touse)
    ei=selectindex(st_data(.,evaluse,touse)); ci=selectindex(st_data(.,canduse,touse))
    radial=st_matrix(radialname); totalslack=st_matrix(totalslackname)
    slackweights=st_matrix(slackweightname)
    ProjectedX=st_matrix(projectedxname); ProjectedY=st_matrix(projectedyname)
    selected=st_matrix(selectedname)
    n=rows(X); ne=rows(ei); nc=rows(ci); m=cols(X); s=cols(Y); nv=n+m+s
    Lmin=J(ne,nc,.); Lmax=J(ne,nc,.); auxstat=J(ne,nc,.); Llex=J(ne,n,.)
    RtsBounds=J(ne,2,.); RtsCode=J(ne,1,.); RtsStatus=J(ne,1,.)
    RtsMargin=J(ne,1,.); RtsStable=J(ne,1,.)
    U0Bounds=J(ne,2,.); U0Code=J(ne,1,.); U0Status=J(ne,1,.)
    U0Agree=J(ne,1,.)
    Lsparse=J(ne,n,.); SparseStatus=J(ne,1,.); sparsefailures=0; sparsenodes=0; sparselps=0
    Lbroad=J(ne,n,.); BroadStatus=J(ne,1,.); broadfailures=0
    Lrelative=J(ne,n,.); Srelative=J(ne,nv,.)
    RelativeStatus=J(ne,1,.); RelativeWitnesses=J(ne,1,.)
    relativefailures=0
    Lbasic=J(ne,n,.); Sbasic=J(ne,nv,.); BasicStatus=J(ne,1,.); basicfailures=0
    Lpolicy=J(ne,n,.); Spolicy=J(ne,nv,.); PolicyStatus=J(ne,1,.)
    PolicyObjective=J(ne,1,.); policyfailures=0
    if (policyweight!="") policyweights=st_data(.,policyweight,touse)
    else if (policymatrix!="") policyweights=st_matrix(policymatrix)
    failures=0; lexfailures=0; auxsolves=0
    opttol=(useropttol>0 ? useropttol : max((100*tol,peertol)))
    if (lexsense!="") {
        lexpriority=(lexorder=="" ? (1::n) : st_data(.,lexorder,touse))
        lexidx=order((lexpriority,(1::n)),(1,2))
    }

    for (a=1; a<=ne; a++) {
        broadacc=J(1,nv,0); broadcount=0; broadstatus=0
        k=ei[a]
        x0=X[k,.]'; y0=Y[k,.]'
        Aeq=(X',I(m),J(m,s,0)) \ (Y',J(s,m,0),-I(s))
        beq=(orientation=="input" ? radial[k]*x0 : x0) \
             (orientation=="input" ? y0 : radial[k]*y0)
        Aub=J(0,nv,0); bub=J(0,1,0)
        if (rts=="vrs") {
            Aeq=Aeq \ (J(1,n,1),J(1,m+s,0)); beq=beq \ 1
        }
        else if (rts=="drs") {
            Aub=(J(1,n,1),J(1,m+s,0)); bub=1
        }
        else if (rts=="irs") {
            Aub=(J(1,n,-1),J(1,m+s,0)); bub=-1
        }
        localopttol=opttol*(1+abs(totalslack[k]))
        Aub=Aub \ (J(1,n,0),slackweights) \
                   (J(1,n,0),-slackweights)
        bub=bub \ (totalslack[k]+localopttol) \
                 (-(totalslack[k]-localopttol))
        lb=J(1,nv,0); ub=J(1,nv,.)
        if (rts=="crs") {
            c=(J(1,n,1),J(1,m+s,0))
            outmin=_dea_lpcontract_dispatch(c,"min",Aeq,beq,Aub,bub,lb,ub,
                engine,tol,scaling,0,trace,fallback,objscale,
                rowscalemin,rowscalemax)
            outmax=_dea_lpcontract_dispatch(c,"max",Aeq,beq,Aub,bub,lb,ub,
                engine,tol,scaling,0,trace,fallback,objscale,
                rowscalemin,rowscalemax)
            auxsolves=auxsolves+2
            if (outmin.status==0 & outmax.status==0) {
                RtsBounds[a,.]=(max((0,outmin.objective)), max((0,outmax.objective)))
                RtsStatus[a]=0
                if (RtsBounds[a,1]<=1+localopttol & RtsBounds[a,2]>=1-localopttol) RtsCode[a]=0
                else if (RtsBounds[a,2]<1-localopttol) RtsCode[a]=1
                else if (RtsBounds[a,1]>1+localopttol) RtsCode[a]=-1
                else RtsStatus[a]=6
                if (RtsStatus[a]==0) {
                    if (RtsCode[a]==1) RtsMargin[a]=1-RtsBounds[a,2]
                    else if (RtsCode[a]==-1) RtsMargin[a]=RtsBounds[a,1]-1
                    else RtsMargin[a]=min((1-RtsBounds[a,1],RtsBounds[a,2]-1))
                    RtsStable[a]=(RtsMargin[a]>localopttol)
                }
            }
            else {
                failures++
                RtsStatus[a]=100*outmin.status+outmax.status
            }
        }
        if (interceptbounds!="") {
            nv0=m+s+2
            c0=(J(1,m+s,0),1,-1)
            lb0=J(1,nv0,0); ub0=J(1,nv0,.)
            if (orientation=="input") {
                Aeq0=(ProjectedX[k,.],J(1,s+2,0)); beq0=1
                Aub0=(-X,Y,J(n,1,-1),J(n,1,1))
                objrow=(J(1,m,0),ProjectedY[k,.],-1,1)
            }
            else {
                Aeq0=(J(1,m,0),ProjectedY[k,.],0,0); beq0=1
                Aub0=(-X,Y,J(n,1,-1),J(n,1,1))
                objrow=(ProjectedX[k,.],J(1,s,0),1,-1)
            }
            Aub0=Aub0\objrow\(-objrow)
            bub0=J(n,1,0)\(1+opttol)\(-(1-opttol))
            outu0min=_dea_lpcontract_dispatch(c0,"min",Aeq0,beq0,Aub0,bub0,
                lb0,ub0,engine,tol,scaling,0,trace,fallback,objscale,
                rowscalemin,rowscalemax)
            outu0max=_dea_lpcontract_dispatch(c0,"max",Aeq0,beq0,Aub0,bub0,
                lb0,ub0,engine,tol,scaling,0,trace,fallback,objscale,
                rowscalemin,rowscalemax)
            auxsolves=auxsolves+2
            if (outu0min.status!=0 & engine=="ipm") {
                outu0min=_dea_lpcontract_dispatch(c0,"min",Aeq0,beq0,Aub0,bub0,
                    lb0,ub0,"simplex",tol,scaling,0,trace,fallback,objscale,
                    rowscalemin,rowscalemax)
                auxsolves++
            }
            if (outu0max.status!=0 & engine=="ipm") {
                outu0max=_dea_lpcontract_dispatch(c0,"max",Aeq0,beq0,Aub0,bub0,
                    lb0,ub0,"simplex",tol,scaling,0,trace,fallback,objscale,
                    rowscalemin,rowscalemax)
                auxsolves++
            }
            if (outu0min.status==0 & outu0max.status==0) {
                U0Bounds[a,.]=(outu0min.objective,outu0max.objective)
                U0Status[a]=0
                if (U0Bounds[a,1]<=opttol & U0Bounds[a,2]>=-opttol) U0Code[a]=0
                else if (U0Bounds[a,2]<-opttol) U0Code[a]=1
                else if (U0Bounds[a,1]>opttol) U0Code[a]=-1
                else U0Status[a]=6
                if (U0Status[a]==0 & RtsStatus[a]==0) U0Agree[a]=(U0Code[a]==RtsCode[a])
            }
            else if (outu0min.status==0 & outu0max.status==2) {
                U0Bounds[a,1]=outu0min.objective
                U0Status[a]=0
                U0Code[a]=(outu0min.objective>opttol ? -1 : 0)
                if (RtsStatus[a]==0) U0Agree[a]=(U0Code[a]==RtsCode[a])
            }
            else if (outu0min.status==2 & outu0max.status==0) {
                U0Bounds[a,2]=outu0max.objective
                U0Status[a]=0
                U0Code[a]=(outu0max.objective < -opttol ? 1 : 0)
                if (RtsStatus[a]==0) U0Agree[a]=(U0Code[a]==RtsCode[a])
            }
            else if (outu0min.status==2 & outu0max.status==2) {
                U0Status[a]=0; U0Code[a]=0
                if (RtsStatus[a]==0) U0Agree[a]=(RtsCode[a]==0)
            }
            else {
                failures++
                U0Status[a]=100*outu0min.status+outu0max.status
            }
        }
        if (rtsonly=="") for (b=1; b<=nc; b++) {
            j=ci[b]
            c=J(1,nv,0); c[j]=1
            outmin=_dea_lpcontract_dispatch(c,"min",Aeq,beq,Aub,bub,lb,ub,
                engine,tol,scaling,0,trace,fallback,objscale,
                rowscalemin,rowscalemax)
            outmax=_dea_lpcontract_dispatch(c,"max",Aeq,beq,Aub,bub,lb,ub,
                engine,tol,scaling,0,trace,fallback,objscale,
                rowscalemin,rowscalemax)
            auxsolves=auxsolves+2
            if (outmin.status==0 & outmax.status==0) {
                Lmin[a,b]=max((0,outmin.objective))
                Lmax[a,b]=max((0,outmax.objective))
                auxstat[a,b]=0
                if (broad!="" & outmax.objective>peertol) {
                    broadacc=broadacc+outmax.x
                    broadcount++
                }
            }
            else {
                failures++
                auxstat[a,b]=100*outmin.status+outmax.status
                if (broad!="") broadstatus=5
            }
        }
        if (broad!="") {
            if (broadstatus==0 & broadcount>0) Lbroad[a,.]=broadacc[1..n]/broadcount
            else {
                if (broadstatus==0) broadstatus=1
                broadfailures++
            }
            BroadStatus[a]=broadstatus
        }
        if (relativeinterior!="") {
            relativestatus=0; relativecount=0; riacc=J(1,nv,0)
            c=J(1,nv,0)
            outbase=_dea_lpcontract_dispatch(c,"min",Aeq,beq,Aub,bub,lb,ub,
                engine,tol,scaling,0,trace,fallback,objscale,
                rowscalemin,rowscalemax)
            auxsolves++
            if (outbase.status!=0) relativestatus=5
            else {
                riacc=outbase.x; relativecount=1
                /* A bounded box around a feasible base point does not remove
                   any relative-interior direction: any witness can be moved
                   sufficiently close to the base by convex combination. */
                riub=colmax((J(1,nv,1)\(outbase.x:+1)))
                for (j=1; j<=nv; j++) {
                    c=J(1,nv,0); c[j]=1
                    outri=_dea_lpcontract_dispatch(c,"max",Aeq,beq,Aub,bub,lb,riub,
                        engine,tol,scaling,0,trace,fallback,objscale,
                        rowscalemin,rowscalemax)
                    auxsolves++
                    if (outri.status!=0) {
                        relativestatus=5
                        break
                    }
                    if (outri.objective>tol) {
                        riacc=riacc+outri.x
                        relativecount++
                    }
                }
                if (relativestatus==0) {
                    for (p=1; p<=rows(Aub); p++) {
                        c=Aub[p,.]
                        outri=_dea_lpcontract_dispatch(c,"min",Aeq,beq,Aub,bub,lb,riub,
                            engine,tol,scaling,0,trace,fallback,objscale,
                            rowscalemin,rowscalemax)
                        auxsolves++
                        if (outri.status!=0) {
                            relativestatus=5
                            break
                        }
                        relativeslack=bub[p]-outri.objective
                        if (relativeslack>tol) {
                            riacc=riacc+outri.x
                            relativecount++
                        }
                    }
                }
            }
            if (relativestatus==0 & relativecount>0) {
                Srelative[a,.]=riacc/relativecount
                Lrelative[a,.]=Srelative[a,1..n]
            }
            else relativefailures++
            RelativeStatus[a]=relativestatus
            RelativeWitnesses[a]=relativecount
        }
        if (basic!="") {
            /* Distinct order weights define a transparent third objective.
               The released revised-simplex adapter returns a basic feasible
               optimum; this is basic recovery, not IPM basis crossover. */
            c=((n..1)/n,J(1,m+s,0))
            outbasic=_dea_lpcontract_dispatch(c,"max",Aeq,beq,Aub,bub,lb,ub,
                "simplex",tol,scaling,0,trace,fallback,objscale,
                rowscalemin,rowscalemax)
            auxsolves++
            BasicStatus[a]=outbasic.status
            if (outbasic.status==0) {
                Sbasic[a,.]=outbasic.x
                Lbasic[a,.]=outbasic.x[1..n]
            }
            else basicfailures++
        }
        if (policyweight!="" | policymatrix!="") {
            c=((policyweight!="" ? policyweights' :
                (rows(policyweights)==1 ? policyweights : policyweights[k,.])),
                J(1,m+s,0))
            outpolicy=_dea_lpcontract_dispatch(c,policysense,Aeq,beq,Aub,bub,lb,ub,
                engine,tol,scaling,0,trace,fallback,objscale,
                rowscalemin,rowscalemax)
            auxsolves++
            PolicyStatus[a]=outpolicy.status
            if (outpolicy.status==0) {
                Spolicy[a,.]=outpolicy.x
                Lpolicy[a,.]=outpolicy.x[1..n]
                PolicyObjective[a]=outpolicy.objective
            }
            else policyfailures++
        }
        if (sparse!="") {
            bigm=(rts=="vrs" | rts=="drs" ? J(1,n,1) : Lmax[a,.])
            sparsestatus=0
            if (cols(bigm)!=n | any(bigm:>=.)) sparsestatus=5
            if (sparsestatus==0) {
                bigm=colmax((J(1,n,0)\bigm))
                milpc=(J(1,nv,0),J(1,n,1))
                milpAeq=(Aeq,J(rows(Aeq),n,0))
                milpAub=(Aub,J(rows(Aub),n,0)) \
                    (I(n),J(n,m+s,0),-diag(bigm))
                milpbub=bub\J(n,1,0)
                milplb=J(1,nv+n,0); milpub=(J(1,nv,.),J(1,n,1))
                ints=(nv+1)..(nv+n)
                st_matrix("__dea_ref_milpc",milpc)
                st_matrix("__dea_ref_milpaeq",milpAeq)
                st_matrix("__dea_ref_milpbeq",beq)
                st_matrix("__dea_ref_milpaub",milpAub)
                st_matrix("__dea_ref_milpbub",milpbub)
                st_matrix("__dea_ref_milplb",milplb)
                st_matrix("__dea_ref_milpub",milpub)
                st_matrix("__dea_ref_milpints",ints)
                if (trace!="") printf("{txt}sparse MILP DMU %g: variables=%g equalities=%g inequalities=%g\n",
                    k,cols(milpc),rows(milpAeq),rows(milpAub))
                milp_bb_solve("__dea_ref_milpc","min","__dea_ref_milpaeq",
                    "__dea_ref_milpbeq","__dea_ref_milpaub","__dea_ref_milpbub",
                    "__dea_ref_milplb","__dea_ref_milpub","__dea_ref_milpints",
                    maxnodes,inttol,tol,engine,scaling,"",trace)
                sparsestatus=st_numscalar("r(status)")
                sparsenodes=sparsenodes+st_numscalar("r(nodes)")
                sparselps=sparselps+st_numscalar("r(lp_solves)")
                if (sparsestatus==0) {
                    milpx=st_matrix("r(x)")
                    Lsparse[a,.]=milpx[1,1..n]
                }
            }
            if (sparsestatus!=0) sparsefailures++
            SparseStatus[a]=sparsestatus
        }
        if (lexsense!="") {
            for (p=1; p<=n; p++) {
                j=lexidx[p]
                c=J(1,nv,0); c[j]=1
                outlex=_dea_lpcontract_dispatch(c,lexsense,Aeq,beq,Aub,bub,lb,ub,
                    engine,tol,scaling,0,trace,fallback,objscale,
                    rowscalemin,rowscalemax)
                auxsolves++
                /* Every lexicographic subproblem contains the preceding
                   accepted solution within its fixing band. Under auto,
                   verify any simplex terminal classification with IPM. */
                if (engine=="auto" & outlex.status!=0) {
                    outlex=_dea_lpcontract_dispatch(c,lexsense,Aeq,beq,Aub,bub,lb,ub,
                        "ipm",tol,scaling,0,trace,fallback,objscale,
                        rowscalemin,rowscalemax)
                    auxsolves++
                }
                if (outlex.status!=0) {
                    lexfailures++
                    break
                }
                /* Sequential fixing accumulates primal residuals across up to
                   n solves. Use a hundredfold band relative to the optimal-set
                   base tolerance; classification still uses peertol. */
                fixtol=100*opttol*(1+abs(outlex.objective))
                /* A one-sided fixing row is sufficient: later problems
                   already remain on the same feasible optimal set. It also
                   keeps the accepted primal point feasible despite small
                   objective/parameter reporting differences. */
                if (lexsense=="max") {
                    Aub=Aub \ (-I(nv)[j,.])
                    bub=bub \ (-(outlex.x[j]-fixtol))
                }
                else {
                    Aub=Aub \ (I(nv)[j,.])
                    bub=bub \ (outlex.x[j]+fixtol)
                }
            }
            if (p>n & outlex.status==0) Llex[a,.]=outlex.x[1..n]
        }
    }
    possible=(Lmax:>peertol) :& (Lmax:<.)
    necessary=(Lmin:>peertol) :& (Lmin:<.)
    classes=possible+necessary
    for (a=1; a<=ne; a++) {
        for (b=1; b<=nc; b++) if (Lmax[a,b]>=.) classes[a,b]=.
    }
    np=sum(possible); nn=sum(necessary); no=np-nn
    ni=sum((classes:==0))
    AuditSupport=J(ne,6,.)
    AuditL1=J(ne,5,.)
    AuditIdentification=J(ne,4,.)
    if (audit!="") {
        AuditSupport[,1]=rowsum(selected[ei,.]:>peertol)
        AuditSupport[,2]=rowsum(Lbasic:>peertol)
        AuditSupport[,3]=rowsum(Lsparse:>peertol)
        AuditSupport[,4]=rowsum(possible)
        AuditSupport[,5]=rowsum(possible)
        if (policyweight!="" | policymatrix!="") {
            AuditSupport[,6]=rowsum(Lpolicy:>peertol)
        }
        AuditL1[,1]=rowsum(abs(Lbasic:-selected[ei,.]))
        AuditL1[,2]=rowsum(abs(Lsparse:-selected[ei,.]))
        AuditL1[,3]=rowsum(abs(Lbroad:-selected[ei,.]))
        AuditL1[,4]=rowsum(abs(Lrelative:-selected[ei,.]))
        if (policyweight!="" | policymatrix!="") {
            AuditL1[,5]=rowsum(abs(Lpolicy:-selected[ei,.]))
        }
        AuditIdentification[,1]=rowsum(possible)
        AuditIdentification[,2]=rowsum(possible :& (selected[ei,.]:<=peertol))
        AuditIdentification[,3]=rowsum(Lmax:-Lmin)/n
        AuditIdentification[,4]=rowmax(Lmax:-Lmin)
    }

    st_numscalar("r(aux_failures)",failures)
    st_numscalar("r(lex_failures)",lexfailures)
    st_numscalar("r(sparse_failures)",sparsefailures)
    st_numscalar("r(sparse_nodes)",sparsenodes)
    st_numscalar("r(broad_failures)",broadfailures)
    st_numscalar("r(relative_failures)",relativefailures)
    st_numscalar("r(basic_failures)",basicfailures)
    st_numscalar("r(policy_failures)",policyfailures)
    st_numscalar("r(aux_lp_solves)",auxsolves)
    st_numscalar("r(sparse_lp_solves)",sparselps)
    st_numscalar("r(n_possible)",np); st_numscalar("r(n_necessary)",nn)
    st_numscalar("r(n_optional)",no); st_numscalar("r(n_impossible)",ni)
    st_numscalar("r(opt_tolerance)",opttol)
    st_matrix("r(lambda_min)",Lmin); st_matrix("r(lambda_max)",Lmax)
    st_matrix("r(rts_sum_bounds)",RtsBounds)
    st_matrix("r(rts_robust_code)",RtsCode)
    st_matrix("r(rts_bound_status)",RtsStatus)
    st_matrix("r(rts_margin)",RtsMargin)
    st_matrix("r(rts_stable)",RtsStable)
    st_matrix("r(u0_bounds)",U0Bounds)
    st_matrix("r(u0_rts_code)",U0Code)
    st_matrix("r(u0_bound_status)",U0Status)
    st_matrix("r(u0_rts_agree)",U0Agree)
    st_matrix("r(possible)",possible); st_matrix("r(necessary)",necessary)
    st_matrix("r(reference_status)",classes); st_matrix("r(aux_status)",auxstat)
    st_matrix("r(evaluated_index)",ei); st_matrix("r(candidate_index)",ci)
    st_matrix("r(lambda_selected_subset)",selected[ei,ci])
    st_matrix("r(lambda_lexicographic)",Llex)
    st_matrix("r(lambda_sparse)",Lsparse)
    st_matrix("r(sparse_support)",rowsum(Lsparse:>peertol))
    st_matrix("r(sparse_status)",SparseStatus)
    st_matrix("r(lambda_broad)",Lbroad)
    st_matrix("r(broad_support)",rowsum(possible))
    st_matrix("r(broad_status)",BroadStatus)
    st_matrix("r(lambda_relative)",Lrelative)
    st_matrix("r(relative_solution)",Srelative)
    st_matrix("r(relative_status)",RelativeStatus)
    st_matrix("r(relative_witnesses)",RelativeWitnesses)
    st_matrix("r(lambda_basic)",Lbasic)
    st_matrix("r(basic_solution)",Sbasic)
    st_matrix("r(basic_status)",BasicStatus)
    st_matrix("r(lambda_policy)",Lpolicy)
    st_matrix("r(policy_solution)",Spolicy)
    st_matrix("r(policy_status)",PolicyStatus)
    st_matrix("r(policy_objective)",PolicyObjective)
    st_matrix("r(audit_support)",AuditSupport)
    st_matrix("r(audit_l1_selected)",AuditL1)
    st_matrix("r(audit_identification)",AuditIdentification)
}
end
