*! version 2.1.1 20AUG2026
*! Authors: Choonjoo Lee(Corresponding), Byung-Ki Jung
*! Program: Technology Forecasting using Data Envelopment Analysis
/**
 * HISTORY:
 * ---------------------------------------------------------------------------
 * 2018-01-03(WED): TFDEA Code by B.K. Jung using Lim et.al.(2014) eq.1-10.
 * 2019-08-20(TUE): TFDEA program feasibility in Stata presented 
 *                          at 2019 Oceania Stata Conference.
 * 2023-01.04(WED): Publicly released at 
 *            https://sourceforge.net/p/deas/code/HEAD/tree/trunk/
 * 2025-01-01(WED): Revise Release Date Matrix, Implement Forecasting 
 *            Formula, Error Handling, Result Combination, Enhanced comments 
 *            and structure for clarity, and Help file.
 * 2026-02-24(MON): FIX: ROC now only for SOA products (theta_tk=1) that 
 *            became inefficient (Inman et al. 2006 methodology).
 *            FIX: ROC formula -> gamma = theta_tf^{1/(t_eff - t_k)} (Eq.10).
 *            FIX: t_eff indexing bug in _spec (sumselambda[i]->[n4]).
 *            FIX: Orientation-aware forecasting formula in _spec.
 *            ADD: t_eff_roc column in ROC result matrix.
 * 2026-02-28(FRI): FIX: Full compatibility with dea.ado v2.0 interface.
 *            dea.ado v2.0 returns theta via gen() and lambda via r(lambda),
 *            not r(dearslt). Super-efficiency via model(supereff), not
 *            separate dea_supereff command. rown(dmu) -> rown(`dmu').
 *            Orientation/RTS mapping: OUT->output, IN->input, VRS->vrs, CRS->crs.
 * 2026-05-17(SUN): FIX: Robust parsing of multiple input variables before "=".
 * 2026-07-20(MON): RELEASE 2.0: enforce tk()/dmu(), honor if/in, preserve
 *            caller data, support arbitrary string DMU labels, and gate
 *            intermediate diagnostics behind detail.
 * 2026-08-15(SAT): ADD edate(weighted|min|max).  The min/max choices run a
 *            lexicographic second LP that fixes the optimal efficiency and
 *            minimizes/maximizes the lambda-weighted reference date.
 *            ADD rocfrontier(final|first). Defaults now follow the archived
 *            original R TFDEA example: final training frontier + min date.
 * 2026-08-20(THU): ADD experimental continuation option. When a radial VRS
 *            forecasting LP has no solution, find the closest convex
 *            training-frontier reference using normalized input/output gaps
 *            and translate that gap to forecast time with the observed ROC.
 *            This is an operational proposal motivated by Xue and Harker's
 *            (2002) infeasibility interpretation, not their proposed model
 *            and not a generally established always-feasible DEA model.
 * ---------------------------------------------------------------------------
 */

********************************************************************************
* Dependency Graph and Call Relationships
********************************************************************************
/*
    Dependency Graph:
    tfdea (Main Program)
    +-- tfdea_roc (Rate of Change Calculation)
    |   +-- dea v2.0 (Core DEA Analysis via gen() + r(lambda))
    |   +-- Outputs: `rocrslt` (ROC Results Matrix)
    +-- tfdea_spec (Technology Forecasting)
    |   +-- dea v2.0 model(supereff) (Super-Efficiency Analysis)
    |   +-- Outputs: `specrslt` (Forecasting Results Matrix)
    +-- Outputs: Displayed and Returned Matrices
        +-- `rocrslt` (ROC Results)
        +-- `specrslt` (Forecasting Results)
        +-- Scalars: `mean_roc`, `num_roc_dmu`
    */
	
/*
    Call Relationships:
    1. Main Program (`tfdea`) handles user input, options, and controls execution flow.
    2. Calls `tfdea_roc` for computing rate of change.
    3. Calls `tfdea_spec` for technology forecasting if MODEL is set to SPECIFICATION.
    4. Displays results from `rocrslt` and `specrslt`.
    */

capture program drop tfdea
program define tfdea, rclass
	#del ;
	syntax anything(name=inoutvars equalok) [if] [in]
	, TK(varname) DMU(varname)
	[	   
		RTS(string) 	// return to scale; default : VRS
		ORT(string) 	// orientation; default : OUT
		TF(string)	 	// cutoff year; omit for the last observed year
		TIME(string)	// future time periods to forecast
		EDATE(string)   // effective-date selection: weighted|min|max
		ROCFRONTIER(string) // ROC evaluation frontier: final|first
		CONTinuation       // experimental always-feasible VRS forecast fallback
		DETAIL
	]
	;
	#del cr
	
	// Parse input/output variables
	local invars ""
	gettoken word inoutvars : inoutvars, parse("=")
	while "`word'" != "=" {
		local invars `invars' `word'
		gettoken word inoutvars : inoutvars, parse("=")
	}
	local outvars `inoutvars'
	
	// Error handling
	if ("`invars'" == "" | "`outvars'" == "") {
    display as error "Error: Input and output variables must be specified."
    exit 1
		}
	
	foreach var in `invars' `outvars' {
    capture confirm variable `var'
    if _rc != 0 {
        display as error "Error: Variable `var' not found."
        exit 1
		    }
		}
		
	confirm variable `dmu'
	local varlist `invars' `outvars' `tk'
	marksample touse
	markout `touse' `invars' `outvars' `tk'
	quietly replace `touse' = 0 if missing(`dmu')
	quietly count if `touse'
	if r(N) < 2 {
		display as error "TFDEA requires at least two usable observations"
		exit 2001
	}
	if "`time'" == "" local time 0
	local _time_num = real("`time'")
	if missing(`_time_num') | `_time_num' < 0 {
		display as error "time() must be nonnegative"
		exit 198
	}

	// Default returns to scale setting (case-insensitive)
	if ("`rts'" == "") {
    local rts = "VRS"  // Default to VRS
	} 
	else if (lower("`rts'") == "vrs" | lower("`rts'") == "bcc") {
    local rts = "VRS"
	} 
	else if (lower("`rts'") == "crs" | lower("`rts'") == "ccr") {
    local rts = "CRS"
	} 
	else {
    display as error "Error: Invalid RTS value. Accepted values are VRS/BCC or CRS/CCR."
    exit 1
	}

	// Default orientation setting (case-insensitive)
	if ("`ort'" == "") {
    local ort = "OUT"  // Default to OUT
	}
	else if (lower("`ort'") == "in" | lower("`ort'") == "input") {
    local ort = "IN"
	}
	else if (lower("`ort'") == "out" | lower("`ort'") == "output") {
    local ort = "OUT"
	}
	else {
    display as error "Error: Invalid ORT value. Accepted values are IN/INPUT or OUT/OUTPUT."
    exit 1
	}

	// Map orientation and RTS to dea.ado v2.0 format
	if ("`ort'" == "OUT") local ort_dea "output"
	else local ort_dea "input"
	
	if ("`rts'" == "VRS") local rts_dea "vrs"
	else local rts_dea "crs"

	// Effective-date selection among alternate optimal lambda solutions.
	if "`edate'" == "" local edate "min"
	local edate = lower("`edate'")
	if !inlist("`edate'", "weighted", "min", "max") {
		display as error "edate() must be weighted, min, or max"
		exit 198
	}
	if "`rocfrontier'" == "" local rocfrontier "final"
	local rocfrontier = lower("`rocfrontier'")
	if !inlist("`rocfrontier'", "final", "first") {
		display as error "rocfrontier() must be final or first"
		exit 198
	}

	// Default forecasting model setting
	if "`tf'" == "" {
		local model = "NOFORECAST"
	} 
	else {
		local model = "SPECIFICATION"
	}

	// Debugging: Display the selected model
	if "`detail'" != "" display as text "Selected model: `model'"

	// end of syntax and macro setting	

***************************************************************************
* _roc : Calculate Rate of Change(ROC)
***************************************************************************
preserve
	quietly keep if `touse'
	//  Declare temporary variables and matrix names
		tempname dataset tktemp tkindex fdataset ftkindex dearslt alltheta tkeff
	tempname tkineff thetafinal rsltroc thetainitial roc specrslt elaprslt theta
		tempname I sumlambda refmat teff_roc contlambda

	
** step 1. Create dataset matrix for ROC computation
	tempvar _tfdea_rowkey
	quietly gen str20 `_tfdea_rowkey' = "d" + string(_n)
	set more off
	sort `tk' `_tfdea_rowkey'
	if "`tf'" != "" local tf = real("`tf'")
	else local tf = `tk'[c(N)]
	local foretime = `tf' + `_time_num'
	
	// Prepare for DEA iterations by creating tk index
	mkmat `tk' `invars' `outvars' if (`tk' <= `tf'), matrix(`dataset') rown(`_tfdea_rowkey')
	mkmat `tk' if(`tk' <= `tf'), matrix(`tktemp') rown(`_tfdea_rowkey')
		matrix colnames `tktemp' = tk
		local n1 = rowsof(`tktemp')	// Number of DMUs for ROC computation
		
	// Create tkindex matrix for iterations
	forvalues i = 1/`n1' {	
		if `tktemp'[`i',1] != `tktemp'[`i'+1,1] {
			matrix `tkindex' = nullmat(`tkindex') \ `tktemp'[`i',1]
		}
	}
		matrix colnames `tkindex' = tk
		local n2 = rowsof(`tkindex')	// Number of unique time periods
		
	// Debugging: Verify dataset matrix
	if "`detail'" != "" matrix list `dataset'
		
** step 2. Compute DEA efficiency for each year that DMU appeared	
	local i = 1
	local tki = `tkindex'[`i',1]
	matrix `alltheta' = J(`n1',`n2',.)	// Efficiency scores for each iteration
	while (`tki' <= `tf'){
		if "`detail'" != "" display "Running DEA iteration", `i', "of", `n2'
		
		// Subset data for the current iteration
		quietly clear
		quietly svmat `dataset', names(col)
		quietly keep if (`tk' <= `tki')
		
		// Run DEA (dea.ado v2.0 interface: gen() + r(lambda))
		cap drop _tfdea_theta
		quietly dea `invars' = `outvars', rts(`rts_dea') ort(`ort_dea') gen(_tfdea_theta)
		matrix `refmat' = r(lambda)
		
		local j = c(N)
		
		// Populate efficiency matrix (alltheta)
		// Convert Shephard to Farrell for output-oriented: if theta < 1, store 1/theta
		forvalues row = 1/`j' {
			local th_val = _tfdea_theta[`row']
			if `th_val' >= . {
				matrix `alltheta'[`row',`i'] = 0
			}
			else if (`th_val' < 1) {
				matrix `alltheta'[`row',`i'] = 1 / `th_val'
			}
			else {
				matrix `alltheta'[`row',`i'] = `th_val'
			}
		}
		
		cap drop _tfdea_theta
		
		// Increment iteration
		local i = `i'+1
		local tki = `tkindex'[`i',1]
	}
	
	// Debugging: Verify efficiency matrix
	if "`detail'" != "" matrix list `alltheta'
	
	// Efficiency score (theta) at current time(`tf')
	matrix `thetafinal' = `alltheta'[1...,`n2']
		matrix colnames `thetafinal' = theta_tf 
		
	// Initial efficiency score (theta) when DMU first appeared(`tk')
	matrix `thetainitial' = J(`n1',1,.)
		matrix colnames `thetainitial' = theta_tk
	
	// Assign initial theta (matrix-based optimization)
	forvalues i = 1/`n1' {
    local first_col = 1
    while (`first_col' <= `n2' & `alltheta'[`i', `first_col'] == .) {
        local first_col = `first_col' + 1
		}
    if (`first_col' <= `n2') {
        matrix `thetainitial'[`i', 1] = `alltheta'[`i', `first_col']
		}
	}
	
	// Identify tk_ineff (last SOA year)
	restore, preserve
	
	// Convert matrices to dataset variables
	quietly clear
	quietly svmat `thetainitial', names(col)
	quietly svmat `thetafinal', names(col)
	
	// Generate and compute tk_ineff: the year DMU first became inefficient
	quietly gen tk_ineff = .
	forvalues i = 1/`n1' {
		local first_ineff_yr = .
		forvalues c = 1/`n2' {
			if (`alltheta'[`i',`c'] != . & `alltheta'[`i',`c'] > 1 + 1e-6) {
				local first_ineff_yr = `tkindex'[`c',1]
				continue, break
			}
		}
		if (`first_ineff_yr' < .) {
			quietly replace tk_ineff = `first_ineff_yr' in `i'
		}
	}
	
	// Convert tk_ineff back to a matrix
	mkmat tk_ineff, matrix(`tkineff')
		matrix colnames `tkineff' = tk_ineff

	restore, preserve


** step 3. Computing ROC (Inman et al. 2006 methodology)
	// ---------------------------------------------------------------
	// dea.ado v2.0 model(supereff) does not return lambda matrix.
	// Workaround: For each superseded SOA DMU A, re-run standard DEA
	// at the iteration where A first became inefficient.
	// A's lambda at that iteration references frontier DMUs (incl. new ones)
	// → t_eff = sum(lambda_j * year_j) / sum(lambda_j)
	// → gamma = theta_A^{1/(t_eff - t_A)}
	// ---------------------------------------------------------------

	// Lambda-based tkeff from last iteration (for result matrix)
	matrix `I' = J(`n1',1,1)
	matrix `tkeff' = `refmat' * `tktemp'
	matrix `sumlambda' = `refmat' * `I'
	forvalues i = 1/`n1' {
		if (`sumlambda'[`i',1] != 0) {
			matrix `tkeff'[`i',1] = `tkeff'[`i',1] / `sumlambda'[`i',1]
		}
		else {
			matrix `tkeff'[`i',1] = .
		}
	}
	matrix colnames `tkeff' = tk_eff

	// RoC result vectors
	matrix `roc' = J(`n1', 1, .)
	matrix colnames `roc' = ROC
	matrix `teff_roc' = J(`n1', 1, .)
	matrix colnames `teff_roc' = t_eff_roc

	tempname iter_lambda

	// For each superseded SOA DMU, select the requested evaluation frontier.
	forvalues k = 1/`n1' {
		// Must be SOA at release (theta_tk ≈ 1)
		if (abs(`thetainitial'[`k',1] - 1) > 1e-6) continue
		
		local eval_c = 0
		if "`rocfrontier'" == "first" {
			// First sequential frontier on which the release-efficient DMU is superseded.
			forvalues c = 1/`n2' {
				if (`alltheta'[`k',`c'] != . & `alltheta'[`k',`c'] > 1 + 1e-6) {
					local eval_c = `c'
					continue, break
				}
			}
		}
		else {
			// Archived R TFDEA rule: evaluate at the final training frontier.
			if (`alltheta'[`k',`n2'] != . & `alltheta'[`k',`n2'] > 1 + 1e-6) {
				local eval_c = `n2'
			}
		}
		if (`eval_c' == 0) continue
		
		local year_c = `tkindex'[`eval_c', 1]
		local theta_k = `alltheta'[`k', `eval_c']
		local tk_k = `tktemp'[`k', 1]
		
		if "`detail'" != "" display "  DMU row `k' (year `tk_k') evaluated at `year_c': theta=" %9.6f `theta_k'
		
		// Use Mata LP solver to get theta AND lambda for DMU k
		// Reference set: all DMUs up to year_c
		mata: _tfdea_roc_solve("`dataset'", "`tktemp'", ///
			"`invars'", "`outvars'", `k', `eval_c', ///
			"`ort_dea'", "`rts_dea'", "`edate'")
		
		// Retrieve results from Mata
		local mata_theta = scalar(_tfdea_mata_theta)
		local mata_teff = scalar(_tfdea_mata_teff)
		
		if (`mata_theta' >= . | `mata_teff' >= .) {
			if "`detail'" != "" display "    -> SKIP: Mata solver failed"
			continue
		}
		
		local elapsed = `mata_teff' - `tk_k'
		
		if "`detail'" != "" display "    t_eff=" %7.2f `mata_teff' " elapsed=" %6.2f `elapsed'
		
		if (`elapsed' < 0.01) {
			if "`detail'" != "" display "    -> SKIP: elapsed <= 0 (t_eff not after release)"
			continue
		}
		
		// gamma = theta^{1/elapsed}
		local gamma_k = `theta_k'^(1/`elapsed')
		
		matrix `roc'[`k',1] = `gamma_k'
		matrix `teff_roc'[`k',1] = `mata_teff'
		
		if "`detail'" != "" display "    -> gamma=" %8.6f `gamma_k'
	}

	// Calculate mean ROC
	local roc_sum = 0
	local numroc = 0

	forvalues i = 1/`n1' {
		if (`roc'[`i',1] != .) {
			local roc_sum = `roc_sum' + `roc'[`i',1]
			local numroc = `numroc' + 1
		}
	}

	// Debugging: Display intermediate values
	if "`detail'" != "" display as text "ROC Sum: `roc_sum'"
	if "`detail'" != "" display as text "Number of valid ROC values (numroc): `numroc'"

	// Get mean ROC
	if (`numroc' > 0) {
		local mroc = `roc_sum' / `numroc'
		if "`detail'" != "" display as text "Mean ROC (mroc): `mroc'"
	} 
	else {
    display as error "Error: No valid ROC values to calculate mean ROC."
    exit 1
	}

	// Return mean_roc and num_roc_dmu
	return scalar mean_roc = `mroc'

	if ("`numroc'" != "" & `numroc' > 0) {
		return scalar num_roc_dmu = `numroc'
	} 
	else {
    display as error "Error: num_roc_dmu is undefined or invalid."
    exit 1
	}

	// Combine results
	if (rowsof(`tktemp') == rowsof(`roc')) {
    matrix `rsltroc' = `tktemp', `thetainitial', `thetafinal', `tkeff', `teff_roc', `roc'
    if "`detail'" != "" display as text "Result matrix created successfully."
	} 
	else {
    display as error "Error: Dimension mismatch when combining matrices."
    exit 1
	}

	// Debugging: Verify the result matrix
	if "`detail'" != "" matrix list `rsltroc'

	// Restore and preserve dataset state
	restore, preserve	
	
** step 4. Forecast Technology by MODEL 
	if ("`model'" == "SPECIFICATION"){
		// Call specification model
		_spec, invars(`invars') outvars(`outvars') dataset(`dataset') ///
				rts(`rts') ort(`ort') rts_dea(`rts_dea') ort_dea(`ort_dea') ///
				tk(`tk') tf(`tf') mroc(`mroc') dmu(`dmu') edate(`edate') ///
				`continuation' `detail'
		
		// Save results from _spec
		matrix `specrslt' = r(specrslt)
		if "`continuation'" != "" matrix `contlambda' = r(contlambda)
		
		// Debugging: Verify specification results
		if "`detail'" != "" matrix list `specrslt'
	}
	
	// Combine results with alltheta
	if "`detail'" != "" matrix list `tktemp'
	if "`detail'" != "" matrix list `alltheta'
	
	if (rowsof(`tktemp') == rowsof(`alltheta')) {
    matrix `alltheta' = `tktemp', `alltheta'
	} 
	else {
    display as error "Error: Dimension mismatch between tktemp and alltheta."
    exit 1
	}
	
	// Restore dataset to original state
	restore

** step 5. display results
	display as text "{hline `c(linesize)'}"
	display as result ""
	display as result "DEA Options: RTS(`rts') ORT(`ort') ROCFRONTIER(`rocfrontier') EDATE(`edate')"
	display as result "{hline `c(linesize)'}"
	
	display as result "Rate of Change (ROC) calculation Results: "
	matrix list `rsltroc', noblank nohalf noheader f(%9.6g)
	
	display as result ""
	display as result "Annual Rate of Change (AROC) is: " `mroc' "(Valid DMUs: `numroc')"
	
	// Model-Specific Results
	if ("`model'" == "SPECIFICATION")  {
		display as result ""
		display as result "TFDEA Forecasting Results: "
		matrix list `specrslt', noblank nohalf noheader f(%9.6g)	
		}
	else if ("`model'" == "NOFORECAST") {
    display as result "No forecasting performed. Specify tf for SPECIFICATION"
	}
	
	// Return results for further analysis
	return scalar mean_roc = `mroc'
	return scalar num_roc_dmu = `numroc'
	return local edate "`edate'"
	return local rocfrontier "`rocfrontier'"
	return local continuation "`continuation'"
	return matrix rocrslt = `rsltroc'
	return matrix theta_all = `alltheta'
	
	if ("`model'" == "SPECIFICATION") {
		return matrix forecasting = `specrslt'
		if "`continuation'" != "" return matrix continuation_lambda = `contlambda'
	}
	return matrix lambda = `refmat'
	
end

************************************************************************
* _SPEC : technology forecasting 
************************************************************************
program define _spec, rclass
    syntax, invars(varlist) outvars(varlist) dataset(name) ///
            rts(string) ort(string) rts_dea(string) ort_dea(string) ///
            tk(string) tf(string) mroc(string) ///
            dmu(string) edate(string) [continuation detail]
    
    preserve

    tempname ftkindex fdataset ndataset supertheta superteff tfeff
    tempname tksuptemp tfexp rsltspec contgap contused contlambda tmplambda

    // Parse input parameters
    local tf = real("`tf'")
    local mroc = real("`mroc'")

    // Input validation
    if missing(`tf') {
        display as error "Error: Forecasting time (tf) is invalid or missing."
        exit 1
    }
    if missing(`mroc') {
        display as error "Error: Mean rate of change (mroc) is invalid or missing."
        exit 1
    }

    // Sort dataset by tk and dmu
    sort `tk' `dmu'

    // Create matrices for forecasting
    mkmat `tk' `invars' `outvars' if (`tk' > `tf'), matrix(`fdataset')
    mkmat `tk' if (`tk' > `tf'), matrix(`ftkindex')
    local n3 = rowsof(`ftkindex')

    // Dataset validation
    if (`n3' == 0) {
        display as error "Error: No observations available for forecasting. Check `tk' and `tf'."
        exit 1
    }

    // Initialize result matrices
    matrix `supertheta' = J(`n3', 1, .)
    matrix colnames `supertheta' = SE_tf
    matrix `tfeff' = J(`n3', 1, .)
    matrix colnames `tfeff' = tf_eff
    matrix `tfexp' = J(`n3', 1, .)
    matrix colnames `tfexp' = tf_exp
    matrix `contgap' = J(`n3', 1, 0)
    matrix colnames `contgap' = cont_gap
    matrix `contused' = J(`n3', 1, 0)
    matrix colnames `contused' = cont_used
    local use_cont = ("`continuation'" != "")
    if (`use_cont') matrix `contlambda' = J(`n3', rowsof(`dataset'), .)

    // Super-efficiency computation using Mata LP solver
    forvalues i = 1/`n3' {
        if "`detail'" != "" display "DEA_SUPEREFF iteration ", `i', " of ", `n3'

        // Call Mata solver: super-efficiency with lambda for t_eff
        mata: _tfdea_supereff_solve("`dataset'", "`fdataset'", ///
            "`invars'", "`outvars'", `i', ///
            "`ort_dea'", "`rts_dea'", "`edate'", `use_cont')
        
        // Retrieve results
        local se_val = scalar(_tfdea_mata_se)
        local teff_val = scalar(_tfdea_mata_se_teff)
        local gap_val = scalar(_tfdea_mata_cont_gap)
        local used_val = scalar(_tfdea_mata_cont_used)
        if (`use_cont') {
            matrix `tmplambda' = _tfdea_mata_cont_lambda
            matrix `contlambda'[`i', 1] = `tmplambda'
        }
        
        if (`se_val' >= .) {
            matrix `supertheta'[`i', 1] = .
        }
        else {
            matrix `supertheta'[`i', 1] = `se_val'
        }
        
        if (`teff_val' >= .) {
            matrix `tfeff'[`i', 1] = .
        }
        else {
            matrix `tfeff'[`i', 1] = `teff_val'
        }
        if (`gap_val' < .) matrix `contgap'[`i', 1] = `gap_val'
        if (`used_val' < .) matrix `contused'[`i', 1] = `used_val'

        // Compute forecasted release date (Inman et al. 2006, Eq. 13)
        if (`supertheta'[`i', 1] != . & `mroc' > 1 & `tfeff'[`i', 1] != .) {
            if ("`ort'" == "OUT") {
                // Output-oriented: SE < 1 means technology exceeds frontier
                // t_exp = t_eff + ln(1/SE) / ln(gamma_bar)
                if (`supertheta'[`i', 1] < 1 - 1e-6) {
                    matrix `tfexp'[`i', 1] = `tfeff'[`i', 1] + ///
                        ln(1 / `supertheta'[`i', 1]) / ln(`mroc')
                }
                else {
                    // SE >= 1: product does not exceed training frontier
                    matrix `tfexp'[`i', 1] = `tfeff'[`i', 1]
                }
            }
            else {
                // Input-oriented: SE > 1 means technology exceeds frontier
                if (`supertheta'[`i', 1] > 1 + 1e-6) {
                    matrix `tfexp'[`i', 1] = `tfeff'[`i', 1] + ///
                        ln(`supertheta'[`i', 1]) / abs(ln(`mroc'))
                }
                else {
                    matrix `tfexp'[`i', 1] = `tfeff'[`i', 1]
                }
            }
        }
		else {
            matrix `tfexp'[`i', 1] = .
        }
	}

    // Combine results into a result matrix
    if (`use_cont') {
        matrix `rsltspec' = `ftkindex', `supertheta', `tfeff', `tfexp', ///
            `contgap', `contused'
    }
    else matrix `rsltspec' = `ftkindex', `supertheta', `tfeff', `tfexp'

    restore

    // Return results for further analysis
    return matrix specrslt = `rsltspec'
    if (`use_cont') return matrix contlambda = `contlambda'
end

********************************************************************************
* Copyright and Licensing Information
********************************************************************************

/**
 * Program: tfdea.ado (Technology Forecasting using Data Envelopment Analysis)
 * Version: 2.1.1 (20AUG2026)
 * Copyright (c) 2018-2026 by Choonjoo Lee, Rep. of Korea. All Rights Reserved.
 * This product includes software developed at the Lab. of Choonjoo Lee.
 * Contributor: Byung-Ki Jung.
 *
 * License:
 * Licensed under the Academic Free License version 3.0.
 * SPDX-License-Identifier: AFL-3.0
 * See LICENSE and LICENSE_CODE.txt in the repository root.
 * Preserve the author and contributor notices when redistributing this code.
 * Notice clarified 22oct2026: use the unmodified AFL-3.0 terms, without
 * an additional commercial-use restriction. Computational code is unchanged.
 *
 * Warranty Disclaimer:
 * - The warranty provisions and limitations in AFL-3.0 Sections 7 and 8 apply.
 *
 * Citation:
 * If you use this program in your research, please cite it as follows:
 *  "Lee, C.J. & Jung, B.K. (2026). "TFDEA: A Stata tfdea.ado Program."
 *
 * ---------------------------------------------------------------------------
 * For inquiries, please contact: Choonjoo Lee <bloom.rampike@gmail.com>
 * ---------------------------------------------------------------------------
 */

// ============================================================================
// Mata helper: Solve DEA for a single DMU and return theta + t_eff via lambda
// Based on _sdea_solve_dea_with_lambda() from sdea.ado
// ============================================================================
mata:
void _tfdea_roc_solve(string scalar dataset_name, string scalar tktemp_name,
                       string scalar invars_str, string scalar outvars_str,
                       real scalar dmu_row, real scalar iter_col,
                       string scalar orientation, string scalar rts,
                       string scalar edate)
{
    real matrix dataset, X_ref, Y_ref
    real colvector tktemp, time_ref
    real rowvector x0, y0, lambda_sol, params
    real scalar n, m, s, ninp, nout, theta, t_eff
    real scalar n_ref, i, lambda_sum
    string rowvector inv_tokens, out_tokens
    
    // Read matrices from Stata
    dataset = st_matrix(dataset_name)
    tktemp = st_matrix(tktemp_name)
    
    inv_tokens = tokens(invars_str)
    out_tokens = tokens(outvars_str)
    ninp = cols(inv_tokens)
    nout = cols(out_tokens)
    
    // Get year cutoff from unique years at iteration iter_col
    real colvector all_years, unique_years
    all_years = dataset[., 1]
    unique_years = uniqrows(all_years)
    
    if (iter_col > rows(unique_years)) {
        st_numscalar("_tfdea_mata_theta", .)
        st_numscalar("_tfdea_mata_teff", .)
        return
    }
    
    real scalar year_cutoff
    year_cutoff = unique_years[iter_col]
    
    // Select reference set: rows where year <= year_cutoff
    real colvector sel
    sel = (all_years :<= year_cutoff)
    
    X_ref = select(dataset[., 2..(1+ninp)], sel)
    Y_ref = select(dataset[., (2+ninp)..(1+ninp+nout)], sel)
    time_ref = select(all_years, sel)
    n_ref = rows(X_ref)
    
    // Target DMU
    x0 = dataset[dmu_row, 2..(1+ninp)]
    y0 = dataset[dmu_row, (2+ninp)..(1+ninp+nout)]
    
    m = ninp
    s = nout
    
    // Solve DEA LP
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real rowvector c, lb, ub
    class LinearProgram scalar lp, lp2
    
    lambda_sol = J(1, n_ref, 0)
    theta = .
    
    lp = LinearProgram()
    lp.setTrace("off")
    
    if (orientation == "output") {
        c = (J(1, n_ref, 0), -1)
        A_ineq = (X_ref', J(m, 1, 0))
        b_ineq = x0'
        A_ineq = A_ineq \ (-Y_ref', y0')
        b_ineq = b_ineq \ J(s, 1, 0)
        lb = (J(1, n_ref, 0), 1)
        ub = (J(1, n_ref, .), .)
        
        if (rts == "vrs") {
            A_eq = (J(1, n_ref, 1), 0)
            b_eq = 1
        }
        else {
            A_eq = J(0, n_ref+1, 0)
            b_eq = J(0, 1, 0)
        }
        
        lp.setCoefficients(c)
        lp.setInequality(A_ineq, b_ineq)
        if (rows(A_eq) > 0) lp.setEquality(A_eq, b_eq)
        lp.setBounds(lb, ub)
        lp.setMaxOrMin("min")
        (void) lp.optimize()
        
        if (lp.converged() == 1 & lp.errorcode() == 0) {
            theta = -lp.value()
            params = lp.parameters()
            lambda_sol = params[1..n_ref]
            if (edate == "min" | edate == "max") {
                lp2 = LinearProgram()
                lp2.setTrace("off")
                lp2.setCoefficients((time_ref', 0))
                lp2.setInequality(A_ineq, b_ineq)
                A_eq = A_eq \ (J(1, n_ref, 0), 1)
                b_eq = b_eq \ theta
                lp2.setEquality(A_eq, b_eq)
                lp2.setBounds(lb, ub)
                lp2.setMaxOrMin(edate)
                (void) lp2.optimize()
                if (lp2.converged() == 1 & lp2.errorcode() == 0) {
                    params = lp2.parameters()
                    lambda_sol = params[1..n_ref]
                }
            }
        }
    }
    else {
        c = (1, J(1, n_ref, 0))
        A_ineq = (-x0', X_ref')
        b_ineq = J(m, 1, 0)
        A_ineq = A_ineq \ (J(s, 1, 0), -Y_ref')
        b_ineq = b_ineq \ (-y0')
        lb = (1e-10, J(1, n_ref, 0))
        ub = (1, J(1, n_ref, .))
        
        if (rts == "vrs") {
            A_eq = (0, J(1, n_ref, 1))
            b_eq = 1
        }
        else {
            A_eq = J(0, n_ref+1, 0)
            b_eq = J(0, 1, 0)
        }
        
        lp.setCoefficients(c)
        lp.setInequality(A_ineq, b_ineq)
        if (rows(A_eq) > 0) lp.setEquality(A_eq, b_eq)
        lp.setBounds(lb, ub)
        lp.setMaxOrMin("min")
        (void) lp.optimize()
        
        if (lp.converged() == 1 & lp.errorcode() == 0) {
            theta = lp.value()
            params = lp.parameters()
            lambda_sol = params[2..cols(params)]
            if (edate == "min" | edate == "max") {
                lp2 = LinearProgram()
                lp2.setTrace("off")
                lp2.setCoefficients((0, time_ref'))
                lp2.setInequality(A_ineq, b_ineq)
                A_eq = A_eq \ (1, J(1, n_ref, 0))
                b_eq = b_eq \ theta
                lp2.setEquality(A_eq, b_eq)
                lp2.setBounds(lb, ub)
                lp2.setMaxOrMin(edate)
                (void) lp2.optimize()
                if (lp2.converged() == 1 & lp2.errorcode() == 0) {
                    params = lp2.parameters()
                    lambda_sol = params[2..cols(params)]
                }
            }
        }
    }
    
    // Compute t_eff from lambda
    t_eff = .
    if (theta != .) {
        lambda_sum = sum(lambda_sol)
        if (lambda_sum > 1e-10) {
            t_eff = sum(lambda_sol :* time_ref') / lambda_sum
        }
    }
    
    st_numscalar("_tfdea_mata_theta", theta)
    st_numscalar("_tfdea_mata_teff", t_eff)
}

// ============================================================================
// Mata helper: Super-efficiency LP for forecasting (_spec)
// Based on _sdea_solve_supereff_tfdea() from sdea.ado
// Excludes test DMU from reference set, returns SE + lambda for t_eff
// ============================================================================
void _tfdea_supereff_solve(string scalar dataset_name, string scalar fdataset_name,
                            string scalar invars_str, string scalar outvars_str,
                            real scalar test_idx,
                            string scalar orientation, string scalar rts,
                            string scalar edate, real scalar use_cont)
{
    real matrix dataset, fdataset, X_ref, Y_ref
    real rowvector x0, y0, lambda_sol, params
    real colvector time_ref
    real scalar n_ref, m, s, se, t_eff, lambda_sum, cont_gap, cont_used
    real scalar ninp, nout
    string rowvector inv_tokens, out_tokens
    
    // Read matrices
    dataset = st_matrix(dataset_name)
    fdataset = st_matrix(fdataset_name)
    
    inv_tokens = tokens(invars_str)
    out_tokens = tokens(outvars_str)
    ninp = cols(inv_tokens)
    nout = cols(out_tokens)
    
    // Reference set = training data (year <= tf) — dataset matrix
    // Test DMU = fdataset[test_idx, .]
    X_ref = dataset[., 2..(1+ninp)]
    Y_ref = dataset[., (2+ninp)..(1+ninp+nout)]
    time_ref = dataset[., 1]
    n_ref = rows(X_ref)
    
    // Target DMU (test observation)
    x0 = fdataset[test_idx, 2..(1+ninp)]
    y0 = fdataset[test_idx, (2+ninp)..(1+ninp+nout)]
    
    m = ninp
    s = nout
    
    // Solve super-efficiency LP (test DMU NOT in reference set)
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real rowvector c, lb, ub
    class LinearProgram scalar lp, lp2
    
    lambda_sol = J(1, n_ref, 0)
    se = .
    cont_gap = .
    cont_used = 0
    
    lp = LinearProgram()
    lp.setTrace("off")
    
    if (orientation == "output") {
        // Output-oriented super-eff: min(-φ), φ can be < 1 (advanced)
        c = (J(1, n_ref, 0), -1)
        A_ineq = (X_ref', J(m, 1, 0))
        b_ineq = x0'
        A_ineq = A_ineq \ (-Y_ref', y0')
        b_ineq = b_ineq \ J(s, 1, 0)
        lb = (J(1, n_ref, 0), 1e-6)
        ub = (J(1, n_ref, .), .)
        
        if (rts == "vrs") {
            A_eq = (J(1, n_ref, 1), 0)
            b_eq = 1
        }
        else {
            A_eq = J(0, n_ref+1, 0)
            b_eq = J(0, 1, 0)
        }
        
        lp.setCoefficients(c)
        lp.setInequality(A_ineq, b_ineq)
        if (rows(A_eq) > 0) lp.setEquality(A_eq, b_eq)
        lp.setBounds(lb, ub)
        lp.setMaxOrMin("min")
        (void) lp.optimize()
        
        if (lp.converged() == 1 & lp.errorcode() == 0) {
            se = -lp.value()
            params = lp.parameters()
            if (cols(params) >= n_ref) {
                lambda_sol = params[1..n_ref]
            }
            if (edate == "min" | edate == "max") {
                lp2 = LinearProgram()
                lp2.setTrace("off")
                lp2.setCoefficients((time_ref', 0))
                lp2.setInequality(A_ineq, b_ineq)
                A_eq = A_eq \ (J(1, n_ref, 0), 1)
                b_eq = b_eq \ se
                lp2.setEquality(A_eq, b_eq)
                lp2.setBounds(lb, ub)
                lp2.setMaxOrMin(edate)
                (void) lp2.optimize()
                if (lp2.converged() == 1 & lp2.errorcode() == 0) {
                    params = lp2.parameters()
                    lambda_sol = params[1..n_ref]
                }
            }
        }
    }
    else {
        // Input-oriented super-eff: min θ, θ can be > 1 (advanced)
        c = (1, J(1, n_ref, 0))
        A_ineq = (-x0', X_ref')
        b_ineq = J(m, 1, 0)
        A_ineq = A_ineq \ (J(s, 1, 0), -Y_ref')
        b_ineq = b_ineq \ (-y0')
        lb = (1e-6, J(1, n_ref, 0))
        ub = (., J(1, n_ref, .))
        
        if (rts == "vrs") {
            A_eq = (0, J(1, n_ref, 1))
            b_eq = 1
        }
        else {
            A_eq = J(0, n_ref+1, 0)
            b_eq = J(0, 1, 0)
        }
        
        lp.setCoefficients(c)
        lp.setInequality(A_ineq, b_ineq)
        if (rows(A_eq) > 0) lp.setEquality(A_eq, b_eq)
        lp.setBounds(lb, ub)
        lp.setMaxOrMin("min")
        (void) lp.optimize()
        
        if (lp.converged() == 1 & lp.errorcode() == 0) {
            se = lp.value()
            params = lp.parameters()
            lambda_sol = params[2..cols(params)]
            if (edate == "min" | edate == "max") {
                lp2 = LinearProgram()
                lp2.setTrace("off")
                lp2.setCoefficients((0, time_ref'))
                lp2.setInequality(A_ineq, b_ineq)
                A_eq = A_eq \ (1, J(1, n_ref, 0))
                b_eq = b_eq \ se
                lp2.setEquality(A_eq, b_eq)
                lp2.setBounds(lb, ub)
                lp2.setMaxOrMin(edate)
                (void) lp2.optimize()
                if (lp2.converged() == 1 & lp2.errorcode() == 0) {
                    params = lp2.parameters()
                    lambda_sol = params[2..cols(params)]
                }
            }
        }
    }
    
    // Experimental VRS continuation, used only when the radial BCC LP has no
    // solution.  Retain the convex training frontier and minimize the sum of
    // normalized input excesses and output shortfalls needed to reach it.
    if (se == . & use_cont & rts == "vrs") {
        real matrix Ac, Aec
        real colvector bc, bec, xscale, yscale
        real rowvector cc, lbc, ubc, pc
        real scalar j, gap_opt
        class LinearProgram scalar lpc, lpc2

        xscale = J(m, 1, 1)
        yscale = J(s, 1, 1)
        for (j=1; j<=m; j++) {
            xscale[j] = max((abs(x0[j]), mean(abs(X_ref[.,j])), 1e-10))
        }
        for (j=1; j<=s; j++) {
            yscale[j] = max((abs(y0[j]), mean(abs(Y_ref[.,j])), 1e-10))
        }

        // Variables are lambda, input excess, and output shortfall.
        cc = (J(1, n_ref, 0), (1:/xscale)', (1:/yscale)')
        Ac = (X_ref', -I(m), J(m, s, 0))
        bc = x0'
        Ac = Ac \ (-Y_ref', J(s, m, 0), -I(s))
        bc = bc \ (-y0')
        Aec = (J(1, n_ref, 1), J(1, m+s, 0))
        bec = 1
        lbc = J(1, n_ref+m+s, 0)
        ubc = J(1, n_ref+m+s, .)

        lpc = LinearProgram()
        lpc.setTrace("off")
        lpc.setCoefficients(cc)
        lpc.setInequality(Ac, bc)
        lpc.setEquality(Aec, bec)
        lpc.setBounds(lbc, ubc)
        lpc.setMaxOrMin("min")
        (void) lpc.optimize()

        if (lpc.converged() == 1 & lpc.errorcode() == 0) {
            gap_opt = lpc.value()
            pc = lpc.parameters()
            lambda_sol = pc[1..n_ref]

            // Among equally close frontier facets, apply the same effective-
            // date convention as the ordinary forecasting LP.
            if (edate == "min" | edate == "max") {
                lpc2 = LinearProgram()
                lpc2.setTrace("off")
                lpc2.setCoefficients((time_ref', J(1, m+s, 0)))
                lpc2.setInequality(Ac, bc)
                Aec = Aec \ cc
                bec = bec \ gap_opt
                lpc2.setEquality(Aec, bec)
                lpc2.setBounds(lbc, ubc)
                lpc2.setMaxOrMin(edate)
                (void) lpc2.optimize()
                if (lpc2.converged() == 1 & lpc2.errorcode() == 0) {
                    pc = lpc2.parameters()
                    lambda_sol = pc[1..n_ref]
                }
            }

            cont_gap = max((0, gap_opt))
            cont_used = 1
            // The existing TFDEA equation then converts ln(1+gap) using the
            // empirically estimated mean ROC, preserving its time scale.
            if (orientation == "output") se = 1 / (1 + cont_gap)
            else se = 1 + cont_gap
        }
    }

    // Compute t_eff from lambda (references training frontier)
    t_eff = .
    if (se != .) {
        lambda_sum = sum(lambda_sol)
        if (lambda_sum > 1e-10) {
            t_eff = sum(lambda_sol :* time_ref') / lambda_sum
        }
    }
    
    st_numscalar("_tfdea_mata_se", se)
    st_numscalar("_tfdea_mata_se_teff", t_eff)
    st_numscalar("_tfdea_mata_cont_gap", cont_gap)
    st_numscalar("_tfdea_mata_cont_used", cont_used)
    st_matrix("_tfdea_mata_cont_lambda", lambda_sol)
}
end
