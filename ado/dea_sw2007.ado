*! version 1.1.0  01may2026
*! dea_sw2007.ado — Simar & Wilson (2007) Algorithm 2
*!                  Full Two-Stage DEA with Bootstrap Inference
*! Requires: Stata 16+, _dea_simplex.ado (ships with dea.ado)
*! Optional: simarwilson (Badunenko & Tauchmann 2019, SJ 19(4)) for
*!           backend(simarwilson). Install via: ssc install simarwilson
*!
*! BACKEND OPTION (added v1.1.0)
*!   backend(internal)     — default; self-contained Mata implementation,
*!                            supports panel data and rts(crs|vrs|drs|irs).
*!   backend(simarwilson)  — dispatch to community simarwilson; cross-
*!                            sectional only, rts(crs|vrs|nirs) only.
*!                            Use for numerical-equivalence verification.
*!
*! ALGORITHM (Simar & Wilson 2007, Journal of Econometrics 136: 31-64):
*!
*!  Step 1  Compute original DEA efficiency scores θ_i* for i=1..n
*!          (year-by-year for panel data)
*!
*!  Step 2  Transform: ỹ_i = 1/θ_i* (left-truncated at τ=1)
*!          Regress ỹ_i on Z_i (environmental variables) using
*!          truncated normal MLE → β̂, σ̂
*!
*!  Step 3  Bootstrap loop b = 1..B:
*!   3a.  For each observation i:
*!        Draw ε_i ~ N(0, σ̂²) left-truncated at (τ - Z_i'β̂)
*!        Compute pseudo-inverse-efficiency: ỹ_i^b = Z_i'β̂ + ε_i
*!        Pseudo-efficiency: θ_i^b = 1/ỹ_i^b
*!   3b.  For each year t (panel) or all obs (cross-section):
*!        Generate pseudo-data:
*!          Output-oriented: Y_i^b = Y_i * (θ_i^b/θ_i*)
*!          Input-oriented:  X_i^b = X_i * (θ_i*/θ_i^b)
*!        Re-solve DEA against pseudo-frontier → θ_i^{**b}
*!   3c.  Transform θ_i^{**b} → ỹ_i^{**b} = 1/θ_i^{**b}
*!        Re-run truncated MLE: ỹ^{**b} ~ Z → β̂^b, σ̂^b
*!
*!  Step 4  Bias-corrected estimate:
*!          β̂^bc = 2β̂ - mean_b(β̂^b)
*!          Pivot CI: [β̂ - q_{1-α/2}, β̂ - q_{α/2}]
*!          where q = quantiles of (β̂^b - β̂)
*!
*! Panel support: DEA re-solved separately per time period (Step 3b).
*! Pre-computed DEA scores accepted via scores() option (Step 1 skipped).
*!
*! AUTHORS:
*! --------
*!   Choonjoo Lee (corresponding author)
*!     Korea National Defense University
*!     Seoul, Korea
*!     bloom.rampike@gmail.com; bloom.rampike@gmail.com
*!
*! REFERENCES:
*! -----------
*!   Badunenko, O., & Tauchmann, H. (2019). Simar and Wilson two-stage
*!     efficiency analysis for Stata. Stata Journal, 19(4), 950-988.
*!
*!   Simar, L., & Wilson, P. W. (1998). Sensitivity analysis of efficiency
*!     scores: How to bootstrap in nonparametric frontier models.
*!     Management Science, 44, 49-61.
*!
*!   Simar, L., & Wilson, P. W. (2007). Estimation and inference in
*!     two-stage, semi-parametric models of production processes.
*!     Journal of Econometrics, 136, 31-64.
*!
*!   McDonald, J. (2009). Using least squares and tobit in second stage
*!     DEA efficiency analyses. European Journal of Operational Research,
*!     197, 792-798.

// ============================================================================
// DEPENDENCY: load _dea_simplex.ado
// ============================================================================
capture findfile _dea_simplex.ado
if _rc == 0 {
    quietly run `"`r(fn)'"'
}


// ============================================================================
// MAIN PROGRAM
// ============================================================================

capture program drop dea_sw2007
program define dea_sw2007, eclass
    version 16.0

    // -------------------------------------------------------------------------
    // Syntax
    // -------------------------------------------------------------------------
    // Form A: dea_sw2007 x1 x2 = y1 y2 [if] [in], envars(Z) [options]
    // Form B: dea_sw2007 [if] [in], inputs(x1 x2) outputs(y1 y2) envars(Z) [options]
    //
    // scores(varname)        — pre-computed DEA efficiency scores (skip Step 1)
    // panel(varname)         — panel identifier (re-solve DEA per period)
    // time(varname)          — time identifier (required with panel())
    // backend(internal|simarwilson)
    //                        — internal (default): full Mata implementation,
    //                          supports panel and all rts options.
    //                        — simarwilson: dispatch to community
    //                          simarwilson; cross-sectional + rts(crs|vrs|nirs)
    //                          only.  Used primarily for verification.

    // Detect Form A
    local raw0 `"`0'"'
    if strpos(`"`raw0'"', "=") > 0 {
        local _io_inputs ""
        gettoken _w 0 : 0, parse("=,")
        while "`_w'" != "=" & "`_w'" != "" {
            if "`_w'" == "," {
                di as error "Syntax error: comma before '='"
                error 198
            }
            local _io_inputs `_io_inputs' `_w'
            gettoken _w 0 : 0, parse("=,")
        }
        local _cpos = strpos(`"`0'"', ",")
        if `_cpos' > 0 {
            local _pre  = substr(`"`0'"', 1, `_cpos')
            local _post = substr(`"`0'"', `_cpos'+1, .)
            local 0 `"`_pre' inputs(`_io_inputs') `_post'"'
        }
        else {
            local 0 `"`0', inputs(`_io_inputs')"'
        }
    }

    syntax [varlist(numeric)] [if] [in] ,  ///
        ENVars(varlist numeric)              ///
        [                                   ///
        INputs(varlist numeric)             ///
        OUtputs(varlist numeric)            ///
        SCores(varname numeric)             ///
        PAnel(name)                         ///
        TIme(varname)                       ///
        RTS(string)                         ///
        ORIentation(string) ORT(string)     ///
        REPS(integer 2000)                  ///
        BCREPS(integer 1000)                ///
        Level(integer 95)                   ///
        SEED(integer 0)                     ///
        TAU(real 1.0)                       ///
        GENerate(name)                      ///
        DMU(varname)                        ///
        SAVing(string)                      ///
        REPLACE                             ///
        NOBCorrected                        ///
        TOLerance(real 1e-7)                ///
        TRACE                               ///
        BACKend(string)                     ///
        ]

    // Resolve varlist → outputs (Form B)
    if "`outputs'" == "" {
        if "`varlist'" != "" local outputs `varlist'
    }

    // Orientation
    if "`ort'" != "" & "`orientation'" != "" {
        di as error "Specify only one of: ort() or orientation()"
        exit 198
    }
    if "`ort'" != "" local orientation "`ort'"
    if "`orientation'" == "" local orientation "output"
    local orientation = lower("`orientation'")
    if inlist("`orientation'", "i", "in") local orientation "input"
    if inlist("`orientation'", "o", "out") local orientation "output"

    // RTS default
    if "`rts'" == "" local rts "vrs"
    local rts = lower("`rts'")
    if inlist("`rts'", "bcc") local rts "vrs"
    if inlist("`rts'", "ccr") local rts "crs"

    // Validate
    if !inlist("`rts'", "crs", "vrs", "drs", "irs", "nirs") {
        di as error "rts() must be crs, vrs, nirs, drs, or irs"
        exit 198
    }
    if !inlist("`orientation'", "input", "output") {
        di as error "orientation() must be input or output"
        exit 198
    }
    if `reps' <= 0 {
        di as error "reps() must be a positive integer"
        exit 198
    }
    if `bcreps' <= 0 {
        di as error "bcreps() must be a positive integer"
        exit 198
    }
    if `level' <= 0 | `level' >= 100 {
        di as error "level() must be between 0 and 100"
        exit 198
    }
    if `tau' <= 0 {
        di as error "tau() must be positive"
        exit 198
    }

    // Backend selection (B3 dual-backend support)
    if "`backend'" == "" local backend "internal"
    local backend = lower("`backend'")
    if !inlist("`backend'", "internal", "simarwilson") {
        di as error "backend() must be {it:internal} or {it:simarwilson}"
        exit 198
    }
    if "`backend'" == "internal" & "`rts'" == "nirs" {
        di as error "backend(internal) supports rts(crs|vrs|drs|irs), not rts(nirs)"
        exit 198
    }

    // Check: either scores() or inputs()+outputs() must be specified
    if "`scores'" == "" {
        if "`inputs'" == "" | "`outputs'" == "" {
            di as error "Specify either scores() or both inputs() and outputs()"
            exit 198
        }
    }
    else if "`backend'" == "internal" & ("`inputs'" == "" | "`outputs'" == "") {
        di as error "backend(internal) with scores() still requires inputs() and outputs()"
        di as error "  Algorithm 2 must re-solve DEA on pseudo-frontiers."
        di as error "  Use backend(simarwilson) for score-only Algorithm 1 inference."
        exit 198
    }

    // Panel: both panel() and time() required together
    if ("`panel'" != "" & "`time'" == "") | ("`panel'" == "" & "`time'" != "") {
        di as error "panel() and time() must both be specified"
        exit 198
    }
    if "`panel'" != "" {
        confirm variable `panel'
    }

    // Generate prefix
    if "`generate'" == "" local generate "sw07"
    if "`replace'" != "" {
        capture drop `generate'_eff
        capture drop `generate'_eff_bc
    }
    else {
        capture confirm new variable `generate'_eff
        if _rc {
            di as error "variable `generate'_eff already exists; specify replace or another generate() prefix"
            exit 110
        }
        if "`nobcorrected'" == "" {
            capture confirm new variable `generate'_eff_bc
            if _rc {
                di as error "variable `generate'_eff_bc already exists; specify replace or another generate() prefix"
                exit 110
            }
        }
    }

    // Seed
    if `seed' != 0 set seed `seed'

    // Mark sample
    tempvar touse
    marksample touse, novarlist
    markout `touse' `envars'
    if "`scores'" != "" markout `touse' `scores'
    if "`inputs'"  != "" markout `touse' `inputs'
    if "`outputs'" != "" markout `touse' `outputs'
    if "`panel'"   != "" {
        capture confirm string variable `panel'
        if _rc {
            markout `touse' `panel'
        }
        else {
            qui replace `touse' = 0 if missing(`panel')
        }
    }
    if "`time'"    != "" markout `touse' `time'

    qui count if `touse'
    local n_obs = r(N)

    // Ensure simplex engine is loaded
    capture findfile _dea_simplex.ado
    if _rc == 0 quietly run `"`r(fn)'"'

    // -------------------------------------------------------------------------
    // Header
    // -------------------------------------------------------------------------
    di as text ""
    di as text "{hline 72}"
    di as text "Simar-Wilson (2007) Two-Stage DEA Bootstrap"
    di as text "{hline 72}"
    if "`scores'" != "" {
        di as text "DEA scores:      " as result "Pre-computed (`scores')"
    }
    else {
        di as text "Inputs:          " as result "`inputs'"
        di as text "Outputs:         " as result "`outputs'"
        di as text "RTS:             " as result "`rts'"
        di as text "Orientation:     " as result "`orientation'"
    }
    di as text "Env. variables:  " as result "`envars'"
    if "`panel'" != "" {
        di as text "Panel var:       " as result "`panel'"
        di as text "Time var:        " as result "`time'"
    }
    di as text "Observations:    " as result `n_obs'
    di as text "Bootstrap reps:  " as result `reps'
    di as text "Confidence:      " as result `level' "%"
    di as text "Truncation (τ):  " as result `tau'
    di as text "Backend:         " as result "`backend'"
    di as text "{hline 72}"

    // -------------------------------------------------------------------------
    // Backend dispatch: simarwilson (Badunenko & Tauchmann 2019, SJ 19(4))
    // Internal Mata path resumes at Step 1 below.
    // -------------------------------------------------------------------------
    if "`backend'" == "simarwilson" {
        if "`panel'" != "" {
            di as error ///
                "backend(simarwilson) does not support panel data" _n ///
                "  Use backend(internal) for panel/time."
            exit 198
        }
        if !inlist("`rts'", "crs", "vrs", "nirs") {
            di as error ///
                "backend(simarwilson) supports only rts(crs|vrs|nirs)" _n ///
                "  For rts(`rts'), use backend(internal)."
            exit 198
        }
        capture which simarwilson
        if _rc {
            di as error ///
                "backend(simarwilson) requires the community package " ///
                "{cmd:simarwilson}" _n ///
                "  Install via: {stata ssc install simarwilson:ssc install simarwilson}" _n ///
                "  Reference:   Badunenko & Tauchmann (2019), Stata Journal 19(4)."
            exit 199
        }

        local sw_base = cond("`orientation'"=="output", "output", "input")
        if `seed' != 0 set seed `seed'

        di as text _n ///
            "Dispatching to simarwilson (Badunenko & Tauchmann 2019, SJ 19(4))..."
        di as text "{hline 72}"

        if "`scores'" != "" {
            di as text "Note: pre-computed scores → simarwilson algorithm(1)"
            simarwilson `scores' `envars' if `touse',                       ///
                algorithm(1) reps(`reps') level(`level') noprint
        }
        else {
            simarwilson (`outputs' = `inputs') `envars' if `touse',         ///
                algorithm(2) rts(`rts') base(`sw_base')                     ///
                reps(`reps') bcreps(`bcreps') level(`level') noprint
        }

        // Mark backend in e() for downstream verification scripts
        ereturn local sw07_backend "simarwilson"
        ereturn local sw07_wrapper "dea_sw2007"

        exit 0
    }

    // -------------------------------------------------------------------------
    // Step 1: DEA efficiency scores
    // -------------------------------------------------------------------------
    tempvar eff_orig
    qui gen double `eff_orig' = . if `touse'

    if "`scores'" != "" {
        // Use pre-computed scores
        di as text _n "Step 1: Using pre-computed DEA scores (`scores')"
        qui replace `eff_orig' = `scores' if `touse'
    }
    else {
        // Compute DEA
        di as text _n "Step 1: Computing DEA efficiency scores..."

        if "`panel'" != "" {
            // Year-by-year DEA for panel data
            qui levelsof `time' if `touse', local(periods)
            foreach t of local periods {
                tempvar touse_t
                qui gen byte `touse_t' = (`touse' & `time' == `t')
                qui count if `touse_t'
                if r(N) >= 3 {
                    mata: _sw07_dea_year("`inputs'", "`outputs'",  ///
                        "`touse_t'", "`eff_orig'",                  ///
                        "`orientation'", "`rts'", `tolerance')
                }
                drop `touse_t'
            }
            di as text "  Completed year-by-year DEA"
        }
        else {
            // Single cross-section DEA
            mata: _sw07_dea_year("`inputs'", "`outputs'", "`touse'", ///
                "`eff_orig'", "`orientation'", "`rts'", `tolerance')
        }
    }

    // Clip scores to valid range and check
    qui replace `eff_orig' = min(`eff_orig', 1) if `touse'
    qui replace `eff_orig' = max(`eff_orig', 1e-6) if `touse' & `eff_orig' < 1e-6
    qui count if `eff_orig' == . & `touse'
    if r(N) > 0 {
        di as text "  Warning: " r(N) " missing DEA scores; these will be excluded"
        qui replace `touse' = 0 if `eff_orig' == .
    }
    qui count if `touse'
    local n_obs = r(N)

    qui sum `eff_orig' if `touse'
    di as text "  Mean efficiency: " as result %8.4f r(mean)
    qui count if `eff_orig' >= (1 - `tolerance') & `touse'
    local n_eff = r(N)
    di as text "  Fraction efficient (=1): " as result %5.1f `=100*`n_eff'/`n_obs'' "%"
    di as text "  N efficient: " as result `n_eff' " of " `n_obs'

    // -------------------------------------------------------------------------
    // Step 2: Initial truncated MLE
    // -------------------------------------------------------------------------
    di as text _n "Step 2: Initial truncated regression (truncation at τ = `tau')"

    // Transform: ỹ = 1/θ (left-truncated at τ=1 for eff ∈ (0,1])
    tempvar eff_inv
    qui gen double `eff_inv' = `tau' / `eff_orig' if `touse'
    // Exclude efficient DMUs (at truncation boundary) per SW(2007)
    tempvar touse_trunc
    qui gen byte `touse_trunc' = (`touse' & `eff_inv' > `tau')

    qui count if `touse_trunc'
    local n_trunc = r(N)
    di as text "  Inefficient obs used: " as result `n_trunc' " of " `n_obs'

    // Run initial truncated regression
    qui truncreg `eff_inv' `envars' if `touse_trunc', ll(`tau')
    tempname beta0 V0
    matrix `beta0' = e(b)
    matrix `V0'    = e(V)
    local sigma0   = e(sigma)
    local k_env    = e(k_eq_model)  // number of beta parameters (excl. sigma)
    local npar     = colsof(`beta0')
    local k_beta   = `npar' - 1     // last element of e(b) is sigma

    // Extract variable names for display
    local bnames : colnames `beta0'

    di as text "  Initial β̂:"
    matrix list `beta0', noheader format(%9.4f)
    di as text "  Initial σ̂: " as result %9.4f `sigma0'

    // -------------------------------------------------------------------------
    // Step 3: Bootstrap loop
    // -------------------------------------------------------------------------
    local dea_reps = cond("`nobcorrected'" == "", max(`reps', `bcreps'), `reps')
    di as text _n "Step 3: Bootstrap (" `reps' " inference reps"
    if "`nobcorrected'" == "" & `bcreps' != `reps' {
        di as text "                 " `bcreps' " score bias-correction reps"
    }
    di as text "                 " `dea_reps' " DEA pseudo-frontier reps total)..."
    di as text "  " _continue

    // Matrix to accumulate β̂^b across replications: reps × k_beta
    tempname B_boots S_boots
    matrix `B_boots' = J(`reps', `k_beta', .)
    matrix `S_boots' = J(`reps', 1, .)

    // Extract initial beta (excluding sigma) as Mata column vector
    // beta0 in Stata is: [b_1, b_2, ..., b_k, /sigma]
    // We need b_1..b_k for the mean-function Z*beta
    tempvar xb_hat
    qui predict double `xb_hat' if `touse', xb   // from the truncreg above

    // Bootstrap loop
    tempvar eff_pseudo eff_boot
    qui gen double `eff_pseudo' = . if `touse'
    qui gen double `eff_boot'   = . if `touse'
    tempvar eff_boot_sum eff_boot_cnt
    qui gen double `eff_boot_sum' = 0 if `touse'
    qui gen double `eff_boot_cnt' = 0 if `touse'

    forvalues b = 1/`dea_reps' {
        qui replace `eff_boot' = . if `touse'

        // --- 3a: Draw pseudo-efficiencies from truncated normal ---
        // Draw ε_i ~ N(0, σ̂²) left-truncated at (τ - Z_i'β̂)
        // ỹ_i^b = Z_i'β̂ + ε_i
        // Equivalently: draw ỹ_i^b ~ TN(Z_i'β̂, σ̂², left at τ)
        // via inverse CDF:
        //   U_i ~ Uniform(0,1)
        //   ỹ_i^b = Z_i'β̂ + σ̂ * Φ^{-1}(Φ((τ - Z_i'β̂)/σ̂) + U_i*(1 - Φ((τ-Z_i'β̂)/σ̂)))
        mata: _sw07_rnd_trunc("`xb_hat'", `sigma0', `tau', "`touse'", "`eff_pseudo'")

        // Pseudo-efficiency: θ_i^b = τ/ỹ_i^b (since ỹ_i^b = τ/θ_i^b)
        // → θ_i^b = τ / ỹ_i^b

        // --- 3b: Re-solve DEA against pseudo-frontier ---
        if "`panel'" != "" {
            qui levelsof `time' if `touse', local(periods)
            foreach t of local periods {
                tempvar touse_t
                qui gen byte `touse_t' = (`touse' & `time' == `t')
                qui count if `touse_t'
                if r(N) >= 3 {
                    mata: _sw07_dea_pseudo("`inputs'", "`outputs'",     ///
                        "`eff_orig'", "`eff_pseudo'", "`tau'",          ///
                        "`touse_t'", "`eff_boot'",                      ///
                        "`orientation'", "`rts'", `tolerance')
                }
                drop `touse_t'
            }
        }
        else {
            mata: _sw07_dea_pseudo("`inputs'", "`outputs'",             ///
                "`eff_orig'", "`eff_pseudo'", "`tau'",                  ///
                "`touse'", "`eff_boot'",                                 ///
                "`orientation'", "`rts'", `tolerance')
        }
        if "`nobcorrected'" == "" & `b' <= `bcreps' {
            qui replace `eff_boot_sum' = `eff_boot_sum' + `eff_boot' ///
                if `touse' & !missing(`eff_boot')
            qui replace `eff_boot_cnt' = `eff_boot_cnt' + 1 ///
                if `touse' & !missing(`eff_boot')
        }

        // --- 3c: Re-run truncated regression on bootstrap DEA scores ---
        if `b' <= `reps' {
            tempvar eff_boot_inv touse_boot
            qui gen double `eff_boot_inv' = `tau' / `eff_boot' ///
                if `touse' & `eff_boot' > 0
            qui gen byte `touse_boot' = (`touse' & `eff_boot_inv' > `tau' ///
                & !missing(`eff_boot_inv'))

            capture qui truncreg `eff_boot_inv' `envars' if `touse_boot', ll(`tau')
            if _rc == 0 {
                // Store β̂^b (all except last element /sigma)
                tempname b_boot
                matrix `b_boot' = e(b)
                forvalues j = 1/`k_beta' {
                    matrix `B_boots'[`b', `j'] = `b_boot'[1, `j']
                }
                matrix `S_boots'[`b', 1] = e(sigma)
            }

            drop `eff_boot_inv' `touse_boot'
        }

        // Progress dot every 100 reps
        if mod(`b', 100) == 0 {
            di as text "." _continue
            if mod(`b', 1000) == 0 di as text " `b'" _continue
        }
    }
    di as text " done."

    // -------------------------------------------------------------------------
    // Step 4: Bias correction and pivot CI
    // -------------------------------------------------------------------------
    di as text _n "Step 4: Computing bias-corrected estimates and CIs..."

    local alpha = (100 - `level') / 100

    // beta0 vector (point estimates, excluding sigma)
    tempname beta0_vec beta_bc CI_lo CI_hi bias_vec
    matrix `beta0_vec' = J(`k_beta', 1, .)
    matrix `beta_bc'   = J(`k_beta', 1, .)
    matrix `CI_lo'     = J(`k_beta', 1, .)
    matrix `CI_hi'     = J(`k_beta', 1, .)
    matrix `bias_vec'  = J(`k_beta', 1, .)

    forvalues j = 1/`k_beta' {
        matrix `beta0_vec'[`j', 1] = `beta0'[1, `j']
    }

    // Compute column-wise: bias, BC estimate, pivot CI
    mata: _sw07_postprocess("`B_boots'", "`beta0_vec'",       ///
        "`beta_bc'", "`CI_lo'", "`CI_hi'", "`bias_vec'",      ///
        `alpha', `reps', `k_beta')
    capture scalar drop SW07_sigmab_mean SW07_sigmab_sd SW07_boot_success
    capture matrix drop SW07_sigmab_mean SW07_sigmab_sd SW07_boot_success
    mata: _sw07_col_stats("`S_boots'", "SW07_sigmab_mean", "SW07_sigmab_sd", "SW07_boot_success")
    local sigmab_mean = scalar(SW07_sigmab_mean)
    local sigmab_sd = scalar(SW07_sigmab_sd)
    local boot_success = scalar(SW07_boot_success)

    // -------------------------------------------------------------------------
    // Display results
    // -------------------------------------------------------------------------
    di as text ""
    di as text "{hline 72}"
    di as text "Simar-Wilson (2007) Second-Stage Results"
    di as text "{hline 72}"
    di as text "Dependent: 1/θ (inverse DEA efficiency, left-truncated at τ=`tau')"
    di as text ""

    local bnames_short : word 1 of `bnames'
    local coln : word count `bnames'

    // Header
    di as text %20s "Variable" "  " %10s "β̂" "  " %10s "β̂_bc" ///
        "  " %10s "CI_lo" "  " %10s "CI_hi" "  " %10s "Bias"
    di as text "{hline 72}"

    forvalues j = 1/`k_beta' {
        local vname : word `j' of `bnames'
        // Remove /sigma if present in last position
        if "`vname'" == "/sigma" continue
        di as text %20s "`vname'" "  "     ///
            as result %10.4f `beta0_vec'[`j',1] "  " ///
            %10.4f `beta_bc'[`j',1] "  "           ///
            %10.4f `CI_lo'[`j',1] "  "             ///
            %10.4f `CI_hi'[`j',1] "  "             ///
            %10.4f `bias_vec'[`j',1]
    }
    di as text "{hline 72}"

    // sigma
    di as text %20s "σ̂ (initial)" "  " as result %10.4f `sigma0'
    di as text %20s "σ̂ (bootstrap)" "  " as result %10.4f `sigmab_mean' ///
        "  [SD: " %7.4f `sigmab_sd' ", successful reps: " %9.0f `boot_success' "]"
    di as text "{hline 72}"
    di as text "Note: Positive β̂ indicates lower efficiency (larger 1/θ)."
    di as text "      Negate β̂ to interpret as effect on θ directly."
    di as text "{hline 72}"

    // -------------------------------------------------------------------------
    // Store bias-corrected efficiency scores
    // -------------------------------------------------------------------------
    qui gen double `generate'_eff = `eff_orig' if `touse'
    label var `generate'_eff "DEA efficiency (original)"
    if "`nobcorrected'" == "" {
        qui gen double `generate'_eff_bc = ///
            2*`eff_orig' - (`eff_boot_sum' / `eff_boot_cnt') ///
            if `touse' & `eff_boot_cnt' > 0
        label var `generate'_eff_bc "DEA efficiency, bootstrap bias-corrected"
    }

    // -------------------------------------------------------------------------
    // Store in e()
    // -------------------------------------------------------------------------
    // Reconstruct full e(b) and e(V) for compatibility
    mata: SW07_beta_bc_hold    = st_matrix("`beta_bc'")
    mata: SW07_ci_lo_hold      = st_matrix("`CI_lo'")
    mata: SW07_ci_hi_hold      = st_matrix("`CI_hi'")
    mata: SW07_b_boot_hold     = st_matrix("`B_boots'")
    mata: SW07_sigma_boot_hold = st_matrix("`S_boots'")
    local boot_success_val = `boot_success'

    ereturn clear
    ereturn post `beta0' `V0', depname(`eff_inv') obs(`n_obs')
    mata: st_matrix("SW07_beta_bc",    SW07_beta_bc_hold)
    mata: st_matrix("SW07_b_bc",       SW07_beta_bc_hold)
    mata: st_matrix("SW07_ci_lo",      SW07_ci_lo_hold)
    mata: st_matrix("SW07_ci_hi",      SW07_ci_hi_hold)
    mata: st_matrix("SW07_b_boot",     SW07_b_boot_hold)
    mata: st_matrix("SW07_sigma_boot", SW07_sigma_boot_hold)
    ereturn matrix beta_bc   = SW07_beta_bc
    ereturn matrix b_bc      = SW07_b_bc
    ereturn matrix ci_lo     = SW07_ci_lo
    ereturn matrix ci_hi     = SW07_ci_hi
    ereturn matrix b_boot    = SW07_b_boot
    ereturn matrix sigma_boot = SW07_sigma_boot
    ereturn scalar sigma     = `sigma0'
    ereturn scalar sigma0    = `sigma0'
    ereturn scalar boot_success = `boot_success_val'
    ereturn scalar reps      = `reps'
    ereturn scalar bcreps    = `bcreps'
    ereturn scalar level     = `level'
    ereturn scalar tau       = `tau'
    ereturn scalar N_trunc   = `n_trunc'
    ereturn scalar N         = `n_obs'
    ereturn local  envars    "`envars'"
    ereturn local  rts       "`rts'"
    ereturn local  orientation "`orientation'"
    ereturn local  backend   "internal"
    ereturn local  cmd       "dea_sw2007"
    capture matrix drop SW07_beta_bc SW07_b_bc SW07_ci_lo SW07_ci_hi SW07_b_boot SW07_sigma_boot
    capture scalar drop SW07_sigmab_mean SW07_sigmab_sd SW07_boot_success
    mata: mata drop SW07_beta_bc_hold SW07_ci_lo_hold SW07_ci_hi_hold
    mata: mata drop SW07_b_boot_hold SW07_sigma_boot_hold

    // LaTeX-ready summary
    di as text _n "Variables generated:"
    di as text "  `generate'_eff  — original DEA efficiency scores"
    if "`nobcorrected'" == "" {
        di as text "  `generate'_eff_bc — bootstrap bias-corrected DEA efficiency scores"
    }
    di as text "  Use e(b_bc), e(ci_lo), e(ci_hi) for bias-corrected β and CIs"

end


// ============================================================================
// MATA FUNCTIONS
// ============================================================================

mata:
mata set matastrict off


// ----------------------------------------------------------------------------
// _sw07_dea_year()
// Run DEA for observations in touse_var, store scores in eff_var.
// Uses the Two-Phase Revised Simplex engine from _dea_simplex.ado.
// ----------------------------------------------------------------------------
void _sw07_dea_year(
    string scalar inputs_str,
    string scalar outputs_str,
    string scalar touse_var,
    string scalar eff_var,
    string scalar orientation,
    string scalar rts,
    real scalar   tol
)
{
    real matrix X, Y
    real colvector eff, obs_idx
    real rowvector Xm, Ym
    real scalar n, i, k, th
    struct RSM_Prob  scalar prob
    struct DEA_Info  scalar info
    struct RSM_Result scalar res

    obs_idx = selectindex(st_data(., touse_var))
    n = rows(obs_idx)
    if (n < 3) return

    X = st_data(obs_idx, tokens(inputs_str))
    Y = st_data(obs_idx, tokens(outputs_str))

    // Normalise (units invariant, improves solver stability)
    Xm = mean(X); Ym = mean(Y)
    for (k=1; k<=cols(X); k++) if (Xm[k]<1e-12) Xm[k]=1
    for (k=1; k<=cols(Y); k++) if (Ym[k]<1e-12) Ym[k]=1
    X = X :/ Xm
    Y = Y :/ Ym

    eff = J(n, 1, .)
    for (i = 1; i <= n; i++) {
        _dea_lp_stage1(prob, info, X, Y, X[i,]', Y[i,]', orientation, rts)
        res = _rsm_two_phase(prob, tol, "")
        if (res.status == 0) {
            th = _dea_extract_theta(res, info)
            if (orientation == "output") {
                eff[i] = (th > tol ? 1/th : .)
            }
            else {
                eff[i] = th
            }
        }
    }

    // Store back
    for (i = 1; i <= n; i++) {
        st_store(obs_idx[i], st_varindex(eff_var), eff[i])
    }
}


// ----------------------------------------------------------------------------
// _sw07_rnd_trunc()
// Draw ỹ_i ~ TN(μ_i, σ², left-truncated at τ) using inverse CDF method.
// μ_i comes from xb_var (Z_i'β̂ from truncreg).
// Result stored in out_var.
// ----------------------------------------------------------------------------
void _sw07_rnd_trunc(
    string scalar xb_var,
    real scalar   sigma,
    real scalar   tau,
    string scalar touse_var,
    string scalar out_var
)
{
    real colvector obs_idx, mu, u, z_tau, Phi_tau, arg, draws
    real scalar n, i

    obs_idx = selectindex(st_data(., touse_var))
    n = rows(obs_idx)
    if (n == 0) return

    mu    = st_data(obs_idx, xb_var)          // Z_i'β̂
    u     = runiform(n, 1)                    // U_i ~ Uniform(0,1)
    z_tau = (tau :- mu) :/ sigma              // (τ - μ_i)/σ
    Phi_tau = normal(z_tau)                   // Φ((τ - μ_i)/σ)

    // Inverse CDF of TN(μ, σ², left at τ):
    //   F^{-1}(u) = μ + σ * Φ^{-1}(Φ_tau + u * (1 - Φ_tau))
    arg   = Phi_tau :+ u :* (1 :- Phi_tau)
    // Clip to avoid Φ^{-1}(0) or Φ^{-1}(1)
    arg   = rowmax((arg, J(n,1,1e-12)))
    arg   = rowmin((arg, J(n,1,1-1e-12)))
    draws = mu :+ sigma :* invnormal(arg)

    // Safety: ensure draws > tau
    draws = rowmax((draws, J(n,1,tau+1e-8)))

    for (i = 1; i <= n; i++) {
        st_store(obs_idx[i], st_varindex(out_var), draws[i])
    }
}


// ----------------------------------------------------------------------------
// _sw07_dea_pseudo()
// Re-solve DEA against pseudo-frontier for bootstrap iteration b.
//
// Pseudo-data construction (SW 2007 eq. 3-4):
//   Output-oriented: Y_i^b = Y_i * (θ_i^b / θ_i*)
//   Input-oriented:  X_i^b = X_i * (θ_i* / θ_i^b)
//
// Then evaluate original DMU i against pseudo-frontier → θ_i^{**b}
// ----------------------------------------------------------------------------
void _sw07_dea_pseudo(
    string scalar inputs_str,
    string scalar outputs_str,
    string scalar eff_orig_var,
    string scalar eff_pseudo_var,
    string scalar tau_str,
    string scalar touse_var,
    string scalar eff_boot_var,
    string scalar orientation,
    string scalar rts,
    real scalar   tol
)
{
    real matrix X, Y, X_pseudo, Y_pseudo
    real colvector obs_idx, eff_orig, eff_pseudo, eff_pseudo_raw
    real colvector eff_boot, ratio
    real rowvector Xm, Ym
    real scalar n, i, k, tau_val, th
    struct RSM_Prob  scalar prob
    struct DEA_Info  scalar info
    struct RSM_Result scalar res

    tau_val = strtoreal(tau_str)
    obs_idx = selectindex(st_data(., touse_var))
    n = rows(obs_idx)
    if (n < 3) return

    X = st_data(obs_idx, tokens(inputs_str))
    Y = st_data(obs_idx, tokens(outputs_str))

    eff_orig   = st_data(obs_idx, eff_orig_var)     // θ_i*
    // eff_pseudo = ỹ_i^b (pseudo inverse-efficiency, > tau)
    // Convert to pseudo-efficiency: θ_i^b = tau / ỹ_i^b
    eff_pseudo_raw = st_data(obs_idx, eff_pseudo_var)
    eff_pseudo = tau_val :/ eff_pseudo_raw           // θ_i^b

    // Guard: avoid division by zero
    eff_orig   = editmissing(eff_orig,   1e-6)
    eff_pseudo = editmissing(eff_pseudo, 1e-6)
    eff_orig   = rowmax((eff_orig,   J(n,1,1e-6)))
    eff_pseudo = rowmax((eff_pseudo, J(n,1,1e-6)))

    // Normalise
    Xm = mean(X); Ym = mean(Y)
    for (k=1; k<=cols(X); k++) if (Xm[k]<1e-12) Xm[k]=1
    for (k=1; k<=cols(Y); k++) if (Ym[k]<1e-12) Ym[k]=1
    X = X :/ Xm
    Y = Y :/ Ym

    // Pseudo-data
    if (orientation == "output") {
        ratio    = eff_pseudo :/ eff_orig             // θ^b / θ*
        Y_pseudo = Y :* ratio                         // Y^b = Y * (θ^b/θ*)
        X_pseudo = X
    }
    else {
        ratio    = eff_orig :/ eff_pseudo             // θ* / θ^b
        X_pseudo = X :* ratio                         // X^b = X * (θ*/θ^b)
        Y_pseudo = Y
    }

    // Re-solve DEA: evaluate original DMU against pseudo-frontier
    eff_boot = J(n, 1, .)
    for (i = 1; i <= n; i++) {
        _dea_lp_stage1(prob, info, X_pseudo, Y_pseudo,
                       X[i,]', Y[i,]', orientation, rts)
        res = _rsm_two_phase(prob, tol, "")
        if (res.status == 0) {
            th = _dea_extract_theta(res, info)
            if (orientation == "output") {
                eff_boot[i] = (th > tol ? 1/th : .)
            }
            else {
                eff_boot[i] = th
            }
        }
    }

    // Store
    for (i = 1; i <= n; i++) {
        st_store(obs_idx[i], st_varindex(eff_boot_var), eff_boot[i])
    }
}


// ----------------------------------------------------------------------------
// _sw07_postprocess()
// Bias correction and pivot CI for bootstrap coefficient matrix.
//
// B_boots: reps × k matrix of bootstrap β̂^b
// beta0:   k × 1 point estimates β̂
//
// Output:
//   beta_bc: bias-corrected β̂^bc = 2β̂ - mean_b(β̂^b)
//   ci_lo, ci_hi: pivot CI [β̂ - q_{1-α/2}, β̂ - q_{α/2}]
//   bias_vec: mean_b(β̂^b) - β̂
// ----------------------------------------------------------------------------
void _sw07_postprocess(
    string scalar B_boots_name,
    string scalar beta0_name,
    string scalar beta_bc_name,
    string scalar ci_lo_name,
    string scalar ci_hi_name,
    string scalar bias_name,
    real scalar   alpha,
    real scalar   reps,
    real scalar   k_beta
)
{
    real matrix B_boots
    real colvector beta0, beta_bc, ci_lo, ci_hi, bias_vec
    real colvector boot_j, dev_j, idx
    real scalar j, mean_j, lo_idx, hi_idx, lo_pct, hi_pct

    B_boots  = st_matrix(B_boots_name)   // reps × k
    beta0    = st_matrix(beta0_name)     // k × 1

    beta_bc  = J(k_beta, 1, .)
    ci_lo    = J(k_beta, 1, .)
    ci_hi    = J(k_beta, 1, .)
    bias_vec = J(k_beta, 1, .)

    lo_pct = alpha / 2
    hi_pct = 1 - alpha / 2

    for (j = 1; j <= k_beta; j++) {
        idx = selectindex(B_boots[.,j] :< .)
        if (rows(idx) == 0) continue
        boot_j = B_boots[idx,j]
        if (rows(boot_j) < 10) continue

        mean_j       = mean(boot_j)
        bias_vec[j]  = mean_j - beta0[j]
        beta_bc[j]   = 2*beta0[j] - mean_j

        // Pivot CI deviations: d^b = β̂^b - β̂
        dev_j   = sort(boot_j :- beta0[j], 1)  // sort ascending
        lo_idx  = max((1,          ceil( rows(dev_j) * lo_pct)))
        hi_idx  = min((rows(dev_j), floor(rows(dev_j) * hi_pct)))
        // Pivot: CI = [β̂ - q_{1-α/2}, β̂ - q_{α/2}]
        ci_lo[j] = beta0[j] - dev_j[hi_idx]
        ci_hi[j] = beta0[j] - dev_j[lo_idx]
    }

    st_matrix(beta_bc_name, beta_bc)
    st_matrix(ci_lo_name,   ci_lo)
    st_matrix(ci_hi_name,   ci_hi)
    st_matrix(bias_name,    bias_vec)
}


// ----------------------------------------------------------------------------
// _sw07_col_stats()
// Store mean, standard deviation, and nonmissing count for a one-column matrix.
// ----------------------------------------------------------------------------
void _sw07_col_stats(
    string scalar mat_name,
    string scalar mean_name,
    string scalar sd_name,
    string scalar n_name
)
{
    real colvector x, idx
    real scalar n, i, mu, ss, sd

    x = st_matrix(mat_name)
    idx = selectindex(x :< .)
    if (rows(idx) > 0) x = x[idx]
    else x = J(0, 1, .)
    n = rows(x)
    mu = .
    sd = .
    if (n > 0) {
        mu = 0
        for (i = 1; i <= n; i++) mu = mu + x[i]
        mu = mu / n
    }
    if (n > 1) {
        ss = 0
        for (i = 1; i <= n; i++) ss = ss + (x[i] - mu)^2
        sd = sqrt(ss / (n - 1))
    }

    st_numscalar(mean_name, mu)
    st_numscalar(sd_name,   sd)
    st_numscalar(n_name,    n)
}


end  // end mata
