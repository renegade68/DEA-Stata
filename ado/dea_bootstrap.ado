*! version 3.0.5  01nov2026
*! 2026-11-01: accept the minimum option abbreviations advertised in help.
*! 2026-11-01: share bootstrap LPs; preserve boundary support in smoothing; repair CI quantiles.
*! 2026-11-01: honor if/in for evaluation and reference samples; reject empty samples.
*! 2026-11-01: preflight output names; explicit replacevars; protect model data.
*! dea_bootstrap.ado — Bootstrap DEA (Simar & Wilson 1998 Algorithm 1)
*! Requires: Stata 16+, _dea_simplex.ado (ships with dea.ado)
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
*!   Banker, R. D., Charnes, A., & Cooper, W. W. (1984). Some models
*!     for estimating technical and scale inefficiencies in data
*!     envelopment analysis. Management Science, 30(9), 1078-1092.
*!     doi:10.1287/mnsc.30.9.1078
*!
*!   Charnes, A., Cooper, W. W., & Rhodes, E. (1978). Measuring the
*!     efficiency of decision making units. European Journal of
*!     Operational Research, 2(6), 429-444.
*!     doi:10.1016/0377-2217(78)90138-8
*!
*!   Simar, L., & Wilson, P. W. (1998). Sensitivity analysis of
*!     efficiency scores: How to bootstrap in nonparametric frontier
*!     models. Management Science, 44(1), 49-61.
*!     doi:10.1287/mnsc.44.1.49
*!
*! Changelog:
*!   3.0.1  21jul2026  Restore documented r(mean_se) stored result.
*!   3.0.0  26apr2026  Switch LP solver to _dea_simplex (Two-Phase Revised
*!                     Simplex with Bland's rule, Harris ratio test, and
*!                     primal→dual fallback).  Stata 16+ LinearProgram had
*!                     persistent false-infeasibility failures on the WHO
*!                     dataset.  The simplex engine is the same solver used
*!                     by dea.ado and has no interior-point convergence issues.
*!   2.3.0  26apr2026  Fixed VRS u0 lower bound (lb[u0]=. → 0 in LinearProgram)
*!   2.2.0  26apr2026  CRITICAL: Fixed pseudo-data ratio direction.
*!                     Input-oriented formula corrected to (theta/beta)*x.
*!                     Verified against simarwilson (Tauchmann & Badunenko).
*!   2.1.0  26apr2026  Switch LP to dual (multiplier) form (m+s+1 vars).
*!   2.0.0  26apr2026  Simar-Wilson smooth bootstrap (KDE with reflection).
*!   1.1.0  11feb2026  Unified syntax; ort() alias
*!   1.0.0  06feb2026  Initial packaged version; based on earlier teaching/seminar versions

// ============================================================================
// DEPENDENCY: load _dea_simplex.ado (Two-Phase Revised Simplex engine)
// Must run BEFORE the mata: block so that RSM_Prob / DEA_Info structs
// are defined when dea_bootstrap's Mata functions are compiled.
// ============================================================================
capture findfile _dea_simplex.ado
if _rc == 0 {
    quietly run `"`r(fn)'"'
}

/**
 * BOOTSTRAP DEA  (Simar & Wilson 1998, Algorithm 1)
 * ==================================================
 *
 * Stage 1  Compute original DEA efficiency scores  theta_i*
 *
 * Stage 2  For b = 1..B:
 *   2a. Draw bootstrap scores  beta_i^b  from the KDE of {theta_i*}
 *       (Gaussian kernel, Silverman bandwidth, reflection at boundary 1,
 *        variance-preserving rescale).
 *   2b. Pseudo-data:
 *         Input-oriented:  x_i^b = (theta_i* / beta_i^b) * x_i,  y_i^b = y_i
 *         Output-oriented: x_i^b = x_i,  y_i^b = (beta_i^b / theta_i*) * y_i
 *   2c. Evaluate original DMU i against pseudo-frontier:
 *         theta_i^{**b} = DEA(x_i, y_i | X^b, Y^b)
 *
 * Stage 3  Bias correction and pivot CI:
 *   theta_i^bc  = 2*theta_i* - mean_b(theta_i^{**b})
 *   CI: [theta_i* - q_{1-a/2}, theta_i* - q_{a/2}]
 *       where q = quantiles of (theta_i^{**b} - theta_i*)
 *
 * Reference:
 *   Simar, L. & Wilson, P.W. (1998). Sensitivity analysis of efficiency
 *   scores: How to bootstrap in nonparametric frontier models.
 *   Management Science, 44(1), 49-61.
 */

capture program drop dea_bootstrap
program define dea_bootstrap, rclass
    version 16.0

    // =========================================================================
    // Unified Syntax Parser
    //   Form A: dea_bootstrap x1 x2 = y1 y2 [if] [in], [options]
    //   Form B: dea_bootstrap [varlist] [if] [in], inputs() outputs() [options]
    // =========================================================================
    _dea_parse_io `0'
    local 0 `"`r(args)'"'

    syntax [varlist(numeric default=none)] [if] [in], ///
        [Inputs(varlist numeric)          ///
         Outputs(varlist numeric)         ///
         DMU(varname)                     ///
         Rts(string)                      ///
         ORT(string) ORientation(string)  ///
         REPS(integer 200)                ///
         Level(real 95)                   ///
         GENerate(string)                 ///
         REPLACEVARS                      ///
         SEED(integer 0)                  ///
         TOLerance(real 1e-8)             ///
         DETail                           ///
         DEBUG]

    // Resolve outputs
    if "`outputs'" == "" {
        if "`varlist'" != "" local outputs `varlist'
        else {
            di as error "outputs() required"
            exit 198
        }
    }
    if "`inputs'" == "" {
        di as error "inputs() required"
        exit 198
    }

    // Reconcile ort/orientation
    if "`ort'" != "" & "`orientation'" != "" {
        di as error "Specify only one of: ort() or orientation()"
        exit 198
    }
    if "`ort'" != "" local orientation "`ort'"

    // Defaults
    if "`rts'"         == "" local rts         "crs"
    if "`orientation'" == "" local orientation "input"
    if "`generate'"    == "" local generate    "_bs_eff"

    if inlist("`orientation'", "i", "in")  local orientation "input"
    if inlist("`orientation'", "o", "out") local orientation "output"

    // Validate
    if !inlist("`rts'", "crs", "vrs", "irs", "drs") {
        di as error "rts() must be crs, vrs, irs, or drs"
        exit 198
    }
    if !inlist("`orientation'", "input", "output") {
        di as error "orientation() must be input or output"
        exit 198
    }
    if `reps' < 50 {
        di as error "reps() should be at least 50"
        exit 198
    }
    if `level' <= 0 | `level' >= 100 {
        di as error "level() must be between 0 and 100"
        exit 198
    }

    // Ensure simplex engine is loaded (always reload to override any mlib version)
    capture findfile _dea_simplex.ado
    if _rc == 0 {
        quietly run `"`r(fn)'"'
    }

    // Validate every output before any existing variable can be replaced.
    _dea_output_check, names(`generate' `generate'_bc `generate'_se `generate'_lo `generate'_hi) ///
        protect(`inputs' `outputs' `dmu') `replacevars'
    local _replace_outputs "`r(existing)'"

    // Mark sample
    marksample touse, novarlist
    markout `touse' `inputs' `outputs'
    if "`dmu'" != "" qui replace `touse' = 0 if missing(`dmu')

    qui count if `touse'
    local n = r(N)
    if `n' == 0 {
        di as error "No observations in estimation sample"
        exit 2000
    }

    if `seed' != 0 set seed `seed'

    // Header
    di as text ""
    di as text "{hline 70}"
    di as text "Bootstrap DEA  (Simar & Wilson 1998)"
    di as text "{hline 70}"
    di as text "Orientation:      " as result "`orientation'"
    di as text "RTS:              " as result "`rts'"
    di as text "DMUs:             " as result "`n'"
    di as text "Bootstrap reps:   " as result "`reps'"
    di as text "Confidence level: " as result "`level'%"
    di as text "LP solver:        " as result "Two-Phase Revised Simplex"
    di as text "{hline 70}"

    // The sample may refer to an existing output; mark it before replacement.
    foreach _ov of local _replace_outputs {
        quietly drop `_ov'
    }

    // Output variables
    qui gen double `generate'     = . if `touse'
    qui gen double `generate'_bc  = . if `touse'
    qui gen double `generate'_se  = . if `touse'
    qui gen double `generate'_lo  = . if `touse'
    qui gen double `generate'_hi  = . if `touse'

    if "`dmu'" == "" {
        tempvar dmu_temp
        qui gen `dmu_temp' = _n if `touse'
        local dmu "`dmu_temp'"
    }

    di as text ""
    di as text "Running bootstrap (Simar-Wilson smooth)... " _continue

    _dea_bootstrap_core, inputs(`inputs') outputs(`outputs') ///
        dmu(`dmu') touse(`touse') eff(`generate') bc(`generate'_bc) ///
        se(`generate'_se) lo(`generate'_lo) hi(`generate'_hi) ///
        orientation(`orientation') rts(`rts') reps(`reps') level(`level') ///
        tolerance(`tolerance') `detail' `debug'

    di as text "done."

    // Summary
    qui summ `generate'    if `touse'
    local mean_orig = r(mean)
    qui summ `generate'_bc if `touse'
    local mean_bc = r(mean)
    qui summ `generate'_se if `touse'
    local mean_se = r(mean)

    di as text ""
    di as text "{hline 70}"
    di as text "Summary"
    di as text "{hline 70}"
    di as text "Mean original efficiency:      " as result %9.4f `mean_orig'
    di as text "Mean bias-corrected:           " as result %9.4f `mean_bc'
    di as text "Mean bootstrap SE:             " as result %9.4f `mean_se'
    di as text "{hline 70}"

    return local orientation "`orientation'"
    return local rts "`rts'"
    return scalar n      = `n'
    return scalar reps   = `reps'
    return scalar level  = `level'
    return scalar mean_eff = `mean_orig'
    return scalar mean_bc  = `mean_bc'
    return scalar mean_se  = `mean_se'
end
