*! version 3.0.1  21jul2026
*! dea_bootstrap.ado — Bootstrap DEA (Simar & Wilson 1998 Algorithm 1)
*! Requires: Stata 16+, _dea_simplex.ado (ships with dea.ado)
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

    syntax [varlist(numeric)] [if] [in], ///
        [Inputs(varlist numeric)          ///
         Outputs(varlist numeric)         ///
         DMU(varname)                     ///
         RTS(string)                      ///
         ORT(string) ORIentation(string)  ///
         REPS(integer 200)                ///
         Level(real 95)                   ///
         GENerate(string)                 ///
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

    if `seed' != 0 set seed `seed'

    // Mark sample
    tempvar touse
    mark `touse'
    markout `touse' `inputs' `outputs'
    if "`dmu'" != "" qui replace `touse' = 0 if missing(`dmu')

    qui count if `touse'
    local n = r(N)

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

    // Output variables
    cap drop `generate'
    cap drop `generate'_bc
    cap drop `generate'_se
    cap drop `generate'_lo
    cap drop `generate'_hi
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

    mata: _dea_bootstrap_sw("`inputs'", "`outputs'", ///
        "`dmu'", "`touse'",                          ///
        "`generate'", "`generate'_bc", "`generate'_se", ///
        "`generate'_lo", "`generate'_hi",            ///
        "`orientation'", "`rts'", `reps', `level', `tolerance', ///
        "`detail'" != "", "`debug'" != "")

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


// ============================================================================
// MATA
// ============================================================================

mata:


// ----------------------------------------------------------------------------
// _sw_bandwidth()
// Silverman (1986) rule-of-thumb bandwidth for KDE of DEA efficiency scores.
// ----------------------------------------------------------------------------
real scalar _sw_bandwidth(real colvector theta)
{
    real colvector v
    real scalar n, sigma

    v     = select(theta, theta :!= .)
    n     = rows(v)
    sigma = sqrt(variance(v))
    if (n < 2 | sigma == 0) return(0.01)
    return(1.06 * sigma * n^(-0.2))
}


// ----------------------------------------------------------------------------
// _sw_sample()
// Draw n bootstrap efficiency scores from the KDE of theta_star.
// Gaussian kernel, Silverman bandwidth, reflection at boundary 1,
// variance-preserving rescale (Efron-Tibshirani).
// ----------------------------------------------------------------------------
real colvector _sw_sample(real colvector theta_star, real scalar h)
{
    real colvector valid, beta
    real scalar n, n_v, theta_bar, sigma2, scale, i, j, b_raw, b_s
    real scalar max_tries

    valid     = select(theta_star, theta_star :!= .)
    n_v       = rows(valid)
    n         = rows(theta_star)
    theta_bar = mean(valid)
    sigma2    = variance(valid)
    scale     = sqrt(1 + h^2 / max((sigma2, 1e-12)))
    max_tries = 200

    beta = J(n, 1, .)

    for (i = 1; i <= n; i++) {
        if (theta_star[i] == .) continue

        b_raw = .
        for (j = 1; j <= max_tries; j++) {
            b_raw = valid[ceil(max((1, n_v * runiform(1,1))))] + h * rnormal(1,1,0,1)
            if (b_raw > 1) b_raw = 2 - b_raw   // reflect at 1
            if (b_raw > 0) break
            b_raw = .
        }
        if (b_raw == .) b_raw = theta_bar   // fallback

        // Variance-preserving rescale: pulls toward theta_bar
        b_s = theta_bar + (b_raw - theta_bar) / scale

        if (b_s <= 0) b_s = 1e-4
        if (b_s >  1) b_s = 1

        beta[i] = b_s
    }

    return(beta)
}


// ----------------------------------------------------------------------------
// _dea_bs_solve()
// Solve one DEA LP using the Two-Phase Revised Simplex engine from
// _dea_simplex.ado (Bland's rule, Harris ratio test, eta-vector updates,
// primal → multiplier fallback).
//
// Evaluates DMU (x0, y0) against reference set (X_ref, Y_ref).
// Returns Farrell efficiency in (0, ∞):
//   Input-oriented:  theta* in (0,1] in standard DEA; can exceed 1 vs pseudo-
//                    frontier (expected: theta^{**b} >= theta* on average)
//   Output-oriented: 1/phi* in (0,1]  (phi* = output expansion factor >= 1)
// ----------------------------------------------------------------------------
real scalar _dea_bs_solve(real rowvector x0,
                           real rowvector y0,
                           real matrix    X_ref,
                           real matrix    Y_ref,
                           string scalar  orientation,
                           string scalar  rts,
                           real scalar    tol)
{
    struct RSM_Prob  scalar prob, prob_m
    struct DEA_Info  scalar info
    struct MULT_Info scalar minfo
    struct RSM_Result scalar res, res_m
    real scalar theta_star
    real colvector x0c, y0c

    x0c = x0'   // 1×m row → m×1 column
    y0c = y0'   // 1×s row → s×1 column

    // --- Primal envelopment LP (Stage 1: efficiency only) ---
    _dea_lp_stage1(prob, info, X_ref, Y_ref, x0c, y0c, orientation, rts)
    res = _rsm_two_phase(prob, tol, "")

    if (res.status == 0) {
        theta_star = _dea_extract_theta(res, info)
        // input:  theta_star = theta  (min theta)
        // output: theta_star = phi    (= -objval of min(-phi))
    }
    else {
        // --- Multiplier (dual) LP fallback ---
        _dea_mult_build_lp(prob_m, minfo, X_ref, Y_ref, x0c, y0c,
                            orientation, rts)
        res_m = _rsm_two_phase(prob_m, tol, "")
        if (res_m.status != 0) return(.)

        // input:  min -(u'y0 + u0)  → theta = -objval
        // output: min  v'x0 - u0    → phi   =  objval
        theta_star = (orientation == "input" ? -res_m.objval : res_m.objval)
    }

    if (theta_star == . || theta_star < 0) return(.)

    // Output-oriented: convert phi → Farrell (= 1/phi)
    if (orientation == "output") return(theta_star > tol ? 1/theta_star : .)
    return(theta_star)
}


// ----------------------------------------------------------------------------
// _dea_bootstrap_sw()
// Main bootstrap routine implementing Simar & Wilson (1998) Algorithm 1.
// ----------------------------------------------------------------------------
void _dea_bootstrap_sw(string scalar inputs,
                        string scalar outputs,
                        string scalar dmu_var,
                        string scalar touse,
                        string scalar eff_var,
                        string scalar eff_bc_var,
                        string scalar eff_se_var,
                        string scalar eff_lo_var,
                        string scalar eff_hi_var,
                        string scalar orientation,
                        string scalar rts,
                        real scalar   B,
                        real scalar   level,
                        real scalar   tol,
                        real scalar   show_detail,
                        real scalar   show_debug)
{
    // ── Declarations ─────────────────────────────────────────────────────────
    real matrix    X, Y, eff_boot, X_pseudo, Y_pseudo
    real colvector obs_idx, eff_orig, beta_b, ratio
    real colvector eff_bc, eff_se, eff_lo, eff_hi
    string colvector dmu_labels
    real rowvector   boot_i, dev_i, X_scale, Y_scale
    real scalar n, m, s, i, b, k
    real scalar h, alpha, lo_pct, hi_pct, lo_idx, hi_idx
    real scalar mean_boot, bias_i, lw
    string scalar lw_fmt, lw_sep
    string rowvector input_names, output_names

    // ── Load data ─────────────────────────────────────────────────────────────
    input_names  = tokens(inputs)
    output_names = tokens(outputs)
    X            = st_data(., input_names,  touse)
    Y            = st_data(., output_names, touse)
    obs_idx      = selectindex(st_data(., touse))

    if (dmu_var != "" & substr(st_vartype(dmu_var),1,3) == "str") {
        dmu_labels = st_sdata(., dmu_var, touse)
    }
    else if (dmu_var != "") {
        dmu_labels = strofreal(st_data(., dmu_var, touse))
    }
    else {
        dmu_labels = strofreal((1::rows(X)))
    }

    n = rows(X)
    m = cols(X)
    s = cols(Y)

    if (show_debug) printf("\n  n=%g, m=%g, s=%g, B=%g\n", n, m, s, B)

    // ── Normalise inputs and outputs ──────────────────────────────────────────
    // DEA efficiency scores are units-invariant, so normalising does not
    // change results.  Dividing each column by its mean reduces scale
    // differences that slow the simplex solver.
    X_scale = mean(X)
    Y_scale = mean(Y)
    for (k = 1; k <= m; k++) if (X_scale[k] < 1e-12) X_scale[k] = 1
    for (k = 1; k <= s; k++) if (Y_scale[k] < 1e-12) Y_scale[k] = 1
    X = X :/ X_scale
    Y = Y :/ Y_scale

    if (show_debug) {
        printf("  X column means (scale): "); X_scale
        printf("  Y column means (scale): "); Y_scale
    }

    // ── Stage 1: Original efficiency scores ───────────────────────────────────
    eff_orig = J(n, 1, .)
    for (i = 1; i <= n; i++) {
        eff_orig[i] = _dea_bs_solve(X[i,.], Y[i,.], X, Y,
                                     orientation, rts, tol)
    }

    if (show_debug) {
        printf("  Original scores: min=%g  mean=%g  max=%g\n",
               min(eff_orig),
               mean(select(eff_orig, eff_orig:!=.)),
               max(eff_orig))
    }

    // ── Stage 2: Simar-Wilson bootstrap ───────────────────────────────────────
    h = _sw_bandwidth(eff_orig)
    if (show_debug) printf("  KDE bandwidth h = %g\n", h)

    eff_boot = J(n, B, .)

    for (b = 1; b <= B; b++) {
        // 2a. Bootstrap scores from KDE with reflection
        beta_b = _sw_sample(eff_orig, h)

        // 2b. Pseudo-data construction (S&W 1998 Alg. 1 Step 2b)
        //
        // Input-oriented:  x_j^b = (theta_j* / beta_j) * x_j  [ratio = theta/beta]
        // Output-oriented: y_j^b = (beta_j / theta_j*) * y_j  [ratio = beta/theta]
        //
        // When beta_j < theta_j* (most draws), ratio > 1:
        //   Inputs enlarged  → pseudo-frontier easier → theta^{**b} >= theta*
        //   => bias = mean(theta^{**b}) - theta* > 0  (upward bias in DEA) ✓
        //   => theta^bc = 2*theta* - mean(theta^{**b}) < theta* ✓
        //
        // Verified against simarwilson (Tauchmann & Badunenko 2019):
        //   X_star = X_ref * (theta_orig / theta_boot)  for input-oriented ✓
        if (orientation == "input") {
            ratio    = eff_orig :/ beta_b       // theta / beta  (n×1)
            ratio    = editmissing(ratio, 1)    // if eff_orig missing → scale 1
            X_pseudo = ratio :* X
            Y_pseudo = Y
        }
        else {
            ratio    = beta_b :/ eff_orig       // beta / theta  (n×1)
            ratio    = editmissing(ratio, 1)
            X_pseudo = X
            Y_pseudo = ratio :* Y
        }

        // 2c. Evaluate original DMUs against pseudo-frontier
        for (i = 1; i <= n; i++) {
            if (eff_orig[i] != .) {
                eff_boot[i, b] = _dea_bs_solve(X[i,.], Y[i,.],
                                                X_pseudo, Y_pseudo,
                                                orientation, rts, tol)
            }
        }

        if (mod(b, max((1, floor(B/10)))) == 0) {
            printf(".")
            displayflush()
        }
    }

    // ── Stage 3: Bias correction and CI ──────────────────────────────────────
    eff_bc = J(n, 1, .)
    eff_se = J(n, 1, .)
    eff_lo = J(n, 1, .)
    eff_hi = J(n, 1, .)

    alpha  = (100 - level) / 100
    lo_pct = alpha / 2
    hi_pct = 1 - alpha / 2

    for (i = 1; i <= n; i++) {
        if (eff_orig[i] == .) continue

        boot_i = eff_boot[i, .]
        boot_i = select(boot_i, boot_i :!= .)
        if (cols(boot_i) == 0) continue

        mean_boot  = mean(boot_i')
        bias_i     = mean_boot - eff_orig[i]

        // Bias-corrected estimate: theta^bc = 2*theta* - mean(theta^{**b})
        eff_bc[i] = 2 * eff_orig[i] - mean_boot
        eff_se[i] = sqrt(variance(boot_i'))

        // Pivot CI (S&W 1998 Alg. 1 Step 3):
        //   deviations d^{(b)} = theta^{**b} - theta*
        //   CI = [theta* - q_{1-a/2}, theta* - q_{a/2}]
        dev_i  = boot_i :- eff_orig[i]
        dev_i  = sort(dev_i', 1)'     // sort ascending

        lo_idx = max((1,          ceil( cols(dev_i) * lo_pct)))
        hi_idx = min((cols(dev_i), floor(cols(dev_i) * hi_pct)))

        eff_lo[i] = eff_orig[i] - dev_i[hi_idx]
        eff_hi[i] = eff_orig[i] - dev_i[lo_idx]
    }

    // ── Display ───────────────────────────────────────────────────────────────
    lw     = max((3, max(strlen(dmu_labels))))
    lw_fmt = "%" + strofreal(lw) + "s"
    lw_sep = "-" * lw

    printf("\n\n")
    printf("Simar-Wilson Bootstrap DEA Results (%g%% CI)\n", level)
    printf("{hline 80}\n")
    printf(lw_fmt + " | %10s %10s %10s | %10s %10s\n",
           "DMU", "Original", "Bias-Corr", "Std.Err", "CI Lo", "CI Hi")
    printf(lw_fmt + " | %10s %10s %10s | %10s %10s\n",
           lw_sep, "----------", "----------", "----------",
           "----------", "----------")

    for (i = 1; i <= n; i++) {
        printf(lw_fmt + " | %10.4f %10.4f %10.4f | %10.4f %10.4f\n",
               dmu_labels[i],
               eff_orig[i], eff_bc[i], eff_se[i],
               eff_lo[i],   eff_hi[i])
    }
    printf("{hline 80}\n")

    // ── Store in Stata ────────────────────────────────────────────────────────
    for (i = 1; i <= n; i++) {
        st_store(obs_idx[i], st_varindex(eff_var),    eff_orig[i])
        st_store(obs_idx[i], st_varindex(eff_bc_var), eff_bc[i])
        st_store(obs_idx[i], st_varindex(eff_se_var), eff_se[i])
        st_store(obs_idx[i], st_varindex(eff_lo_var), eff_lo[i])
        st_store(obs_idx[i], st_varindex(eff_hi_var), eff_hi[i])
    }
}

end
