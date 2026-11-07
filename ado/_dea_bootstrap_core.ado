*! version 1.0.0  01nov2026
*! Shared smoothed bootstrap DEA calculation; called after output/sample validation.
*! Bootstrap LPs use the package simplex solver, including when initial scores
*! are supplied by core dea's IPM engine. Author: Choonjoo Lee.

quietly findfile _dea_simplex.ado
quietly run `"`r(fn)'"'

program define _dea_bootstrap_core
    version 16.0
    syntax , Inputs(string) Outputs(string) DMU(varname) TOUSE(varname) ///
        EFF(varname) BC(varname) SE(varname) LO(varname) HI(varname) ///
        ORIentation(string) RTS(string) REPS(integer) Level(real) TOLerance(real) ///
        [INITIAL(varname) DETail DEBUG]

    mata: _dea_bootstrap_sw("`inputs'", "`outputs'", "`dmu'", "`touse'", ///
        "`eff'", "`bc'", "`se'", "`lo'", "`hi'", "`orientation'", "`rts'", ///
        `reps', `level', `tolerance', "`detail'"!="", "`debug'"!="", "`initial'")
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
// variance correction on the mirrored distribution BEFORE reflection.
// ----------------------------------------------------------------------------
real colvector _sw_sample(real colvector theta_star, real scalar h)
{
    real colvector valid, reflected, beta
    real scalar n, n_v, sigma2, scale, i, j, center, candidate, max_tries

    valid = select(theta_star, theta_star :!= .)
    n_v = rows(valid)
    n = rows(theta_star)
    reflected = valid \ (2 :- valid)
    sigma2 = variance(reflected)
    scale = sqrt(1 + h^2 / max((sigma2, 1e-12)))
    max_tries = 200
    beta = J(n, 1, .)

    for (i = 1; i <= n; i++) {
        if (theta_star[i] == .) continue
        candidate = .
        for (j = 1; j <= max_tries; j++) {
            center = reflected[ceil(max((1, 2*n_v*runiform(1,1))))]
            // The mirrored empirical distribution has mean 1. Correct its
            // variance before folding at that boundary, preserving support at 1.
            candidate = 1 + (center + h*rnormal(1,1,0,1) - 1)/scale
            if (candidate > 1) candidate = 2 - candidate
            if (candidate > 0) break
            candidate = .
        }
        if (candidate == .) candidate = theta_star[i]
        beta[i] = candidate
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
                        real scalar   show_debug,
                        string scalar initial_var)
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
    if (initial_var != "") {
        // Core dea retains its selected engine's point estimates and references.
        eff_orig = st_data(., initial_var, touse)
    }
    else {
        eff_orig = J(n, 1, .)
        for (i = 1; i <= n; i++) {
            eff_orig[i] = _dea_bs_solve(X[i,.], Y[i,.], X, Y,
                                         orientation, rts, tol)
        }
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

        // Empirical inverse-CDF quantiles: use ceil for BOTH probabilities.
        // Mixing ceil/floor can reverse endpoints for odd B and narrow levels.
        lo_idx = max((1,          ceil(cols(dev_i) * lo_pct)))
        hi_idx = min((cols(dev_i), ceil(cols(dev_i) * hi_pct)))

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
