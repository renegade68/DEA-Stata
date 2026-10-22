*! _dea_simplex.ado
*! Two-Phase Revised Simplex LP Solver for DEA
*! Version 2.0.6 - 2026-08-05
*!
*! Changelog:
*!   2.0.6  05aug2026  Scale TE input/output columns before simplex LP solves;
*!                      restore slacks to source units on output.  Prevents
*!                      false VRS optima with heterogeneous financial magnitudes.
*!   2.0.5  20jul2026  Align input multiplier u0 sign and DRS/IRS bounds with Cooper
*!   2.0.4  20jul2026  Return the signed multiplier intercept to dea
*!   2.0.2  10feb2026  primal-dual fallback: when envelopment (primal) Stage 1
*!                      fails, automatically retry via multiplier (dual) model;
*!                      by strong duality both yield the same θ*
*!   2.0.1  08feb2026  supereff: add lambda_var for reference set extraction
*!   2.0.0  07feb2026  Initial packaged version; based on earlier teaching/seminar versions
*!
*!*!
*! AUTHORS:
*! --------
*!   Choonjoo Lee (corresponding author)
*!     Korea National Defense University
*!     Seoul, Korea
*!     bloom.rampike@gmail.com; bloom.rampike@gmail.com
*!
*! CONTRIBUTORS TO PRIOR VERSIONS:
*! --------------------------------
*!   Yong-Bae Ji, Kyongrok Lee, Byung-Ihn Lim
*!     (contributors to dea.ado v1.x, dea_supereff.ado, dea_sbm.ado)
*!
*! ORIGINAL SOFTWARE:
*! ------------------
*!   Ji, Y.-B., & Lee, C. (2010). Data envelopment analysis.
*!     The Stata Journal, 10(2), 267-280.
*!     https://doi.org/10.1177/1536867X1001000207
*!
*!   Lee, C., & Ji, Y.-B. (2009). Data envelopment analysis in Stata.
*!     Presented at the 2009 Stata Conference, Washington, DC.
*!     https://ideas.repec.org/p/boc/dcon09/4.html
*!
*! NOTE: This Simplex engine is a complete rewrite of ldea.do (the
*!   original Mata LP solver from Ji & Lee, 2010). No code is copied;
*!   the two-phase structure is a standard LP technique. See below
*!   for specific algorithmic improvements over ldea.do.
*!
*! ALGORITHM REFERENCES:
*! ---------------------
*! 1. Dantzig, G. B. (1963). Linear Programming and Extensions.
*!    Princeton University Press.
*!    - Two-Phase Simplex Method, Revised Simplex Method
*!
*! 2. Chen, D.-S., Batson, R. G., & Dang, Y. (2010).
*!    Applied Integer Programming. Wiley. Section 9.8.
*!    - Revised simplex tableau: B^{-1} maintenance, pivot column/row
*!      computation, and systematic update procedure
*!
*! 3. Bazaraa, M. S., Jarvis, J. J., & Sherali, H. D. (2010).
*!    Linear Programming and Network Flows (4th ed.). Wiley.
*!    - Product-form (eta-vector) updates, degeneracy / cycling theory
*!
*! 4. Cooper, W. W., Seiford, L. M., & Tone, K. (2007).
*!    Data Envelopment Analysis: A Comprehensive Text (2nd ed.). Springer.
*!    - DEA envelopment LP formulations (CRS/VRS/IRS/DRS)
*!    - Cost, revenue, profit efficiency (Ch.8)
*!
*! 5. Bland, R. G. (1977). New finite pivoting rules for the simplex method.
*!    Mathematics of Operations Research, 2(2), 103-107.
*!    - Anti-cycling pivot rule
*!
*! 6. Harris, P. M. J. (1973). Pivot selection methods of the Devex LP
*!    code. Mathematical Programming, 5(1), 1-28.
*!    - Harris ratio test (two-pass leaving variable selection)
*!
*! DEA MODEL REFERENCES (see dea.ado header for complete list):
*! ------------------------------------------------------------
*!   Charnes, Cooper, & Rhodes (1978) — CCR (CRS)
*!   Banker, Charnes, & Cooper (1984) — BCC (VRS)
*!   Andersen & Petersen (1993) — super-efficiency
*!   Chambers, Chung, & Färe (1996) — DDF
*!   Tone (2001) — SBM
*!   Tone (2002) — Super-SBM
*!   Charnes & Cooper (1962) — Charnes-Cooper transformation
*!
*! IMPROVEMENTS OVER ldea.do:
*! -------------------------
*! - Vectorized pricing: y=c_B*B^{-1}, then rc=c_N-y*N (single matmul)
*!   (ldea.do computes reduced costs in element-wise loop)
*! - Proper Bland's rule with early activation on degeneracy detection
*!   (ldea.do activates Bland's only after 150+ stalled iterations)
*! - Harris ratio test with largest-pivot tie-breaking for stability
*!   (ldea.do uses simple min-ratio without Harris tolerance)
*! - Eta-vector (product form) basis updates: O(m^2) per iteration
*!   (ldea.do recomputes full B^{-1} via lusolve/svsolve EACH iteration)
*! - Periodic basis refactorization for numerical hygiene
*! - Clean Phase1->Phase2 transition with artificial purging
*! - No BoundCond struct overhead; clean standard-form LP throughout
*!
*! INTERFACE:
*! ----------
*! Compatible with dea.ado. Provides two entry points:
*!   _dea_simplex_envelopment()  - efficiency scores only (stores in var)
*!   _dea_simplex_full()         - full results in r(dea_results) matrix

version 16.0

// Clean up existing Mata structs and functions
capture mata: mata drop RSM_Prob()
capture mata: mata drop RSM_Basis()
capture mata: mata drop RSM_Result()
capture mata: mata drop DEA_Info()
capture mata: mata drop _dea_simplex_envelopment()
capture mata: mata drop _dea_simplex_full()
capture mata: mata drop _dea_solve_dmu()
capture mata: mata drop _dea_lp_stage1()
capture mata: mata drop _dea_lp_stage2()
capture mata: mata drop _rsm_two_phase()
capture mata: mata drop _rsm_solve()
capture mata: mata drop _rsm_enter()
capture mata: mata drop _rsm_leave()
capture mata: mata drop _rsm_pivot()
capture mata: mata drop _rsm_refac()
capture mata: mata drop _rsm_purge_artif()
capture mata: mata drop _rsm_drop_artif()
capture mata: mata drop _rsm_inv()
capture mata: mata drop _dea_extract()
capture mata: mata drop _dea_extract_theta()
capture mata: mata drop _dea_tol_pivot()
capture mata: mata drop _dea_tol_feas()
capture mata: mata drop _dea_tol_opt()
capture mata: mata drop _dea_tol_zero()
capture mata: mata drop _dea_refac_freq()
capture mata: mata drop _dea_harris_tol()
capture mata: mata drop MULT_Info()
capture mata: mata drop _dea_simplex_mult_envelopment()
capture mata: mata drop _dea_simplex_mult_full()
capture mata: mata drop _dea_simplex_multiplier()
capture mata: mata drop _dea_mult_solve_dmu()
capture mata: mata drop _dea_mult_build_lp()
capture mata: mata drop _dea_mult_extract()
capture mata: mata drop _rsm_std_form()
capture mata: mata drop _rsm_solve_lp()
capture mata: mata drop _dea_simplex_cost()
capture mata: mata drop _dea_simplex_revenue()
capture mata: mata drop _dea_simplex_profit()
capture mata: mata drop _dea_simplex_ddf()
capture mata: mata drop _dea_simplex_ddf_undesirable()
capture mata: mata drop _dea_simplex_sbm()
capture mata: mata drop _dea_simplex_additive()
capture mata: mata drop _dea_simplex_supereff()
capture mata: mata drop _dea_simplex_supersbm()
capture mata: mata drop _dea_simplex_allocative()

mata:
mata set matastrict on

// ============================================================================
// STRUCT DEFINITIONS
// ============================================================================

/**
 * RSM_Prob: Standard-form LP problem
 *   min  c'x   s.t.  Ax = b,  x >= 0
 * b >= 0 is guaranteed for DEA problems.
 */
struct RSM_Prob {
    real matrix    A           // m x n constraint matrix
    real colvector b           // m x 1 RHS vector (>= 0)
    real rowvector c           // 1 x n objective coefficients
    real scalar    m_con       // number of constraints
    real scalar    n_var       // total variables including artificials
    real scalar    n_orig      // original variables (before artificials)
    real rowvector is_artif    // 1 x n: 1 if artificial, 0 otherwise
    real colvector init_basic  // m x 1: initial basic variable index per row
}

/**
 * RSM_Basis: Revised simplex basis state
 * Stores B^{-1} explicitly; updated via eta-vectors between refactorizations.
 */
struct RSM_Basis {
    real colvector basic       // m x 1: column indices of basic variables
    real matrix    Binv        // m x m: current basis inverse
    real scalar    n_upd       // eta updates since last refactorization
}

/**
 * RSM_Result: LP solution result
 */
struct RSM_Result {
    real scalar    status      // 0=optimal, 1=infeasible, 2=unbounded, 3=maxiter
    real scalar    objval      // objective function value
    real colvector x           // n x 1 primal solution
    real scalar    iters       // total simplex iterations
    string scalar  msg         // human-readable status
}

/**
 * DEA_Info: DEA problem metadata for result extraction
 */
struct DEA_Info {
    real scalar    n_dmu       // number of DMUs
    real scalar    n_inp       // number of inputs
    real scalar    n_out       // number of outputs
    real scalar    idx_theta   // column index of theta/phi (0 if absent)
    real rowvector idx_lam     // column indices of lambda_1..lambda_n
    real rowvector idx_si      // column indices of input slacks s_1^-..s_m^-
    real rowvector idx_so      // column indices of output slacks s_1^+..s_s^+
    real scalar    idx_sc      // column index of convexity slack (0 if none)
    string scalar  ort         // "input" or "output"
    string scalar  rts         // "crs", "vrs", "irs", "drs"
}

/**
 * MULT_Info: Multiplier (dual) model metadata
 *
 * LP variables: [v_1..v_m, u_1..u_s, (u0+), (u0-), s_1..s_n]
 *   v_i  = input weights  (dual of input constraint)
 *   u_r  = output weights (dual of output constraint)
 *   u0   = free variable  (dual of convexity constraint)
 *          split u0 = u0+ - u0-  for standard form
 *   s_j  = DMU constraint slack
 */
struct MULT_Info {
    real scalar    n_dmu       // number of DMUs
    real scalar    n_inp       // number of inputs (m)
    real scalar    n_out       // number of outputs (s)
    string scalar  ort         // "input" or "output"
    string scalar  rts         // "crs", "vrs", "irs", "drs"
    real rowvector idx_v       // input weight variable indices
    real rowvector idx_u       // output weight variable indices
    real scalar    idx_u0p     // u0+ index (0 if absent)
    real scalar    idx_u0m     // u0- index (0 if absent)
    real rowvector idx_s       // DMU slack variable indices
}


// ============================================================================
// NUMERICAL CONSTANTS (as functions to avoid Mata global issues)
// ============================================================================
real scalar _dea_tol_pivot()  return(1e-10)   // minimum pivot magnitude
real scalar _dea_tol_feas()   return(1e-8)    // feasibility tolerance
real scalar _dea_tol_opt()    return(1e-8)    // optimality tolerance
real scalar _dea_tol_zero()   return(1e-12)   // zero threshold
real scalar _dea_refac_freq() return(40)      // refactorize interval
real scalar _dea_harris_tol() return(1e-7)    // Harris ratio tolerance


// ============================================================================
// MAIN ENTRY: _dea_simplex_envelopment()
// ============================================================================
/**
 * Solves DEA envelopment model for all DMUs. Stores efficiency in eff_var.
 * Called from dea.ado as drop-in alternative to the IPM solver.
 *
 * @param inputs         Space-separated input variable names
 * @param outputs        Space-separated output variable names
 * @param eff_var        Stata variable to store efficiency scores
 * @param touse          Stata variable marking sample observations
 * @param orientation    "input" or "output"
 * @param rts            "crs", "vrs", "irs", "drs"
 * @param tol            Numerical tolerance
 * @param lambda_vars    Space-separated variable names for lambda weights ("" if none)
 * @param slack_in_vars  Space-separated variable names for input slacks ("" if none)
 * @param slack_out_vars Space-separated variable names for output slacks ("" if none)
 */
void _dea_simplex_envelopment(string scalar inputs,
                               string scalar outputs,
                               string scalar eff_var,
                               string scalar touse,
                               string scalar orientation,
                               string scalar rts,
                               real scalar tol,
                               string scalar lambda_vars,
                               | string scalar slack_in_vars,
                                 string scalar slack_out_vars)
{
    real matrix X, Y
    real rowvector x_scale, y_scale
    real colvector obs_indices, result_k
    real scalar n_dmu, n_inp, n_out, k, obs_k, j
    string rowvector lv, siv, sov
    string scalar _lv_j, _sv_j
    real scalar has_lambda, has_slacks, _val_j

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))
    n_dmu = rows(X)
    n_inp = cols(X)
    n_out = cols(Y)

    // DEA radial scores and lambdas are invariant to positive changes of units.
    // Scaling each input/output column to a maximum of one prevents the revised
    // simplex absolute tolerances from declaring false optima when financial
    // statement columns have very different magnitudes.  Slacks are converted
    // back to the source units before storage below.
    x_scale = J(1, n_inp, 1)
    y_scale = J(1, n_out, 1)
    for (j = 1; j <= n_inp; j++) {
        x_scale[j] = max(abs(X[., j]))
        if (x_scale[j] <= 0 | x_scale[j] >= .) x_scale[j] = 1
        X[., j] = X[., j] / x_scale[j]
    }
    for (j = 1; j <= n_out; j++) {
        y_scale[j] = max(abs(Y[., j]))
        if (y_scale[j] <= 0 | y_scale[j] >= .) y_scale[j] = 1
        Y[., j] = Y[., j] / y_scale[j]
    }

    has_lambda = (lambda_vars != "")
    if (has_lambda) lv = tokens(lambda_vars)

    has_slacks = (args() >= 9 ? (slack_in_vars != "" & slack_out_vars != "") : 0)
    if (has_slacks) {
        siv = tokens(slack_in_vars)
        sov = tokens(slack_out_vars)
    }

    for (k = 1; k <= n_dmu; k++) {
        obs_k = obs_indices[k]
        result_k = _dea_solve_dmu(X, Y, k, orientation, rts, tol, "")
        st_store(obs_k, st_varindex(eff_var), (result_k[1]))
        // Lambda: result_k[2..n+1]
        if (has_lambda) {
            for (j = 1; j <= min((n_dmu, cols(lv))); j++) {
                if (rows(result_k) >= j + 1) {
                    _lv_j = lv[j]
                    _val_j = (result_k[j + 1])
                    st_store(obs_k, st_varindex(_lv_j), _val_j)
                }
            }
        }
        // Slacks: result_k[n+2..n+1+m] = input slacks, result_k[n+2+m..n+1+m+s] = output slacks
        if (has_slacks) {
            for (j = 1; j <= min((n_inp, cols(siv))); j++) {
                if (rows(result_k) >= n_dmu + 1 + j) {
                    _sv_j = siv[j]
                    _val_j = (result_k[n_dmu + 1 + j]) * x_scale[j]
                    st_store(obs_k, st_varindex(_sv_j), _val_j)
                }
            }
            for (j = 1; j <= min((n_out, cols(sov))); j++) {
                if (rows(result_k) >= n_dmu + 1 + n_inp + j) {
                    _sv_j = sov[j]
                    _val_j = (result_k[n_dmu + 1 + n_inp + j]) * y_scale[j]
                    st_store(obs_k, st_varindex(_sv_j), _val_j)
                }
            }
        }
    }
}


// ============================================================================
// FULL ENTRY: _dea_simplex_full()
// ============================================================================
/**
 * Full DEA solution: efficiency + lambdas + slacks for all DMUs.
 * Result stored in r(dea_results): n x (1 + n + m + s [+1]) matrix.
 *
 * Row format: [theta, lambda_1..n, sinp_1..m, sout_1..s, (sconv)]
 */
void _dea_simplex_full(string scalar inputs,
                        string scalar outputs,
                        string scalar touse,
                        string scalar orientation,
                        string scalar rts,
                        real scalar tol,
                        | string scalar trace)
{
    real matrix X, Y, RESULTS
    real rowvector x_scale, y_scale
    real colvector result_k
    real scalar n_dmu, n_inp, n_out, k, nc, nr_k, j
    string scalar do_trace

    do_trace = (args() >= 7 ? trace : "")

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    n_dmu = rows(X)
    n_inp = cols(X)
    n_out = cols(Y)

    x_scale = J(1, n_inp, 1)
    y_scale = J(1, n_out, 1)
    for (j = 1; j <= n_inp; j++) {
        x_scale[j] = max(abs(X[., j]))
        if (x_scale[j] <= 0 | x_scale[j] >= .) x_scale[j] = 1
        X[., j] = X[., j] / x_scale[j]
    }
    for (j = 1; j <= n_out; j++) {
        y_scale[j] = max(abs(Y[., j]))
        if (y_scale[j] <= 0 | y_scale[j] >= .) y_scale[j] = 1
        Y[., j] = Y[., j] / y_scale[j]
    }

    nc = 1 + n_dmu + n_inp + n_out
    if (rts == "irs" || rts == "drs") nc = nc + 1
    RESULTS = J(n_dmu, nc, 0)

    for (k = 1; k <= n_dmu; k++) {
        if (do_trace == "trace") {
            printf("\n{txt}=== DMU %g / %g ===\n", k, n_dmu)
        }
        result_k = _dea_solve_dmu(X, Y, k, orientation, rts, tol, do_trace)
        nr_k = rows(result_k)
        if (nr_k >= nc) {
            RESULTS[k, .] = result_k[1..nc]'
        }
        else if (nr_k > 0) {
            RESULTS[k, 1..nr_k] = result_k'
        }
        // Return stage-2 slacks in the source variables' units.
        if (nr_k >= 1 + n_dmu + n_inp + n_out) {
            for (j = 1; j <= n_inp; j++) {
                RESULTS[k, 1 + n_dmu + j] = RESULTS[k, 1 + n_dmu + j] * x_scale[j]
            }
            for (j = 1; j <= n_out; j++) {
                RESULTS[k, 1 + n_dmu + n_inp + j] = RESULTS[k, 1 + n_dmu + n_inp + j] * y_scale[j]
            }
        }
    }

    st_matrix("r(dea_results)", RESULTS)
}


// ============================================================================
// SOLVE SINGLE DMU: Two-stage DEA procedure
// ============================================================================
/**
 * _dea_solve_dmu()
 *
 * Two-stage DEA solution for DMU k with primal-dual fallback:
 *   Stage 1: Optimize efficiency (min theta or max phi) via envelopment
 *   Fallback: If Stage 1 fails, re-solve via multiplier (dual) model
 *   Stage 2: Maximize total slacks with efficiency fixed at theta*
 *
 * By strong duality the envelopment (primal) and multiplier (dual)
 * formulations yield the same optimal θ*.  When the primal LP encounters
 * numerical difficulties the dual formulation may still succeed because
 * it has a fundamentally different constraint structure (n DMU constraints
 * vs. m+s technology constraints).
 *
 * @return Column vector [theta, lam_1..n, si_1..m, so_1..s, (sc)]
 */
real colvector _dea_solve_dmu(real matrix X,
                               real matrix Y,
                               real scalar k,
                               string scalar ort,
                               string scalar rts,
                               real scalar tol,
                               string scalar trace)
{
    real scalar n_dmu, n_inp, n_out, nc, theta_star
    real colvector x0, y0, result

    struct RSM_Prob scalar prob1, prob2
    struct DEA_Info scalar info1, info2
    struct RSM_Result scalar res1, res2

    // Primal-dual fallback variables
    struct RSM_Prob scalar prob_m
    struct MULT_Info scalar minfo
    struct RSM_Result scalar res_m

    n_dmu = rows(X)
    n_inp = cols(X)
    n_out = cols(Y)
    x0 = X[k, .]'
    y0 = Y[k, .]'
    nc = 1 + n_dmu + n_inp + n_out
    if (rts == "irs" || rts == "drs") nc = nc + 1

    // =================================================================
    // STAGE 1: Optimize efficiency
    // =================================================================
    _dea_lp_stage1(prob1, info1, X, Y, x0, y0, ort, rts)

    res1 = _rsm_two_phase(prob1, tol, trace)

    if (res1.status != 0) {
        // ─────────────────────────────────────────────────────────────
        // PRIMAL-DUAL FALLBACK
        // Envelopment (primal) failed → try Multiplier (dual) model.
        // By strong duality, both formulations yield the same θ*.
        // If multiplier succeeds, use its θ* and attempt Stage 2
        // to recover lambdas and slacks.
        // ─────────────────────────────────────────────────────────────
        if (trace == "trace") {
            printf("{err}DMU %g Stage1 (envelopment) FAIL: %s\n", k, res1.msg)
            printf("{txt}DMU %g: Attempting multiplier (dual) fallback...\n", k)
        }

        _dea_mult_build_lp(prob_m, minfo, X, Y, x0, y0, ort, rts)
        res_m = _rsm_two_phase(prob_m, tol, trace)

        if (res_m.status != 0) {
            // Both primal and dual failed
            if (trace == "trace") {
                printf("{err}DMU %g: Multiplier fallback also FAILED: %s\n",
                       k, res_m.msg)
            }
            return(J(nc, 1, .))
        }

        // Extract θ* from multiplier objective
        //   Input:  min -(u'y₀+u₀), so θ* = -objval
        //   Output: min  (v'x₀+u₀), so φ* =  objval
        if (ort == "input") {
            theta_star = -res_m.objval
        }
        else {
            theta_star = res_m.objval
        }

        if (trace == "trace") {
            printf("{txt}DMU %g: Multiplier fallback OK, theta* = %12.8f\n",
                   k, theta_star)
        }

        // Attempt Stage 2 with θ* from multiplier to recover λ and slacks
        _dea_lp_stage2(prob2, info2, X, Y, x0, y0, ort, rts, theta_star)
        res2 = _rsm_two_phase(prob2, tol, trace)

        if (res2.status != 0) {
            // Stage 2 failed: return θ* with zero lambdas/slacks
            if (trace == "trace") {
                printf("{txt}DMU %g: Stage2 after fallback failed; " +
                       "returning theta only\n", k)
            }
            result = J(nc, 1, 0)
            if (ort == "output" && theta_star != 0) {
                result[1] = 1 / theta_star
            }
            else {
                result[1] = theta_star
            }
            return(result)
        }

        if (trace == "trace") {
            printf("{txt}DMU %g: Stage2 after fallback OK, " +
                   "sum_slacks = %12.8f (iter=%g)\n",
                   k, -res2.objval, res2.iters)
        }
        result = _dea_extract(res2, info2, theta_star, 2)
        return(result)
    }

    theta_star = _dea_extract_theta(res1, info1)

    if (trace == "trace") {
        printf("{txt}DMU %g Stage1: theta* = %12.8f (iter=%g)\n",
               k, theta_star, res1.iters)
    }

    // =================================================================
    // STAGE 2: Maximize slacks with theta fixed
    // =================================================================
    _dea_lp_stage2(prob2, info2, X, Y, x0, y0, ort, rts, theta_star)

    res2 = _rsm_two_phase(prob2, tol, trace)

    if (res2.status != 0) {
        if (trace == "trace") {
            printf("{err}DMU %g Stage2 FAIL: %s (falling back to Stage1)\n",
                   k, res2.msg)
        }
        result = _dea_extract(res1, info1, theta_star, 1)
    }
    else {
        if (trace == "trace") {
            printf("{txt}DMU %g Stage2: sum_slacks = %12.8f (iter=%g)\n",
                   k, -res2.objval, res2.iters)
        }
        result = _dea_extract(res2, info2, theta_star, 2)
    }

    return(result)
}


// ============================================================================
// BUILD LP: STAGE 1 (Optimize efficiency)
// ============================================================================
/**
 * _dea_lp_stage1()
 *
 * Builds the standard-form LP for DEA Stage 1.
 *
 * INPUT-ORIENTED (min theta):
 *   Variables: [theta, lambda_1..n, si_1..m, so_1..s, (sc), artif...]
 *   Constraints:
 *     Row i=1..m:   -x0_i*theta + SUM_j(X_ij*lam_j) + si_i       = 0
 *     Row m+r=1..s: SUM_j(Y_rj*lam_j)               - so_r + a_r = y0_r
 *     [VRS: SUM(lam)=1 +artif;  IRS: SUM(lam)-sc=1 +artif;  DRS: SUM(lam)+sc=1]
 *
 * OUTPUT-ORIENTED (min -phi, i.e., max phi):
 *   Variables: [phi, lambda_1..n, si_1..m, so_1..s, (sc), artif...]
 *   Constraints:
 *     Row i=1..m:   SUM_j(X_ij*lam_j) + si_i                     = x0_i
 *     Row m+r=1..s: -y0_r*phi + SUM_j(Y_rj*lam_j) - so_r  + a_r = 0
 *     [convexity constraint as above]
 *
 * Artificials are placed where natural slack has wrong sign (-1) or
 * where there is no slack (VRS equality, IRS surplus).
 * init_basic[i] records which variable starts as basic for row i.
 */
void _dea_lp_stage1(struct RSM_Prob scalar prob,
                     struct DEA_Info scalar info,
                     real matrix X,
                     real matrix Y,
                     real colvector x0,
                     real colvector y0,
                     string scalar ort,
                     string scalar rts)
{
    real scalar n_dmu, n_inp, n_out, nvar, ncon, nart, has_conv
    real scalar i, acol
    real matrix A
    real colvector b, ib
    real rowvector c, is_art

    n_dmu = rows(X)
    n_inp = cols(X)
    n_out = cols(Y)
    has_conv = (rts == "vrs" || rts == "irs" || rts == "drs")

    ncon = n_inp + n_out + (has_conv ? 1 : 0)
    nvar = 1 + n_dmu + n_inp + n_out + (rts == "irs" || rts == "drs" ? 1 : 0)

    // Artificial count: output constraints always need them (so has -1 coeff)
    nart = n_out
    if (rts == "vrs" || rts == "irs") nart = nart + 1

    // Store DEA index mapping
    info.n_dmu = n_dmu
    info.n_inp = n_inp
    info.n_out = n_out
    info.ort   = ort
    info.rts   = rts
    info.idx_theta = 1
    info.idx_lam   = (2..n_dmu + 1)
    info.idx_si    = (n_dmu + 2..n_dmu + n_inp + 1)
    info.idx_so    = (n_dmu + n_inp + 2..n_dmu + n_inp + n_out + 1)
    info.idx_sc    = (rts == "irs" || rts == "drs" ? n_dmu + n_inp + n_out + 2 : 0)

    // Allocate LP
    prob.m_con  = ncon
    prob.n_orig = nvar
    prob.n_var  = nvar + nart
    A      = J(ncon, prob.n_var, 0)
    b      = J(ncon, 1, 0)
    c      = J(1, prob.n_var, 0)
    is_art = J(1, prob.n_var, 0)
    ib     = J(ncon, 1, 0)

    // ---- Fill constraint matrix ----
    if (ort == "input") {
        // Input constraints: -x0_i*theta + X'_i*lam + si_i = 0
        for (i = 1; i <= n_inp; i++) {
            A[i, 1]                = -x0[i]
            A[i, info.idx_lam]     = X[., i]'
            A[i, info.idx_si[i]]   = 1
            b[i] = 0
            ib[i] = info.idx_si[i]                   // slack basic at 0
        }
        // Output constraints: Y'_r*lam - so_r + a_r = y0_r
        for (i = 1; i <= n_out; i++) {
            A[n_inp + i, info.idx_lam]     = Y[., i]'
            A[n_inp + i, info.idx_so[i]]   = -1
            b[n_inp + i] = y0[i]
            acol = nvar + i
            A[n_inp + i, acol] = 1
            is_art[acol] = 1
            ib[n_inp + i] = acol                     // artificial basic at y0_r
        }
        // Objective: min theta
        c[1] = 1

    }
    else {
        // OUTPUT-ORIENTED: min(-phi) to maximize phi
        // Input constraints: X'_i*lam + si_i = x0_i
        for (i = 1; i <= n_inp; i++) {
            A[i, info.idx_lam]     = X[., i]'
            A[i, info.idx_si[i]]   = 1
            b[i] = x0[i]
            ib[i] = info.idx_si[i]
        }
        // Output constraints: -y0_r*phi + Y'_r*lam - so_r + a_r = 0
        for (i = 1; i <= n_out; i++) {
            A[n_inp + i, 1]                = -y0[i]
            A[n_inp + i, info.idx_lam]     = Y[., i]'
            A[n_inp + i, info.idx_so[i]]   = -1
            b[n_inp + i] = 0
            acol = nvar + i
            A[n_inp + i, acol] = 1
            is_art[acol] = 1
            ib[n_inp + i] = acol
        }
        // Objective: min -phi
        c[1] = -1
    }

    // ---- Convexity constraint ----
    if (has_conv) {
        i = n_inp + n_out + 1
        A[i, info.idx_lam] = J(1, n_dmu, 1)

        if (rts == "vrs") {
            b[i] = 1
            acol = nvar + n_out + 1
            A[i, acol] = 1
            is_art[acol] = 1
            ib[i] = acol
        }
        else if (rts == "irs") {
            A[i, info.idx_sc] = -1
            b[i] = 1
            acol = nvar + n_out + 1
            A[i, acol] = 1
            is_art[acol] = 1
            ib[i] = acol
        }
        else if (rts == "drs") {
            A[i, info.idx_sc] = 1
            b[i] = 1
            ib[i] = info.idx_sc
        }
    }

    // Ensure b >= 0
    for (i = 1; i <= ncon; i++) {
        if (b[i] < 0) {
            A[i, .] = -A[i, .]
            b[i] = -b[i]
        }
    }

    prob.A = A
    prob.b = b
    prob.c = c
    prob.is_artif = is_art
    prob.init_basic = ib
}


// ============================================================================
// BUILD LP: STAGE 2 (Maximize slacks, theta fixed)
// ============================================================================
/**
 * _dea_lp_stage2()
 *
 * min -(SUM si + SUM so)   [maximize total slacks]
 * with theta = theta* substituted into RHS.
 * theta is NOT a variable in Stage 2.
 * Variables: [lam_1..n, si_1..m, so_1..s, (sc)]
 */
void _dea_lp_stage2(struct RSM_Prob scalar prob,
                     struct DEA_Info scalar info,
                     real matrix X,
                     real matrix Y,
                     real colvector x0,
                     real colvector y0,
                     string scalar ort,
                     string scalar rts,
                     real scalar theta_star)
{
    real scalar n_dmu, n_inp, n_out, nvar, ncon, nart, has_conv
    real scalar i, acol
    real matrix A
    real colvector b, ib
    real rowvector c, is_art

    n_dmu = rows(X)
    n_inp = cols(X)
    n_out = cols(Y)
    has_conv = (rts == "vrs" || rts == "irs" || rts == "drs")
    ncon = n_inp + n_out + (has_conv ? 1 : 0)

    // No theta variable
    nvar = n_dmu + n_inp + n_out + (rts == "irs" || rts == "drs" ? 1 : 0)
    nart = n_out
    if (rts == "vrs" || rts == "irs") nart = nart + 1

    // Set metadata (was missing → caused J():3300 in _dea_extract)
    info.n_dmu = n_dmu
    info.n_inp = n_inp
    info.n_out = n_out
    info.ort   = ort
    info.rts   = rts

    // Update index mapping (no theta → indices shift by -1)
    info.idx_theta = 0
    info.idx_lam   = (1..n_dmu)
    info.idx_si    = (n_dmu + 1..n_dmu + n_inp)
    info.idx_so    = (n_dmu + n_inp + 1..n_dmu + n_inp + n_out)
    info.idx_sc    = (rts == "irs" || rts == "drs" ? n_dmu + n_inp + n_out + 1 : 0)

    prob.m_con  = ncon
    prob.n_orig = nvar
    prob.n_var  = nvar + nart
    A      = J(ncon, prob.n_var, 0)
    b      = J(ncon, 1, 0)
    c      = J(1, prob.n_var, 0)
    is_art = J(1, prob.n_var, 0)
    ib     = J(ncon, 1, 0)

    if (ort == "input") {
        // Input: X'_i*lam + si_i = theta* * x0_i
        for (i = 1; i <= n_inp; i++) {
            A[i, info.idx_lam]     = X[., i]'
            A[i, info.idx_si[i]]   = 1
            b[i] = theta_star * x0[i]
            ib[i] = info.idx_si[i]
        }
        // Output: Y'_r*lam - so_r + a_r = y0_r
        for (i = 1; i <= n_out; i++) {
            A[n_inp + i, info.idx_lam]     = Y[., i]'
            A[n_inp + i, info.idx_so[i]]   = -1
            b[n_inp + i] = y0[i]
            acol = nvar + i
            A[n_inp + i, acol] = 1
            is_art[acol] = 1
            ib[n_inp + i] = acol
        }
    }
    else {
        // Output-oriented, phi fixed:
        // Input: X'_i*lam + si_i = x0_i
        for (i = 1; i <= n_inp; i++) {
            A[i, info.idx_lam]     = X[., i]'
            A[i, info.idx_si[i]]   = 1
            b[i] = x0[i]
            ib[i] = info.idx_si[i]
        }
        // Output: Y'_r*lam - so_r + a_r = theta_star * y0_r
        for (i = 1; i <= n_out; i++) {
            A[n_inp + i, info.idx_lam]     = Y[., i]'
            A[n_inp + i, info.idx_so[i]]   = -1
            b[n_inp + i] = theta_star * y0[i]
            acol = nvar + i
            A[n_inp + i, acol] = 1
            is_art[acol] = 1
            ib[n_inp + i] = acol
        }
    }

    // Convexity constraint
    if (has_conv) {
        i = n_inp + n_out + 1
        A[i, info.idx_lam] = J(1, n_dmu, 1)
        if (rts == "vrs") {
            b[i] = 1
            acol = nvar + n_out + 1
            A[i, acol] = 1
            is_art[acol] = 1
            ib[i] = acol
        }
        else if (rts == "irs") {
            A[i, info.idx_sc] = -1
            b[i] = 1
            acol = nvar + n_out + 1
            A[i, acol] = 1
            is_art[acol] = 1
            ib[i] = acol
        }
        else if (rts == "drs") {
            A[i, info.idx_sc] = 1
            b[i] = 1
            ib[i] = info.idx_sc
        }
    }

    // Objective: min -(SUM si + SUM so)
    c[info.idx_si] = J(1, n_inp, -1)
    c[info.idx_so] = J(1, n_out, -1)

    for (i = 1; i <= ncon; i++) {
        if (b[i] < 0) {
            A[i, .] = -A[i, .]
            b[i] = -b[i]
        }
    }

    prob.A = A
    prob.b = b
    prob.c = c
    prob.is_artif = is_art
    prob.init_basic = ib
}


// ============================================================================
// TWO-PHASE REVISED SIMPLEX METHOD
// ============================================================================
/**
 * _rsm_two_phase()
 *
 * Phase 1: min SUM(artificials) → find initial BFS
 * Phase 2: Optimize original objective with artificial columns removed
 *
 * Reference: Chen et al. (2010), Section 9.3 & 9.8:
 *   "The two-phase procedure is the only method utilized in practice."
 */
struct RSM_Result scalar _rsm_two_phase(struct RSM_Prob scalar prob,
                                         real scalar tol,
                                         string scalar trace)
{
    real scalar i, max_iter, p1_iter
    real colvector xb
    real rowvector c_orig
    real scalar art_sum
    real matrix B0

    struct RSM_Basis scalar basis
    struct RSM_Result scalar result

    max_iter = max((10000, 20 * prob.m_con))

    // Save Phase 2 objective
    c_orig = prob.c

    // =================================================================
    // PHASE 1: Minimize sum of artificials
    // =================================================================
    prob.c = prob.is_artif

    // Initialize basis
    basis.basic = prob.init_basic
    basis.n_upd = 0

    // Initial B^{-1} (should be identity by construction)
    B0 = prob.A[., basis.basic]
    if (mreldif(B0, I(prob.m_con)) < _dea_tol_feas()) {
        basis.Binv = I(prob.m_con)
    }
    else {
        basis.Binv = _rsm_inv(B0)
        if (hasmissing(basis.Binv)) {
            result.status = 1
            result.msg = "Singular initial basis"
            result.objval = .
            result.x = J(prob.n_var, 1, .)
            result.iters = 0
            prob.c = c_orig
            return(result)
        }
    }

    if (trace == "trace") {
        printf("{txt}Phase1: %g con, %g var, %g artif\n",
               prob.m_con, prob.n_var, sum(prob.is_artif))
    }

    // Solve Phase 1
    result = _rsm_solve(prob, basis, tol, max_iter, 1, trace)
    p1_iter = result.iters

    if (result.status == 2) {
        result.status = 1
        result.msg = "Phase1 unbounded (bug)"
        prob.c = c_orig
        return(result)
    }

    // Check artificial sum
    xb = basis.Binv * prob.b
    art_sum = 0
    for (i = 1; i <= prob.m_con; i++) {
        if (basis.basic[i] > prob.n_orig) {
            art_sum = art_sum + max((0, xb[i]))
        }
    }

    if (trace == "trace") {
        printf("{txt}Phase1 done: art_sum=%e iter=%g\n", art_sum, p1_iter)
    }

    if (art_sum > _dea_tol_feas() * 100) {
        result.status = 1
        result.msg = sprintf("Infeasible (art_sum=%e)", art_sum)
        prob.c = c_orig
        return(result)
    }

    // =================================================================
    // TRANSITION: Purge artificials, drop columns, restore objective
    // =================================================================
    _rsm_purge_artif(prob, basis, tol)
    _rsm_drop_artif(prob, basis)
    prob.c = c_orig[1..prob.n_var]

    // Refactorize for clean Phase 2 start
    _rsm_refac(basis, prob)

    // =================================================================
    // PHASE 2: Optimize original objective
    // =================================================================
    if (trace == "trace") {
        printf("{txt}Phase2: %g con, %g var\n", prob.m_con, prob.n_var)
    }

    result = _rsm_solve(prob, basis, tol, max_iter - p1_iter, 2, trace)
    result.iters = result.iters + p1_iter

    // Extract full solution
    if (result.status == 0) {
        result.x = J(prob.n_var, 1, 0)
        xb = basis.Binv * prob.b
        for (i = 1; i <= prob.m_con; i++) {
            if (abs(xb[i]) > _dea_tol_zero()) {
                result.x[basis.basic[i]] = max((0, xb[i]))
            }
        }
        result.objval = prob.c * result.x
    }

    return(result)
}


// ============================================================================
// CORE REVISED SIMPLEX LOOP
// ============================================================================
/**
 * _rsm_solve()
 *
 * Implements the revised simplex algorithm per Chen et al. (2010) Sec 9.8:
 *
 *   Step 0 (Init): B^{-1}=I, u^T = c_B^T * B^{-1}
 *   Step 1 (Pivot Col): a_bar_j = B^{-1}*a_j,  c_bar_j = u^T*a_bar_j - c_j
 *   Step 2 (Optimality): c_bar >= 0 for all j in N → optimal
 *   Step 3 (Pivot Row): r = argmin{ b_bar_i / a_bar_ik : a_bar_ik > 0 }
 *   Step 4 (Pivoting): Update B^{-1} via eta-vector
 *
 * Degeneracy handling: consecutive degenerate pivots → Bland's rule
 */
struct RSM_Result scalar _rsm_solve(struct RSM_Prob scalar prob,
                                     struct RSM_Basis scalar basis,
                                     real scalar tol,
                                     real scalar max_iter,
                                     real scalar phase,
                                     string scalar trace)
{
    struct RSM_Result scalar result
    real scalar iter, entering, pivot_row
    real scalar use_bland, consec_degen, stall_count
    real scalar obj_val, prev_obj
    real colvector xb, d
    real scalar i

    result.status = 3
    result.objval = .
    result.x = J(prob.n_var, 1, .)
    result.iters = 0
    result.msg = "Maximum iterations"

    use_bland = 0
    consec_degen = 0
    stall_count = 0
    prev_obj = .

    for (iter = 1; iter <= max_iter; iter++) {
        result.iters = iter

        // Periodic refactorization
        if (basis.n_upd > 0 && mod(basis.n_upd, _dea_refac_freq()) == 0) {
            _rsm_refac(basis, prob)
        }

        // Current BFS: xb = B^{-1} * b
        xb = basis.Binv * prob.b
        xb = xb :* (abs(xb) :> _dea_tol_zero())

        // Objective: z = c_B' * xb
        obj_val = 0
        for (i = 1; i <= prob.m_con; i++) {
            obj_val = obj_val + prob.c[basis.basic[i]] * xb[i]
        }

        // Degeneracy / stall detection
        if (prev_obj < . && abs(obj_val - prev_obj) < _dea_tol_zero()) {
            stall_count = stall_count + 1
            if (stall_count > prob.m_con && !use_bland) {
                use_bland = 1
                if (trace == "trace") {
                    printf("{txt}  iter%g: Bland ON (stall=%g)\n", iter, stall_count)
                }
            }
        }
        else {
            stall_count = 0
        }
        prev_obj = obj_val

        // =============================================================
        // PRICING: select entering variable (Steps 1+2 of Chen et al.)
        // =============================================================
        entering = _rsm_enter(prob, basis, use_bland, tol)

        if (entering == 0) {
            // Optimal
            result.status = 0
            result.objval = obj_val
            result.msg = sprintf("Optimal (P%g iter%g)", phase, iter)
            result.x = J(prob.n_var, 1, 0)
            for (i = 1; i <= prob.m_con; i++) {
                result.x[basis.basic[i]] = max((0, xb[i]))
            }
            if (trace == "trace") {
                printf("{txt}  P%g optimal iter%g: obj=%12.8e\n",
                       phase, iter, obj_val)
            }
            return(result)
        }

        // Pivot column: d = B^{-1} * a_entering  (Step 1 cont.)
        d = basis.Binv * prob.A[., entering]

        // =============================================================
        // RATIO TEST: select pivot row (Step 3)
        // =============================================================
        pivot_row = _rsm_leave(xb, d, basis, use_bland, tol)

        if (pivot_row == 0) {
            result.status = 2
            result.objval = obj_val
            result.msg = sprintf("Unbounded (P%g iter%g)", phase, iter)
            return(result)
        }

        // Track degeneracy
        if (abs(xb[pivot_row]) < _dea_tol_feas()) {
            consec_degen = consec_degen + 1
            if (consec_degen > 2 * prob.m_con && !use_bland) {
                use_bland = 1
                if (trace == "trace") {
                    printf("{txt}  iter%g: Bland ON (degen=%g)\n",
                           iter, consec_degen)
                }
            }
        }
        else {
            consec_degen = 0
        }

        // =============================================================
        // PIVOT: update basis (Step 4)
        // =============================================================
        _rsm_pivot(basis, d, pivot_row)
        basis.basic[pivot_row] = entering

        // Progress trace
        if (trace == "trace" && mod(iter, 500) == 0) {
            printf("{txt}  P%g iter%g: obj=%12.8e deg=%g\n",
                   phase, iter, obj_val, consec_degen)
        }
    }

    return(result)
}


// ============================================================================
// PRICING: Select entering variable
// ============================================================================
/**
 * _rsm_enter()
 *
 * Computes reduced costs and selects entering variable.
 *
 * Key vectorization (improvement over ldea.do element-wise loop):
 *   y = c_B * B^{-1}                     (simplex multiplier, 1xm)
 *   rc_all = c_N - y * N                 (all non-basic reduced costs, 1xn_nb)
 *
 * For minimization, entering has most negative rc (Dantzig) or first
 * negative rc by index (Bland).
 *
 * @return Column index of entering variable (0 if optimal)
 */
real scalar _rsm_enter(struct RSM_Prob scalar prob,
                        struct RSM_Basis scalar basis,
                        real scalar use_bland,
                        real scalar tol)
{
    real rowvector y, c_B, rc_all
    real scalar i, j, n_nb, rc_j, best_rc, opt_tol
    real colvector nb, is_basic

    opt_tol = _dea_tol_opt()

    // Simplex multiplier: y = c_B * B^{-1}  [VECTORIZED]
    c_B = J(1, prob.m_con, 0)
    for (i = 1; i <= prob.m_con; i++) {
        c_B[i] = prob.c[basis.basic[i]]
    }
    y = c_B * basis.Binv

    // Build non-basic index set
    is_basic = J(prob.n_var, 1, 0)
    for (i = 1; i <= prob.m_con; i++) {
        is_basic[basis.basic[i]] = 1
    }
    nb = select((1::prob.n_var), !is_basic)
    n_nb = rows(nb)
    if (n_nb == 0) return(0)

    if (use_bland) {
        // BLAND'S RULE: first negative rc by original variable index
        nb = sort(nb, 1)
        for (i = 1; i <= n_nb; i++) {
            j = nb[i]
            rc_j = prob.c[j] - y * prob.A[., j]
            if (rc_j < -opt_tol) return(j)
        }
        return(0)
    }
    else {
        // DANTZIG'S RULE: most negative rc [VECTORIZED]
        rc_all = prob.c[nb'] - y * prob.A[., nb]    // 1 x n_nb
        best_rc = min(rc_all)
        if (best_rc >= -opt_tol) return(0)

        for (i = 1; i <= n_nb; i++) {
            if (abs(rc_all[i] - best_rc) < _dea_tol_zero()) {
                return(nb[i])
            }
        }
        return(0)
    }
}


// ============================================================================
// RATIO TEST: Select leaving variable (pivot row)
// ============================================================================
/**
 * _rsm_leave()
 *
 * Harris two-pass ratio test:
 *   Pass 1: min_ratio = min{ xb_i / d_i : d_i > tol_pivot }
 *   Pass 2: Among rows with ratio <= min_ratio + harris_tol,
 *           pick largest |d_i| (stable) or smallest basis index (Bland)
 *
 * @return Pivot row (1..m), 0 if unbounded
 */
real scalar _rsm_leave(real colvector xb,
                        real colvector d,
                        struct RSM_Basis scalar basis,
                        real scalar use_bland,
                        real scalar tol)
{
    real scalar m, i, pivot_row
    real scalar min_ratio, harris_bound
    real scalar best_d, min_idx
    real colvector ratios
    real scalar pvt_tol

    m = rows(d)
    pvt_tol = _dea_tol_pivot()

    if (max(d) <= pvt_tol) return(0)

    // Pass 1: compute ratios
    ratios = J(m, 1, .)
    for (i = 1; i <= m; i++) {
        if (d[i] > pvt_tol) {
            if (xb[i] >= -_dea_tol_feas()) {
                ratios[i] = max((0, xb[i])) / d[i]
            }
        }
    }

    min_ratio = min(ratios)
    if (min_ratio == .) return(0)

    // Pass 2: select among eligible
    harris_bound = min_ratio + _dea_harris_tol()

    if (use_bland) {
        pivot_row = 0
        min_idx = .
        for (i = 1; i <= m; i++) {
            if (ratios[i] < . && ratios[i] <= harris_bound) {
                if (basis.basic[i] < min_idx) {
                    min_idx = basis.basic[i]
                    pivot_row = i
                }
            }
        }
    }
    else {
        pivot_row = 0
        best_d = 0
        for (i = 1; i <= m; i++) {
            if (ratios[i] < . && ratios[i] <= harris_bound) {
                if (d[i] > best_d) {
                    best_d = d[i]
                    pivot_row = i
                }
            }
        }
    }

    return(pivot_row)
}


// ============================================================================
// PIVOT: Eta-vector basis update
// ============================================================================
/**
 * _rsm_pivot()
 *
 * Updates B^{-1} via eta-vector (product form of inverse):
 *
 *   For i != p:  B_inv[i,.] -= (d[i]/d[p]) * B_inv[p,.]
 *   For i =  p:  B_inv[p,.] /= d[p]
 *
 * This is O(m^2) per pivot vs O(m^3) for full refactorization.
 * ldea.do recomputes full lusolve(B, I(m)) at EVERY iteration,
 * which is O(m^3) each time -- the eta update is a key improvement.
 */
void _rsm_pivot(struct RSM_Basis scalar basis,
                 real colvector d,
                 real scalar pivot_row)
{
    real scalar m, i, pv
    real rowvector prow

    m = rows(d)
    pv = d[pivot_row]

    if (abs(pv) < _dea_tol_pivot()) return

    // Save pivot row before modification
    prow = basis.Binv[pivot_row, .]

    // Update non-pivot rows
    for (i = 1; i <= m; i++) {
        if (i != pivot_row && abs(d[i]) > _dea_tol_zero()) {
            basis.Binv[i, .] = basis.Binv[i, .] - (d[i] / pv) * prow
        }
    }

    // Scale pivot row
    basis.Binv[pivot_row, .] = prow / pv

    basis.n_upd = basis.n_upd + 1
}


// ============================================================================
// REFACTORIZE: Recompute B^{-1} from scratch
// ============================================================================
/**
 * Corrects accumulated round-off from eta updates.
 * Strategy: LU first, SVD fallback if ill-conditioned.
 */
void _rsm_refac(struct RSM_Basis scalar basis,
                 struct RSM_Prob scalar prob)
{
    real matrix B, Binv_new

    B = prob.A[., basis.basic]
    Binv_new = _rsm_inv(B)

    if (!hasmissing(Binv_new)) {
        basis.Binv = Binv_new
        basis.n_upd = 0
    }
}


// ============================================================================
// SAFE MATRIX INVERSE
// ============================================================================
/**
 * LU → SVD fallback strategy for inverting B.
 */
real matrix _rsm_inv(real matrix B)
{
    real scalar m, cnd
    real matrix Binv

    m = rows(B)
    if (m == 0) return(J(0, 0, .))

    if (mreldif(B, I(m)) < _dea_tol_zero()) return(I(m))

    cnd = cond(B)

    if (cnd < 1e10) {
        Binv = lusolve(B, I(m))
        if (!hasmissing(Binv)) return(Binv)
    }

    Binv = svsolve(B, I(m), 1e-14)
    if (!hasmissing(Binv)) return(Binv)

    return(J(m, m, .))
}


// ============================================================================
// PURGE ARTIFICIALS FROM BASIS
// ============================================================================
/**
 * After Phase 1, artificial vars may remain in basis at value 0.
 * Try to pivot each out using a non-artificial non-basic variable.
 */
void _rsm_purge_artif(struct RSM_Prob scalar prob,
                       struct RSM_Basis scalar basis,
                       real scalar tol)
{
    real scalar i, j, ii, found_basic
    real colvector d
    real scalar best_j, best_val, val

    for (i = 1; i <= prob.m_con; i++) {
        if (basis.basic[i] <= prob.n_orig) continue

        best_j = 0
        best_val = 0

        for (j = 1; j <= prob.n_orig; j++) {
            // Check if j is already basic
            found_basic = 0
            for (ii = 1; ii <= prob.m_con; ii++) {
                if (basis.basic[ii] == j) {
                    found_basic = 1
                    break
                }
            }
            if (found_basic) continue

            // Check if j can pivot into row i
            d = basis.Binv * prob.A[., j]
            val = abs(d[i])
            if (val > _dea_tol_pivot() && val > best_val) {
                best_val = val
                best_j = j
            }
        }

        if (best_j > 0) {
            d = basis.Binv * prob.A[., best_j]
            _rsm_pivot(basis, d, i)
            basis.basic[i] = best_j
        }
    }
}


// ============================================================================
// DROP ARTIFICIAL COLUMNS
// ============================================================================
/**
 * Removes artificial columns from problem after Phase 1.
 */
void _rsm_drop_artif(struct RSM_Prob scalar prob,
                      struct RSM_Basis scalar basis)
{
    real scalar new_n, i, j, found
    real colvector used

    if (prob.n_var == prob.n_orig) return

    new_n = prob.n_orig
    prob.A = prob.A[., 1..new_n]
    prob.c = prob.c[1..new_n]
    prob.is_artif = J(1, new_n, 0)
    prob.n_var = new_n

    // Fix any remaining artificial in basis (shouldn't happen after purge)
    for (i = 1; i <= prob.m_con; i++) {
        if (basis.basic[i] > new_n) {
            used = J(new_n, 1, 0)
            for (j = 1; j <= prob.m_con; j++) {
                if (basis.basic[j] <= new_n) {
                    used[basis.basic[j]] = 1
                }
            }
            found = 0
            for (j = 1; j <= new_n; j++) {
                if (!used[j]) {
                    basis.basic[i] = j
                    found = 1
                    break
                }
            }
            if (!found) basis.basic[i] = 1
        }
    }
}


// ============================================================================
// EXTRACT THETA FROM LP RESULT
// ============================================================================
real scalar _dea_extract_theta(struct RSM_Result scalar res,
                                struct DEA_Info scalar info)
{
    real scalar theta

    if (res.status != 0) return(.)

    if (info.ort == "input") {
        // Input-oriented: theta = x[idx_theta]
        if (info.idx_theta > 0 && info.idx_theta <= rows(res.x)) {
            theta = res.x[info.idx_theta]
        }
        else theta = .
    }
    else {
        // Output-oriented: objective was min(-phi), so phi = -objval
        theta = -res.objval
    }

    return(theta)
}


// ============================================================================
// EXTRACT FULL DEA RESULT
// ============================================================================
/**
 * Builds result vector [theta, lam_1..n, si_1..m, so_1..s, (sc)]
 * from LP solution.
 */
real colvector _dea_extract(struct RSM_Result scalar res,
                             struct DEA_Info scalar info,
                             real scalar theta_star,
                             real scalar stage)
{
    real scalar nc, i
    real colvector result

    nc = 1 + info.n_dmu + info.n_inp + info.n_out
    if (info.rts == "irs" || info.rts == "drs") nc = nc + 1
    result = J(nc, 1, 0)

    if (res.status != 0) return(J(nc, 1, .))

    // For output-oriented, theta_star = phi (expansion factor >= 1)
    // DEA convention: efficiency = 1/phi for output-oriented
    if (info.ort == "output" && theta_star != 0) {
        result[1] = 1 / theta_star
    }
    else {
        result[1] = theta_star
    }

    // Lambda
    for (i = 1; i <= info.n_dmu; i++) {
        if (info.idx_lam[i] <= rows(res.x)) {
            result[1 + i] = max((0, res.x[info.idx_lam[i]]))
        }
    }

    // Input slacks
    for (i = 1; i <= info.n_inp; i++) {
        if (info.idx_si[i] <= rows(res.x)) {
            result[1 + info.n_dmu + i] = max((0, res.x[info.idx_si[i]]))
        }
    }

    // Output slacks
    for (i = 1; i <= info.n_out; i++) {
        if (info.idx_so[i] <= rows(res.x)) {
            result[1 + info.n_dmu + info.n_inp + i] = max((0, res.x[info.idx_so[i]]))
        }
    }

    // Convexity slack
    if (info.idx_sc > 0 && info.idx_sc <= rows(res.x)) {
        result[nc] = max((0, res.x[info.idx_sc]))
    }

    // Clean near-zeros
    for (i = 1; i <= nc; i++) {
        if (abs(result[i]) < _dea_tol_feas()) result[i] = 0
    }

    return(result)
}



// ############################################################################
// #                                                                          #
// #               MULTIPLIER (DUAL) MODEL                                    #
// #                                                                          #
// ############################################################################
//
// The multiplier model is the LP dual of the envelopment model.
// It determines optimal input/output weights (v*, u*) directly.
//
// INPUT-ORIENTED (all RTS):
//   max  u'y₀ + u₀          (u₀=0 for CRS)
//   s.t. u'y_j - v'x_j + u₀ ≤ 0   for all j=1..n
//        v'x₀ = 1           (normalization)
//        v, u ≥ 0
//        u₀:  0 (CRS), free (VRS), ≤0 (IRS), ≥0 (DRS)
//
// OUTPUT-ORIENTED (all RTS):
//   min  v'x₀ + u₀          (u₀=0 for CRS)
//   s.t. v'x_j - u'y_j + u₀ ≥ 0   for all j=1..n
//        u'y₀ = 1           (normalization)
//        v, u ≥ 0
//        u₀:  0 (CRS), free (VRS), ≤0 (IRS), ≥0 (DRS)
//
// By strong duality:
//   Input:  θ* = u*'y₀ + u₀* = envelopment θ*
//   Output: φ* = v*'x₀ + u₀* = envelopment φ*
//
// Standard form conversion (both orientations):
//   DMU row j:  -v'x_j + u'y_j + [u₀ terms] + s_j = 0
//   (sign of u₀ terms differs between input/output orientation)


// ============================================================================
// MULTIPLIER ENTRY: efficiency scores only
// ============================================================================
/**
 * _dea_simplex_mult_envelopment()
 *
 * Multiplier model counterpart of _dea_simplex_envelopment().
 * Solves for optimal weights and stores efficiency score in eff_var.
 */
void _dea_simplex_mult_envelopment(string scalar inputs,
                                    string scalar outputs,
                                    string scalar eff_var,
                                    string scalar touse,
                                    string scalar orientation,
                                    string scalar rts,
                                    real scalar tol)
{
    real matrix X, Y
    real colvector obs_indices, result_k
    real scalar n_dmu, k, obs_k

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))
    n_dmu = rows(X)

    for (k = 1; k <= n_dmu; k++) {
        obs_k = obs_indices[k]
        result_k = _dea_mult_solve_dmu(X, Y, k, orientation, rts, tol, "")
        st_store(obs_k, st_varindex(eff_var), result_k[1])
    }
}


// ============================================================================
// MULTIPLIER ENTRY: full results (efficiency + weights)
// ============================================================================
/**
 * _dea_simplex_mult_full()
 *
 * Stores results in r(dea_weights): n × (1 + m + s + 1) matrix.
 * Row format: [efficiency, v_1..v_m, u_1..u_s, u₀]
 */
void _dea_simplex_mult_full(string scalar inputs,
                             string scalar outputs,
                             string scalar touse,
                             string scalar orientation,
                             string scalar rts,
                             real scalar tol,
                             | string scalar trace)
{
    real matrix X, Y, RESULTS
    real colvector result_k
    real scalar n_dmu, n_inp, n_out, nc, k, nr_k
    string scalar do_trace

    do_trace = (args() >= 7 ? trace : "")

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    n_dmu = rows(X)
    n_inp = cols(X)
    n_out = cols(Y)

    // [efficiency, v_1..v_m, u_1..u_s, u₀]
    nc = 1 + n_inp + n_out + 1
    RESULTS = J(n_dmu, nc, 0)

    for (k = 1; k <= n_dmu; k++) {
        if (do_trace == "trace") {
            printf("\n{txt}=== DMU %g / %g (Multiplier) ===\n", k, n_dmu)
        }
        result_k = _dea_mult_solve_dmu(X, Y, k, orientation, rts, tol, do_trace)
        nr_k = rows(result_k)
        if (nr_k >= nc) {
            RESULTS[k, .] = result_k[1..nc]'
        }
        else if (nr_k > 0) {
            RESULTS[k, 1..nr_k] = result_k'
        }
    }

    st_matrix("r(dea_weights)", RESULTS)
}


// ============================================================================
// SOLVE SINGLE DMU: Multiplier model (single-stage LP)
// ============================================================================
/**
 * _dea_mult_solve_dmu()
 *
 * Unlike envelopment (two-stage: θ then slacks), the multiplier model
 * is a single LP that directly yields efficiency and optimal weights.
 *
 * @return Column vector [efficiency, v_1..v_m, u_1..u_s, u₀]
 */
real colvector _dea_mult_solve_dmu(real matrix X,
                                    real matrix Y,
                                    real scalar k,
                                    string scalar ort,
                                    string scalar rts,
                                    real scalar tol,
                                    string scalar trace)
{
    real scalar n_inp, n_out, nc
    real colvector x0, y0, result

    struct RSM_Prob scalar prob
    struct MULT_Info scalar minfo
    struct RSM_Result scalar res

    n_inp = cols(X)
    n_out = cols(Y)
    x0 = X[k, .]'
    y0 = Y[k, .]'
    nc = 1 + n_inp + n_out + 1

    _dea_mult_build_lp(prob, minfo, X, Y, x0, y0, ort, rts)

    res = _rsm_two_phase(prob, tol, trace)

    if (res.status != 0) {
        if (trace == "trace") {
            printf("{err}DMU %g Multiplier FAIL: %s\n", k, res.msg)
        }
        return(J(nc, 1, .))
    }

    result = _dea_mult_extract(res, minfo)

    if (trace == "trace") {
        printf("{txt}DMU %g Mult: eff=%12.8f (iter=%g)\n",
               k, result[1], res.iters)
    }

    return(result)
}

/**
 * _dea_simplex_multiplier()
 *
 * Solves the multiplier form for all DMUs and optionally stores the
 * stage-2 optimal input/output weights in user-visible variables.
 */
void _dea_simplex_multiplier(string scalar inputs,
                              string scalar outputs,
                              string scalar eff_var,
                              string scalar v_var,
                              string scalar u_var,
                              string scalar u0_var,
                              string scalar touse,
                              string scalar orientation,
                              string scalar rts,
                              real scalar tol,
                              string scalar trace)
{
    real matrix X, Y
    real colvector obs_indices, result
    real scalar n, m, s, k, obs_k, j, has_weights
    string rowvector v_varlist, u_varlist

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))

    n = rows(X)
    m = cols(X)
    s = cols(Y)

    has_weights = (v_var != "" & u_var != "")
    if (has_weights) {
        v_varlist = tokens(v_var)
        u_varlist = tokens(u_var)
    }

    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        result = _dea_mult_solve_dmu(X, Y, k, orientation, rts, tol, trace)

        st_store(obs_k, st_varindex(eff_var), result[1])

        if (has_weights) {
            for (j = 1; j <= m; j++) {
                st_store(obs_k, st_varindex(v_varlist[j]), result[1+j])
            }
            for (j = 1; j <= s; j++) {
                st_store(obs_k, st_varindex(u_varlist[j]), result[1+m+j])
            }
            if (u0_var != "") {
                st_store(obs_k, st_varindex(u0_var), result[2+m+s])
            }
        }
    }
}


// ============================================================================
// BUILD LP: Multiplier model
// ============================================================================
/**
 * _dea_mult_build_lp()
 *
 * Constructs the standard-form LP for the multiplier (dual) DEA model.
 *
 * Variables: [v(m), u(s), u0+(0/1), u0-(0/1), s(n)]
 * Constraints: n+1  (n DMU + 1 normalization)
 * Artificials: 1    (for normalization equality)
 *
 * The u₀ free variable is handled via splitting: u₀ = u₀⁺ - u₀⁻
 *   CRS:  no u₀ variable
 *   VRS:  u₀ free  → both u₀⁺, u₀⁻ present
 *   IRS:  u₀ ≤ 0   → only u₀⁻ (u₀ = -u₀⁻)
 *   DRS:  u₀ ≥ 0   → only u₀⁺ (u₀ = u₀⁺)
 */
void _dea_mult_build_lp(struct RSM_Prob scalar prob,
                          struct MULT_Info scalar minfo,
                          real matrix X,
                          real matrix Y,
                          real colvector x0,
                          real colvector y0,
                          string scalar ort,
                          string scalar rts)
{
    real scalar n_dmu, n_inp, n_out, nvar, ncon
    real scalar has_u0p, has_u0m, col_off
    real scalar i, j, acol
    real matrix A
    real colvector b, ib
    real rowvector c, is_art

    n_dmu = rows(X)
    n_inp = cols(X)
    n_out = cols(Y)

    // Metadata
    minfo.n_dmu = n_dmu
    minfo.n_inp = n_inp
    minfo.n_out = n_out
    minfo.ort   = ort
    minfo.rts   = rts

    // u₀ split
    if (ort == "input") {
        /* The LP stores w=-u0 in the input inequalities. */
        has_u0p = (rts == "vrs" || rts == "irs")
        has_u0m = (rts == "vrs" || rts == "drs")
    }
    else {
        has_u0p = (rts == "vrs" || rts == "drs")
        has_u0m = (rts == "vrs" || rts == "irs")
    }

    // ---- Variable layout ----
    // [v_1..v_m | u_1..u_s | u0+(?) | u0-(?) | s_1..s_n]
    minfo.idx_v = (1..n_inp)
    minfo.idx_u = (n_inp + 1..n_inp + n_out)

    col_off = n_inp + n_out
    if (has_u0p) {
        col_off = col_off + 1
        minfo.idx_u0p = col_off
    }
    else {
        minfo.idx_u0p = 0
    }
    if (has_u0m) {
        col_off = col_off + 1
        minfo.idx_u0m = col_off
    }
    else {
        minfo.idx_u0m = 0
    }

    minfo.idx_s = (col_off + 1..col_off + n_dmu)
    nvar = col_off + n_dmu
    ncon = n_dmu + 1              // n DMU rows + 1 normalization

    // Allocate LP (1 artificial for normalization)
    prob.m_con  = ncon
    prob.n_orig = nvar
    prob.n_var  = nvar + 1        // +1 artificial
    A      = J(ncon, prob.n_var, 0)
    b      = J(ncon, 1, 0)
    c      = J(1, prob.n_var, 0)
    is_art = J(1, prob.n_var, 0)
    ib     = J(ncon, 1, 0)

    // ---- DMU constraints (rows 1..n) ----
    //
    // Both orientations share the same basic structure after
    // converting to ≤ 0 then to standard form (Ax = b, x ≥ 0):
    //
    //   -v'x_j + u'y_j + [u₀ sign terms] + s_j = 0
    //
    // The SIGN of u₀ in the constraint differs:
    //
    //   Input:  original  u'y_j - v'x_j + u₀ ≤ 0
    //           std form  -v'x_j + u'y_j + u0+ - u0- + s_j = 0
    //
    //   Output: original  v'x_j - u'y_j + u₀ ≥ 0
    //           negate     -v'x_j + u'y_j - u₀ ≤ 0
    //           std form  -v'x_j + u'y_j - u0+ + u0- + s_j = 0

    for (j = 1; j <= n_dmu; j++) {
        // v columns: -x_ij
        A[j, minfo.idx_v] = -X[j, .]
        // u columns: +y_rj
        A[j, minfo.idx_u] = Y[j, .]
        // u₀ terms
        if (ort == "input") {
            if (has_u0p) A[j, minfo.idx_u0p] =  1
            if (has_u0m) A[j, minfo.idx_u0m] = -1
        }
        else {
            if (has_u0p) A[j, minfo.idx_u0p] = -1
            if (has_u0m) A[j, minfo.idx_u0m] =  1
        }
        // slack
        A[j, minfo.idx_s[j]] = 1
        b[j] = 0
        ib[j] = minfo.idx_s[j]         // basic at 0 (degenerate start)
    }

    // ---- Normalization constraint (row n+1) ----
    //   Input:  v'x₀ = 1
    //   Output: u'y₀ = 1
    if (ort == "input") {
        A[ncon, minfo.idx_v] = x0'
    }
    else {
        A[ncon, minfo.idx_u] = y0'
    }
    b[ncon] = 1
    acol = nvar + 1
    A[ncon, acol] = 1
    is_art[acol] = 1
    ib[ncon] = acol                     // artificial basic at 1

    // ---- Objective ----
    //   Input:  min -(u'y₀ + u₀) = min(-Σu_r·y₀_r - u₀⁺ + u₀⁻)
    //   Output: min  (v'x₀ + u₀) = min( Σv_i·x₀_i + u₀⁺ - u₀⁻)
    if (ort == "input") {
        c[minfo.idx_u] = -y0'
        if (has_u0p) c[minfo.idx_u0p] = -1
        if (has_u0m) c[minfo.idx_u0m] =  1
    }
    else {
        c[minfo.idx_v] = x0'
        if (has_u0p) c[minfo.idx_u0p] =  1
        if (has_u0m) c[minfo.idx_u0m] = -1
    }

    prob.A = A
    prob.b = b
    prob.c = c
    prob.is_artif = is_art
    prob.init_basic = ib
}


// ============================================================================
// EXTRACT MULTIPLIER RESULTS
// ============================================================================
/**
 * _dea_mult_extract()
 *
 * Extracts efficiency and optimal weights from LP solution.
 *
 * @return [efficiency, v_1..v_m, u_1..u_s, u₀]
 *
 * Efficiency:
 *   Input:  θ* = -(objective) = -objval  (∈ (0,1])
 *   Output: efficiency = 1/φ* = 1/objval (∈ (0,1])
 */
real colvector _dea_mult_extract(struct RSM_Result scalar res,
                                  struct MULT_Info scalar minfo)
{
    real scalar nc, i, u0_val, eff
    real colvector result

    nc = 1 + minfo.n_inp + minfo.n_out + 1
    result = J(nc, 1, 0)

    if (res.status != 0) return(J(nc, 1, .))

    // Input weights v*
    for (i = 1; i <= minfo.n_inp; i++) {
        result[1 + i] = max((0, res.x[minfo.idx_v[i]]))
    }

    // Output weights u*
    for (i = 1; i <= minfo.n_out; i++) {
        result[1 + minfo.n_inp + i] = max((0, res.x[minfo.idx_u[i]]))
    }

    // Reconstruct u₀ = u₀⁺ - u₀⁻
    u0_val = 0
    if (minfo.idx_u0p > 0) u0_val = u0_val + max((0, res.x[minfo.idx_u0p]))
    if (minfo.idx_u0m > 0) u0_val = u0_val - max((0, res.x[minfo.idx_u0m]))
    if (minfo.ort == "input") u0_val = -u0_val
    result[nc] = u0_val

    // Efficiency
    if (minfo.ort == "input") {
        eff = -res.objval          // θ* = max(u'y₀+u₀) = -(min objective)
    }
    else {
        // φ* = min(v'x₀+u₀) = objval;  efficiency = 1/φ*
        if (abs(res.objval) > _dea_tol_zero()) {
            eff = 1 / res.objval
        }
        else {
            eff = .
        }
    }
    result[1] = eff

    // Clean near-zeros
    for (i = 1; i <= nc; i++) {
        if (abs(result[i]) < _dea_tol_feas()) result[i] = 0
    }

    return(result)
}


// ############################################################################
// #                                                                          #
// #               GENERAL-PURPOSE LP SOLVER INTERFACE                        #
// #                                                                          #
// ############################################################################
//
// These functions allow any LP in inequality/equality form to be solved
// by the Two-Phase Revised Simplex core.
//
// Usage pattern (same as IPM):
//   1. Build LP in natural form:  min c'x, A_ineq*x <= b_ineq, A_eq*x = b_eq, x >= 0
//   2. Call _rsm_std_form() to convert to standard form RSM_Prob
//   3. Call _rsm_two_phase() to solve
//   4. Extract solution from result.x[1..n_orig_vars]


// ============================================================================
// STANDARD FORM CONVERTER
// ============================================================================
/**
 * _rsm_std_form()
 *
 * Converts LP from inequality/equality form to standard form (Ax=b, x>=0).
 *
 * Input:  min c'x  s.t.  A_ineq*x <= b_ineq,  A_eq*x = b_eq,  x >= 0
 *
 * Output (RSM_Prob):
 *   Variables: [x(nv), s_ineq(ni), a_eq(ne)]
 *   Rows:      [inequality rows (with slacks), equality rows (with artificials)]
 *   n_orig   = nv + ni  (survive into Phase 2)
 *   n_var    = nv + ni + ne
 *
 * Any row with b[i]<0 is multiplied by -1 to ensure b>=0.
 *
 * @param A_iq    ni × nv  inequality constraint matrix (may be 0×nv)
 * @param b_iq    ni × 1   inequality RHS
 * @param A_eq    ne × nv  equality constraint matrix (may be 0×nv)
 * @param b_eq    ne × 1   equality RHS
 * @param c       1  × nv  objective coefficients
 */
void _rsm_std_form(struct RSM_Prob scalar prob,
                    real matrix A_iq,
                    real colvector b_iq,
                    real matrix A_eq,
                    real colvector b_eq,
                    real rowvector c)
{
    real scalar nv, ni, ne, ncon
    real scalar n_neg_iq, n_art, n_orig, n_total
    real scalar i, row, art_col
    real matrix A
    real colvector b, ib
    real rowvector c_std, is_art

    nv = cols(c)
    ni = rows(A_iq)
    ne = rows(A_eq)
    ncon = ni + ne

    if (ncon == 0) {
        prob.m_con = 0
        prob.n_orig = nv
        prob.n_var = nv
        prob.A = J(0, nv, 0)
        prob.b = J(0, 1, 0)
        prob.c = c
        prob.is_artif = J(1, nv, 0)
        prob.init_basic = J(0, 1, 0)
        return
    }

    // Count inequality rows with b < 0 (need artificials)
    n_neg_iq = 0
    for (i = 1; i <= ni; i++) {
        if (b_iq[i] < 0) n_neg_iq = n_neg_iq + 1
    }

    // Variable layout:
    //   [x(nv), slack_ineq(ni), artificial(n_neg_iq + ne)]
    //
    // - b >= 0 inequality: slack +1, basic at b (no artificial)
    // - b <  0 inequality: flip row, surplus -1, artificial +1, basic at |b|
    // - equality:          artificial +1, basic at |b|
    //
    n_art   = n_neg_iq + ne
    n_orig  = nv + ni               // survive into Phase 2
    n_total = n_orig + n_art

    A      = J(ncon, n_total, 0)
    b      = J(ncon, 1, 0)
    c_std  = J(1, n_total, 0)
    is_art = J(1, n_total, 0)
    ib     = J(ncon, 1, 0)

    c_std[1..nv] = c

    // Track next artificial column
    art_col = n_orig

    // ---- Inequality rows ----
    row = 0
    for (i = 1; i <= ni; i++) {
        row = row + 1
        A[row, 1..nv] = A_iq[i, .]
        b[row] = b_iq[i]

        if (b[row] >= 0) {
            // Standard: A_row * x + slack = b,  slack basic at b
            A[row, nv + i] = 1
            ib[row] = nv + i
        }
        else {
            // b < 0: flip entire original part → -A_row * x = -b > 0
            A[row, 1..nv] = -A[row, 1..nv]
            b[row] = -b[row]

            // Surplus variable (flipped slack): coefficient -1
            A[row, nv + i] = -1

            // Artificial variable: coefficient +1, basic at |b|
            art_col = art_col + 1
            A[row, art_col] = 1
            is_art[art_col] = 1
            ib[row] = art_col
        }
    }

    // ---- Equality rows ----
    for (i = 1; i <= ne; i++) {
        row = row + 1
        A[row, 1..nv] = A_eq[i, .]
        b[row] = b_eq[i]

        if (b[row] < 0) {
            A[row, 1..nv] = -A[row, 1..nv]
            b[row] = -b[row]
        }

        // Artificial variable
        art_col = art_col + 1
        A[row, art_col] = 1
        is_art[art_col] = 1
        ib[row] = art_col
    }

    prob.m_con     = ncon
    prob.n_orig    = n_orig
    prob.n_var     = n_total
    prob.A         = A
    prob.b         = b
    prob.c         = c_std
    prob.is_artif  = is_art
    prob.init_basic = ib
}


// ============================================================================
// GENERAL LP SOLVER
// ============================================================================
/**
 * _rsm_solve_lp()
 *
 * One-call solver: build standard form → two-phase RSM → return result.
 * Solution vector in result.x has length n_orig_vars (original x only).
 *
 * @return RSM_Result with status, objval, x[1..nv] (original variables only)
 */
struct RSM_Result scalar _rsm_solve_lp(real matrix A_iq,
                                        real colvector b_iq,
                                        real matrix A_eq,
                                        real colvector b_eq,
                                        real rowvector c,
                                        real scalar tol)
{
    struct RSM_Prob scalar prob
    struct RSM_Result scalar res
    real scalar nv

    nv = cols(c)

    _rsm_std_form(prob, A_iq, b_iq, A_eq, b_eq, c)
    res = _rsm_two_phase(prob, tol, "")

    // Trim solution to original variables only
    if (res.status == 0 && rows(res.x) >= nv) {
        res.x = res.x[1..nv]
    }
    else if (res.status == 0) {
        res.x = J(nv, 1, 0)
    }

    return(res)
}


// ############################################################################
// #                                                                          #
// #               EXTENDED DEA MODELS                                        #
// #                                                                          #
// ############################################################################


// ============================================================================
// COST EFFICIENCY MODEL
// ============================================================================
/**
 * _dea_simplex_cost()
 *
 * min w'x  s.t.  X'λ ≤ x,  Y'λ ≥ y₀,  [convexity],  λ,x ≥ 0
 * Cost Efficiency = min_cost / observed_cost
 *
 * Interface matches _dea_ipm_cost() exactly.
 */
void _dea_simplex_cost(string scalar inputs,
                        string scalar outputs,
                        string scalar eff_var,
                        string scalar w_var,
                        string scalar touse,
                        string scalar rts)
{
    real matrix X, Y, W
    real colvector obs_indices, x0, y0, w0
    real rowvector c
    real matrix A_iq, A_eq
    real colvector b_iq, b_eq
    real scalar n, m, s, k, obs_k
    real scalar min_cost, obs_cost, cost_eff

    struct RSM_Result scalar res

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    W = st_data(., tokens(w_var), touse)
    obs_indices = selectindex(st_data(., touse))

    n = rows(X); m = cols(X); s = cols(Y)

    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, .]'; y0 = Y[k, .]'; w0 = W[k, .]'
        cost_eff = .

        obs_cost = sum(w0 :* x0)

        // Variables: [λ(n), x(m)]
        // min w'x
        c = J(1, n, 0), w0'

        // X'λ - x ≤ 0;   -Y'λ ≤ -y0
        A_iq = (X', -I(m)) \ (-Y', J(s, m, 0))
        b_iq = J(m, 1, 0) \ (-y0)

        A_eq = J(0, n+m, 0)
        b_eq = J(0, 1, 0)

        if (rts == "vrs") {
            A_eq = (J(1, n, 1), J(1, m, 0))
            b_eq = 1
        }
        else if (rts == "drs") {
            A_iq = A_iq \ (J(1, n, 1), J(1, m, 0))
            b_iq = b_iq \ 1
        }
        else if (rts == "irs") {
            A_iq = A_iq \ (J(1, n, -1), J(1, m, 0))
            b_iq = b_iq \ (-1)
        }

        res = _rsm_solve_lp(A_iq, b_iq, A_eq, b_eq, c, 1e-8)

        if (res.status == 0) {
            min_cost = res.objval
            if (obs_cost > 0) cost_eff = min_cost / obs_cost
        }

        st_store(obs_k, st_varindex(eff_var), cost_eff)
    }
}


// ============================================================================
// REVENUE EFFICIENCY MODEL
// ============================================================================
/**
 * _dea_simplex_revenue()
 *
 * max p'y  s.t.  X'λ ≤ x₀,  Y'λ ≥ y,  [convexity],  λ,y ≥ 0
 * Revenue Efficiency = observed_rev / max_rev
 *
 * Interface matches _dea_ipm_revenue() exactly.
 */
void _dea_simplex_revenue(string scalar inputs,
                           string scalar outputs,
                           string scalar eff_var,
                           string scalar p_var,
                           string scalar touse,
                           string scalar rts)
{
    real matrix X, Y, P
    real colvector obs_indices, x0, y0, p0
    real rowvector c
    real matrix A_iq, A_eq
    real colvector b_iq, b_eq
    real scalar n, m, s, k, obs_k
    real scalar max_rev, obs_rev, rev_eff

    struct RSM_Result scalar res

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    P = st_data(., tokens(p_var), touse)
    obs_indices = selectindex(st_data(., touse))

    n = rows(X); m = cols(X); s = cols(Y)

    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, .]'; y0 = Y[k, .]'; p0 = P[k, .]'
        rev_eff = .

        obs_rev = sum(p0 :* y0)

        // Variables: [λ(n), y(s)]
        // max p'y → min -p'y
        c = J(1, n, 0), -p0'

        // X'λ ≤ x₀;   -Y'λ + y ≤ 0
        A_iq = (X', J(m, s, 0)) \ (-Y', I(s))
        b_iq = x0 \ J(s, 1, 0)

        A_eq = J(0, n+s, 0)
        b_eq = J(0, 1, 0)

        if (rts == "vrs") {
            A_eq = (J(1, n, 1), J(1, s, 0))
            b_eq = 1
        }
        else if (rts == "drs") {
            A_iq = A_iq \ (J(1, n, 1), J(1, s, 0))
            b_iq = b_iq \ 1
        }
        else if (rts == "irs") {
            A_iq = A_iq \ (J(1, n, -1), J(1, s, 0))
            b_iq = b_iq \ (-1)
        }

        res = _rsm_solve_lp(A_iq, b_iq, A_eq, b_eq, c, 1e-8)

        if (res.status == 0) {
            max_rev = -res.objval
            if (max_rev > 0) rev_eff = obs_rev / max_rev
        }

        st_store(obs_k, st_varindex(eff_var), rev_eff)
    }
}


// ============================================================================
// PROFIT EFFICIENCY MODEL
// ============================================================================
/**
 * _dea_simplex_profit()
 *
 * Profit efficiency (Cooper, Seiford & Tone, Ch.8 eq. 8.8).
 *
 * Variables: λ only (x = Xλ, y = Yλ determined by technology).
 * max (pY - cX)'λ  s.t.  X'λ ≤ x_o,  -Y'λ ≤ -y_o,  [convexity],  λ ≥ 0
 * Profit Efficiency = (p'y_o - c'x_o) / (p'y* - c'x*)
 */
void _dea_simplex_profit(string scalar inputs,
                          string scalar outputs,
                          string scalar eff_var,
                          string scalar w_var,
                          string scalar p_var,
                          string scalar touse,
                          string scalar rts)
{
    real matrix X, Y, W, P
    real colvector obs_indices, x0, y0, w0, p0
    real rowvector c
    real matrix A_iq, A_eq
    real colvector b_iq, b_eq
    real scalar n, m, s, k, obs_k
    real scalar max_profit, obs_profit, profit_eff

    struct RSM_Result scalar res

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    W = st_data(., tokens(w_var), touse)
    P = st_data(., tokens(p_var), touse)
    obs_indices = selectindex(st_data(., touse))

    n = rows(X); m = cols(X); s = cols(Y)

    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, .]'; y0 = Y[k, .]'
        w0 = W[k, .]'; p0 = P[k, .]'
        profit_eff = .

        obs_profit = sum(p0 :* y0) - sum(w0 :* x0)

        // Variables: λ(n)
        // min (Xw0 - Yp0)'  (negated max objective)
        c = (X * w0 - Y * p0)'

        // X'λ ≤ x0;   -Y'λ ≤ -y0
        A_iq = X' \ (-Y')
        b_iq = x0 \ (-y0)

        A_eq = J(0, n, 0)
        b_eq = J(0, 1, 0)

        if (rts == "vrs") {
            A_eq = J(1, n, 1)
            b_eq = 1
        }
        else if (rts == "drs") {
            A_iq = A_iq \ J(1, n, 1)
            b_iq = b_iq \ 1
        }
        else if (rts == "irs") {
            A_iq = A_iq \ J(1, n, -1)
            b_iq = b_iq \ (-1)
        }

        res = _rsm_solve_lp(A_iq, b_iq, A_eq, b_eq, c, 1e-8)

        if (res.status == 0) {
            max_profit = -res.objval
            if (max_profit > 0) profit_eff = obs_profit / max_profit
        }

        st_store(obs_k, st_varindex(eff_var), profit_eff)
    }
}


// ============================================================================
// DIRECTIONAL DISTANCE FUNCTION (DDF)
// ============================================================================
/**
 * _dea_simplex_ddf()
 *
 * max β  s.t.  X'λ + β·gx ≤ x₀,  -Y'λ + β·gy ≤ -y₀,
 *              [convexity],  λ,β ≥ 0
 *
 * Interface matches _dea_ipm_ddf() exactly.
 */
void _dea_simplex_ddf(string scalar inputs,
                       string scalar outputs,
                       string scalar eff_var,
                       string scalar gx_var,
                       string scalar gy_var,
                       string scalar touse,
                       string scalar rts)
{
    real matrix X, Y, GX, GY
    real colvector obs_indices, x0, y0, gx0, gy0
    real rowvector c, sol
    real matrix A_iq, A_eq
    real colvector b_iq, b_eq
    real scalar n, m, s, k, obs_k, beta

    struct RSM_Result scalar res

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))

    n = rows(X); m = cols(X); s = cols(Y)

    if (gx_var != "") GX = st_data(., tokens(gx_var), touse)
    else GX = X
    if (gy_var != "") GY = st_data(., tokens(gy_var), touse)
    else GY = Y

    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, .]'; y0 = Y[k, .]'
        gx0 = GX[k, .]'; gy0 = GY[k, .]'
        beta = .

        // Variables: [λ(n), β]
        // max β → min -β
        c = J(1, n, 0), -1

        // X'λ + β·gx ≤ x₀;   -Y'λ - β·gy ≤ -y₀
        A_iq = (X', gx0) \ (-Y', gy0)
        b_iq = x0 \ (-y0)

        A_eq = J(0, n+1, 0)
        b_eq = J(0, 1, 0)

        if (rts == "vrs") {
            A_eq = (J(1, n, 1), 0)
            b_eq = 1
        }
        else if (rts == "drs") {
            A_iq = A_iq \ (J(1, n, 1), 0)
            b_iq = b_iq \ 1
        }
        else if (rts == "irs") {
            A_iq = A_iq \ (J(1, n, -1), 0)
            b_iq = b_iq \ (-1)
        }

        res = _rsm_solve_lp(A_iq, b_iq, A_eq, b_eq, c, 1e-8)

        if (res.status == 0) {
            sol = res.x'
            beta = sol[n+1]
        }

        st_store(obs_k, st_varindex(eff_var), beta)
    }
}


// ============================================================================
// DDF WITH UNDESIRABLE OUTPUTS
// ============================================================================
/**
 * _dea_simplex_ddf_undesirable()
 *
 * max β  s.t.  X'λ + β·gx ≤ x₀,   (inputs — no contraction if gx=0)
 *              Y'λ ≥ y₀ + β·gy,    (good outputs expand)
 *              B'λ + β·gb = b₀,    (bad outputs: equality, weak disposability)
 *              [convexity],  λ,β ≥ 0
 *
 * direction: "yb"  → gx=0,  gy=y₀, gb=b₀  (default)
 *            "xyb" → gx=x₀, gy=y₀, gb=b₀
 *            "x"   → gx=x₀, gy=0,  gb=0
 */
void _dea_simplex_ddf_undesirable(string scalar inputs,
                                   string scalar outputs,
                                   string scalar bad_outputs,
                                   string scalar eff_var,
                                   string scalar direction,
                                   string scalar touse,
                                   string scalar rts)
{
    real matrix X, Y, B
    real colvector obs_indices, x0, y0, b0, gx0, gy0, gb0
    real rowvector c, lb, ub, sol
    real matrix A_iq, A_eq
    real colvector b_iq, b_eq
    real scalar n, m, s, q, k, obs_k, beta, rc
    class LinearProgram scalar S
    string scalar dir

    X = st_data(., tokens(inputs),      touse)
    Y = st_data(., tokens(outputs),     touse)
    B = st_data(., tokens(bad_outputs), touse)
    obs_indices = selectindex(st_data(., touse))

    n = rows(X)
    m = cols(X)
    s = cols(Y)
    q = cols(B)

    dir = (direction == "" ? "yb" : direction)

    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, .]'
        y0 = Y[k, .]'
        b0 = B[k, .]'
        beta = .

        if (dir == "xyb") {
            gx0 = x0
            gy0 = y0
            gb0 = b0
        }
        else if (dir == "x") {
            gx0 = x0
            gy0 = J(s, 1, 0)
            gb0 = J(q, 1, 0)
        }
        else {
            gx0 = J(m, 1, 0)
            gy0 = y0
            gb0 = b0
        }

        c  = J(1, n, 0), -1
        lb = J(1, n + 1, 0)
        ub = J(1, n + 1, .)

        A_iq = (X', gx0) \ (-Y', -gy0)
        b_iq = x0 \ (-y0)

        A_eq = (B', gb0)
        b_eq = b0

        if (rts == "vrs") {
            A_eq = A_eq \ (J(1, n, 1), 0)
            b_eq = b_eq \ 1
        }
        else if (rts == "drs") {
            A_iq = A_iq \ (J(1, n, 1), 0)
            b_iq = b_iq \ 1
        }
        else if (rts == "irs") {
            A_iq = A_iq \ (J(1, n, -1), 0)
            b_iq = b_iq \ (-1)
        }

        S = LinearProgram()
        S.setMaxOrMin("min")
        S.setCoefficients(c)
        S.setInequality(A_iq, b_iq)
        S.setEquality(A_eq, b_eq)
        S.setBounds(lb, ub)
        S.setTrace("off")

        rc = S.optimize()
        if (S.converged() && S.errorcode() == 0) {
            sol = S.parameters()
            beta = sol[n+1]
        }

        st_store(obs_k, st_varindex(eff_var), beta)
    }
}


// ============================================================================
// ADDITIVE MODEL
// ============================================================================
/**
 * _dea_simplex_additive()
 *
 * max Σs⁻ + Σs⁺  s.t.  X'λ + s⁻ = x₀,  Y'λ - s⁺ = y₀,
 *                        [convexity],  λ,s⁻,s⁺ ≥ 0
 *
 * Score = total slack (0 = efficient).
 * Interface matches _dea_ipm_additive() exactly.
 */
void _dea_simplex_additive(string scalar inputs,
                            string scalar outputs,
                            string scalar eff_var,
                            string scalar slack_in_var,
                            string scalar slack_out_var,
                            string scalar touse,
                            string scalar rts)
{
    real matrix X, Y
    real colvector obs_indices, x0, y0
    real rowvector c, sol, slack_in, slack_out
    real matrix A_iq, A_eq
    real colvector b_iq, b_eq
    real scalar n, m, s, k, obs_k, add_score, j

    struct RSM_Result scalar res

    string rowvector si_varlist, so_varlist
    string scalar _sv_j
    real scalar has_slacks, _sval_j

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))

    n = rows(X); m = cols(X); s = cols(Y)

    has_slacks = (slack_in_var != "" & slack_out_var != "")
    if (has_slacks) {
        si_varlist = tokens(slack_in_var)
        so_varlist = tokens(slack_out_var)
    }

    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, .]'; y0 = Y[k, .]'
        add_score = .
        slack_in = J(1, m, 0); slack_out = J(1, s, 0)

        // Variables: [λ(n), s⁻(m), s⁺(s)]
        // max Σs⁻+Σs⁺ → min -Σs⁻-Σs⁺
        c = J(1, n, 0), J(1, m, -1), J(1, s, -1)

        // Equalities: X'λ+s⁻=x₀,  Y'λ-s⁺=y₀
        A_eq = (X', I(m), J(m, s, 0)) \ (Y', J(s, m, 0), -I(s))
        b_eq = x0 \ y0

        A_iq = J(0, n+m+s, 0)
        b_iq = J(0, 1, 0)

        if (rts == "vrs") {
            A_eq = A_eq \ (J(1, n, 1), J(1, m+s, 0))
            b_eq = b_eq \ 1
        }
        else if (rts == "drs") {
            A_iq = (J(1, n, 1), J(1, m+s, 0))
            b_iq = 1
        }
        else if (rts == "irs") {
            A_iq = (J(1, n, -1), J(1, m+s, 0))
            b_iq = -1
        }

        res = _rsm_solve_lp(A_iq, b_iq, A_eq, b_eq, c, 1e-8)

        if (res.status == 0) {
            sol = res.x'
            slack_in = sol[n+1..n+m]
            slack_out = sol[n+m+1..n+m+s]
            add_score = sum(slack_in) + sum(slack_out)
        }

        st_store(obs_k, st_varindex(eff_var), add_score)

        if (has_slacks) {
            for (j = 1; j <= m; j++) {
                _sv_j = si_varlist[j]
                _sval_j = (slack_in[j])
                st_store(obs_k, st_varindex(_sv_j), _sval_j)
            }
            for (j = 1; j <= s; j++) {
                _sv_j = so_varlist[j]
                _sval_j = (slack_out[j])
                st_store(obs_k, st_varindex(_sv_j), _sval_j)
            }
        }
    }
}


// ============================================================================
// SBM (SLACKS-BASED MEASURE) MODEL
// ============================================================================
/**
 * _dea_simplex_sbm()
 *
 * Tone (1997, Saitama Univ. Research Report; published EJOR 2001) SBM.
 * Non-oriented SBM: Charnes-Cooper transformation.
 * Oriented SBM (input/output): direct LP formulation following the
 *   orientation framework of Ali, Lerme & Seiford (1995, EJOR 80, 462-473).
 *
 * Transformed LP (all equalities):
 *   min t - (1/m)·Σ(S⁻ᵢ/x₀ᵢ)
 *   s.t. t + (1/s)·Σ(S⁺ᵣ/y₀ᵣ) = 1       (normalizing)
 *        X'Λ + S⁻ - t·x₀ = 0              (input)
 *        Y'Λ - S⁺ - t·y₀ = 0              (output)
 *        [convexity: ΣΛ = t  or  ΣΛ ≤/≥ t]
 *        Λ, S⁻, S⁺, t ≥ 0
 *
 * Variables: [Λ(n), S⁻(m), S⁺(s), t]
 * Interface matches _dea_ipm_sbm() exactly.
 */
void _dea_simplex_sbm(string scalar inputs,
                       string scalar outputs,
                       string scalar eff_var,
                       string scalar slack_in_var,
                       string scalar slack_out_var,
                       string scalar touse,
                       string scalar rts,
                       string scalar orientation)
{
    real matrix X, Y
    real colvector obs_indices, x0, y0
    real rowvector c, sol, coef_sin, coef_sout, slack_in, slack_out
    real matrix A_iq, A_eq
    real colvector b_iq, b_eq
    real scalar n, m, s, k, obs_k, sbm_eff, t_val, z_val, j, nv

    struct RSM_Result scalar res

    string rowvector si_varlist, so_varlist
    string scalar _sv_j
    real scalar has_slacks, _sval_j
    string scalar ort

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))

    n = rows(X); m = cols(X); s = cols(Y)

    ort = strlower(orientation)
    if (ort == "" | ort == "inout") ort = "inout"
    else if (ort == "in" | ort == "input") ort = "in"
    else ort = "out"

    has_slacks = (slack_in_var != "" & slack_out_var != "")
    if (has_slacks) {
        si_varlist = tokens(slack_in_var)
        so_varlist = tokens(slack_out_var)
    }

    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, .]'; y0 = Y[k, .]'
        sbm_eff = .
        slack_in = J(1, m, 0); slack_out = J(1, s, 0)

        // Coefficients
        coef_sin = J(1, m, 0)
        for (j = 1; j <= m; j++) {
            if (x0[j] > 0) coef_sin[j] = -1 / (m * x0[j])
        }
        coef_sout = J(1, s, 0)
        for (j = 1; j <= s; j++) {
            if (y0[j] > 0) coef_sout[j] = 1 / (s * y0[j])
        }

        if (ort == "inout") {
            // ── Non-oriented: Charnes-Cooper transformation ─────────────────
            nv = n + m + s + 1

            // Objective: min t + Σ(coef_sin · S⁻)
            c = J(1, n, 0), coef_sin, J(1, s, 0), 1

            // (1) Normalizing: t + Σ(coef_sout · S⁺) = 1
            A_eq = (J(1, n, 0), J(1, m, 0), coef_sout, 1)
            b_eq = 1

            // (2) Input: X'Λ + S⁻ - t·x₀ = 0
            A_eq = A_eq \ (X', I(m), J(m, s, 0), -x0)
            b_eq = b_eq \ J(m, 1, 0)

            // (3) Output: Y'Λ - S⁺ - t·y₀ = 0
            A_eq = A_eq \ (Y', J(s, m, 0), -I(s), -y0)
            b_eq = b_eq \ J(s, 1, 0)

            A_iq = J(0, nv, 0)
            b_iq = J(0, 1, 0)

            if (rts == "vrs") {
                A_eq = A_eq \ (J(1, n, 1), J(1, m+s, 0), -1)
                b_eq = b_eq \ 0
            }
            else if (rts == "drs") {
                A_iq = (J(1, n, 1), J(1, m+s, 0), -1)
                b_iq = 0
            }
            else if (rts == "irs") {
                A_iq = (J(1, n, -1), J(1, m+s, 0), 1)
                b_iq = 0
            }

            res = _rsm_solve_lp(A_iq, b_iq, A_eq, b_eq, c, 1e-8)

            if (res.status == 0) {
                sol = res.x'
                t_val = sol[nv]
                if (t_val > 0) {
                    slack_in  = sol[n+1..n+m]   / t_val
                    slack_out = sol[n+m+1..n+m+s] / t_val
                    sbm_eff   = res.objval
                }
            }
        }
        else {
            // ── Oriented: direct LP, variables [λ(n), s⁻(m), s⁺(s)] ────────
            nv = n + m + s

            if (ort == "in") {
                // min -(1/m)·Σ(s⁻ᵢ/x₀ᵢ);  ρ* = 1 + objval
                c = J(1, n, 0), coef_sin, J(1, s, 0)
            }
            else {
                // out: min -(1/s)·Σ(s⁺ᵣ/y₀ᵣ);  ρ* = 1/(1-objval)
                c = J(1, n, 0), J(1, m, 0), -coef_sout
            }

            // Input balance: X'λ + s⁻ = x₀
            A_eq = (X', I(m), J(m, s, 0))
            b_eq = x0

            // Output balance: Y'λ - s⁺ = y₀
            A_eq = A_eq \ (Y', J(s, m, 0), -I(s))
            b_eq = b_eq \ y0

            A_iq = J(0, nv, 0)
            b_iq = J(0, 1, 0)

            if (rts == "vrs") {
                A_eq = A_eq \ (J(1, n, 1), J(1, m+s, 0))
                b_eq = b_eq \ 1
            }
            else if (rts == "drs") {
                A_iq = (J(1, n, 1), J(1, m+s, 0))
                b_iq = 1
            }
            else if (rts == "irs") {
                A_iq = (J(1, n, -1), J(1, m+s, 0))
                b_iq = -1
            }

            res = _rsm_solve_lp(A_iq, b_iq, A_eq, b_eq, c, 1e-8)

            if (res.status == 0) {
                sol = res.x'
                slack_in  = sol[n+1..n+m]
                slack_out = sol[n+m+1..n+m+s]
                if (ort == "in") {
                    sbm_eff = 1 + res.objval
                }
                else {
                    z_val = -res.objval
                    if (1 + z_val > 0) sbm_eff = 1 / (1 + z_val)
                }
            }
        }

        st_store(obs_k, st_varindex(eff_var), sbm_eff)

        if (has_slacks) {
            for (j = 1; j <= m; j++) {
                _sv_j = si_varlist[j]
                _sval_j = (slack_in[j])
                st_store(obs_k, st_varindex(_sv_j), _sval_j)
            }
            for (j = 1; j <= s; j++) {
                _sv_j = so_varlist[j]
                _sval_j = (slack_out[j])
                st_store(obs_k, st_varindex(_sv_j), _sval_j)
            }
        }
    }
}


// ============================================================================
// SUPER-EFFICIENCY MODEL
// ============================================================================
/**
 * _dea_simplex_supereff()
 *
 * Same LP as TE but DMU k is excluded from reference set.
 * Input-oriented:  theta > 1 → super-efficient (stored as theta)
 * Output-oriented: phi < 1 → super-efficient (stored as 1/phi > 1)
 *
 * Arguments:
 *   lambda_var   - string: lambda variable names (n, space-separated), "" if none
 *                  LP solves over n-1 lambdas (excluding DMU k);
 *                  results are mapped back to n positions with lambda_k = 0.
 *
 * Interface matches _dea_ipm_supereff() exactly.
 */
void _dea_simplex_supereff(string scalar inputs,
                            string scalar outputs,
                            string scalar eff_var,
                            string scalar lambda_var,
                            string scalar touse,
                            string scalar orientation,
                            string scalar rts)
{
    real matrix X, Y, X_ref, Y_ref
    real colvector obs_indices, x0, y0
    real rowvector c, lam_ref, lam_k
    real matrix A_iq, A_eq
    real colvector b_iq, b_eq
    real scalar n, m, s, k, obs_k, n_ref, j, eff_score, phi_raw, has_lambda
    string scalar _lv_j
    real scalar _lval_j

    struct RSM_Result scalar res

    // Lambda variable handling
    string rowvector lam_varlist
    has_lambda = (lambda_var != "")
    if (has_lambda) lam_varlist = tokens(lambda_var)

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))

    n = rows(X); m = cols(X); s = cols(Y)

    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, .]'; y0 = Y[k, .]'
        eff_score = .
        lam_k = J(1, n, 0)  // Initialize: all zeros (k-th is always 0)
        n_ref = n - 1

        // Build reference set excluding DMU k
        if (k == 1) {
            X_ref = X[2..n, ]; Y_ref = Y[2..n, ]
        }
        else if (k == n) {
            X_ref = X[1..n-1, ]; Y_ref = Y[1..n-1, ]
        }
        else {
            X_ref = X[1..k-1, ] \ X[k+1..n, ]
            Y_ref = Y[1..k-1, ] \ Y[k+1..n, ]
        }

        if (orientation == "input") {
            // min θ  s.t. X_ref'λ - θ·x₀ ≤ 0,  -Y_ref'λ ≤ -y₀
            // Variables: [λ(1..n_ref), θ]
            c = J(1, n_ref, 0), 1

            A_iq = (X_ref', -x0) \ (-Y_ref', J(s, 1, 0))
            b_iq = J(m, 1, 0) \ (-y0)

            A_eq = J(0, n_ref+1, 0)
            b_eq = J(0, 1, 0)

            if (rts == "vrs") {
                A_eq = (J(1, n_ref, 1), 0)
                b_eq = 1
            }
            else if (rts == "drs") {
                A_iq = A_iq \ (J(1, n_ref, 1), 0)
                b_iq = b_iq \ 1
            }
            else if (rts == "irs") {
                A_iq = A_iq \ (J(1, n_ref, -1), 0)
                b_iq = b_iq \ (-1)
            }

            res = _rsm_solve_lp(A_iq, b_iq, A_eq, b_eq, c, 1e-8)

            if (res.status == 0) {
                eff_score = res.objval
                lam_ref = res.x[1..n_ref]'
                
                // Map n-1 lambdas back to n positions (k-th = 0)
                if (k == 1) {
                    lam_k[2..n] = lam_ref
                }
                else if (k == n) {
                    lam_k[1..n-1] = lam_ref
                }
                else {
                    lam_k[1..k-1] = lam_ref[1..k-1]
                    lam_k[k+1..n] = lam_ref[k..n_ref]
                }
            }
        }
        else {
            // max φ → min -φ
            // X_ref'λ ≤ x₀,  -Y_ref'λ + φ·y₀ ≤ 0
            // Variables: [λ(1..n_ref), φ]
            c = J(1, n_ref, 0), -1

            A_iq = (X_ref', J(m, 1, 0)) \ (-Y_ref', y0)
            b_iq = x0 \ J(s, 1, 0)

            A_eq = J(0, n_ref+1, 0)
            b_eq = J(0, 1, 0)

            if (rts == "vrs") {
                A_eq = (J(1, n_ref, 1), 0)
                b_eq = 1
            }
            else if (rts == "drs") {
                A_iq = A_iq \ (J(1, n_ref, 1), 0)
                b_iq = b_iq \ 1
            }
            else if (rts == "irs") {
                A_iq = A_iq \ (J(1, n_ref, -1), 0)
                b_iq = b_iq \ (-1)
            }

            res = _rsm_solve_lp(A_iq, b_iq, A_eq, b_eq, c, 1e-8)

            if (res.status == 0) {
                phi_raw = -res.objval
                if (phi_raw > 0) eff_score = 1 / phi_raw
                lam_ref = res.x[1..n_ref]'
                
                // Map n-1 lambdas back to n positions (k-th = 0)
                if (k == 1) {
                    lam_k[2..n] = lam_ref
                }
                else if (k == n) {
                    lam_k[1..n-1] = lam_ref
                }
                else {
                    lam_k[1..k-1] = lam_ref[1..k-1]
                    lam_k[k+1..n] = lam_ref[k..n_ref]
                }
            }
        }

        // Store efficiency score
        st_store(obs_k, st_varindex(eff_var), eff_score)
        
        // Store lambda weights (mapped to original n DMU positions)
        if (has_lambda) {
            for (j = 1; j <= n; j++) {
                _lv_j = lam_varlist[j]
                _lval_j = (lam_k[j])
                st_store(obs_k, st_varindex(_lv_j), _lval_j)
            }
        }
    }
}


// ============================================================================
// SUPER-EFFICIENCY SBM MODEL
// ============================================================================
/**
 * _dea_simplex_supersbm()
 *
 * Tone (2002) Super-SBM: same as SBM but DMU k excluded from reference set.
 * Uses Charnes-Cooper transformation.
 *
 * Interface matches _dea_ipm_supersbm() exactly.
 */
void _dea_simplex_supersbm(string scalar inputs,
                            string scalar outputs,
                            string scalar eff_var,
                            string scalar touse,
                            string scalar rts)
{
    real matrix X, Y, X_ref, Y_ref
    real colvector obs_indices, x0, y0
    real rowvector c, coef_xbar, coef_ybar
    real matrix A_iq, A_eq
    real colvector b_iq, b_eq
    real scalar n, m, s, k, obs_k, n_ref, j
    real scalar sbm_eff, sbm_score, tol

    struct RSM_Result scalar res

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))

    n = rows(X); m = cols(X); s = cols(Y)
    tol = 1e-6

    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, .]'; y0 = Y[k, .]'
        sbm_eff = .

        // ==============================================================
        // Step 1: Standard SBM to identify efficient DMUs
        //   Charnes-Cooper linearization of Tone (2001):
        //   Variables: [lambda_hat(n), s_minus(m), s_plus(s), t]
        //   min  t - (1/m) * sum(s_minus_i / x_i0)
        //   s.t. t + (1/s) * sum(s_plus_r / y_r0) = 1
        //        X' * lambda_hat + s_minus - t * x_0 = 0
        //        Y' * lambda_hat - s_plus  - t * y_0 = 0
        //        VRS: sum(lambda_hat) = t
        //        all >= 0
        // ==============================================================
        {
            real rowvector c_sbm, coef_sm, coef_sp
            real matrix A_eq_sbm, A_iq_sbm
            real colvector b_eq_sbm, b_iq_sbm
            struct RSM_Result scalar res_sbm

            // Coefficients: -(1/m)/x_{i0} for s_minus; (1/s)/y_{r0} for s_plus
            coef_sm = J(1, m, 0)
            for (j = 1; j <= m; j++) {
                if (x0[j] > 0) coef_sm[j] = -1 / (m * x0[j])
            }
            coef_sp = J(1, s, 0)
            for (j = 1; j <= s; j++) {
                if (y0[j] > 0) coef_sp[j] = 1 / (s * y0[j])
            }

            // Objective: min t - (1/m) * sum(s_minus_i / x_i0)
            c_sbm = J(1, n, 0), coef_sm, J(1, s, 0), 1

            // Eq 1: t + (1/s) * sum(s_plus_r / y_r0) = 1
            A_eq_sbm = (J(1, n, 0), J(1, m, 0), coef_sp, 1)
            b_eq_sbm = 1

            // Eq 2: X'*lambda + s_minus - t*x_0 = 0
            A_eq_sbm = A_eq_sbm \ (X', I(m), J(m, s, 0), -x0)
            b_eq_sbm = b_eq_sbm \ J(m, 1, 0)

            // Eq 3: Y'*lambda - s_plus - t*y_0 = 0
            A_eq_sbm = A_eq_sbm \ (Y', J(s, m, 0), -I(s), -y0)
            b_eq_sbm = b_eq_sbm \ J(s, 1, 0)

            // RTS
            if (rts == "vrs") {
                A_eq_sbm = A_eq_sbm \ (J(1, n, 1), J(1, m+s, 0), -1)
                b_eq_sbm = b_eq_sbm \ 0
            }

            A_iq_sbm = J(0, n+m+s+1, 0)
            b_iq_sbm = J(0, 1, 0)
            if (rts == "drs") {
                A_iq_sbm = (J(1, n, 1), J(1, m+s, 0), -1)
                b_iq_sbm = 0
            }
            else if (rts == "irs") {
                A_iq_sbm = (J(1, n, -1), J(1, m+s, 0), 1)
                b_iq_sbm = 0
            }

            res_sbm = _rsm_solve_lp(A_iq_sbm, b_iq_sbm, A_eq_sbm, b_eq_sbm,
                                     c_sbm, 1e-8)

            sbm_score = (res_sbm.status == 0 ? res_sbm.objval : .)
        }

        // If SBM-inefficient, return SBM score directly
        if (sbm_score == . || sbm_score < 1 - tol) {
            st_store(obs_k, st_varindex(eff_var), sbm_score)
            continue
        }

        // ==============================================================
        // Step 2: Super-SBM for efficient DMUs — Tone (2002) eq. 14
        //   Charnes-Cooper linearization:
        //   Variables: [Lambda(n_ref), x_tilde(m), y_tilde(s), t]
        //
        //   min  (1/m) * sum(x_tilde_i / x_i0)
        //   s.t. (1/s) * sum(y_tilde_r / y_r0) = 1
        //        x_tilde >= X_ref' * Lambda
        //        y_tilde <= Y_ref' * Lambda
        //        x_tilde >= t * x_0              <-- was missing
        //        y_tilde <= t * y_0
        //        VRS: sum(Lambda) = t
        //        Lambda >= 0, y_tilde >= 0, t > 0
        // ==============================================================

        n_ref = n - 1
        if (n_ref < 1) {
            st_store(obs_k, st_varindex(eff_var), .)
            continue
        }

        if (k == 1) {
            X_ref = X[2..n, ]; Y_ref = Y[2..n, ]
        }
        else if (k == n) {
            X_ref = X[1..n-1, ]; Y_ref = Y[1..n-1, ]
        }
        else {
            X_ref = X[1..k-1, ] \ X[k+1..n, ]
            Y_ref = Y[1..k-1, ] \ Y[k+1..n, ]
        }

        coef_xbar = J(1, m, 0)
        for (j = 1; j <= m; j++) {
            if (x0[j] > 0) coef_xbar[j] = 1 / (m * x0[j])
        }
        coef_ybar = J(1, s, 0)
        for (j = 1; j <= s; j++) {
            if (y0[j] > 0) coef_ybar[j] = 1 / (s * y0[j])
        }

        // Objective: min (1/m) * sum(x_tilde_i / x_i0)   [FIX: no +t]
        c = J(1, n_ref, 0), coef_xbar, J(1, s, 0), 0

        // Normalizing: (1/s) * sum(y_tilde_r / y_r0) = 1 [FIX: no t]
        A_eq = (J(1, n_ref, 0), J(1, m, 0), coef_ybar, 0)
        b_eq = 1

        // Inequality 1: x_tilde >= X_ref' * Lambda
        A_iq = (X_ref', -I(m), J(m, s, 0), J(m, 1, 0))
        b_iq = J(m, 1, 0)

        // Inequality 2: y_tilde <= Y_ref' * Lambda
        A_iq = A_iq \ (-Y_ref', J(s, m, 0), I(s), J(s, 1, 0))
        b_iq = b_iq \ J(s, 1, 0)

        // Inequality 3: x_tilde >= t * x_0  [FIX: was missing]
        A_iq = A_iq \ (J(m, n_ref, 0), -I(m), J(m, s, 0), x0)
        b_iq = b_iq \ J(m, 1, 0)

        // Inequality 4: y_tilde <= t * y_0
        A_iq = A_iq \ (J(s, n_ref, 0), J(s, m, 0), I(s), -y0)
        b_iq = b_iq \ J(s, 1, 0)

        // RTS
        if (rts == "vrs") {
            A_eq = A_eq \ (J(1, n_ref, 1), J(1, m+s, 0), -1)
            b_eq = b_eq \ 0
        }
        else if (rts == "drs") {
            A_iq = A_iq \ (J(1, n_ref, 1), J(1, m+s, 0), -1)
            b_iq = b_iq \ 0
        }
        else if (rts == "irs") {
            A_iq = A_iq \ (J(1, n_ref, -1), J(1, m+s, 0), 1)
            b_iq = b_iq \ 0
        }

        res = _rsm_solve_lp(A_iq, b_iq, A_eq, b_eq, c, 1e-8)

        if (res.status == 0) {
            sbm_eff = res.objval
        }

        st_store(obs_k, st_varindex(eff_var), sbm_eff)
    }
}


// ============================================================================
// ALLOCATIVE EFFICIENCY MODEL
// ============================================================================
/**
 * _dea_simplex_allocative()
 *
 * Allocative Efficiency = Economic Efficiency / Technical Efficiency
 *   Input-oriented:  AE = Cost_Eff / TE
 *   Output-oriented: AE = Revenue_Eff / TE
 *
 * This is a composite model that calls TE + Cost/Revenue internally.
 * Interface matches _dea_ipm_allocative() (via _dea_ipm_allocative program).
 */
void _dea_simplex_allocative(string scalar inputs,
                              string scalar outputs,
                              string scalar eff_var,
                              string scalar touse,
                              string scalar orientation,
                              string scalar rts,
                              string scalar w_var,
                              string scalar p_var)
{
    real matrix X, Y
    real colvector obs_indices, x0, y0
    real scalar n, m, s, k, obs_k
    real scalar te_k, econ_k, alloc_k
    real colvector te_result

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))

    n = rows(X); m = cols(X); s = cols(Y)

    // Prices
    real matrix W, P
    real colvector w0, p0
    real rowvector c_econ
    real matrix A_iq, A_eq
    real colvector b_iq, b_eq
    real scalar obs_val, opt_val

    struct RSM_Result scalar res

    if (orientation == "input" & w_var != "") {
        W = st_data(., tokens(w_var), touse)
    }
    if (orientation == "output" & p_var != "") {
        P = st_data(., tokens(p_var), touse)
    }

    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, .]'; y0 = Y[k, .]'
        alloc_k = .

        // --- Step 1: Technical Efficiency ---
        te_result = _dea_solve_dmu(X, Y, k, orientation, rts, 1e-8, "")
        te_k = te_result[1]

        if (te_k == . || te_k == 0) {
            st_store(obs_k, st_varindex(eff_var), .)
            continue
        }

        // --- Step 2: Economic Efficiency ---
        if (orientation == "input") {
            w0 = W[k, .]'
            obs_val = sum(w0 :* x0)

            // Cost minimization
            c_econ = J(1, n, 0), w0'
            A_iq = (X', -I(m)) \ (-Y', J(s, m, 0))
            b_iq = J(m, 1, 0) \ (-y0)
            A_eq = J(0, n+m, 0); b_eq = J(0, 1, 0)

            if (rts == "vrs") {
                A_eq = (J(1, n, 1), J(1, m, 0)); b_eq = 1
            }
            else if (rts == "drs") {
                A_iq = A_iq \ (J(1, n, 1), J(1, m, 0)); b_iq = b_iq \ 1
            }
            else if (rts == "irs") {
                A_iq = A_iq \ (J(1, n, -1), J(1, m, 0)); b_iq = b_iq \ (-1)
            }

            res = _rsm_solve_lp(A_iq, b_iq, A_eq, b_eq, c_econ, 1e-8)

            if (res.status == 0 && obs_val > 0) {
                econ_k = res.objval / obs_val
                alloc_k = econ_k / te_k
            }
        }
        else {
            p0 = P[k, .]'
            obs_val = sum(p0 :* y0)

            // Revenue maximization (min -p'y)
            c_econ = J(1, n, 0), -p0'
            A_iq = (X', J(m, s, 0)) \ (-Y', I(s))
            b_iq = x0 \ J(s, 1, 0)
            A_eq = J(0, n+s, 0); b_eq = J(0, 1, 0)

            if (rts == "vrs") {
                A_eq = (J(1, n, 1), J(1, s, 0)); b_eq = 1
            }
            else if (rts == "drs") {
                A_iq = A_iq \ (J(1, n, 1), J(1, s, 0)); b_iq = b_iq \ 1
            }
            else if (rts == "irs") {
                A_iq = A_iq \ (J(1, n, -1), J(1, s, 0)); b_iq = b_iq \ (-1)
            }

            res = _rsm_solve_lp(A_iq, b_iq, A_eq, b_eq, c_econ, 1e-8)

            if (res.status == 0 && obs_val > 0) {
                opt_val = -res.objval
                econ_k = obs_val / opt_val
                alloc_k = econ_k / te_k
            }
        }

        st_store(obs_k, st_varindex(eff_var), alloc_k)
    }
}


end
