*! _dea_simplex.ado
*! SBM: zero-term means, unit rescaling, and oriented super-SBM.
*! Forward non-TE tolerance and validate original LP constraints before returning results.
*! Reject Phase-I and final real-variable solutions outside the shared residual tolerance.
*! Two-Phase Revised Simplex LP Solver for DEA
*! Version 2.0.25 - 2026-11-05
*!
*! Changelog:
*!   2.0.25 05nov2026  Bounded certified stage-link retry from original rows;
*!                      per-DMU old/new radial diagnostics and matched loaders.
*!   2.0.24 05nov2026  Bounded exact original-row selection and empty-basis proof;
*!                      certified row recovery after all earlier retries fail.
*!   2.0.23 04nov2026  Check successful bases in original coordinates; refine only
*!                      verified defects, with exact-zero rows and fast bounds.
*!   2.0.22 04nov2026  Recover exact binary-power duplicate rows/columns only
*!                      after earlier failures; restore and check original units.
*!   2.0.21 04nov2026  Retry blocked Phase-1 pivots with an original-column basis
*!                      only after earlier paths fail; retain original-unit checks.
*!   2.0.20 04nov2026  Bound basis-solve errors in certified prices; stop repeated
*!                      Bland states and retain existing restart priority.
*!   2.0.19 03nov2026  Guard failure-only basis exchanges and reconstruct real-column
*!                      bases without dropping original rows; extend wide selection.
*!   2.0.18 03nov2026  Preserve accepted results; refine failure-only artificial-free
*!                      retries and bound their primal-dual quality checks.
*!   2.0.17 03nov2026  Check both inverse identities with compensated residuals;
*!                      conservative fast acceptance and exact binary scaling.
*!   2.0.12 02nov2026  Return projection residual diagnostics; preserve scores and checked references.
*!   2.0.11 02nov2026  Report stage status/provenance; retain only the radial
*!                      score when requested stage-2 projections fail.
*!   2.0.10 02nov2026  Validate basis inverses; propagate numerical failures;
*!                      strict Bland ties, feasibility-bounded Harris test,
*!                      redundant-row handling and one equilibrated restart.
*!   2.0.9  01nov2026  Check original LP residuals after both phases.
*!   2.0.8  01nov2026  Honor stage(1); preserve missing reference/slack results
*!                      when only a dual score is available.
*!   2.0.7  01nov2026  Correct good-output beta sign in undesirable-output DDF.
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
*!     bloom.rampike@gmail.com
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
// Load the shared acceptance criterion before compiling solver functions.
capture mata: assert(_dea_lpcheck_ok(J(1,6,0),1e-8))
if _rc {
    quietly findfile _dea_lpcheck.ado
    quietly run "`r(fn)'"
}

capture mata: assert(_dea_nonte_lp_ready()==1)
if _rc {
    quietly findfile _dea_nonte_lp.ado
    quietly run "`r(fn)'"
}

capture mata: assert(_dea_teaudit_ready()==1)
if _rc {
    quietly findfile _dea_teaudit.ado
    quietly run "`r(fn)'"
}

capture mata: assert(_dea_recovery_ready()==1)
if _rc {
    quietly findfile _dea_recovery.ado
    quietly run "`r(fn)'"
}
capture mata: mata drop _dea_simplex_audit_ready()
capture mata: mata drop RSM_Prob()
capture mata: mata drop RSM_Basis()
capture mata: mata drop RSM_Result()
capture mata: mata drop DEA_Info()
capture mata: mata drop _dea_simplex_envelopment()
capture mata: mata drop _dea_simplex_full()
capture mata: mata drop _dea_link_reset()
capture mata: mata drop _dea_link_read()
capture mata: mata drop _dea_link_store()
capture mata: mata drop _dea_link_small()
capture mata: mata drop _dea_link_residual()
capture mata: mata drop _dea_link_basis()
capture mata: mata drop _dea_link_apply()
capture mata: mata drop _dea_link_stage2()
capture mata: mata drop _dea_solve_dmu()
capture mata: mata drop _dea_lp_stage1()
capture mata: mata drop _dea_lp_stage2()
capture mata: mata drop _rsm_two_phase()
capture mata: mata drop _rsm_duplicate_ratio()
capture mata: mata drop _rsm_duplicate_map()
capture mata: mata drop _rsm_duplicate_problem()
capture mata: mata drop _rsm_duplicate_accept()
capture mata: mata drop _rsm_duplicate_retry()
capture mata: mata drop _rsm_quality_fast()
capture mata: mata drop _rsm_affine_quality()
capture mata: mata drop _rsm_affine_refine()
capture mata: mata drop _rsm_success_fastread()
capture mata: mata drop _rsm_xgrow()
capture mata: mata drop _rsm_xadd()
capture mata: mata drop _rsm_xmul()
capture mata: mata drop _rsm_xnormalize()
capture mata: mata drop _rsm_rowmap()
capture mata: mata drop _rsm_row_empty()
capture mata: mata drop _rsm_row_accept()
capture mata: mata drop _rsm_row_retry()
capture mata: mata drop _rsm_success_reset()
capture mata: mata drop _rsm_success_read()
capture mata: mata drop _rsm_success_counts_read()
capture mata: mata drop _rsm_success_check()
capture mata: mata drop _rsm_success_apply()
capture mata: mata drop _rsm_two_phase_core()
capture mata: mata drop _rsm_two_phase_once()
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
capture mata: mata drop _dea_nonte_simplex()
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

capture mata: mata drop _rsm_eta_two_phase()
capture mata: mata drop _rsm_factor_residual()
capture mata: mata drop _rsm_inverse_residual()
capture mata: mata drop _rsm_inverse_fast()
capture mata: mata drop _rsm_factor_solve()
capture mata: mata drop _rsm_retry_refine()
capture mata: mata drop _rsm_q_sum()
capture mata: mata drop _rsm_q_add()
capture mata: mata drop _rsm_q_product()
capture mata: mata drop _rsm_q_mulup()
capture mata: mata drop _rsm_q_dot()
capture mata: mata drop _rsm_q_compare()
capture mata: mata drop _rsm_q_scaled()
capture mata: mata drop _rsm_cert_divup()
capture mata: mata drop _rsm_cert_dual_error()
capture mata: mata drop _rsm_cert_price_bounds()
capture mata: mata drop _rsm_cert_enter()
capture mata: mata drop _rsm_precision_two_phase()
capture mata: mata drop _rsm_precision_once()
capture mata: mata drop _rsm_precision_solve()
capture mata: mata drop _rsm_precision_accept()
capture mata: mata drop _rsm_retry_quality()
capture mata: mata drop _rsm_retry_two_phase()
capture mata: mata drop _rsm_retry_once()
capture mata: mata drop _rsm_retry_realbasis()
capture mata: mata drop _rsm_retry_columns()
capture mata: mata drop _rsm_retry_solve()
capture mata: mata drop _rsm_retry_leave()
capture mata: mata drop _rsm_retry_exchange()
capture mata: mata drop _rsm_fp_two_phase()
capture mata: mata drop _rsm_fp_once()
capture mata: mata drop _rsm_fp_solve()
capture mata: mata drop _rsm_fp_enter()
capture mata: mata drop _rsm_recovery_reason()
capture mata: mata drop _rsm_audit_reset()
capture mata: mata drop _rsm_audit_emit()
capture mata: mata drop _rsm_audit_read()
capture mata: mata drop _rsm_audit_counts()
capture mata: mata drop _rsm_audit_enable()
mata:
mata set matastrict on
real scalar _dea_simplex_audit_ready() return(16)

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
    real colvector basiscols   // actual standard-form Phase-2 columns, including zero basics
    real scalar    status      // 0=optimal, 1=infeasible, 2=unbounded, 3=numerical/maxiter
    real scalar    objval      // objective function value
    real colvector x           // n x 1 primal solution
    real scalar    iters       // total simplex iterations
    real rowvector recovery    // common LP recovery contract, 10 fields
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
 * @param stage          1: radial optimum only; 2: also maximize slacks (default)
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
                                 string scalar slack_out_vars,
                                 real scalar stage,
                                 string scalar stage1_var, string scalar stage2_var,
                                 string scalar score_source_var, string scalar ref_stage_var,
                                 string scalar diagnostics_mat, string scalar linkage_mat)
{
    external real scalar _dea_rec_on
    real matrix X, Y, XR, YR, D, L
    real rowvector x_scale, y_scale, status, linkmeta
    real colvector obs_indices, result_k
    real scalar n_dmu, n_inp, n_out, k, obs_k, j
    string rowvector lv, siv, sov
    string scalar _lv_j, _sv_j
    real scalar has_lambda, has_slacks, _val_j, use_stage, raw_radial, used_radial

    use_stage = (args() >= 11 ? stage : 2)
    if (args()<12) stage1_var=""
    if (args()<13) stage2_var=""
    if (args()<14) score_source_var=""
    if (args()<15) ref_stage_var=""
    if (args()<16) diagnostics_mat=""
    if (args()<17) linkage_mat=""

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    if (diagnostics_mat!="") { XR=X; YR=Y; }
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

    D=J((diagnostics_mat!="" ? n_dmu : 0),23,.)
    L=J((linkage_mat!="" ? n_dmu : 0),30,.)
    for (k = 1; k <= n_dmu; k++) {
        obs_k = obs_indices[k]
        if (_dea_rec_on) _dea_rec_context(obs_k,1,rts)
        result_k = _dea_solve_dmu(X, Y, k, orientation, rts, tol, "", use_stage, status, raw_radial, linkmeta)
        if (linkage_mat!="") {
            // No primal stage-1 result was audited on the dual-fallback path.
            // Missing certificates stay missing; reason 12 means unavailable.
            if (missing(linkmeta[24])) {
                linkmeta[2]=(orientation=="output");linkmeta[3]=use_stage
                linkmeta[15]=0;linkmeta[16]=raw_radial;linkmeta[17]=raw_radial
                linkmeta[23]=status[2];linkmeta[24]=12;linkmeta[27]=tol
                linkmeta[28]=(rts=="crs" ? 1 : (rts=="vrs" ? 2 : (rts=="drs" ? 3 : 4)))
            }
            L[k,.]=(linkmeta,tol);L[k,1]=obs_k
        }
        if (_dea_rec_on) _dea_rec_group(obs_k,1,rts,result_k[1]<.,status[4]>0 & !hasmissing(result_k[2..(n_dmu+1)]))
        if (diagnostics_mat!="") {
            used_radial=(status[2]==0 ? raw_radial : .)
            D[k,.]=_dea_teaudit(X,Y,XR,YR,k,orientation,rts,x_scale,y_scale,
                result_k[2..(n_dmu+1)]',
                result_k[(n_dmu+2)..(n_dmu+1+n_inp)]',
                result_k[(n_dmu+2+n_inp)..(n_dmu+1+n_inp+n_out)]',
                raw_radial,used_radial,status[2],tol)
        }
        if (stage1_var!="") st_store(obs_k,st_varindex(stage1_var),status[1])
        if (stage2_var!="") st_store(obs_k,st_varindex(stage2_var),status[2])
        if (score_source_var!="") st_store(obs_k,st_varindex(score_source_var),status[3])
        if (ref_stage_var!="") st_store(obs_k,st_varindex(ref_stage_var),status[4])
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
    _dea_teaudit_store(diagnostics_mat,D)
    if (linkage_mat!="") {
        st_matrix(linkage_mat,L)
        st_matrixcolstripe(linkage_mat,(J(30,1,""),tokens("obs output stage policy gate old_basis negative old_upper old_error retry_tol retry_status new_basis new_lower new_error applied radial_before radial_after relative_change residual_lower residual_upper residual_after safe_radial stage2_status reason old_iters retry_iters effective_tol rts stage2_quality requested_tol" )'))
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
                        | string scalar trace,
                          real scalar stage)
{
    real matrix X, Y, RESULTS
    real rowvector x_scale, y_scale
    real colvector result_k
    real scalar n_dmu, n_inp, n_out, k, nc, nr_k, j, use_stage
    string scalar do_trace

    do_trace = (args() >= 7 ? trace : "")
    use_stage = (args() >= 8 ? stage : 2)

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
        result_k = _dea_solve_dmu(X, Y, k, orientation, rts, tol, do_trace, use_stage)
        nr_k = rows(result_k)
        if (nr_k >= nc) {
            RESULTS[k, .] = result_k[1..nc]'
        }
        else if (nr_k > 0) {
            RESULTS[k, 1..nr_k] = result_k'
        }
        // Return the selected stage's slacks in the source variables' units.
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
 * One- or two-stage DEA solution for DMU k with primal-dual fallback:
 *   Stage 1: Optimize efficiency (min theta or max phi) via envelopment
 *   Fallback: If Stage 1 fails, re-solve via multiplier (dual) model
 *   Stage 2 (only if requested): Maximize slacks with efficiency fixed at theta*
 * Omitted stage defaults to 2 for existing internal Mata callers.
 *
 * By strong duality the envelopment (primal) and multiplier (dual)
 * formulations yield the same optimal θ*.  When the primal LP encounters
 * numerical difficulties the dual formulation may still succeed because
 * it has a fundamentally different constraint structure (n DMU constraints
 * vs. m+s technology constraints).
 *
 * @return Column vector [theta, lam_1..n, si_1..m, so_1..s, (sc)]
 */
// T39 isolated DEA-only candidate. No public option or common-LP policy change.
void _dea_link_reset(real scalar mode,real scalar record,|real scalar fast)
{
 external real scalar _dea_link_mode,_dea_link_record,_dea_link_fast
 external real matrix _dea_link_records
 if(args()<3) fast=1
 _dea_link_fast=fast
 _dea_link_mode=mode;_dea_link_record=record;_dea_link_records=J(0,29,.)
}
// Sufficient small-residual screen on a bounded normal-arithmetic domain.
// A wide margin separates this shortcut from the later EFT trigger threshold.
// It only retains an existing result; it cannot authorize a correction.
real scalar _dea_link_small(struct RSM_Prob scalar P,real colvector x)
{
 real scalar m,n,eta
 real matrix A
 real colvector v,res,den
 m=P.m_con;n=P.n_orig
 if(m<1 | m>32 | n<1 | m!=floor(m) | n!=floor(n) | n+2*m>4096 | m*n>16384 |
  rows(P.A)!=m | cols(P.A)<n | rows(P.b)!=m | cols(P.b)!=1 |
  rows(x)!=n | cols(x)!=1 | hasmissing(x)) return(0)
 A=P.A[,1..n]
 if(hasmissing(A) | hasmissing(P.b)) return(0)
 v=(vec(A)\P.b\x);v=select(abs(v),v:!=0)
 if(rows(v)>0) if(min(v)<3.8725919148493183e-121 | max(v)>2.5822498780869086e120) return(0)
 res=P.b-A*x;den=1:+abs(P.b)+abs(A)*abs(x)
 if(hasmissing((res,den)) | min(den)<=0) return(0)
 eta=64*(m+n+4)*2.2204460492503131e-16
 return(max(abs(res):/den)<eta/4)
}
real matrix _dea_link_read()
{
 external real matrix _dea_link_records
 return(_dea_link_records)
}
void _dea_link_store(real rowvector meta,real scalar s2,|real scalar quality)
{
 external real scalar _dea_link_record
 external real matrix _dea_link_records
 meta[23]=s2
 if(args()>=3) meta[29]=quality
 if(_dea_link_record==1 & rows(_dea_link_records)<10000) _dea_link_records=_dea_link_records\meta
}
// Enclose the maximum row residual, scaled by 1+|b_i|+sum |A_ij*x_j|.
// The lower bound only triggers inspection. It never certifies feasibility.
real rowvector _dea_link_residual(struct RSM_Prob scalar P,real colvector x)
{
 real scalar i,m,n,lo,hi,denlo,denhi,maxlo,maxhi
 real rowvector Q,D
 m=P.m_con;n=P.n_orig;maxlo=maxhi=0
 if(m<1 | m>32 | n<1 | m!=floor(m) | n!=floor(n) | n+2*m>4096 | m*n>16384 |
  rows(P.A)!=m | cols(P.A)<n | rows(P.b)!=m | cols(P.b)!=1 |
  rows(x)!=n | cols(x)!=1 | hasmissing(x)) return((.,.))
 for(i=1;i<=m;i++) {
  Q=_rsm_q_dot(P.A[i,1..n],x,P.b[i])
  D=_rsm_q_dot(-abs(P.A[i,1..n]),abs(x),0)
  if(hasmissing((Q,D))) return((.,.))
  denlo=_rsm_q_add(_rsm_q_add(1,abs(P.b[i]),0),max((0,_rsm_q_add(D[1],-D[2],0))),0)
  denhi=_rsm_q_add(_rsm_q_add(1,abs(P.b[i]),1),_rsm_q_add(D[1],D[2],1),1)
  lo=max((0,_rsm_q_add(abs(Q[1]),-Q[2],0)))
  hi=_rsm_q_add(abs(Q[1]),Q[2],1)
  if(hasmissing((denlo,denhi,lo,hi)) | denlo<=0) return((.,.))
  // Upward enclosure of hi/denlo. Downward division for lo is verified below.
  hi=_rsm_cert_divup(hi,denlo)
  lo=lo/denhi
  if(lo>0) lo=_rsm_q_add(lo,-abs(lo)*4.4408920985006262e-16,0)
  if(hasmissing((lo,hi))) return((.,.))
  maxlo=max((maxlo,lo));maxhi=max((maxhi,hi))
 }
 return((maxlo,maxhi))
}
// Uses the existing inverse-residual enclosure to certify exact B^-1*b signs.
// Only a complete original-row basis is supported in this isolated candidate.
real scalar _dea_link_basis(struct RSM_Prob scalar P,struct RSM_Result scalar R,
 real matrix XP,real matrix bounds,real rowvector detail)
{
 real scalar m,n,i,delta,radius
 real rowvector Q
 real matrix B
 m=P.m_con;n=P.n_orig;XP=bounds=J(0,2,.);detail=J(1,4,.)
 if(m<1 | m>32 | n<1 | m!=floor(m) | n!=floor(n) | n+2*m>4096 | m*n>16384 |
  rows(P.A)!=m | cols(P.A)<n | rows(P.b)!=m | cols(P.b)!=1 |
  rows(R.basiscols)!=m | cols(R.basiscols)!=1 | hasmissing(R.basiscols)) return(0)
 for(i=1;i<=m;i++) if(R.basiscols[i]<1 | R.basiscols[i]>n |
  R.basiscols[i]!=floor(R.basiscols[i]) | sum(R.basiscols:==R.basiscols[i])!=1) return(0)
 B=P.A[,R.basiscols];XP=_rsm_retry_refine(B,P.b)
 if(hasmissing(XP)) return(0)
 delta=_rsm_cert_dual_error(B,P.b,XP,detail)
 if(missing(delta)) return(0)
 bounds=J(m,2,.)
 for(i=1;i<=m;i++) {
  Q=_rsm_q_dot((-1,-1),(XP[i,1]\XP[i,2]),0)
  radius=_rsm_q_add(Q[2],delta,1)
  bounds[i,1]=_rsm_q_add(Q[1],-radius,0)
  bounds[i,2]=_rsm_q_add(Q[1],radius,1)
 }
 return(!hasmissing(bounds))
}
// meta: DMU,orientation,stage,mode,gate,old_basis_available,negative_count,
// min_old_upper,old_error,retry_tol,retry_status,new_basis_available,
// min_new_lower,new_error,accepted,old_radial,new_radial,relative_change,
// old_residual_lower,old_residual_upper,new_residual_upper,safe_radial,
// stage2_status,reason,old_iterations,retry_iterations,effective_stage2_tol,RTS,
// stage2_quality_outcome (repaired stage1 only).
// Modes 0=off,1=diagnose/retain,2=bounded repair.
// Reasons 0=unchanged,1=range,2=small residual,3=uncertified,4=no proven negative,
// 5=audit only,6=retry failed,7=new basis uncertified,8=radial/budget,
// 9=quality/residual,10=repaired,11=fast small-residual screen.
// No uncertain result is silently rejected.
void _dea_link_apply(struct RSM_Prob scalar P,struct DEA_Info scalar I,
 struct RSM_Result scalar R,real scalar k,real scalar stage,real scalar tol,
 real scalar linktol,real rowvector meta)
{
 external real scalar _dea_link_mode,_dea_link_fast,_dea_rec_on
 struct RSM_Prob scalar trial
 struct RSM_Result scalar T
 real matrix XP,bounds
 real rowvector detail,res,check,q
 real colvector x
 real scalar mode,gate,oldrad,newrad,safe,i,pos,accepted
 mode=_dea_link_mode;linktol=tol;oldrad=_dea_extract_theta(R,I)
 meta=(k,I.ort=="output",stage,mode,0,0,0,.,.,.,.,0,.,.,0,oldrad,oldrad,0,.,.,.,.,.,0,R.iters,0,tol,(I.rts=="crs" ? 1 : (I.rts=="vrs" ? 2 : (I.rts=="drs" ? 3 : 4))),.)
 if(mode==0) return
 if(R.status!=0 | rows(R.x)!=P.n_orig | hasmissing(R.x) | missing(oldrad) | oldrad<=0) {meta[24]=1;return;}
 if(_dea_link_fast==1) if(_dea_link_small(P,R.x)) {meta[24]=11;return;}
 res=_dea_link_residual(P,R.x);meta[19..20]=res
 if(hasmissing(res)) {meta[24]=1;return;}
 gate=64*(P.m_con+P.n_orig+4)*2.2204460492503131e-16
 if(res[1]<=gate) {meta[24]=2;return;}
 meta[5]=1
 if(!_dea_link_basis(P,R,XP,bounds,detail)) {meta[24]=3;return;}
 meta[6]=1;meta[7]=sum(bounds[,2]:<0);meta[8]=min(bounds[,2]);meta[9]=detail[2]
 if(meta[7]==0) {meta[24]=4;return;}
 if(mode!=2) {meta[24]=5;return;}
 linktol=min((tol,1e-10));meta[10]=linktol
 trial=P
 // Exactly one fresh top-level retry; existing solver-internal caps stay intact.
 T=_rsm_two_phase(trial,linktol,"");meta[11]=T.status;meta[26]=T.iters
 if(_dea_rec_on) _dea_rec_status(T.status)
 if(T.status!=0) {linktol=tol;meta[24]=6;return;}
 if(!_dea_link_basis(P,T,XP,bounds,detail)) {linktol=tol;meta[24]=7;return;}
 meta[12]=1;meta[13]=min(bounds[,1]);meta[14]=detail[2]
 if(min(bounds[,1])<0) {linktol=tol;meta[24]=7;return;}
 pos=0
 for(i=1;i<=rows(T.basiscols);i++) if(T.basiscols[i]==I.idx_theta) pos=i
 if(pos==0) {linktol=tol;meta[24]=8;return;}
 newrad=_dea_extract_theta(T,I)
 safe=(I.ort=="output" ? min((newrad,bounds[pos,1])) : max((newrad,bounds[pos,2])))
 meta[22]=safe
 if(missing(safe) | safe<=0 | abs(safe-oldrad)/(1+abs(oldrad))>tol) {linktol=tol;meta[24]=8;return;}
 x=T.x;x[I.idx_theta]=safe
 check=_dea_lpcheck(x',P.A[,1..P.n_orig],P.b,J(0,P.n_orig,0),J(0,1,0),J(1,P.n_orig,0),J(1,P.n_orig,.))
 res=_dea_link_residual(P,x);meta[21]=res[2]
 T.x=x;T.objval=P.c[1..P.n_orig]*x
 // Original LP primal/dual quality is checked again after the radial change.
 q=_rsm_success_check(P,T,T.basiscols,linktol)
 if(hasmissing(res) | res[2]>linktol | !_dea_lpcheck_ok(check,linktol) |
  q[1]!=1 | hasmissing(T.x) | min(T.x)<0 | missing(T.objval)) {linktol=tol;meta[24]=9;return;}
 T.iters=R.iters+T.iters
 R=T;R.msg="Optimal after verified DEA stage-link retry"
 meta[15]=1;meta[17]=_dea_extract_theta(R,I);meta[18]=abs(meta[17]-oldrad)/(1+abs(oldrad));meta[24]=10;meta[27]=linktol
}
// A repaired radial result may expose a Stage-2 reference only after the
// fixed-target original LP passes both primal residual and primal/dual quality.
real scalar _dea_link_stage2(struct RSM_Prob scalar P,struct RSM_Result scalar R,
 real scalar tol)
{
 real rowvector q,check
 if(R.status!=0) return(.)
 q=_rsm_success_check(P,R,R.basiscols,tol)
 check=_dea_lpcheck(R.x',P.A[,1..P.n_orig],P.b,J(0,P.n_orig,0),J(0,1,0),J(1,P.n_orig,0),J(1,P.n_orig,.))
 if(missing(q[1]) | q[1]<=0 | !_dea_lpcheck_ok(check,tol) | hasmissing(R.x) |
  min(R.x)<0 | missing(R.objval)) {
  R.status=3;R.objval=.;R.x=J(P.n_orig,1,.)
  R.msg="DEA stage-link fixed-target certificate unavailable"
 }
 return(q[1])
}

real colvector _dea_solve_dmu(real matrix X,
                               real matrix Y,
                               real scalar k,
                               string scalar ort,
                               string scalar rts,
                               real scalar tol,
                               string scalar trace,
                               | real scalar stage, real rowvector status, real scalar raw_radial, real rowvector linkmeta)
{
    external real scalar _dea_rec_on
    real scalar n_dmu, n_inp, n_out, nc, theta_star, use_stage
    real colvector x0, y0, result
    real scalar linktol, linkquality

    struct RSM_Prob scalar prob1, prob2, link_original, link_second
    struct DEA_Info scalar info1, info2
    struct RSM_Result scalar res1, res2

    // Primal-dual fallback variables
    struct RSM_Prob scalar prob_m
    struct MULT_Info scalar minfo
    struct RSM_Result scalar res_m

    // Optional output: radial availability, stage-2 code, score source, reference stage.
    status = (1, ., 0, 0)
    raw_radial = .
    linkmeta=J(1,29,.)
    use_stage = (args() >= 8 ? stage : 2)
    if (use_stage != 1 & use_stage != 2) _error(3300)

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
    link_original=prob1

    if (_dea_rec_on) _dea_rec_phase(1,1)

    res1 = _rsm_two_phase(prob1, tol, trace)

    if (_dea_rec_on) _dea_rec_status(res1.status)

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
        if (_dea_rec_on) _dea_rec_phase(1,2)
        res_m = _rsm_two_phase(prob_m, tol, trace)
        if (_dea_rec_on) _dea_rec_status(res_m.status)

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

        raw_radial=theta_star
        status[1]=0; status[3]=2
        // The dual optimum supplies a score, not a primal reference solution.
        // A single-stage request must not run a slack-maximization LP to recover it.
        if (use_stage == 1) {
            result = J(nc, 1, .)
            if (ort == "output" && theta_star != 0) result[1] = 1 / theta_star
            else result[1] = theta_star
            return(result)
        }

        // Attempt Stage 2 with θ* from multiplier to recover λ and slacks
        _dea_lp_stage2(prob2, info2, X, Y, x0, y0, ort, rts, theta_star)
        if (_dea_rec_on) _dea_rec_phase(2,1)
        res2 = _rsm_two_phase(prob2, tol, trace)
        if (_dea_rec_on) _dea_rec_status(res2.status)
        status[2]=res2.status

        if (res2.status != 0) {
            // Stage 2 failed: retain the dual score, but no reference was recovered.
            if (trace == "trace") {
                printf("{txt}DMU %g: Stage2 after fallback failed; " +
                       "returning theta only\n", k)
            }
            result = J(nc, 1, .)
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
        status[4]=2
        result = _dea_extract(res2, info2, theta_star, 2, prob2, tol)
        return(result)
    }

    _dea_link_apply(link_original,info1,res1,k,use_stage,tol,linktol,linkmeta)
    theta_star = _dea_extract_theta(res1, info1)
    raw_radial=theta_star

    if (trace == "trace") {
        printf("{txt}DMU %g Stage1: theta* = %12.8f (iter=%g)\n",
               k, theta_star, res1.iters)
    }

    status[1]=0; status[3]=1
    if (use_stage == 1) {
        status[4]=1
        _dea_link_store(linkmeta,.)
        return(_dea_extract(res1, info1, theta_star, 1, link_original, linktol))
    }

    // =================================================================
    // STAGE 2: Maximize slacks with theta fixed
    // =================================================================
    _dea_lp_stage2(prob2, info2, X, Y, x0, y0, ort, rts, theta_star)
    linkquality=.
    if(linkmeta[15]==1) link_second=prob2

    if (_dea_rec_on) _dea_rec_phase(2,1)

    res2 = _rsm_two_phase(prob2, linktol, trace)
    if(linkmeta[15]==1) linkquality=_dea_link_stage2(link_second,res2,linktol)

    if (_dea_rec_on) _dea_rec_status(res2.status)
    status[2]=res2.status

    if (res2.status != 0) {
        if (trace == "trace") {
            printf("{err}DMU %g Stage2 FAIL: %s (returning radial score only)\n",
                   k, res2.msg)
        }
        // Preserve the radial score independently of reference availability.
        result = _dea_extract(res1, info1, theta_star, 1, link_original, linktol)
        result[2..nc] = J(nc-1,1,.)
    }
    else {
        if (trace == "trace") {
            printf("{txt}DMU %g Stage2: sum_slacks = %12.8f (iter=%g)\n",
                   k, -res2.objval, res2.iters)
        }
        status[4]=2
        result = _dea_extract(res2, info2, theta_star, 2, prob2, linktol)
    }

    _dea_link_store(linkmeta,status[2],linkquality)
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
struct RSM_Result scalar _rsm_eta_two_phase(struct RSM_Prob scalar prob,
                                         real scalar tol,
                                         string scalar trace,| real colvector final_basic)
{
    struct RSM_Prob scalar original, retry
    struct RSM_Result scalar first, second
    real scalar i, scale
    real rowvector check, column_scale

    original=prob
    first=_rsm_two_phase_once(prob,tol,trace,0,final_basic)
    if (first.status!=3) return(first)

    // Bounded numerical recovery, inspired by Benchmarking's fresh-basis
    // retry and equilibration. Rebuild from the untouched LP, never from
    // a partially pivoted state. Positive row/objective scaling is equivalent;
    // variable units are mapped back and the user's radial target is unchanged.
    retry=original
    for (i=1;i<=retry.m_con;i++) {
        scale=max((abs(retry.A[i,1..retry.n_orig]),abs(retry.b[i])))
        if (scale>0 & scale<.) {
            retry.A[i,.]=retry.A[i,.]/scale
            retry.b[i]=retry.b[i]/scale
        }
    }
    // Artificial columns are temporary Phase-I coordinates. Restore their
    // unit columns so Phase I minimizes violations in the equilibrated rows,
    // rather than retaining weights from the original measurement units.
    for (i=retry.n_orig+1;i<=retry.n_var;i++) {
        if (retry.is_artif[i]) retry.A[.,i]=(retry.A[.,i]:!=0)
    }
    column_scale=J(1,retry.n_orig,1)
    for (i=1;i<=retry.n_orig;i++) {
        scale=max(abs(retry.A[.,i]))
        if (scale>0 & scale<.) {
            column_scale[i]=1/scale
            retry.A[.,i]=retry.A[.,i]/scale
            retry.c[i]=retry.c[i]/scale
        }
    }
    scale=max(abs(retry.c))
    if (scale>0 & scale<.) retry.c=retry.c/scale
    if (hasmissing(retry.A) | hasmissing(retry.b) | hasmissing(retry.c) | hasmissing(column_scale)) return(first)
    if (trace=="trace") printf("{txt}RSM restart: fresh basis, row/column/objective equilibration, Bland\n")
    second=_rsm_two_phase_once(retry,tol,trace,1,final_basic)
    second.iters=first.iters+second.iters
    if (second.status==0) {
        second.x=second.x:*column_scale'
        check=_dea_lpcheck(second.x',original.A[.,1..original.n_orig],original.b,
            J(0,original.n_orig,0),J(0,1,0),J(1,original.n_orig,0),J(1,original.n_orig,.))
        second.objval=original.c[1..original.n_orig]*second.x
        if (_dea_lpcheck_ok(check,tol) & second.objval<.) {
            second.msg="Optimal after equilibrated restart"
            prob=original
            return(second)
        }
    }
    // A numerical restart is not an infeasibility/unboundedness certificate.
    first.iters=second.iters
    return(first)
}

struct RSM_Result scalar _rsm_two_phase_once(struct RSM_Prob scalar prob,
                                         real scalar tol,
                                         string scalar trace,
                                         real scalar restart,| real colvector final_basic)
{
    real scalar i, max_iter, p1_iter
    real colvector xb
    real rowvector c_orig
    real scalar art_sum
    real matrix B0, original_A
    real colvector original_b, feasible_x
    real rowvector check

    struct RSM_Basis scalar basis
    struct RSM_Result scalar result

    final_basic=J(0,1,.)
    max_iter = max((10000, 20 * prob.m_con))

    // Save Phase 2 objective
    c_orig = prob.c
    original_A = prob.A[.,1..prob.n_orig]
    original_b = prob.b

    // =================================================================
    // PHASE 1: Minimize sum of artificials
    // =================================================================
    prob.c = prob.is_artif

    // Initialize basis
    basis.basic = prob.init_basic
    basis.n_upd = 0

    // Initial B^{-1} (should be identity by construction)
    B0 = prob.A[., basis.basic]
    if (all(B0 :== I(prob.m_con))) {
        basis.Binv = I(prob.m_con)
    }
    else {
        basis.Binv = _rsm_inv(B0)
        if (hasmissing(basis.Binv)) {
            result.status = 3
            result.msg = "Invalid initial basis inverse"
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
    result = _rsm_solve(prob, basis, tol, max_iter, 1, trace, restart)
    p1_iter = result.iters

    if (result.status == 2) {
        result.status = 3
        result.msg = "Phase1 unbounded; no feasibility certificate"
        prob.c = c_orig
        return(result)
    }

    if (result.status != 0) {
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

    // A tiny artificial objective alone is not a feasible original LP.
    feasible_x=J(prob.n_orig,1,0)
    for (i=1;i<=prob.m_con;i++) {
        if (basis.basic[i]<=prob.n_orig) feasible_x[basis.basic[i]]=max((0,xb[i]))
    }
    check=_dea_lpcheck(feasible_x',original_A,original_b,
        J(0,prob.n_orig,0),J(0,1,0),J(1,prob.n_orig,0),J(1,prob.n_orig,.))
    if (!_dea_lpcheck_ok(check,tol)) {
        result.status=3; result.msg="Residual check failed"
        result.x=feasible_x; result.objval=.
        prob.c=c_orig
        return(result)
    }

    // =================================================================
    // TRANSITION: Purge artificials, drop columns, restore objective
    // =================================================================
    if (!_rsm_purge_artif(prob, basis, tol)) {
        result.status=3; result.msg="Artificial-basis pivot failed"
        result.objval=.; result.x=feasible_x
        return(result)
    }
    if (!_rsm_drop_artif(prob, basis)) {
        result.status=3; result.msg="Artificial-basis transition failed"
        result.objval=.; result.x=feasible_x
        return(result)
    }
    prob.c = c_orig[1..prob.n_var]

    // Refactorize for clean Phase 2 start
    if (!_rsm_refac(basis, prob)) {
        result.status=3; result.msg="Invalid Phase2 basis inverse"
        result.objval=.; result.x=feasible_x
        return(result)
    }

    // =================================================================
    // PHASE 2: Optimize original objective
    // =================================================================
    if (trace == "trace") {
        printf("{txt}Phase2: %g con, %g var\n", prob.m_con, prob.n_var)
    }

    result = _rsm_solve(prob, basis, tol, max_iter - p1_iter, 2, trace, restart)
    result.iters = result.iters + p1_iter

    // Extract full solution
    if (result.status == 0) {
        result.x = J(prob.n_var, 1, 0)
        xb = basis.Binv * prob.b
        for (i = 1; i <= prob.m_con; i++) {
            result.x[basis.basic[i]] = max((0, xb[i]))
        }
        result.objval = prob.c * result.x
        check=_dea_lpcheck(result.x',original_A,original_b,
            J(0,prob.n_var,0),J(0,1,0),J(1,prob.n_var,0),J(1,prob.n_var,.))
        if (!_dea_lpcheck_ok(check,tol) | result.objval>=.) {
            result.status=3; result.msg="Residual check failed"; result.objval=.
        }
    }

    if(result.status==0) final_basic=basis.basic
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
                                     string scalar trace,
                                     real scalar restart)
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

    use_bland = restart
    consec_degen = 0
    stall_count = 0
    prev_obj = .

    for (iter = 1; iter <= max_iter; iter++) {
        result.iters = iter

        // Periodic refactorization
        if (basis.n_upd > 0 && mod(basis.n_upd, _dea_refac_freq()) == 0) {
            if (!_rsm_refac(basis, prob)) {
                result.msg="Basis refactorization failed"
                return(result)
            }
        }

        // Current BFS: xb = B^{-1} * b
        xb = basis.Binv * prob.b
        if (hasmissing(xb)) { result.msg="Nonfinite basic solution"; return(result); }
        if (min(xb)<-min((_dea_tol_feas(),tol))) {
            if (!_rsm_refac(basis,prob)) {
                result.msg="Basis refactorization failed"; return(result)
            }
            xb=basis.Binv*prob.b
            if (hasmissing(xb) | min(xb)<-min((_dea_tol_feas(),tol))) {
                result.msg="Primal basis lost feasibility"; return(result)
            }
        }

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

        // Confirm termination using a freshly factored basis.
        if (entering==0 & basis.n_upd>0) {
            if (!_rsm_refac(basis,prob)) {
                result.msg="Basis refactorization failed"; return(result)
            }
            xb=basis.Binv*prob.b
            if (hasmissing(xb) | min(xb)<-min((_dea_tol_feas(),tol))) {
                result.msg="Primal basis lost feasibility"; return(result)
            }
            entering=_rsm_enter(prob,basis,use_bland,tol)
            obj_val=prob.c[basis.basic']*xb
        }

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
        if (hasmissing(d)) { result.msg="Nonfinite pivot column"; return(result); }

        // =============================================================
        // RATIO TEST: select pivot row (Step 3)
        // =============================================================
        pivot_row = _rsm_leave(xb, d, basis, use_bland, tol)

        if (pivot_row == 0) {
            if (basis.n_upd>0) {
                if (!_rsm_refac(basis,prob)) {
                    result.msg="Basis refactorization failed"; return(result)
                }
                continue
            }
            if (rows(d)>0 && max(d)>0) {
                result.msg="No stable positive pivot"; return(result)
            }
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
        if (trace=="trace") printf("{txt}  pivot P%g i%g enter%g leave%g step%e pivot%e Bland%g\n",
            phase,iter,entering,basis.basic[pivot_row],max((0,xb[pivot_row]))/d[pivot_row],d[pivot_row],use_bland)
        if (!_rsm_pivot(basis, d, pivot_row)) {
            result.msg="Basis pivot failed"; return(result)
        }
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

    // Price small objectives in their own units. A large coefficient must
    // not loosen the absolute cutoff for other, small-cost columns.
    opt_tol = min((_dea_tol_opt(),tol))*min((1,max(abs(prob.c))))

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
            if (rc_all[i] == best_rc) {
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
 *   Pass 1: relaxed bound = min{ max(0,xb_i+feas_tol) / d_i }.
 *   Pass 2: choose the largest stable pivot whose exact ratio is within it.
 * Bland uses only exact minimum-ratio ties, never the Harris candidate set.
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
    real scalar pvt_tol, feas_tol

    m = rows(d)
    pvt_tol = _dea_tol_pivot()
    feas_tol = min((_dea_tol_feas(),tol))

    if (max(d) <= pvt_tol) return(0)

    // Pass 1: compute ratios
    ratios = J(m, 1, .)
    harris_bound=.
    for (i = 1; i <= m; i++) {
        if (d[i] > pvt_tol) {
            if (xb[i] >= -feas_tol) {
                ratios[i] = max((0, xb[i])) / d[i]
                harris_bound=min((harris_bound,max((0,xb[i]+feas_tol))/d[i]))
            }
        }
    }

    min_ratio = min(ratios)
    if (min_ratio == .) return(0)

    // Pass 2: select among eligible
    if (use_bland) {
        pivot_row = 0
        min_idx = .
        for (i = 1; i <= m; i++) {
            if (ratios[i] < . && ratios[i] == min_ratio) {
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
real scalar _rsm_pivot(struct RSM_Basis scalar basis,
                 real colvector d,
                 real scalar pivot_row)
{
    real scalar m, i, pv
    real rowvector prow

    m = rows(d)
    pv = d[pivot_row]

    if (hasmissing(d) | abs(pv) <= _dea_tol_pivot()) return(0)

    // Save pivot row before modification
    prow = basis.Binv[pivot_row, .]

    // Update non-pivot rows
    for (i = 1; i <= m; i++) {
        if (i != pivot_row && d[i] != 0) {
            basis.Binv[i, .] = basis.Binv[i, .] - (d[i] / pv) * prow
        }
    }

    // Scale pivot row
    basis.Binv[pivot_row, .] = prow / pv

    basis.n_upd = basis.n_upd + 1
    return(!hasmissing(basis.Binv))
}


// ============================================================================
// REFACTORIZE: Recompute B^{-1} from scratch
// ============================================================================
/**
 * Corrects accumulated round-off from eta updates.
 * Strategy: LU first, SVD fallback if ill-conditioned.
 */
real scalar _rsm_refac(struct RSM_Basis scalar basis,
                 struct RSM_Prob scalar prob)
{
    real matrix B, Binv_new

    B = prob.A[., basis.basic]
    Binv_new = _rsm_inv(B)

    if (!hasmissing(Binv_new)) {
        basis.Binv = Binv_new
        basis.n_upd = 0
        return(1)
    }
    basis.Binv=J(prob.m_con,prob.m_con,.)
    return(0)
}


// ============================================================================
// SAFE MATRIX INVERSE
// ============================================================================
/**
 * Equilibrated LU, bounded residual correction, then checked SVD fallback.
 * A finite pseudo-inverse is insufficient: both inverse identities must hold.
 */
real matrix _rsm_inv(real matrix B)
{
    real scalar m, i, j, k, method, err
    real matrix Binv, C, CI, correction, left, right
    real colvector rs
    real rowvector cs

    m = rows(B)
    if (m == 0) return(J(0, 0, .))

    if (cols(B)!=m | hasmissing(B)) return(J(m,m,.))
    if (all(B :== I(m))) return(I(m))
    rs=rowmax(abs(B))
    if (any(rs:==0)) return(J(m,m,.))
    C=B:/rs
    cs=colmax(abs(C))
    if (any(cs:==0)) return(J(m,m,.))
    C=C:/cs
    for (method=1;method<=2;method++) {
        CI=(method==1 ? lusolve(C,I(m)) : svsolve(C,I(m),1e-14))
        if (hasmissing(CI)) continue
        for (k=0;k<=2;k++) {
            Binv=CI
            for (i=1;i<=m;i++) for (j=1;j<=m;j++) Binv[i,j]=CI[i,j]/cs[i]/rs[j]
            if (!hasmissing(Binv)) {
                if (_rsm_inverse_fast(B,Binv)) return(Binv)
                left=_rsm_inverse_residual(B,Binv,I(m))
                right=_rsm_inverse_residual(Binv,B,I(m))
                if (!hasmissing(left) & !hasmissing(right)) {
                    err=max((max(abs(left)),max(abs(right))))
                    if (err<. & err<=1e-8) return(Binv)
                }
            }
            if (k==2) break
            correction=lusolve(C,I(m)-C*CI)
            if (hasmissing(correction)) break
            CI=CI+correction
        }
    }

    return(J(m, m, .))
}


// ============================================================================
// PURGE ARTIFICIALS FROM BASIS
// ============================================================================
/**
 * After Phase 1, artificial vars may remain in basis at value 0.
 * Try to pivot each out using a non-artificial non-basic variable.
 */
real scalar _rsm_purge_artif(struct RSM_Prob scalar prob,
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
            if (!_rsm_pivot(basis, d, i)) return(0)
            basis.basic[i] = best_j
        }
    }
    return(1)
}


// ============================================================================
// DROP ARTIFICIAL COLUMNS
// ============================================================================
/**
 * Removes artificial columns from problem after Phase 1.
 */
real scalar _rsm_drop_artif(struct RSM_Prob scalar prob,
                      struct RSM_Basis scalar basis)
{
    real scalar new_n, i
    real colvector keep, transformed_b
    real matrix transformed_A, error_scale

    if (prob.n_var == prob.n_orig) return(1)

    new_n = prob.n_orig
    if (any(basis.basic:>new_n)) {
        // A remaining zero artificial represents a dependent tableau row.
        // Remove only rows with no real coefficient beyond multiplication
        // roundoff. Never fill a basis slot with an arbitrary unused column.
        transformed_A=basis.Binv*prob.A[.,1..new_n]
        transformed_b=basis.Binv*prob.b
        error_scale=abs(basis.Binv)*abs(prob.A[.,1..new_n])
        for (i=1;i<=prob.m_con;i++) {
            if (basis.basic[i]<=new_n) continue
            if (any(abs(transformed_A[i,.]):>64*epsilon(1)*prob.m_con:*error_scale[i,.])) return(0)
        }
        keep=select((1::prob.m_con),basis.basic:<=new_n)
        prob.A=transformed_A[keep,.]
        prob.b=transformed_b[keep]
        basis.basic=basis.basic[keep]
        prob.m_con=rows(keep)
        basis.Binv=I(prob.m_con)
        basis.n_upd=0
    }
    prob.A = prob.A[., 1..new_n]
    prob.c = prob.c[1..new_n]
    prob.is_artif = J(1, new_n, 0)
    prob.n_var = new_n

    return(1)
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
                             real scalar stage,
                             | struct RSM_Prob prob, real scalar tol)
{
    real scalar nc, i
    real colvector result, raw, candidate

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

    raw=result

    // Lambda
    for (i = 1; i <= info.n_dmu; i++) {
        if (info.idx_lam[i] <= rows(res.x)) {
            raw[1 + i] = res.x[info.idx_lam[i]]
            result[1 + i] = max((0, raw[1 + i]))
        }
    }

    // Input slacks
    for (i = 1; i <= info.n_inp; i++) {
        if (info.idx_si[i] <= rows(res.x)) {
            raw[1 + info.n_dmu + i] = res.x[info.idx_si[i]]
            result[1 + info.n_dmu + i] = max((0, raw[1 + info.n_dmu + i]))
        }
    }

    // Output slacks
    for (i = 1; i <= info.n_out; i++) {
        if (info.idx_so[i] <= rows(res.x)) {
            raw[1 + info.n_dmu + info.n_inp + i] = res.x[info.idx_so[i]]
            result[1 + info.n_dmu + info.n_inp + i] = max((0, raw[1 + info.n_dmu + info.n_inp + i]))
        }
    }

    // Convexity slack
    if (info.idx_sc > 0 && info.idx_sc <= rows(res.x)) {
        raw[nc] = res.x[info.idx_sc]
        result[nc] = max((0, raw[nc]))
    }

    // Never round the radial score to zero. Retain the legacy presentation
    // cleanup only if it preserves the original LP's requested residual test.
    for (i = 2; i <= nc; i++) {
        if (abs(result[i]) < _dea_tol_feas()) result[i] = 0
    }

    if (args()>=6) {
        candidate=res.x[1..prob.n_orig]
        candidate[info.idx_lam]=result[2..(1+info.n_dmu)]
        candidate[info.idx_si]=result[(2+info.n_dmu)..(1+info.n_dmu+info.n_inp)]
        candidate[info.idx_so]=result[(2+info.n_dmu+info.n_inp)..(1+info.n_dmu+info.n_inp+info.n_out)]
        if (info.idx_sc>0) candidate[info.idx_sc]=result[nc]
        if (!_dea_lpcheck_ok(_dea_lpcheck(candidate',prob.A[.,1..prob.n_orig],prob.b,
            J(0,prob.n_orig,0),J(0,1,0),J(1,prob.n_orig,0),J(1,prob.n_orig,.)),tol)) {
            // The solver already checked res.x. Return that accepted vector,
            // not a rounded vector that breaks its constraints. No new LP solve.
            result=raw
        }
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

    // Preserve the solver-accepted multiplier solution. A small positive
    // weight can be essential to normalization when observations are large.

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
                                        real scalar tol,| real scalar check_success)
{
    struct RSM_Prob scalar prob
    struct RSM_Result scalar res
    real scalar nv

    if(args()<7) check_success=1
    nv = cols(c)

    _rsm_std_form(prob, A_iq, b_iq, A_eq, b_eq, c)
    res = _rsm_two_phase(prob, tol, "",check_success)

    // Trim solution to original variables only
    if (res.status == 0 && rows(res.x) >= nv) {
        res.x = res.x[1..nv]
    }
    else if (res.status == 0) {
        res.x = J(nv, 1, 0)
    }

    return(res)
}

// Non-TE models also validate the submitted rows before exporting a score.
struct RSM_Result scalar _dea_nonte_simplex(real matrix A_iq,
                                        real colvector b_iq,
                                        real matrix A_eq,
                                        real colvector b_eq,
                                        real rowvector c,
                                        real scalar tol)
{
    external real scalar _dea_rec_on
    struct RSM_Result scalar res
    real scalar nv

    nv = cols(c)

    res = _rsm_solve_lp(A_iq,b_iq,A_eq,b_eq,c,tol)

    // Trim solution to original variables only
    if (res.status == 0 && rows(res.x) >= nv) {
        res.x = res.x[1..nv]
    }
    else if (res.status == 0) {
        res.status = 4
        res.objval = .
        res.x = J(nv, 1, .)
    }

    if (res.status == 0) {
        if (!_dea_lpcheck_ok(_dea_lpcheck(res.x',A_eq,b_eq,A_iq,b_iq,
            J(1,nv,0),J(1,nv,.)),tol) | res.objval>=. | c*res.x>=.) {
            res.status=4
            res.objval=.
            res.x=J(nv,1,.)
        }
    }
    if (_dea_rec_on) _dea_rec_status(res.status)
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
                        string scalar rts,
    | real scalar tol)
{
    real matrix X, Y, W
    real colvector obs_indices, x0, y0, w0
    real rowvector c
    real matrix A_iq, A_eq
    real colvector b_iq, b_eq
    real scalar n, m, s, k, obs_k
    real scalar min_cost, obs_cost, cost_eff

    struct RSM_Result scalar res

    if (args()<7) tol=1e-7

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

        res = _dea_nonte_simplex(A_iq, b_iq, A_eq, b_eq, c, tol)

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
                           string scalar rts,
    | real scalar tol)
{
    real matrix X, Y, P
    real colvector obs_indices, x0, y0, p0
    real rowvector c
    real matrix A_iq, A_eq
    real colvector b_iq, b_eq
    real scalar n, m, s, k, obs_k
    real scalar max_rev, obs_rev, rev_eff

    struct RSM_Result scalar res

    if (args()<7) tol=1e-7

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

        res = _dea_nonte_simplex(A_iq, b_iq, A_eq, b_eq, c, tol)

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
                          string scalar rts,
    | real scalar tol)
{
    real matrix X, Y, W, P
    real colvector obs_indices, x0, y0, w0, p0
    real rowvector c
    real matrix A_iq, A_eq
    real colvector b_iq, b_eq
    real scalar n, m, s, k, obs_k
    real scalar max_profit, obs_profit, profit_eff

    struct RSM_Result scalar res

    if (args()<8) tol=1e-7

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

        res = _dea_nonte_simplex(A_iq, b_iq, A_eq, b_eq, c, tol)

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
                       string scalar rts,
    | real scalar tol)
{
    real matrix X, Y, GX, GY
    real colvector obs_indices, x0, y0, gx0, gy0
    real rowvector c, sol
    real matrix A_iq, A_eq
    real colvector b_iq, b_eq
    real scalar n, m, s, k, obs_k, beta

    struct RSM_Result scalar res

    if (args()<8) tol=1e-7

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

        res = _dea_nonte_simplex(A_iq, b_iq, A_eq, b_eq, c, tol)

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
 *              B'λ + β·gb = b₀,    (bad outputs: equality, no reduction slack)
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
                                   string scalar rts,
    | real scalar tol)
{
    real matrix X, Y, B
    real colvector obs_indices, x0, y0, b0, gx0, gy0, gb0
    real rowvector c, lb, ub, sol, xscale, yscale, bscale
    real matrix A_iq, A_eq
    real colvector b_iq, b_eq
    real scalar n, m, s, q, k, obs_k, beta, rc
    class LinearProgram scalar S
    string scalar dir

    if (args()<8) tol=1e-7

    X = st_data(., tokens(inputs),      touse)
    Y = st_data(., tokens(outputs),     touse)
    B = st_data(., tokens(bad_outputs), touse)
    obs_indices = selectindex(st_data(., touse))

    n = rows(X)
    m = cols(X)
    s = cols(Y)
    q = cols(B)

    dir = (direction == "" ? "yb" : direction)

    // Scale the data and their observed directions together; beta and lambda
    // stay in their original dimensionless units. Keep zero columns unchanged.
    xscale = colmax(abs(X)); xscale = xscale :+ (xscale:==0)
    yscale = colmax(abs(Y)); yscale = yscale :+ (yscale:==0)
    bscale = colmax(abs(B)); bscale = bscale :+ (bscale:==0)
    X = X :/ xscale
    Y = Y :/ yscale
    B = B :/ bscale

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

        // Good outputs expand: Y'lambda >= y0 + beta*gy0.
        A_iq = (X', gx0) \ (-Y', gy0)
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

        rc = _dea_nonte_optimize(S,tol,1e-12)
        if (rc == 0) {
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
                            string scalar rts,
    | real scalar tol)
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

    if (args()<8) tol=1e-7

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

        res = _dea_nonte_simplex(A_iq, b_iq, A_eq, b_eq, c, tol)

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
                       string scalar orientation,
    | real scalar tol)
{
    external real scalar _dea_rec_on
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

    if (args()<9) tol=1e-7

    real rowvector x_scale, y_scale
    real scalar mp, sp

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    // Column rescaling preserves the SBM ratios and balances the LP.
    x_scale = colmax(X); y_scale = colmax(Y)
    for (j=1; j<=cols(X); j++) if (x_scale[j]==0) x_scale[j]=1
    for (j=1; j<=cols(Y); j++) if (y_scale[j]==0) y_scale[j]=1
    X = X :/ x_scale; Y = Y :/ y_scale
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
        if (_dea_rec_on) _dea_rec_context(obs_k,2,rts)
        x0 = X[k, .]'; y0 = Y[k, .]'
        sbm_eff = .
        slack_in = J(1, m, 0); slack_out = J(1, s, 0)

        // Coefficients
        // Tone (2001), section 6: omitted zero terms also reduce the mean.
        // Zero outputs use the structural-zero interpretation (no perturbation).
        mp = sum(x0:>0); sp = sum(y0:>0)
        if (mp==0 | sp==0) {
            st_store(obs_k, st_varindex(eff_var), .)
            continue
        }
        coef_sin = J(1, m, 0)
        for (j = 1; j <= m; j++) {
            if (x0[j] > 0) coef_sin[j] = -1 / (mp * x0[j])
        }
        coef_sout = J(1, s, 0)
        for (j = 1; j <= s; j++) {
            if (y0[j] > 0) coef_sout[j] = 1 / (sp * y0[j])
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

            res = _dea_nonte_simplex(A_iq, b_iq, A_eq, b_eq, c, tol)

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

            res = _dea_nonte_simplex(A_iq, b_iq, A_eq, b_eq, c, tol)

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

        // Public/internal slack requests remain in the original data units.
        slack_in = slack_in :* x_scale
        slack_out = slack_out :* y_scale
        if (_dea_rec_on) _dea_rec_group(obs_k,2,rts,sbm_eff<.,.)
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
                            string scalar rts,
    | real scalar tol)
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

    if (args()<8) tol=1e-7

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

            res = _dea_nonte_simplex(A_iq, b_iq, A_eq, b_eq, c, tol)

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

            res = _dea_nonte_simplex(A_iq, b_iq, A_eq, b_eq, c, tol)

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
 * Tone (2002) Super-SBM uses projected quantities, with DMU k excluded.
 * Ordinary SBM scores below 1-tol are retained in the chosen orientation.
 * Non-oriented uses Charnes-Cooper; oriented models fix the opposite side.
 *
 * Interface matches _dea_ipm_supersbm() exactly.
 */
void _dea_simplex_supersbm(string scalar inputs,
                            string scalar outputs,
                            string scalar eff_var,
                            string scalar touse,
                            string scalar rts,
    | real scalar tol, string scalar orientation)
{
    external real scalar _dea_rec_on
    real matrix X, Y, X_ref, Y_ref
    real colvector obs_indices, x0, y0
    real rowvector c, coef_xbar, coef_ybar
    real matrix A_iq, A_eq
    real colvector b_iq, b_eq
    real scalar n, m, s, k, obs_k, n_ref, j
    real scalar sbm_eff, sbm_score

    struct RSM_Result scalar res

    real colvector ordinary_scores
    real rowvector x_scale, y_scale
    string scalar ort

    if (args()<6) tol=1e-7

    if (args()<7) orientation=""
    ort = strlower(orientation)
    if (ort=="in") ort="input"
    if (ort=="out") ort="output"
    if (ort=="inout") ort=""
    // Classify and retain ordinary SBM scores using the requested orientation.
    _dea_simplex_sbm(inputs, outputs, eff_var, "", "", touse, rts, ort, tol)
    ordinary_scores = st_data(., eff_var, touse)

    // No peer exists even if rounding puts the ordinary score just below 1.
    if (rows(ordinary_scores)<2) {
        st_store(., eff_var, touse, J(rows(ordinary_scores),1,.))
        return
    }

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    x_scale = colmax(X); y_scale = colmax(Y)
    X = X :/ x_scale; Y = Y :/ y_scale
    obs_indices = selectindex(st_data(., touse))

    n = rows(X); m = cols(X); s = cols(Y)

    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        if (_dea_rec_on) _dea_rec_context(obs_k,3,rts)
        x0 = X[k, .]'; y0 = Y[k, .]'
        sbm_eff = .

        sbm_score = ordinary_scores[k]

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
        //        x_tilde >= t * x_0
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

        // Non-oriented: fractional model. Oriented: t=1 and fix the
        // opposite-side projection. Output score is the reciprocal maximum.
        c = J(1, n_ref, 0), coef_xbar, J(1, s, 0), 0
        if (ort=="") {
            A_eq = (J(1, n_ref+m, 0), coef_ybar, 0)
            b_eq = 1
        }
        else {
            A_eq = (J(1, n_ref+m+s, 0), 1)
            b_eq = 1
            if (ort=="input") {
                A_eq = A_eq \ (J(s,n_ref+m,0), I(s), -y0)
                b_eq = b_eq \ J(s,1,0)
            }
            else {
                c = J(1,n_ref+m,0), -coef_ybar, 0
                A_eq = A_eq \ (J(m,n_ref,0), I(m), J(m,s,0), -x0)
                b_eq = b_eq \ J(m,1,0)
            }
        }

        // Inequality 1: x_tilde >= X_ref' * Lambda
        A_iq = (X_ref', -I(m), J(m, s, 0), J(m, 1, 0))
        b_iq = J(m, 1, 0)

        // Inequality 2: y_tilde <= Y_ref' * Lambda
        A_iq = A_iq \ (-Y_ref', J(s, m, 0), I(s), J(s, 1, 0))
        b_iq = b_iq \ J(s, 1, 0)

        // Inequality 3: x_tilde >= t * x_0
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

        res = _dea_nonte_simplex(A_iq, b_iq, A_eq, b_eq, c, tol)

        if (res.status == 0) {
            if (ort!="output") sbm_eff = res.objval
            else if (res.objval<0) sbm_eff = -1/res.objval
        }

        if (_dea_rec_on) _dea_rec_group(obs_k,3,rts,sbm_eff<.,.)
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
                              string scalar p_var,
    | real scalar tol)
{
    real matrix X, Y
    real colvector obs_indices, x0, y0
    real scalar n, m, s, k, obs_k
    real scalar te_k, econ_k, alloc_k
    real colvector te_result

    if (args()<9) tol=1e-7

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
        te_result = _dea_solve_dmu(X, Y, k, orientation, rts, tol, "")
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

            res = _dea_nonte_simplex(A_iq, b_iq, A_eq, b_eq, c_econ, tol)

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

            res = _dea_nonte_simplex(A_iq, b_iq, A_eq, b_eq, c_econ, tol)

            if (res.status == 0 && obs_val > 0) {
                opt_val = -res.objval
                econ_k = obs_val / opt_val
                alloc_k = econ_k / te_k
            }
        }

        st_store(obs_k, st_varindex(eff_var), alloc_k)
    }
}


// Conservative fast acceptance only. Zero means use compensated verification.
// See tests/INVERSE_VALIDATION.txt for the bound and arithmetic assumptions.
real scalar _rsm_inverse_fast(real matrix B, real matrix X)
{
    real scalar m,a,b,radius,err
    real matrix AB,AX,left,right
    m=rows(B)
    if(m==0 | m>4096 | cols(B)!=m | rows(X)!=m | cols(X)!=m) return(0)
    if(hasmissing(B) | hasmissing(X)) return(0)
    AB=abs(B);AX=abs(X);a=max(AB);b=max(AX)
    if(a==0 | b==0 | a>1.6069380442589903e60 | b>1.6069380442589903e60) return(0)
    if(any((AB:>0):&(AB:<6.223015277861142e-61)) | any((AX:>0):&(AX:<6.223015277861142e-61))) return(0)
    radius=64*epsilon(1)*m*m*(a*b+1)
    if(missing(radius) | radius>1e-8/4) return(0)
    left=B*X-I(m);right=X*B-I(m)
    if(hasmissing(left) | hasmissing(right)) return(0)
    err=max((max(abs(left)),max(abs(right))))
    return(err<=1e-8/2)
}

real matrix _rsm_inverse_residual(real matrix B,real matrix X,real matrix D)
{
    real scalar i,j,k,a,b,p,e,t,ah,al,bh,bl,shift
    real rowvector terms
    real matrix R
    if(cols(B)!=rows(X) | rows(B)!=rows(D) | cols(X)!=cols(D)) return(J(rows(D),cols(D),.))
    R=J(rows(D),cols(D),.)
    for(i=1;i<=rows(D);i++) for(j=1;j<=cols(D);j++) {
        terms=J(1,2*cols(B)+1,0);terms[1]=D[i,j]
        for(k=1;k<=cols(B);k++) {
            a=B[i,k];b=X[k,j]
            if(missing(a) | missing(b)) return(J(rows(D),cols(D),.))
            if(a==0 | b==0) continue
            // Reciprocal binary shifts preserve each product while moving
            // extreme operands into the existing splitter's safe range.
            // Exact binary64 literal for 2^256: Mata exponentiation is approximate.
            // Reject if eight shifts cannot do this without missing values.
            for(shift=0;shift<8;shift++) {
                if(abs(a)>=3.4395525670743494e-136 & abs(a)<=2.9073548971824276e+135 & abs(b)>=3.4395525670743494e-136 & abs(b)<=2.9073548971824276e+135) break
                if(abs(a)>2.9073548971824276e+135 | abs(b)<3.4395525670743494e-136) {a=a/1.157920892373162e77;b=b*1.157920892373162e77;}
                else {a=a*1.157920892373162e77;b=b/1.157920892373162e77;}
                if(missing(a) | missing(b) | a==0 | b==0) return(J(rows(D),cols(D),.))
            }
            if(shift==8) return(J(rows(D),cols(D),.))
            p=a*b
            t=134217729*a;ah=t-(t-a);al=a-ah
            t=134217729*b;bh=t-(t-b);bl=b-bh
            e=((ah*bh-p)+ah*bl+al*bh)+al*bl
            terms[2*k]=-p;terms[2*k+1]=-e
        }
        R[i,j]=quadcross(J(cols(terms),1,1),terms')
    }
    return(R)
}


real matrix _rsm_factor_residual(real matrix B,real matrix X,real matrix D)
{
    real scalar i,j,k,a,b,p,e,t,ah,al,bh,bl
    real rowvector terms
    real matrix R
    if(cols(B)!=rows(X) | rows(B)!=rows(D) | cols(X)!=cols(D)) return(J(rows(D),cols(D),.))
    R=J(rows(D),cols(D),.)
    for(i=1;i<=rows(D);i++) for(j=1;j<=cols(D);j++) {
        terms=J(1,2*cols(B)+1,0);terms[1]=D[i,j]
        for(k=1;k<=cols(B);k++) {
            a=B[i,k];b=X[k,j]
            if(missing(a) | missing(b)) return(J(rows(D),cols(D),.))
            if(a==0 | b==0) continue
            if(abs(a)<3.4395525670743494e-136 | abs(a)>2.9073548971824276e+135 | abs(b)<3.4395525670743494e-136 | abs(b)>2.9073548971824276e+135) return(J(rows(D),cols(D),.))
            p=a*b
            t=134217729*a;ah=t-(t-a);al=a-ah
            t=134217729*b;bh=t-(t-b);bl=b-bh
            e=((ah*bh-p)+ah*bl+al*bh)+al*bl
            terms[2*k]=-p;terms[2*k+1]=-e
        }
        R[i,j]=quadcross(J(cols(terms),1,1),terms')
    }
    return(R)
}

real matrix _rsm_factor_solve(real matrix B,real matrix D)
{
    real scalar i,iteration,old_err,err
    real matrix C,X,correction,R,best,comparison_scale
    real colvector rs
    real rowvector cs
    if(rows(B)!=cols(B) | rows(B)!=rows(D) | hasmissing(B) | hasmissing(D)) return(J(rows(B),cols(D),.))
    if(rows(B)==0) return(J(0,cols(D),.))
    rs=rowmax(abs(B));if(any(rs:==0)) return(J(rows(B),cols(D),.))
    C=B:/rs;cs=colmax(abs(C));if(any(cs:==0)) return(J(rows(B),cols(D),.))
    C=C:/cs;X=lusolve(C,D:/rs)
    if(hasmissing(X)) return(X)
    X=X:/cs';best=X;old_err=.
    comparison_scale=abs(D)+abs(B)*(abs(X):+1):+1e-300
    for(iteration=0;iteration<=4;iteration++) {
        R=_rsm_factor_residual(B,X,D)
        if(hasmissing(R)) return(J(rows(B),cols(D),.))
        err=max(abs(R):/comparison_scale)
        if(err<old_err) {best=X;old_err=err;}
        else break
        if(iteration==4 | err==0) break
        correction=lusolve(C,R:/rs)
        if(hasmissing(correction)) break
        X=X+correction:/cs'
    }
    return(best)
}

struct RSM_Result scalar _rsm_fp_two_phase(struct RSM_Prob scalar prob,
                                         real scalar tol,
                                         string scalar trace,| real colvector final_basic)
{
    struct RSM_Prob scalar original, retry
    struct RSM_Result scalar first, second
    real scalar i, scale
    real rowvector check, column_scale

    original=prob
    first=_rsm_fp_once(prob,tol,trace,0,final_basic)
    if (first.status!=3) return(first)

    // Bounded numerical recovery, inspired by Benchmarking's fresh-basis
    // retry and equilibration. Rebuild from the untouched LP, never from
    // a partially pivoted state. Positive row/objective scaling is equivalent;
    // variable units are mapped back and the user's radial target is unchanged.
    retry=original
    for (i=1;i<=retry.m_con;i++) {
        scale=max((abs(retry.A[i,1..retry.n_orig]),abs(retry.b[i])))
        if (scale>0 & scale<.) {
            retry.A[i,.]=retry.A[i,.]/scale
            retry.b[i]=retry.b[i]/scale
        }
    }
    // Artificial columns are temporary Phase-I coordinates. Restore their
    // unit columns so Phase I minimizes violations in the equilibrated rows,
    // rather than retaining weights from the original measurement units.
    for (i=retry.n_orig+1;i<=retry.n_var;i++) {
        if (retry.is_artif[i]) retry.A[.,i]=(retry.A[.,i]:!=0)
    }
    column_scale=J(1,retry.n_orig,1)
    for (i=1;i<=retry.n_orig;i++) {
        scale=max(abs(retry.A[.,i]))
        if (scale>0 & scale<.) {
            column_scale[i]=1/scale
            retry.A[.,i]=retry.A[.,i]/scale
            retry.c[i]=retry.c[i]/scale
        }
    }
    scale=max(abs(retry.c))
    if (scale>0 & scale<.) retry.c=retry.c/scale
    if (hasmissing(retry.A) | hasmissing(retry.b) | hasmissing(retry.c) | hasmissing(column_scale)) return(first)
    if (trace=="trace") printf("{txt}RSM restart: fresh basis, row/column/objective equilibration, Bland\n")
    second=_rsm_fp_once(retry,tol,trace,1,final_basic)
    second.iters=first.iters+second.iters
    if (second.status==0) {
        second.x=second.x:*column_scale'
        check=_dea_lpcheck(second.x',original.A[.,1..original.n_orig],original.b,
            J(0,original.n_orig,0),J(0,1,0),J(1,original.n_orig,0),J(1,original.n_orig,.))
        second.objval=original.c[1..original.n_orig]*second.x
        if (_dea_lpcheck_ok(check,tol) & second.objval<.) {
            second.msg="Optimal after equilibrated restart"
            prob=original
            return(second)
        }
    }
    // A numerical restart is not an infeasibility/unboundedness certificate.
    first.iters=second.iters
    return(first)
}

struct RSM_Result scalar _rsm_fp_once(struct RSM_Prob scalar prob,
                                         real scalar tol,
                                         string scalar trace,
                                         real scalar restart,| real colvector final_basic)
{
    real scalar i, max_iter, p1_iter
    real colvector xb
    real rowvector c_orig
    real scalar art_sum
    real matrix B0, original_A
    real colvector original_b, feasible_x
    real rowvector check

    struct RSM_Basis scalar basis
    struct RSM_Result scalar result

    final_basic=J(0,1,.)
    max_iter = max((10000, 20 * prob.m_con))

    // Save Phase 2 objective
    c_orig = prob.c
    original_A = prob.A[.,1..prob.n_orig]
    original_b = prob.b

    // =================================================================
    // PHASE 1: Minimize sum of artificials
    // =================================================================
    prob.c = prob.is_artif

    // Initialize basis
    basis.basic = prob.init_basic
    basis.n_upd = 0

    // Initial B^{-1} (should be identity by construction)
    B0 = prob.A[., basis.basic]
    if (all(B0 :== I(prob.m_con))) {
        basis.Binv = I(prob.m_con)
    }
    else {
        basis.Binv = _rsm_inv(B0)
        if (hasmissing(basis.Binv)) {
            result.status = 3
            result.msg = "Invalid initial basis inverse"
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
    result = _rsm_fp_solve(prob, basis, tol, max_iter, 1, trace, restart)
    p1_iter = result.iters

    if (result.status == 2) {
        result.status = 3
        result.msg = "Phase1 unbounded; no feasibility certificate"
        prob.c = c_orig
        return(result)
    }

    if (result.status != 0) {
        prob.c = c_orig
        return(result)
    }

    // Check artificial sum
    xb = _rsm_factor_solve(prob.A[,basis.basic],prob.b)
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

    // A tiny artificial objective alone is not a feasible original LP.
    feasible_x=J(prob.n_orig,1,0)
    for (i=1;i<=prob.m_con;i++) {
        if (basis.basic[i]<=prob.n_orig) feasible_x[basis.basic[i]]=max((0,xb[i]))
    }
    check=_dea_lpcheck(feasible_x',original_A,original_b,
        J(0,prob.n_orig,0),J(0,1,0),J(1,prob.n_orig,0),J(1,prob.n_orig,.))
    if (!_dea_lpcheck_ok(check,tol)) {
        result.status=3; result.msg="Residual check failed"
        result.x=feasible_x; result.objval=.
        prob.c=c_orig
        return(result)
    }

    // =================================================================
    // TRANSITION: Purge artificials, drop columns, restore objective
    // =================================================================
    if (!_rsm_purge_artif(prob, basis, tol)) {
        result.status=3; result.msg="Artificial-basis pivot failed"
        result.objval=.; result.x=feasible_x
        return(result)
    }
    if (!_rsm_drop_artif(prob, basis)) {
        result.status=3; result.msg="Artificial-basis transition failed"
        result.objval=.; result.x=feasible_x
        return(result)
    }
    prob.c = c_orig[1..prob.n_var]

    // Refactorize for clean Phase 2 start
    if (!_rsm_refac(basis, prob)) {
        result.status=3; result.msg="Invalid Phase2 basis inverse"
        result.objval=.; result.x=feasible_x
        return(result)
    }

    // =================================================================
    // PHASE 2: Optimize original objective
    // =================================================================
    if (trace == "trace") {
        printf("{txt}Phase2: %g con, %g var\n", prob.m_con, prob.n_var)
    }

    result = _rsm_fp_solve(prob, basis, tol, max_iter - p1_iter, 2, trace, restart)
    result.iters = result.iters + p1_iter

    // Extract full solution
    if (result.status == 0) {
        result.x = J(prob.n_var, 1, 0)
        xb = _rsm_factor_solve(prob.A[,basis.basic],prob.b)
        for (i = 1; i <= prob.m_con; i++) {
            result.x[basis.basic[i]] = max((0, xb[i]))
        }
        result.objval = prob.c * result.x
        check=_dea_lpcheck(result.x',original_A,original_b,
            J(0,prob.n_var,0),J(0,1,0),J(1,prob.n_var,0),J(1,prob.n_var,.))
        if (!_dea_lpcheck_ok(check,tol) | result.objval>=.) {
            result.status=3; result.msg="Residual check failed"; result.objval=.
        }
    }

    if(result.status==0) final_basic=basis.basic
    return(result)
}

struct RSM_Result scalar _rsm_fp_solve(struct RSM_Prob scalar prob,
    struct RSM_Basis scalar basis,real scalar tol,real scalar max_iter,
    real scalar phase,string scalar trace,real scalar restart)
{
    struct RSM_Result scalar result
    real scalar iter,entering,pivot_row,use_bland,stall_count,consec_degen,prev_obj,obj_val,i
    real colvector xb,d
    result.status=3;result.objval=.;result.x=J(prob.n_var,1,.);result.iters=0
    result.msg="Maximum iterations";use_bland=restart;stall_count=consec_degen=0;prev_obj=.
    for(iter=1;iter<=max_iter;iter++) {
        result.iters=iter
        xb=_rsm_factor_solve(prob.A[,basis.basic],prob.b)
        if(hasmissing(xb)) {result.msg="Factorized basic solution unavailable";return(result);}
        if(min(xb)<-min((_dea_tol_feas(),tol))) {result.msg="Factorized primal basis lost feasibility";return(result);}
        obj_val=prob.c[basis.basic']*xb
        if(prev_obj<. & abs(obj_val-prev_obj)<_dea_tol_zero()) {
            stall_count++
            if(stall_count>prob.m_con) use_bland=1
        }
        else stall_count=0
        prev_obj=obj_val
        entering=_rsm_fp_enter(prob,basis,use_bland,tol)
        if(entering<0) {result.msg="Factorized price unavailable";return(result);}
        if(entering==0) {
            // Original Phase-1 transition still requires a checked inverse.
            // The failure-only retry does not change dependent-row removal.
            if(phase==1) {
                basis.Binv=_rsm_inv(prob.A[,basis.basic]);basis.n_upd=0
                if(hasmissing(basis.Binv)) {result.msg="Phase1 transition inverse unavailable";return(result);}
            }
            result.status=0;result.objval=obj_val;result.msg=sprintf("Factorized optimal P%g i%g",phase,iter)
            result.x=J(prob.n_var,1,0)
            for(i=1;i<=prob.m_con;i++) result.x[basis.basic[i]]=max((0,xb[i]))
            return(result)
        }
        d=_rsm_factor_solve(prob.A[,basis.basic],prob.A[,entering])
        if(hasmissing(d)) {result.msg="Factorized pivot unavailable";return(result);}
        pivot_row=_rsm_leave(xb,d,basis,use_bland,tol)
        if(pivot_row==0) {
            if(rows(d)>0 & max(d)>0) {result.msg="No stable factorized pivot";return(result);}
            result.status=2;result.objval=obj_val;result.msg="Factorized unbounded";return(result)
        }
        if(abs(xb[pivot_row])<_dea_tol_feas()) {
            consec_degen++
            if(consec_degen>2*prob.m_con) use_bland=1
        }
        else consec_degen=0
        if(trace=="trace") printf("FS P%g i%g enter%g leave%g step%e pivot%e\\n",phase,iter,entering,basis.basic[pivot_row],max((0,xb[pivot_row]))/d[pivot_row],d[pivot_row])
        basis.basic[pivot_row]=entering
    }
    return(result)
}

real scalar _rsm_fp_enter(struct RSM_Prob scalar prob,
                        struct RSM_Basis scalar basis,
                        real scalar use_bland,
                        real scalar tol)
{
    real rowvector y, c_B, rc_all
    real scalar i, j, n_nb, rc_j, best_rc, opt_tol
    real colvector nb, is_basic

    // Price small objectives in their own units. A large coefficient must
    // not loosen the absolute cutoff for other, small-cost columns.
    opt_tol = min((_dea_tol_opt(),tol))*min((1,max(abs(prob.c))))

    // Simplex multiplier: y = c_B * B^{-1}  [VECTORIZED]
    c_B = J(1, prob.m_con, 0)
    for (i = 1; i <= prob.m_con; i++) {
        c_B[i] = prob.c[basis.basic[i]]
    }
    y = _rsm_factor_solve(prob.A[,basis.basic]',c_B')'
    if(hasmissing(y)) return(-1)

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
            rc_j = _rsm_factor_residual(y,prob.A[.,j],prob.c[j])
            if(missing(rc_j)) return(-1)
            if (rc_j < -opt_tol) return(j)
        }
        return(0)
    }
    else {
        // DANTZIG'S RULE: most negative rc [VECTORIZED]
        rc_all = _rsm_factor_residual(y,prob.A[.,nb],prob.c[nb'])
        if(hasmissing(rc_all)) return(-1)    // 1 x n_nb
        best_rc = min(rc_all)
        if (best_rc >= -opt_tol) return(0)

        for (i = 1; i <= n_nb; i++) {
            if (rc_all[i] == best_rc) {
                return(nb[i])
            }
        }
        return(0)
    }
}


struct RSM_Result scalar _rsm_retry_two_phase(struct RSM_Prob scalar prob,
    real scalar tol,string scalar trace,real colvector certbasis)
{
    real scalar cycle_seen
    return(_rsm_precision_two_phase(prob,tol,trace,certbasis,0,cycle_seen))
}


// Propose a complete real-column basis, preserving every original row.
// QR selects columns only; the unchanged original equations approve the trial.
// This is a bounded two-proposal reconstruction, not artificial-variable pivoting.

// Thin pivoted Householder selection. Store m-by-n work only.
// Recompute scaled trailing norms, avoiding cancellation-prone norm downdates.
// This proposes column indices; unchanged original-basis checks approve them.
real colvector _rsm_retry_columns(real matrix C,real rowvector fixed)
{
    real scalar m,n,k,j,nfixed,best,scale,beta,tau,alpha
    real matrix W,T
    real colvector v,order
    real rowvector p,norms,scales,safe
    m=rows(C);n=cols(C)
    if(m<1 | n<m | cols(fixed)!=n | hasmissing(C) | hasmissing(fixed)) return(J(0,1,.))
    p=1..n;nfixed=0
    // Match fixed-column front-loading by swaps, as in pivoted QR.
    for(j=1;j<=n;j++) if(fixed[j]) {
        nfixed++
        best=p[j];p[j]=p[nfixed];p[nfixed]=best
    }
    if(nfixed>m) return(J(0,1,.))
    W=C[,p]
    for(k=1;k<=m;k++) {
        if(k>nfixed) {
            T=W[k..m,k..n]
            scales=colmax(abs(T));safe=scales
            for(j=1;j<=cols(safe);j++) if(safe[j]==0) safe[j]=1
            norms=scales:*sqrt(quadcross(J(rows(T),1,1),(T:/safe):^2))
            if(hasmissing(norms) | max(norms)==0) return(J(0,1,.))
            order=select((1::cols(T)),norms':==max(norms));best=k+order[1]-1
            if(best!=k) {
                v=W[,k];W[,k]=W[,best];W[,best]=v
                j=p[k];p[k]=p[best];p[best]=j
            }
        }
        v=W[k..m,k];scale=max(abs(v))
        if(scale==0 | scale>=.) return(J(0,1,.))
        v=v/scale;alpha=v[1]
        beta=(alpha>=0 ? -1 : 1)*sqrt(quadcross(v,v))
        if(beta==0 | beta>=.) return(J(0,1,.))
        tau=(beta-alpha)/beta;v[1]=alpha-beta;v=v/v[1]
        T=W[k..m,k..n]
        W[k..m,k..n]=T-v*(tau*quadcross(v,T))
        W[k,k]=beta*scale
        if(k<m) W[k+1..m,k]=J(m-k,1,0)
        if(hasmissing(W)) return(J(0,1,.))
    }
    return(p[1..m]')
}

real scalar _rsm_retry_realbasis(struct RSM_Prob scalar prob,
    struct RSM_Basis scalar basis,real scalar tol)
{
    real scalar m,n,attempt,i,feas
    real matrix A,C,trial_B,H,tau,R,selection
    real colvector rs,trial_basic,trial_x,full_x
    real rowvector cs,p,check
    m=prob.m_con;n=prob.n_orig;feas=min((_dea_tol_feas(),tol))
    if(m<1 | n<m | rows(basis.basic)!=m | tol<=0 | tol>=.) return(0)
    A=prob.A[,1..n]
    if(hasmissing(A) | hasmissing(prob.b)) return(0)
    rs=rowmax(abs(A));if(any(rs:==0)) return(0)
    C=A:/rs;cs=colmax(abs(C))
    for(i=1;i<=n;i++) if(cs[i]==0) cs[i]=1
    C=C:/cs;if(hasmissing(C)) return(0)
    for(attempt=1;attempt<=2;attempt++) {
        p=J(1,n,0)
        if(attempt==1) for(i=1;i<=m;i++) if(basis.basic[i]<=n) p[basis.basic[i]]=1
        // Retain the previously validated selector for bounded small cases.
        // Wide problems use thin Householder work: no n-by-n allocation.
        if(n<=128) {
            selection=(C\J(n-m,n,0))
            hqrdp(selection,H,tau,R,p)
            if(hasmissing(R) | hasmissing(p)) continue
            trial_basic=p[1..m]'
        }
        else trial_basic=_rsm_retry_columns(C,p)
        if(rows(trial_basic)!=m) continue
        trial_B=A[,trial_basic]
        trial_x=_rsm_factor_solve(trial_B,prob.b)
        if(hasmissing(trial_x)) continue
        if(min(trial_x)<-feas) continue
        full_x=J(n,1,0)
        for(i=1;i<=m;i++) full_x[trial_basic[i]]=max((0,trial_x[i]))
        check=_dea_lpcheck(full_x',A,prob.b,
            J(0,n,0),J(0,1,0),J(1,n,0),J(1,n,.))
        if(!_dea_lpcheck_ok(check,tol)) continue
        basis.basic=trial_basic;basis.n_upd=0
        // No inverse is used after this artificial-free transition.
        basis.Binv=J(m,m,.)
        return(1)
    }
    return(0)
}

struct RSM_Result scalar _rsm_retry_once(struct RSM_Prob scalar prob,
    real scalar tol,string scalar trace,real scalar restart,real colvector certbasis)
{
    return(_rsm_precision_once(prob,tol,trace,restart,certbasis,0))
}


// Keep small positive directions in the failure-only feasible step bound.
// Only stable pivots may enter the candidate set; bounded arithmetic rejects
// uncertain comparisons. This screens stored xb/d, not the exact basis solve.
real scalar _rsm_retry_leave(real colvector xb,real colvector d,
    struct RSM_Basis scalar basis,real scalar use_bland,real scalar tol,
    real colvector excluded)
{
    real scalar m,i,j,theta,feas,pvt,ok,upper,pivot,minratio,best,minindex
    real colvector ratios
    real rowvector bound
    m=rows(d);feas=min((_dea_tol_feas(),tol));pvt=_dea_tol_pivot()
    if(m<1 | rows(xb)!=m | rows(basis.basic)!=m | rows(excluded)!=m |
        hasmissing(xb) | hasmissing(d) | hasmissing(basis.basic) |
        hasmissing(excluded) | tol<=0 | tol>=.) return(0)
    if(min(xb)<-feas) return(0)
    ratios=J(m,1,.)
    for(i=1;i<=m;i++) {
        if(excluded[i] | d[i]<=pvt) continue
        theta=max((0,xb[i]))/d[i]
        if(theta>=.) continue
        ok=1
        for(j=1;j<=m;j++) {
            if(d[j]<=0) continue
            // Upper bound for d[j]*theta-xb[j]-feas, without division gates.
            bound=_rsm_q_dot((-d[j],1,1),(theta\xb[j]\feas),0)
            if(hasmissing(bound)) {ok=0;break;}
            upper=_rsm_q_add(bound[1],bound[2],1)
            if(upper>=. | upper>0) {ok=0;break;}
        }
        if(ok) ratios[i]=theta
    }
    minratio=min(ratios);if(minratio>=.) return(0)
    pivot=0;best=0;minindex=.
    for(i=1;i<=m;i++) {
        if(ratios[i]>=.) continue
        if(use_bland) {
            if(ratios[i]==minratio & basis.basic[i]<minindex) {pivot=i;minindex=basis.basic[i];}
        }
        else if(d[i]>best) {pivot=i;best=d[i];}
    }
    return(pivot)
}
// Screen an actual replacement basis before committing the exchange.
// A rejected candidate remains a step-limiting row for all later candidates.
real scalar _rsm_retry_exchange(struct RSM_Prob scalar prob,
    struct RSM_Basis scalar basis,real colvector xb,real colvector d,
    real scalar entering,real scalar use_bland,real scalar tol)
{
    real scalar m,pivot,attempt,feas
    real colvector excluded,trial_basic,trial_x
    real matrix trial_B
    real rowvector check
    m=prob.m_con;feas=min((_dea_tol_feas(),tol));excluded=J(m,1,0)
    if(m<1 | rows(xb)!=m | rows(d)!=m | rows(basis.basic)!=m |
        entering<1 | entering>prob.n_var | entering!=floor(entering) |
        any(basis.basic:==entering)) return(0)
    for(attempt=1;attempt<=m;attempt++) {
        pivot=_rsm_retry_leave(xb,d,basis,use_bland,tol,excluded)
        if(pivot==0) return(0)
        trial_basic=basis.basic;trial_basic[pivot]=entering
        trial_B=prob.A[,trial_basic]
        trial_x=_rsm_factor_solve(trial_B,prob.b)
        if(!hasmissing(trial_x)) {
            if(min(trial_x)>=-feas) {
                check=_dea_lpcheck(trial_x',trial_B,prob.b,
                    J(0,m,0),J(0,1,0),J(1,m,0),J(1,m,.))
                if(_dea_lpcheck_ok(check,tol)) return(pivot)
            }
        }
        excluded[pivot]=1
    }
    return(0)
}

struct RSM_Result scalar _rsm_retry_solve(struct RSM_Prob scalar prob,
    struct RSM_Basis scalar basis,real scalar tol,real scalar max_iter,
    real scalar phase,string scalar trace,real scalar restart)
{
    return(_rsm_precision_solve(prob,basis,tol,max_iter,phase,trace,restart,0))
}


// Original-basis refinement with a separately stored low part.
real matrix _rsm_retry_refine(real matrix B,real colvector D)
{
    real scalar m,iteration,err,old_err
    real colvector hi,lo,R,correction,S,V,E,T,H,L,comparison
    real matrix best
    m=rows(B)
    if(m==0 | cols(B)!=m | rows(D)!=m | cols(D)!=1 | hasmissing(B) | hasmissing(D)) return(J(m,2,.))
    hi=_rsm_factor_solve(B,D);lo=J(m,1,0)
    if(hasmissing(hi)) return(J(m,2,.))
    comparison=abs(D)+abs(B)*(abs(hi):+1):+1e-300
    if(hasmissing(comparison) | any(comparison:<=0)) return(J(m,2,.))
    best=(hi,lo);old_err=.
    for(iteration=0;iteration<=8;iteration++) {
        R=_rsm_inverse_residual((B,B),(hi\lo),D)
        if(hasmissing(R)) return(J(m,2,.))
        err=max(abs(R):/comparison)
        if(err<old_err) {best=(hi,lo);old_err=err;}
        else break
        if(iteration==8 | err==0) break
        correction=_rsm_factor_solve(B,R)
        if(hasmissing(correction)) break
        // TwoSum twice; do not round the accumulated low component into hi.
        S=hi+correction;V=S-hi
        E=(hi-(S-V))+(correction-V)
        T=lo+E;H=S+T;V=H-S
        L=(S-(H-V))+(T-V)
        if(hasmissing(H) | hasmissing(L)) break
        hi=H;lo=L
    }
    return(best)
}
// Bounded arithmetic on a deliberately finite exponent domain.
// See tests/PRECISION_RECOVERY.txt for assumptions, limits and exact checks.
real rowvector _rsm_q_sum(real scalar a,real scalar b)
{
    real scalar s,v,e
    s=a+b;v=s-a;e=(a-(s-v))+(b-v)
    return((s,e))
}
real scalar _rsm_q_add(real scalar a,real scalar b,real scalar upper)
{
    real rowvector p
    real scalar s
    if(missing(a) | missing(b) | max((abs(a),abs(b)))>8.452712498170644e270) return(.)
    if((a!=0 & abs(a)<1.0020841800044864e-292) | (b!=0 & abs(b)<1.0020841800044864e-292)) return(.)
    p=_rsm_q_sum(a,b);s=p[1]
    if(missing(s) | (s!=0 & abs(s)<1.0020841800044864e-292)) return(.)
    if(upper & p[2]>0) s=s+abs(s)*2.2204460492503131e-16
    if(!upper & p[2]<0) s=s-abs(s)*2.2204460492503131e-16
    return(s)
}
real rowvector _rsm_q_product(real scalar a,real scalar b)
{
    real scalar p,t,ah,al,bh,bl,e
    if(missing(a) | missing(b)) return((.,.))
    if(a==0 | b==0) return((0,0))
    // Exact powers 2^-400 and 2^400. Products and splitting stay normal.
    if(min((abs(a),abs(b)))<3.8725919148493183e-121 | max((abs(a),abs(b)))>2.5822498780869086e120) return((.,.))
    p=a*b
    t=134217729*a;ah=t-(t-a);al=a-ah
    t=134217729*b;bh=t-(t-b);bl=b-bh
    e=al*bl-(((p-ah*bh)-al*bh)-ah*bl)
    return((p,e))
}
real scalar _rsm_q_mulup(real scalar a,real scalar b)
{
    real rowvector p
    if(a<0 | b<0) return(.)
    p=_rsm_q_product(a,b)
    if(hasmissing(p)) return(.)
    if(p[2]>0) return(p[1]+abs(p[1])*2.2204460492503131e-16)
    return(p[1])
}
// Enclose D-a*b as (center,radius), retaining actual discarded EFT tails.
real rowvector _rsm_q_dot(real rowvector a,real colvector b,real scalar D)
{
    real scalar i,hi,lo,radius
    real rowvector p,s,t,u,v
    if(cols(a)!=rows(b) | cols(a)>4096 | missing(D) | hasmissing(a) | hasmissing(b)) return((.,.))
    if(D!=0 & (abs(D)<3.8725919148493183e-121 | abs(D)>2.5822498780869086e120)) return((.,.))
    hi=D;lo=radius=0
    for(i=1;i<=cols(a);i++) {
        p=_rsm_q_product(a[i],b[i]);if(hasmissing(p)) return((.,.))
        s=_rsm_q_sum(hi,-p[1]);t=_rsm_q_sum(lo,-p[2])
        u=_rsm_q_sum(s[2],t[1]);v=_rsm_q_sum(s[1],u[1])
        hi=v[1];lo=v[2]
        radius=_rsm_q_add(radius,_rsm_q_add(abs(t[2]),abs(u[2]),1),1)
        if(missing(radius)) return((.,.))
    }
    v=_rsm_q_sum(hi,lo);radius=_rsm_q_add(radius,abs(v[2]),1)
    return((v[1],radius))
}
// Upper enclosure of numerator-cutoff*denominator. No rounded division gate.
real scalar _rsm_q_compare(real scalar numerator,real scalar cutoff,real scalar denominator)
{
    real rowvector p
    p=_rsm_q_dot(cutoff,denominator,numerator)
    if(hasmissing(p)) return(.)
    return(_rsm_q_add(p[1],p[2],1))
}
// Preserve 1+abs(objective) as an expression, avoiding denominator rounding.
real scalar _rsm_q_scaled(real scalar numerator,real scalar cutoff,real rowvector O)
{
    real rowvector p
    if(abs(O[1])<O[2]) return(_rsm_q_compare(numerator,cutoff,1))
    p=_rsm_q_dot((cutoff,cutoff,-cutoff),(1\abs(O[1])\O[2]),numerator)
    if(hasmissing(p)) return(.)
    return(_rsm_q_add(p[1],p[2],1))
}
real rowvector _rsm_retry_quality(struct RSM_Prob scalar original,
    real colvector x,real rowvector yh,real rowvector yl,real scalar tol)
{
    real matrix A
    real rowvector c,O,Y,G,P,R,cmp
    real scalar i,j,m,n,denlo,gapu,pu,cu,du,ayu,rpu,rdu,lower,cutoff,price,cs,ok
    m=original.m_con;n=original.n_orig
    if(m<1 | n<1 | n+2*m>4096 | rows(x)!=n | cols(yh)!=m | cols(yl)!=m | missing(tol) | tol<=0 | hasmissing(x) | hasmissing(yh) | hasmissing(yl)) return(J(1,16,.))
    A=original.A[,1..n];c=original.c[1..n]
    if(hasmissing(A) | hasmissing(c) | hasmissing(original.b)) return(J(1,16,.))
    O=_rsm_q_dot(-c,x,0);Y=_rsm_q_dot((-yh,-yl),(original.b\original.b),0)
    G=_rsm_q_dot((c,-yh,-yl),(x\original.b\original.b),0)
    if(hasmissing((O,Y,G))) return(J(1,16,.))
    lower=_rsm_q_add(abs(O[1]),-O[2],0)
    if(missing(lower)) return(J(1,16,.))
    denlo=_rsm_q_add(1,max((0,lower)),0)
    gapu=_rsm_q_add(abs(G[1]),G[2],1)
    pu=cu=du=0
    for(i=1;i<=m;i++) {
        P=_rsm_q_dot(A[i,.],x,original.b[i]);if(hasmissing(P)) return(J(1,16,.))
        ayu=_rsm_q_add(abs(yh[i]),abs(yl[i]),1)
        rpu=_rsm_q_add(abs(P[1]),P[2],1)
        pu=_rsm_q_add(pu,_rsm_q_mulup(ayu,rpu),1)
    }
    for(j=1;j<=n;j++) {
        R=_rsm_q_dot((yh,yl),(A[,j]\A[,j]),c[j]);if(hasmissing(R)) return(J(1,16,.))
        rdu=_rsm_q_add(abs(R[1]),R[2],1)
        cu=_rsm_q_add(cu,_rsm_q_mulup(abs(x[j]),rdu),1)
        lower=_rsm_q_add(R[1],-R[2],0)
        if(missing(lower)) return(J(1,16,.))
        du=max((du,-lower))
    }
    cutoff=min((1e-8,tol));price=min((_dea_tol_opt(),tol));cs=min((1,max(abs(c))))
    if(hasmissing((denlo,gapu,pu,cu,du))) return(J(1,16,.))
    cmp=(_rsm_q_scaled(gapu,cutoff,O),_rsm_q_scaled(pu,cutoff,O),_rsm_q_scaled(cu,cutoff,O),_rsm_q_compare(du,price,cs))
    ok=(!hasmissing(cmp) & max(cmp)<=0)
    return((ok,O[1],Y[1],gapu/denlo,pu/denlo,cu/denlo,du,price*cs,denlo,gapu,pu,cu,cmp))
}

// Optional common LP recovery diagnostics. Recording is off by default.
real scalar _rsm_recovery_reason(string scalar msg)
{
    if(strpos(msg,"initial basis")) return(1)
    if(strpos(msg,"basic solution")) return(2)
    if(strpos(msg,"price unavailable")) return(3)
    if(strpos(msg,"pivot unavailable") | strpos(msg,"stable factorized pivot")) return(4)
    if(strpos(msg,"Phase1 transition inverse")) return(5)
    if(strpos(msg,"Artificial-basis pivot")) return(6)
    if(strpos(msg,"Artificial-basis transition")) return(7)
    if(strpos(msg,"Phase2 basis inverse")) return(8)
    if(strpos(msg,"Residual check")) return(9)
    if(strpos(msg,"Maximum iterations")) return(10)
    if(strpos(msg,"lost feasibility")) return(11)
    if(strpos(msg,"refactorization")) return(12)
    return(0)
}
void _rsm_audit_reset(real scalar enabled)
{
    external real scalar _rsm_recording
    external real matrix _rsm_records
    external real rowvector _rsm_counts
    _rsm_recording=enabled
    _rsm_records=J(0,12,.)
    _rsm_counts=J(1,8,0)
}
void _rsm_audit_emit(real rowvector row,real scalar m,real scalar n)
{
    external real scalar _dea_rec_on
    external real scalar _rsm_recording
    external real matrix _rsm_records
    external real rowvector _rsm_counts
    if (_dea_rec_on) _dea_rec_emit(row,m,n)
    if(_rsm_recording!=1) return
    _rsm_counts[1]=_rsm_counts[1]+1
    _rsm_counts[2]=_rsm_counts[2]+row[1]
    _rsm_counts[3]=_rsm_counts[3]+row[2]
    _rsm_counts[4]=_rsm_counts[4]+row[5]
    _rsm_counts[5]=_rsm_counts[5]+row[6]
    _rsm_counts[6]=_rsm_counts[6]+row[7]
    if(rows(_rsm_records)<10000) _rsm_records=_rsm_records\(row,m,n)
    else _rsm_counts[7]=_rsm_counts[7]+1
    _rsm_counts[8]=rows(_rsm_records)
}
real matrix _rsm_audit_read()
{
    external real matrix _rsm_records
    return(_rsm_records)
}
real rowvector _rsm_audit_counts()
{
    external real rowvector _rsm_counts
    return(_rsm_counts)
}
void _rsm_audit_enable(real scalar enabled)
{
    external real scalar _rsm_recording
    _rsm_recording=enabled
}

// Certified recovery still requires the final original-unit primal/dual screen.

// Bounds use the exact real arithmetic interpretation
// of stored binary64 inputs. If validation fails, return unavailable.
real scalar _rsm_cert_divup(real scalar num,real scalar den)
{
 real scalar q,i,upper
 real rowvector R
 if(missing(num) | missing(den) | num<0 | den<=0) return(.)
 if(num==0) return(0)
 q=num/den
 for(i=1;i<=3;i++) {
  if(missing(q) | q<=0) return(.)
  R=_rsm_q_dot(den,q,num)
  if(hasmissing(R)) return(.)
  upper=_rsm_q_add(R[1],R[2],1)
  if(upper<. & upper<=0) return(q)
  q=_rsm_q_add(q,abs(q)*2.2204460492503131e-16,1)
 }
 return(.)
}
// M*y=D. For approximate inverse Z and E=I-Z*M, rho>=||E||inf<1
// certifies ||y-(yh+yl)||inf <= || |Z| |D-M*(yh+yl)| ||inf/(1-rho).
real scalar _rsm_cert_dual_error(real matrix M,real colvector D,
 real matrix Y,real rowvector detail)
{
 real scalar m,i,j,rho,rowsum,den,num,delta,upper
 real matrix Z
 real colvector residual
 real rowvector Q
 detail=J(1,4,.)
 m=rows(M)
 if(m<1 | cols(M)!=m | rows(D)!=m | cols(D)!=1 |
  rows(Y)!=m | cols(Y)!=2 | 2*m>4096 | hasmissing(M) |
  hasmissing(D) | hasmissing(Y)) return(.)
 Z=_rsm_factor_solve(M,I(m))
 if(hasmissing(Z)) return(.)
 rho=0
 for(i=1;i<=m;i++) {
  rowsum=0
  for(j=1;j<=m;j++) {
   Q=_rsm_q_dot(Z[i,.],M[,j],i==j)
   if(hasmissing(Q)) return(.)
   upper=_rsm_q_add(abs(Q[1]),Q[2],1)
   rowsum=_rsm_q_add(rowsum,upper,1)
   if(missing(rowsum)) return(.)
  }
  rho=max((rho,rowsum))
 }
 detail[1]=rho
 if(rho>=1) return(.)
 den=_rsm_q_add(1,-rho,0)
 if(missing(den) | den<=0) return(.)
 residual=J(m,1,.)
 for(i=1;i<=m;i++) {
  Q=_rsm_q_dot((M[i,.],M[i,.]),(Y[,1]\Y[,2]),D[i])
  if(hasmissing(Q)) return(.)
  residual[i]=_rsm_q_add(abs(Q[1]),Q[2],1)
 }
 if(hasmissing(residual)) return(.)
 num=0
 for(i=1;i<=m;i++) {
  Q=_rsm_q_dot(-abs(Z[i,.]),residual,0)
  if(hasmissing(Q)) return(.)
  upper=_rsm_q_add(Q[1],Q[2],1)
  if(missing(upper)) return(.)
  num=max((num,upper))
 }
 delta=_rsm_cert_divup(num,den)
 detail=(rho,delta,den,num)
 return(delta)
}
real matrix _rsm_cert_price_bounds(struct RSM_Prob scalar prob,
 struct RSM_Basis scalar basis,real matrix Y,real rowvector detail)
{
 real scalar m,n,j,k,delta,error,radius,lower,upper
 real matrix out
 real rowvector Q,R
 m=prob.m_con;n=prob.n_var;out=J(n,4,.)
 Y=J(m,2,.);detail=J(1,4,.)
 if(m<1 | n<m | rows(basis.basic)!=m | hasmissing(basis.basic) |
  hasmissing(prob.A) | hasmissing(prob.c) | n+2*m>4096) return(out)
 for(k=1;k<=m;k++) {
  if(basis.basic[k]<1 | basis.basic[k]>n | basis.basic[k]!=floor(basis.basic[k]) |
    sum(basis.basic:==basis.basic[k])!=1) return(out)
 }
 Y=_rsm_retry_refine(prob.A[,basis.basic]',prob.c[basis.basic']')
 if(hasmissing(Y)) return(out)
 delta=_rsm_cert_dual_error(prob.A[,basis.basic]',prob.c[basis.basic']',Y,detail)
 if(missing(delta)) return(out)
 for(j=1;j<=n;j++) {
  Q=_rsm_q_dot((Y[,1]',Y[,2]'),(prob.A[,j]\prob.A[,j]),prob.c[j])
  R=_rsm_q_dot(J(1,m,-delta),abs(prob.A[,j]),0)
  if(hasmissing(Q) | hasmissing(R)) return(J(n,4,.))
  error=_rsm_q_add(R[1],R[2],1)
  radius=_rsm_q_add(Q[2],error,1)
  lower=_rsm_q_add(Q[1],-radius,0);upper=_rsm_q_add(Q[1],radius,1)
  if(hasmissing((radius,lower,upper))) return(J(n,4,.))
  out[j,.]=(Q[1],radius,lower,upper)
 }
 return(out)
}
real scalar _rsm_cert_enter(struct RSM_Prob scalar prob,
 struct RSM_Basis scalar basis,real scalar bland,real scalar tol)
{
 real matrix bounds,Y
 real rowvector detail
 real colvector isbasic
 real scalar j,cutoff,best,enter,uncertain
 if(missing(tol) | tol<=0) return(-1)
 bounds=_rsm_cert_price_bounds(prob,basis,Y,detail)
 if(hasmissing(bounds)) return(-1)
 cutoff=min((_dea_tol_opt(),tol))*min((1,max(abs(prob.c))))
 isbasic=J(prob.n_var,1,0);isbasic[basis.basic]=J(rows(basis.basic),1,1)
 best=0;enter=0;uncertain=0
 for(j=1;j<=prob.n_var;j++) {
  if(isbasic[j]) continue
  if(bounds[j,4]<-cutoff) {
   if(bland) return(j)
   if(enter==0 | bounds[j,1]<best) {enter=j;best=bounds[j,1];}
  }
  else if(bounds[j,3]<-cutoff) {
   if(bland) return(-1)
   uncertain=1
  }
 }
 if(enter==0 & uncertain) return(-1)
 return(enter)
}


struct RSM_Result scalar _rsm_precision_two_phase(struct RSM_Prob scalar prob,
                                         real scalar tol,
                                         string scalar trace, real colvector certbasis, real scalar certified, real scalar cycle_seen,
                                         | real scalar allow_rebuild, real scalar pivot_seen)
{
    struct RSM_Prob scalar original, retry
    struct RSM_Result scalar first, second
    real colvector second_basis
    real scalar i, scale, second_pivot
    real rowvector check, column_scale

    if(args()<7) allow_rebuild=0
    pivot_seen=0
    cycle_seen=0
    original=prob
    first=_rsm_precision_once(prob,tol,trace,0,certbasis,certified,allow_rebuild,pivot_seen)
    cycle_seen=(first.msg=="Repeated Bland pivot state")
    if (first.status!=3) return(first)

    // Bounded numerical recovery, inspired by Benchmarking's fresh-basis
    // retry and equilibration. Rebuild from the untouched LP, never from
    // a partially pivoted state. Positive row/objective scaling is equivalent;
    // variable units are mapped back and the user's radial target is unchanged.
    retry=original
    for (i=1;i<=retry.m_con;i++) {
        scale=max((abs(retry.A[i,1..retry.n_orig]),abs(retry.b[i])))
        if (scale>0 & scale<.) {
            retry.A[i,.]=retry.A[i,.]/scale
            retry.b[i]=retry.b[i]/scale
        }
    }
    // Artificial columns are temporary Phase-I coordinates. Restore their
    // unit columns so Phase I minimizes violations in the equilibrated rows,
    // rather than retaining weights from the original measurement units.
    for (i=retry.n_orig+1;i<=retry.n_var;i++) {
        if (retry.is_artif[i]) retry.A[.,i]=(retry.A[.,i]:!=0)
    }
    column_scale=J(1,retry.n_orig,1)
    for (i=1;i<=retry.n_orig;i++) {
        scale=max(abs(retry.A[.,i]))
        if (scale>0 & scale<.) {
            column_scale[i]=1/scale
            retry.A[.,i]=retry.A[.,i]/scale
            retry.c[i]=retry.c[i]/scale
        }
    }
    scale=max(abs(retry.c))
    if (scale>0 & scale<.) retry.c=retry.c/scale
    if (hasmissing(retry.A) | hasmissing(retry.b) | hasmissing(retry.c) | hasmissing(column_scale)) return(first)
    if (trace=="trace") printf("{txt}RSM restart: fresh basis, row/column/objective equilibration, Bland\n")
    second=_rsm_precision_once(retry,tol,trace,1,second_basis,certified,allow_rebuild,second_pivot)
    pivot_seen=pivot_seen | second_pivot
    cycle_seen=cycle_seen | (second.msg=="Repeated Bland pivot state")
    second.iters=first.iters+second.iters
    if (second.status==0) {
        second.x=second.x:*column_scale'
        check=_dea_lpcheck(second.x',original.A[.,1..original.n_orig],original.b,
            J(0,original.n_orig,0),J(0,1,0),J(1,original.n_orig,0),J(1,original.n_orig,.))
        second.objval=original.c[1..original.n_orig]*second.x
        if (_dea_lpcheck_ok(check,tol) & second.objval<.) {
            certbasis=second_basis
            second.msg="Optimal after equilibrated transition restart"
            prob=original
            return(second)
        }
    }
    // A numerical restart is not an infeasibility/unboundedness certificate.
    first.iters=second.iters
    return(first)
}

struct RSM_Result scalar _rsm_precision_once(struct RSM_Prob scalar prob,
                                         real scalar tol,
                                         string scalar trace,
                                         real scalar restart, real colvector certbasis, real scalar certified,
                                         | real scalar allow_rebuild, real scalar pivot_seen)
{
    real scalar i, max_iter, p1_iter
    real colvector xb
    real rowvector c_orig
    real scalar art_sum
    real matrix B0, original_A
    real colvector original_b, feasible_x
    real rowvector check

    struct RSM_Basis scalar basis
    struct RSM_Result scalar result

    if(args()<7) allow_rebuild=0
    pivot_seen=0
    certbasis=J(prob.m_con,1,.)
    max_iter = max((10000, 20 * prob.m_con))

    // Save Phase 2 objective
    c_orig = prob.c
    original_A = prob.A[.,1..prob.n_orig]
    original_b = prob.b

    // =================================================================
    // PHASE 1: Minimize sum of artificials
    // =================================================================
    prob.c = prob.is_artif

    // Initialize basis
    basis.basic = prob.init_basic
    basis.n_upd = 0

    // Initial B^{-1} (should be identity by construction)
    B0 = prob.A[., basis.basic]
    if (all(B0 :== I(prob.m_con))) {
        basis.Binv = I(prob.m_con)
    }
    else {
        basis.Binv = _rsm_inv(B0)
        if (hasmissing(basis.Binv)) {
            result.status = 3
            result.msg = "Invalid initial basis inverse"
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
    result = _rsm_precision_solve(prob, basis, tol, max_iter, 1, trace, restart, certified)
    p1_iter = result.iters

    if (result.status == 2) {
        result.status = 3
        result.msg = "Phase1 unbounded; no feasibility certificate"
        prob.c = c_orig
        return(result)
    }

    // Enable only on a final numerical retry, after every preceding
    // acceptance path failed. Keep all original rows and require a checked
    // artificial-free basis; the caller still checks quality in original units.
    pivot_seen=(result.status==3 & result.msg=="No stable factorized pivot")
    if(allow_rebuild & pivot_seen & p1_iter<max_iter) {
        if(_rsm_retry_realbasis(prob,basis,tol)) {
            result.status=0
            result.msg="Phase1 original-column reconstruction"
        }
    }

    if (result.status != 0) {
        prob.c = c_orig
        return(result)
    }

    // Check artificial sum
    xb = _rsm_factor_solve(prob.A[,basis.basic],prob.b)
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

    // A tiny artificial objective alone is not a feasible original LP.
    feasible_x=J(prob.n_orig,1,0)
    for (i=1;i<=prob.m_con;i++) {
        if (basis.basic[i]<=prob.n_orig) feasible_x[basis.basic[i]]=max((0,xb[i]))
    }
    check=_dea_lpcheck(feasible_x',original_A,original_b,
        J(0,prob.n_orig,0),J(0,1,0),J(1,prob.n_orig,0),J(1,prob.n_orig,.))
    if (!_dea_lpcheck_ok(check,tol)) {
        result.status=3; result.msg="Residual check failed"
        result.x=feasible_x; result.objval=.
        prob.c=c_orig
        return(result)
    }

    // =================================================================
    // TRANSITION: Purge artificials, drop columns, restore objective
    // =================================================================
    if (any(basis.basic:>prob.n_orig)) {
        if (!_rsm_retry_realbasis(prob,basis,tol)) {
            result.status=3;result.msg="Transition retry requires artificial-free basis";result.objval=.
            return(result)
        }
    }
    if (!_rsm_purge_artif(prob, basis, tol)) {
        result.status=3; result.msg="Artificial-basis pivot failed"
        result.objval=.; result.x=feasible_x
        return(result)
    }
    if (!_rsm_drop_artif(prob, basis)) {
        result.status=3; result.msg="Artificial-basis transition failed"
        result.objval=.; result.x=feasible_x
        return(result)
    }
    prob.c = c_orig[1..prob.n_var]

    // Failure-only retry: no basic artificial and factorized Phase 2 uses no inverse.
    // No artificial-row removal or pivot tolerance change is permitted here.

    // =================================================================
    // PHASE 2: Optimize original objective
    // =================================================================
    if (trace == "trace") {
        printf("{txt}Phase2: %g con, %g var\n", prob.m_con, prob.n_var)
    }

    result = _rsm_precision_solve(prob, basis, tol, max_iter - p1_iter, 2, trace, restart, certified)
    result.iters = result.iters + p1_iter

    // Extract full solution
    if (result.status == 0) {
        result.x = J(prob.n_var, 1, 0)
        xb = _rsm_factor_solve(prob.A[,basis.basic],prob.b)
        for (i = 1; i <= prob.m_con; i++) {
            result.x[basis.basic[i]] = max((0, xb[i]))
        }
        certbasis=basis.basic
        result.objval = prob.c * result.x
        check=_dea_lpcheck(result.x',original_A,original_b,
            J(0,prob.n_var,0),J(0,1,0),J(1,prob.n_var,0),J(1,prob.n_var,.))
        if (!_dea_lpcheck_ok(check,tol) | result.objval>=.) {
            result.status=3; result.msg="Residual check failed"; result.objval=.
        }
    }

    return(result)
}

struct RSM_Result scalar _rsm_precision_solve(struct RSM_Prob scalar prob,
    struct RSM_Basis scalar basis,real scalar tol,real scalar max_iter,
    real scalar phase,string scalar trace,real scalar restart,real scalar certified)
{
    struct RSM_Result scalar result
    real scalar iter,entering,pivot_row,use_bland,stall_count,consec_degen,prev_obj,obj_val,i
    real colvector xb,d
    real matrix seen
    real rowvector state
    real scalar k,hit
    result.status=3;result.objval=.;result.x=J(prob.n_var,1,.);result.iters=0
    result.msg="Maximum iterations";use_bland=restart;stall_count=consec_degen=0;prev_obj=.
    seen=J(0,3*prob.m_con+4,0)
    for(iter=1;iter<=max_iter;iter++) {
        result.iters=iter
        xb=_rsm_factor_solve(prob.A[,basis.basic],prob.b)
        if(hasmissing(xb)) {result.msg="Factorized basic solution unavailable";return(result);}
        if(min(xb)<-min((_dea_tol_feas(),tol))) {result.msg="Factorized primal basis lost feasibility";return(result);}
        obj_val=prob.c[basis.basic']*xb
        if(prev_obj<. & abs(obj_val-prev_obj)<_dea_tol_zero()) {
            stall_count++
            if(stall_count>prob.m_con) use_bland=1
        }
        else stall_count=0
        prev_obj=obj_val
        if(certified) entering=_rsm_cert_enter(prob,basis,use_bland,tol)
        else entering=_rsm_fp_enter(prob,basis,use_bland,tol)
        if(entering<0) {result.msg="Factorized price unavailable";return(result);}
        if(entering==0) {
            // Original Phase-1 transition still requires a checked inverse.
            // The failure-only retry does not change dependent-row removal.
            if(phase==1 & any(basis.basic:>prob.n_orig)) {
                basis.Binv=_rsm_inv(prob.A[,basis.basic]);basis.n_upd=0
                if(hasmissing(basis.Binv)) {result.msg="Phase1 transition inverse unavailable";return(result);}
            }
            result.status=0;result.objval=obj_val;result.msg=sprintf("Factorized optimal P%g i%g",phase,iter)
            result.x=J(prob.n_var,1,0)
            for(i=1;i<=prob.m_con;i++) result.x[basis.basic[i]]=max((0,xb[i]))
            return(result)
        }
        d=_rsm_factor_solve(prob.A[,basis.basic],prob.A[,entering])
        if(hasmissing(d)) {result.msg="Factorized pivot unavailable";return(result);}
        pivot_row=_rsm_retry_exchange(prob,basis,xb,d,entering,use_bland,tol)
        if(pivot_row==0) {
            if(rows(d)>0 & max(d)>0) {result.msg="No stable factorized pivot";return(result);}
            result.status=2;result.objval=obj_val;result.msg="Factorized unbounded";return(result)
        }
        // A repeated deterministic Bland pivot state cannot yield a later
        // success within this fixed phase/restart. End only this path so the
        // existing subsequent restart retains priority and its chosen weights.
        if(!certified & phase==2 & use_bland==1) {
            state=(basis.basic',xb',d',obj_val,entering,pivot_row,use_bland)
            hit=0
            for(k=1;k<=rows(seen);k++) if(all(seen[k,.]:==state)) {hit=1;break;}
            if(hit) {result.msg="Repeated Bland pivot state";return(result);}
            seen=seen\state
            if(rows(seen)>8) seen=seen[2..rows(seen),.]
        }
        if(abs(xb[pivot_row])<_dea_tol_feas()) {
            consec_degen++
            if(consec_degen>2*prob.m_con) use_bland=1
        }
        else consec_degen=0
        if(trace=="trace") printf("FS P%g i%g enter%g leave%g step%e pivot%e\\n",phase,iter,entering,basis.basic[pivot_row],max((0,xb[pivot_row]))/d[pivot_row],d[pivot_row])
        basis.basic[pivot_row]=entering
    }
    return(result)
}

real scalar _rsm_precision_accept(struct RSM_Prob scalar original,
    struct RSM_Result scalar candidate,real colvector certbasis,real scalar tol)
{
    real scalar valid_basis,i,accept
    real matrix B,XP,YP
    real rowvector check,quality
    if(candidate.status!=0 | rows(candidate.x)!=original.n_orig) return(0)
    valid_basis=(rows(certbasis)==original.m_con & !hasmissing(certbasis))
    if(valid_basis) {
        for(i=1;i<=rows(certbasis);i++) {
            if(certbasis[i]<1 | certbasis[i]>original.n_orig | certbasis[i]!=floor(certbasis[i])) valid_basis=0
            if(sum(certbasis:==certbasis[i])!=1) valid_basis=0
        }
    }
    XP=YP=J(original.m_con,2,.)
    if(valid_basis) {
        B=original.A[,certbasis]
        XP=_rsm_retry_refine(B,original.b)
        YP=_rsm_retry_refine(B',original.c[certbasis']')
        if(!hasmissing(XP) & !hasmissing(YP)) {
            candidate.x=J(original.n_orig,1,0)
            for(i=1;i<=rows(certbasis);i++) candidate.x[certbasis[i]]=max((0,XP[i,1]+XP[i,2]))
        }
    }
    check=_dea_lpcheck(candidate.x',original.A[,1..original.n_orig],original.b,
        J(0,original.n_orig,0),J(0,1,0),J(1,original.n_orig,0),J(1,original.n_orig,.))
    candidate.objval=original.c[1..original.n_orig]*candidate.x
    quality=_rsm_retry_quality(original,candidate.x,YP[,1]',YP[,2]',tol)
    accept=(_dea_lpcheck_ok(check,tol) & candidate.objval<. & valid_basis & !hasmissing(XP) & !hasmissing(YP) & !hasmissing(quality) & quality[1]==1)
    return(accept)
}


// Exact binary-power duplicates after all earlier numerical paths fail.
// No approximate rank test, no tolerance-based deletion, no integer aggregation.
real scalar _rsm_duplicate_ratio(real rowvector a,real rowvector b)
{
    real scalar i,k,ratio,exponent,scale
    real rowvector product
    if(cols(a)<1 | cols(a)!=cols(b) | hasmissing(a) | hasmissing(b)) return(.)
    k=0
    for(i=1;i<=cols(b);i++) if(b[i]!=0) {k=i;break;}
    if(k==0) return(.) // Zero rows/columns are outside this narrow transform.
    ratio=a[k]/b[k]
    if(missing(ratio) | ratio<=0) return(.)
    exponent=round(ln(ratio)/ln(2))
    if(exponent< -400 | exponent>400) return(.)
    // Mata exponentiation is approximate at large integer exponents.
    // Binary shifts by2 or.5 stay exact and normal over this bounded range.
    scale=1
    for(i=1;i<=abs(exponent);i++) scale=scale*(exponent<0 ? .5 : 2)
    if(ratio!=scale) return(.)
    for(i=1;i<=cols(a);i++) {
        // EFT proves equality of the actual binary64 values, not rounded ratios.
        product=_rsm_q_product(b[i],scale)
        if(hasmissing(product) | product[2]!=0 | product[1]!=a[i]) return(.)
    }
    return(scale)
}

real scalar _rsm_duplicate_map(struct RSM_Prob scalar original,
    real colvector kr,real colvector kc,real colvector representative,
    real scalar work,| real scalar budget)
{
    real scalar m,n,i,j,hit,ratio
    real matrix A
    real rowvector c
    if(args()<6) budget=1000000
    kr=kc=J(0,1,.);representative=J(0,1,.);work=0
    m=original.m_con;n=original.n_orig
    if(missing(m) | missing(n) | m<1 | n<1 | m!=floor(m) | n!=floor(n) |
       n+2*m>4096 | budget<0 | missing(budget) | rows(original.A)!=m |
       cols(original.A)<n | rows(original.b)!=m | cols(original.b)!=1 |
       cols(original.c)<n) return(0)
    A=original.A[,1..n];c=original.c[1..n]
    if(hasmissing(A) | hasmissing(original.b) | hasmissing(c)) return(0)
    representative=J(n,1,.)
    for(i=1;i<=m;i++) {
        hit=0
        for(j=1;j<=rows(kr);j++) {
            if(work+n+1>budget) return(0)
            work=work+n+1
            ratio=_rsm_duplicate_ratio((A[i,.],original.b[i]),(A[kr[j],.],original.b[kr[j]]))
            if(ratio<.) {hit=1;break;}
        }
        if(!hit) kr=kr\i
    }
    for(i=1;i<=n;i++) {
        hit=0
        for(j=1;j<=rows(kc);j++) {
            if(work+m+1>budget) return(0)
            work=work+m+1
            ratio=_rsm_duplicate_ratio((A[,i]' ,c[i]),(A[,kc[j]]',c[kc[j]]))
            if(ratio<.) {hit=kc[j];break;}
        }
        if(!hit) {kc=kc\i;hit=i;}
        representative[i]=hit
    }
    return(rows(kr)<m | rows(kc)<n)
}

// Preserve original slack/artificial initialization and all bound rows.
// The core receives a continuous nonnegative standard-form LP. Finite bounds
// have already become rows; MILP branching always sees restored coordinates.
real scalar _rsm_duplicate_problem(struct RSM_Prob scalar original,
    real colvector kr,real colvector kc,real colvector representative,
    struct RSM_Prob scalar reduced)
{
    real scalar i,j,k,n
    real colvector ac,allcols,found
    n=original.n_orig;ac=J(0,1,.)
    for(j=n+1;j<=original.n_var;j++) if(any(original.A[kr,j]:!=0)) ac=ac\j
    allcols=kc\ac
    reduced.m_con=rows(kr);reduced.n_orig=rows(kc);reduced.n_var=rows(allcols)
    reduced.A=original.A[kr,allcols];reduced.b=original.b[kr]
    reduced.c=original.c[allcols'];reduced.is_artif=original.is_artif[allcols']
    reduced.init_basic=J(rows(kr),1,.)
    for(i=1;i<=rows(kr);i++) {
        k=original.init_basic[kr[i]]
        if(k<1 | k>original.n_var | missing(k) | k!=floor(k)) return(0)
        if(k<=n) k=representative[k]
        found=selectindex(allcols:==k)
        if(rows(found)!=1) return(0)
        reduced.init_basic[i]=found[1]
    }
    return(1)
}

real scalar _rsm_duplicate_accept(struct RSM_Prob scalar original,
    struct RSM_Prob scalar reduced,struct RSM_Result scalar candidate,
    real colvector kr,real colvector kc,real scalar tol)
{
    real scalar i,m,n
    real colvector basic,x
    real matrix B,XP,YP,Y
    real rowvector check,quality
    m=reduced.m_con;n=reduced.n_orig
    if(candidate.status!=0 | rows(candidate.x)!=n | hasmissing(candidate.x)) return(0)
    // Require a support-sized nonsingular basis for this bounded recovery.
    // Degenerate/nonbasic approximate supports are conservatively left unresolved.
    basic=selectindex(candidate.x:>0)
    if(rows(basic)!=m) return(0)
    B=reduced.A[,basic]
    XP=_rsm_retry_refine(B,reduced.b)
    YP=_rsm_retry_refine(B',reduced.c[basic']')
    if(hasmissing(XP) | hasmissing(YP)) return(0)
    x=J(original.n_orig,1,0);Y=J(original.m_con,2,0)
    for(i=1;i<=m;i++) x[kc[basic[i]]]=max((0,XP[i,1]+XP[i,2]))
    Y[kr,.]=YP
    check=_dea_lpcheck(x',original.A[,1..original.n_orig],original.b,
        J(0,original.n_orig,0),J(0,1,0),J(1,original.n_orig,0),J(1,original.n_orig,.))
    quality=_rsm_retry_quality(original,x,Y[,1]',Y[,2]',tol)
    if(!_dea_lpcheck_ok(check,tol) | hasmissing(quality) | quality[1]!=1) return(0)
    candidate.basiscols=kc[basic]
    candidate.x=x;candidate.objval=original.c[1..original.n_orig]*x
    return(candidate.objval<.)
}

struct RSM_Result scalar _rsm_duplicate_retry(struct RSM_Prob scalar original,
    real scalar tol,string scalar trace,| real scalar check_success)
{
    struct RSM_Prob scalar reduced,saved
    struct RSM_Result scalar result
    real colvector kr,kc,representative
    real scalar work
    if(args()<4) check_success=1
    result.status=3;result.objval=.;result.x=J(original.n_orig,1,.)
    result.iters=0;result.msg="Exact duplicate recovery unavailable"
    if(!_rsm_duplicate_map(original,kr,kc,representative,work)) return(result)
    if(!_rsm_duplicate_problem(original,kr,kc,representative,reduced)) return(result)
    saved=reduced
    // One reduced run, with duplicate retry disabled and no nested audit emit.
    // Existing per-path iteration limits remain; map comparisons have a fixed cap.
    result=_rsm_two_phase_core(reduced,tol,trace,0,check_success,0)
    if(_rsm_duplicate_accept(original,saved,result,kr,kc,tol)) {
        result.msg="Optimal after exact duplicate recovery"
        return(result)
    }
    result.status=3;result.objval=.;result.x=J(original.n_orig,1,.)
    result.msg="Exact duplicate recovery unavailable"
    return(result)
}

// T31 sufficient screen: conventional finite normal binary64 arithmetic only.
// eta=64*(m+n+4)*eps dominates dot, norm, and positive aggregation roundoff.
// Rejection/unavailability falls through to the existing EFT certificate.
// T35 isolated candidate: bounded exact expansions, not a tolerance rank test.
// Components are accumulated low-to-high by error-free TwoSum/TwoProduct.
// No component is discarded; capacity/range/work exhaustion means unavailable.
real rowvector _rsm_xgrow(real rowvector e,real scalar b,
 real scalar work,real scalar err,real scalar budget)
{
 real scalar i,Q
 real rowvector h,s
 if(err) return(J(1,0,.))
 if(b==0) return(e)
 h=J(1,0,.);Q=b
 for(i=1;i<=cols(e);i++) {
  if(work+1>budget) {err=4;return(J(1,0,.));}
  work=work+1
  if(hasmissing((Q,e[i])) | max(abs((Q,e[i])))>8.452712498170644e270) {err=2;return(J(1,0,.));}
  s=_rsm_q_sum(Q,e[i])
  if(hasmissing(s)) {err=2;return(J(1,0,.));}
  if(any((s:!=0):&(abs(s):<9.3326361850321888e-302))) {err=2;return(J(1,0,.));}
  if(s[2]!=0) h=(h,s[2])
  Q=s[1]
 }
 if(Q!=0 | cols(h)==0) h=(h,Q)
 if(cols(h)>32) {err=3;return(J(1,0,.));}
 return(h)
}
real rowvector _rsm_xadd(real rowvector a,real rowvector b,
 real scalar work,real scalar err,real scalar budget)
{
 real scalar i
 real rowvector out
 out=a
 for(i=1;i<=cols(b);i++) {
  out=_rsm_xgrow(out,b[i],work,err,budget)
  if(err) return(J(1,0,.))
 }
 return(out)
}
real rowvector _rsm_xmul(real rowvector a,real rowvector b,
 real scalar work,real scalar err,real scalar budget)
{
 real scalar i,j
 real rowvector out,p
 if(err) return(J(1,0,.))
 out=0
 for(i=1;i<=cols(a);i++) for(j=1;j<=cols(b);j++) {
  if(a[i]==0 | b[j]==0) continue
  if(work+1>budget) {err=4;return(J(1,0,.));}
  work=work+1;p=_rsm_q_product(a[i],b[j])
  if(hasmissing(p)) {err=2;return(J(1,0,.));}
  out=_rsm_xgrow(out,p[2],work,err,budget)
  out=_rsm_xgrow(out,p[1],work,err,budget)
  if(err) return(J(1,0,.))
 }
 return(out)
}
real scalar _rsm_xnormalize(real matrix T,
 real scalar work,real scalar err,real scalar budget)
{
 real scalar mx,e,i,s,cost
 real matrix V
 if(err | hasmissing(T)) {if(!err) err=2;return(0);}
 mx=max(abs(vec(T)))
 if(mx==0) return(1)
 e=ceil(ln(mx)/ln(2))
 if(missing(e) | abs(e)>900) {err=2;return(0);}
 cost=4*rows(T)*cols(T)+abs(e)
 if(work+cost>budget) {err=4;return(0);}
 work=work+cost;s=1
 for(i=1;i<=abs(e);i++) s=s*(e>0 ? .5 : 2)
 V=T*s
 if(hasmissing(V) | any(vec((V:==0):&(T:!=0))) | any(vec((V:/s):!=T))) {err=5;return(0);}
 if(any((vec(V):!=0):&(abs(vec(V)):<3.8725919148493183e-121)) |
    max(abs(vec(V)))>2.5822498780869086e120) {err=2;return(0);}
 T=V
 return(1)
}

// Fraction-free elimination on [A|b], retaining original row/column indices.
// A new pivot is selected only from an exact nonzero expansion. An omitted
// row has every original coefficient AND its RHS exactly zero after invertible
// row operations. Contradictions and almost-dependent independent rows survive.
// detail: available,reason,rank,work,max expansion length.
real scalar _rsm_rowmap(struct RSM_Prob scalar original,
 real colvector kr,real colvector kc,real rowvector detail,| real scalar budget)
{
 real scalar m,n,i,j,l,col,k,work,err,maxlen
 real matrix piv,T,P,N
 real rowvector a,d,u,v,w
 if(args()<5) budget=1000000
 kr=kc=J(0,1,.);detail=(0,1,0,0,0);work=err=maxlen=0
 m=original.m_con;n=original.n_orig
 if(missing(m) | missing(n) | m<0 | n<1 | m!=floor(m) | n!=floor(n) |
  m>256 | n>256 | n+2*m>4096 | missing(budget) | budget<0 |
  rows(original.A)!=m | cols(original.A)<n | rows(original.b)!=m |
  cols(original.b)!=1 | hasmissing(original.A[,1..n]) | hasmissing(original.b)) return(0)
 piv=J(0,32,0);k=0
 for(i=1;i<=m;i++) {
  T=J(n+1,32,0);T[,1]=(original.A[i,1..n],original.b[i])'
  if(!_rsm_xnormalize(T,work,err,budget)) {detail=(0,10+err,k,work,maxlen);return(0);}
  for(j=1;j<=k;j++) {
   a=select(T[kc[j],.],T[kc[j],.]:!=0)
   if(cols(a)==0) continue
   P=piv[| (j-1)*(n+1)+1,1 \ j*(n+1),32 |]
   d=select(P[kc[j],.],P[kc[j],.]:!=0)
   N=J(n+1,32,0)
   for(l=1;l<=n+1;l++) {
    u=_rsm_xmul(d,select(T[l,.],T[l,.]:!=0),work,err,budget)
    v=_rsm_xmul(a,select(P[l,.],P[l,.]:!=0),work,err,budget)
    w=_rsm_xadd(u,-v,work,err,budget)
    if(err) {detail=(0,10+err,k,work,maxlen);return(0);}
    maxlen=max((maxlen,cols(w)));N[l,1..cols(w)]=w
   }
   if(any(N[kc[j],.]:!=0)) {detail=(0,15,k,work,maxlen);return(0);}
   T=N
   if(!_rsm_xnormalize(T,work,err,budget)) {detail=(0,10+err,k,work,maxlen);return(0);}
  }
  col=0
  for(l=1;l<=n;l++) if(any(T[l,.]:!=0)) {col=l;break;}
  if(col==0) {
   if(any(T[n+1,.]:!=0)) {detail=(0,2,k,work,maxlen);return(0);}
   continue
  }
  if(k>=32) {detail=(0,3,k,work,maxlen);return(0);}
  k=k+1;kr=kr\i;kc=kc\col;piv=piv\T
 }
 detail=(1,0,k,work,maxlen)
 return(1)
}

// Empty-basis proof is coordinatewise: A=0,b=0,c>=0,x>=L and
// c_j>0 => x_j=L_j. It retains alternative optimal zero-cost coordinates.
real rowvector _rsm_row_empty(struct RSM_Prob scalar original,
 struct RSM_Result scalar candidate,real scalar tol,real colvector L)
{
 real scalar m,n,j,optimal
 real colvector x
 real rowvector meta,c
 m=original.m_con;n=original.n_orig
 meta=(-1,2,m,0,n,0,0,0,.,.,.,.,.,.)
 if(rows(candidate.x)!=n | hasmissing(candidate.x) | rows(L)!=n | hasmissing(L) |
  rows(original.A)!=m | cols(original.A)<n | rows(original.b)!=m |
  hasmissing(original.A) | hasmissing(original.b) | cols(original.c)<n) return(meta)
 if(any(vec(original.A[,1..n]):!=0) | any(original.b:!=0)) return(meta)
 c=original.c[1..n]
 if(hasmissing(c) | any(c:<0)) {meta[1]=0;return(meta);}
 optimal=all(candidate.x:>=L)
 for(j=1;j<=n;j++) if(c[j]>0 & candidate.x[j]!=L[j]) optimal=0
 if(optimal & candidate.objval<. & candidate.objval==c*candidate.x) {
  meta[1]=1;meta[2]=0;meta[9]=1;meta[11]=meta[12]=0;return(meta)
 }
 x=candidate.x
 for(j=1;j<=n;j++) if(c[j]>0 | x[j]<L[j]) x[j]=L[j]
 if(c*x>=.) return(meta)
 candidate.x=x;candidate.objval=c*x;candidate.msg="Optimal after empty-basis success refinement"
 meta[1]=2;meta[2]=0;meta[10]=1;meta[13]=meta[14]=0
 return(meta)
}

real scalar _rsm_row_accept(struct RSM_Prob scalar original,
 struct RSM_Result scalar candidate,real colvector kr,real scalar tol)
{
 real rowvector meta,detail
 real colvector checkrows,checkcols
 if(candidate.status!=0 | hasmissing(candidate.x) | rows(candidate.x)!=original.n_orig) return(0)
 if(!_rsm_rowmap(original,checkrows,checkcols,detail)) return(0)
 if(rows(kr)!=rows(checkrows)) return(0)
 if(rows(kr)>0) if(any(kr:!=checkrows)) return(0)
 if(rows(candidate.basiscols)!=rows(kr)) return(0)
 meta=_rsm_success_check(original,candidate,candidate.basiscols,tol)
 return(meta[1]==1 | meta[1]==2)
}

struct RSM_Result scalar _rsm_row_retry(struct RSM_Prob scalar original,
 real scalar tol,string scalar trace,real scalar check_success)
{
 struct RSM_Prob scalar reduced
 struct RSM_Result scalar result
 real colvector kr,kc,allcols
 real rowvector detail
 result.status=3;result.objval=.;result.x=J(original.n_orig,1,.)
 result.iters=0;result.msg="Certified row recovery unavailable"
 if(!_rsm_rowmap(original,kr,kc,detail)) return(result)
 if(rows(kr)>=original.m_con) return(result)
 allcols=(1::original.n_orig)
 if(!_rsm_duplicate_problem(original,kr,allcols,allcols,reduced)) return(result)
 // Exactly one reduced core call; no duplicate/row recursion or audit emit.
 result=_rsm_two_phase_core(reduced,tol,trace,0,check_success,0)
 if(_rsm_row_accept(original,result,kr,tol)) {
  result.msg="Optimal after certified row recovery"
  return(result)
 }
 result.status=3;result.objval=.;result.x=J(original.n_orig,1,.)
 result.msg="Certified row recovery unavailable"
 return(result)
}

real rowvector _rsm_success_fastread()
{
 external real rowvector _rsm_success_fastcounts
 return(_rsm_success_fastcounts)
}
real rowvector _rsm_quality_fast(struct RSM_Prob scalar P,
 real colvector x,real rowvector y,real colvector L,real scalar tol)
{
 real scalar m,n,eta,o,oe,pu,cu,du,gu,cutoff,price,cs,den
 real matrix A
 real rowvector c,r,er,z,cmp
 real colvector p,ep,d
 m=P.m_con;n=P.n_orig
 if(m<1 | n<1 | n+2*m>4096 | rows(x)!=n | rows(L)!=n | cols(y)!=m |
    missing(tol) | tol<=0) return(J(1,16,.))
 A=P.A[,1..n];c=P.c[1..n]
 z=abs((vec(A)',P.b',c,x',y,L'))
 if(hasmissing(z)) return(J(1,16,.))
 z=select(z,z:!=0)
 // Narrower than EFT: triple products, aggregation and error radii stay normal.
 // Outside [2^-100,2^100], use the existing [2^-400,2^400] guarded EFT path.
 if(cols(z)>0) if(min(z)<7.8886090522101181e-31 | max(z)>1.2676506002282294e30) return(J(1,16,.))
 eta=64*(m+n+4)*epsilon(1)
 p=P.b-A*x;ep=eta*(abs(P.b)+abs(A)*abs(x))
 r=c-y*A;er=eta*(abs(c)+abs(y)*abs(A))
 d=abs(x-L)+eta*(abs(x)+abs(L))
 d= d :* (x:!=L) // Equal stored endpoints have exactly zero displacement.
 o=c*x;oe=eta*(abs(c)*abs(x))
 pu=(abs(y)*(abs(p)+ep))*(1+eta)
 cu=((abs(r)+er)*d)*(1+eta)
 du=max((0,max(-r+er)))*(1+eta)
 gu=(pu+cu)*(1+eta) // Gap identity includes the exact lower-bound offset.
 cutoff=min((1e-8,tol));price=min((_dea_tol_opt(),tol));cs=min((1,max(abs(c))))
 if(hasmissing((o,oe,pu,cu,du,gu)) | abs(o)<oe) return(J(1,16,.))
 cmp=(_rsm_q_scaled(gu,cutoff,(o,oe)),_rsm_q_scaled(pu,cutoff,(o,oe)),
      _rsm_q_scaled(cu,cutoff,(o,oe)),_rsm_q_compare(du,price,cs))
 den=1+max((0,abs(o)-oe))
 // The third field is an approximate dual value; acceptance uses only bounds.
 return((!hasmissing(cmp) & max(cmp)<=0,o,o-(r*(x-L)-y*p),gu/den,pu/den,cu/den,du,price*cs,den,gu,pu,cu,cmp))
}

// Exact-coordinate affine certificate: x>=L, A*x=b, min c*x.
// D=b'y+L'(c-A'y). Hence c'x-D=(c-A'y)'(x-L)-y'(b-Ax).
// Bounding both terms avoids rounding a lower-bound objective constant into c'x.
real rowvector _rsm_affine_quality(struct RSM_Prob scalar P,
 real colvector x,real rowvector yh,real rowvector yl,real colvector L,real scalar tol)
{
 real scalar m,n,i,j,pu,cu,du,gu,den,ayu,rpu,rdu,deltau,cutoff,price,cs,gcenter
 real matrix A
 real rowvector c,O,R,Q,D,cmp
 m=P.m_con;n=P.n_orig
 if(m<1 | n<1 | n+2*m>4096 | rows(x)!=n | rows(L)!=n | cols(yh)!=m |
    cols(yl)!=m | missing(tol) | tol<=0 | hasmissing((x,L)) |
    hasmissing((yh,yl))) return(J(1,16,.))
 A=P.A[,1..n];c=P.c[1..n];O=_rsm_q_dot(-c,x,0)
 if(hasmissing(O)) return(J(1,16,.))
 pu=cu=du=gcenter=0
 for(i=1;i<=m;i++) {
  Q=_rsm_q_dot(A[i,.],x,P.b[i]);if(hasmissing(Q)) return(J(1,16,.))
  ayu=_rsm_q_add(abs(yh[i]),abs(yl[i]),1)
  rpu=_rsm_q_add(abs(Q[1]),Q[2],1)
  pu=_rsm_q_add(pu,_rsm_q_mulup(ayu,rpu),1)
  gcenter=gcenter-(yh[i]+yl[i])*Q[1]
 }
 for(j=1;j<=n;j++) {
  R=_rsm_q_dot((yh,yl),(A[,j]\A[,j]),c[j])
  D=_rsm_q_dot(1,L[j],x[j])
  if(hasmissing((R,D))) return(J(1,16,.))
  rdu=_rsm_q_add(abs(R[1]),R[2],1)
  deltau=_rsm_q_add(abs(D[1]),D[2],1)
  cu=_rsm_q_add(cu,_rsm_q_mulup(deltau,rdu),1)
  du=max((du,-_rsm_q_add(R[1],-R[2],0)))
  gcenter=gcenter+R[1]*D[1]
 }
 gu=_rsm_q_add(pu,cu,1)
 den=_rsm_q_add(1,max((0,_rsm_q_add(abs(O[1]),-O[2],0))),0)
 cutoff=min((1e-8,tol));price=min((_dea_tol_opt(),tol));cs=min((1,max(abs(c))))
 if(hasmissing((gu,pu,cu,du,den,gcenter))) return(J(1,16,.))
 cmp=(_rsm_q_scaled(gu,cutoff,O),_rsm_q_scaled(pu,cutoff,O),
      _rsm_q_scaled(cu,cutoff,O),_rsm_q_compare(du,price,cs))
 return((!hasmissing(cmp) & max(cmp)<=0,O[1],O[1]-gcenter,gu/den,pu/den,cu/den,du,price*cs,den,gu,pu,cu,cmp))
}

// Refine actual basic coordinates with every nonbasic coordinate fixed at L.
// Residuals are always evaluated in the original unshifted rows, including RHS.
real matrix _rsm_affine_refine(real matrix A,real colvector b,
 real colvector basic,real colvector L)
{
 real scalar iteration,err,old_err,m,n
 real matrix B,best
 real colvector hi,lo,xh,xl,R,correction,S,V,E,T,H,Q,comparison
 m=rows(A);n=cols(A);B=A[,basic]
 if(rows(basic)!=m | rows(L)!=n) return(J(m,2,.))
 hi=L[basic]+_rsm_factor_solve(B,b-A*L);lo=J(m,1,0)
 if(hasmissing(hi)) return(J(m,2,.))
 best=(hi,lo);old_err=.
 xh=L;xh[basic]=hi
 comparison=abs(b)+abs(A)*(abs(xh):+1):+1e-300
 for(iteration=0;iteration<=8;iteration++) {
  xh=L;xh[basic]=hi;xl=J(n,1,0);xl[basic]=lo
  R=_rsm_inverse_residual((A,A),(xh\xl),b)
  if(hasmissing(R)) return(J(m,2,.))
  err=max(abs(R):/comparison)
  if(err<old_err) {best=(hi,lo);old_err=err;}
  else break
  if(iteration==8 | err==0) break
  correction=_rsm_factor_solve(B,R);if(hasmissing(correction)) break
  S=hi+correction;V=S-hi;E=(hi-(S-V))+(correction-V)
  T=lo+E;H=S+T;V=H-S;Q=(S-(H-V))+(T-V)
  if(hasmissing(H) | hasmissing(Q)) break
  hi=H;lo=Q
 }
 return(best)
}

// T31 research-only accepted-result screen. No public option or version change.
// Actual Phase-2 column indices are passed through; zero basic values are valid.
void _rsm_success_reset(real scalar policy,real scalar recording,| real scalar fast)
{
 external real scalar _rsm_success_policy,_rsm_success_recording
 external real matrix _rsm_success_records
 external real rowvector _rsm_success_counts
 external real scalar _rsm_success_fast
 external real rowvector _rsm_success_fastcounts
 if(args()<3) fast=1
 _rsm_success_fast=fast;_rsm_success_fastcounts=J(1,4,0)
 _rsm_success_policy=policy;_rsm_success_recording=recording
 _rsm_success_records=J(0,14,.);_rsm_success_counts=J(1,8,0)
}
real matrix _rsm_success_read()
{
 external real matrix _rsm_success_records
 return(_rsm_success_records)
}
real rowvector _rsm_success_counts_read()
{
 external real rowvector _rsm_success_counts
 return(_rsm_success_counts)
}

// meta: outcome,reason,m,basis_size,n,map_work,dual_solves,primal_solves,
// old_quality,new_quality,old_gap,old_weighted_primal,new_gap,new_weighted_primal.
// outcome: 1=preserved,2=refined,0=quality failed,-1=unavailable.
real rowvector _rsm_success_check(struct RSM_Prob scalar original,
    struct RSM_Result scalar candidate,real colvector basic,real scalar tol,| real colvector L)
{
 external real scalar _rsm_success_fast
 external real rowvector _rsm_success_fastcounts
 real scalar m,n,k,i,j,hit,work,ratio,affine
 real colvector kr,x
 real matrix B,YP,Y,XP
 real colvector yf
 real rowvector meta,q,check,rowdetail
 real colvector rowcols
 m=original.m_con;n=original.n_orig;k=rows(basic)
 affine=(args()==5)
 if(!affine) L=J(n,1,0)
 meta=(-1,1,m,k,n,0,0,0,.,.,.,.,.,.)
 if(m<0 | n<1 | n+2*m>4096 | missing(tol) | tol<=0) return(meta)
 meta[2]=2
 if(rows(L)!=n | hasmissing(L)) return(meta)
 if(k==0) return(_rsm_row_empty(original,candidate,tol,L))
 if(k<1 | k>m | k>n | hasmissing(basic) | rows(candidate.x)!=n) return(meta)
 for(i=1;i<=k;i++) {
  if(basic[i]<1 | basic[i]>n | basic[i]!=floor(basic[i]) |
     sum(basic:==basic[i])!=1) return(meta)
 }
 kr=(1::m);work=0
 if(k<m) {
  // Exact positive binary-power row/RHS duplicates only. Keep original columns.
  // Otherwise use a bounded exact original-row proof, including the RHS.
  meta[2]=3;kr=J(0,1,.)
  for(i=1;i<=m;i++) {
   if(all(original.A[i,1..n]:==0) & original.b[i]==0) continue
   hit=0
   for(j=1;j<=rows(kr);j++) {
    if(work+n+1>1000000) {meta[6]=work;return(meta);}
    work=work+n+1
    ratio=_rsm_duplicate_ratio((original.A[i,1..n],original.b[i]),
       (original.A[kr[j],1..n],original.b[kr[j]]))
    if(ratio<.) {hit=1;break;}
   }
   if(!hit) kr=kr\i
  }
  meta[6]=work
  if(rows(kr)!=k) {
   if(!_rsm_rowmap(original,kr,rowcols,rowdetail,1000000-work)) {meta[6]=meta[6]+rowdetail[4];return(meta);}
   meta[6]=meta[6]+rowdetail[4]
   if(rows(kr)!=k) return(meta)
  }
 }
 B=original.A[kr,basic]
 meta[2]=4
 check=_dea_lpcheck(candidate.x',original.A[,1..n],original.b,
    J(0,n,0),J(0,1,0),L',J(1,n,.))
 if(_rsm_success_fast==1) {
  _rsm_success_fastcounts[1]=_rsm_success_fastcounts[1]+1
  meta[7]=1;yf=lusolve(B',original.c[basic']')
  if(!hasmissing(yf)) {
   Y=J(m,1,0);Y[kr,.]=yf
   q=_rsm_quality_fast(original,candidate.x,Y',L,tol)
   if(!hasmissing(q)) if(q[1]==1 & _dea_lpcheck_ok(check,tol) & candidate.objval<.) {
    _rsm_success_fastcounts[2]=_rsm_success_fastcounts[2]+1
    meta[1]=1;meta[2]=0;meta[9]=1;meta[11]=q[4];meta[12]=q[5]
    return(meta)
   }
  }
 }
 _rsm_success_fastcounts[3]=_rsm_success_fastcounts[3]+1
 meta[7]=meta[7]+1
 YP=_rsm_retry_refine(B',original.c[basic']')
 if(hasmissing(YP)) return(meta)
 Y=J(m,2,0);Y[kr,.]=YP
 if(affine) q=_rsm_affine_quality(original,candidate.x,Y[,1]',Y[,2]',L,tol)
 else q=_rsm_retry_quality(original,candidate.x,Y[,1]',Y[,2]',tol)
 meta[2]=5
 if(hasmissing(q)) return(meta)
 meta[9]=q[1];meta[11]=q[4];meta[12]=q[5]
 check=_dea_lpcheck(candidate.x',original.A[,1..n],original.b,
    J(0,n,0),J(0,1,0),L',J(1,n,.))
 if(q[1]==1 & _dea_lpcheck_ok(check,tol) & candidate.objval<.) {
  meta[1]=1;meta[2]=0;return(meta)
 }
 // No alteration of a previously certified result. Refinement only after failure.
 meta[1]=0;meta[2]=6;meta[8]=1
 _rsm_success_fastcounts[4]=_rsm_success_fastcounts[4]+1
 if(affine) XP=_rsm_affine_refine(original.A[kr,1..n],original.b[kr],basic,L)
 else XP=_rsm_retry_refine(B,original.b[kr])
 if(hasmissing(XP)) return(meta)
 x=L
 for(i=1;i<=k;i++) x[basic[i]]=max((L[basic[i]],XP[i,1]+XP[i,2]))
 if(affine) q=_rsm_affine_quality(original,x,Y[,1]',Y[,2]',L,tol)
 else q=_rsm_retry_quality(original,x,Y[,1]',Y[,2]',tol)
 meta[2]=7
 if(hasmissing(q)) return(meta)
 meta[10]=q[1];meta[13]=q[4];meta[14]=q[5]
 check=_dea_lpcheck(x',original.A[,1..n],original.b,
    J(0,n,0),J(0,1,0),L',J(1,n,.))
 if(q[1]!=1 | !_dea_lpcheck_ok(check,tol) | original.c[1..n]*x>=.) return(meta)
 candidate.x=x;candidate.objval=original.c[1..n]*x
 candidate.msg="Optimal after original-basis success refinement"
 meta[1]=2;meta[2]=0
 return(meta)
}
void _rsm_success_apply(struct RSM_Prob scalar original,
    struct RSM_Result scalar candidate,real colvector basic,real scalar tol,| real colvector L)
{
 external real scalar _rsm_success_policy,_rsm_success_recording
 external real matrix _rsm_success_records
 external real rowvector _rsm_success_counts
 real rowvector meta
 if(_rsm_success_policy==0 | candidate.status!=0) return
 if(args()==5) meta=_rsm_success_check(original,candidate,basic,tol,L)
 else meta=_rsm_success_check(original,candidate,basic,tol)
 _rsm_success_counts[1]=_rsm_success_counts[1]+1
 if(meta[1]==1) _rsm_success_counts[2]=_rsm_success_counts[2]+1
 if(meta[1]==2) _rsm_success_counts[3]=_rsm_success_counts[3]+1
 if(meta[1]==0) _rsm_success_counts[4]=_rsm_success_counts[4]+1
 if(meta[1]==-1) _rsm_success_counts[5]=_rsm_success_counts[5]+1
 _rsm_success_counts[6]=_rsm_success_counts[6]+meta[7]
 _rsm_success_counts[7]=_rsm_success_counts[7]+meta[8]
 if(_rsm_success_recording==1) {
  if(rows(_rsm_success_records)<10000) _rsm_success_records=_rsm_success_records\meta
  else _rsm_success_counts[8]=_rsm_success_counts[8]+1
 }
 // Isolated comparison policy: 1 retains uncertified old successes; 2 rejects.
 if(_rsm_success_policy==2 & meta[1]<=0) {
  candidate.status=3;candidate.objval=.;candidate.x=J(original.n_orig,1,.)
  candidate.msg="Residual check: accepted-result certificate unavailable"
  if(meta[1]==0) candidate.msg="Residual check: accepted-result quality failed"
  candidate.recovery[2]=0;candidate.recovery[8]=3
 }
}

struct RSM_Result scalar _rsm_two_phase_core(struct RSM_Prob scalar prob,
    real scalar tol,string scalar trace,real scalar allow_duplicates,| real scalar check_success,real scalar allow_rows)
{
    struct RSM_Prob scalar original,retry
    struct RSM_Result scalar first,second,third
    real rowvector check,row
    real scalar certified,cycle_seen,pivot_seen,rebuild_seen
    real colvector certbasis,success_basic
    if(args()<5) check_success=1
    if(args()<6) allow_rows=1
    original=prob
    first=_rsm_eta_two_phase(prob,tol,trace,success_basic)
    row=(0,0,first.status,.,first.iters,0,first.iters,first.status,_rsm_recovery_reason(first.msg),.)
    if(first.status!=3) {
        first.recovery=row
        first.basiscols=success_basic
        if(check_success) _rsm_success_apply(original,first,success_basic,tol)

        return(first)
    }
    retry=original
    second=_rsm_fp_two_phase(retry,tol,trace,success_basic)
    row[1]=1;row[4]=second.status;row[6]=second.iters;row[7]=first.iters+second.iters;row[10]=_rsm_recovery_reason(second.msg)
    second.iters=row[7]
    if(second.status==0 & rows(second.x)==original.n_orig) {
        check=_dea_lpcheck(second.x',original.A[,1..original.n_orig],original.b,
            J(0,original.n_orig,0),J(0,1,0),J(1,original.n_orig,0),J(1,original.n_orig,.))
        second.objval=original.c[1..original.n_orig]*second.x
        if(_dea_lpcheck_ok(check,tol) & second.objval<.) {
            row[2]=1;row[8]=0;second.recovery=row
            second.msg="Optimal after factorized recovery";prob=original

            second.basiscols=success_basic
            if(check_success) _rsm_success_apply(original,second,success_basic,tol)
            return(second)
        }
        row[10]=9
    }
    // Preserve existing restart priority, including alternative optimal weights.
    // Certified pricing is attempted only after an observed repeated pivot state
    // and failure of every preceding acceptance screen.
    cycle_seen=0;rebuild_seen=0
    for(certified=0;certified<=1;certified++) {
        if(certified & !cycle_seen) break
        retry=original
        third=_rsm_precision_two_phase(retry,tol,trace,certbasis,certified,cycle_seen,0,pivot_seen)
        rebuild_seen=rebuild_seen | pivot_seen
        row[6]=row[6]+third.iters;row[7]=row[7]+third.iters
        if(_rsm_precision_accept(original,third,certbasis,tol)) {
            row[2]=1;row[4]=0;row[8]=0;row[10]=0
            third.basiscols=certbasis;third.iters=row[7];third.recovery=row
            if(certified) third.msg="Optimal after certified-price retry"
            else third.msg="Optimal after original-basis precision retry"
            prob=original;
            return(third)
        }
    }
    // Preserve every earlier success, including restart-selected weights.
    // Re-run from the untouched LP only after an observed Phase1 pivot stop.
    if(rebuild_seen) {
        retry=original
        third=_rsm_precision_two_phase(retry,tol,trace,certbasis,0,cycle_seen,1,pivot_seen)
        row[6]=row[6]+third.iters;row[7]=row[7]+third.iters
        if(_rsm_precision_accept(original,third,certbasis,tol)) {
            row[2]=1;row[4]=0;row[8]=0;row[10]=0
            third.basiscols=certbasis;third.iters=row[7];third.recovery=row
            third.msg="Optimal after original-column reconstruction"
            prob=original;
            return(third)
        }
    }
    if(allow_duplicates) {
        third=_rsm_duplicate_retry(original,tol,trace,check_success)
        row[6]=row[6]+third.iters;row[7]=row[7]+third.iters
        if(third.status==0) {
            row[2]=1;row[4]=0;row[8]=0;row[10]=0
            third.iters=row[7];third.recovery=row;prob=original
            return(third)
        }
    }
    if(allow_rows) {
        third=_rsm_row_retry(original,tol,trace,check_success)
        row[6]=row[6]+third.iters;row[7]=row[7]+third.iters
        if(third.status==0) {
            row[2]=1;row[4]=0;row[8]=0;row[10]=0
            third.iters=row[7];third.recovery=row;prob=original
            return(third)
        }
    }
    first.iters=row[7];first.recovery=row

    return(first)
}

struct RSM_Result scalar _rsm_two_phase(struct RSM_Prob scalar prob,
    real scalar tol,string scalar trace,| real scalar check_success)
{
    struct RSM_Result scalar result
    real scalar m,n
    if(args()<4) check_success=1
    m=prob.m_con;n=prob.n_orig
    result=_rsm_two_phase_core(prob,tol,trace,1,check_success)
    _rsm_audit_emit(result.recovery,m,n)
    return(result)
}

_rsm_audit_reset(0)
_rsm_success_reset(1,0)
_dea_link_reset(2,0)

end
