*! version 2.0.6  05aug2026
*! _dea_ipm.ado — IPM Engine for dea v2.0
*! Interior Point Method engine using LinearProgram() class (Stata 16+)
*!
*! AUTHORS:
*! --------
*!   Choonjoo Lee (corresponding author)
*!     Korea National Defense University
*!     Seoul, Korea
*!     bloom.rampike@gmail.com; bloom.rampike@gmail.com
*!
*!
*! NOTE: This IPM engine module is entirely new in v2.0 (no code from
*!   prior versions). It uses Stata 16+ LinearProgram() class with
*!   primal-dual fallback. Prior versions used ldea.do (Revised Simplex).
*!
*! MODEL REFERENCES (see dea.ado header for complete list):
*! --------------------------------------------------------
*!   [2]  Charnes, Cooper, & Rhodes (1978) — CCR (CRS)
*!   [3]  Banker, Charnes, & Cooper (1984) — BCC (VRS)
*!   [4]  Cooper, Seiford, & Tone (2007) — envelopment LP, cost/revenue/profit
*!   [6]  Andersen & Petersen (1993) — super-efficiency
*!   [7]  Chambers, Chung, & Färe (1996) — DDF
*!   [9]  Tone (2001) — SBM
*!   [10] Tone (2002) — Super-SBM
*!   [11] Deprins, Simar, & Tulkens (2006) — FDH
*!   [12] Cooper, Park, & Yu (1999) — IDEA
*!   [13] Färe, Grosskopf, & Lovell (1993) — RTS determination
*!   [15] Charnes & Cooper (1962) — Charnes-Cooper transformation
*! 
*! Changelog:
*!   2.0.6  05aug2026  Scale TE input/output columns before LinearProgram
*!                      stage 1/2 solves and restore slacks to source units.
*!   2.0.5  20jul2026  Align input multiplier u0 sign and DRS/IRS bounds with Cooper
*!   2.0.4  20jul2026  Return a common signed multiplier intercept
*!   2.0.3  20jul2026  Fix DRS/IRS stage-2 scalar min/max Mata syntax
*!   2.0.1  08feb2026  Fix: supereff output-oriented φ lower bound 1 → 1e-8
*!                      (φ≥1 prevented super-efficient DMU scores < 1)
*!   2.0.0  05feb2026  Initial packaged version; based on earlier teaching/seminar versions
*!
*! Features:
*!   - Unified envelopment: CRS, VRS, DRS, IRS
*!   - Stage 2 slack maximization
*!   - Primal → Dual (Multiplier) fallback
*!   - Lambda (reference weights) extraction
*!   - VRS Frontier analysis (SCALE, RTS)

/**
 * MATA FUNCTION INDEX:
 * ---------------------------------------------------------------------------
 * Core Solver:
 *   _dea_ipm_envelopment()     Unified envelopment: all RTS, all orientations
 *                               Returns: theta/phi, lambda, convergence info
 *   _dea_ipm_slacks()          Stage 2: slack maximization (theta fixed)
 *                               Returns: input slacks, output slacks
 *
 * Dual Fallback:
 *   _dea_ipm_dual_input()      Multiplier form, input-oriented
 *   _dea_ipm_dual_output()     Multiplier form, output-oriented
 *
 * VRS Frontier:
 *   _dea_ipm_vrsfrontier()     CRS→VRS→DRS + SCALE + RTS
 *
 * Utilities:
 *   _dea_ipm_count_results()   Count converged/failed DMUs
 * ---------------------------------------------------------------------------
 */


// ============================================================================
// Clear existing functions to allow clean reload
// ============================================================================
capture mata: mata drop _dea_ipm_envelopment()
capture mata: mata drop _dea_ipm_slacks()
capture mata: mata drop _dea_ipm_dual_input()
capture mata: mata drop _dea_ipm_dual_output()
capture mata: mata drop _dea_ipm_vrsfrontier()
capture mata: mata drop _dea_ipm_count_results()
capture mata: mata drop _dea_ipm_supereff()
capture mata: mata drop _dea_ipm_ddf()
capture mata: mata drop _dea_ipm_ddf_undesirable()
capture mata: mata drop _dea_ipm_additive()
capture mata: mata drop _dea_ipm_cost()
capture mata: mata drop _dea_ipm_revenue()
capture mata: mata drop _dea_ipm_profit()
capture mata: mata drop _dea_ipm_sbm()
capture mata: mata drop _dea_ipm_sbm_oriented()
capture mata: mata drop _dea_ipm_multiplier()
capture mata: mata drop _dea_ipm_supersbm()
capture mata: mata drop _dea_fdh()
capture mata: mata drop _dea_idea()

// Helper function for bootstrap CI computation
capture mata: mata drop _dea_bootstrap_ci()
mata:
void _dea_bootstrap_ci(string scalar matname, 
                       real scalar lo_pct, 
                       real scalar hi_pct,
                       string scalar lo_var,
                       string scalar hi_var)
{
    real matrix B
    real scalar n, reps, i, lo_idx, hi_idx
    real rowvector sorted_row
    real colvector ci_lo, ci_hi
    
    B = st_matrix(matname)
    n = rows(B)
    reps = cols(B)
    
    lo_idx = max((1, ceil(reps * lo_pct / 100)))
    hi_idx = min((reps, floor(reps * hi_pct / 100)))
    
    ci_lo = J(n, 1, .)
    ci_hi = J(n, 1, .)
    
    for (i = 1; i <= n; i++) {
        sorted_row = sort(B[i, .]', 1)'
        ci_lo[i] = sorted_row[lo_idx]
        ci_hi[i] = sorted_row[hi_idx]
    }
    
    st_store(., st_varindex(lo_var), ci_lo)
    st_store(., st_varindex(hi_var), ci_hi)
}
end

mata:
mata set matastrict on

// ============================================================================
// 1. UNIFIED ENVELOPMENT SOLVER
// ============================================================================
/**
 * _dea_ipm_envelopment()
 *
 * Solves the DEA envelopment (primal) form for all DMUs with dual fallback.
 *
 * Supports all four RTS assumptions:
 *   CRS: no convexity constraint
 *   VRS: sum(lambda) = 1
 *   DRS: sum(lambda) <= 1  (Non-Increasing RTS)
 *   IRS: sum(lambda) >= 1  (Non-Decreasing RTS)
 *
 * Input-oriented:  min theta  s.t. X*lambda <= theta*x0, Y*lambda >= y0
 * Output-oriented: max phi    s.t. X*lambda <= x0,       Y*lambda >= phi*y0
 *
 * Arguments:
 *   inputs       - string: input variable names
 *   outputs      - string: output variable names
 *   eff_var      - string: efficiency score variable name
 *   lambda_var   - string: lambda variable names (space-separated), "" if none
 *   slack_in_var - string: input slack variable names, "" if none
 *   slack_out_var- string: output slack variable names, "" if none
 *   touse        - string: estimation sample marker
 *   orientation  - string: "input" or "output"
 *   rts          - string: "crs", "vrs", "drs", "irs"
 *   do_stage2    - real:   1 = run Stage 2 slack maximization, 0 = skip
 *   tol          - real:   numerical tolerance
 */
void _dea_ipm_envelopment(string scalar inputs,
                           string scalar outputs,
                           string scalar eff_var,
                           string scalar lambda_var,
                           string scalar slack_in_var,
                           string scalar slack_out_var,
                           string scalar touse,
                           string scalar orientation,
                           string scalar rts,
                           real scalar do_stage2,
                           real scalar tol)
{
    // ----------------------------------------------------------------
    // Data loading
    // ----------------------------------------------------------------
    real matrix X, Y
    real rowvector x_scale, y_scale
    real colvector obs_indices
    real scalar n, m, s, k, obs_k, j
    
    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)

    // Positive column scaling leaves radial DEA scores and lambdas unchanged,
    // while protecting LinearProgram stage 1 and stage 2 from heterogeneous
    // financial-statement magnitudes.  Stored slacks are restored below.
    x_scale = J(1, m, 1)
    y_scale = J(1, s, 1)
    for (j = 1; j <= m; j++) {
        x_scale[j] = max(abs(X[., j]))
        if (x_scale[j] <= 0 | x_scale[j] >= .) x_scale[j] = 1
        X[., j] = X[., j] / x_scale[j]
    }
    for (j = 1; j <= s; j++) {
        y_scale[j] = max(abs(Y[., j]))
        if (y_scale[j] <= 0 | y_scale[j] >= .) y_scale[j] = 1
        Y[., j] = Y[., j] / y_scale[j]
    }
    
    // Lambda variable list (may be empty)
    string rowvector lam_varlist, si_varlist, so_varlist
    real scalar has_lambda, has_slacks
    
    has_lambda = (lambda_var != "")
    has_slacks = (slack_in_var != "" & slack_out_var != "")
    
    if (has_lambda) lam_varlist = tokens(lambda_var)
    if (has_slacks) {
        si_varlist = tokens(slack_in_var)
        so_varlist = tokens(slack_out_var)
    }
    
    // ----------------------------------------------------------------
    // Solve each DMU
    // ----------------------------------------------------------------
    real colvector x0, y0
    real rowvector c, lb, ub, sol, lam_k
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real scalar rc, eff_score, primal_ok
    real scalar phi_raw
    class LinearProgram scalar S
    
    // Stage 2 variables
    real rowvector slacks_in, slacks_out
    
    for (k = 1; k <= n; k++) {
        
        obs_k = obs_indices[k]
        x0 = X[k, ]'   // column vector (m × 1)
        y0 = Y[k, ]'   // column vector (s × 1)
        primal_ok = 0
        eff_score = .
        lam_k = J(1, n, .)
        slacks_in = J(1, m, 0)
        slacks_out = J(1, s, 0)
        
        // ============================================================
        // PRIMAL (ENVELOPMENT)
        // ============================================================
        if (orientation == "input") {
            // ---------------------------------------------------------
            // Input-oriented: min theta
            //   s.t. X'*lambda - theta*x0 <= 0   (m rows)
            //        -Y'*lambda           <= -y0  (s rows)
            //        [convexity constraint]
            //        lambda >= 0, theta >= 0
            // Variables: [lambda(1..n), theta]
            // ---------------------------------------------------------
            
            c = J(1, n, 0), 1    // min theta
            
            A_ineq = (X', -x0) \ (-Y', J(s, 1, 0))
            b_ineq = J(m, 1, 0) \ (-y0)
            
            lb = J(1, n+1, 0)
            ub = J(1, n+1, .)
            
            // Convexity constraint based on RTS
            A_eq = J(0, n+1, .)
            b_eq = J(0, 1, .)
            
            if (rts == "vrs") {
                // sum(lambda) = 1
                A_eq = (J(1, n, 1), 0)
                b_eq = 1
            }
            else if (rts == "drs") {
                // sum(lambda) <= 1 → add as inequality
                A_ineq = A_ineq \ (J(1, n, 1), 0)
                b_ineq = b_ineq \ 1
            }
            else if (rts == "irs") {
                // sum(lambda) >= 1 → -sum(lambda) <= -1
                A_ineq = A_ineq \ (J(1, n, -1), 0)
                b_ineq = b_ineq \ (-1)
            }
            // CRS: no convexity constraint
            
            // Solve
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setInequality(A_ineq, b_ineq)
            if (rows(A_eq) > 0) {
                S.setEquality(A_eq, b_eq)
            }
            S.setBounds(lb, ub)
            S.setTrace("off")
            
            rc = S.optimize()
            
            if (S.converged() && S.errorcode() == 0) {
                sol = S.parameters()
                eff_score = S.value()
                lam_k = sol[1..n]
                primal_ok = 1
            }
        }
        else {
            // ---------------------------------------------------------
            // Output-oriented: max phi  (= min -phi)
            //   s.t. X'*lambda           <= x0   (m rows)
            //        -Y'*lambda + phi*y0 <= 0    (s rows)
            //        [convexity constraint]
            //        lambda >= 0, phi >= 1
            // Variables: [lambda(1..n), phi]
            // ---------------------------------------------------------
            
            c = J(1, n, 0), -1   // min(-phi) = max(phi)
            
            A_ineq = (X', J(m, 1, 0)) \ (-Y', y0)
            b_ineq = x0 \ J(s, 1, 0)
            
            lb = J(1, n, 0), 1   // phi >= 1
            ub = J(1, n+1, .)
            
            // Convexity constraint
            A_eq = J(0, n+1, .)
            b_eq = J(0, 1, .)
            
            if (rts == "vrs") {
                A_eq = (J(1, n, 1), 0)
                b_eq = 1
            }
            else if (rts == "drs") {
                A_ineq = A_ineq \ (J(1, n, 1), 0)
                b_ineq = b_ineq \ 1
            }
            else if (rts == "irs") {
                A_ineq = A_ineq \ (J(1, n, -1), 0)
                b_ineq = b_ineq \ (-1)
            }
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setInequality(A_ineq, b_ineq)
            if (rows(A_eq) > 0) {
                S.setEquality(A_eq, b_eq)
            }
            S.setBounds(lb, ub)
            S.setTrace("off")
            
            rc = S.optimize()
            
            if (S.converged() && S.errorcode() == 0) {
                sol = S.parameters()
                phi_raw = -S.value()    // phi (>= 1)
                eff_score = 1 / phi_raw // efficiency in [0, 1]
                lam_k = sol[1..n]
                primal_ok = 1
            }
        }
        
        // ============================================================
        // DUAL FALLBACK (if primal failed)
        // ============================================================
        if (!primal_ok) {
            if (orientation == "input") {
                eff_score = _dea_ipm_dual_input(X, Y, k, rts)
            }
            else {
                eff_score = _dea_ipm_dual_output(X, Y, k, rts)
            }
            lam_k = J(1, n, .)  // lambda unavailable from dual
        }
        
        // ============================================================
        // STAGE 2: SLACK MAXIMIZATION (if requested and Stage 1 OK)
        // ============================================================
        if (do_stage2 && primal_ok && eff_score != .) {
            _dea_ipm_slacks(X, Y, k, orientation, rts,
                            eff_score, lam_k,
                            slacks_in, slacks_out, tol)
        }
        
        // ============================================================
        // STORE RESULTS
        // ============================================================
        st_store(obs_k, st_varindex(eff_var), eff_score)
        
        if (has_lambda) {
            for (j = 1; j <= n; j++) {
                st_store(obs_k, st_varindex(lam_varlist[j]), lam_k[j])
            }
        }
        
        if (has_slacks && do_stage2) {
            for (j = 1; j <= m; j++) {
                st_store(obs_k, st_varindex(si_varlist[j]), slacks_in[j] * x_scale[j])
            }
            for (j = 1; j <= s; j++) {
                st_store(obs_k, st_varindex(so_varlist[j]), slacks_out[j] * y_scale[j])
            }
        }
    }
    
    // Store convergence counts
    _dea_ipm_count_results(eff_var, touse, n)
}


// ============================================================================
// 2. STAGE 2: SLACK MAXIMIZATION
// ============================================================================
/**
 * _dea_ipm_slacks()
 *
 * Stage 2: Given theta* from Stage 1, maximize total slacks.
 *
 * Input-oriented:
 *   max sum(s_in) + sum(s_out)
 *   s.t. X*lambda + s_in  = theta* * x_k    (inputs)
 *        Y*lambda - s_out = y_k              (outputs)
 *        [convexity constraint]
 *        lambda, s_in, s_out >= 0
 *
 * Output-oriented:
 *   max sum(s_in) + sum(s_out)
 *   s.t. X*lambda + s_in  = x_k             (inputs)
 *        Y*lambda - s_out = phi* * y_k       (outputs)
 *        [convexity constraint]
 *        lambda, s_in, s_out >= 0
 *
 * Variables: [lambda(n), s_in(m), s_out(s)]
 *
 * Modifies slacks_in and slacks_out in place.
 */
void _dea_ipm_slacks(real matrix X,
                      real matrix Y,
                      real scalar k,
                      string scalar orientation,
                      string scalar rts,
                      real scalar theta_star,
                      real rowvector lam_k,
                      real rowvector slacks_in,
                      real rowvector slacks_out,
                      real scalar tol)
{
    real scalar n, m, s, total_vars, j
    real rowvector c, lb, ub, sol
    real matrix A_eq
    real colvector b_eq, x0, y0
    real scalar rc, sum_lam
    class LinearProgram scalar S2
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    total_vars = n + m + s  // lambda(n) + s_in(m) + s_out(s)
    
    x0 = X[k, ]'
    y0 = Y[k, ]'
    
    // Compute sum of lambda from Stage 1 for CRS constraint
    sum_lam = 0
    for (j = 1; j <= n; j++) {
        if (lam_k[j] < . && lam_k[j] > 0) sum_lam = sum_lam + lam_k[j]
    }
    if (sum_lam == 0) sum_lam = 1  // fallback
    
    // Objective: max sum(s_in) + sum(s_out)
    //         = min -sum(s_in) - sum(s_out)
    c = (J(1, n, 0), -J(1, m, 1), -J(1, s, 1))
    
    // Equality constraints
    if (orientation == "input") {
        // X'*lambda + I*s_in = theta* * x_k
        // Y'*lambda - I*s_out = y_k
        A_eq = (X', I(m), J(m, s, 0)) \
               (Y', J(s, m, 0), -I(s))
        b_eq = (theta_star * x0) \ y0
    }
    else {
        // X'*lambda + I*s_in = x_k
        // Y'*lambda - I*s_out = (1/theta_star) * y_k  [phi* = 1/theta_star]
        real scalar phi_star
        phi_star = 1 / theta_star
        
        A_eq = (X', I(m), J(m, s, 0)) \
               (Y', J(s, m, 0), -I(s))
        b_eq = x0 \ (phi_star * y0)
    }
    
    // Convexity constraint based on RTS
    if (rts == "crs") {
        // For CRS, fix sum(lambda) to Stage 1 value to ensure feasibility
        A_eq = A_eq \ (J(1, n, 1), J(1, m + s, 0))
        b_eq = b_eq \ sum_lam
    }
    else if (rts == "vrs") {
        A_eq = A_eq \ (J(1, n, 1), J(1, m + s, 0))
        b_eq = b_eq \ 1
    }
    else if (rts == "drs") {
        // Use sum(lambda) = min(1, sum_lam) for DRS
        A_eq = A_eq \ (J(1, n, 1), J(1, m + s, 0))
        b_eq = b_eq \ min((1, sum_lam))
    }
    else if (rts == "irs") {
        // Use sum(lambda) = max(1, sum_lam) for IRS
        A_eq = A_eq \ (J(1, n, 1), J(1, m + s, 0))
        b_eq = b_eq \ max((1, sum_lam))
    }
    
    // Bounds: all >= 0
    lb = J(1, total_vars, 0)
    ub = J(1, total_vars, .)
    
    // Solve
    S2 = LinearProgram()
    S2.setMaxOrMin("min")
    S2.setCoefficients(c)
    S2.setEquality(A_eq, b_eq)
    S2.setBounds(lb, ub)
    S2.setTrace("off")
    
    rc = S2.optimize()
    
    if (S2.converged() && S2.errorcode() == 0) {
        sol = S2.parameters()
        // Extract slacks
        slacks_in  = sol[(n+1)..(n+m)]
        slacks_out = sol[(n+m+1)..(n+m+s)]
        
        // Clean tiny values
        for (j = 1; j <= m; j++) {
            if (abs(slacks_in[j]) < tol) slacks_in[j] = 0
        }
        for (j = 1; j <= s; j++) {
            if (abs(slacks_out[j]) < tol) slacks_out[j] = 0
        }
    }
    else {
        // Stage 2 failed — leave slacks at 0
        slacks_in  = J(1, m, 0)
        slacks_out = J(1, s, 0)
    }
}


// ============================================================================
// 3. DUAL FALLBACK SOLVERS
// ============================================================================

/**
 * _dea_ipm_dual_input()
 *
 * Multiplier form, input-oriented:
 *   max u'y0 [+ u0]  s.t. v'x0 = 1, u'Y - v'X [+ u0] <= 0
 *   theta = u'y0 [+ u0]
 *
 * Supports CRS/VRS/DRS/IRS:
 *   CRS: no u0 term
 *   VRS: u0 free (unrestricted)
 *   DRS: u0 <= 0  (non-positive)
 *   IRS: u0 >= 0  (non-negative)
 */
real scalar _dea_ipm_dual_input(real matrix X,
                                real matrix Y,
                                real scalar k,
                                string scalar rts)
{
    real scalar n, m_in, s_out, total_vars, eff
    real rowvector x0, y0, c, lb, ub
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    class LinearProgram scalar lp
    
    n = rows(X)
    m_in = cols(X)
    s_out = cols(Y)
    x0 = X[k, .]
    y0 = Y[k, .]
    
    lp = LinearProgram()
    lp.setTrace("off")
    
    if (rts == "crs") {
        // Variables: [u(s), v(m)]
        total_vars = s_out + m_in
        c = (-y0, J(1, m_in, 0))        // min -u'y0 = max u'y0
        A_eq = (J(1, s_out, 0), x0)     // v'x0 = 1
        A_ineq = (Y, -X)                // u'Y - v'X <= 0
        lb = J(1, total_vars, 0)
    }
    else {
        // Variables: [u(s), v(m), u0]
        total_vars = s_out + m_in + 1
        c = (-y0, J(1, m_in, 0), -1)    // min -(u'y0 + u0)
        A_eq = (J(1, s_out, 0), x0, 0)  // v'x0 = 1
        A_ineq = (Y, -X, J(n, 1, 1))   // u'Y - v'X + u0 <= 0
        
        // Bounds on u0 depend on RTS
        if (rts == "vrs") {
            lb = (J(1, s_out + m_in, 0), .)  // u0 free
        }
        else if (rts == "drs") {
            lb = (J(1, s_out + m_in, 0), .)  // u0 free (but u0 <= 0)
            // DRS dual: u0 <= 0, implemented via upper bound
        }
        else if (rts == "irs") {
            lb = (J(1, s_out + m_in, 0), 0)  // u0 >= 0
        }
    }
    b_eq = 1
    b_ineq = J(n, 1, 0)
    ub = J(1, total_vars, .)
    
    // DRS: u0 <= 0
    if (rts == "drs") {
        ub[total_vars] = 0
    }
    
    lp.setCoefficients(c)
    lp.setMaxOrMin("min")
    lp.setEquality(A_eq, b_eq)
    lp.setInequality(A_ineq, b_ineq)
    lp.setBounds(lb, ub)
    
    (void) lp.optimize()
    
    if (lp.converged() == 1 & lp.errorcode() == 0) {
        eff = -lp.value()   // theta = u'y0 [+ u0]
    }
    else {
        eff = .
    }
    
    return(eff)
}


/**
 * _dea_ipm_dual_output()
 *
 * Multiplier form, output-oriented:
 *   min v'x0 [- u0]  s.t. u'y0 = 1, u'Y - v'X [+ u0] <= 0
 *   1/phi = v'x0 [- u0]   →  efficiency = v'x0 [- u0]
 *
 * Returns efficiency in [0,1] scale.
 */
real scalar _dea_ipm_dual_output(real matrix X,
                                  real matrix Y,
                                  real scalar k,
                                  string scalar rts)
{
    real scalar n, m_in, s_out, total_vars, eff
    real rowvector x0, y0, c, lb, ub
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    class LinearProgram scalar lp
    
    n = rows(X)
    m_in = cols(X)
    s_out = cols(Y)
    x0 = X[k, .]
    y0 = Y[k, .]
    
    lp = LinearProgram()
    lp.setTrace("off")
    
    if (rts == "crs") {
        total_vars = s_out + m_in
        c = (J(1, s_out, 0), x0)        // min v'x0
        A_eq = (y0, J(1, m_in, 0))      // u'y0 = 1
        A_ineq = (Y, -X)                // u'Y - v'X <= 0
        lb = J(1, total_vars, 0)
    }
    else {
        total_vars = s_out + m_in + 1
        c = (J(1, s_out, 0), x0, -1)    // min v'x0 - u0
        A_eq = (y0, J(1, m_in, 0), 0)   // u'y0 = 1
        A_ineq = (Y, -X, J(n, 1, 1))   // u'Y - v'X + u0 <= 0
        
        if (rts == "vrs") {
            lb = (J(1, s_out + m_in, 0), .)  // u0 free
        }
        else if (rts == "drs") {
            lb = (J(1, s_out + m_in, 0), .)  // u0 free, capped at 0
        }
        else if (rts == "irs") {
            lb = (J(1, s_out + m_in, 0), 0)  // u0 >= 0
        }
    }
    b_eq = 1
    b_ineq = J(n, 1, 0)
    ub = J(1, total_vars, .)
    
    if (rts == "drs") {
        ub[total_vars] = 0  // u0 <= 0
    }
    
    lp.setCoefficients(c)
    lp.setMaxOrMin("min")
    lp.setEquality(A_eq, b_eq)
    lp.setInequality(A_ineq, b_ineq)
    lp.setBounds(lb, ub)
    
    (void) lp.optimize()
    
    if (lp.converged() == 1 & lp.errorcode() == 0) {
        eff = lp.value()   // = 1/phi (already in [0,1])
        if (eff <= 0) eff = .
    }
    else {
        eff = .
    }
    
    return(eff)
}


// ============================================================================
// 4. VRS FRONTIER ANALYSIS
// ============================================================================
/**
 * _dea_ipm_vrsfrontier()
 *
 * Runs CRS → VRS → DRS sequentially for VRS Frontier analysis.
 * Computes:
 *   CRS_TE  = CRS efficiency
 *   VRS_TE  = VRS efficiency
 *   SCALE   = CRS_TE / VRS_TE
 *   RTS     = based on DRS_TE comparison:
 *             if CRS_TE == VRS_TE → CRS (0)
 *             else if DRS_TE == VRS_TE → DRS (-1)
 *             else → IRS (1)
 *
 * Stores results in 5 Stata variables: crs_te, vrs_te, drs_te, scale, rts_type
 *
 * Alternative RTS method (sum of lambda):
 *   sum(lambda_CRS) < 1 → IRS
 *   sum(lambda_CRS) > 1 → DRS
 *   sum(lambda_CRS) = 1 → CRS
 */
void _dea_ipm_vrsfrontier(string scalar inputs,
                           string scalar outputs,
                           string scalar crs_var,
                           string scalar vrs_var,
                           string scalar drs_var,
                           string scalar scale_var,
                           string scalar rts_var,
                           string scalar lambda_sum_var,
                           string scalar touse,
                           string scalar orientation,
                           real scalar tol)
{
    real matrix X, Y
    real colvector obs_indices
    real scalar n, m, s, k, obs_k
    real scalar crs_eff, vrs_eff, drs_eff, scale_eff, rts_type
    real scalar sum_lambda
    
    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    // Solve each DMU for three RTS assumptions
    for (k = 1; k <= n; k++) {
        
        obs_k = obs_indices[k]
        
        real colvector x0, y0
        real rowvector c, lb, ub, sol, lam_crs
        real matrix A_ineq, A_eq
        real colvector b_ineq, b_eq
        real scalar rc, phi_raw
        class LinearProgram scalar S
        
        x0 = X[k, ]'
        y0 = Y[k, ]'
        
        crs_eff = .
        vrs_eff = .
        drs_eff = .
        sum_lambda = .
        lam_crs = J(1, n, .)
        
        // ==================================================
        // (A) CRS — with lambda extraction for RTS method 2
        // ==================================================
        if (orientation == "input") {
            c = J(1, n, 0), 1
            A_ineq = (X', -x0) \ (-Y', J(s, 1, 0))
            b_ineq = J(m, 1, 0) \ (-y0)
            lb = J(1, n+1, 0)
            ub = J(1, n+1, .)
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setInequality(A_ineq, b_ineq)
            S.setBounds(lb, ub)
            S.setTrace("off")
            rc = S.optimize()
            
            if (S.converged() && S.errorcode() == 0) {
                sol = S.parameters()
                crs_eff = S.value()
                lam_crs = sol[1..n]
                sum_lambda = sum(lam_crs)
            }
            else {
                crs_eff = _dea_ipm_dual_input(X, Y, k, "crs")
            }
        }
        else {
            c = J(1, n, 0), -1
            A_ineq = (X', J(m, 1, 0)) \ (-Y', y0)
            b_ineq = x0 \ J(s, 1, 0)
            lb = J(1, n, 0), 1
            ub = J(1, n+1, .)
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setInequality(A_ineq, b_ineq)
            S.setBounds(lb, ub)
            S.setTrace("off")
            rc = S.optimize()
            
            if (S.converged() && S.errorcode() == 0) {
                sol = S.parameters()
                phi_raw = -S.value()
                crs_eff = 1 / phi_raw
                lam_crs = sol[1..n]
                sum_lambda = sum(lam_crs)
            }
            else {
                crs_eff = _dea_ipm_dual_output(X, Y, k, "crs")
            }
        }
        
        // ==================================================
        // (B) VRS
        // ==================================================
        if (orientation == "input") {
            c = J(1, n, 0), 1
            A_ineq = (X', -x0) \ (-Y', J(s, 1, 0))
            b_ineq = J(m, 1, 0) \ (-y0)
            A_eq = (J(1, n, 1), 0)
            b_eq = 1
            lb = J(1, n+1, 0)
            ub = J(1, n+1, .)
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setInequality(A_ineq, b_ineq)
            S.setEquality(A_eq, b_eq)
            S.setBounds(lb, ub)
            S.setTrace("off")
            rc = S.optimize()
            
            if (S.converged() && S.errorcode() == 0) {
                vrs_eff = S.value()
            }
            else {
                vrs_eff = _dea_ipm_dual_input(X, Y, k, "vrs")
            }
        }
        else {
            c = J(1, n, 0), -1
            A_ineq = (X', J(m, 1, 0)) \ (-Y', y0)
            b_ineq = x0 \ J(s, 1, 0)
            A_eq = (J(1, n, 1), 0)
            b_eq = 1
            lb = J(1, n, 0), 1
            ub = J(1, n+1, .)
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setInequality(A_ineq, b_ineq)
            S.setEquality(A_eq, b_eq)
            S.setBounds(lb, ub)
            S.setTrace("off")
            rc = S.optimize()
            
            if (S.converged() && S.errorcode() == 0) {
                phi_raw = -S.value()
                vrs_eff = 1 / phi_raw
            }
            else {
                vrs_eff = _dea_ipm_dual_output(X, Y, k, "vrs")
            }
        }
        
        // ==================================================
        // (C) DRS (for Färe-Grosskopf-Lovell RTS determination)
        // ==================================================
        if (orientation == "input") {
            c = J(1, n, 0), 1
            A_ineq = (X', -x0) \ (-Y', J(s, 1, 0)) \ (J(1, n, 1), 0)
            b_ineq = J(m, 1, 0) \ (-y0) \ 1
            lb = J(1, n+1, 0)
            ub = J(1, n+1, .)
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setInequality(A_ineq, b_ineq)
            S.setBounds(lb, ub)
            S.setTrace("off")
            rc = S.optimize()
            
            if (S.converged() && S.errorcode() == 0) {
                drs_eff = S.value()
            }
            else {
                drs_eff = _dea_ipm_dual_input(X, Y, k, "drs")
            }
        }
        else {
            c = J(1, n, 0), -1
            A_ineq = (X', J(m, 1, 0)) \ (-Y', y0) \ (J(1, n, 1), 0)
            b_ineq = x0 \ J(s, 1, 0) \ 1
            lb = J(1, n, 0), 1
            ub = J(1, n+1, .)
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setInequality(A_ineq, b_ineq)
            S.setBounds(lb, ub)
            S.setTrace("off")
            rc = S.optimize()
            
            if (S.converged() && S.errorcode() == 0) {
                phi_raw = -S.value()
                drs_eff = 1 / phi_raw
            }
            else {
                drs_eff = _dea_ipm_dual_output(X, Y, k, "drs")
            }
        }
        
        // ==================================================
        // (D) Compute SCALE and RTS
        // ==================================================
        
        // Scale efficiency = CRS_TE / VRS_TE
        if (crs_eff != . & vrs_eff != . & vrs_eff > 0) {
            scale_eff = crs_eff / vrs_eff
        }
        else {
            scale_eff = .
        }
        
        // RTS determination (Färe-Grosskopf-Lovell method):
        //   if CRS_TE == VRS_TE → CRS (0)
        //   else if DRS_TE == VRS_TE → DRS (-1)
        //   else → IRS (1)
        if (crs_eff != . & vrs_eff != .) {
            if (abs(crs_eff - vrs_eff) < tol) {
                rts_type = 0   // CRS
            }
            else if (drs_eff != . & abs(drs_eff - vrs_eff) < tol) {
                rts_type = -1  // DRS
            }
            else {
                rts_type = 1   // IRS
            }
        }
        else {
            rts_type = .
        }
        
        // Store results
        st_store(obs_k, st_varindex(crs_var), crs_eff)
        st_store(obs_k, st_varindex(vrs_var), vrs_eff)
        st_store(obs_k, st_varindex(drs_var), drs_eff)
        st_store(obs_k, st_varindex(scale_var), scale_eff)
        st_store(obs_k, st_varindex(rts_var), rts_type)
        if (lambda_sum_var != "") {
            st_store(obs_k, st_varindex(lambda_sum_var), sum_lambda)
        }
    }
}


// ============================================================================
// 5. UTILITY: COUNT RESULTS
// ============================================================================
void _dea_ipm_count_results(string scalar eff_var,
                             string scalar touse,
                             real scalar n)
{
    real colvector eff_data
    real scalar n_conv
    
    eff_data = st_data(., eff_var, touse)
    n_conv = rows(eff_data) - colmissing(eff_data)
    st_numscalar("r(n_converged)", n_conv)
    st_numscalar("r(n_failed)", n - n_conv)
}


// ============================================================================
// 6. SUPER-EFFICIENCY MODEL
// ============================================================================
/**
 * _dea_ipm_supereff()
 *
 * Computes super-efficiency by excluding the evaluated DMU from its reference set.
 * Input-oriented:  min theta  s.t. X_{-o}*lambda <= theta*x0, Y_{-o}*lambda >= y0
 *   theta > 1 → super-efficient (stored as theta)
 * Output-oriented: max phi    s.t. X_{-o}*lambda <= x0,       Y_{-o}*lambda >= phi*y0
 *   phi < 1 → super-efficient (stored as 1/phi > 1)
 *
 * Arguments:
 *   lambda_var   - string: lambda variable names (n, space-separated), "" if none
 *                  LP solves over n-1 lambdas (excluding DMU k);
 *                  results are mapped back to n positions with lambda_k = 0.
 */
void _dea_ipm_supereff(string scalar inputs,
                        string scalar outputs,
                        string scalar eff_var,
                        string scalar lambda_var,
                        string scalar touse,
                        string scalar orientation,
                        string scalar rts)
{
    real matrix X, Y, X_ref, Y_ref
    real colvector obs_indices, x0, y0
    real rowvector c, lb, ub, sol, lam_ref, lam_k
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real scalar n, m, s, k, obs_k, rc, eff_score, n_ref, j
    real scalar phi_raw, has_lambda
    class LinearProgram scalar S
    
    // Lambda variable handling
    string rowvector lam_varlist
    has_lambda = (lambda_var != "")
    if (has_lambda) lam_varlist = tokens(lambda_var)
    
    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, ]'
        y0 = Y[k, ]'
        eff_score = .
        lam_k = J(1, n, 0)  // Initialize: all zeros (k-th is always 0)
        
        // Create reference set excluding DMU k
        n_ref = n - 1
        if (k == 1) {
            X_ref = X[2..n, ]
            Y_ref = Y[2..n, ]
        }
        else if (k == n) {
            X_ref = X[1..n-1, ]
            Y_ref = Y[1..n-1, ]
        }
        else {
            X_ref = X[1..k-1, ] \ X[k+1..n, ]
            Y_ref = Y[1..k-1, ] \ Y[k+1..n, ]
        }
        
        if (orientation == "input") {
            // min theta s.t. X_ref'*lambda <= theta*x0, Y_ref'*lambda >= y0
            // Variables: [lambda(1..n_ref), theta]
            c = J(1, n_ref, 0), 1
            A_ineq = (X_ref', -x0) \ (-Y_ref', J(s, 1, 0))
            b_ineq = J(m, 1, 0) \ (-y0)
            lb = J(1, n_ref+1, 0)
            ub = J(1, n_ref+1, .)
            
            A_eq = J(0, n_ref+1, .)
            b_eq = J(0, 1, .)
            
            if (rts == "vrs") {
                A_eq = (J(1, n_ref, 1), 0)
                b_eq = 1
            }
            else if (rts == "drs") {
                A_ineq = A_ineq \ (J(1, n_ref, 1), 0)
                b_ineq = b_ineq \ 1
            }
            else if (rts == "irs") {
                A_ineq = A_ineq \ (J(1, n_ref, -1), 0)
                b_ineq = b_ineq \ (-1)
            }
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setInequality(A_ineq, b_ineq)
            if (rows(A_eq) > 0) {
                S.setEquality(A_eq, b_eq)
            }
            S.setBounds(lb, ub)
            S.setTrace("off")
            
            rc = S.optimize()
            
            if (S.converged() && S.errorcode() == 0) {
                sol = S.parameters()
                eff_score = S.value()
                lam_ref = sol[1..n_ref]
                
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
            // Output-oriented: max phi
            // Super-efficiency allows phi < 1 (super-efficient DMUs)
            // Variables: [lambda(1..n_ref), phi]
            c = J(1, n_ref, 0), -1
            A_ineq = (X_ref', J(m, 1, 0)) \ (-Y_ref', y0)
            b_ineq = x0 \ J(s, 1, 0)
            lb = J(1, n_ref, 0), 1e-8
            ub = J(1, n_ref+1, .)
            
            A_eq = J(0, n_ref+1, .)
            b_eq = J(0, 1, .)
            
            if (rts == "vrs") {
                A_eq = (J(1, n_ref, 1), 0)
                b_eq = 1
            }
            else if (rts == "drs") {
                A_ineq = A_ineq \ (J(1, n_ref, 1), 0)
                b_ineq = b_ineq \ 1
            }
            else if (rts == "irs") {
                A_ineq = A_ineq \ (J(1, n_ref, -1), 0)
                b_ineq = b_ineq \ (-1)
            }
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setInequality(A_ineq, b_ineq)
            if (rows(A_eq) > 0) {
                S.setEquality(A_eq, b_eq)
            }
            S.setBounds(lb, ub)
            S.setTrace("off")
            
            rc = S.optimize()
            
            if (S.converged() && S.errorcode() == 0) {
                sol = S.parameters()
                phi_raw = -S.value()
                eff_score = 1 / phi_raw
                lam_ref = sol[1..n_ref]
                
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
                st_store(obs_k, st_varindex(lam_varlist[j]), lam_k[j])
            }
        }
    }
}


// ============================================================================
// 7. DIRECTIONAL DISTANCE FUNCTION (DDF)
// ============================================================================
/**
 * _dea_ipm_ddf()
 *
 * Solves the Directional Distance Function model.
 *
 * max beta s.t. X*lambda <= x0 - beta*gx
 *               Y*lambda >= y0 + beta*gy
 *               [convexity constraint]
 *               lambda >= 0, beta >= 0
 */
void _dea_ipm_ddf(string scalar inputs,
                  string scalar outputs,
                  string scalar eff_var,
                  string scalar gx_var,
                  string scalar gy_var,
                  string scalar touse,
                  string scalar rts)
{
    real matrix X, Y, GX, GY
    real colvector obs_indices, x0, y0, gx0, gy0
    real rowvector c, lb, ub, sol
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real scalar n, m, s, k, obs_k, rc, beta
    class LinearProgram scalar S
    
    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    // Load direction vectors (or use observed values as default)
    if (gx_var != "") {
        GX = st_data(., tokens(gx_var), touse)
    }
    else {
        GX = X  // Default: direction = observed input values
    }
    
    if (gy_var != "") {
        GY = st_data(., tokens(gy_var), touse)
    }
    else {
        GY = Y  // Default: direction = observed output values
    }
    
    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, ]'
        y0 = Y[k, ]'
        gx0 = GX[k, ]'
        gy0 = GY[k, ]'
        beta = .
        
        // max beta => min -beta
        // Variables: [lambda(1..n), beta]
        // Constraints:
        //   X'*lambda + beta*gx0 <= x0  =>  X'*lambda + beta*gx0 - x0 <= 0
        //   Y'*lambda - beta*gy0 >= y0  => -Y'*lambda + beta*gy0 <= -y0
        
        c = J(1, n, 0), -1   // min(-beta) = max(beta)
        
        A_ineq = (X', gx0) \ (-Y', gy0)
        b_ineq = x0 \ (-y0)
        
        lb = J(1, n+1, 0)
        ub = J(1, n+1, .)
        
        A_eq = J(0, n+1, .)
        b_eq = J(0, 1, .)
        
        if (rts == "vrs") {
            A_eq = (J(1, n, 1), 0)
            b_eq = 1
        }
        else if (rts == "drs") {
            A_ineq = A_ineq \ (J(1, n, 1), 0)
            b_ineq = b_ineq \ 1
        }
        else if (rts == "irs") {
            A_ineq = A_ineq \ (J(1, n, -1), 0)
            b_ineq = b_ineq \ (-1)
        }
        
        S = LinearProgram()
        S.setMaxOrMin("min")
        S.setCoefficients(c)
        S.setInequality(A_ineq, b_ineq)
        if (rows(A_eq) > 0) {
            S.setEquality(A_eq, b_eq)
        }
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
// 7b. DDF WITH UNDESIRABLE (BAD) OUTPUTS — Wang et al. (2022)
// ============================================================================
/**
 * _dea_ipm_ddf_undesirable()
 *
 * Directional Distance Function with undesirable outputs (weak disposability).
 * Based on Chung et al. (1997), Färe et al. (2007), Wang et al. (2022).
 *
 * max β
 * s.t.  X'λ + β·gx ≤ x₀     (inputs, if gx ≠ 0)
 *       Y'λ + β·gy ≥ y₀     (good outputs: (-Y',-gy)·[λ;β] ≤ -y₀)
 *       B'λ + β·gb = b₀     (bad outputs: EQUALITY — weak disposability)
 *       [Σλ = 1  VRS]
 *       λ ≥ 0, β ≥ 0
 *
 * Direction vectors (direction()):
 *   "yb"  (default): gx=0,  gy=y₀, gb=b₀  — adjust good/bad outputs only
 *   "xyb":           gx=x₀, gy=y₀, gb=b₀  — adjust all three
 *   "x":             gx=x₀, gy=0,  gb=0   — input-oriented (bad outputs fixed)
 *
 * Good outputs expand in direction gy:
 *   Y'λ - β·gy ≥ y₀, equivalently Y'λ ≥ y₀ + β·gy.
 */
void _dea_ipm_ddf_undesirable(string scalar inputs,
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
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real scalar n, m, s, q, k, obs_k, rc, beta
    class LinearProgram scalar S

    X = st_data(., tokens(inputs),      touse)
    Y = st_data(., tokens(outputs),     touse)
    B = st_data(., tokens(bad_outputs), touse)
    obs_indices = selectindex(st_data(., touse))

    n = rows(X)
    m = cols(X)    // number of inputs
    s = cols(Y)    // number of good outputs
    q = cols(B)    // number of bad outputs

    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0  = X[k, .]'
        y0  = Y[k, .]'
        b0  = B[k, .]'
        beta = .

        // --- Direction vectors ---
        if (direction == "xyb") {
            gx0 = x0
            gy0 = y0
            gb0 = b0
        }
        else if (direction == "x") {
            gx0 = x0
            gy0 = J(s, 1, 0)
            gb0 = J(q, 1, 0)
        }
        else {  // "yb" (default)
            gx0 = J(m, 1, 0)
            gy0 = y0
            gb0 = b0
        }

        // Variables: [λ(n), β];  objective: min -β (= max β)
        c  = (J(1, n, 0), -1)
        lb = J(1, n+1, 0)
        ub = J(1, n+1, .)

        // Inequality constraints
        // X'λ + β·gx ≤ x₀ :  (X', gx0) · [λ;β] ≤ x₀
        // Y'λ - β·gy ≥ y₀ :  (-Y', -gy0) · [λ;β] ≤ -y₀
        A_ineq = (X', gx0) \ (-Y', -gy0)
        b_ineq = x0 \ (-y0)

        // Equality constraints: bad outputs + (optional) convexity
        // B'λ + β·gb = b₀ :  (B', gb0) · [λ;β] = b₀
        A_eq = (B', gb0)
        b_eq = b0

        if (rts == "vrs") {
            A_eq = A_eq \ (J(1, n, 1), 0)
            b_eq = b_eq \ 1
        }
        else if (rts == "drs") {
            A_ineq = A_ineq \ (J(1, n, 1), 0)
            b_ineq = b_ineq \ 1
        }
        else if (rts == "irs") {
            A_ineq = A_ineq \ (J(1, n, -1), 0)
            b_ineq = b_ineq \ (-1)
        }

        S = LinearProgram()
        S.setMaxOrMin("min")
        S.setCoefficients(c)
        S.setInequality(A_ineq, b_ineq)
        S.setEquality(A_eq, b_eq)
        S.setBounds(lb, ub)
        S.setTrace("off")

        rc = S.optimize()

        if (S.converged() && S.errorcode() == 0) {
            sol  = S.parameters()
            beta = sol[n+1]   // β is the last variable
        }

        st_store(obs_k, st_varindex(eff_var), beta)
    }
}


// ============================================================================
// 8. ADDITIVE MODEL
// ============================================================================
/**
 * _dea_ipm_additive()
 *
 * Maximizes the sum of input and output slacks.
 *
 * max sum(s-) + sum(s+)
 * s.t. X*lambda + s- = x0
 *      Y*lambda - s+ = y0
 *      [convexity constraint]
 *      lambda, s-, s+ >= 0
 */
void _dea_ipm_additive(string scalar inputs,
                        string scalar outputs,
                        string scalar eff_var,
                        string scalar slack_in_var,
                        string scalar slack_out_var,
                        string scalar touse,
                        string scalar rts)
{
    real matrix X, Y
    real colvector obs_indices, x0, y0
    real rowvector c, lb, ub, sol
    real matrix A_eq
    real colvector b_eq
    real matrix A_ineq
    real colvector b_ineq
    real scalar n, m, s, k, obs_k, rc, add_score
    real rowvector slack_in, slack_out
    class LinearProgram scalar S
    
    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    string rowvector si_varlist, so_varlist
    real scalar has_slacks, j
    
    has_slacks = (slack_in_var != "" & slack_out_var != "")
    if (has_slacks) {
        si_varlist = tokens(slack_in_var)
        so_varlist = tokens(slack_out_var)
    }
    
    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, ]'
        y0 = Y[k, ]'
        add_score = .
        slack_in = J(1, m, 0)
        slack_out = J(1, s, 0)
        
        // Variables: [lambda(n), s-(m), s+(s)]
        // Total: n + m + s variables
        // Objective: max sum(s-) + sum(s+) => min -sum(s-) - sum(s+)
        
        c = J(1, n, 0), J(1, m, -1), J(1, s, -1)
        
        // Equality constraints:
        // X'*lambda + I*s- = x0  => X'*lambda + I*s- - x0 = 0
        // Y'*lambda - I*s+ = y0  => Y'*lambda - I*s+ - y0 = 0
        
        A_eq = (X', I(m), J(m, s, 0)) \ (Y', J(s, m, 0), -I(s))
        b_eq = x0 \ y0
        
        // RTS constraints
        A_ineq = J(0, n+m+s, .)
        b_ineq = J(0, 1, .)
        
        if (rts == "vrs") {
            A_eq = A_eq \ (J(1, n, 1), J(1, m+s, 0))
            b_eq = b_eq \ 1
        }
        else if (rts == "drs") {
            A_ineq = (J(1, n, 1), J(1, m+s, 0))
            b_ineq = 1
        }
        else if (rts == "irs") {
            A_ineq = (J(1, n, -1), J(1, m+s, 0))
            b_ineq = -1
        }
        
        lb = J(1, n+m+s, 0)
        ub = J(1, n+m+s, .)
        
        S = LinearProgram()
        S.setMaxOrMin("min")
        S.setCoefficients(c)
        S.setEquality(A_eq, b_eq)
        if (rows(A_ineq) > 0) {
            S.setInequality(A_ineq, b_ineq)
        }
        S.setBounds(lb, ub)
        S.setTrace("off")
        
        rc = S.optimize()
        
        if (S.converged() && S.errorcode() == 0) {
            sol = S.parameters()
            slack_in = sol[n+1..n+m]
            slack_out = sol[n+m+1..n+m+s]
            add_score = sum(slack_in) + sum(slack_out)
        }
        
        // Store efficiency (0 if efficient, >0 if inefficient)
        st_store(obs_k, st_varindex(eff_var), add_score)
        
        // Store slacks if requested
        if (has_slacks) {
            for (j = 1; j <= m; j++) {
                st_store(obs_k, st_varindex(si_varlist[j]), slack_in[j])
            }
            for (j = 1; j <= s; j++) {
                st_store(obs_k, st_varindex(so_varlist[j]), slack_out[j])
            }
        }
    }
}


// ============================================================================
// 9. COST EFFICIENCY MODEL
// ============================================================================
/**
 * _dea_ipm_cost()
 *
 * Computes cost efficiency given input prices.
 *
 * min w'*x
 * s.t. X*lambda <= x
 *      Y*lambda >= y0
 *      [convexity constraint]
 *      lambda, x >= 0
 *
 * Cost Efficiency = min_cost / observed_cost = (w'*x*) / (w'*x0)
 */
void _dea_ipm_cost(string scalar inputs,
                   string scalar outputs,
                   string scalar eff_var,
                   string scalar w_var,
                   string scalar touse,
                   string scalar rts)
{
    real matrix X, Y, W
    real colvector obs_indices, x0, y0, w0
    real rowvector c, lb, ub, sol
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real scalar n, m, s, k, obs_k, rc
    real scalar min_cost, obs_cost, cost_eff
    class LinearProgram scalar S
    
    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    W = st_data(., tokens(w_var), touse)
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, ]'
        y0 = Y[k, ]'
        w0 = W[k, ]'
        cost_eff = .
        
        // Observed cost
        obs_cost = sum(w0 :* x0)
        
        // Variables: [lambda(n), x(m)]
        // Objective: min w'*x
        // Constraints:
        //   X'*lambda - x <= 0
        //   -Y'*lambda <= -y0
        //   [convexity]
        
        c = J(1, n, 0), w0'
        
        A_ineq = (X', -I(m)) \ (-Y', J(s, m, 0))
        b_ineq = J(m, 1, 0) \ (-y0)
        
        lb = J(1, n+m, 0)
        ub = J(1, n+m, .)
        
        A_eq = J(0, n+m, .)
        b_eq = J(0, 1, .)
        
        if (rts == "vrs") {
            A_eq = (J(1, n, 1), J(1, m, 0))
            b_eq = 1
        }
        else if (rts == "drs") {
            A_ineq = A_ineq \ (J(1, n, 1), J(1, m, 0))
            b_ineq = b_ineq \ 1
        }
        else if (rts == "irs") {
            A_ineq = A_ineq \ (J(1, n, -1), J(1, m, 0))
            b_ineq = b_ineq \ (-1)
        }
        
        S = LinearProgram()
        S.setMaxOrMin("min")
        S.setCoefficients(c)
        S.setInequality(A_ineq, b_ineq)
        if (rows(A_eq) > 0) {
            S.setEquality(A_eq, b_eq)
        }
        S.setBounds(lb, ub)
        S.setTrace("off")
        
        rc = S.optimize()
        
        if (S.converged() && S.errorcode() == 0) {
            min_cost = S.value()
            if (obs_cost > 0) {
                cost_eff = min_cost / obs_cost
            }
        }
        
        st_store(obs_k, st_varindex(eff_var), cost_eff)
    }
}


// ============================================================================
// 10. REVENUE EFFICIENCY MODEL
// ============================================================================
/**
 * _dea_ipm_revenue()
 *
 * Computes revenue efficiency given output prices.
 *
 * max p'*y
 * s.t. X*lambda <= x0
 *      Y*lambda >= y
 *      [convexity constraint]
 *      lambda, y >= 0
 *
 * Revenue Efficiency = observed_revenue / max_revenue = (p'*y0) / (p'*y*)
 */
void _dea_ipm_revenue(string scalar inputs,
                      string scalar outputs,
                      string scalar eff_var,
                      string scalar p_var,
                      string scalar touse,
                      string scalar rts)
{
    real matrix X, Y, P
    real colvector obs_indices, x0, y0, p0
    real rowvector c, lb, ub, sol
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real scalar n, m, s, k, obs_k, rc
    real scalar max_rev, obs_rev, rev_eff
    class LinearProgram scalar S
    
    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    P = st_data(., tokens(p_var), touse)
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, ]'
        y0 = Y[k, ]'
        p0 = P[k, ]'
        rev_eff = .
        
        // Observed revenue
        obs_rev = sum(p0 :* y0)
        
        // Variables: [lambda(n), y(s)]
        // Objective: max p'*y => min -p'*y
        // Constraints:
        //   X'*lambda <= x0
        //   Y'*lambda - y >= 0 => -Y'*lambda + y <= 0
        //   [convexity]
        
        c = J(1, n, 0), -p0'
        
        A_ineq = (X', J(m, s, 0)) \ (-Y', I(s))
        b_ineq = x0 \ J(s, 1, 0)
        
        lb = J(1, n+s, 0)
        ub = J(1, n+s, .)
        
        A_eq = J(0, n+s, .)
        b_eq = J(0, 1, .)
        
        if (rts == "vrs") {
            A_eq = (J(1, n, 1), J(1, s, 0))
            b_eq = 1
        }
        else if (rts == "drs") {
            A_ineq = A_ineq \ (J(1, n, 1), J(1, s, 0))
            b_ineq = b_ineq \ 1
        }
        else if (rts == "irs") {
            A_ineq = A_ineq \ (J(1, n, -1), J(1, s, 0))
            b_ineq = b_ineq \ (-1)
        }
        
        S = LinearProgram()
        S.setMaxOrMin("min")
        S.setCoefficients(c)
        S.setInequality(A_ineq, b_ineq)
        if (rows(A_eq) > 0) {
            S.setEquality(A_eq, b_eq)
        }
        S.setBounds(lb, ub)
        S.setTrace("off")
        
        rc = S.optimize()
        
        if (S.converged() && S.errorcode() == 0) {
            max_rev = -S.value()
            if (max_rev > 0) {
                rev_eff = obs_rev / max_rev
            }
        }
        
        st_store(obs_k, st_varindex(eff_var), rev_eff)
    }
}


// ============================================================================
// 10a. PROFIT EFFICIENCY MODEL
// ============================================================================
/**
 * _dea_ipm_profit()
 *
 * Computes profit efficiency (Cooper, Seiford & Tone, Ch.8).
 *
 * [Profit] eq. 8.8:
 *   py* - cx* = max  p'(Yλ) - c'(Xλ)
 *               λ
 *   s.t. x = Xλ ≤ x_o       (inputs bounded by observed)
 *        y = Yλ ≥ y_o        (outputs meet observed)
 *        [convexity]
 *        λ ≥ 0
 *
 * Variables: λ(n) only — x and y are determined by λ.
 * The LP is bounded under CRS because Xλ ≤ x_o constrains λ.
 *
 * Profit Efficiency = observed_profit / max_profit  (eq. 8.9)
 *   where observed_profit = p'y_o - c'x_o
 *         max_profit = p'(Yλ*) - c'(Xλ*)
 */
void _dea_ipm_profit(string scalar inputs,
                     string scalar outputs,
                     string scalar eff_var,
                     string scalar w_var,
                     string scalar p_var,
                     string scalar touse,
                     string scalar rts)
{
    real matrix X, Y, W, P
    real colvector obs_indices, x0, y0, w0, p0
    real rowvector c, lb, ub
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real scalar n, m, s, k, obs_k, rc
    real scalar max_profit, obs_profit, profit_eff
    class LinearProgram scalar S
    
    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    W = st_data(., tokens(w_var), touse)
    P = st_data(., tokens(p_var), touse)
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, ]'
        y0 = Y[k, ]'
        w0 = W[k, ]'
        p0 = P[k, ]'
        profit_eff = .
        
        // Observed profit
        obs_profit = sum(p0 :* y0) - sum(w0 :* x0)
        
        // Variables: λ(n) only
        // Objective: max p0'(Y'λ) - w0'(X'λ) = (Yp0 - Xw0)'λ
        //         => min (Xw0 - Yp0)'λ
        c = (X * w0 - Y * p0)'
        
        // Constraints:
        //   X'λ ≤ x0
        //   -Y'λ ≤ -y0
        A_ineq = X' \ (-Y')
        b_ineq = x0 \ (-y0)
        
        lb = J(1, n, 0)
        ub = J(1, n, .)
        
        A_eq = J(0, n, .)
        b_eq = J(0, 1, .)
        
        if (rts == "vrs") {
            A_eq = J(1, n, 1)
            b_eq = 1
        }
        else if (rts == "drs") {
            A_ineq = A_ineq \ J(1, n, 1)
            b_ineq = b_ineq \ 1
        }
        else if (rts == "irs") {
            A_ineq = A_ineq \ J(1, n, -1)
            b_ineq = b_ineq \ (-1)
        }
        
        S = LinearProgram()
        S.setMaxOrMin("min")
        S.setCoefficients(c)
        S.setInequality(A_ineq, b_ineq)
        if (rows(A_eq) > 0) {
            S.setEquality(A_eq, b_eq)
        }
        S.setBounds(lb, ub)
        S.setTrace("off")
        
        rc = S.optimize()
        
        if (S.converged() && S.errorcode() == 0) {
            max_profit = -S.value()
            if (max_profit > 0) {
                profit_eff = obs_profit / max_profit
            }
        }
        
        st_store(obs_k, st_varindex(eff_var), profit_eff)
    }
}


// ============================================================================
// 11. SBM (SLACKS-BASED MEASURE) MODEL
// ============================================================================
/**
 * _dea_ipm_sbm()
 *
 * Computes Slacks-Based Measure efficiency (Tone, 2001).
 * Uses Charnes-Cooper transformation to convert fractional program to LP.
 *
 * Original fractional program:
 *   min rho = [1 - (1/m)*sum(s-_i/x_io)] / [1 + (1/s)*sum(s+_r/y_ro)]
 *   s.t. x_o = X*lambda + s-
 *        y_o = Y*lambda - s+
 *        [convexity]
 *        lambda, s-, s+ >= 0
 *
 * After Charnes-Cooper transformation (t = 1/denominator):
 *   min t - (1/m)*sum(S-_i/x_io)
 *   s.t. t + (1/s)*sum(S+_r/y_ro) = 1
 *        t*x_o = X*Lambda + S-
 *        t*y_o = Y*Lambda - S+
 *        [convexity: sum(Lambda) = t]
 *        Lambda, S-, S+, t >= 0
 */
void _dea_ipm_sbm(string scalar inputs,
                  string scalar outputs,
                  string scalar eff_var,
                  string scalar slack_in_var,
                  string scalar slack_out_var,
                  string scalar touse,
                  string scalar rts)
{
    real matrix X, Y
    real colvector obs_indices, x0, y0
    real rowvector c, lb, ub, sol
    real matrix A_eq
    real colvector b_eq
    real matrix A_ineq
    real colvector b_ineq
    real scalar n, m, s, k, obs_k, rc
    real scalar sbm_eff, t_val
    real rowvector slack_in, slack_out, coef_sin, coef_sout
    class LinearProgram scalar S
    real scalar j
    
    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    string rowvector si_varlist, so_varlist
    real scalar has_slacks
    
    has_slacks = (slack_in_var != "" & slack_out_var != "")
    if (has_slacks) {
        si_varlist = tokens(slack_in_var)
        so_varlist = tokens(slack_out_var)
    }
    
    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, ]'
        y0 = Y[k, ]'
        sbm_eff = .
        slack_in = J(1, m, 0)
        slack_out = J(1, s, 0)
        
        // Variables: [Lambda(n), S-(m), S+(s), t]
        // Total: n + m + s + 1 variables
        
        // Coefficients for input slacks: -(1/m)/x_io
        coef_sin = J(1, m, 0)
        for (j = 1; j <= m; j++) {
            if (x0[j] > 0) {
                coef_sin[j] = -1 / (m * x0[j])
            }
        }
        
        // Coefficients for output slacks in normalizing constraint: (1/s)/y_ro
        coef_sout = J(1, s, 0)
        for (j = 1; j <= s; j++) {
            if (y0[j] > 0) {
                coef_sout[j] = 1 / (s * y0[j])
            }
        }
        
        // Objective: min t + sum(coef_sin * S-)
        c = J(1, n, 0), coef_sin, J(1, s, 0), 1
        
        // Equality constraints:
        // (1) Normalizing: t + sum(coef_sout * S+) = 1
        // (2) Input:  X'*Lambda + I*S- - t*x0 = 0
        // (3) Output: Y'*Lambda - I*S+ - t*y0 = 0
        
        // (1) Normalizing constraint
        A_eq = (J(1, n, 0), J(1, m, 0), coef_sout, 1)
        b_eq = 1
        
        // (2) Input constraints: X'*Lambda + I*S- = t*x0
        //     => X'*Lambda + I*S- - x0*t = 0
        A_eq = A_eq \ (X', I(m), J(m, s, 0), -x0)
        b_eq = b_eq \ J(m, 1, 0)
        
        // (3) Output constraints: Y'*Lambda - I*S+ = t*y0
        //     => Y'*Lambda - I*S+ - y0*t = 0
        A_eq = A_eq \ (Y', J(s, m, 0), -I(s), -y0)
        b_eq = b_eq \ J(s, 1, 0)
        
        // RTS constraint
        A_ineq = J(0, n+m+s+1, .)
        b_ineq = J(0, 1, .)
        
        if (rts == "vrs") {
            // sum(Lambda) = t => sum(Lambda) - t = 0
            A_eq = A_eq \ (J(1, n, 1), J(1, m+s, 0), -1)
            b_eq = b_eq \ 0
        }
        else if (rts == "drs") {
            // sum(Lambda) <= t => sum(Lambda) - t <= 0
            A_ineq = (J(1, n, 1), J(1, m+s, 0), -1)
            b_ineq = 0
        }
        else if (rts == "irs") {
            // sum(Lambda) >= t => -sum(Lambda) + t <= 0
            A_ineq = (J(1, n, -1), J(1, m+s, 0), 1)
            b_ineq = 0
        }
        // CRS: sum(Lambda) is free (no constraint added)
        
        lb = J(1, n+m+s+1, 0)
        ub = J(1, n+m+s+1, .)
        
        S = LinearProgram()
        S.setMaxOrMin("min")
        S.setCoefficients(c)
        S.setEquality(A_eq, b_eq)
        if (rows(A_ineq) > 0) {
            S.setInequality(A_ineq, b_ineq)
        }
        S.setBounds(lb, ub)
        S.setTrace("off")
        
        rc = S.optimize()
        
        if (S.converged() && S.errorcode() == 0) {
            sol = S.parameters()
            t_val = sol[n+m+s+1]
            if (t_val > 0) {
                // Convert back: s- = S-/t, s+ = S+/t
                slack_in = sol[n+1..n+m] / t_val
                slack_out = sol[n+m+1..n+m+s] / t_val
                sbm_eff = S.value()
            }
        }
        
        st_store(obs_k, st_varindex(eff_var), sbm_eff)
        
        if (has_slacks) {
            for (j = 1; j <= m; j++) {
                st_store(obs_k, st_varindex(si_varlist[j]), slack_in[j])
            }
            for (j = 1; j <= s; j++) {
                st_store(obs_k, st_varindex(so_varlist[j]), slack_out[j])
            }
        }
    }
}


// ============================================================================
// 12. MULTIPLIER (DUAL) MODEL
// ============================================================================
/**
 * _dea_ipm_multiplier()
 *
 * Solves the multiplier (dual) form of DEA.
 *
 * Input-oriented (max virtual output ratio):
 *   max u'*y0 - u0  (for VRS)
 *   s.t. v'*x0 = 1
 *        u'*Y - v'*X - u0 <= 0  (for each DMU)
 *        u, v >= epsilon
 *
 * Output-oriented (min virtual input ratio):
 *   min v'*x0 + v0  (for VRS)
 *   s.t. u'*y0 = 1
 *        u'*Y - v'*X - v0 >= 0
 *        u, v >= epsilon
 */
void _dea_ipm_multiplier(string scalar inputs,
                          string scalar outputs,
                          string scalar eff_var,
                          string scalar v_var,
                          string scalar u_var,
                          string scalar u0_var,
                          string scalar touse,
                          string scalar orientation,
                          string scalar rts)
{
    real matrix X, Y
    real colvector obs_indices, x0, y0
    real rowvector c, lb, ub, sol
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real scalar n, m, s, k, obs_k, rc, eff_score, u0_value
    real scalar epsilon, has_u0
    real rowvector v_weights, u_weights
    class LinearProgram scalar S
    
    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    epsilon = 0     // Standard multiplier form permits zero weights
    has_u0 = (rts == "vrs" | rts == "drs" | rts == "irs")
    
    string rowvector v_varlist, u_varlist
    real scalar has_weights, j
    
    has_weights = (v_var != "" & u_var != "")
    if (has_weights) {
        v_varlist = tokens(v_var)
        u_varlist = tokens(u_var)
    }
    
    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, ]'
        y0 = Y[k, ]'
        eff_score = .
        v_weights = J(1, m, 0)
        u_weights = J(1, s, 0)
        u0_value = 0
        
        if (orientation == "input") {
            // Variables: [v(m), u(s), u0] if VRS, else [v(m), u(s)]
            // Objective: max u'*y0 - u0 => min -u'*y0 + u0
            
            if (has_u0) {
                c = J(1, m, 0), -y0', 1
                
                // Normalizing: v'*x0 = 1
                A_eq = (x0', J(1, s, 0), 0)
                b_eq = 1
                
                // Constraints: u'*y_j - v'*x_j - u0 <= 0 for all j
                A_ineq = (-X, Y, J(n, 1, -1))
                b_ineq = J(n, 1, 0)
                
                lb = J(1, m, epsilon), J(1, s, epsilon), .
                ub = J(1, m+s+1, .)
                
                // Cooper signed intercept: DRS u0 >= 0; IRS u0 <= 0
                if (rts == "drs") {
                    lb[m+s+1] = 0
                }
                else if (rts == "irs") {
                    ub[m+s+1] = 0
                }
            }
            else {
                // CRS
                c = J(1, m, 0), -y0'
                
                A_eq = (x0', J(1, s, 0))
                b_eq = 1
                
                A_ineq = (-X, Y)
                b_ineq = J(n, 1, 0)
                
                lb = J(1, m, epsilon), J(1, s, epsilon)
                ub = J(1, m+s, .)
            }
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setEquality(A_eq, b_eq)
            S.setInequality(A_ineq, b_ineq)
            S.setBounds(lb, ub)
            S.setTrace("off")
            
            rc = S.optimize()
            
            if (S.converged() && S.errorcode() == 0) {
                sol = S.parameters()
                eff_score = -S.value()
                v_weights = sol[1..m]
                u_weights = sol[m+1..m+s]
                if (has_u0) u0_value = sol[m+s+1]
            }
        }
        else {
            // Output-oriented
            // Variables: [v(m), u(s), v0] if VRS
            // Objective: min v'*x0 + v0
            
            if (has_u0) {
                c = x0', J(1, s, 0), 1
                
                // Normalizing: u'*y0 = 1
                A_eq = (J(1, m, 0), y0', 0)
                b_eq = 1
                
                // Constraints: v'*x_j - u'*y_j + v0 >= 0 => -v'*x_j + u'*y_j - v0 <= 0
                A_ineq = (-X, Y, J(n, 1, -1))
                b_ineq = J(n, 1, 0)
                
                lb = J(1, m, epsilon), J(1, s, epsilon), .
                ub = J(1, m+s+1, .)
                
                if (rts == "drs") {
                    lb[m+s+1] = 0
                }
                else if (rts == "irs") {
                    ub[m+s+1] = 0
                }
            }
            else {
                c = x0', J(1, s, 0)
                
                A_eq = (J(1, m, 0), y0')
                b_eq = 1
                
                A_ineq = (-X, Y)
                b_ineq = J(n, 1, 0)
                
                lb = J(1, m, epsilon), J(1, s, epsilon)
                ub = J(1, m+s, .)
            }
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setEquality(A_eq, b_eq)
            S.setInequality(A_ineq, b_ineq)
            S.setBounds(lb, ub)
            S.setTrace("off")
            
            rc = S.optimize()
            
            if (S.converged() && S.errorcode() == 0) {
                sol = S.parameters()
                eff_score = 1 / S.value()  // Convert to standard efficiency
                v_weights = sol[1..m]
                u_weights = sol[m+1..m+s]
                if (has_u0) u0_value = sol[m+s+1]
            }
        }
        
        st_store(obs_k, st_varindex(eff_var), eff_score)
        
        if (has_weights) {
            for (j = 1; j <= m; j++) {
                st_store(obs_k, st_varindex(v_varlist[j]), v_weights[j])
            }
            for (j = 1; j <= s; j++) {
                st_store(obs_k, st_varindex(u_varlist[j]), u_weights[j])
            }
            if (u0_var != "") {
                st_store(obs_k, st_varindex(u0_var), u0_value)
            }
        }
    }
}


// ============================================================================
// 13. SUPER-SBM MODEL (Tone, 2002)
// ============================================================================
/**
 * _dea_ipm_supersbm()
 *
 * Computes Super-efficiency SBM by excluding the evaluated DMU from reference set.
 * For efficient DMUs: rho > 1
 * For inefficient DMUs: returns standard SBM score
 *
 * Super-SBM for efficient unit o:
 *   min rho = [1 + (1/m)*sum((x_bar_i - x_io)/x_io)] / [1 - (1/s)*sum((y_bar_r - y_ro)/y_ro)]
 *   s.t. x_bar >= X_{-o}*lambda
 *        y_bar <= Y_{-o}*lambda
 *        y_bar <= y_o
 *        lambda >= 0
 */
void _dea_ipm_supersbm(string scalar inputs,
                        string scalar outputs,
                        string scalar eff_var,
                        string scalar touse,
                        string scalar rts)
{
    real matrix X, Y, X_ref, Y_ref
    real colvector obs_indices, x0, y0
    real rowvector c, lb, ub, coef_xbar, coef_ybar
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real scalar n, m, s, k, obs_k, rc, n_ref, j
    real scalar sbm_eff, sbm_score, tol
    class LinearProgram scalar S
    
    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    tol = 1e-6
    
    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, ]'
        y0 = Y[k, ]'
        sbm_eff = .
        
        // ==============================================================
        // Step 1: Standard SBM to identify efficient DMUs
        //   Charnes-Cooper linearization of Tone (2001):
        //   Variables: [lambda_hat(n), s_minus(m), s_plus(s), t]
        //   min  t - (1/m) * sum(s_minus_i / x_i0)
        //   s.t. t + (1/s) * sum(s_plus_r / y_r0) = 1
        //        X' * lambda_hat + s_minus = t * x_0
        //        Y' * lambda_hat - s_plus  = t * y_0
        //        VRS: sum(lambda_hat) = t
        //        lambda_hat, s_minus, s_plus >= 0, t > 0
        // ==============================================================
        {
            real rowvector c_sbm, lb_sbm, ub_sbm, coef_sm, coef_sp
            real matrix A_eq_sbm, A_ineq_sbm
            real colvector b_eq_sbm, b_ineq_sbm
            real scalar nv_sbm
            class LinearProgram scalar S_sbm
            
            nv_sbm = n + m + s + 1
            
            // Coefficients for s_minus in objective: -(1/m)/x_{i0}
            coef_sm = J(1, m, 0)
            for (j = 1; j <= m; j++) {
                if (x0[j] > 0) coef_sm[j] = -1 / (m * x0[j])
            }
            // Coefficients for s_plus in normalizing: (1/s)/y_{r0}
            coef_sp = J(1, s, 0)
            for (j = 1; j <= s; j++) {
                if (y0[j] > 0) coef_sp[j] = 1 / (s * y0[j])
            }
            
            // Objective: min t - (1/m) * sum(s_minus_i / x_i0)
            c_sbm = J(1, n, 0), coef_sm, J(1, s, 0), 1
            
            // Eq 1: t + (1/s) * sum(s_plus_r / y_r0) = 1
            A_eq_sbm = (J(1, n, 0), J(1, m, 0), coef_sp, 1)
            b_eq_sbm = 1
            
            // Eq 2: X'*lambda_hat + s_minus - t*x_0 = 0
            A_eq_sbm = A_eq_sbm \ (X', I(m), J(m, s, 0), -x0)
            b_eq_sbm = b_eq_sbm \ J(m, 1, 0)
            
            // Eq 3: Y'*lambda_hat - s_plus - t*y_0 = 0
            A_eq_sbm = A_eq_sbm \ (Y', J(s, m, 0), -I(s), -y0)
            b_eq_sbm = b_eq_sbm \ J(s, 1, 0)
            
            // RTS
            if (rts == "vrs") {
                A_eq_sbm = A_eq_sbm \ (J(1, n, 1), J(1, m+s, 0), -1)
                b_eq_sbm = b_eq_sbm \ 0
            }
            
            lb_sbm = J(1, nv_sbm, 0)
            ub_sbm = J(1, nv_sbm, .)
            
            S_sbm = LinearProgram()
            S_sbm.setMaxOrMin("min")
            S_sbm.setCoefficients(c_sbm)
            S_sbm.setEquality(A_eq_sbm, b_eq_sbm)
            
            if (rts == "drs") {
                A_ineq_sbm = (J(1, n, 1), J(1, m+s, 0), -1)
                b_ineq_sbm = 0
                S_sbm.setInequality(A_ineq_sbm, b_ineq_sbm)
            }
            else if (rts == "irs") {
                A_ineq_sbm = (J(1, n, -1), J(1, m+s, 0), 1)
                b_ineq_sbm = 0
                S_sbm.setInequality(A_ineq_sbm, b_ineq_sbm)
            }
            
            S_sbm.setBounds(lb_sbm, ub_sbm)
            S_sbm.setTrace("off")
            
            rc = S_sbm.optimize()
            
            if (S_sbm.converged() && S_sbm.errorcode() == 0) {
                sbm_score = S_sbm.value()
            }
            else {
                sbm_score = .
            }
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
            X_ref = X[2..n, ]
            Y_ref = Y[2..n, ]
        }
        else if (k == n) {
            X_ref = X[1..n-1, ]
            Y_ref = Y[1..n-1, ]
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
        A_ineq = (X_ref', -I(m), J(m, s, 0), J(m, 1, 0))
        b_ineq = J(m, 1, 0)
        
        // Inequality 2: y_tilde <= Y_ref' * Lambda
        A_ineq = A_ineq \ (-Y_ref', J(s, m, 0), I(s), J(s, 1, 0))
        b_ineq = b_ineq \ J(s, 1, 0)
        
        // Inequality 3: x_tilde >= t * x_0  [FIX: was missing]
        A_ineq = A_ineq \ (J(m, n_ref, 0), -I(m), J(m, s, 0), x0)
        b_ineq = b_ineq \ J(m, 1, 0)
        
        // Inequality 4: y_tilde <= t * y_0
        A_ineq = A_ineq \ (J(s, n_ref, 0), J(s, m, 0), I(s), -y0)
        b_ineq = b_ineq \ J(s, 1, 0)
        
        // RTS constraints
        if (rts == "vrs") {
            A_eq = A_eq \ (J(1, n_ref, 1), J(1, m+s, 0), -1)
            b_eq = b_eq \ 0
        }
        else if (rts == "drs") {
            A_ineq = A_ineq \ (J(1, n_ref, 1), J(1, m+s, 0), -1)
            b_ineq = b_ineq \ 0
        }
        else if (rts == "irs") {
            A_ineq = A_ineq \ (J(1, n_ref, -1), J(1, m+s, 0), 1)
            b_ineq = b_ineq \ 0
        }
        
        lb = J(1, n_ref + m + s + 1, 0)
        lb[n_ref + m + s + 1] = 1e-10   // t > 0
        ub = J(1, n_ref + m + s + 1, .)
        
        S = LinearProgram()
        S.setMaxOrMin("min")
        S.setCoefficients(c)
        S.setEquality(A_eq, b_eq)
        S.setInequality(A_ineq, b_ineq)
        S.setBounds(lb, ub)
        S.setTrace("off")
        
        rc = S.optimize()
        
        if (S.converged() && S.errorcode() == 0) {
            sbm_eff = S.value()
        }
        
        st_store(obs_k, st_varindex(eff_var), sbm_eff)
    }
}


// ============================================================================
// 14. FDH (FREE DISPOSAL HULL) - No LP Required
// ============================================================================
/**
 * _dea_fdh()
 *
 * Free Disposal Hull efficiency using dominance comparison.
 * No convexity assumption - only actual observed DMUs as references.
 *
 * Input-oriented:
 *   theta_k = min_j max_i (x_ji / x_ki), over references with y_j >= y_k.
 *   This is the radial Farrell input score relative to the non-convex FDH.
 *
 * Output-oriented:
 *   phi_k = max_j min_r (y_jr / y_kr), over references with x_j <= x_k.
 *   Efficiency is stored as 1/phi, matching dea.ado output-score convention.
 */
void _dea_fdh(string scalar inputs,
              string scalar outputs,
              string scalar eff_var,
              string scalar touse,
              string scalar orientation)
{
    real matrix X, Y
    real colvector obs_indices
    real scalar n, m, s, k, j, i, r, obs_k
    real scalar theta, phi, ratio, min_ratio, max_ratio
    real scalar feasible_ref
    real colvector x_k, y_k, x_j, y_j
    
    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x_k = X[k, ]'
        y_k = Y[k, ]'
        
        if (orientation == "input") {
            // Input-oriented radial FDH:
            // among output-feasible references, minimize the maximum input ratio.
            theta = 1  // default: efficient
            
            for (j = 1; j <= n; j++) {
                x_j = X[j, ]'
                y_j = Y[j, ]'
                
                feasible_ref = 1
                for (r = 1; r <= s; r++) {
                    if (y_j[r] < y_k[r]) {
                        feasible_ref = 0
                        break
                    }
                }
                
                if (feasible_ref) {
                    max_ratio = 0
                    for (i = 1; i <= m; i++) {
                        if (x_k[i] > 0) {
                            ratio = x_j[i] / x_k[i]
                            if (ratio > max_ratio) {
                                max_ratio = ratio
                            }
                        }
                    }
                    if (max_ratio < theta) {
                        theta = max_ratio
                    }
                }
            }
            
            st_store(obs_k, st_varindex(eff_var), theta)
        }
        else {
            // Output-oriented radial FDH:
            // among input-feasible references, maximize the minimum output ratio.
            phi = 1  // default: efficient (phi=1 means eff=1)
            
            for (j = 1; j <= n; j++) {
                x_j = X[j, ]'
                y_j = Y[j, ]'
                
                feasible_ref = 1
                for (i = 1; i <= m; i++) {
                    if (x_j[i] > x_k[i]) {
                        feasible_ref = 0
                        break
                    }
                }
                
                if (feasible_ref) {
                    min_ratio = .
                    for (r = 1; r <= s; r++) {
                        if (y_k[r] > 0) {
                            ratio = y_j[r] / y_k[r]
                            if (min_ratio == . || ratio < min_ratio) {
                                min_ratio = ratio
                            }
                        }
                    }
                    if (min_ratio != . && min_ratio > phi) {
                        phi = min_ratio
                    }
                }
            }
            
            // Output efficiency = 1/phi (or just store phi for expansion factor)
            st_store(obs_k, st_varindex(eff_var), 1/phi)
        }
    }
}


end  // end of mata block


// ============================================================================
// 15. IDEA (IMPRECISE DEA) - Interval/Bounded Data
// ============================================================================
/**
 * IDEA computes efficiency under data uncertainty.
 * 
 * For interval data [xlo, xhi] and [ylo, yhi]:
 *   - Pessimistic: DMU k at worst (high inputs, low outputs)
 *                  Reference set at best (low inputs, high outputs)
 *   - Optimistic:  DMU k at best (low inputs, high outputs)
 *                  Reference set at worst (high inputs, low outputs)
 *
 * Uses IPM solver for LP.
 */

mata:

void _dea_idea(string scalar inputs,
               string scalar outputs,
               string scalar xlo_vars,
               string scalar xhi_vars,
               string scalar ylo_vars,
               string scalar yhi_vars,
               string scalar eff_pess_var,
               string scalar eff_opt_var,
               string scalar touse,
               string scalar orientation,
               string scalar rts,
               real scalar tol)
{
    real matrix X, Y, Xlo, Xhi, Ylo, Yhi
    real matrix X_pess, Y_pess, X_opt, Y_opt
    real colvector obs_indices, x0, y0
    real rowvector c, lb, ub, sol, lam_k
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real scalar n, m, s, k, j, obs_k, rc
    real scalar eff_pess, eff_opt
    class LinearProgram scalar S
    
    // Load base data
    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    // Load interval bounds (if provided)
    if (xlo_vars != "") {
        Xlo = st_data(., tokens(xlo_vars), touse)
    } else {
        Xlo = X
    }
    if (xhi_vars != "") {
        Xhi = st_data(., tokens(xhi_vars), touse)
    } else {
        Xhi = X
    }
    if (ylo_vars != "") {
        Ylo = st_data(., tokens(ylo_vars), touse)
    } else {
        Ylo = Y
    }
    if (yhi_vars != "") {
        Yhi = st_data(., tokens(yhi_vars), touse)
    } else {
        Yhi = Y
    }
    
    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        eff_pess = .
        eff_opt = .
        
        // ================================================================
        // PESSIMISTIC EFFICIENCY
        // DMU k: worst case (Xhi, Ylo)
        // Reference: best case (Xlo, Yhi)
        // ================================================================
        if (orientation == "input") {
            // min theta s.t. Xlo' * lambda - theta * Xhi[k] <= 0
            //                -Yhi' * lambda <= -Ylo[k]
            x0 = Xhi[k, ]'
            y0 = Ylo[k, ]'
            
            c = J(1, n, 0), 1
            A_ineq = (Xlo', -x0) \ (-Yhi', J(s, 1, 0))
            b_ineq = J(m, 1, 0) \ (-y0)
            
            lb = J(1, n+1, 0)
            ub = J(1, n+1, .)
            
            A_eq = J(0, n+1, .)
            b_eq = J(0, 1, .)
            if (rts == "vrs") {
                A_eq = (J(1, n, 1), 0)
                b_eq = 1
            }
            else if (rts == "drs") {
                A_ineq = A_ineq \ (J(1, n, 1), 0)
                b_ineq = b_ineq \ 1
            }
            else if (rts == "irs") {
                A_ineq = A_ineq \ (J(1, n, -1), 0)
                b_ineq = b_ineq \ (-1)
            }
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setInequality(A_ineq, b_ineq)
            if (rows(A_eq) > 0) S.setEquality(A_eq, b_eq)
            S.setBounds(lb, ub)
            S.setTrace("off")
            
            rc = S.optimize()
            if (S.converged() && S.errorcode() == 0) {
                eff_pess = S.value()
            }
        }
        else {
            // Output-oriented pessimistic
            x0 = Xhi[k, ]'
            y0 = Ylo[k, ]'
            
            c = J(1, n, 0), -1
            A_ineq = (Xlo', J(m, 1, 0)) \ (-Yhi', y0)
            b_ineq = x0 \ J(s, 1, 0)
            
            lb = J(1, n, 0), 1
            ub = J(1, n+1, .)
            
            A_eq = J(0, n+1, .)
            b_eq = J(0, 1, .)
            if (rts == "vrs") {
                A_eq = (J(1, n, 1), 0)
                b_eq = 1
            }
            else if (rts == "drs") {
                A_ineq = A_ineq \ (J(1, n, 1), 0)
                b_ineq = b_ineq \ 1
            }
            else if (rts == "irs") {
                A_ineq = A_ineq \ (J(1, n, -1), 0)
                b_ineq = b_ineq \ (-1)
            }
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setInequality(A_ineq, b_ineq)
            if (rows(A_eq) > 0) S.setEquality(A_eq, b_eq)
            S.setBounds(lb, ub)
            S.setTrace("off")
            
            rc = S.optimize()
            if (S.converged() && S.errorcode() == 0) {
                eff_pess = 1 / (-S.value())
            }
        }
        
        // ================================================================
        // OPTIMISTIC EFFICIENCY
        // DMU k: best case (Xlo, Yhi)
        // Reference: worst case (Xhi, Ylo)
        // ================================================================
        if (orientation == "input") {
            x0 = Xlo[k, ]'
            y0 = Yhi[k, ]'
            
            c = J(1, n, 0), 1
            A_ineq = (Xhi', -x0) \ (-Ylo', J(s, 1, 0))
            b_ineq = J(m, 1, 0) \ (-y0)
            
            lb = J(1, n+1, 0)
            ub = J(1, n+1, .)
            
            A_eq = J(0, n+1, .)
            b_eq = J(0, 1, .)
            if (rts == "vrs") {
                A_eq = (J(1, n, 1), 0)
                b_eq = 1
            }
            else if (rts == "drs") {
                A_ineq = A_ineq \ (J(1, n, 1), 0)
                b_ineq = b_ineq \ 1
            }
            else if (rts == "irs") {
                A_ineq = A_ineq \ (J(1, n, -1), 0)
                b_ineq = b_ineq \ (-1)
            }
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setInequality(A_ineq, b_ineq)
            if (rows(A_eq) > 0) S.setEquality(A_eq, b_eq)
            S.setBounds(lb, ub)
            S.setTrace("off")
            
            rc = S.optimize()
            if (S.converged() && S.errorcode() == 0) {
                eff_opt = S.value()
            }
        }
        else {
            // Output-oriented optimistic
            x0 = Xlo[k, ]'
            y0 = Yhi[k, ]'
            
            c = J(1, n, 0), -1
            A_ineq = (Xhi', J(m, 1, 0)) \ (-Ylo', y0)
            b_ineq = x0 \ J(s, 1, 0)
            
            lb = J(1, n, 0), 1
            ub = J(1, n+1, .)
            
            A_eq = J(0, n+1, .)
            b_eq = J(0, 1, .)
            if (rts == "vrs") {
                A_eq = (J(1, n, 1), 0)
                b_eq = 1
            }
            else if (rts == "drs") {
                A_ineq = A_ineq \ (J(1, n, 1), 0)
                b_ineq = b_ineq \ 1
            }
            else if (rts == "irs") {
                A_ineq = A_ineq \ (J(1, n, -1), 0)
                b_ineq = b_ineq \ (-1)
            }
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setInequality(A_ineq, b_ineq)
            if (rows(A_eq) > 0) S.setEquality(A_eq, b_eq)
            S.setBounds(lb, ub)
            S.setTrace("off")
            
            rc = S.optimize()
            if (S.converged() && S.errorcode() == 0) {
                eff_opt = 1 / (-S.value())
            }
        }
        
        // Store results
        st_store(obs_k, st_varindex(eff_pess_var), eff_pess)
        st_store(obs_k, st_varindex(eff_opt_var), eff_opt)
    }
}

// Oriented SBM (Cooper, Seiford, and Tone 2007, eqs. 4.78--4.79).
// Input orientation maximizes normalized input slack and reports 1-z.
// Output orientation maximizes normalized output slack and reports 1/(1+z).
void _dea_ipm_sbm_oriented(string scalar inputs,
                           string scalar outputs,
                           string scalar eff_var,
                           string scalar touse,
                           string scalar rts,
                           string scalar orientation)
{
    real matrix X, Y, A_eq, A_ineq
    real colvector obs_indices, x0, y0, b_eq, b_ineq
    real rowvector c, lb, ub, sol, coef_sin, coef_sout
    real scalar n, m, s, k, j, obs_k, nv, obj, eff
    class LinearProgram scalar L

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))
    n = rows(X); m = cols(X); s = cols(Y); nv = n + m + s

    for (k=1; k<=n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k,.]'; y0 = Y[k,.]'
        coef_sin = J(1,m,0); coef_sout = J(1,s,0)
        for (j=1; j<=m; j++) if (x0[j]>0) coef_sin[j] = 1/(m*x0[j])
        for (j=1; j<=s; j++) if (y0[j]>0) coef_sout[j] = 1/(s*y0[j])

        if (orientation=="input") c = J(1,n,0), -coef_sin, J(1,s,0)
        else c = J(1,n,0), J(1,m,0), -coef_sout

        A_eq = (X', I(m), J(m,s,0))
        b_eq = x0
        A_eq = A_eq \ (Y', J(s,m,0), -I(s))
        b_eq = b_eq \ y0
        A_ineq = J(0,nv,0); b_ineq = J(0,1,0)

        if (rts=="vrs") {
            A_eq = A_eq \ (J(1,n,1), J(1,m+s,0))
            b_eq = b_eq \ 1
        }
        else if (rts=="drs") {
            A_ineq = (J(1,n,1), J(1,m+s,0)); b_ineq = 1
        }
        else if (rts=="irs") {
            A_ineq = (J(1,n,-1), J(1,m+s,0)); b_ineq = -1
        }

        lb = J(1,nv,0); ub = J(1,nv,.)
        L = LinearProgram()
        L.setMaxOrMin("min")
        L.setCoefficients(c)
        L.setEquality(A_eq,b_eq)
        if (rows(A_ineq)>0) L.setInequality(A_ineq,b_ineq)
        L.setBounds(lb,ub); L.setTrace("off")
        L.optimize()
        eff = .
        if (L.converged() && L.errorcode()==0) {
            obj = L.value()
            if (orientation=="input") eff = 1 + obj
            else if (1-obj>0) eff = 1/(1-obj)
        }
        st_store(obs_k,st_varindex(eff_var),eff)
    }
}

end
