*! _dea_spbland.ado
*! Single-Phase Bland's Rule Simplex for DEA (model(te) only)
*! Version 1.0.0 - 29apr2026
*!
*! AUTHORS:
*! --------
*!   Choonjoo Lee (corresponding author)
*!     Korea National Defense University
*!     Seoul, Korea
*!     bloom.rampike@gmail.com; bloom.rampike@gmail.com
*!
*! ALGORITHM:
*! ----------
*! Single-phase primal simplex exploiting DEA structural feasibility:
*! each DMU is always a feasible starting point (lambda_o=1 satisfies
*! all constraints), so Phase 1 is unnecessary.
*!
*!   - lusolve(B, b) per iteration   (full LU re-factorisation each step)
*!   - Bland's rule for both entering AND leaving variable selection
*!     (min-index tie-breaking guarantees finite termination)
*!   - No Harris ratio test, no eta-vector updates (simpler code)
*!
*! This differs from _dea_simplex.ado (Two-Phase Revised Simplex) in that:
*!   - No artificial variables / Phase 1
*!   - B^{-1} maintained by full lusolve each iteration (not eta updates)
*!   - Only model(te) supported; all RTS (crs/vrs/irs/drs); both orientations
*!
*! REFERENCES:
*! -----------
*! Bland, R. G. (1977). New finite pivoting rules for the simplex method.
*!   Mathematics of Operations Research, 2(2), 103-107.
*!
*! Charnes, A., Cooper, W. W., & Rhodes, E. (1978). Measuring the efficiency
*!   of decision making units. European Journal of Operational Research, 2, 429-444.
*!
*! Banker, R. D., Charnes, A., & Cooper, W. W. (1984). Some models for
*!   estimating technical and scale inefficiencies in data envelopment analysis.
*!   Management Science, 30(9), 1078-1092.
*!
*! LP FORMULATION (input-oriented, standard-max form):
*! ---------------------------------------------------
*!   Variables: [theta, lam_1..lam_n, s_1^-..s_M^-, s_1^+..s_N^+, (s_conv)]
*!   Maximize:  -theta
*!   s.t.
*!     Input  m : -x_om*theta + X_m'*lam + s_m^-          = 0   (m=1..M)
*!     Output n :             + Y_n'*lam          - s_n^+  = y_on  (n=1..N)
*!     Convexity:               sum(lam)                    = 1   (VRS only)
*!               :              sum(lam) - s_conv           = 1   (IRS)
*!               :              sum(lam) + s_conv           = 1   (DRS)
*!   Starting basis: {theta(r1), s2^-..sM^-(r2..M),
*!                    s1^+..sN^+(rM+1..M+N), lam_o(rM+N[+1])}
*!
*! LP FORMULATION (output-oriented, standard-max form):
*! ----------------------------------------------------
*!   Maximize:  +phi
*!   s.t.
*!     Input  m :              X_m'*lam + s_m^-          = x_om  (m=1..M)
*!     Output n : -y_on*phi  + Y_n'*lam          -s_n^+  = 0     (n=1..N)
*!     Convexity: (same as input-oriented)
*!   Starting basis: {s1^-..sM^-(r1..M), phi(rM+1),
*!                    s2^+..sN^+(rM+2..M+N), lam_o(rM+N[+1])}
*!   Note: CRS output-oriented uses sN^+ in last row (no lam_o in basis);
*!         phi starts at 0 (degenerate but feasible).

version 16.0

// ============================================================================
// Clean up existing Mata functions
// ============================================================================
capture mata: mata drop _sp_bland()
capture mata: mata drop _dea_spbland_te_dmu()
capture mata: mata drop _dea_spbland_envelopment()

mata:
mata set matastrict on

// ============================================================================
// 1. CORE SIMPLEX: Single-Phase Bland's Rule
// ============================================================================
// Input:
//   A      : m x n constraint matrix (standard form: Ax = b, x >= 0)
//   b      : m x 1 RHS (must be >= 0 for the provided starting basis)
//   c      : 1 x n objective row vector (maximize c'x)
//   basis0 : 1 x m initial basis (column indices of A, 1-based)
// Output:
//   1 x n primal solution vector
//   Returns J(1,n,.) if infeasible, unbounded, or max iterations reached

real rowvector _sp_bland(real matrix  A,
                          real colvector b,
                          real rowvector c,
                          real rowvector basis0)
{
    real scalar    m, n, i, j, k
    real scalar    j_enter, i_leave
    real scalar    ratio, min_ratio, EPS, in_basis, rc_j, iter
    real rowvector basis, c_B, x_sol
    real colvector x_B, d_col, pi_col
    real matrix    B, Bt

    m     = rows(A)
    n     = cols(A)
    EPS   = 1e-9
    basis = basis0   // 1 x m working copy

    for (iter = 1; iter <= 300 * n; iter++) {

        // ------------------------------------------------------------------
        // Build basis matrix B (m x m) and cost-basis c_B (1 x m)
        // ------------------------------------------------------------------
        B   = J(m, m, 0)
        c_B = J(1, m, 0)
        for (j = 1; j <= m; j++) {
            B[., j] = A[., basis[j]]
            c_B[j]  = c[basis[j]]
        }

        // ------------------------------------------------------------------
        // Primal solution: x_B = B \ b
        // ------------------------------------------------------------------
        x_B = lusolve(B, b)

        // ------------------------------------------------------------------
        // Dual multipliers: pi = (B')^{-1} * c_B'
        // Reduced cost of column j: rc_j = c[j] - pi' * A[.,j]
        // ------------------------------------------------------------------
        Bt     = B'
        pi_col = lusolve(Bt, c_B')

        // ------------------------------------------------------------------
        // Entering variable (Bland's rule: smallest nonbasic index with rc>0)
        // ------------------------------------------------------------------
        j_enter = 0
        for (j = 1; j <= n; j++) {
            // Skip if j is already in basis
            in_basis = 0
            for (k = 1; k <= m; k++) {
                if (basis[k] == j) {
                    in_basis = 1
                    break
                }
            }
            if (!in_basis) {
                rc_j = c[j]
                for (i = 1; i <= m; i++) rc_j = rc_j - pi_col[i] * A[i, j]
                if (rc_j > EPS) {
                    j_enter = j
                    break    // Bland's: stop at first (smallest-index) improving col
                }
            }
        }
        if (j_enter == 0) break   // No improving column → optimal

        // ------------------------------------------------------------------
        // Pivot column: d = B \ A[., j_enter]
        // ------------------------------------------------------------------
        d_col = lusolve(B, A[., j_enter])

        // ------------------------------------------------------------------
        // Leaving variable (min-ratio test; Bland's tie-breaking: min basis index)
        // ------------------------------------------------------------------
        i_leave   = 0
        min_ratio = .
        for (i = 1; i <= m; i++) {
            if (d_col[i] > EPS) {
                ratio = x_B[i] / d_col[i]
                if (i_leave == 0 | ratio < min_ratio - EPS) {
                    min_ratio = ratio
                    i_leave   = i
                }
                else if (abs(ratio - min_ratio) <= EPS) {
                    if (basis[i] < basis[i_leave]) i_leave = i
                }
            }
        }
        if (i_leave == 0) return(J(1, n, .))   // Unbounded

        basis[i_leave] = j_enter
    }

    // ------------------------------------------------------------------
    // Extract final primal solution
    // ------------------------------------------------------------------
    B = J(m, m, 0)
    for (j = 1; j <= m; j++) B[., j] = A[., basis[j]]
    x_B   = lusolve(B, b)
    x_sol = J(1, n, 0)
    for (i = 1; i <= m; i++) {
        if (x_B[i] > EPS) x_sol[basis[i]] = x_B[i]
    }
    return(x_sol)
}


// ============================================================================
// 2. DEA-LP BUILDER & SOLVER FOR ONE DMU
// ============================================================================
// Solves the DEA TE problem for one DMU using _sp_bland().
//
// Arguments:
//   X      : n x M input matrix  (all n DMUs)
//   Y      : n x N output matrix (all n DMUs)
//   x_o    : 1 x M input  vector of DMU under evaluation
//   y_o    : 1 x N output vector of DMU under evaluation
//   o_idx  : row index of DMU o in X/Y (1-based; used for starting basis)
//   ort    : "input" or "output"
//   rts    : "crs" | "vrs" | "irs" | "drs"
// Returns:
//   theta* (input) or phi* (output); missing on failure

real scalar _dea_spbland_te_dmu(
    real matrix X, real matrix Y,
    real rowvector x_o, real rowvector y_o,
    real scalar o_idx,
    string scalar ort, string scalar rts)
{
    real scalar    n, M, N, nrows, K, m, nn, j
    real scalar    has_conv, is_irs, is_drs
    real scalar    conv_col, obj_val
    real matrix    Aeq
    real colvector beq
    real rowvector c_obj, basis0, x_opt

    n        = rows(X)
    M        = cols(X)
    N        = cols(Y)
    has_conv = (rts != "crs")
    is_irs   = (rts == "irs")
    is_drs   = (rts == "drs")
    nrows    = M + N + (has_conv ? 1 : 0)

    // Variable layout:
    //   theta/phi : col 1
    //   lambda_j  : col 1+j          (j = 1..n)
    //   sminus_m  : col 1+n+m        (m = 1..M)
    //   splus_n   : col 1+n+M+n      (n = 1..N)
    //   s_conv    : col 1+n+M+N+1    (IRS surplus OR DRS slack; absent for CRS/VRS)
    K        = 1 + n + M + N + ((is_irs | is_drs) ? 1 : 0)
    conv_col = 1 + n + M + N + 1   // only valid when is_irs or is_drs

    Aeq   = J(nrows, K, 0)
    beq   = J(nrows, 1, 0)
    c_obj = J(1, K, 0)

    // ------------------------------------------------------------------
    // Build LP matrices
    // ------------------------------------------------------------------
    if (ort == "input") {
        // maximize -theta
        c_obj[1] = -1

        // Input rows m = 1..M:  -x_om*theta + X_m'*lam + sminus_m = 0
        for (m = 1; m <= M; m++) {
            Aeq[m, 1]       = -x_o[m]
            for (j = 1; j <= n; j++) Aeq[m, 1+j] = X[j, m]
            Aeq[m, 1+n+m]  = 1
            beq[m]          = 0
        }
        // Output rows n = 1..N:  Y_n'*lam - splus_n = y_on
        for (nn = 1; nn <= N; nn++) {
            for (j = 1; j <= n; j++) Aeq[M+nn, 1+j] = Y[j, nn]
            Aeq[M+nn, 1+n+M+nn] = -1
            beq[M+nn]            = y_o[nn]
        }

    } else {   // output-oriented
        // maximize phi
        c_obj[1] = 1

        // Input rows m = 1..M:  X_m'*lam + sminus_m = x_om
        for (m = 1; m <= M; m++) {
            for (j = 1; j <= n; j++) Aeq[m, 1+j] = X[j, m]
            Aeq[m, 1+n+m] = 1
            beq[m]         = x_o[m]
        }
        // Output rows n = 1..N:  -y_on*phi + Y_n'*lam - splus_n = 0
        for (nn = 1; nn <= N; nn++) {
            Aeq[M+nn, 1]         = -y_o[nn]
            for (j = 1; j <= n; j++) Aeq[M+nn, 1+j] = Y[j, nn]
            Aeq[M+nn, 1+n+M+nn] = -1
            beq[M+nn]            = 0
        }
    }

    // Convexity row (VRS / IRS / DRS)
    if (has_conv) {
        for (j = 1; j <= n; j++) Aeq[M+N+1, 1+j] = 1
        beq[M+N+1] = 1
        if (is_irs) Aeq[M+N+1, conv_col] = -1   // Slam >= 1:  Slam - s_conv = 1
        if (is_drs) Aeq[M+N+1, conv_col] =  1   // Slam <= 1:  Slam + s_conv = 1
    }

    // ------------------------------------------------------------------
    // Build starting basis (1 x nrows)
    // ------------------------------------------------------------------
    // Strategy: exploit the fact that lambda_o = 1, theta/phi = 1 is always
    // feasible (each DMU is on its own "frontier" when evaluated alone).
    // The starting BFS is degenerate (many slacks = 0) but valid.
    // ------------------------------------------------------------------
    basis0 = J(1, nrows, 0)

    if (ort == "input") {
        // Row 1               → theta  (col 1)
        // Rows 2..M           → sminus_2..sminus_M
        // Rows M+1..M+N[−1]  → splus_1..splus_N[−1]
        // Row M+N[+1]         → lambda_o
        basis0[1] = 1                                   // theta
        for (m = 2; m <= M; m++) basis0[m] = 1+n+m     // sminus_m

        if (rts == "crs") {
            // M+N constraints: last output row gets lambda_o
            for (nn = 1; nn <= N-1; nn++) basis0[M+nn] = 1+n+M+nn
            basis0[M+N] = 1 + o_idx                    // lambda_o
        } else {
            // M+N+1 constraints: all output rows get splus, convexity gets lambda_o
            for (nn = 1; nn <= N; nn++) basis0[M+nn] = 1+n+M+nn
            basis0[M+N+1] = 1 + o_idx                  // lambda_o
        }

    } else {   // output-oriented
        // Rows 1..M           → sminus_1..sminus_M
        // Row M+1             → phi  (col 1)
        // Rows M+2..M+N[−1]  → splus_2..splus_N[−1]
        // Row M+N[+1]         → lambda_o  (VRS/IRS/DRS) OR splus_N  (CRS)
        for (m = 1; m <= M; m++) basis0[m] = 1+n+m    // sminus_m
        basis0[M+1] = 1                                 // phi

        if (rts == "crs") {
            // CRS: no convexity row; lambda_o cannot be in basis (linearly
            // dependent with phi in output rows without the convexity anchor).
            // Use splus_N in the last row instead. Starting BFS: phi=0.
            for (nn = 2; nn <= N; nn++) basis0[M+nn] = 1+n+M+nn
            // basis0[M+N] already set to splus_N by the loop above (nn=N)
        } else {
            // VRS/IRS/DRS: convexity row pins lambda_o separately from phi
            for (nn = 2; nn <= N; nn++) basis0[M+nn] = 1+n+M+nn
            basis0[M+N+1] = 1 + o_idx                  // lambda_o
        }
    }

    // ------------------------------------------------------------------
    // Edge-case guard: if any basis entry is still 0, abort
    // (Happens if M=0 or N=0, which DEA doesn't allow)
    // ------------------------------------------------------------------
    for (j = 1; j <= nrows; j++) {
        if (basis0[j] == 0) return(.)
    }

    // ------------------------------------------------------------------
    // Solve
    // ------------------------------------------------------------------
    x_opt = _sp_bland(Aeq, beq, c_obj, basis0)
    if (missing(x_opt[1])) return(.)

    if (ort == "input") {
        obj_val = 0
        for (j = 1; j <= K; j++) obj_val = obj_val + c_obj[j]*x_opt[j]
        return(-obj_val)         // theta* = -max(-theta), in [0,1]
    } else {
        if (x_opt[1] == 0) return(.)    // degenerate / infeasible guard
        return(1 / x_opt[1])   // TE = 1/phi*, in [0,1] (Farrell convention)
    }
}


// ============================================================================
// 3. MAIN ENTRY: Loop over all DMUs
// ============================================================================
// Arguments (Stata variable names as strings):
//   inputs   : space-separated input variable names
//   outputs  : space-separated output variable names
//   effvar   : efficiency score variable (already exists, fill in values)
//   tousevar : 0/1 indicator for sample inclusion
//   ort      : "input" or "output"
//   rts      : "crs" | "vrs" | "irs" | "drs"
//   tol      : convergence tolerance (used only for missing-value threshold)
//   trace    : "" or "trace" (verbose iteration output)

void _dea_spbland_envelopment(
    string scalar inputs,   string scalar outputs,
    string scalar effvar,   string scalar tousevar,
    string scalar ort,      string scalar rts,
    real scalar   tol,      string scalar trace)
{
    real scalar      n, M, N, i, o, te_val, st_touse, st_eff
    real matrix      X, Y
    real rowvector   x_o, y_o, inp_idx, out_idx
    real colvector   touse_col, idx
    string rowvector inp_names, out_names

    // ------------------------------------------------------------------
    // Load data from Stata
    // ------------------------------------------------------------------
    inp_names = tokens(inputs)
    out_names = tokens(outputs)
    M         = length(inp_names)
    N         = length(out_names)

    st_touse  = st_varindex(tousevar)
    touse_col = st_data(., st_touse)

    // Build index list of DMUs in sample
    n   = 0
    idx = J(rows(touse_col), 1, 0)
    for (i = 1; i <= rows(touse_col); i++) {
        if (touse_col[i] != 0) {
            n++
            idx[n] = i
        }
    }
    idx = idx[1..n]   // trim

    if (n == 0) return

    // Collect X and Y matrices (n x M, n x N)
    inp_idx = J(1, M, 0)
    out_idx = J(1, N, 0)
    for (i = 1; i <= M; i++) inp_idx[i] = st_varindex(inp_names[i])
    for (i = 1; i <= N; i++) out_idx[i] = st_varindex(out_names[i])

    X = st_data(idx, inp_idx)   // n x M
    Y = st_data(idx, out_idx)   // n x N

    // Efficiency variable index
    st_eff = st_varindex(effvar)

    // ------------------------------------------------------------------
    // Solve for each DMU
    // ------------------------------------------------------------------
    for (o = 1; o <= n; o++) {
        x_o = X[o, .]   // 1 x M
        y_o = Y[o, .]   // 1 x N

        if (trace != "") {
            printf("  SP-Bland: DMU %g of %g\n", o, n)
            displayflush()
        }

        te_val = _dea_spbland_te_dmu(X, Y, x_o, y_o, o, ort, rts)

        // Store result
        st_store(idx[o], st_eff, te_val)
    }
}

end
