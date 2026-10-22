*! version 1.4.1  21jul2026
*! dea_crosseff.ado — Cross-Efficiency DEA
*! Requires: Stata 16+ (LinearProgram)
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
*!   Doyle, J., & Green, R. (1994). Efficiency and cross-efficiency in
*!     DEA: Derivations, meanings and uses. Journal of the Operational
*!     Research Society, 45(5), 567-578. doi:10.1057/jors.1994.84
*!
*!   Lim, S., & Zhu, J. (2015). DEA cross-efficiency evaluation under
*!     variable returns to scale. Journal of the Operational Research
*!     Society, 66(3), 476-487. doi:10.1057/jors.2014.13
*!
*!   Sexton, T. R., Silkman, R. H., & Hogan, A. J. (1986). Data
*!     envelopment analysis: Critique and extensions. New Directions
*!     for Program Evaluation, 1986(32), 73-105. doi:10.1002/ev.1441
*!
*! Changelog:
*!   1.4.0  25apr2026  Implement method(benevolent/aggressive): 2-stage LP
*!   1.3.1  25apr2026  VRS CE peer-eval: use u*y/v*x only (no w), avoids negative CE
*!   1.3.0  25apr2026  Add VRS and output orientation support
*!   1.2.2  25apr2026  Fix string DMU: remove [varlist] from syntax (Stata defaults to
*!                     all dataset vars causing r(109)); Form A transform puts outputs
*!                     into outputs() option; markout excludes dmu variable
*!   1.1.0  11feb2026  Unified syntax: supports both
*!                        dea_crosseff x = y, [options]  (Form A)
*!                        dea_crosseff, inputs(x) outputs(y) [options] (Form B)
*!   1.0.0  06feb2026  Initial packaged version; based on earlier teaching/seminar versions

/**
 * CROSS-EFFICIENCY DEA
 * ====================
 * 
 * Cross-efficiency evaluates each DMU using the optimal weights
 * of all other DMUs, providing a peer evaluation mechanism.
 * 
 * Key concepts:
 *   - Self-appraisal: DMU evaluated with its own optimal weights
 *   - Peer-appraisal: DMU evaluated with other DMUs' weights
 *   - Cross-efficiency score: Average of peer appraisals
 * 
 * Advantages:
 *   - Provides unique ranking (unlike standard DEA)
 *   - Reduces unrealistic weight flexibility
 *   - Identifies "maverick" DMUs with atypical weights
 * 
 * Methods:
 *   - Aggressive: Secondary goal minimizes cross-efficiency of others
 *   - Benevolent: Secondary goal maximizes cross-efficiency of others
 *   - Arbitrary: Uses primary optimal weights only
 */

capture program drop dea_crosseff
program define dea_crosseff, rclass
    version 16.0

    // =========================================================================
    // Unified Syntax Parser
    //   Form A: dea_crosseff x1 x2 = y1 y2 [if] [in], [options]
    //     → transforms to: , inputs(x1 x2) outputs(y1 y2) [options]
    //   Form B: dea_crosseff [if] [in], inputs() outputs() [options]
    //
    // No [varlist] in syntax: avoids r(109) from Stata defaulting to all
    // dataset vars (including string dmu) when no varlist is specified.
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
        // `0' now contains: "outputs [, options]"
        local _cpos = strpos(`"`0'"', ",")
        if `_cpos' > 0 {
            local _io_outputs = strtrim(substr(`"`0'"', 1, `_cpos'-1))
            local _io_rest    = substr(`"`0'"', `_cpos'+1, .)
            local 0 `", inputs(`_io_inputs') outputs(`_io_outputs') `_io_rest'"'
        }
        else {
            local 0 `", inputs(`_io_inputs') outputs(`0')"'
        }
    }

    // No [varlist]: Stata defaults to ALL dataset vars when varlist is omitted,
    // causing r(109) if any string var (e.g. dmu) exists. Omit varlist entirely;
    // inputs/outputs handled via Inputs()/Outputs() options and Form A transform.
    syntax [if] [in], ///
        [DMU(varname)           ///
         Inputs(varlist numeric) ///
         Outputs(varlist numeric) ///
         RTS(string) ///
         ORT(string) ORIentation(string) ///
         METHOD(string) ///
         GENerate(string) ///
         TOLerance(real 1e-7) ///
         MATrix ///
         DETail ///
         DEBUG]

    // Validate required options (Form A already placed both into Inputs/Outputs options)
    if "`outputs'" == "" {
        di as error "outputs() required, or use: dea_crosseff inputs = outputs, [options]"
        exit 198
    }
    if "`inputs'" == "" {
        di as error "inputs() required, or use: dea_crosseff inputs = outputs, [options]"
        exit 198
    }

    // Validate dmu is a real variable (if provided)
    if "`dmu'" != "" {
        capture confirm variable `dmu'
        if _rc != 0 {
            di as error "dmu(`dmu') is not a variable in the dataset"
            exit 111
        }
    }
    
    // Reconcile ort() / orientation()
    if "`ort'" != "" & "`orientation'" != "" {
        di as error "Specify only one of: ort() or orientation()"
        exit 198
    }
    if "`ort'" != "" local orientation "`ort'"
    
    // Defaults
    if "`rts'" == "" local rts "crs"
    if "`orientation'" == "" local orientation "input"
    if "`method'" == "" local method "arbitrary"
    if "`generate'" == "" local generate "_ce_eff"
    
    // Abbreviations
    if inlist("`orientation'", "i", "in")  local orientation "input"
    if inlist("`orientation'", "o", "out") local orientation "output"
    if "`rts'" == "c" local rts "crs"
    if "`rts'" == "v" local rts "vrs"

    // Validate
    if !inlist("`rts'", "crs", "vrs") {
        di as error "rts() must be crs or vrs"
        exit 198
    }
    if !inlist("`orientation'", "input", "output") {
        di as error "orientation() must be input or output"
        exit 198
    }
    if !inlist("`method'", "arbitrary", "benevolent", "aggressive") {
        di as error "method() must be arbitrary, benevolent, or aggressive"
        exit 198
    }
    
    // Mark sample — do NOT include dmu in markout (string DMU causes r(109))
    marksample touse
    markout `touse' `inputs' `outputs'
    if "`dmu'" != "" {
        qui replace `touse' = 0 if missing(`dmu')   // handles both string("") and numeric(.)
    }
    
    // Count
    qui count if `touse'
    local n = r(N)
    
    // Display header
    di as text ""
    di as text "{hline 70}"
    di as text "Cross-Efficiency DEA"
    di as text "{hline 70}"
    di as text "Method:        " as result "`method'"
    di as text "RTS:           " as result "`rts'"
    di as text "Orientation:   " as result "`orientation'"
    di as text "DMUs:          " as result "`n'"
    di as text "Inputs:        " as result "`inputs'"
    di as text "Outputs:       " as result "`outputs'"
    di as text "{hline 70}"
    
    // Create output variables
    cap drop `generate'
    cap drop `generate'_self
    cap drop `generate'_maverick
    qui gen double `generate' = . if `touse'
    qui gen double `generate'_self = . if `touse'
    qui gen double `generate'_maverick = . if `touse'
    
    // Create DMU identifier if not provided
    if "`dmu'" == "" {
        tempvar dmu_temp
        qui gen `dmu_temp' = _n if `touse'
        local dmu "`dmu_temp'"
    }
    
    // Call Mata function.  The tempname keeps the optional matrix alive long
    // enough for the public r(ce_matrix) contract to be restored below.
    tempname ce_matrix_return
    mata: _dea_crosseff_compute("`inputs'", "`outputs'", ///
        "`dmu'", "`touse'", ///
        "`generate'", "`generate'_self", "`generate'_maverick", ///
        "`method'", "`rts'", "`orientation'", `tolerance', ///
        "`matrix'" != "", "`detail'" != "", "`debug'" != "", ///
        "`ce_matrix_return'")
    
    // Summary statistics
    qui summ `generate' if `touse'
    local mean_ce = r(mean)
    local min_ce = r(min)
    local max_ce = r(max)
    
    qui summ `generate'_self if `touse'
    local mean_self = r(mean)
    
    di as text ""
    di as text "{hline 70}"
    di as text "Summary"
    di as text "{hline 70}"
    di as text "Mean self-appraisal:      " as result %9.4f `mean_self'
    di as text "Mean cross-efficiency:    " as result %9.4f `mean_ce'
    di as text "Min cross-efficiency:     " as result %9.4f `min_ce'
    di as text "Max cross-efficiency:     " as result %9.4f `max_ce'
    di as text "{hline 70}"
    
    // Return values
    return local method "`method'"
    return local rts "`rts'"
    return local orientation "`orientation'"
    return scalar n = `n'
    return scalar mean_ce = `mean_ce'
    return scalar mean_self = `mean_self'
    if "`matrix'" != "" {
        return matrix ce_matrix = `ce_matrix_return'
    }
end


// ============================================================================
// MATA FUNCTIONS
// ============================================================================

capture mata: mata drop _dea_dmu_read()
capture mata: mata drop _dea_crosseff_compute()
capture mata: mata drop _dea_crosseff_multiplier()
capture mata: mata drop _dea_rank()

mata:
mata set matastrict on

// ─────────────────────────────────────────────────────────────────────────────
// _dea_dmu_read()
//
// Reads a DMU variable that may be string or numeric.
// Returns: dmu_ids (numeric, 1-indexed for strings; original values for numeric)
//          dmu_labels (string representations for display)
// ─────────────────────────────────────────────────────────────────────────────
void _dea_dmu_read(string scalar  dmu_var,
                   string scalar  touse,
                   real colvector dmu_ids,
                   string colvector dmu_labels)
{
    real scalar n, i, j
    string colvector dmu_str, uniq

    if (substr(st_vartype(dmu_var), 1, 3) == "str") {
        dmu_str    = st_sdata(., dmu_var, touse)
        n          = rows(dmu_str)
        uniq       = sort(uniqrows(dmu_str), 1)
        dmu_ids    = J(n, 1, .)
        dmu_labels = dmu_str
        for (i = 1; i <= n; i++) {
            for (j = 1; j <= rows(uniq); j++) {
                if (dmu_str[i] == uniq[j]) {
                    dmu_ids[i] = j
                    break
                }
            }
        }
    }
    else {
        dmu_ids    = st_data(., dmu_var, touse)
        dmu_labels = strofreal(dmu_ids)
    }
}


/**
 * _dea_crosseff_compute()
 *
 * Computes Cross-Efficiency scores. Supports string and numeric DMU identifiers.
 */
void _dea_crosseff_compute(string scalar inputs,
                            string scalar outputs,
                            string scalar dmu_var,
                            string scalar touse,
                            string scalar eff_ce_var,
                            string scalar eff_self_var,
                            string scalar eff_maverick_var,
                            string scalar method,
                            string scalar rts,
                            string scalar orientation,
                            real scalar tol,
                            real scalar show_matrix,
                            real scalar show_detail,
                            real scalar show_debug,
                            string scalar ce_matrix_name)
{
    // Load data
    real matrix X, Y
    real colvector dmu_ids, obs_indices
    string colvector dmu_labels
    string rowvector input_names, output_names

    input_names  = tokens(inputs)
    output_names = tokens(outputs)

    X           = st_data(., input_names,  touse)
    Y           = st_data(., output_names, touse)
    obs_indices = selectindex(st_data(., touse))

    // Read DMU — handles string and numeric
    _dea_dmu_read(dmu_var, touse, dmu_ids, dmu_labels)
    
    real scalar n, m, s
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    if (show_debug) {
        printf("  n=%g, m=%g (inputs), s=%g (outputs)\n", n, m, s)
    }
    
    // Step 1: Compute optimal weights for each DMU using multiplier model
    // For input-oriented CRS: max sum(u*y0) s.t. sum(v*x0)=1, sum(u*yj)-sum(v*xj)<=0
    
    real matrix V_opt, U_opt   // Optimal weights: V (inputs), U (outputs)
    real colvector W_opt       // VRS free intercept (0 for CRS)
    V_opt = J(n, m, .)
    U_opt = J(n, s, .)
    W_opt = J(n, 1, 0)

    real colvector self_eff
    self_eff = J(n, 1, .)
    
    real scalar i, j
    real rowvector x0, y0, v_i, u_i
    real scalar eff_i, w_i
    
    for (i = 1; i <= n; i++) {
        x0 = X[i, .]
        y0 = Y[i, .]
        
        // Solve multiplier model (Stage 1 + optional Stage 2 for benevolent/aggressive)
        _dea_crosseff_multiplier(x0, y0, X, Y, rts, orientation, method, tol, v_i, u_i, w_i, eff_i)

        V_opt[i, .] = v_i
        U_opt[i, .] = u_i
        W_opt[i]    = w_i
        self_eff[i] = eff_i
        
        if (show_debug) {
            printf("  DMU %g: self-eff = %9.4f\n", dmu_ids[i], eff_i)
        }
    }
    
    // Step 2: Compute cross-efficiency matrix
    // CE[i,j] = efficiency of DMU j using weights of DMU i
    real matrix CE
    CE = J(n, n, .)
    
    real rowvector v_i2, u_i2
    real scalar num, denom

    // Peer-evaluation CE formula: CE[k,j] = (u_k*y_j) / (v_k*x_j).
    // Under VRS, the free intercept w remains in the self-appraisal multiplier
    // LP. It is excluded from peer ratios because intercept-driven peer scores
    // can become negative under VRS, which makes the cross-efficiency average
    // hard to interpret.
    for (i = 1; i <= n; i++) {
        v_i2 = V_opt[i, .]
        u_i2 = U_opt[i, .]

        for (j = 1; j <= n; j++) {
            num   = u_i2 * Y[j, .]'
            denom = v_i2 * X[j, .]'
            CE[i, j] = (abs(denom) > tol ? num / denom : .)
        }
    }
    
    // Step 3: Compute cross-efficiency scores (column means excluding self)
    real colvector ce_scores
    ce_scores = J(n, 1, .)
    
    real scalar sum_ce, count_ce
    
    for (j = 1; j <= n; j++) {
        sum_ce = 0
        count_ce = 0
        
        for (i = 1; i <= n; i++) {
            if (i != j && CE[i, j] != .) {
                sum_ce = sum_ce + CE[i, j]
                count_ce++
            }
        }
        
        if (count_ce > 0) {
            ce_scores[j] = sum_ce / count_ce
        }
    }
    
    // Step 4: Compute maverick index (self-appraisal vs cross-efficiency)
    real colvector maverick
    maverick = self_eff :/ ce_scores
    
    // Compute max label width for aligned display
    real scalar lw
    string scalar fmt_dmu
    lw = 10
    for (i = 1; i <= n; i++) {
        if (strlen(dmu_labels[i]) > lw) lw = strlen(dmu_labels[i])
    }
    fmt_dmu = "%" + strofreal(lw) + "s"

    // Display results — use string labels
    printf("\n")
    printf("Cross-Efficiency Results\n")
    printf("{hline 70}\n")
    printf(fmt_dmu + " | %10s %10s %10s %10s\n",
           "DMU", "Self", "Cross-Eff", "Maverick", "Rank")
    printf(fmt_dmu + " | %10s %10s %10s %10s\n",
           "-":*lw, "----------", "----------", "----------", "----------")

    // Rank by cross-efficiency
    real colvector ranks
    ranks = _dea_rank(ce_scores)

    for (i = 1; i <= n; i++) {
        printf(fmt_dmu + " | %10.4f %10.4f %10.4f %10.0f\n",
               dmu_labels[i], self_eff[i], ce_scores[i], maverick[i], ranks[i])
    }
    printf("{hline 70}\n")

    // Show cross-efficiency matrix if requested
    if (show_matrix) {
        printf("\nCross-Efficiency Matrix (row i evaluates column j)\n")
        printf("{hline 70}\n")

        // Header
        printf(fmt_dmu + " |", "Evaluator")
        for (j = 1; j <= n; j++) {
            printf(" %7s", dmu_labels[j])
        }
        printf("\n")

        printf(fmt_dmu + " |", "-":*lw)
        for (j = 1; j <= n; j++) {
            printf(" -------")
        }
        printf("\n")

        // Data
        for (i = 1; i <= n; i++) {
            printf(fmt_dmu + " |", dmu_labels[i])
            for (j = 1; j <= n; j++) {
                printf(" %7.4f", CE[i, j])
            }
            printf("\n")
        }
        printf("{hline 70}\n")
        st_matrix(ce_matrix_name, CE)
    }
    
    // Store results in Stata
    for (i = 1; i <= n; i++) {
        st_store(obs_indices[i], st_varindex(eff_ce_var), ce_scores[i])
        st_store(obs_indices[i], st_varindex(eff_self_var), self_eff[i])
        st_store(obs_indices[i], st_varindex(eff_maverick_var), maverick[i])
    }
}


/**
 * _dea_crosseff_multiplier()
 *
 * Solves the multiplier (dual) DEA model to get optimal weights.
 * 
 * For input-oriented CRS:
 *   max  sum(u * y0)
 *   s.t. sum(v * x0) = 1
 *        sum(u * yj) - sum(v * xj) <= 0  for all j
 *        u, v >= epsilon
 */
void _dea_crosseff_multiplier(real rowvector x0,
                               real rowvector y0,
                               real matrix X,
                               real matrix Y,
                               string scalar rts,
                               string scalar orientation,
                               string scalar method,
                               real scalar tol,
                               real rowvector v_opt,
                               real rowvector u_opt,
                               real scalar w_opt,
                               real scalar eff)
{
    real scalar n, m, s
    n = rows(X)
    m = cols(X)
    s = cols(Y)

    real scalar is_vrs, is_output, nvars, eps
    real rowvector c, lb, ub
    real matrix    A_ineq, A_eq
    real colvector b_ineq, b_eq
    class LinearProgram scalar S
    real rowvector sol

    eps       = 1e-6
    is_vrs    = (rts == "vrs")
    is_output = (orientation == "output")
    nvars     = m + s + is_vrs        // add 1 free variable w if VRS

    // ── Objective (Stage 1: self-efficiency) ─────────────────────────
    if (!is_output) {
        // Input: min −u·y₀ [− w]
        c = (J(1, m, 0), -y0)
        if (is_vrs) c = (c, -1)
    }
    else {
        // Output: min v·x₀ [+ w]
        c = (x0, J(1, s, 0))
        if (is_vrs) c = (c, 1)
    }

    // ── Equality constraint ──────────────────────────────────────────
    if (!is_output) {
        // Input: v·x₀ = 1
        A_eq = !is_vrs ? (x0, J(1, s, 0)) : (x0, J(1, s, 0), 0)
    }
    else {
        // Output: u·y₀ = 1
        A_eq = !is_vrs ? (J(1, m, 0), y0) : (J(1, m, 0), y0, 0)
    }
    b_eq = 1

    // ── Inequality constraints: [−X, Y {, ±1}]·z ≤ 0 ────────────────
    if (!is_vrs) {
        A_ineq = (-X, Y)                          // CRS: −v·xⱼ + u·yⱼ ≤ 0
    }
    else if (!is_output) {
        A_ineq = (-X, Y, J(n, 1,  1))            // Input VRS:  +w term (u·y+w ≤ v·x)
    }
    else {
        A_ineq = (-X, Y, J(n, 1, -1))            // Output VRS: −w term (u·y−w ≤ v·x)
    }
    b_ineq = J(n, 1, 0)

    // ── Bounds: v, u ≥ ε; w free ────────────────────────────────────
    lb = !is_vrs ? J(1, nvars, eps) : (J(1, m+s, eps), .)
    ub = J(1, nvars, .)

    // ── Stage 1: Solve self-efficiency ───────────────────────────────
    S = LinearProgram()
    S.setMaxOrMin("min")
    S.setCoefficients(c)
    S.setEquality(A_eq, b_eq)
    S.setInequality(A_ineq, b_ineq)
    S.setBounds(lb, ub)
    S.setTrace("off")

    (void) S.optimize()

    if (S.converged() && S.errorcode() == 0) {
        sol   = S.parameters()
        v_opt = sol[1..m]
        u_opt = sol[(m+1)..(m+s)]
        w_opt = (is_vrs ? sol[nvars] : 0)
        if (!is_output) {
            eff = -S.value()           // Input:  −(−u·y₀ [−w]) = u·y₀ [+w]
        }
        else {
            eff = (abs(S.value()) > tol ? 1 / S.value() : .)   // Output: 1/(v·x₀ [+w])
        }
    }
    else {
        v_opt = J(1, m, .)
        u_opt = J(1, s, .)
        w_opt = .
        eff   = .
        return
    }

    // ── Stage 2: Benevolent / Aggressive secondary LP ────────────────
    // Keep eff from Stage 1; re-optimize weights to max/min peer CE.
    // For input-oriented: fix u·y₀[+w]=eff*, then max/min u·Σ_{j≠k}y_j
    // For output-oriented: fix v·x₀[+w]=1/eff*, then min/max v·Σ_{j≠k}x_j
    if (method != "arbitrary" && eff != .) {
        real scalar eff_prime
        real rowvector c2, A_eq2_row
        real matrix    A_eq2
        real colvector b_eq2
        class LinearProgram scalar S2
        real rowvector sol2

        eff_prime = eff
        if (!is_output & eff_prime > 1 - 100*tol & eff_prime < 1 + 100*tol) eff_prime = 1

        // Second equality: lock self-efficiency at eff_prime
        if (!is_output) {
            // u·y₀ [+w] = eff_prime
            A_eq2_row = !is_vrs ? (J(1,m,0), y0) : (J(1,m,0), y0, 1)
            b_eq2 = (1 \ eff_prime)
        }
        else {
            // v·x₀ [+w] = 1/eff_prime
            A_eq2_row = !is_vrs ? (x0, J(1,s,0)) : (x0, J(1,s,0), 1)
            b_eq2 = (1 \ (1/eff_prime))
        }
        A_eq2 = (A_eq \ A_eq2_row)

        // Secondary objective
        real rowvector sum_y, sum_x
        sum_y = colsum(Y) - y0    // Σ_{j≠k} y_j  (1×s)
        sum_x = colsum(X) - x0    // Σ_{j≠k} x_j  (1×m)

        if (!is_output) {
            // Sexton et al. secondary goal uses the peer constraint slack
            // d_j = v*x_j - u*y_j.  Benevolent minimizes sum d_j; aggressive
            // maximizes it, while preserving the Stage-1 self-efficiency.
            c2 = !is_vrs ? (sum_x, -sum_y) : (sum_x, -sum_y, 0)
            if (method == "aggressive") c2 = -c2
        }
        else {
            // Benevolent: min v·sum_x  (minimise peer input → higher CE for peers)
            // Aggressive: max v·sum_x  → min −v·sum_x
            c2 = !is_vrs ? (sum_x, J(1,s,0)) : (sum_x, J(1,s,0), 0)
            if (method == "aggressive") c2 = -c2
        }

        S2 = LinearProgram()
        S2.setMaxOrMin("min")
        S2.setCoefficients(c2)
        S2.setEquality(A_eq2, b_eq2)
        S2.setInequality(A_ineq, b_ineq)
        S2.setBounds(lb, ub)
        S2.setTrace("off")

        (void) S2.optimize()

        if (S2.converged() && S2.errorcode() == 0) {
            sol2  = S2.parameters()
            v_opt = sol2[1..m]
            u_opt = sol2[(m+1)..(m+s)]
            w_opt = (is_vrs ? sol2[nvars] : 0)
            // eff unchanged — keeps Stage 1 self-efficiency
        }
        // If Stage 2 fails: keep Stage 1 solution (already in v_opt, u_opt, w_opt)
    }
}


/**
 * _dea_rank()
 *
 * Returns ranks (1 = highest value)
 */
real colvector _dea_rank(real colvector x)
{
    real scalar n, i, j, rank
    real colvector ranks
    
    n = rows(x)
    ranks = J(n, 1, .)
    
    for (i = 1; i <= n; i++) {
        if (x[i] == .) {
            ranks[i] = .
            continue
        }
        
        rank = 1
        for (j = 1; j <= n; j++) {
            if (j != i && x[j] != . && x[j] > x[i]) {
                rank++
            }
        }
        ranks[i] = rank
    }
    
    return(ranks)
}

end
