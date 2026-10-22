*! version 2.2.0  03may2026
*! dea_crossgame.ado — Cross-Efficiency Game DEA
*! Includes: Aggressive, Benevolent, Nash Equilibrium strategies
*! Requires: Stata 16+ (LinearProgram)
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
*!   Doyle, J., & Green, R. (1994). Efficiency and cross-efficiency in
*!     DEA: Derivations, meanings and uses. Journal of the Operational
*!     Research Society, 45(5), 567-578. doi:10.1057/jors.1994.84
*!
*!   Liang, L., Wu, J., Cook, W. D., & Zhu, J. (2008). The DEA game
*!     cross-efficiency model and its Nash equilibrium. Operations
*!     Research, 56(5), 1278-1288. doi:10.1287/opre.1070.0487
*!
*!   Wu, J., Liang, L., & Yang, F. (2009). Achievement and benchmarking
*!     of countries at the Summer Olympics using cross efficiency
*!     evaluation method. European Journal of Operational Research,
*!     197(2), 722-730. doi:10.1016/j.ejor.2008.06.030
*!
*! Changelog:
*!   2.2.0  03may2026  Implement Liang et al. (2008) model (5) for strategy(nash)
*!   2.1.0  28apr2026  Add Form A syntax (x1 x2 = y1 y2)

/**
 * CROSS-EFFICIENCY GAME DEA
 * =========================
 * 
 * Models DMU competition/cooperation through cross-efficiency evaluation.
 * 
 * Strategies:
 *   - Aggressive: Minimize others' cross-efficiency (competitive)
 *   - Benevolent: Maximize others' cross-efficiency (cooperative)
 *   - Nash: Find equilibrium where no DMU can unilaterally improve
 * 
 * References:
 *   - Doyle & Green (1994): Cross-efficiency derivations, meanings, and uses
 *   - Liang et al. (2008): DEA game cross-efficiency model and Nash equilibrium
 *   - Wu et al. (2009): Cross-efficiency benchmarking application
 */

capture program drop dea_crossgame
program define dea_crossgame, rclass
    version 16.0

    // =========================================================================
    // Unified Syntax Parser (Form A / Form B)
    //   Form A: dea_crossgame x1 x2 = y1 y2 [if] [in], [options]
    //   Form B: dea_crossgame [if] [in], inputs(x1 x2) outputs(y1 y2) [options]
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
            local _io_outputs = strtrim(substr(`"`0'"', 1, `_cpos'-1))
            local _io_rest    = substr(`"`0'"', `_cpos'+1, .)
            local 0 `", inputs(`_io_inputs') outputs(`_io_outputs') `_io_rest'"'
        }
        else {
            local 0 `", inputs(`_io_inputs') outputs(`0')"'
        }
    }

    syntax [if] [in], ///
        Inputs(varlist) ///
        Outputs(varlist) ///
        [DMU(varname) ///
         STRAtegy(string) ///
         MAXiter(integer 100) ///
         GENerate(string) ///
         TOLerance(real 1e-6) ///
         MATrix ///
         DETail ///
         DEBUG]
    
    // Defaults
    if "`strategy'" == "" local strategy "all"
    if "`generate'" == "" local generate "_cg"
    
    // Validate strategy
    if !inlist("`strategy'", "aggressive", "benevolent", "nash", "all") {
        di as error "strategy() must be aggressive, benevolent, nash, or all"
        exit 198
    }
    
    // Mark sample — do NOT include dmu in markout (string DMU causes r(109))
    tempvar touse
    mark `touse'
    markout `touse' `inputs' `outputs'

    // Count
    qui count if `touse'
    local n = r(N)

    // Display header
    di as text ""
    di as text "{hline 70}"
    di as text "Cross-Efficiency Game DEA"
    di as text "{hline 70}"
    di as text "Strategy:      " as result "`strategy'"
    di as text "DMUs:          " as result "`n'"
    di as text "Max iterations:" as result "`maxiter'" as text " (for Nash)"
    di as text "{hline 70}"
    
    // Create DMU identifier if not provided
    if "`dmu'" == "" {
        tempvar dmu_temp
        qui gen `dmu_temp' = _n if `touse'
        local dmu "`dmu_temp'"
    }
    
    // Create output variables
    if "`strategy'" == "all" | "`strategy'" == "aggressive" {
        cap drop `generate'_agg
        qui gen double `generate'_agg = . if `touse'
    }
    if "`strategy'" == "all" | "`strategy'" == "benevolent" {
        cap drop `generate'_ben
        qui gen double `generate'_ben = . if `touse'
    }
    if "`strategy'" == "all" | "`strategy'" == "nash" {
        cap drop `generate'_nash
        qui gen double `generate'_nash = . if `touse'
    }
    
    // Call Mata function
    mata: _dea_crossgame_compute("`inputs'", "`outputs'", ///
        "`dmu'", "`touse'", "`generate'", ///
        "`strategy'", `maxiter', `tolerance', ///
        "`matrix'" != "", "`detail'" != "", "`debug'" != "")
    
    di as text ""
    di as text "{hline 70}"
    
    // Return values
    return local strategy "`strategy'"
    return scalar n = `n'
end


// ============================================================================
// MATA FUNCTIONS
// ============================================================================

mata:
mata set matastrict on

// ─────────────────────────────────────────────────────────────────────────────
// _dea_dmu_read() — handles string and numeric DMU variables
// Returns: dmu_ids (numeric), dmu_labels (string for display)
// ─────────────────────────────────────────────────────────────────────────────
void _dea_cg_dmu_read(string scalar  dmu_var,
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
 * _dea_crossgame_compute()
 *
 * Main computation function for cross-efficiency game.
 */
void _dea_crossgame_compute(string scalar inputs,
                             string scalar outputs,
                             string scalar dmu_var,
                             string scalar touse,
                             string scalar gen_prefix,
                             string scalar strategy,
                             real scalar maxiter,
                             real scalar tol,
                             real scalar show_matrix,
                             real scalar show_detail,
                             real scalar show_debug)
{
    // Load data
    real matrix X, Y
    real colvector dmu_ids, obs_indices
    string colvector dmu_labels

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    _dea_cg_dmu_read(dmu_var, touse, dmu_ids, dmu_labels)
    obs_indices = selectindex(st_data(., touse))
    
    real scalar n, m, s
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    // Step 1: Compute self-appraisal (standard DEA efficiency)
    real colvector self_eff
    real matrix V_self, U_self  // Optimal weights
    self_eff = J(n, 1, .)
    V_self = J(n, m, .)
    U_self = J(n, s, .)
    
    real scalar i
    real rowvector v_i, u_i
    real scalar eff_i
    
    printf("\nStep 1: Computing self-appraisal scores...\n")
    
    for (i = 1; i <= n; i++) {
        _dea_cg_multiplier(X[i, .], Y[i, .], X, Y, tol, v_i, u_i, eff_i)
        self_eff[i] = eff_i
        V_self[i, .] = v_i
        U_self[i, .] = u_i
    }
    
    // Compute and display results based on strategy
    real colvector ce_agg, ce_ben, ce_nash
    real matrix CE_agg, CE_ben, CE_nash
    
    if (strategy == "all" | strategy == "aggressive") {
        printf("\nStep 2a: Computing aggressive cross-efficiency...\n")
        _dea_cg_aggressive(X, Y, self_eff, tol, ce_agg, CE_agg, show_debug)
        
        // Store results
        for (i = 1; i <= n; i++) {
            st_store(obs_indices[i], st_varindex(gen_prefix + "_agg"), ce_agg[i])
        }
    }
    
    if (strategy == "all" | strategy == "benevolent") {
        printf("\nStep 2b: Computing benevolent cross-efficiency...\n")
        _dea_cg_benevolent(X, Y, self_eff, tol, ce_ben, CE_ben, show_debug)
        
        // Store results
        for (i = 1; i <= n; i++) {
            st_store(obs_indices[i], st_varindex(gen_prefix + "_ben"), ce_ben[i])
        }
    }
    
    if (strategy == "all" | strategy == "nash") {
        printf("\nStep 2c: Computing Nash equilibrium cross-efficiency...\n")
        _dea_cg_nash(X, Y, self_eff, maxiter, tol, ce_nash, CE_nash, show_debug)
        
        // Store results
        for (i = 1; i <= n; i++) {
            st_store(obs_indices[i], st_varindex(gen_prefix + "_nash"), ce_nash[i])
        }
    }
    
    // Display results
    printf("\n")
    printf("Cross-Efficiency Game Results\n")
    printf("{hline 70}\n")
    
    if (strategy == "all") {
        printf("%10s | %10s | %10s %10s %10s\n", 
               "DMU", "Self", "Aggressive", "Benevolent", "Nash")
        printf("%10s | %10s | %10s %10s %10s\n", 
               "----------", "----------", "----------", "----------", "----------")
        
        for (i = 1; i <= n; i++) {
            printf("%10s | %10.4f | %10.4f %10.4f %10.4f\n",
                   dmu_labels[i], self_eff[i], ce_agg[i], ce_ben[i], ce_nash[i])
        }

        printf("{hline 70}\n")
        printf("%10s | %10.4f | %10.4f %10.4f %10.4f\n",
               "Mean", mean(self_eff), mean(ce_agg), mean(ce_ben), mean(ce_nash))
    }
    else {
        real colvector ce_result
        string scalar strat_name
        
        if (strategy == "aggressive") {
            ce_result = ce_agg
            strat_name = "Aggressive"
        }
        else if (strategy == "benevolent") {
            ce_result = ce_ben
            strat_name = "Benevolent"
        }
        else {
            ce_result = ce_nash
            strat_name = "Nash"
        }
        
        printf("%10s | %10s | %10s\n", "DMU", "Self", strat_name)
        printf("%10s | %10s | %10s\n", "----------", "----------", "----------")
        
        for (i = 1; i <= n; i++) {
            printf("%10s | %10.4f | %10.4f\n",
                   dmu_labels[i], self_eff[i], ce_result[i])
        }
        
        printf("{hline 50}\n")
        printf("%10s | %10.4f | %10.4f\n", "Mean", mean(self_eff), mean(ce_result))
    }
    printf("{hline 70}\n")
    
    // Display interpretation
    printf("\nInterpretation:\n")
    printf("  Self:       Standard DEA efficiency (self-appraisal)\n")
    printf("  Aggressive: Competitive strategy (minimize others' scores)\n")
    printf("  Benevolent: Cooperative strategy (maximize others' scores)\n")
    printf("  Nash:       Equilibrium (no unilateral improvement possible)\n")
}


/**
 * _dea_cg_multiplier()
 *
 * Solves multiplier DEA model to get optimal weights.
 */
void _dea_cg_multiplier(real rowvector x0,
                         real rowvector y0,
                         real matrix X,
                         real matrix Y,
                         real scalar tol,
                         real rowvector v_opt,
                         real rowvector u_opt,
                         real scalar eff)
{
    real scalar n, m, s
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    real rowvector c, lb, ub
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    
    c = (J(1, m, 0), -y0)
    A_eq = (x0, J(1, s, 0))
    b_eq = 1
    A_ineq = (-X, Y)
    b_ineq = J(n, 1, 0)
    
    real scalar eps
    eps = 1e-6
    lb = J(1, m + s, eps)
    ub = J(1, m + s, .)
    
    class LinearProgram scalar S
    S = LinearProgram()
    S.setMaxOrMin("min")
    S.setCoefficients(c)
    S.setEquality(A_eq, b_eq)
    S.setInequality(A_ineq, b_ineq)
    S.setBounds(lb, ub)
    S.setTrace("off")
    
    (void) S.optimize()
    
    if (S.converged() && S.errorcode() == 0) {
        real rowvector sol
        sol = S.parameters()
        v_opt = sol[1..m]
        u_opt = sol[(m+1)..(m+s)]
        eff = -S.value()
    }
    else {
        v_opt = J(1, m, .)
        u_opt = J(1, s, .)
        eff = .
    }
}


/**
 * _dea_cg_aggressive()
 *
 * Aggressive strategy: minimize average cross-efficiency of others
 * while maintaining own efficiency.
 */
void _dea_cg_aggressive(real matrix X,
                         real matrix Y,
                         real colvector self_eff,
                         real scalar tol,
                         real colvector ce_scores,
                         real matrix CE,
                         real scalar show_debug)
{
    real scalar n, m, s, i, j
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    // Optimal weights under aggressive strategy
    real matrix V_agg, U_agg
    V_agg = J(n, m, .)
    U_agg = J(n, s, .)
    
    for (i = 1; i <= n; i++) {
        real rowvector x0, y0
        x0 = X[i, .]
        y0 = Y[i, .]
        
        // Secondary goal: minimize sum of cross-efficiencies
        // min Σ_{j≠i} u*y_j (with v*x normalized to 1)
        // Subject to: own efficiency = self_eff[i]
        
        real rowvector c, lb, ub
        real matrix A_ineq, A_eq
        real colvector b_ineq, b_eq
        
        // Objective: minimize cross-efficiency to others
        // c = (0...0 for v, Σ_{j≠i} y_j for u)
        real rowvector sum_y_others
        sum_y_others = colsum(Y) - Y[i, .]
        c = (J(1, m, 0), sum_y_others)
        
        // Constraints
        // v*x0 = 1 (normalization)
        A_eq = (x0, J(1, s, 0))
        b_eq = 1
        
        // u*y0 = self_eff[i] (maintain own efficiency)
        A_eq = A_eq \ (J(1, m, 0), y0)
        b_eq = b_eq \ self_eff[i]
        
        // u*y_j - v*x_j <= 0 for all j (efficiency <= 1)
        A_ineq = (-X, Y)
        b_ineq = J(n, 1, 0)
        
        real scalar eps
        eps = 1e-6
        lb = J(1, m + s, eps)
        ub = J(1, m + s, .)
        
        class LinearProgram scalar S
        S = LinearProgram()
        S.setMaxOrMin("min")
        S.setCoefficients(c)
        S.setEquality(A_eq, b_eq)
        S.setInequality(A_ineq, b_ineq)
        S.setBounds(lb, ub)
        S.setTrace("off")
        
        (void) S.optimize()
        
        if (S.converged() && S.errorcode() == 0) {
            real rowvector sol
            sol = S.parameters()
            V_agg[i, .] = sol[1..m]
            U_agg[i, .] = sol[(m+1)..(m+s)]
        }
        else {
            // Fallback to self-appraisal weights
            real rowvector v_i, u_i
            real scalar eff_i
            _dea_cg_multiplier(x0, y0, X, Y, tol, v_i, u_i, eff_i)
            V_agg[i, .] = v_i
            U_agg[i, .] = u_i
        }
    }
    
    // Compute cross-efficiency matrix
    CE = J(n, n, .)
    for (i = 1; i <= n; i++) {
        for (j = 1; j <= n; j++) {
            real scalar num, denom
            num = U_agg[i, .] * Y[j, .]'
            denom = V_agg[i, .] * X[j, .]'
            if (denom > tol) {
                CE[i, j] = num / denom
            }
        }
    }
    
    // Compute cross-efficiency scores (average excluding self)
    ce_scores = J(n, 1, .)
    for (j = 1; j <= n; j++) {
        real scalar sum_ce, count_ce
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
}


/**
 * _dea_cg_benevolent()
 *
 * Benevolent strategy: maximize average cross-efficiency of others.
 */
void _dea_cg_benevolent(real matrix X,
                         real matrix Y,
                         real colvector self_eff,
                         real scalar tol,
                         real colvector ce_scores,
                         real matrix CE,
                         real scalar show_debug)
{
    real scalar n, m, s, i, j
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    real matrix V_ben, U_ben
    V_ben = J(n, m, .)
    U_ben = J(n, s, .)
    
    for (i = 1; i <= n; i++) {
        real rowvector x0, y0
        x0 = X[i, .]
        y0 = Y[i, .]
        
        // Objective: maximize cross-efficiency to others = min negative
        // c = (0...0 for v, -Σ_{j≠i} y_j for u)
        real rowvector sum_y_others
        sum_y_others = colsum(Y) - Y[i, .]
        real rowvector c
        c = (J(1, m, 0), -sum_y_others)
        
        // Constraints
        real matrix A_eq, A_ineq
        real colvector b_eq, b_ineq
        
        A_eq = (x0, J(1, s, 0))
        b_eq = 1
        
        A_eq = A_eq \ (J(1, m, 0), y0)
        b_eq = b_eq \ self_eff[i]
        
        A_ineq = (-X, Y)
        b_ineq = J(n, 1, 0)
        
        real scalar eps
        eps = 1e-6
        real rowvector lb, ub
        lb = J(1, m + s, eps)
        ub = J(1, m + s, .)
        
        class LinearProgram scalar S
        S = LinearProgram()
        S.setMaxOrMin("min")
        S.setCoefficients(c)
        S.setEquality(A_eq, b_eq)
        S.setInequality(A_ineq, b_ineq)
        S.setBounds(lb, ub)
        S.setTrace("off")
        
        (void) S.optimize()
        
        if (S.converged() && S.errorcode() == 0) {
            real rowvector sol
            sol = S.parameters()
            V_ben[i, .] = sol[1..m]
            U_ben[i, .] = sol[(m+1)..(m+s)]
        }
        else {
            real rowvector v_i, u_i
            real scalar eff_i
            _dea_cg_multiplier(x0, y0, X, Y, tol, v_i, u_i, eff_i)
            V_ben[i, .] = v_i
            U_ben[i, .] = u_i
        }
    }
    
    // Compute cross-efficiency matrix
    CE = J(n, n, .)
    for (i = 1; i <= n; i++) {
        for (j = 1; j <= n; j++) {
            real scalar num, denom
            num = U_ben[i, .] * Y[j, .]'
            denom = V_ben[i, .] * X[j, .]'
            if (denom > tol) {
                CE[i, j] = num / denom
            }
        }
    }
    
    // Compute cross-efficiency scores
    ce_scores = J(n, 1, .)
    for (j = 1; j <= n; j++) {
        real scalar sum_ce, count_ce
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
}


/**
 * _dea_cg_nash()
 *
 * Nash equilibrium: iteratively update until no DMU can improve.
 */
void _dea_cg_nash(real matrix X,
                   real matrix Y,
                   real colvector self_eff,
                   real scalar maxiter,
                   real scalar tol,
                   real colvector ce_scores,
                   real matrix CE,
                   real scalar show_debug)
{
    real scalar n, m, s, d, i, j, iter
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    // Liang et al. (2008), model (5): initialize theta with ordinary
    // cross-efficiency scores and iterate best responses that preserve each
    // DMU d's current score theta[d].
    real matrix V_self, U_self, CE0, CE_new
    real colvector theta_prev, theta_new
    V_self = J(n, m, .)
    U_self = J(n, s, .)
    
    for (i = 1; i <= n; i++) {
        real rowvector v_i, u_i
        real scalar eff_i
        _dea_cg_multiplier(X[i, .], Y[i, .], X, Y, tol, v_i, u_i, eff_i)
        V_self[i, .] = v_i
        U_self[i, .] = u_i
    }

    CE0 = J(n, n, .)
    for (i = 1; i <= n; i++) {
        for (j = 1; j <= n; j++) {
            real scalar num0, den0
            num0 = U_self[i, .] * Y[j, .]'
            den0 = V_self[i, .] * X[j, .]'
            if (den0 > tol) CE0[i, j] = num0 / den0
        }
    }

    theta_prev = J(n, 1, .)
    for (j = 1; j <= n; j++) {
        theta_prev[j] = mean(CE0[., j])
    }

    real scalar converged, max_change, status
    converged = 0
    
    for (iter = 1; iter <= maxiter; iter++) {
        CE_new = J(n, n, .)

        for (d = 1; d <= n; d++) {
            for (j = 1; j <= n; j++) {
                real rowvector c, lb, ub
                real matrix A_ineq, A_eq
                real colvector b_ineq, b_eq
                class LinearProgram scalar S

                // Maximize DMU j's efficiency, implemented as minimize -u*y_j.
                c = (J(1, m, 0), -Y[j, .])

                // Normalization: v*x_j = 1.
                A_eq = (X[j, .], J(1, s, 0))
                b_eq = 1

                // Frontier: u*y_l - v*x_l <= 0 for all l.
                A_ineq = (-X, Y)
                b_ineq = J(n, 1, 0)

                // Game constraint: theta_d * v*x_d - u*y_d <= 0.
                A_ineq = A_ineq \ (theta_prev[d] * X[d, .], -Y[d, .])
                b_ineq = b_ineq \ 0

                lb = J(1, m + s, 0)
                ub = J(1, m + s, .)

                S = LinearProgram()
                S.setMaxOrMin("min")
                S.setCoefficients(c)
                S.setEquality(A_eq, b_eq)
                S.setInequality(A_ineq, b_ineq)
                S.setBounds(lb, ub)
                S.setTrace("off")

                status = S.optimize()

                if (S.converged() && S.errorcode() == 0) {
                    CE_new[d, j] = -S.value()
                    if (CE_new[d, j] < 0 & abs(CE_new[d, j]) < tol) {
                        CE_new[d, j] = 0
                    }
                }
            }
        }

        theta_new = J(n, 1, .)
        for (j = 1; j <= n; j++) theta_new[j] = mean(CE_new[., j])

        max_change = max(abs(theta_new :- theta_prev))

        if (show_debug) {
            printf("  Iteration %g: max change = %g\n", iter, max_change)
        }

        if (max_change < tol) {
            converged = 1
            break
        }

        theta_prev = theta_new
    }
    
    if (converged) {
        printf("  Game cross-efficiency converged in %g iterations.\n", iter)
    }
    else {
        printf("  Warning: game cross-efficiency did not converge in %g iterations.\n", maxiter)
    }

    CE = CE_new
    ce_scores = theta_new
}

end
