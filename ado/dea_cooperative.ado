*! version 1.1.1  02may2026
*! dea_cooperative.ado — Cooperative Game DEA (Shapley Value)
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
*!   Nakabayashi, K., & Tone, K. (2006). Egoist's dilemma: A DEA game.
*!     Omega, 34(2), 135-148. doi:10.1016/j.omega.2004.08.003
*!
*!   Lozano, S. (2012). Information sharing in DEA: A cooperative game
*!     theory approach. European Journal of Operational Research,
*!     222(3), 558-565. doi:10.1016/j.ejor.2012.05.014
*!
*!   Li, Y., & Liang, L. (2010). A Shapley value index on the importance
*!     of variables in DEA models. Expert Systems with Applications,
*!     37(9), 6287-6292. doi:10.1016/j.eswa.2010.02.093
*!
*! Changelog:
*!   1.1.1  02may2026  Return r(v_grand) and display coalition values with detail
*!   1.1.0  28apr2026  Add Form A syntax (x1 x2 = y1 y2)

/**
 * COOPERATIVE GAME DEA
 * ====================
 * 
 * Uses cooperative game theory to fairly allocate efficiency gains
 * among DMUs that form coalitions.
 * 
 * Key concepts:
 *   - Characteristic function v(S): Efficiency of coalition S
 *   - Shapley value: Fair allocation based on marginal contributions
 *   - Core: Stable allocations where no coalition wants to deviate
 * 
 * Shapley value formula:
 *   φ_i = Σ_{S⊆N\{i}} [|S|!(n-|S|-1)!/n!] × [v(S∪{i}) - v(S)]
 * 
 * References:
 *   - Nakabayashi & Tone (2006): Egoist's dilemma as a DEA game
 *   - Lozano (2012): Information sharing in DEA as a cooperative game
 *   - Li & Liang (2010): Shapley value index for DEA variables
 */

capture program drop dea_cooperative
program define dea_cooperative, rclass
    version 16.0

    // =========================================================================
    // Unified Syntax Parser (Form A / Form B)
    //   Form A: dea_cooperative x1 x2 = y1 y2 [if] [in], [options]
    //   Form B: dea_cooperative [if] [in], inputs(x1 x2) outputs(y1 y2) [options]
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
         RTS(string) ///
         GENerate(string) ///
         TOLerance(real 1e-7) ///
         DETail ///
         DEBUG]
    
    // Defaults
    if "`rts'" == "" local rts "crs"
    if "`generate'" == "" local generate "_coop"
    
    // Validate
    if !inlist("`rts'", "crs", "vrs") {
        di as error "rts() must be crs or vrs"
        exit 198
    }
    
    // Mark sample — do NOT include dmu in markout (string DMU causes r(109))
    tempvar touse
    mark `touse'
    markout `touse' `inputs' `outputs'

    // Count
    qui count if `touse'
    local n = r(N)
    
    // Warning for large n (Shapley is O(2^n))
    if `n' > 15 {
        di as error "Warning: n=`n' is large. Shapley value computation is O(2^n)."
        di as error "Consider using a subset or approximation."
    }
    if `n' > 20 {
        di as error "n=`n' is too large for exact Shapley computation."
        exit 198
    }
    
    // Display header
    di as text ""
    di as text "{hline 70}"
    di as text "Cooperative Game DEA (Shapley Value)"
    di as text "{hline 70}"
    di as text "RTS:           " as result "`rts'"
    di as text "DMUs:          " as result "`n'"
    di as text "Coalitions:    " as result "2^`n' = `=2^`n''"
    di as text "{hline 70}"
    
    // Create output variables
    cap drop `generate'_indiv
    cap drop `generate'_shapley
    cap drop `generate'_contrib
    cap drop `generate'_rank
    qui gen double `generate'_indiv = . if `touse'
    qui gen double `generate'_shapley = . if `touse'
    qui gen double `generate'_contrib = . if `touse'
    qui gen double `generate'_rank = . if `touse'
    
    // Create DMU identifier if not provided
    if "`dmu'" == "" {
        tempvar dmu_temp
        qui gen `dmu_temp' = _n if `touse'
        local dmu "`dmu_temp'"
    }
    
    // Call Mata function
    mata: _dea_cooperative_compute("`inputs'", "`outputs'", ///
        "`dmu'", "`touse'", "`generate'", ///
        "`rts'", `tolerance', ///
        "`detail'" != "", "`debug'" != "")
    tempname v_grand
    scalar `v_grand' = __dea_coop_v_grand
    capture scalar drop __dea_coop_v_grand
    
    di as text ""
    di as text "{hline 70}"
    
    // Return values
    return local rts "`rts'"
    return scalar n = `n'
    return scalar v_grand = `v_grand'
end


// ============================================================================
// MATA FUNCTIONS
// ============================================================================

mata:
mata set matastrict on

// ─────────────────────────────────────────────────────────────────────────────
// _dea_coop_dmu_read() — handles string and numeric DMU variables
// ─────────────────────────────────────────────────────────────────────────────
void _dea_coop_dmu_read(string scalar  dmu_var,
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


void _dea_cooperative_compute(string scalar inputs,
                               string scalar outputs,
                               string scalar dmu_var,
                               string scalar touse,
                               string scalar gen_prefix,
                               string scalar rts,
                               real scalar tol,
                               real scalar show_detail,
                               real scalar show_debug)
{
    // Load data
    real matrix X, Y
    real colvector dmu_ids, obs_indices
    string colvector dmu_labels

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    _dea_coop_dmu_read(dmu_var, touse, dmu_ids, dmu_labels)
    obs_indices = selectindex(st_data(., touse))
    
    real scalar n, m, s, i
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    // Step 1: Compute individual efficiencies
    printf("\nStep 1: Computing individual efficiencies...\n")
    real colvector eff_indiv
    eff_indiv = J(n, 1, .)
    
    for (i = 1; i <= n; i++) {
        eff_indiv[i] = _dea_coop_efficiency(X[i, .], Y[i, .], X, Y, rts, tol)
    }
    
    // Step 2: Compute characteristic function for all coalitions
    printf("Step 2: Computing coalition values (2^%g = %g coalitions)...\n", 
           n, 2^n)
    real colvector v  // v[k] = value of coalition k (binary representation)
    v = J(2^n, 1, 0)
    
    real scalar k, j, count
    real matrix X_coal, Y_coal
    real colvector members
    
    for (k = 1; k <= 2^n - 1; k++) {
        // Get coalition members from binary representation
        members = _dea_coop_get_members(k, n)
        
        if (rows(members) > 0) {
            // Coalition value = total efficiency of the coalition members.
            X_coal = X[members, .]
            Y_coal = Y[members, .]
            
            // Aggregate efficiency: can the coalition achieve more together?
            v[k+1] = _dea_coop_coalition_value(X_coal, Y_coal, X, Y, rts, tol)
        }
        
        if (show_debug && mod(k, 100) == 0) {
            printf("  Coalition %g/%g...\n", k, 2^n - 1)
        }
    }
    st_numscalar("__dea_coop_v_grand", v[2^n])

    if (show_detail) {
        printf("\nCoalition values\n")
        printf("{hline 40}\n")
        printf("%10s | %12s\n", "Coalition", "v(S)")
        printf("%10s | %12s\n", "----------", "------------")
        for (k = 1; k <= 2^n - 1; k++) {
            printf("%10.0f | %12.4f\n", k, v[k+1])
        }
        printf("{hline 40}\n")
    }
    
    // Step 3: Compute Shapley values
    printf("Step 3: Computing Shapley values...\n")
    real colvector shapley
    shapley = _dea_coop_shapley(v, n)
    
    // Step 4: Compute marginal contribution (relative to individual)
    real colvector contrib
    contrib = J(n, 1, .)
    for (i = 1; i <= n; i++) {
        contrib[i] = shapley[i] - eff_indiv[i]
    }
    
    // Rank by Shapley value
    real colvector ranks
    ranks = _dea_coop_rank(shapley)
    
    // Display results
    printf("\n")
    printf("Cooperative Game DEA Results (Shapley Value)\n")
    printf("{hline 70}\n")
    printf("%10s | %10s %12s %12s %6s\n", 
           "DMU", "Individual", "Shapley", "Contribution", "Rank")
    printf("%10s | %10s %12s %12s %6s\n", 
           "----------", "----------", "------------", "------------", "------")
    
    for (i = 1; i <= n; i++) {
        printf("%10s | %10.4f %12.4f %12.4f %6.0f\n",
               dmu_labels[i], eff_indiv[i], shapley[i], contrib[i], ranks[i])
    }
    
    printf("{hline 70}\n")
    printf("%10s | %10.4f %12.4f %12.4f\n",
           "Mean", mean(eff_indiv), mean(shapley), mean(contrib))
    printf("%10s | %10.4f %12.4f\n",
           "Sum", sum(eff_indiv), sum(shapley))
    printf("{hline 70}\n")
    
    printf("\nInterpretation:\n")
    printf("  Individual:   Standard DEA efficiency (alone)\n")
    printf("  Shapley:      Fair share of coalition value\n")
    printf("  Contribution: Shapley allocation minus individual efficiency\n")
    printf("  Rank:         Ranking by Shapley value (1 = highest)\n")
    
    // Store results
    for (i = 1; i <= n; i++) {
        st_store(obs_indices[i], st_varindex(gen_prefix + "_indiv"), eff_indiv[i])
        st_store(obs_indices[i], st_varindex(gen_prefix + "_shapley"), shapley[i])
        st_store(obs_indices[i], st_varindex(gen_prefix + "_contrib"), contrib[i])
        st_store(obs_indices[i], st_varindex(gen_prefix + "_rank"), ranks[i])
    }
}


/**
 * _dea_coop_efficiency()
 *
 * Standard DEA efficiency.
 */
real scalar _dea_coop_efficiency(real rowvector x0,
                                  real rowvector y0,
                                  real matrix X,
                                  real matrix Y,
                                  string scalar rts,
                                  real scalar tol)
{
    real scalar n, m, s
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    real rowvector c, lb, ub
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    
    c = (J(1, n, 0), 1)
    A_ineq = (X', -x0') \ (-Y', J(s, 1, 0))
    b_ineq = J(m, 1, 0) \ (-y0')
    lb = J(1, n + 1, 0)
    ub = J(1, n + 1, .)
    
    A_eq = J(0, n + 1, .)
    b_eq = J(0, 1, .)
    
    if (rts == "vrs") {
        A_eq = (J(1, n, 1), 0)
        b_eq = 1
    }
    
    class LinearProgram scalar S
    S = LinearProgram()
    S.setMaxOrMin("min")
    S.setCoefficients(c)
    S.setInequality(A_ineq, b_ineq)
    if (rows(A_eq) > 0) {
        S.setEquality(A_eq, b_eq)
    }
    S.setBounds(lb, ub)
    S.setTrace("off")
    
    (void) S.optimize()
    
    if (S.converged() && S.errorcode() == 0) {
        return(S.value())
    }
    return(.)
}


/**
 * _dea_coop_get_members()
 *
 * Get coalition members from binary representation.
 * k = 5 (binary 101) with n=3 means DMUs {1, 3}
 */
real colvector _dea_coop_get_members(real scalar k, real scalar n)
{
    real colvector members
    real scalar i, bit, count
    
    members = J(n, 1, 0)
    count = 0
    
    for (i = 1; i <= n; i++) {
        bit = mod(floor(k / 2^(i-1)), 2)
        if (bit == 1) {
            count++
            members[count] = i
        }
    }
    
    if (count > 0) {
        return(members[1..count])
    }
    return(J(0, 1, .))
}


/**
 * _dea_coop_coalition_value()
 *
 * Compute the value (total efficiency) of a coalition.
 * This is the total efficiency of coalition members under the reference frontier.
 */
real scalar _dea_coop_coalition_value(real matrix X_coal,
                                       real matrix Y_coal,
                                       real matrix X_all,
                                       real matrix Y_all,
                                       string scalar rts,
                                       real scalar tol)
{
    real scalar n_coal, i
    real scalar total_eff
    
    n_coal = rows(X_coal)
    total_eff = 0
    
    // Evaluate each coalition member against the full frontier
    for (i = 1; i <= n_coal; i++) {
        real scalar eff_i
        eff_i = _dea_coop_efficiency(X_coal[i, .], Y_coal[i, .], 
                                      X_all, Y_all, rts, tol)
        if (eff_i != .) {
            total_eff = total_eff + eff_i
        }
    }
    
    return(total_eff)
}


/**
 * _dea_coop_shapley()
 *
 * Compute Shapley values from characteristic function.
 * 
 * φ_i = Σ_{S⊆N\{i}} [|S|!(n-|S|-1)!/n!] × [v(S∪{i}) - v(S)]
 */
real colvector _dea_coop_shapley(real colvector v, real scalar n)
{
    real colvector shapley
    shapley = J(n, 1, 0)
    
    real scalar i, k, s_size
    real colvector members
    real scalar v_with, v_without, weight
    
    for (i = 1; i <= n; i++) {
        // Sum over all coalitions not containing i
        for (k = 0; k <= 2^n - 1; k++) {
            members = _dea_coop_get_members(k, n)
            
            // Check if i is not in coalition k
            if (!_dea_coop_in_set(i, members)) {
                s_size = rows(members)
                if (s_size == 0) s_size = 0
                
                // Coalition with i
                real scalar k_with_i
                k_with_i = k + 2^(i-1)
                
                v_without = v[k + 1]
                v_with = v[k_with_i + 1]
                
                // Shapley weight
                weight = (factorial(s_size) * factorial(n - s_size - 1)) / factorial(n)
                
                shapley[i] = shapley[i] + weight * (v_with - v_without)
            }
        }
    }
    
    return(shapley)
}


/**
 * _dea_coop_in_set()
 *
 * Check if element i is in the set.
 */
real scalar _dea_coop_in_set(real scalar i, real colvector members)
{
    real scalar j
    
    if (rows(members) == 0) return(0)
    
    for (j = 1; j <= rows(members); j++) {
        if (members[j] == i) return(1)
    }
    return(0)
}


/**
 * _dea_coop_rank()
 *
 * Rank by value (1 = highest).
 */
real colvector _dea_coop_rank(real colvector x)
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
