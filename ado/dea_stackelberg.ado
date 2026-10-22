*! version 1.1.1  02may2026
*! dea_stackelberg.ado — Stackelberg Game DEA
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
*!   Liang, L., Yang, F., Cook, W. D., & Zhu, J. (2006). DEA models for
*!     supply chain efficiency evaluation. Annals of Operations Research,
*!     145(1), 35-49. doi:10.1007/s10479-006-0026-7
*!
*!   Li, Y., Chen, Y., Liang, L., & Xie, J. (2012). DEA models for
*!     extended two-stage network structures. Omega, 40(5), 611-618.
*!     doi:10.1016/j.omega.2011.11.007
*!
*! Changelog:
*!   1.1.1  02may2026  Tighten help alignment and support string DMU identifiers
*!   1.1.0  28apr2026  Remove spurious [varlist] from syntax

/**
 * STACKELBERG GAME DEA
 * ====================
 * 
 * Computes a leader-follower-style two-stage DEA decomposition.
 * 
 * Applications:
 *   - Supply chain: manufacturer (leader) vs retailer (follower)
 *   - Two-stage production: Stage 1 (leader) vs Stage 2 (follower)
 *   - Principal-agent relationships
 * 
 * Approach:
 *   1. Compute independent Stage 1 and Stage 2 efficiencies
 *   2. Report the selected leader stage and follower stage efficiencies
 *   3. Form an overall score as the product of leader and follower scores
 * 
 * References:
 *   - Liang et al. (2006): DEA models for supply chain efficiency evaluation
 *   - Li et al. (2012): DEA models for extended two-stage network structures
 */

capture program drop dea_stackelberg
program define dea_stackelberg, rclass
    version 16.0
    
    syntax [if] [in], ///
        Inputs(varlist numeric) ///
        INTERmediates(varlist numeric) ///
        Outputs(varlist numeric) ///
        [DMU(varname) ///
         LEADer(string) ///
         RTS(string) ///
         GENerate(string) ///
         TOLerance(real 1e-7) ///
         DETail ///
         DEBUG]
    
    // Defaults
    if "`leader'" == "" local leader "stage1"
    local leader = lower("`leader'")
    if "`rts'" == "" local rts "crs"
    local rts = lower("`rts'")
    if "`generate'" == "" local generate "_stack"
    
    // Validate
    if !inlist("`rts'", "crs", "vrs") {
        di as error "rts() must be crs or vrs"
        exit 198
    }
    if !inlist("`leader'", "stage1", "stage2") {
        di as error "leader() must be stage1 or stage2"
        exit 198
    }
    
    // Mark sample
    tempvar touse
    marksample touse
    markout `touse' `inputs' `intermediates' `outputs'
    if "`dmu'" != "" {
        qui replace `touse' = 0 if missing(`dmu')
    }
    
    // Count
    qui count if `touse'
    local n = r(N)
    if `n' == 0 {
        di as error "no observations"
        exit 2000
    }
    
    // Display header
    di as text ""
    di as text "{hline 70}"
    di as text "Stackelberg Game DEA (Two-Stage)"
    di as text "{hline 70}"
    di as text "Leader:        " as result "Stage `=substr("`leader'", 6, 1)'"
    di as text "Follower:      " as result "Stage `=cond("`leader'"=="stage1", "2", "1")'"
    di as text "RTS:           " as result "`rts'"
    di as text "DMUs:          " as result "`n'"
    di as text "{hline 70}"
    
    // Create output variables
    cap drop `generate'_leader
    cap drop `generate'_follower
    cap drop `generate'_overall
    cap drop `generate'_indep1
    cap drop `generate'_indep2
    qui gen double `generate'_leader = . if `touse'
    qui gen double `generate'_follower = . if `touse'
    qui gen double `generate'_overall = . if `touse'
    qui gen double `generate'_indep1 = . if `touse'
    qui gen double `generate'_indep2 = . if `touse'
    
    // Create DMU identifier if not provided
    if "`dmu'" == "" {
        tempvar dmu_temp
        qui gen `dmu_temp' = _n if `touse'
        local dmu "`dmu_temp'"
    }
    
    // Call Mata function
    mata: _dea_stackelberg_compute("`inputs'", "`intermediates'", "`outputs'", ///
        "`dmu'", "`touse'", "`generate'", ///
        "`leader'", "`rts'", `tolerance', ///
        "`detail'" != "", "`debug'" != "")
    
    di as text ""
    di as text "{hline 70}"
    
    // Return values
    return local leader "`leader'"
    return local rts "`rts'"
    return scalar n = `n'
end


// ============================================================================
// MATA FUNCTIONS
// ============================================================================

capture mata: mata drop _dea_stack_dmu_read()
capture mata: mata drop _dea_stackelberg_compute()
capture mata: mata drop _dea_stack_stage_eff()
capture mata: mata drop _dea_stack_s1_leads()
capture mata: mata drop _dea_stack_s2_leads()

mata:
mata set matastrict on

void _dea_stack_dmu_read(string scalar dmu_var,
                          string scalar touse,
                          string colvector dmu_labels)
{
    real colvector dmu_num

    if (substr(st_vartype(dmu_var), 1, 3) == "str") {
        dmu_labels = st_sdata(., dmu_var, touse)
    }
    else {
        dmu_num = st_data(., dmu_var, touse)
        dmu_labels = strofreal(dmu_num)
    }
}

void _dea_stackelberg_compute(string scalar inputs,
                               string scalar intermediates,
                               string scalar outputs,
                               string scalar dmu_var,
                               string scalar touse,
                               string scalar gen_prefix,
                               string scalar leader,
                               string scalar rts,
                               real scalar tol,
                               real scalar show_detail,
                               real scalar show_debug)
{
    // Load data
    real matrix X, Z, Y
    real colvector obs_indices
    string colvector dmu_labels
    
    X = st_data(., tokens(inputs), touse)
    Z = st_data(., tokens(intermediates), touse)
    Y = st_data(., tokens(outputs), touse)
    _dea_stack_dmu_read(dmu_var, touse, dmu_labels)
    obs_indices = selectindex(st_data(., touse))
    
    real scalar n, m, d, s, i
    n = rows(X)
    m = cols(X)
    d = cols(Z)
    s = cols(Y)
    
    // Step 1: Compute independent efficiencies (benchmark)
    printf("\nStep 1: Computing independent stage efficiencies...\n")
    real colvector eff_indep1, eff_indep2
    eff_indep1 = J(n, 1, .)
    eff_indep2 = J(n, 1, .)
    
    for (i = 1; i <= n; i++) {
        eff_indep1[i] = _dea_stack_stage_eff(X[i, .], Z[i, .], X, Z, rts, tol)
        eff_indep2[i] = _dea_stack_stage_eff(Z[i, .], Y[i, .], Z, Y, rts, tol)
    }
    
    // Step 2: Stackelberg game
    printf("Step 2: Computing leader-follower ordering...\n")
    real colvector eff_leader, eff_follower, eff_overall
    eff_leader = J(n, 1, .)
    eff_follower = J(n, 1, .)
    eff_overall = J(n, 1, .)
    
    real rowvector result
    for (i = 1; i <= n; i++) {
        if (leader == "stage1") {
            // Stage 1 is reported as leader; Stage 2 as follower.
            result = _dea_stack_s1_leads(X[i, .], Z[i, .], Y[i, .], 
                                          X, Z, Y, rts, tol)
        }
        else {
            // Stage 2 is reported as leader; Stage 1 as follower.
            result = _dea_stack_s2_leads(X[i, .], Z[i, .], Y[i, .], 
                                          X, Z, Y, rts, tol)
        }
        
        eff_leader[i] = result[1]
        eff_follower[i] = result[2]
        
        // Overall efficiency = product
        if (eff_leader[i] != . && eff_follower[i] != .) {
            eff_overall[i] = eff_leader[i] * eff_follower[i]
        }
    }
    
    // Display results
    printf("\n")
    printf("Stackelberg DEA Results\n")
    printf("{hline 80}\n")
    printf("%10s | %10s %10s | %10s %10s | %10s\n", 
           "DMU", "Indep S1", "Indep S2", "Leader", "Follower", "Overall")
    printf("%10s | %10s %10s | %10s %10s | %10s\n", 
           "----------", "----------", "----------", "----------", "----------", "----------")
    
    for (i = 1; i <= n; i++) {
        printf("%10s | %10.4f %10.4f | %10.4f %10.4f | %10.4f\n",
               dmu_labels[i], eff_indep1[i], eff_indep2[i], 
               eff_leader[i], eff_follower[i], eff_overall[i])
    }
    
    printf("{hline 80}\n")
    printf("%10s | %10.4f %10.4f | %10.4f %10.4f | %10.4f\n",
           "Mean", mean(eff_indep1), mean(eff_indep2), 
           mean(eff_leader), mean(eff_follower), mean(eff_overall))
    printf("{hline 80}\n")
    
    printf("\nInterpretation:\n")
    if (leader == "stage1") {
        printf("  Stage 1 is reported as leader; Stage 2 as follower\n")
    }
    else {
        printf("  Stage 2 is reported as leader; Stage 1 as follower\n")
    }
    printf("  Overall = Leader × Follower\n")
    
    // Store results
    for (i = 1; i <= n; i++) {
        st_store(obs_indices[i], st_varindex(gen_prefix + "_leader"), eff_leader[i])
        st_store(obs_indices[i], st_varindex(gen_prefix + "_follower"), eff_follower[i])
        st_store(obs_indices[i], st_varindex(gen_prefix + "_overall"), eff_overall[i])
        st_store(obs_indices[i], st_varindex(gen_prefix + "_indep1"), eff_indep1[i])
        st_store(obs_indices[i], st_varindex(gen_prefix + "_indep2"), eff_indep2[i])
    }
}


real scalar _dea_stack_stage_eff(real rowvector x0,
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
 * _dea_stack_s1_leads()
 *
 * Stage 1 as leader: report Stage 1 efficiency first, then Stage 2.
 * Returns: (eff_leader, eff_follower)
 */
real rowvector _dea_stack_s1_leads(real rowvector x0,
                                    real rowvector z0,
                                    real rowvector y0,
                                    real matrix X,
                                    real matrix Z,
                                    real matrix Y,
                                    string scalar rts,
                                    real scalar tol)
{
    real scalar eff_leader, eff_follower
    
    eff_leader = _dea_stack_stage_eff(x0, z0, X, Z, rts, tol)
    eff_follower = _dea_stack_stage_eff(z0, y0, Z, Y, rts, tol)
    
    return((eff_leader, eff_follower))
}


/**
 * _dea_stack_s2_leads()
 *
 * Stage 2 as leader: report Stage 2 efficiency first, then Stage 1.
 * Returns: (eff_leader, eff_follower)
 */
real rowvector _dea_stack_s2_leads(real rowvector x0,
                                    real rowvector z0,
                                    real rowvector y0,
                                    real matrix X,
                                    real matrix Z,
                                    real matrix Y,
                                    string scalar rts,
                                    real scalar tol)
{
    real scalar eff_leader, eff_follower
    
    eff_leader = _dea_stack_stage_eff(z0, y0, Z, Y, rts, tol)
    eff_follower = _dea_stack_stage_eff(x0, z0, X, Z, rts, tol)
    
    return((eff_leader, eff_follower))
}

end
