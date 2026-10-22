*! version 1.1.2  03may2026
*! dea_network.ado — Network DEA (Two-Stage)
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
*!   Fare, R., & Grosskopf, S. (2000). Network DEA. Socio-Economic
*!     Planning Sciences, 34(1), 35-49. doi:10.1016/S0038-0121(99)00012-9
*!
*!   Kao, C., & Hwang, S.-N. (2008). Efficiency decomposition in
*!     two-stage data envelopment analysis: An application to non-life
*!     insurance companies in Taiwan. European Journal of Operational
*!     Research, 185(1), 418-429. doi:10.1016/j.ejor.2006.11.041
*!
*!   Cooper, W. W., Seiford, L. M., & Tone, K. (2007). Data Envelopment
*!     Analysis: A Comprehensive Text with Models, Applications, References
*!     and DEA-Solver Software (2nd ed.). New York: Springer.
*!     doi:10.1007/978-0-387-45283-8
*!
*! Changelog:
*!   1.1.2  03may2026  Implement Kao-Hwang CRS/input shared intermediate-weight
*!                     relational multiplier model for method(multiplicative)
*!   1.1.1  02may2026  Fix string DMU handling and pass orientation to stage LPs
*!   1.1.0  11feb2026  Unified syntax: dea_network x = y | inputs()/outputs()
*!   1.0.0  06feb2026  Initial packaged version; based on earlier teaching/seminar versions

/**
 * NETWORK DEA (TWO-STAGE)
 * =======================
 * 
 * Network DEA models multi-stage production processes where
 * intermediate products connect stages.
 * 
 * Two-Stage Model:
 *   Stage 1: X → Z (inputs to intermediates)
 *   Stage 2: Z → Y (intermediates to outputs)
 * 
 * Overall efficiency can be decomposed:
 *   - Multiplicative: E_overall = E_stage1 × E_stage2
 *   - Additive: E_overall = w1*E_stage1 + w2*E_stage2
 * 
 * This implementation uses the relational model (Kao & Hwang, 2008)
 * which ensures consistency between stages.
 */

capture program drop dea_network
program define dea_network, rclass
    version 16.0
    
    // =========================================================================
    // Unified Syntax Parser
    //   Form A: dea_network x1 = y1 y2, intermediates(z) [options]
    //   Form B: dea_network [varlist], inputs() intermediates() outputs() [opts]
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
            local _pre = substr(`"`0'"', 1, `_cpos')
            local _post = substr(`"`0'"', `_cpos'+1, .)
            local 0 `"`_pre' inputs(`_io_inputs') `_post'"'
        }
        else {
            local 0 `"`0', inputs(`_io_inputs')"'
        }
    }
    
    syntax [varlist(numeric)] [if] [in], ///
        INTERmediates(varlist numeric) ///
        [Inputs(varlist numeric) ///
         Outputs(varlist numeric) ///
         DMU(string) ///
         RTS(string) ///
         ORT(string) ORIentation(string) ///
         METHOD(string) ///
         GENerate(string) ///
         TOLerance(real 1e-7) ///
         DETail ///
         DEBUG]
    
    // Resolve outputs
    if "`outputs'" == "" {
        if "`varlist'" != "" {
            local outputs `varlist'
        }
        else {
            di as error "outputs() required, or use: dea_network inputs = outputs, intermediates() [options]"
            exit 198
        }
    }
    if "`inputs'" == "" {
        di as error "inputs() required, or use: dea_network inputs = outputs, intermediates() [options]"
        exit 198
    }
    if "`dmu'" != "" {
        confirm variable `dmu'
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
    if "`method'" == "" local method "multiplicative"
    if "`generate'" == "" local generate "_net_eff"
    
    // Abbreviations for orientation
    if inlist("`orientation'", "i", "in")  local orientation "input"
    if inlist("`orientation'", "o", "out") local orientation "output"
    
    // Validate
    if !inlist("`rts'", "crs", "vrs") {
        di as error "rts() must be crs or vrs"
        exit 198
    }
    if !inlist("`orientation'", "input", "output") {
        di as error "orientation() must be input or output"
        exit 198
    }
    if !inlist("`method'", "multiplicative", "additive", "independent") {
        di as error "method() must be multiplicative, additive, or independent"
        exit 198
    }
    
    // Mark sample
    marksample touse, novarlist
    markout `touse' `inputs' `intermediates' `outputs'
    if "`dmu'" != "" {
        qui replace `touse' = 0 if missing(`dmu')
    }
    
    // Count
    qui count if `touse'
    local n = r(N)
    
    // Display header
    di as text ""
    di as text "{hline 70}"
    di as text "Network DEA (Two-Stage Model)"
    di as text "{hline 70}"
    di as text "Method:        " as result "`method'"
    di as text "RTS:           " as result "`rts'"
    di as text "DMUs:          " as result "`n'"
    di as text "Inputs:        " as result "`inputs'"
    di as text "Intermediates: " as result "`intermediates'"
    di as text "Outputs:       " as result "`outputs'"
    di as text "{hline 70}"
    
    // Create output variables
    cap drop `generate'_overall
    cap drop `generate'_stage1
    cap drop `generate'_stage2
    qui gen double `generate'_overall = . if `touse'
    qui gen double `generate'_stage1 = . if `touse'
    qui gen double `generate'_stage2 = . if `touse'
    
    // Create DMU identifier if not provided
    if "`dmu'" == "" {
        tempvar dmu_temp
        qui gen `dmu_temp' = _n if `touse'
        local dmu "`dmu_temp'"
    }
    
    // Call Mata function
    mata: _dea_network_compute("`inputs'", "`intermediates'", "`outputs'", ///
        "`dmu'", "`touse'", ///
        "`generate'_overall", "`generate'_stage1", "`generate'_stage2", ///
        "`method'", "`rts'", "`orientation'", `tolerance', ///
        "`detail'" != "", "`debug'" != "")
    
    // Summary statistics
    qui summ `generate'_overall if `touse'
    local mean_overall = r(mean)
    qui summ `generate'_stage1 if `touse'
    local mean_stage1 = r(mean)
    qui summ `generate'_stage2 if `touse'
    local mean_stage2 = r(mean)
    
    di as text ""
    di as text "{hline 70}"
    di as text "Summary"
    di as text "{hline 70}"
    di as text "Mean Stage 1 efficiency:  " as result %9.4f `mean_stage1'
    di as text "Mean Stage 2 efficiency:  " as result %9.4f `mean_stage2'
    di as text "Mean Overall efficiency:  " as result %9.4f `mean_overall'
    di as text "{hline 70}"
    
    // Return values
    return local method "`method'"
    return local rts "`rts'"
    return local orientation "`orientation'"
    return scalar n = `n'
    return scalar mean_overall = `mean_overall'
    return scalar mean_stage1 = `mean_stage1'
    return scalar mean_stage2 = `mean_stage2'
end


// ============================================================================
// MATA FUNCTIONS
// ============================================================================

capture mata: mata drop _dea_network_dmu_read()
capture mata: mata drop _dea_network_compute()
capture mata: mata drop _dea_network_stage_eff()
capture mata: mata drop _dea_network_relational_kh()
capture mata: mata drop _dea_network_multiplicative()
capture mata: mata drop _dea_network_additive()

mata:
mata set matastrict on

// ─────────────────────────────────────────────────────────────────────────────
// _dea_network_dmu_read() — handles string and numeric DMU variables
// ─────────────────────────────────────────────────────────────────────────────
void _dea_network_dmu_read(string scalar  dmu_var,
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
 * _dea_network_compute()
 *
 * Computes Network DEA efficiency scores.
 */
void _dea_network_compute(string scalar inputs,
                           string scalar intermediates,
                           string scalar outputs,
                           string scalar dmu_var,
                           string scalar touse,
                           string scalar eff_overall_var,
                           string scalar eff_stage1_var,
                           string scalar eff_stage2_var,
                           string scalar method,
                           string scalar rts,
                           string scalar orientation,
                           real scalar tol,
                           real scalar show_detail,
                           real scalar show_debug)
{
    // Load data
    real matrix X, Z, Y
    real colvector dmu_ids, obs_indices
    string colvector dmu_labels
    string rowvector input_names, inter_names, output_names
    
    input_names = tokens(inputs)
    inter_names = tokens(intermediates)
    output_names = tokens(outputs)
    
    X = st_data(., input_names, touse)
    Z = st_data(., inter_names, touse)
    Y = st_data(., output_names, touse)
    _dea_network_dmu_read(dmu_var, touse, dmu_ids, dmu_labels)
    obs_indices = selectindex(st_data(., touse))
    
    real scalar n, m, d, s, i, input_oriented, lw
    string scalar fmt_dmu
    n = rows(X)
    m = cols(X)
    d = cols(Z)
    s = cols(Y)
    
    if (show_debug) {
        printf("  n=%g, m=%g (inputs), d=%g (intermediates), s=%g (outputs)\n", 
               n, m, d, s)
    }
    
    // Store efficiency scores
    real colvector eff_overall, eff_stage1, eff_stage2
    real rowvector kh_eff
    eff_overall = J(n, 1, .)
    eff_stage1 = J(n, 1, .)
    eff_stage2 = J(n, 1, .)
    
    input_oriented = (orientation == "input")
    
    if (method == "independent") {
        // Independent model: compute each stage separately
        for (i = 1; i <= n; i++) {
            // Stage 1: X → Z
            eff_stage1[i] = _dea_network_stage_eff(X[i, .], Z[i, .], 
                                                    X, Z, rts, tol, input_oriented)
            
            // Stage 2: Z → Y
            eff_stage2[i] = _dea_network_stage_eff(Z[i, .], Y[i, .], 
                                                    Z, Y, rts, tol, input_oriented)
            
            // Overall = product
            eff_overall[i] = eff_stage1[i] * eff_stage2[i]
        }
    }
    else if (method == "multiplicative") {
        // Multiplicative relational model (Kao & Hwang, 2008)
        for (i = 1; i <= n; i++) {
            if (rts == "crs" & input_oriented) {
                kh_eff = _dea_network_relational_kh(X[i, .], Z[i, .], Y[i, .],
                                                    X, Z, Y, tol)
                eff_overall[i] = kh_eff[1]
                eff_stage1[i]  = kh_eff[2]
                eff_stage2[i]  = kh_eff[3]
            }
            else {
                // VRS or output-oriented network decompositions are not the
                // Kao-Hwang CCR multiplier model; retain the legacy fallback.
                eff_stage1[i] = _dea_network_stage_eff(X[i, .], Z[i, .],
                                                        X, Z, rts, tol, input_oriented)
                eff_stage2[i] = _dea_network_stage_eff(Z[i, .], Y[i, .],
                                                        Z, Y, rts, tol, input_oriented)
                eff_overall[i] = eff_stage1[i] * eff_stage2[i]
            }
        }
    }
    else {
        // Additive model
        for (i = 1; i <= n; i++) {
            // Stage 1: X → Z
            eff_stage1[i] = _dea_network_stage_eff(X[i, .], Z[i, .], 
                                                    X, Z, rts, tol, input_oriented)
            
            // Stage 2: Z → Y
            eff_stage2[i] = _dea_network_stage_eff(Z[i, .], Y[i, .], 
                                                    Z, Y, rts, tol, input_oriented)
            
            // Additive decomposition with equal weights
            eff_overall[i] = 0.5 * eff_stage1[i] + 0.5 * eff_stage2[i]
        }
    }
    
    // Display results
    printf("\n")
    printf("Network DEA Results\n")
    printf("{hline 60}\n")
    lw = 10
    for (i = 1; i <= n; i++) {
        if (strlen(dmu_labels[i]) > lw) lw = strlen(dmu_labels[i])
    }
    fmt_dmu = "%" + strofreal(lw) + "s"
    printf(fmt_dmu + " | %10s %10s %10s\n", "DMU", "Stage 1", "Stage 2", "Overall")
    printf(fmt_dmu + " | %10s %10s %10s\n", "-":*lw, "----------", "----------", "----------")
    
    for (i = 1; i <= n; i++) {
        printf(fmt_dmu + " | %10.4f %10.4f %10.4f\n", 
               dmu_labels[i], eff_stage1[i], eff_stage2[i], eff_overall[i])
    }
    printf("{hline 60}\n")
    
    // Store results in Stata
    for (i = 1; i <= n; i++) {
        st_store(obs_indices[i], st_varindex(eff_overall_var), eff_overall[i])
        st_store(obs_indices[i], st_varindex(eff_stage1_var), eff_stage1[i])
        st_store(obs_indices[i], st_varindex(eff_stage2_var), eff_stage2[i])
    }
}

/**
 * _dea_network_relational_kh()
 *
 * Kao-Hwang two-stage relational multiplier model under CRS/input orientation.
 *
 * Variables are ordered as v (inputs), a (intermediates), u (outputs).
 * The same intermediate weights a are used as Stage 1 output weights and
 * Stage 2 input weights:
 *
 *   max u*y0
 *   s.t. u*yj - v*xj <= 0             overall system
 *        a*zj - v*xj <= 0             stage 1
 *        u*yj - a*zj <= 0             stage 2
 *        v*x0 = 1
 *        v, a, u >= 0
 *
 * Then E = u*y0, E1 = a*z0, and E2 = E/E1.
 */
real rowvector _dea_network_relational_kh(real rowvector x0,
                                          real rowvector z0,
                                          real rowvector y0,
                                          real matrix X,
                                          real matrix Z,
                                          real matrix Y,
                                          real scalar tol)
{
    real scalar n, m, d, s, nvars, val, mid
    real rowvector c, lb, ub, sol
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    class LinearProgram scalar S, S2

    n = rows(X)
    m = cols(X)
    d = cols(Z)
    s = cols(Y)
    nvars = m + d + s

    // Minimize -u*y0 to maximize overall efficiency.
    c = (J(1, m, 0), J(1, d, 0), -y0)

    A_ineq = ( -X,  J(n, d, 0),  Y ) \ ///
             ( -X,  Z,           J(n, s, 0) ) \ ///
             ( J(n, m, 0), -Z,   Y )
    b_ineq = J(3*n, 1, 0)

    A_eq = (x0, J(1, d, 0), J(1, s, 0))
    b_eq = 1

    lb = J(1, nvars, 0)
    ub = J(1, nvars, .)

    S = LinearProgram()
    S.setMaxOrMin("min")
    S.setCoefficients(c)
    S.setInequality(A_ineq, b_ineq)
    S.setEquality(A_eq, b_eq)
    S.setBounds(lb, ub)
    S.setTrace("off")

    (void) S.optimize()

    if (S.converged() && S.errorcode() == 0) {
        sol = S.parameters()
        val = -S.value()

        if (val < 0 & abs(val) < tol) val = 0
        if (val > 1 & val < 1 + 100*tol) val = 1

        // Kao-Hwang decomposition: keep the optimal overall efficiency and
        // choose the largest first-stage efficiency among alternative optima.
        c = (J(1, m, 0), -z0, J(1, s, 0))
        A_eq = ( (x0, J(1, d, 0), J(1, s, 0)) \
                 (J(1, m, 0), J(1, d, 0), y0) )
        b_eq = (1 \ val)

        S2 = LinearProgram()
        S2.setMaxOrMin("min")
        S2.setCoefficients(c)
        S2.setInequality(A_ineq, b_ineq)
        S2.setEquality(A_eq, b_eq)
        S2.setBounds(lb, ub)
        S2.setTrace("off")

        (void) S2.optimize()

        if (S2.converged() && S2.errorcode() == 0) {
            sol = S2.parameters()
        }

        mid = sol[(m+1)..(m+d)] * z0'
        if (mid < 0 & abs(mid) < tol) mid = 0
        if (mid > 1 & mid < 1 + 100*tol) mid = 1

        real scalar stage2
        stage2 = (abs(mid) > tol ? val / mid : .)
        if (stage2 > 1 & stage2 < 1 + 100*tol) stage2 = 1
        return((val, mid, stage2))
    }
    else {
        return((., ., .))
    }
}


/**
 * _dea_network_stage_eff()
 *
 * Computes standard DEA efficiency for a single stage.
 * Used in independent model.
 */
real scalar _dea_network_stage_eff(real rowvector x0,
                                    real rowvector y0,
                                    real matrix X_ref,
                                    real matrix Y_ref,
                                    string scalar rts,
                                    real scalar tol,
                                    real scalar input_oriented)
{
    real scalar n, m, s, eff
    real rowvector c, lb, ub
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    class LinearProgram scalar S
    
    n = rows(X_ref)
    m = cols(X_ref)
    s = cols(Y_ref)
    
    if (input_oriented) {
        // min theta s.t. X'lambda <= theta*x0, Y'lambda >= y0
        c = (J(1, n, 0), 1)
        
        A_ineq = (X_ref', -x0') \ (-Y_ref', J(s, 1, 0))
        b_ineq = J(m, 1, 0) \ (-y0')
        
        lb = J(1, n + 1, 0)
        ub = J(1, n + 1, .)
        
        A_eq = J(0, n + 1, .)
        b_eq = J(0, 1, .)
        
        if (rts == "vrs") {
            A_eq = (J(1, n, 1), 0)
            b_eq = 1
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
        
        (void) S.optimize()
        
        if (S.converged() && S.errorcode() == 0) {
            eff = S.value()
        }
        else {
            eff = .
        }
    }
    else {
        // Output-oriented: max phi
        c = (J(1, n, 0), -1)
        
        A_ineq = (X_ref', J(m, 1, 0)) \ (-Y_ref', y0')
        b_ineq = x0' \ J(s, 1, 0)
        
        lb = (J(1, n, 0), 1)
        ub = J(1, n + 1, .)
        
        A_eq = J(0, n + 1, .)
        b_eq = J(0, 1, .)
        
        if (rts == "vrs") {
            A_eq = (J(1, n, 1), 0)
            b_eq = 1
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
        
        (void) S.optimize()
        
        if (S.converged() && S.errorcode() == 0) {
            real scalar phi
            phi = -S.value()
            eff = 1 / phi
        }
        else {
            eff = .
        }
    }
    
    return(eff)
}


/**
 * _dea_network_multiplicative()
 *
 * Multiplicative relational model (Kao & Hwang, 2008).
 * 
 * For input-oriented CRS:
 *   E_overall = E_1 × E_2
 *   
 * We first compute overall efficiency, then decompose.
 */
void _dea_network_multiplicative(real rowvector x0,
                                  real rowvector z0,
                                  real rowvector y0,
                                  real matrix X,
                                  real matrix Z,
                                  real matrix Y,
                                  string scalar rts,
                                  real scalar tol,
                                  real scalar eff_overall,
                                  real scalar eff_stage1,
                                  real scalar eff_stage2)
{
    real scalar n, m, d, s
    n = rows(X)
    m = cols(X)
    d = cols(Z)
    s = cols(Y)
    
    // For multiplicative model, first solve for overall efficiency
    // using system model where X → Z → Y
    
    // Stage 1 efficiency: X → Z
    eff_stage1 = _dea_network_stage_eff(x0, z0, X, Z, rts, tol, 1)
    
    // Stage 2 efficiency: Z → Y  
    eff_stage2 = _dea_network_stage_eff(z0, y0, Z, Y, rts, tol, 1)
    
    // Overall efficiency (multiplicative decomposition)
    eff_overall = eff_stage1 * eff_stage2
}


/**
 * _dea_network_additive()
 *
 * Additive model with equal weights.
 *   E_overall = 0.5 * E_1 + 0.5 * E_2
 */
void _dea_network_additive(real rowvector x0,
                            real rowvector z0,
                            real rowvector y0,
                            real matrix X,
                            real matrix Z,
                            real matrix Y,
                            string scalar rts,
                            real scalar tol,
                            real scalar eff_overall,
                            real scalar eff_stage1,
                            real scalar eff_stage2)
{
    // Compute stage efficiencies independently
    eff_stage1 = _dea_network_stage_eff(x0, z0, X, Z, rts, tol, 1)
    eff_stage2 = _dea_network_stage_eff(z0, y0, Z, Y, rts, tol, 1)
    
    // Additive decomposition with equal weights
    eff_overall = 0.5 * eff_stage1 + 0.5 * eff_stage2
}

end
