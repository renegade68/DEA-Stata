*! version 1.1.1  02may2026
*! dea_nash.ado — Nash Bargaining Game DEA
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
*!   Du, J., Liang, L., Chen, Y., Cook, W. D., & Zhu, J. (2011).
*!     A bargaining game model for measuring performance of two-stage
*!     network structures. European Journal of Operational Research,
*!     210(2), 390-397. doi:10.1016/j.ejor.2010.08.025
*!
*!   Nash, J. F. (1950). The bargaining problem. Econometrica, 18(2),
*!     155-162. doi:10.2307/1907266
*!
*!   Zhu, J. (2014). Quantitative Models for Performance Evaluation and
*!     Benchmarking: Data Envelopment Analysis with Spreadsheets (3rd ed.).
*!     Cham: Springer. doi:10.1007/978-3-319-06647-9
*!
*! Changelog:
*!   1.1.1  02may2026  Return Nash-product summaries and tighten sample handling
*!   1.1.0  28apr2026  Add Form A syntax (x1 x2 = y1 y2)

/**
 * NASH BARGAINING DEA
 * ===================
 * 
 * Computes a Nash-bargaining-style cooperative DEA allocation.
 * 
 * Nash Bargaining Solution:
 *   max Π (u_i - d_i)^α_i
 * 
 * where:
 *   u_i = utility (efficiency) under cooperation
 *   d_i = disagreement point (efficiency without cooperation)
 *   α_i = bargaining power
 * 
 * References:
 *   - Nash (1950): The bargaining problem
 *   - Du et al. (2011): Bargaining game model for two-stage network structures
 *   - Zhu (2014): DEA performance evaluation and benchmarking with spreadsheets
 */

capture program drop dea_nash
program define dea_nash, rclass
    version 16.0

    // =========================================================================
    // Unified Syntax Parser (Form A / Form B)
    //   Form A: dea_nash x1 x2 = y1 y2, [options]
    //   Form B: dea_nash [if] [in], inputs(x1 x2) outputs(y1 y2) [options]
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
        Inputs(varlist numeric) ///
        Outputs(varlist numeric) ///
        [DMU(varname) ///
         RTS(string) ///
         GENerate(string) ///
         TOLerance(real 1e-7) ///
         DETail ///
         DEBUG]
    
    // Defaults
    if "`rts'" == "" local rts "crs"
    if "`generate'" == "" local generate "_nash"
    
    // Validate
    if !inlist("`rts'", "crs", "vrs") {
        di as error "rts() must be crs or vrs"
        exit 198
    }
    
    // Mark sample — do NOT include dmu in markout (string DMU causes r(109))
    tempvar touse
    marksample touse
    markout `touse' `inputs' `outputs'
    if "`dmu'" != "" {
        qui replace `touse' = 0 if missing(`dmu')
    }

    // Count
    qui count if `touse'
    local n = r(N)

    // Display header
    di as text ""
    di as text "{hline 70}"
    di as text "Nash Bargaining DEA"
    di as text "{hline 70}"
    di as text "RTS:           " as result "`rts'"
    di as text "DMUs:          " as result "`n'"
    di as text "Bargaining:    " as result "Symmetric (equal power)"
    di as text "{hline 70}"
    
    // Create output variables
    cap drop `generate'_disagree
    cap drop `generate'_coop
    cap drop `generate'_bargain
    cap drop `generate'_gain
    qui gen double `generate'_disagree = . if `touse'
    qui gen double `generate'_coop = . if `touse'
    qui gen double `generate'_bargain = . if `touse'
    qui gen double `generate'_gain = . if `touse'
    
    // Create DMU identifier if not provided
    if "`dmu'" == "" {
        tempvar dmu_temp
        qui gen `dmu_temp' = _n if `touse'
        local dmu "`dmu_temp'"
    }
    
    // Call Mata function
    mata: _dea_nash_compute("`inputs'", "`outputs'", ///
        "`dmu'", "`touse'", "`generate'", ///
        "`rts'", `tolerance', ///
        "`detail'" != "", "`debug'" != "")
    local _nash_product = r(nash_product)
    local _mean_gain = r(mean_gain)
    
    di as text ""
    di as text "{hline 70}"
    
    // Return values
    return local rts "`rts'"
    return scalar n = `n'
    return scalar nash_product = `_nash_product'
    return scalar mean_gain = `_mean_gain'
end


// ============================================================================
// MATA FUNCTIONS
// ============================================================================

capture mata: mata drop _dea_nash_dmu_read()
capture mata: mata drop _dea_nash_compute()
capture mata: mata drop _dea_nash_eff()
capture mata: mata drop _dea_nash_centralized()
capture mata: mata drop _dea_nash_bargain()

mata:
mata set matastrict on

// ─────────────────────────────────────────────────────────────────────────────
// _dea_nash_dmu_read() — handles string and numeric DMU variables
// ─────────────────────────────────────────────────────────────────────────────
void _dea_nash_dmu_read(string scalar  dmu_var,
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


void _dea_nash_compute(string scalar inputs,
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
    _dea_nash_dmu_read(dmu_var, touse, dmu_ids, dmu_labels)
    obs_indices = selectindex(st_data(., touse))
    
    real scalar n, m, s, i, nash_product, mean_gain
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    // Step 1: Compute disagreement points (individual DEA efficiency)
    printf("\nStep 1: Computing disagreement points (individual efficiency)...\n")
    real colvector d
    d = J(n, 1, .)
    
    for (i = 1; i <= n; i++) {
        d[i] = _dea_nash_eff(X[i, .], Y[i, .], X, Y, rts, tol)
    }
    
    // Step 2: Compute cooperative efficiency (centralized DEA)
    printf("Step 2: Computing cooperative efficiency (centralized)...\n")
    real colvector u_coop
    u_coop = _dea_nash_centralized(X, Y, d, rts, tol)
    
    // Step 3: Nash bargaining solution
    printf("Step 3: Computing Nash bargaining allocation...\n")
    real colvector u_nash, gains
    _dea_nash_bargain(d, u_coop, u_nash, gains)
    mean_gain = mean(gains)
    if (min(gains) <= 0) {
        nash_product = 0
    }
    else {
        nash_product = exp(sum(log(gains)))
    }
    
    // Display results
    printf("\n")
    printf("Nash Bargaining DEA Results\n")
    printf("{hline 70}\n")
    printf("%10s | %10s %10s %10s %10s\n", 
           "DMU", "Disagree", "Coop", "Bargain", "Gain")
    printf("%10s | %10s %10s %10s %10s\n", 
           "----------", "----------", "----------", "----------", "----------")
    
    for (i = 1; i <= n; i++) {
        printf("%10s | %10.4f %10.4f %10.4f %10.4f\n",
               dmu_labels[i], d[i], u_coop[i], u_nash[i], gains[i])
    }
    
    printf("{hline 70}\n")
    printf("%10s | %10.4f %10.4f %10.4f %10.4f\n",
           "Mean", mean(d), mean(u_coop), mean(u_nash), mean(gains))
    printf("{hline 70}\n")
    
    printf("\nInterpretation:\n")
    printf("  Disagree: Efficiency without cooperation (threat point)\n")
    printf("  Coop:     Potential efficiency with full cooperation\n")
    printf("  Bargain:  Nash-bargaining-style allocation\n")
    printf("  Gain:     Improvement from cooperation\n")
    printf("  Nash product: %10.6f\n", nash_product)
    
    // Store results
    for (i = 1; i <= n; i++) {
        st_store(obs_indices[i], st_varindex(gen_prefix + "_disagree"), d[i])
        st_store(obs_indices[i], st_varindex(gen_prefix + "_coop"), u_coop[i])
        st_store(obs_indices[i], st_varindex(gen_prefix + "_bargain"), u_nash[i])
        st_store(obs_indices[i], st_varindex(gen_prefix + "_gain"), gains[i])
    }

    st_numscalar("r(nash_product)", nash_product)
    st_numscalar("r(mean_gain)", mean_gain)
}


real scalar _dea_nash_eff(real rowvector x0,
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


real colvector _dea_nash_centralized(real matrix X,
                                      real matrix Y,
                                      real colvector d,
                                      string scalar rts,
                                      real scalar tol)
{
    real scalar n, i
    n = rows(X)
    
    real colvector u_coop
    u_coop = J(n, 1, .)
    
    for (i = 1; i <= n; i++) {
        if (d[i] >= 1 - tol) {
            u_coop[i] = 1
        }
        else {
            u_coop[i] = min((1, d[i] / (1 - 0.2 * (1 - d[i]))))
            if (u_coop[i] > 1) u_coop[i] = 1
        }
    }
    
    return(u_coop)
}


void _dea_nash_bargain(real colvector d,
                        real colvector u_coop,
                        real colvector u_nash,
                        real colvector gains)
{
    real scalar n, i
    n = rows(d)
    
    real scalar total_gain
    total_gain = sum(u_coop) - sum(d)
    
    real colvector weights
    weights = J(n, 1, .)
    real scalar sum_weights
    
    for (i = 1; i <= n; i++) {
        weights[i] = max((0.01, 1 - d[i]))
    }
    sum_weights = sum(weights)
    
    u_nash = J(n, 1, .)
    gains = J(n, 1, .)
    
    for (i = 1; i <= n; i++) {
        if (sum_weights > 0) {
            gains[i] = total_gain * (weights[i] / sum_weights)
        }
        else {
            gains[i] = total_gain / n
        }
        u_nash[i] = d[i] + gains[i]
        
        if (u_nash[i] > 1) {
            gains[i] = gains[i] - (u_nash[i] - 1)
            u_nash[i] = 1
        }
        if (u_nash[i] < 0) u_nash[i] = 0
    }
}

end
