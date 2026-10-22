*! version 1.0.1  03may2026
*! dea_fuzzy.ado — Fuzzy DEA via α-cut approach
*! Requires: Stata 16.0 or later (LinearProgram class)
*!
*! AUTHORS:
*! --------
*!   Choonjoo Lee (corresponding author)
*!     Korea National Defense University
*!     Seoul, Korea
*!     bloom.rampike@gmail.com; bloom.rampike@gmail.com
*!
*! METHOD:
*! -------
*! Implements fuzzy DEA using the α-cut (possibility) approach.
*! For each α ∈ [0,1], triangular fuzzy inputs/outputs are sliced into
*! crisp intervals [x_L(α), x_U(α)], and two parametric LPs are solved
*! per DMU:
*!
*!   Pessimistic LP: DMU k uses unfavorable bounds (high inputs, low outputs);
*!                   peer frontier uses favorable bounds, with DMU k itself
*!                   entering the reference set at its same unfavorable bounds.
*!                   → lower bound on efficiency  θ_L(α)
*!
*!   Optimistic LP:  DMU k uses favorable bounds (low inputs, high outputs);
*!                   peer frontier uses unfavorable bounds, with DMU k itself
*!                   entering the reference set at its same favorable bounds.
*!                   → upper bound on efficiency  θ_U(α)
*!
*! α-cut bounds for triangular fuzzy number (x_lo, x_center, x_hi):
*!   x_L(α) = x_lo    + α * (x_center − x_lo)      [lo at α=0, center at α=1]
*!   x_U(α) = x_hi    − α * (x_hi     − x_center)  [hi at α=0, center at α=1]
*!
*! At α = 1: x_L = x_U = x_center  →  collapses to standard (crisp) DEA.
*! At α = 0: full fuzzy spread      →  maximum parametric uncertainty.
*!
*! For output-oriented models, φ (≥ 1) is converted to 1/φ ∈ (0,1] for
*! display, consistent with dea.ado conventions for IDEA output-oriented.
*!
*! REFERENCES:
*! -----------
*! Kao, C., & Liu, S.-T. (2000). Fuzzy efficiency measures in data
*!   envelopment analysis. Fuzzy Sets and Systems, 113(3), 427-437.
*!   → α-cut LP formulation (pessimistic/optimistic)
*!
*! Lertworasirikul, S., Fang, S.-C., Joines, J. A., & Nuttle, H. L. W.
*!   (2003). Fuzzy data envelopment analysis (DEA): a possibility approach.
*!   Fuzzy Sets and Systems, 139(2), 379-394.
*!   → Possibility measure DEA; pessimistic = lower α-cut efficiency
*!
*! Guo, P., & Tanaka, H. (2001). Fuzzy DEA: a perceptual evaluation method.
*!   Fuzzy Sets and Systems, 119(1), 149-160.
*!   → Fuzzy efficiency intervals and membership functions
*!
*! Cooper, W. W., Park, K. S., & Yu, G. (1999). IDEA and AR-IDEA: Models
*!   for dealing with imprecise data in DEA. Management Science, 45(4),
*!   597-607.
*!   → IDEA LP structure (pessimistic/optimistic); basis for this extension
*!
*! Charnes, A., Cooper, W. W., & Rhodes, E. (1978). Measuring the efficiency
*!   of decision making units. European Journal of Operational Research,
*!   2, 429-444.
*!
*! Banker, R. D., Charnes, A., & Cooper, W. W. (1984). Some models for
*!   estimating technical and scale inefficiencies in data envelopment
*!   analysis. Management Science, 30, 1078-1092.
*!
*! Ji, Y.-B., & Lee, C. (2010). Data envelopment analysis. Stata Journal,
*!   10(2), 267-280.
*!
*! Changelog:
*!   1.0.1  03may2026  Align best/worst reference-set convention with
*!                      deaR::modelfuzzy_kaoliu (evaluated DMU enters the
*!                      reference set at the same alpha-cut realization).
*!   1.0.0  07apr2026  Initial packaged version; based on earlier teaching/seminar versions


// ============================================================================
// Compile Mata function when this file is (re)loaded
// ============================================================================
capture mata: mata drop _dea_fuzzy_solve()

mata:
mata set matastrict on

// ----------------------------------------------------------------------------
// _dea_fuzzy_solve()
//
// Solves pessimistic and optimistic Fuzzy DEA LPs for all DMUs at one
// α-cut level.
//
// Arguments:
//   inputs      - space-separated center input variable names
//   outputs     - space-separated center output variable names
//   xlo_vars    - space-separated lower-bound input variable names
//   xhi_vars    - space-separated upper-bound input variable names
//   ylo_vars    - space-separated lower-bound output variable names
//   yhi_vars    - space-separated upper-bound output variable names
//   alpha       - α-cut level ∈ [0,1]
//   orientation - "input" or "output"
//   rts         - "crs", "vrs", "drs", or "irs"
//   tol         - numerical tolerance (unused; reserved for future use)
//   pess_var    - variable name to store pessimistic efficiency scores
//   opt_var     - variable name to store optimistic efficiency scores
//   touse       - sample marker variable name
//
// Pessimistic input-oriented LP:
//   min θ  s.t.  X_L(α)' λ  ≤  θ · x_{k,U}(α)    [frontier best, DMU worst]
//               −Y_U(α)' λ  ≤ −y_{k,L}(α)
//                         λ  ≥  0
//
// Optimistic input-oriented LP:
//   min θ  s.t.  X_U(α)' λ  ≤  θ · x_{k,L}(α)    [frontier worst, DMU best]
//               −Y_L(α)' λ  ≤ −y_{k,U}(α)
//                         λ  ≥  0
//
// Output-oriented: analogously; result stored as 1/φ ∈ (0,1].
// ----------------------------------------------------------------------------

void _dea_fuzzy_solve(
    string scalar inputs,
    string scalar outputs,
    string scalar xlo_vars,
    string scalar xhi_vars,
    string scalar ylo_vars,
    string scalar yhi_vars,
    real scalar   alpha,
    string scalar orientation,
    string scalar rts,
    real scalar   tol,
    string scalar pess_var,
    string scalar opt_var,
    string scalar touse)
{
    // ----------------------------------------------------------------
    // Load data
    // ----------------------------------------------------------------
    real matrix X_c, X_lo, X_hi, Y_c, Y_lo, Y_hi
    real matrix X_L, X_U, Y_L, Y_U, X_ref, Y_ref
    real colvector obs_indices
    real scalar n, m, s, k, obs_k, rc

    X_c  = st_data(., tokens(inputs),   touse)
    X_lo = st_data(., tokens(xlo_vars), touse)
    X_hi = st_data(., tokens(xhi_vars), touse)
    Y_c  = st_data(., tokens(outputs),  touse)
    Y_lo = st_data(., tokens(ylo_vars), touse)
    Y_hi = st_data(., tokens(yhi_vars), touse)

    obs_indices = selectindex(st_data(., touse))
    n = rows(X_c)
    m = cols(X_c)
    s = cols(Y_c)

    // ----------------------------------------------------------------
    // Compute α-cut bounds
    //   x_L(α) = x_lo + α*(x_center − x_lo)   [lo → center as α: 0→1]
    //   x_U(α) = x_hi − α*(x_hi − x_center)   [hi → center as α: 0→1]
    // ----------------------------------------------------------------
    X_L = X_lo + alpha :* (X_c  - X_lo)   // n × m
    X_U = X_hi - alpha :* (X_hi - X_c )   // n × m
    Y_L = Y_lo + alpha :* (Y_c  - Y_lo)   // n × s
    Y_U = Y_hi - alpha :* (Y_hi - Y_c )   // n × s

    // ----------------------------------------------------------------
    // Result vectors
    // ----------------------------------------------------------------
    real colvector eff_pess, eff_opt
    eff_pess = J(n, 1, .)
    eff_opt  = J(n, 1, .)

    // ----------------------------------------------------------------
    // LP building blocks
    // ----------------------------------------------------------------
    real rowvector c, lb, ub
    real matrix    A_ineq, A_eq
    real colvector b_ineq, b_eq, x0, y0
    real scalar    phi_val
    class LinearProgram scalar S

    for (k = 1; k <= n; k++) {

        obs_k = obs_indices[k]

        if (orientation == "input") {
            // ============================================================
            // INPUT-ORIENTED
            //
            // Variables: [λ₁, …, λₙ, θ]    (n+1 variables)
            //
            // Pessimistic:
            //   DMU k worst → x0 = X_U[k], y0 = Y_L[k]
            //   Peer frontier best; DMU k enters reference at worst case
            //   min θ  s.t.  X_L'λ − θ·X_U[k] ≤ 0   (m rows)
            //                −Y_U'λ             ≤ −Y_L[k]   (s rows)
            //                λ ≥ 0, θ ≥ 0
            //
            // Optimistic:
            //   DMU k best  → x0 = X_L[k], y0 = Y_U[k]
            //   Peer frontier worst; DMU k enters reference at best case
            //   min θ  s.t.  X_U'λ − θ·X_L[k] ≤ 0   (m rows)
            //                −Y_L'λ             ≤ −Y_U[k]   (s rows)
            //                λ ≥ 0, θ ≥ 0
            // ============================================================

            c  = J(1, n, 0), 1
            lb = J(1, n+1, 0)
            ub = J(1, n+1, .)

            // --- PESSIMISTIC ---
            x0 = X_U[k, ]'
            y0 = Y_L[k, ]'
            X_ref = X_L
            Y_ref = Y_U
            X_ref[k, ] = X_U[k, ]
            Y_ref[k, ] = Y_L[k, ]

            A_ineq = (X_ref', -x0) \ (-Y_ref', J(s, 1, 0))
            b_ineq = J(m, 1, 0) \ (-y0)

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
            if (S.converged() & S.errorcode() == 0) {
                eff_pess[k] = S.value()
            }

            // --- OPTIMISTIC ---
            x0 = X_L[k, ]'
            y0 = Y_U[k, ]'
            X_ref = X_U
            Y_ref = Y_L
            X_ref[k, ] = X_L[k, ]
            Y_ref[k, ] = Y_U[k, ]

            A_ineq = (X_ref', -x0) \ (-Y_ref', J(s, 1, 0))
            b_ineq = J(m, 1, 0) \ (-y0)

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
            if (S.converged() & S.errorcode() == 0) {
                eff_opt[k] = S.value()
            }
        }
        else {
            // ============================================================
            // OUTPUT-ORIENTED
            //
            // Variables: [λ₁, …, λₙ, φ]    (n+1 variables)
            //
            // Standard output-oriented LP (min −φ form):
            //   c       = [0,…,0, −1]
            //   A_ineq  = (X_ref', J(m,1,0)) row-block
            //             \ (−Y_ref', y0)     row-block
            //   b_ineq  = x0 \ J(s,1,0)
            //   lb      = J(1,n,0), 1     (φ ≥ 1)
            //   result  stored as 1/φ ∈ (0,1]  (consistent with IDEA)
            //
            // Pessimistic:
            //   DMU k worst → x0 = X_U[k], y0 = Y_L[k]
            //   Peer frontier best; DMU k enters reference at worst case
            //   → large φ (less efficient); stored as 1/φ_pess (small)
            //
            // Optimistic:
            //   DMU k best  → x0 = X_L[k], y0 = Y_U[k]
            //   Peer frontier worst; DMU k enters reference at best case
            //   → small φ (more efficient); stored as 1/φ_opt (large)
            // ============================================================

            c  = J(1, n, 0), -1
            lb = J(1, n, 0), 1
            ub = J(1, n+1, .)

            // --- PESSIMISTIC ---
            x0 = X_U[k, ]'
            y0 = Y_L[k, ]'
            X_ref = X_L
            Y_ref = Y_U
            X_ref[k, ] = X_U[k, ]
            Y_ref[k, ] = Y_L[k, ]

            A_ineq = (X_ref', J(m, 1, 0)) \ (-Y_ref', y0)
            b_ineq = x0 \ J(s, 1, 0)

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
            if (S.converged() & S.errorcode() == 0) {
                phi_val = -S.value()                 // φ ≥ 1
                if (phi_val > 0) eff_pess[k] = 1 / phi_val   // 1/φ ∈ (0,1]
            }

            // --- OPTIMISTIC ---
            x0 = X_L[k, ]'
            y0 = Y_U[k, ]'
            X_ref = X_U
            Y_ref = Y_L
            X_ref[k, ] = X_L[k, ]
            Y_ref[k, ] = Y_U[k, ]

            A_ineq = (X_ref', J(m, 1, 0)) \ (-Y_ref', y0)
            b_ineq = x0 \ J(s, 1, 0)

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
            if (S.converged() & S.errorcode() == 0) {
                phi_val = -S.value()
                if (phi_val > 0) eff_opt[k] = 1 / phi_val
            }
        }
    }

    // ----------------------------------------------------------------
    // Store results back to Stata
    // ----------------------------------------------------------------
    st_store(obs_indices, st_varindex(pess_var), eff_pess)
    st_store(obs_indices, st_varindex(opt_var),  eff_opt)
}

end   // end mata


// ============================================================================
// Main Stata program
// ============================================================================
capture program drop dea_fuzzy
program define dea_fuzzy, rclass
    version 16.0

    // =========================================================================
    // Syntax
    // =========================================================================
    syntax [if] [in] , ///
        INputs(varlist numeric)  ///
        OUTputs(varlist numeric) ///
        [XLO(varlist numeric)    ///
         XHI(varlist numeric)    ///
         YLO(varlist numeric)    ///
         YHI(varlist numeric)    ///
         ALPha(numlist >=0 <=1 sort) ///
         RTS(string)             ///
         ORT(string)             ///
         ORIentation(string)     ///
         DMU(varname)            ///
         GENerate(string)        ///
         TOLerance(real 1e-7)    ///
         DETail                  ///
         DEBUG]

    // =========================================================================
    // Defaults
    // =========================================================================
    if "`alpha'"       == "" local alpha "0 0.25 0.5 0.75 1"
    if "`rts'"         == "" local rts "crs"
    if "`generate'"    == "" local generate "fdea_eff"

    // Resolve orientation (ort() is an alias)
    if "`orientation'" == "" & "`ort'" != "" local orientation "`ort'"
    if "`orientation'" == "" local orientation "input"

    // Normalize RTS
    local rts = lower("`rts'")
    if "`rts'" == "ccr"  local rts "crs"
    if "`rts'" == "bcc"  local rts "vrs"
    if "`rts'" == "ndrs" local rts "drs"
    if "`rts'" == "nirs" local rts "irs"
    if !inlist("`rts'", "crs", "vrs", "drs", "irs") {
        di as error "rts() must be one of: crs ccr vrs bcc drs irs"
        exit 198
    }

    // Normalize orientation
    local orientation = lower("`orientation'")
    if substr("`orientation'", 1, 1) == "i" local orientation "input"
    if substr("`orientation'", 1, 1) == "o" local orientation "output"
    if !inlist("`orientation'", "input", "output") {
        di as error "orientation() must be 'input' or 'output'"
        exit 198
    }

    // Count inputs / outputs
    local m_in  : word count `inputs'
    local m_out : word count `outputs'

    // =========================================================================
    // Validate bounds: if provided, must match count of inputs/outputs
    // =========================================================================
    if "`xlo'" != "" {
        local n_xlo : word count `xlo'
        if `n_xlo' != `m_in' {
            di as error "xlo() must contain `m_in' variable(s) (same as inputs())"
            exit 198
        }
    }
    if "`xhi'" != "" {
        local n_xhi : word count `xhi'
        if `n_xhi' != `m_in' {
            di as error "xhi() must contain `m_in' variable(s) (same as inputs())"
            exit 198
        }
    }
    if "`ylo'" != "" {
        local n_ylo : word count `ylo'
        if `n_ylo' != `m_out' {
            di as error "ylo() must contain `m_out' variable(s) (same as outputs())"
            exit 198
        }
    }
    if "`yhi'" != "" {
        local n_yhi : word count `yhi'
        if `n_yhi' != `m_out' {
            di as error "yhi() must contain `m_out' variable(s) (same as outputs())"
            exit 198
        }
    }

    // Warn if no fuzzy bounds supplied (degenerates to standard DEA)
    if "`xlo'`xhi'`ylo'`yhi'" == "" {
        di as text "Note: No fuzzy bounds specified; all variables treated as crisp."
        di as text "      Results will be identical to standard DEA at all α-levels."
    }

    // =========================================================================
    // Build sample
    // =========================================================================
    tempvar touse
    mark    `touse'
    markout `touse' `inputs' `outputs' `xlo' `xhi' `ylo' `yhi'

    if "`dmu'" == "" {
        cap drop _dmu
        qui gen str _dmu = "dmu" + string(_n) if `touse'
        local dmu "_dmu"
    }

    qui count if `touse'
    if r(N) == 0 {
        di as error "No observations in estimation sample"
        exit 2000
    }
    local n_dmus = r(N)

    // =========================================================================
    // Resolve bound variable lists
    // (use center values where bounds not specified → crisp for that variable)
    // =========================================================================
    local xlo_vars "`xlo'"
    local xhi_vars "`xhi'"
    local ylo_vars "`ylo'"
    local yhi_vars "`yhi'"

    if "`xlo_vars'" == "" local xlo_vars "`inputs'"
    if "`xhi_vars'" == "" local xhi_vars "`inputs'"
    if "`ylo_vars'" == "" local ylo_vars "`outputs'"
    if "`yhi_vars'" == "" local yhi_vars "`outputs'"

    // =========================================================================
    // Parse alpha list
    // =========================================================================
    local n_alpha : word count `alpha'

    // =========================================================================
    // Orientation label for display
    // =========================================================================
    if "`orientation'" == "input" {
        local ort_label "Input-oriented"
        local eff_label "θ (input efficiency)"
    }
    else {
        local ort_label "Output-oriented"
        local eff_label "1/φ (output efficiency)"
    }

    // RTS label
    if "`rts'" == "crs" local rts_label "CRS"
    if "`rts'" == "vrs" local rts_label "VRS"
    if "`rts'" == "drs" local rts_label "DRS"
    if "`rts'" == "irs" local rts_label "IRS"

    // =========================================================================
    // Header
    // =========================================================================
    di as text "{hline `c(linesize)'}"
    di as text "dea_fuzzy v1.0 — Fuzzy DEA (α-cut / possibility approach)"
    di as text "  `ort_label', `rts_label'"
    di as text "  DMUs: `n_dmus'  |  Inputs: `m_in'  |  Outputs: `m_out'"
    di as text "  α-levels: `alpha'"
    di as text "{hline `c(linesize)'}"

    // =========================================================================
    // Loop over α-levels: create temp vars, call Mata solver
    // =========================================================================
    local ia = 0
    foreach a of numlist `alpha' {
        local ++ia

        // Suffix: round(alpha*100) → 0, 25, 50, 75, 100
        local asuf = string(round(`a' * 100))

        // Temporary variables to hold results for this α level
        tempvar _pess_`ia' _opt_`ia'
        qui gen double `_pess_`ia'' = . if `touse'
        qui gen double `_opt_`ia''  = . if `touse'

        // Call Mata solver
        // quietly suppresses per-LP infeasibility warnings from LinearProgram();
        // a single summary note is printed below if any scores are missing.
        quietly {
            mata: _dea_fuzzy_solve(                ///
                "`inputs'",   "`outputs'",         ///
                "`xlo_vars'", "`xhi_vars'",        ///
                "`ylo_vars'", "`yhi_vars'",        ///
                `a',          "`orientation'",     ///
                "`rts'",      `tolerance',         ///
                "`_pess_`ia''", "`_opt_`ia''",     ///
                "`touse'")
        }

        // Report infeasible LP count (missing scores) for this α-level
        qui count if (`_pess_`ia'' == . | `_opt_`ia'' == .) & `touse'
        if r(N) > 0 {
            di as text "  Note (α=`a'): `r(N)' LP(s) infeasible; " ///
                "corresponding scores set to missing."
        }

        if "`debug'" != "" {
            di as text "  [α=`a'] Mata solver completed."
        }
    }

    // =========================================================================
    // Save results to permanent named variables
    // =========================================================================
    local ia = 0
    foreach a of numlist `alpha' {
        local ++ia
        local asuf = string(round(`a' * 100))

        local vp "`generate'_pess_a`asuf'"
        local vo "`generate'_opt_a`asuf'"

        // Drop existing (silently)
        cap drop `vp' `vo'

        qui gen double `vp' = `_pess_`ia'' if `touse'
        qui gen double `vo' = `_opt_`ia''  if `touse'

        label var `vp' "Fuzzy DEA pessimistic efficiency (α=`a', `rts_label', `ort_label')"
        label var `vo' "Fuzzy DEA optimistic efficiency  (α=`a', `rts_label', `ort_label')"
    }

    // =========================================================================
    // Display: summary table across α-levels
    // =========================================================================
    di as text ""
    di as text "Fuzzy DEA Efficiency Summary  [`eff_label']"
    di as text "{hline 68}"
    di as text "  α      " ///
        _col(12) "Pessimistic" ///
        _col(38) "Optimistic"
    di as text "         " ///
        _col(12) "Mean       [Min ,  Max ]" ///
        _col(38) "Mean       [Min ,  Max ]"
    di as text "{hline 68}"

    local ia = 0
    foreach a of numlist `alpha' {
        local ++ia
        local asuf = string(round(`a' * 100))
        local vp "`generate'_pess_a`asuf'"
        local vo "`generate'_opt_a`asuf'"

        qui sum `vp' if `touse'
        local mp  = string(r(mean), "%8.4f")
        local lop = string(r(min),  "%6.4f")
        local hip = string(r(max),  "%6.4f")

        qui sum `vo' if `touse'
        local mo  = string(r(mean), "%8.4f")
        local loo = string(r(min),  "%6.4f")
        local hio = string(r(max),  "%6.4f")

        di as result "  `a'" ///
            _col(12) "`mp'  [`lop', `hip']" ///
            _col(38) "`mo'  [`loo', `hio']"
    }

    di as text "{hline 68}"
    di as text "Variables created: `generate'_pess_a* (pessimistic)"
    di as text "                   `generate'_opt_a*  (optimistic)"

    // =========================================================================
    // DMU-level detail (if requested)
    // =========================================================================
    if "`detail'" != "" {
        di as text ""
        di as text "DMU-level Fuzzy Efficiency  [`eff_label']"

        // Build header
        local hdr1 "DMU"
        local hdr2 "---"
        foreach a of numlist `alpha' {
            local asuf = string(round(`a' * 100))
            local hdr1 "`hdr1'  α=`a' (pess/opt)"
            local hdr2 "`hdr2'  ----------------"
        }
        di as text "{hline `c(linesize)'}"
        di as text "`hdr1'"
        di as text "`hdr2'"
        di as text "{hline `c(linesize)'}"

        // Print one row per DMU
        // Loop by observation index to correctly handle DMU names with spaces.
        forval k = 1/`=_N' {
            if `touse'[`k'] == 1 {
                local d_name = `dmu'[`k']
                local row `"`d_name'"'

                local ia = 0
                foreach a of numlist `alpha' {
                    local ++ia
                    local asuf = string(round(`a' * 100))
                    local vp "`generate'_pess_a`asuf'"
                    local vo "`generate'_opt_a`asuf'"

                    local ep_val = `vp'[`k']
                    local eo_val = `vo'[`k']

                    if `ep_val' == . local ep "     ."
                    else             local ep = string(`ep_val', "%6.4f")

                    if `eo_val' == . local eo "     ."
                    else             local eo = string(`eo_val', "%6.4f")

                    local row `"`row'  `ep'/`eo'"'
                }
                di as result `"`row'"'
            }
        }
        di as text "{hline `c(linesize)'}"
    }

    di as text "{hline `c(linesize)'}"

    // =========================================================================
    // Return values
    // =========================================================================
    return local cmd         "dea_fuzzy"
    return local orientation "`orientation'"
    return local rts         "`rts'"
    return local alpha       "`alpha'"
    return local generate    "`generate'"
    return scalar N          = `n_dmus'
    return scalar m_in       = `m_in'
    return scalar m_out      = `m_out'

end
