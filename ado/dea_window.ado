*! version 1.3.2  20jul2026
*! dea_window.ado — Window Analysis DEA
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
*!   Charnes, A., Cooper, W. W., Golany, B., Seiford, L., & Stutz, J.
*!     (1985). Foundations of data envelopment analysis for
*!     Pareto-Koopmans efficient empirical production functions.
*!     Journal of Econometrics, 30(1-2), 91-107.
*!     doi:10.1016/0304-4076(85)90133-2
*!
*!   Cooper, W. W., Seiford, L. M., & Tone, K. (2007). Data Envelopment
*!     Analysis: A Comprehensive Text with Models, Applications, References
*!     and DEA-Solver Software (2nd ed.). New York: Springer.
*!     doi:10.1007/978-0-387-45283-8
*!
*! Changelog:
*!   1.3.2  20jul2026  Fix inputs()/outputs() sample marking by explicitly
*!                      populating the parser varlist before marksample
*!   1.3.0  02may2026  Implement generate(): save window scores as stub_w#_p#
*!                      plus stub_mean and stub_sd
*!   1.2.3  02may2026  Warn that generate() is accepted but not saved to dataset
*!   1.2.2  28apr2026  Fix r(mean_eff)/r(sd_eff)/r(n_obs): propagate Mata
*!                      st_numscalar results via explicit return scalar statements
*!   1.2.1  28apr2026  dmu() optional: auto-generate sequential IDs if omitted
*!   1.2.0  25apr2026  Fix string DMU: markout excludes dmu var; Mata uses
*!                      _dea_dmu_read_win() for string/numeric dispatch; display
*!                      shows original DMU labels (same fix as dea_crosseff 1.2.2)
*!   1.1.0  11feb2026  Unified syntax: dea_window x = y | inputs()/outputs()
*!                      ort() alias for orientation(); abbreviations i/o/in/out
*!   1.0.0  06feb2026  Initial packaged version; based on earlier teaching/seminar versions

/**
 * WINDOW ANALYSIS DEA
 * ===================
 * 
 * Window Analysis tracks DMU performance over time using moving windows.
 * 
 * Key features:
 *   - Pools data from multiple periods within each window
 *   - Each DMU-period is evaluated against the pooled reference set
 *   - Allows tracking performance trends and stability
 * 
 * Example with 3 periods and window width = 2:
 *   Window 1: periods 1-2 (6 DMU-period observations)
 *   Window 2: periods 2-3 (6 DMU-period observations)
 * 
 * Each DMU appears multiple times in each window (once per period).
 */

capture program drop dea_window
program define dea_window, rclass
    version 16.0
    
    // =========================================================================
    // Unified Syntax Parser
    //   Form A: dea_window x1 x2 = y1 y2 [if] [in], dmu() time() [options]
    //   Form B: dea_window [varlist] [if] [in], inputs() outputs() dmu() time() [options]
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
    
    syntax [varlist] [if] [in], ///
        TIME(varname) ///
        [DMU(varname)] ///
        [Inputs(varlist numeric) ///
         Outputs(varlist numeric) ///
         RTS(string) ///
         ORT(string) ORIentation(string) ///
         Width(integer 0) ///
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
            di as error "outputs() required, or use: dea_window inputs = outputs, dmu() time() [options]"
            exit 198
        }
    }
    if "`inputs'" == "" {
        di as error "inputs() required, or use: dea_window inputs = outputs, dmu() time() [options]"
        exit 198
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
    if "`generate'" == "" local generate "_win"
    
    // Abbreviations
    if inlist("`orientation'", "i", "in")  local orientation "input"
    if inlist("`orientation'", "o", "out") local orientation "output"
    if "`rts'" == "c" local rts "crs"
    if "`rts'" == "v" local rts "vrs"
    
    // Validate
    if !inlist("`rts'", "crs", "vrs", "drs", "irs") {
        di as error "rts() must be crs, vrs, drs, or irs"
        exit 198
    }
    if !inlist("`orientation'", "input", "output") {
        di as error "orientation() must be input or output"
        exit 198
    }
    
    // Mark sample — do NOT include dmu in markout (string DMU causes r(109))
    // marksample reads local varlist.  Populate it explicitly because the
    // inputs()/outputs() syntax form has no main varlist.
    local varlist `inputs' `outputs' `time'
    marksample touse
    markout `touse' `inputs' `outputs' `time'

    // Auto-generate DMU identifier if not specified
    if "`dmu'" == "" {
        tempvar _dmu_auto
        qui bysort `time' (`touse'): gen long `_dmu_auto' = _n if `touse'
        local dmu "`_dmu_auto'"
    }

    qui replace `touse' = 0 if missing(`dmu')   // handles string("") and numeric(.)

    // Get unique time periods
    qui levelsof `time' if `touse', local(periods)
    local n_periods : word count `periods'
    
    if `n_periods' < 2 {
        di as error "Window analysis requires at least 2 time periods"
        exit 198
    }
    
    // Default window width
    if `width' == 0 {
        local width = min(3, `n_periods')
    }
    if `width' < 2 {
        di as error "Window width must be at least 2"
        exit 198
    }
    if `width' > `n_periods' {
        di as error "Window width cannot exceed number of periods"
        exit 198
    }
    
    // Get unique DMUs
    qui levelsof `dmu' if `touse', local(dmus)
    local n_dmus : word count `dmus'
    
    // Display header
    di as text ""
    di as text "{hline 70}"
    di as text "Window Analysis DEA"
    di as text "{hline 70}"
    di as text "Orientation:   " as result "`orientation'"
    di as text "RTS:           " as result "`rts'"
    di as text "Periods:       " as result "`n_periods'"
    di as text "DMUs:          " as result "`n_dmus'"
    di as text "Window width:  " as result "`width'"
    di as text "{hline 70}"
    
    // Calculate number of windows
    local n_windows = `n_periods' - `width' + 1
    di as text "Number of windows: " as result "`n_windows'"
    di as text "{hline 70}"

    // Prepare generated variables. Scores are stored by window and by the
    // observation's position within that window.
    forvalues w = 1/`n_windows' {
        forvalues p = 1/`width' {
            capture drop `generate'_w`w'_p`p'
            qui gen double `generate'_w`w'_p`p' = .
        }
    }
    capture drop `generate'_mean
    capture drop `generate'_sd
    qui gen double `generate'_mean = .
    qui gen double `generate'_sd = .
    
    // Create period list as numeric
    local period_list ""
    foreach p of local periods {
        local period_list "`period_list' `p'"
    }
    
    // Call Mata function
    mata: _dea_window_compute("`inputs'", "`outputs'", ///
        "`dmu'", "`time'", "`touse'", ///
        "`period_list'", `n_periods', `width', ///
        "`orientation'", "`rts'", "`generate'", `tolerance', ///
        "`detail'" != "", "`debug'" != "")
    
    di as text ""
    di as text "{hline 70}"

    // Return values — capture Mata-set r() scalars before rclass posts results
    local _mean_eff = r(mean_eff)
    local _sd_eff   = r(sd_eff)
    local _n_obs    = r(n_obs)

    return local orientation "`orientation'"
    return local rts "`rts'"
    return scalar n_periods = `n_periods'
    return scalar n_dmus = `n_dmus'
    return scalar width = `width'
    return scalar n_windows = `n_windows'
    return scalar mean_eff = `_mean_eff'
    return scalar sd_eff   = `_sd_eff'
    return scalar n_obs    = `_n_obs'
end


// ============================================================================
// MATA FUNCTIONS
// ============================================================================

capture mata: mata drop _dea_dmu_read_win()
capture mata: mata drop _dea_window_compute()
capture mata: mata drop _dea_window_efficiency()

mata:

// ─────────────────────────────────────────────────────────────────────────────
// _dea_dmu_read_win()
// Reads a DMU variable that may be string or numeric.
// Returns: dmu_ids    (numeric: 1-indexed for strings, original for numeric)
//          dmu_labels (string: original labels for display)
// ─────────────────────────────────────────────────────────────────────────────
void _dea_dmu_read_win(string scalar   dmu_var,
                        string scalar   touse,
                        real colvector  dmu_ids,
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
 * _dea_window_compute()
 *
 * Computes Window Analysis DEA for all windows.
 */
void _dea_window_compute(string scalar inputs,
                          string scalar outputs,
                          string scalar dmu_var,
                          string scalar time_var,
                          string scalar touse,
                          string scalar period_list,
                          real scalar n_periods,
                          real scalar width,
                          string scalar orientation,
                          string scalar rts,
                          string scalar gen_prefix,
                          real scalar tol,
                          real scalar show_detail,
                          real scalar show_debug)
{
    // ── All declarations at top ───────────────────────────────────────────────
    real matrix X_all, Y_all, X_win, Y_win, results
    real colvector time_all, dmu_ids, unique_ids, in_window, idx_window, obs_all
    real colvector dmu_win, time_win, row_vals, all_eff
    string colvector dmu_labels, unique_labels
    string rowvector input_names, output_names
    real rowvector periods, x0, y0
    real scalar n_all, m, s, n_dmus, n_windows, n_win
    real scalar w, p, i, k, pos, ii, p_in_window, obs_row
    real scalar start_period, end_period, dmu_id, period_id, eff
    real scalar lw, row_mean, row_sd
    string scalar lw_fmt, lw_sep, dmu_lab, score_var

    // ── Load data ─────────────────────────────────────────────────────────────
    input_names  = tokens(inputs)
    output_names = tokens(outputs)
    periods      = strtoreal(tokens(period_list))

    X_all    = st_data(., input_names, touse)
    Y_all    = st_data(., output_names, touse)
    time_all = st_data(., time_var, touse)
    obs_all  = selectindex(st_data(., touse))
    _dea_dmu_read_win(dmu_var, touse, dmu_ids, dmu_labels)

    n_all = rows(X_all)
    m     = cols(X_all)
    s     = cols(Y_all)

    // ── Unique DMUs and labels ────────────────────────────────────────────────
    unique_ids    = sort(uniqrows(dmu_ids), 1)
    n_dmus        = rows(unique_ids)
    unique_labels = J(n_dmus, 1, "")
    for (ii = 1; ii <= rows(dmu_ids); ii++) {
        for (pos = 1; pos <= n_dmus; pos++) {
            if (dmu_ids[ii] == unique_ids[pos]) {
                unique_labels[pos] = dmu_labels[ii]
                break
            }
        }
    }

    // ── Window setup ──────────────────────────────────────────────────────────
    n_windows = n_periods - width + 1
    results   = J(n_dmus * n_periods, 2 + n_windows, .)

    // Fill DMU and period columns
    k = 0
    for (i = 1; i <= n_dmus; i++) {
        for (p = 1; p <= n_periods; p++) {
            k++
            results[k, 1] = unique_ids[i]
            results[k, 2] = periods[p]
        }
    }

    // ── Process each window ───────────────────────────────────────────────────
    for (w = 1; w <= n_windows; w++) {
        start_period = periods[w]
        end_period   = periods[w + width - 1]

        if (show_debug) {
            printf("\nWindow %g: periods %g to %g\n", w, start_period, end_period)
        }

        in_window  = (time_all :>= start_period) :& (time_all :<= end_period)
        idx_window = selectindex(in_window)

        X_win    = X_all[idx_window, .]
        Y_win    = Y_all[idx_window, .]
        dmu_win  = dmu_ids[idx_window]
        time_win = time_all[idx_window]
        n_win    = rows(X_win)

        if (show_debug) {
            printf("  Observations in window: %g\n", n_win)
        }

        for (i = 1; i <= n_win; i++) {
            dmu_id    = dmu_win[i]
            period_id = time_win[i]
            x0 = X_win[i, .]
            y0 = Y_win[i, .]

            eff = _dea_window_efficiency(x0, y0, X_win, Y_win,
                                          orientation, rts, tol)
            p_in_window = selectindex(periods :== period_id)[1] - w + 1
            obs_row = obs_all[idx_window[i]]
            score_var = gen_prefix + "_w" + strofreal(w) + "_p" + strofreal(p_in_window)
            st_store(obs_row, st_varindex(score_var), eff)

            for (k = 1; k <= rows(results); k++) {
                if (results[k, 1] == dmu_id && results[k, 2] == period_id) {
                    results[k, 2 + w] = eff
                    break
                }
            }
        }
    }

    // ── Display ───────────────────────────────────────────────────────────────
    lw = 10
    for (i = 1; i <= n_dmus; i++) {
        if (strlen(unique_labels[i]) > lw) lw = strlen(unique_labels[i])
    }
    lw     = lw + 1
    lw_fmt = "%" + strofreal(lw) + "s"
    lw_sep = "-" * lw

    printf("\n")
    printf("Efficiency Scores by Window\n")
    printf("{hline 70}\n")

    printf(lw_fmt + " %8s |", "DMU", "Period")
    for (w = 1; w <= n_windows; w++) printf(" Window%g", w)
    printf(" |   Mean     SD\n")

    printf(lw_fmt + " %8s |", lw_sep, "--------")
    for (w = 1; w <= n_windows; w++) printf(" -------")
    printf(" | ------- ------\n")

    for (k = 1; k <= rows(results); k++) {
        row_vals = J(0, 1, .)
        for (w = 1; w <= n_windows; w++) {
            if (results[k, 2 + w] != .) row_vals = row_vals \ results[k, 2 + w]
        }

        if (rows(row_vals) > 0) {
            row_mean = mean(row_vals)
            row_sd   = (rows(row_vals) > 1 ? sqrt(variance(row_vals)) : 0)

            dmu_lab = strofreal(results[k, 1])
            for (pos = 1; pos <= n_dmus; pos++) {
                if (unique_ids[pos] == results[k, 1]) {
                    dmu_lab = unique_labels[pos]
                    break
                }
            }

            printf(lw_fmt + " %8.0f |", dmu_lab, results[k, 2])
            for (w = 1; w <= n_windows; w++) {
                if (results[k, 2 + w] != .) printf(" %7.4f", results[k, 2 + w])
                else                         printf(" %7s", ".")
            }
            printf(" | %7.4f %6.4f\n", row_mean, row_sd)

            for (ii = 1; ii <= rows(dmu_ids); ii++) {
                if (dmu_ids[ii] == results[k, 1] && time_all[ii] == results[k, 2]) {
                    st_store(obs_all[ii], st_varindex(gen_prefix + "_mean"), row_mean)
                    st_store(obs_all[ii], st_varindex(gen_prefix + "_sd"), row_sd)
                    break
                }
            }
        }
    }
    printf("{hline 70}\n")

    printf("\nSummary by DMU (Mean across all windows)\n")
    printf("{hline 50}\n")
    printf(lw_fmt + " | %10s %10s\n", "DMU", "Mean", "SD")
    printf(lw_fmt + " | %10s %10s\n", lw_sep, "----------", "----------")

    for (i = 1; i <= n_dmus; i++) {
        row_vals = J(0, 1, .)
        for (k = 1; k <= rows(results); k++) {
            if (results[k, 1] == unique_ids[i]) {
                for (w = 1; w <= n_windows; w++) {
                    if (results[k, 2 + w] != .) row_vals = row_vals \ results[k, 2 + w]
                }
            }
        }
        if (rows(row_vals) > 0) {
            row_mean = mean(row_vals)
            row_sd   = (rows(row_vals) > 1 ? sqrt(variance(row_vals)) : 0)
            printf(lw_fmt + " | %10.4f %10.4f\n", unique_labels[i], row_mean, row_sd)
        }
    }
    printf("{hline 50}\n")

    // ── Store results ─────────────────────────────────────────────────────────
    all_eff = J(0, 1, .)
    for (k = 1; k <= rows(results); k++) {
        for (w = 1; w <= n_windows; w++) {
            if (results[k, 2 + w] != .) all_eff = all_eff \ results[k, 2 + w]
        }
    }

    st_numscalar("r(mean_eff)", mean(all_eff))
    st_numscalar("r(sd_eff)",   sqrt(variance(all_eff)))
    st_numscalar("r(n_obs)",    rows(all_eff))
}


/**
 * _dea_window_efficiency()
 *
 * Computes DEA efficiency score for a single DMU against window reference set.
 */
real scalar _dea_window_efficiency(real rowvector x0,
                                    real rowvector y0,
                                    real matrix X_ref,
                                    real matrix Y_ref,
                                    string scalar orientation,
                                    string scalar rts,
                                    real scalar tol)
{
    real scalar n, m, s, eff
    real rowvector c, lb, ub
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    class LinearProgram scalar S

    n = rows(X_ref)
    m = cols(X_ref)
    s = cols(Y_ref)

    if (orientation == "input") {
        c      = (J(1, n, 0), 1)
        A_ineq = (X_ref', -x0') \ (-Y_ref', J(s, 1, 0))
        b_ineq = J(m, 1, 0) \ (-y0')
        lb     = J(1, n + 1, 0)
        ub     = J(1, n + 1, .)
        A_eq   = J(0, n + 1, .)
        b_eq   = J(0, 1, .)
        if (rts == "vrs") {
            A_eq = (J(1, n, 1), 0)
            b_eq = 1
        }
        else if (rts == "drs") {
            A_ineq = A_ineq \ (J(1,n,1),0)
            b_ineq = b_ineq \ 1
        }
        else if (rts == "irs") {
            A_ineq = A_ineq \ (J(1,n,-1),0)
            b_ineq = b_ineq \ (-1)
        }
    }
    else {
        c      = (J(1, n, 0), -1)
        A_ineq = (X_ref', J(m, 1, 0)) \ (-Y_ref', y0')
        b_ineq = x0' \ J(s, 1, 0)
        lb     = (J(1, n, 0), 1)
        ub     = J(1, n + 1, .)
        A_eq   = J(0, n + 1, .)
        b_eq   = J(0, 1, .)
        if (rts == "vrs") {
            A_eq = (J(1, n, 1), 0)
            b_eq = 1
        }
        else if (rts == "drs") {
            A_ineq = A_ineq \ (J(1,n,1),0)
            b_ineq = b_ineq \ 1
        }
        else if (rts == "irs") {
            A_ineq = A_ineq \ (J(1,n,-1),0)
            b_ineq = b_ineq \ (-1)
        }
    }

    S = LinearProgram()
    S.setMaxOrMin("min")
    S.setCoefficients(c)
    S.setInequality(A_ineq, b_ineq)
    if (rows(A_eq) > 0) S.setEquality(A_eq, b_eq)
    S.setBounds(lb, ub)
    S.setTrace("off")
    (void) S.optimize()

    if (S.converged() && S.errorcode() == 0) {
        eff = (orientation == "input" ? S.value() : 1 / (-S.value()))
    }
    else {
        eff = .
    }

    return(eff)
}

end
