*! version 1.3.1  28apr2026
*! dea_malmquist.ado — Malmquist Productivity Index
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
*!   Caves, D. W., Christensen, L. R., & Diewert, W. E. (1982). The
*!     economic theory of index numbers and the measurement of input,
*!     output, and productivity. Econometrica, 50(6), 1393-1414.
*!     doi:10.2307/1913388
*!
*!   Fare, R., Grosskopf, S., Norris, M., & Zhang, Z. (1994).
*!     Productivity growth, technical progress, and efficiency change in
*!     industrialized countries. American Economic Review, 84(1), 66-83.
*!
*! Changelog:
*!   1.3.0  27apr2026  Add global option (global/pooled-frontier Malmquist);
*!                      implement generate() to save EC/TC/MPI per DMU to dataset;
*!                      add r(malmquist_type) return macro
*!   1.2.0  25apr2026  Fix string DMU: markout excludes dmu var; Mata uses
*!                      _dea_dmu_read() for string/numeric dispatch; display
*!                      shows original DMU labels (same fix as dea_crosseff 1.2.2)
*!   1.1.0  11feb2026  Unified syntax: supports both
*!                        dea_malmquist x = y, dmu() time() [options]  (Form A)
*!                        dea_malmquist, inputs(x) outputs(y) dmu() time() [options] (Form B)
*!                      ort() alias for orientation(); abbreviations i/o/in/out
*!   1.0.0  06feb2026  Initial packaged version; based on earlier teaching/seminar versions

/**
 * MALMQUIST PRODUCTIVITY INDEX
 * ============================
 *
 * Adjacent-period MPI (default):
 *   MPI = EC × TC
 *   EC = D^{t+1}(x^{t+1}, y^{t+1}) / D^t(x^t, y^t)       [catch-up]
 *   TC = sqrt[ D^t(x^{t+1}, y^{t+1})/D^{t+1}(x^{t+1}, y^{t+1})
 *            × D^t(x^t, y^t)/D^{t+1}(x^t, y^t) ]           [frontier shift]
 *
 * Global (pooled-frontier) MPI (option: global):
 *   Uses a single frontier estimated from all periods pooled.
 *   TFPCH = D^G(x^{t+1}, y^{t+1}) / D^G(x^t, y^t)         [TFP change]
 *   EFFCH = D^{t+1}(x^{t+1}, y^{t+1}) / D^t(x^t, y^t)     [efficiency change]
 *   TECHCH = TFPCH / EFFCH                                   [technical change]
 *   where D^G is distance to the global (all-period pooled) frontier.
 *   This index satisfies the circularity property; cross-period infeasibility
 *   problems under VRS are avoided because the global frontier is larger.
 *
 * Interpretation (both methods):
 *   Index > 1: improvement     Index = 1: no change     Index < 1: decline
 */

capture program drop dea_malmquist
program define dea_malmquist, rclass
    version 16.0

    // =========================================================================
    // Unified Syntax Parser
    //   Form A: dea_malmquist x1 x2 = y1 y2 [if] [in], dmu() time() [options]
    //   Form B: dea_malmquist [varlist] [if] [in], inputs() outputs() dmu() time() [options]
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

    // Note: [varlist] has no (numeric) constraint so string dmu in the
    // positional slot does not trigger r(109).
    syntax [varlist] [if] [in], ///
        TIME(varname) ///
        [DMU(varname)] ///
        [Inputs(varlist numeric) ///
         Outputs(varlist numeric) ///
         RTS(string) ///
         ORT(string) ORIentation(string) ///
         GLOBal ///
         GENerate(string) ///
         TOLerance(real 1e-7) ///
         DETail ///
         DEBUG]

    // Resolve outputs: positional varlist or outputs() option
    if "`outputs'" == "" {
        if "`varlist'" != "" {
            local outputs `varlist'
        }
        else {
            di as error "outputs() required, or use: dea_malmquist inputs = outputs, dmu() time() [options]"
            exit 198
        }
    }
    if "`inputs'" == "" {
        di as error "inputs() required, or use: dea_malmquist inputs = outputs, dmu() time() [options]"
        exit 198
    }

    // Reconcile ort() / orientation()
    if "`ort'" != "" & "`orientation'" != "" {
        di as error "Specify only one of: ort() or orientation()"
        exit 198
    }
    if "`ort'" != "" local orientation "`ort'"

    // Defaults
    if "`rts'"         == "" local rts "crs"
    if "`orientation'" == "" local orientation "input"
    if "`generate'"    == "" local generate "_malm"

    // Abbreviations for orientation
    if inlist("`orientation'", "i", "in")  local orientation "input"
    if inlist("`orientation'", "o", "out") local orientation "output"

    // Validate
    if !inlist("`rts'", "crs", "vrs") {
        di as error "Malmquist requires rts(crs) or rts(vrs)"
        exit 198
    }
    if !inlist("`orientation'", "input", "output") {
        di as error "orientation() must be input or output"
        exit 198
    }

    // Mark sample — do NOT include dmu in markout (string DMU causes r(109))
    marksample touse
    markout `touse' `inputs' `outputs' `time'

    // Auto-generate DMU identifier if not specified
    // Assign sequential numbers within each time period (assumes consistent
    // ordering across periods, same as dea auto-numbering convention)
    if "`dmu'" == "" {
        tempvar _dmu_auto
        qui bysort `time' (`touse'): gen long `_dmu_auto' = _n if `touse'
        local dmu "`_dmu_auto'"
    }

    qui replace `touse' = 0 if missing(`dmu')

    // Get unique time periods
    qui levelsof `time' if `touse', local(periods)
    local n_periods : word count `periods'

    if `n_periods' < 2 {
        di as error "Malmquist requires at least 2 time periods"
        exit 198
    }

    // Determine method label
    if "`global'" != "" {
        local malm_type "global"
        local malm_label "Global (Pooled-Frontier)"
    }
    else {
        local malm_type "adjacent"
        local malm_label "Adjacent Period"
    }

    // Prepare generate variables (drop if exist, then create)
    local gen_ec  "`generate'_ec"
    local gen_tc  "`generate'_tc"
    local gen_mpi "`generate'_mpi"
    foreach v in `gen_ec' `gen_tc' `gen_mpi' {
        capture drop `v'
        qui gen double `v' = .
    }

    // Display header
    di as text ""
    di as text "{hline 70}"
    di as text "Malmquist Productivity Index — `malm_label'"
    di as text "{hline 70}"
    di as text "Orientation: " as result "`orientation'"
    di as text "RTS:         " as result "`rts'"
    di as text "Periods:     " as result "`n_periods'"
    di as text "{hline 70}"

    // Get consecutive period pairs
    local prev_period ""
    local pair_count = 0

    foreach t of local periods {
        if "`prev_period'" != "" {
            local ++pair_count
            local period1_`pair_count' "`prev_period'"
            local period2_`pair_count' "`t'"
        }
        local prev_period "`t'"
    }

    // Process each consecutive pair
    forvalues p = 1/`pair_count' {
        local t0 "`period1_`p''"
        local t1 "`period2_`p''"

        di as text ""
        di as text "Period `t0' → `t1'"
        di as text "{hline 50}"

        if "`global'" != "" {
            mata: _dea_malmquist_global_compute("`inputs'", "`outputs'", ///
                "`dmu'", "`time'", "`touse'", ///
                `t0', `t1', "`orientation'", "`rts'", `tolerance', ///
                "`gen_ec'", "`gen_tc'", "`gen_mpi'", ///
                "`detail'" != "", "`debug'" != "")
        }
        else {
            mata: _dea_malmquist_compute("`inputs'", "`outputs'", ///
                "`dmu'", "`time'", "`touse'", ///
                `t0', `t1', "`orientation'", "`rts'", `tolerance', ///
                "`gen_ec'", "`gen_tc'", "`gen_mpi'", ///
                "`detail'" != "", "`debug'" != "")
        }
    }

    di as text ""
    di as text "{hline 70}"
    di as text "Index > 1: improvement   = 1: no change   < 1: decline"
    if "`global'" == "" {
        di as text "EC > 1: catch-up to frontier   TC > 1: frontier shift outward"
    }
    else {
        di as text "EFFCH > 1: catch-up   TECHCH > 1: frontier shift   TFPCH = EFFCH × TECHCH"
    }
    di as text "{hline 70}"

    // Capture Mata-stored r() scalars
    local _mean_ec  = r(mean_ec)
    local _mean_tc  = r(mean_tc)
    local _mean_mpi = r(mean_mpi)
    local _n_dmus   = r(n_dmus)

    // Return values
    return local malmquist_type "`malm_type'"
    return local orientation    "`orientation'"
    return local rts            "`rts'"
    return scalar n_periods  = `n_periods'
    return scalar mean_ec    = `_mean_ec'
    return scalar mean_tc    = `_mean_tc'
    return scalar mean_mpi   = `_mean_mpi'
    return scalar n_dmus     = `_n_dmus'
end


// ============================================================================
// MATA FUNCTIONS
// ============================================================================

mata:

// ─────────────────────────────────────────────────────────────────────────────
// _dea_dmu_read_malm()  — read DMU variable (string or numeric)
// ─────────────────────────────────────────────────────────────────────────────
void _dea_dmu_read_malm(string scalar  dmu_var,
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


// ─────────────────────────────────────────────────────────────────────────────
// _dea_malmquist_save()  — store EC/TC/MPI for one DMU into dataset
// Finds the dataset row where dmu==dmu_id and time==t1 (within touse)
// ─────────────────────────────────────────────────────────────────────────────
void _dea_malmquist_save(real colvector obs_all,
                         real colvector dmu_ids_all,
                         real colvector time_all,
                         real scalar    dmu_id,
                         real scalar    t1,
                         real scalar    ec,
                         real scalar    tc,
                         real scalar    mpi,
                         string scalar  var_ec,
                         string scalar  var_tc,
                         string scalar  var_mpi)
{
    real colvector match
    real scalar    row

    match = selectindex((dmu_ids_all :== dmu_id) :& (time_all :== t1))
    if (rows(match) == 0) return
    row = obs_all[match[1]]
    if (var_ec  != "") st_store(row, st_varindex(var_ec),  ec)
    if (var_tc  != "") st_store(row, st_varindex(var_tc),  tc)
    if (var_mpi != "") st_store(row, st_varindex(var_mpi), mpi)
}


// ─────────────────────────────────────────────────────────────────────────────
// _dea_malmquist_compute()  — Adjacent-period Malmquist
// ─────────────────────────────────────────────────────────────────────────────
void _dea_malmquist_compute(string scalar inputs,
                             string scalar outputs,
                             string scalar dmu_var,
                             string scalar time_var,
                             string scalar touse,
                             real scalar t0,
                             real scalar t1,
                             string scalar orientation,
                             string scalar rts,
                             real scalar tol,
                             string scalar var_ec,
                             string scalar var_tc,
                             string scalar var_mpi,
                             real scalar show_detail,
                             real scalar show_debug)
{
    real matrix X_all, Y_all, X_t0, Y_t0, X_t1, Y_t1
    real colvector time_all, dmu_ids, idx_t0, idx_t1
    real colvector dmu_t0, dmu_t1, common_ids, obs_all
    real colvector D_tt, D_t1t1, D_t_t1, D_t1_t, EC, TC, MPI
    string colvector dmu_labels, lab_t0, lab_t1, common_labels
    string rowvector input_names, output_names
    real rowvector x_t0, y_t0, x_t1, y_t1
    real scalar n0, n1, m, s, n_common, i, j, d
    real scalar dmu_id, idx0, idx1, mean_ec, mean_tc, mean_mpi

    input_names  = tokens(inputs)
    output_names = tokens(outputs)

    X_all    = st_data(., input_names,  touse)
    Y_all    = st_data(., output_names, touse)
    time_all = st_data(., time_var,     touse)
    obs_all  = selectindex(st_data(., touse))
    _dea_dmu_read_malm(dmu_var, touse, dmu_ids, dmu_labels)

    idx_t0 = selectindex(time_all :== t0)
    idx_t1 = selectindex(time_all :== t1)

    X_t0   = X_all[idx_t0, .];  Y_t0   = Y_all[idx_t0, .]
    dmu_t0 = dmu_ids[idx_t0];   lab_t0 = dmu_labels[idx_t0]

    X_t1   = X_all[idx_t1, .];  Y_t1   = Y_all[idx_t1, .]
    dmu_t1 = dmu_ids[idx_t1];   lab_t1 = dmu_labels[idx_t1]

    n0 = rows(X_t0);  n1 = rows(X_t1)
    m  = cols(X_t0);  s  = cols(Y_t0)

    // Find common DMUs
    common_ids    = J(0, 1, .)
    common_labels = J(0, 1, "")
    for (i = 1; i <= n0; i++) {
        d = dmu_t0[i]
        for (j = 1; j <= n1; j++) {
            if (dmu_t1[j] == d) {
                common_ids    = common_ids    \ d
                common_labels = common_labels \ lab_t0[i]
                break
            }
        }
    }
    n_common = rows(common_ids)
    if (n_common == 0) { printf("  Warning: No common DMUs between periods\n"); return; }

    D_tt   = J(n_common, 1, .)
    D_t1t1 = J(n_common, 1, .)
    D_t_t1 = J(n_common, 1, .)
    D_t1_t = J(n_common, 1, .)

    for (i = 1; i <= n_common; i++) {
        dmu_id = common_ids[i]
        idx0 = selectindex(dmu_t0 :== dmu_id)[1]
        idx1 = selectindex(dmu_t1 :== dmu_id)[1]
        x_t0 = X_t0[idx0, .];  y_t0 = Y_t0[idx0, .]
        x_t1 = X_t1[idx1, .];  y_t1 = Y_t1[idx1, .]
        D_tt[i]   = _dea_malm_eff(x_t0, y_t0, X_t0, Y_t0, orientation, rts, tol)
        D_t1t1[i] = _dea_malm_eff(x_t1, y_t1, X_t1, Y_t1, orientation, rts, tol)
        D_t_t1[i] = _dea_malm_eff(x_t1, y_t1, X_t0, Y_t0, orientation, rts, tol)
        D_t1_t[i] = _dea_malm_eff(x_t0, y_t0, X_t1, Y_t1, orientation, rts, tol)
    }

    EC  = D_t1t1 :/ D_tt
    TC  = sqrt((D_t_t1 :/ D_t1t1) :* (D_tt :/ D_t1_t))
    MPI = EC :* TC

    // Display
    printf("\n")
    printf("%15s | %10s %10s %10s %10s | %9s %9s %9s\n",
           "DMU", "D_tt", "D_t1t1", "D_t_t1", "D_t1_t", "EC", "TC", "MPI")
    printf("%s\n", "-"*88)
    for (i = 1; i <= n_common; i++) {
        printf("%15s | %10.4f %10.4f %10.4f %10.4f | %9.4f %9.4f %9.4f\n",
               common_labels[i], D_tt[i], D_t1t1[i], D_t_t1[i], D_t1_t[i],
               EC[i], TC[i], MPI[i])
        // Save to dataset
        _dea_malmquist_save(obs_all, dmu_ids, time_all,
                            common_ids[i], t1,
                            EC[i], TC[i], MPI[i],
                            var_ec, var_tc, var_mpi)
    }
    printf("%s\n", "-"*88)

    mean_ec  = mean(EC)
    mean_tc  = mean(TC)
    mean_mpi = mean(MPI)
    printf("%15s | %10s %10s %10s %10s | %9.4f %9.4f %9.4f\n",
           "Mean", "", "", "", "", mean_ec, mean_tc, mean_mpi)

    st_numscalar("r(mean_ec)",  mean_ec)
    st_numscalar("r(mean_tc)",  mean_tc)
    st_numscalar("r(mean_mpi)", mean_mpi)
    st_numscalar("r(n_dmus)",   n_common)
}


// ─────────────────────────────────────────────────────────────────────────────
// _dea_malmquist_global_compute()  — Global (pooled-frontier) Malmquist
// ─────────────────────────────────────────────────────────────────────────────
void _dea_malmquist_global_compute(string scalar inputs,
                                    string scalar outputs,
                                    string scalar dmu_var,
                                    string scalar time_var,
                                    string scalar touse,
                                    real scalar t0,
                                    real scalar t1,
                                    string scalar orientation,
                                    string scalar rts,
                                    real scalar tol,
                                    string scalar var_ec,
                                    string scalar var_tc,
                                    string scalar var_mpi,
                                    real scalar show_detail,
                                    real scalar show_debug)
{
    real matrix X_all, Y_all, X_t0, Y_t0, X_t1, Y_t1
    real colvector time_all, dmu_ids, idx_t0, idx_t1
    real colvector dmu_t0, dmu_t1, common_ids, obs_all
    real colvector TFPCH, EFFCH, TECHCH
    string colvector dmu_labels, lab_t0, common_labels
    string rowvector input_names, output_names
    real rowvector x_t0, y_t0, x_t1, y_t1
    real scalar n0, n1, n_common, i, j, d
    real scalar dmu_id, idx0, idx1
    real scalar d_t_g, d_t1_g, d_tt, d_t1t1
    real scalar mean_ec, mean_tc, mean_mpi

    input_names  = tokens(inputs)
    output_names = tokens(outputs)

    X_all    = st_data(., input_names,  touse)
    Y_all    = st_data(., output_names, touse)
    time_all = st_data(., time_var,     touse)
    obs_all  = selectindex(st_data(., touse))
    _dea_dmu_read_malm(dmu_var, touse, dmu_ids, dmu_labels)

    idx_t0 = selectindex(time_all :== t0)
    idx_t1 = selectindex(time_all :== t1)

    X_t0   = X_all[idx_t0, .];  Y_t0   = Y_all[idx_t0, .]
    X_t1   = X_all[idx_t1, .];  Y_t1   = Y_all[idx_t1, .]
    dmu_t0 = dmu_ids[idx_t0]
    dmu_t1 = dmu_ids[idx_t1]
    lab_t0 = dmu_labels[idx_t0]

    n0 = rows(X_t0);  n1 = rows(X_t1)

    // Find common DMUs
    common_ids    = J(0, 1, .)
    common_labels = J(0, 1, "")
    for (i = 1; i <= n0; i++) {
        d = dmu_t0[i]
        for (j = 1; j <= n1; j++) {
            if (dmu_t1[j] == d) {
                common_ids    = common_ids    \ d
                common_labels = common_labels \ lab_t0[i]
                break
            }
        }
    }
    n_common = rows(common_ids)
    if (n_common == 0) { printf("  Warning: No common DMUs between periods\n"); return; }

    TFPCH  = J(n_common, 1, .)
    EFFCH  = J(n_common, 1, .)
    TECHCH = J(n_common, 1, .)

    for (i = 1; i <= n_common; i++) {
        dmu_id = common_ids[i]
        idx0 = selectindex(dmu_t0 :== dmu_id)[1]
        idx1 = selectindex(dmu_t1 :== dmu_id)[1]
        x_t0 = X_t0[idx0, .];  y_t0 = Y_t0[idx0, .]
        x_t1 = X_t1[idx1, .];  y_t1 = Y_t1[idx1, .]

        // Distance to global (pooled) frontier
        d_t_g  = _dea_malm_eff(x_t0, y_t0, X_all, Y_all, orientation, rts, tol)
        d_t1_g = _dea_malm_eff(x_t1, y_t1, X_all, Y_all, orientation, rts, tol)

        // Distance to own-period frontiers
        d_tt   = _dea_malm_eff(x_t0, y_t0, X_t0, Y_t0, orientation, rts, tol)
        d_t1t1 = _dea_malm_eff(x_t1, y_t1, X_t1, Y_t1, orientation, rts, tol)

        if (d_t_g > 0 & d_t1_g > 0 & d_tt > 0 & d_t1t1 > 0) {
            TFPCH[i]  = d_t1_g  / d_t_g
            EFFCH[i]  = d_t1t1  / d_tt
            TECHCH[i] = TFPCH[i] / EFFCH[i]
        }
    }

    // Display (column labels: TFPCH = MPI for global; EFFCH = EC; TECHCH = TC)
    printf("\n")
    printf("%15s | %9s %9s %9s\n", "DMU", "EFFCH", "TECHCH", "TFPCH")
    printf("%s\n", "-"*50)
    for (i = 1; i <= n_common; i++) {
        printf("%15s | %9.4f %9.4f %9.4f\n",
               common_labels[i], EFFCH[i], TECHCH[i], TFPCH[i])
        _dea_malmquist_save(obs_all, dmu_ids, time_all,
                            common_ids[i], t1,
                            EFFCH[i], TECHCH[i], TFPCH[i],
                            var_ec, var_tc, var_mpi)
    }
    printf("%s\n", "-"*50)

    mean_ec  = mean(EFFCH)
    mean_tc  = mean(TECHCH)
    mean_mpi = mean(TFPCH)
    printf("%15s | %9.4f %9.4f %9.4f\n", "Mean", mean_ec, mean_tc, mean_mpi)

    st_numscalar("r(mean_ec)",  mean_ec)
    st_numscalar("r(mean_tc)",  mean_tc)
    st_numscalar("r(mean_mpi)", mean_mpi)
    st_numscalar("r(n_dmus)",   n_common)
}


// ─────────────────────────────────────────────────────────────────────────────
// _dea_malm_eff()  — DEA distance for one DMU against a reference set
// ─────────────────────────────────────────────────────────────────────────────
real scalar _dea_malm_eff(real rowvector x0,
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
        if (rts == "vrs") { A_eq = (J(1, n, 1), 0); b_eq = 1; }
    }
    else {
        c      = (J(1, n, 0), -1)
        A_ineq = (X_ref', J(m, 1, 0)) \ (-Y_ref', y0')
        b_ineq = x0' \ J(s, 1, 0)
        // Cross-period Malmquist distances can evaluate observations outside
        // the reference frontier, so the output expansion factor may be < 1.
        lb     = J(1, n + 1, 0)
        ub     = J(1, n + 1, .)
        A_eq   = J(0, n + 1, .)
        b_eq   = J(0, 1, .)
        if (rts == "vrs") { A_eq = (J(1, n, 1), 0); b_eq = 1; }
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
        eff = (orientation == "input") ? S.value() : 1 / (-S.value())
    }
    else {
        eff = .
    }

    return(eff)
}

end
