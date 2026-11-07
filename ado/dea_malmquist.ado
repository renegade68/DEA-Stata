*! version 1.3.7  02nov2026
*! 2026-11-02: honor tolerance, validate radial LPs, and report panel failures.
*! 2026-11-01: number selected rows in stable within-period order when dmu() is omitted.
*! 2026-11-01: reject repeated selected DMU-period keys before replacing results.
*! 2026-11-01: accept the minimum option abbreviations advertised in help.
*! 2026-11-01: use shared Form A/B parser; preserve if/in expressions.
*! 2026-11-01: preflight output names; explicit replacevars; protect model data.
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
 * E denotes theta for input orientation and 1/phi for output orientation.
 * The D notation in the formulas below refers to E, not the input Shephard
 * distance (which is 1/theta). Adjacent-period MPI (default):
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
program define dea_malmquist, rclass sortpreserve
    version 16.0

    // =========================================================================
    // Unified Syntax Parser
    //   Form A: dea_malmquist x1 x2 = y1 y2 [if] [in], dmu() time() [options]
    //   Form B: dea_malmquist [varlist] [if] [in], inputs() outputs() dmu() time() [options]
    // =========================================================================
    _dea_parse_io `0'
    local 0 `"`r(args)'"'

    syntax [varlist(numeric default=none)] [if] [in], ///
        TIME(varname) ///
        [DMU(varname)] ///
        [Inputs(varlist numeric) ///
         Outputs(varlist numeric) ///
         Rts(string) ///
         ORT(string) ORientation(string) ///
         GLOBal ///
         GENerate(string) ///
         REPLACEVARS      ///
         TOLerance(real 1e-7) ///
         DETail ///
         DEBUG]

    if missing(`tolerance') | `tolerance'<=0 {
        di as error "tolerance() must be positive and finite"
        exit 198
    }
    capture mata: assert(_dea_nonte_lp_ready()==1)
    if _rc {
        quietly findfile _dea_nonte_lp.ado
        quietly run "`r(fn)'"
    }


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

    confirm numeric variable `time'

    // Validate every output before any existing variable can be replaced.
    _dea_output_check, names(`generate'_ec `generate'_tc `generate'_mpi) ///
        protect(`inputs' `outputs' `dmu' `time') `replacevars'
    local _replace_outputs "`r(existing)'"

    // Mark sample — do NOT include dmu in markout (string DMU causes r(109))
    marksample touse, novarlist
    markout `touse' `inputs' `outputs' `time'

    // Positional IDs use only selected observations in their incoming order.
    // Restore the row order before solving so if/in and output rows stay aligned.
    local _dmu_method "variable"
    if "`dmu'" == "" {
        tempvar _dmu_auto _panel_order
        qui gen double `_panel_order' = _n
        qui bysort `time' (`_panel_order'): gen long `_dmu_auto' = sum(`touse') if `touse'
        qui sort `_panel_order'
        local dmu "`_dmu_auto'"
        local _dmu_method "position"
        di as text "Note: dmu() omitted; periods are matched by selected within-period row order."
    }

    qui replace `touse' = 0 if missing(`dmu')

    // Validate only the final selected sample, before any output replacement.
    foreach _side in inputs outputs {
        tempvar _positive_side
        quietly generate byte `_positive_side'=0 if `touse'
        foreach _v of local `_side' {
            quietly count if `touse' & `_v'<0
            if r(N) {
                di as error "DEA panel quantities must be nonnegative"
                exit 459
            }
            quietly replace `_positive_side'=1 if `touse' & `_v'>0
        }
        quietly count if `touse' & `_positive_side'==0
        if r(N) {
            di as error "each selected DMU-period requires a positive input and output"
            exit 459
        }
    }

    _dea_panel_check, dmu(`dmu') time(`time') touse(`touse')

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

    // The sample may refer to an existing output; mark it before replacement.
    foreach _ov of local _replace_outputs {
        quietly drop `_ov'
    }

    // Create the validated result columns
    local gen_ec  "`generate'_ec"
    local gen_tc  "`generate'_tc"
    local gen_mpi "`generate'_mpi"
    foreach v in `gen_ec' `gen_tc' `gen_mpi' {
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

    tempname _pair_stats
    matrix `_pair_stats'=J(`pair_count',10,.)
    matrix colnames `_pair_stats'=period0 period1 N n_failed n_failed_ec n_failed_tc n_failed_mpi mean_ec mean_tc mean_mpi
    local _failed_total=0

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
        matrix `_pair_stats'[`p',1]=(`t0',`t1',r(n_dmus),r(n_failed),r(n_failed_ec),r(n_failed_tc),r(n_failed_mpi),r(mean_ec),r(mean_tc),r(mean_mpi))
        local _failed_total=`_failed_total'+r(n_failed)

    }

    di as text ""
    di as text "{hline 70}"
    di as text "Index > 1: improvement   = 1: no change   < 1: decline"
    if "`global'" == "" {
        di as text "EC > 1: catch-up to frontier   TC > 1: frontier shift outward"
    }
    else {
        di as text "EFFCH > 1: catch-up   TECHCH > 1: best-practice gap improvement   TFPCH = EFFCH × TECHCH"
    }
    di as text "{hline 70}"

    // Capture Mata-stored r() scalars
    local _mean_ec  = r(mean_ec)
    local _mean_tc  = r(mean_tc)
    local _mean_mpi = r(mean_mpi)
    local _n_dmus   = r(n_dmus)

    local _n_failed=r(n_failed)
    local _n_failed_ec=r(n_failed_ec)
    local _n_failed_tc=r(n_failed_tc)
    local _n_failed_mpi=r(n_failed_mpi)

    // Return values
    return local malmquist_type "`malm_type'"
    return local dmu_method "`_dmu_method'"
    return local orientation    "`orientation'"
    return local rts            "`rts'"
    return matrix pair_summary=`_pair_stats'
    return scalar n_pairs=`pair_count'
    return scalar n_failed=`_n_failed'
    return scalar n_failed_ec=`_n_failed_ec'
    return scalar n_failed_tc=`_n_failed_tc'
    return scalar n_failed_mpi=`_n_failed_mpi'
    return scalar n_failed_total=`_failed_total'
    return scalar tolerance=`tolerance'
    return scalar n_periods  = `n_periods'
    return scalar mean_ec    = `_mean_ec'
    return scalar mean_tc    = `_mean_tc'
    return scalar mean_mpi   = `_mean_mpi'
    return scalar n_dmus     = `_n_dmus'
end


// ============================================================================
// MATA FUNCTIONS
// ============================================================================

capture mata: mata drop _dea_dmu_read_malm()
capture mata: mata drop _dea_malmquist_save()
capture mata: mata drop _dea_malmquist_compute()
capture mata: mata drop _dea_malmquist_global_compute()
capture mata: mata drop _dea_malm_eff()
capture mata: mata drop _dea_malm_mean()
capture mata: mata drop _dea_malm_stats()
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
        dmu_labels = strtrim(strofreal(dmu_ids,"%21.17g"))
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
    if (n_common == 0) {
        printf("  Warning: No common DMUs between periods\n")
        _dea_malm_stats(J(0,1,.),J(0,1,.),J(0,1,.))
        return
    }

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
           "DMU", "E_tt", "E_t1t1", "E_t_t1", "E_t1_t", "EC", "TC", "MPI")
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

    mean_ec  = _dea_malm_mean(EC)
    mean_tc  = _dea_malm_mean(TC)
    mean_mpi = _dea_malm_mean(MPI)
    printf("%15s | %10s %10s %10s %10s | %9.4f %9.4f %9.4f\n",
           "Mean", "", "", "", "", mean_ec, mean_tc, mean_mpi)

    _dea_malm_stats(EC,TC,MPI)
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
    if (n_common == 0) {
        printf("  Warning: No common DMUs between periods\n")
        _dea_malm_stats(J(0,1,.),J(0,1,.),J(0,1,.))
        return
    }

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

        // Retain every computable component when another distance is unavailable.
        if (d_t_g>0 & d_t_g<. & d_t1_g>0 & d_t1_g<.) TFPCH[i]=d_t1_g/d_t_g
        if (d_tt>0 & d_tt<. & d_t1t1>0 & d_t1t1<.) EFFCH[i]=d_t1t1/d_tt
        if (TFPCH[i]<. & EFFCH[i]>0 & EFFCH[i]<.) TECHCH[i]=TFPCH[i]/EFFCH[i]
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

    mean_ec  = _dea_malm_mean(EFFCH)
    mean_tc  = _dea_malm_mean(TECHCH)
    mean_mpi = _dea_malm_mean(TFPCH)
    printf("%15s | %9.4f %9.4f %9.4f\n", "Mean", mean_ec, mean_tc, mean_mpi)

    _dea_malm_stats(EFFCH,TECHCH,TFPCH)
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

    real matrix X, Y
    real rowvector x, y, sx, sy
    real scalar j, rc

    // Include the evaluated row when balancing cross-period constraints.
    // Work on copies: reference arrays and the caller's observations are unchanged.
    sx=colmax(X_ref \ x0); sy=colmax(Y_ref \ y0)
    for (j=1; j<=cols(sx); j++) if (sx[j]==0) sx[j]=1
    for (j=1; j<=cols(sy); j++) if (sy[j]==0) sy[j]=1
    X=X_ref:/sx; Y=Y_ref:/sy; x=x0:/sx; y=y0:/sy

    n = rows(X_ref)
    m = cols(X_ref)
    s = cols(Y_ref)

    if (orientation == "input") {
        c      = (J(1, n, 0), 1)
        A_ineq = (X', -x') \ (-Y', J(s, 1, 0))
        b_ineq = J(m, 1, 0) \ (-y')
        lb     = J(1, n + 1, 0)
        ub     = J(1, n + 1, .)
        A_eq   = J(0, n + 1, .)
        b_eq   = J(0, 1, .)
        if (rts == "vrs") { A_eq = (J(1, n, 1), 0); b_eq = 1; }
    }
    else {
        c      = (J(1, n, 0), -1)
        A_ineq = (X', J(m, 1, 0)) \ (-Y', y')
        b_ineq = x' \ J(s, 1, 0)
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
    rc = _dea_nonte_optimize(S,tol)

    if (rc == 0 & S.value()<.) {
        eff = (orientation == "input") ? S.value() : 1 / (-S.value())
    }
    else {
        eff = .
    }

    if (eff<=0 | eff>=.) eff=.
    return(eff)
}


real scalar _dea_malm_mean(real colvector x)
{
    real colvector good
    good=select(x,x:<.)
    return(rows(good)>0 ? mean(good) : .)
}

void _dea_malm_stats(real colvector ec, real colvector tc, real colvector mpi)
{
    real scalar nf
    nf=sum((ec:>=.):|(tc:>=.):|(mpi:>=.))
    st_numscalar("r(n_dmus)",rows(ec))
    st_numscalar("r(mean_ec)",_dea_malm_mean(ec))
    st_numscalar("r(mean_tc)",_dea_malm_mean(tc))
    st_numscalar("r(mean_mpi)",_dea_malm_mean(mpi))
    st_numscalar("r(n_failed)",nf)
    st_numscalar("r(n_failed_ec)",sum(ec:>=.))
    st_numscalar("r(n_failed_tc)",sum(tc:>=.))
    st_numscalar("r(n_failed_mpi)",sum(mpi:>=.))
    if (nf>0) printf("  Note: %g matched DMU(s) have unavailable components; available components are retained.\n",nf)
}

end
