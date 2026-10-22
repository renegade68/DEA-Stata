*! version 2.0.28  05aug2026
*! DEA: Data Envelopment Analysis
*! Requires: Stata 16.0 or later
*! Program: dea — an update
*! Purpose: Data Envelopment Analysis (DEA) — Unified Command
*!          Supports DEA, IDEA, FDH with IPM and Simplex engines (full model parity)
*!
*! AUTHORS:
*! --------
*!   Choonjoo Lee (corresponding author)
*!     Korea National Defense University
*!     Seoul, Korea
*!     bloom.rampike@gmail.com
*!
*! CONTRIBUTORS TO PRIOR VERSIONS:
*! --------------------------------
*!   Yong-Bae Ji, Kyongrok Lee, Byung-Ihn Lim
*!     (contributors to dea.ado v1.x, dea_supereff.ado, dea_sbm.ado)
*!
*! ORIGINAL SOFTWARE:
*! ------------------
*!   Ji, Y.-B., & Lee, C. (2010). Data envelopment analysis.
*!     The Stata Journal, 10(2), 267-280.
*!     https://doi.org/10.1177/1536867X1001000207
*!
*!   Lee, C., & Ji, Y.-B. (2009). Data envelopment analysis in Stata.
*!     Presented at the 2009 Stata Conference, Washington, DC.
*!     https://ideas.repec.org/p/boc/dcon09/4.html
*!
*! PROVENANCE:
*! -----------
*! This program (dea.ado v2.0) is a major rewrite of the original dea.ado
*! v1.x package (Ji & Lee, 2010) with two new LP engine modules:
*!
*! (a) Inherited from dea.ado v1.x (Lee, Ji, Lee K., & Lim):
*!     - Syntax form 1: "input = output" parsing pattern
*!     - RTS/orientation normalization logic (CCR→CRS, BCC→VRS, etc.)
*!     - VRS frontier analysis structure (CRS→VRS→DRS sequential)
*!     - FDH and IDEA dispatch concepts
*!
*! (b) Substantially rewritten in v2.0:
*!     - _dea_simplex.ado: Two-Phase Revised Simplex with vectorized
*!       pricing, Bland's rule, Harris ratio test, eta-vector updates
*!       (replaces ldea.do which used per-iteration lusolve B^{-1})
*!     - _dea_ipm.ado: Interior Point Method engine via Stata 16+
*!       LinearProgram() class (entirely new)
*!     - SBM/Super-SBM via Charnes-Cooper transformation
*!       (dea_sbm.ado used direct objective modification)
*!     - Super-efficiency via DMU exclusion from constraint matrix
*!       (dea_supereff.ado used column-zeroing approach)
*!     - New models: cost, revenue, profit, DDF, multiplier,
*!       allocative, SBM variants, undesirable outputs, congestion,
*!       and bootstrap (Simar-Wilson)
*!     - Unified engine dispatch with auto-fallback
*!     - Modern st_data()/st_store() data interface
*!
*! REFERENCES — DEA Foundational Models:
*! --------------------------------------
*! [1]  Farrell, M. J. (1957). The measurement of productive efficiency.
*!        Journal of the Royal Statistical Society: Series A, 120(3), 253-290.
*!
*! [2]  Charnes, A., Cooper, W. W., & Rhodes, E. (1978). Measuring the
*!        efficiency of decision making units. European Journal of
*!        Operational Research, 2(6), 429-444.
*!        → CCR model: model(te), rts(crs)
*!
*! [3]  Banker, R. D., Charnes, A., & Cooper, W. W. (1984). Some models
*!        for estimating technical and scale inefficiencies in data
*!        envelopment analysis. Management Science, 30(9), 1078-1092.
*!        → BCC model: model(te), rts(vrs); scale efficiency
*!
*! [4]  Cooper, W. W., Seiford, L. M., & Tone, K. (2007). Data
*!        Envelopment Analysis: A Comprehensive Text with Models,
*!        Applications, References and DEA-Solver Software (2nd ed.).
*!        Springer.
*!        → Envelopment/multiplier LP formulations; cost, revenue,
*!          profit efficiency (Ch.8, eq.8.8-8.9)
*!
*! REFERENCES — Extended Models:
*! -----------------------------
*! [5]  Charnes, A., Cooper, W. W., Golany, B., Seiford, L., & Stutz, J.
*!        (1985). Foundations of data envelopment analysis for Pareto-
*!        Koopmans efficient empirical production functions. Journal of
*!        Econometrics, 30(1-2), 91-107.
*!        → Additive model: model(additive)
*!
*! [6]  Andersen, P., & Petersen, N. C. (1993). A procedure for ranking
*!        efficient units in data envelopment analysis. Management
*!        Science, 39(10), 1261-1264.
*!        → Radial super-efficiency: model(supereff)
*!
*! [7]  Chambers, R. G., Chung, Y., & Färe, R. (1996). Benefit and
*!        distance functions. Journal of Economic Theory, 70(2), 407-419.
*!        → Directional distance function: model(ddf)
*!
*! [8]  Chambers, R. G., Chung, Y., & Färe, R. (1998). Profit,
*!        directional distance functions, and Nerlovian efficiency.
*!        Journal of Optimization Theory and Applications, 98(2), 351-364.
*!        → DDF and profit efficiency relationship
*!
*! [8b] Ali, A.I., Lerme, C.S., & Seiford, L.M. (1995). Components of
*!        efficiency evaluation in data envelopment analysis. European
*!        Journal of Operational Research, 80(3), 462-473.
*!        → Framework: envelopment surface, orientation, pricing mechanism
*!           (conceptual basis for non-oriented/oriented SBM distinction)
*!
*! [9]  Tone, K. (2001). A slacks-based measure of efficiency in data
*!        envelopment analysis. European Journal of Operational Research,
*!        130(3), 498-509.
*!        → SBM: model(sbm) — oriented and non-oriented
*!
*! [10] Tone, K. (2002). A slacks-based measure of super-efficiency in
*!        data envelopment analysis. European Journal of Operational
*!        Research, 143(1), 32-41.
*!        → Super-SBM: model(supersbm)
*!
*! [10b] Tone, K. (2008). Variations on the theme of slacks-based measure
*!        of efficiency in DEA. Discussion Paper 08-14, National Graduate
*!        Institute for Policy Studies (GRIPS), Tokyo.
*!        [Published: EJOR 200(3), 901-907, 2010]
*!        → SBM variations: Variation I–IV (facet-based)
*!
*! [10c] Tone, K. (2004). Dealing with undesirable outputs in DEA: A
*!        Slacks-Based Measure (SBM) approach. GRIPS Discussion Paper
*!        I-2003-0005. doi:10.24545/00000955
*!        → SBM with undesirable outputs: model(sbmundesirable)
*!
*! [11] Deprins, D., Simar, L., & Tulkens, H. (2006). Measuring labor-
*!        efficiency in post offices. In P. Chander, J. Drèze, C. K.
*!        Lovell, & J. Mintz (Eds.), Public Goods, Environmental
*!        Externalities and Fiscal Competition. Springer, Boston, MA.
*!        doi:10.1007/978-0-387-25534-7_16
*!        → Free Disposal Hull: model(fdh)
*!
*! [12] Cooper, W. W., Park, K. S., & Yu, G. (1999). IDEA and AR-IDEA:
*!        Models for dealing with imprecise data in DEA. Management
*!        Science, 45(4), 597-607.
*!        → Imprecise DEA: model(imprecise)
*!
*! [12b] Cooper, W. W., Deng, H., Seiford, L. M., & Zhu, J. (2004).
*!        Congestion: Its identification and management with DEA. In
*!        W. W. Cooper, L. M. Seiford, & J. Zhu (Eds.), Handbook on
*!        Data Envelopment Analysis (pp. 177-201). Kluwer Academic
*!        Publishers. doi:10.1007/1-4020-7798-X_7
*!        → Input congestion: model(congestion)
*!
*! REFERENCES — Decomposition & Scale:
*! ------------------------------------
*! [13] Färe, R., Grosskopf, S., & Lovell, C. A. K. (1993). Production
*!        Frontiers. Cambridge University Press.
*!        → RTS determination via CRS/VRS/DRS comparison
*!
*! [14] Färe, R., & Lovell, C. A. K. (1978). Measuring the technical
*!        efficiency of production. Journal of Economic Theory, 19(1),
*!        150-162.
*!        → Allocative efficiency decomposition: model(allocative)
*!
*! [15] Charnes, A., & Cooper, W. W. (1962). Programming with linear
*!        fractional functionals. Naval Research Logistics Quarterly,
*!        9(3-4), 181-186.
*!        → Charnes-Cooper transformation (SBM, Super-SBM linearization)
*!
*! REFERENCES — Statistical Inference:
*! ------------------------------------
*! [16] Simar, L., & Wilson, P. W. (1998). Sensitivity analysis of
*!        efficiency scores: How to bootstrap in nonparametric frontier
*!        models. Management Science, 44(1), 49-61.
*!        → Bootstrap Algorithm 1 (bias correction)
*!
*! [17] Simar, L., & Wilson, P. W. (2007). Estimation and inference in
*!        two-stage, semi-parametric models of production processes.
*!        Journal of Econometrics, 136(1), 31-64.
*!        → Bootstrap Algorithm 2 (two-stage with environmental variables)
*!
*! REFERENCES — LP Algorithms (see also _dea_simplex.ado, _dea_ipm.ado):
*! ----------------------------------------------------------------------
*! [18] Dantzig, G. B. (1963). Linear Programming and Extensions.
*!        Princeton University Press.
*!
*! [19] Bland, R. G. (1977). New finite pivoting rules for the simplex
*!        method. Mathematics of Operations Research, 2(2), 103-107.
*!
*! [20] Harris, P. M. J. (1973). Pivot selection methods of the Devex LP
*!        code. Mathematical Programming, 5(1), 1-28.
*!        → Harris ratio test
*!
*! [21] Chen, D.-S., Batson, R. G., & Dang, Y. (2010). Applied Integer
*!        Programming. Wiley.
*!        → Revised simplex: B^{-1} maintenance (Section 9.8)
*!
*! [22] Bazaraa, M. S., Jarvis, J. J., & Sherali, H. D. (2010). Linear
*!        Programming and Network Flows (4th ed.). Wiley.
*!        → Eta-vector (product form) updates, degeneracy theory

/**
 * HISTORY:
 * ---------------------------------------------------------------------------
 * 2009-07-30: DEA Primer at 2009 DC Stata Conference
 * 2013-01-04: Public release at SourceForge
 * 2024-12-28: dea.ado v1.2.2 (Simplex engine, enhanced comments)
 * 2026-01-26: dea_lp.ado v6.1.0 (IPM engine, dual fallback)
 * 2026-02-05: dea.ado v2.0.0 — Unified command integrating both engines
 * 2026-02-08: dea.ado v2.0.1 — supereff: reference set (lambda) support
 * 2026-02-10: dea.ado v2.0.2 — fix: lambda tempvar scoping for reference sets
 * 2026-02-13: dea.ado v2.0.3 — engine(auto): simplex-first with IPM fallback
 * 2026-02-14: dea.ado v2.0.4 — add model(profit): profit efficiency
 * 2026-05-02: dea.ado v2.0.7 — allow bootstrap after simplex-generated scores
 * 2026-05-03: dea.ado v2.0.8 — fix SBM-undesirable score recovery after Charnes-Cooper linearization
 * 2026-07-19: dea.ado v2.0.9 — preserve VRS returns and validate optional matrices before copying
 * 2026-07-19: dea.ado v2.0.10 — restore documented top-level parser returns
 * 2026-07-19: dea.ado v2.0.11 — reject unsupported explicit orientation for SBM
 * 2026-07-19: dea.ado v2.0.12 — correct additive ranking and display labels
 * 2026-07-19: dea.ado v2.0.13 — add simplex/IPM parity for oriented SBM
 * 2026-07-19: dea.ado v2.0.14 — validate economic-model prices and ratio profit domain
 * 2026-07-20: dea.ado v2.0.16 — apply parser-returned if/in conditions to the estimation sample
 * 2026-07-20: dea.ado v2.0.17 — distinguish Form B equality tests in if qualifiers from Form A '=' syntax
 * 2026-07-20: dea.ado v2.0.18 — fix IPM DRS/IRS stage-2 scalar min/max syntax
 * 2026-07-20: dea.ado v2.0.19 — add developer-only canonical shadow comparison
 * 2026-07-20: dea.ado v2.0.20 — return shadow fallback counts for promotion reports
 * 2026-07-20: dea.ado v2.0.21 — use row scaling in stage-1 shadow comparisons
 * 2026-07-20: dea.ado v2.0.22 — preselect IPM for heterogeneous-scale radial auto calls
 * 2026-07-20: dea.ado v2.0.23 — expose the signed VRS multiplier intercept
 * 2026-07-20: dea.ado v2.0.24 — correct desirable-output multiplier rows
 * 2026-07-20: dea.ado v2.0.25 — align input u0 sign and DRS/IRS bounds with Cooper
 * 2026-07-21: dea.ado v2.0.26 — disclose the effective specification and defaults
 * 2026-07-21: dea.ado v2.0.27 — package-level stored-result contract repair
 * 2026-08-05: dea.ado v2.0.28 — stabilize radial TE scaling and restore numeric dmu() VRS classification
 * ---------------------------------------------------------------------------
 *
 * ARCHITECTURE:
 * ---------------------------------------------------------------------------
 * dea.ado v2.0 (this file)
 *   Entry point: syntax parsing, version detection, engine dispatch, results
 *
 * _dea_ipm.ado (IPM engine, Stata 16+)
 *   LinearProgram() class, Primal→Dual fallback, all RTS, Stage 1+2
 *   Mata routines used by radial, non-radial, imprecise, FDH,
 *   undesirable-output, and congestion-related dispatchers
 *
 * _dea_simplex.ado (Simplex engine, Mata)
 *   Two-Phase Revised Simplex, Bland's rule, Harris ratio test
 *   Package-supplied LP routines for the standard DEA model family
 *   General-purpose LP converter for non-radial models
 *
 * dea.sthlp (help file)
 * ---------------------------------------------------------------------------
 *
 * MODULE DISPATCH MAP:
 * ---------------------------------------------------------------------------
 *   dea (entry point)
 *   ├─ _dea_parse_equal      → syntax form: inputs = outputs
 *   ├─ _dea_parse_standard   → syntax form: outputs, inputs()
 *   ├─ _dea_normalize_options→ option normalization and engine selection
 *   ├─ _dea_validate_model_opts
 *   ├─ _dea_model_dispatch   → standard model dispatcher → IPM or Simplex
 *   │   ├─ model(te/cost/revenue/profit/ddf/sbm/additive)
 *   │   ├─ model(supereff/multiplier/allocative/supersbm)
 *   │   └─ model(undesirable)
 *   ├─ _dea_fdh_dispatch     → FDH dominance comparison
 *   ├─ _dea_imprecise_dispatch
 *   ├─ _dea_congestion_dispatch
 *   ├─ _dea_sbmundesirable_dispatch
 *   ├─ _dea_vrsfrontier      → VRS frontier analysis (SCALE, RTS)
 *   └─ _dea_display          → console output + reference sets
 * ---------------------------------------------------------------------------
 */


********************************************************************************
* 1. MAIN ENTRY POINT
********************************************************************************

program define dea, rclass
    version 16.0

    // =========================================================================
    // 1.0 Load IPM engine Mata functions (if Stata 16+ and not already loaded)
    // =========================================================================
    if `c(stata_version)' >= 16 {
        // Check core function AND undesirable model function (added in v2.0)
        local _ipm_needs_load 0
        capture mata: mata which _dea_ipm_envelopment()
        if _rc != 0 local _ipm_needs_load 1
        capture mata: mata which _dea_ipm_ddf_undesirable()
        if _rc != 0 local _ipm_needs_load 1
        if `_ipm_needs_load' {
            capture findfile _dea_ipm.ado
            if _rc == 0 {
                qui run "`r(fn)'"
            }
            else {
                di as error "_dea_ipm.ado not found. Please install it in your ado path."
                exit 601
            }
        }
    }

    // =========================================================================
    // 1.1 Preserve original command line (syntax consumes `0')
    // =========================================================================
    local cmdline `"dea `0'"'
    local raw0 `"`0'"'
    local _no_option_call = (strpos(`"`raw0'"', ",") == 0)

    // =========================================================================
    // 1.2 Detect syntax form and parse
    // =========================================================================
    //   Form A (legacy):   dea x1 x2 = y1 y2, rts() ort() ...
    //   Form B (extended): dea y1 y2, inputs(x1 x2) model() ...
    //
    // Form A's separator must precede any if/in qualifier.  Equality tests
    // such as "if group==1" belong to Form B and are not syntax separators.
    local _eqpos = strpos(`"`raw0'"', "=")
    local _ifpos = strpos(lower(`" `raw0'"'), " if ")
    local _inpos = strpos(lower(`" `raw0'"'), " in ")
    local _comma = strpos(`"`raw0'"', ",")
    local _boundary = length(`"`raw0'"') + 1
    foreach p in _ifpos _inpos _comma {
        if ``p'' > 0 & ``p'' < `_boundary' local _boundary = ``p''
    }
    local has_equal = (`_eqpos' > 0 & `_eqpos' < `_boundary')

    if `has_equal' {
        _dea_parse_equal `raw0'
    }
    else {
        _dea_parse_standard `raw0'
    }

    // At this point, the following locals are set by the parser:
    //   inputs, outputs, dmu, rts, orientation, model, stage,
    //   method, tolerance, generate, trace, saving, replace, debug,
    //   fdh, reference, saveref, saveweights,
    //   w, p, gx, gy, ddfbeta,
    //   xlo, xhi, ylo, yhi, xord, yord,
    //   bounded, ordinal,
    //   bootstrap, envars, reps1, reps2, level,
    //   ifcond, incond, m_in, m_out,
    //   engine, engine_label
    // These are returned via r() from the parser and captured here.

    // Capture all parsed locals from r()
    local inputs    "`r(inputs)'"
    local outputs   "`r(outputs)'"
    local dmu       "`r(dmu)'"
    local rts       "`r(rts)'"
    local orientation "`r(orientation)'"
    local model     "`r(model)'"
    local stage     "`r(stage)'"
    local method    "`r(method)'"
    local tolerance "`r(tolerance)'"
    local generate  "`r(generate)'"
    local trace     "`r(trace)'"
    local saving    "`r(saving)'"
    local replace   "`r(replace)'"
    local debug     "`r(debug)'"
    local contractcheck "`r(contractcheck)'"
    local fdh       "`r(fdh)'"
    local reference "`r(reference)'"
    local saveref   "`r(saveref)'"
    local saveweights "`r(saveweights)'"
    local w         "`r(w)'"
    local p         "`r(p)'"
    local gx        "`r(gx)'"
    local gy        "`r(gy)'"
    local ddfbeta   "`r(ddfbeta)'"
    local bad       "`r(bad)'"
    local direction "`r(direction)'"
    local xlo       "`r(xlo)'"
    local xhi       "`r(xhi)'"
    local ylo       "`r(ylo)'"
    local yhi       "`r(yhi)'"
    local xord      "`r(xord)'"
    local yord      "`r(yord)'"
    local bounded   "`r(bounded)'"
    local ordinal   "`r(ordinal)'"
    local bootstrap "`r(bootstrap)'"
    local envars    "`r(envars)'"
    local reps1     "`r(reps1)'"
    local reps2     "`r(reps2)'"
    local level     "`r(level)'"
    local m_in      "`r(m_in)'"
    local m_out     "`r(m_out)'"
    local engine    "`r(engine)'"
    local engine_label "`r(engine_label)'"
    local auto_fallback "`r(auto_fallback)'"
    local ifcond    "`r(ifcond)'"
    local incond    "`r(incond)'"
    local user_set_ort = r(user_set_ort)
    if "`user_set_ort'" == "" | "`user_set_ort'" == "." local user_set_ort 0
    // sbm_ort: 0=non-oriented, 1=input-oriented, 2=output-oriented
    // Integer avoids Stata keyword conflicts (e.g. "input" is a Stata command)
    local sbm_ort_num = 0
    if `user_set_ort' {
        if "`orientation'" == "input"  local sbm_ort_num 1
        if "`orientation'" == "output" local sbm_ort_num 2
    }

    // =========================================================================
    // 1.3 Build estimation sample (must be done in main program)
    // =========================================================================
    tempvar touse
    mark `touse' `ifcond' `incond'
    local _markout_vars `inputs' `outputs'
    if "`w'" != "" local _markout_vars `_markout_vars' `w'
    if "`p'" != "" local _markout_vars `_markout_vars' `p'
    if "`bounded'" != "" {
        local _markout_vars : list _markout_vars - bounded
    }
    markout `touse' `_markout_vars'

    // Economic models require strictly positive prices.  Ratio profit
    // efficiency additionally requires strictly positive observed profit;
    // negative-profit/Nerlovian measures have a separate contract.
    if inlist("`model'", "cost", "revenue", "profit", "allocative") {
        foreach _price of varlist `w' `p' {
            quietly count if `_price' <= 0 & `touse'
            if r(N) > 0 {
                di as error "economic-efficiency price variables must be strictly positive"
                exit 459
            }
        }
    }
    if "`model'" == "profit" {
        tempvar _observed_profit
        quietly generate double `_observed_profit' = 0 if `touse'
        forvalues _j = 1/`m_out' {
            local _yj : word `_j' of `outputs'
            local _pj : word `_j' of `p'
            quietly replace `_observed_profit' = `_observed_profit' + `_pj'*`_yj' if `touse'
        }
        forvalues _i = 1/`m_in' {
            local _xi : word `_i' of `inputs'
            local _wi : word `_i' of `w'
            quietly replace `_observed_profit' = `_observed_profit' - `_wi'*`_xi' if `touse'
        }
        quietly count if `_observed_profit' <= 0 & `touse'
        if r(N) > 0 {
            di as error "model(profit) ratio requires strictly positive observed profit for every evaluated DMU"
            di as error "use a separately specified directional/Nerlovian measure for zero or negative profit"
            exit 459
        }
    }
    
    // Handle DMU variable
    local _generated_dmu 0
    if "`dmu'" == "" {
        cap drop _dmu
        qui gen str _dmu = "dmu" + string(_n) if `touse'
        local dmu "_dmu"
        local _generated_dmu 1
    }
    
    // Count DMUs
    qui count if `touse'
    if r(N) == 0 {
        di as error "No observations found"
        exit 2000
    }
    local n_dmus = r(N)

    // Under engine(auto), avoid sending a radial DEA problem with strongly
    // heterogeneous dimensional scales through the unscaled public simplex
    // path.  Explicit engine(simplex) remains an unconditional user choice.
    local _auto_scale_ratio 1
    local _auto_preemptive_ipm 0
    local _auto_scale_threshold 1e8
    if "`auto_fallback'" != "" & "`model'" == "te" & "`fdh'" == "" {
        local _auto_scale_min .
        local _auto_scale_max 0
        foreach _scale_var of varlist `inputs' `outputs' {
            quietly summarize `_scale_var' if `touse', meanonly
            local _var_scale = max(abs(r(min)), abs(r(max)))
            if `_var_scale' > 0 {
                if missing(`_auto_scale_min') | `_var_scale' < `_auto_scale_min' {
                    local _auto_scale_min = `_var_scale'
                }
                if `_var_scale' > `_auto_scale_max' {
                    local _auto_scale_max = `_var_scale'
                }
            }
        }
        if !missing(`_auto_scale_min') & `_auto_scale_min' > 0 {
            local _auto_scale_ratio = `_auto_scale_max' / `_auto_scale_min'
        }
        if `_auto_scale_ratio' > `_auto_scale_threshold' {
            local engine "ipm"
            local engine_label "Interior Point Method (auto-scale safeguard)"
            local _auto_preemptive_ipm 1
        }
    }

    // =========================================================================
    // 1.4 Display header
    // =========================================================================
    di as text "{hline `c(linesize)'}"
    di as text "dea v2.0.28 — Data Envelopment Analysis"
    di as text "Engine: `engine_label'"
    di as text "Specification: model(`model') rts(`rts') orientation(`orientation') stage(`stage')"
    if `_no_option_call' {
        di as text "Defaults applied: model(te) rts(crs) orientation(input) stage(2)" ///
            " engine(simplex) tolerance(1e-7)"
    }
    if `_generated_dmu' {
        di as text "DMU identifiers: generated from observation order (_dmu)"
    }
    if "`debug'" != "" {
        di as text "DMUs: `n_dmus', Inputs: `m_in', Outputs: `m_out'"
        di as text "Inputs: `inputs'"
        di as text "Outputs: `outputs'"
        di as text "Model: `rts' `orientation' (`model')"
        di as text "Stage: `stage'"
    }

    // =========================================================================
    // 1.5 Dispatch to engine
    // =========================================================================

    // FDH: Branch & Bound — no LP needed, engine-independent
    if "`fdh'" != "" | "`model'" == "fdh" {
        _dea_fdh_dispatch, inputs(`inputs') outputs(`outputs') ///
            dmu(`dmu') orientation(`orientation') stage(`stage') ///
            touse(`touse') generate(`generate') tolerance(`tolerance') ///
            `trace' `debug' `reference' `saveref' ///
            saving(`saving') `replace'
    }
    // IDEA (Imprecise DEA) — with bounded/ordinal/interval data
    else if "`bounded'" != "" | "`ordinal'" != "" | ///
            "`xlo'" != "" | "`xhi'" != "" | "`ylo'" != "" | "`yhi'" != "" | ///
            "`model'" == "imprecise" {
        _dea_imprecise_dispatch, inputs(`inputs') outputs(`outputs') ///
            dmu(`dmu') rts(`rts') orientation(`orientation') ///
            stage(`stage') engine(`engine') touse(`touse') ///
            generate(`generate') tolerance(`tolerance') ///
            bounded(`bounded') ordinal(`ordinal') ///
            xlo(`xlo') xhi(`xhi') ylo(`ylo') yhi(`yhi') ///
            xord(`xord') yord(`yord') ///
            `trace' `debug' saving(`saving') `replace'
    }
    // Congestion model (Cooper, Deng, Seiford & Zhu, 2004)
    else if "`model'" == "congestion" {
        _dea_congestion_dispatch, inputs(`inputs') outputs(`outputs') ///
            dmu(`dmu') engine(`engine') touse(`touse') ///
            generate(`generate') tolerance(`tolerance') `debug'
    }
    // SBM-Undesirable (Tone, 2004) — SBM with bad outputs
    else if "`model'" == "sbmundesirable" {
        _dea_sbmundesirable_dispatch, inputs(`inputs') outputs(`outputs') ///
            bad(`bad') dmu(`dmu') rts(`rts') touse(`touse') ///
            generate(`generate') tolerance(`tolerance') `debug'
    }
    // Standard DEA model(te) and extended models
    else {
        _dea_model_dispatch, inputs(`inputs') outputs(`outputs') ///
            dmu(`dmu') rts(`rts') orientation(`orientation') ///
            model(`model') stage(`stage') engine(`engine') ///
            touse(`touse') generate(`generate') tolerance(`tolerance') ///
            method(`method') ///
            w(`w') p(`p') gx(`gx') gy(`gy') ddfbeta(`ddfbeta') ///
            bad(`bad') direction(`direction') ///
            `reference' `saveref' `saveweights' ///
            sbmort(`sbm_ort_num') ///
            `trace' `debug' saving(`saving') `replace'
    }

    // =========================================================================
    // 1.5a Auto fallback: simplex failed → retry all with IPM
    // =========================================================================
    // Save r() from _dea_model_dispatch before subsequent calls clear it
    local _r_N_eff    = r(n_efficient)
    local _r_mean_eff = r(mean_efficiency)
    local _r_min_eff  = r(min_eff)
    local _r_max_eff  = r(max_eff)
    local _r_n_failed = r(n_failed)
    tempname _r_lambda_mat _r_slack_in_mat _r_slack_out_mat _r_vrsfrontier_mat _r_multiplier_u0_mat
    local _r_has_lambda   0
    local _r_has_slack_in 0
    local _r_has_slack_out 0
    local _r_has_vrsfrontier 0
    local _r_has_multiplier_u0 0
    // Test existence before copying.  Direct assignment from a nonexistent
    // r(name) creates a 1 x 1 missing matrix and falsely looks successful.
    capture quietly matrix list r(lambda)
    if _rc == 0 {
        matrix `_r_lambda_mat' = r(lambda)
        local _r_has_lambda 1
    }
    capture quietly matrix list r(slack_in)
    if _rc == 0 {
        matrix `_r_slack_in_mat' = r(slack_in)
        local _r_has_slack_in 1
    }
    capture quietly matrix list r(slack_out)
    if _rc == 0 {
        matrix `_r_slack_out_mat' = r(slack_out)
        local _r_has_slack_out 1
    }
    capture quietly matrix list r(multiplier_u0)
    if _rc == 0 {
        matrix `_r_multiplier_u0_mat' = r(multiplier_u0)
        local _r_has_multiplier_u0 1
    }

    if "`auto_fallback'" != "" & "`fdh'" == "" & ///
       !inlist("`model'", "congestion", "sbmundesirable") {
        qui count if `generate' >= . & `touse'
        local n_failed = r(N)
        if `n_failed' > 0 {
            di as text ""
            di as text "Note: `n_failed' DMU(s) returned missing with simplex engine."
            di as text "      Retrying all `n_dmus' DMUs with IPM engine..."
            di as text ""
            local engine "ipm"
            local engine_label "Interior Point Method (auto-fallback)"
            // Reset efficiency variable
            qui replace `generate' = . if `touse'
            if "`bounded'" != "" | "`ordinal'" != "" | ///
                    "`xlo'" != "" | "`xhi'" != "" | "`ylo'" != "" | "`yhi'" != "" | ///
                    "`model'" == "imprecise" {
                _dea_imprecise_dispatch, inputs(`inputs') outputs(`outputs') ///
                    dmu(`dmu') rts(`rts') orientation(`orientation') ///
                    stage(`stage') engine(`engine') touse(`touse') ///
                    generate(`generate') tolerance(`tolerance') ///
                    bounded(`bounded') ordinal(`ordinal') ///
                    xlo(`xlo') xhi(`xhi') ylo(`ylo') yhi(`yhi') ///
                    xord(`xord') yord(`yord') ///
                    `trace' `debug' saving(`saving') `replace'
            }
            else {
                _dea_model_dispatch, inputs(`inputs') outputs(`outputs') ///
                    dmu(`dmu') rts(`rts') orientation(`orientation') ///
                    model(`model') stage(`stage') engine(`engine') ///
                    touse(`touse') generate(`generate') tolerance(`tolerance') ///
                    method(`method') ///
                    w(`w') p(`p') gx(`gx') gy(`gy') ddfbeta(`ddfbeta') ///
                    bad(`bad') direction(`direction') ///
                    `reference' `saveref' `saveweights' ///
                    sbmort(`sbm_ort_num') ///
                    `trace' `debug' saving(`saving') `replace'
            }
        }
    }

    // =========================================================================
    // 1.5 VRS Frontier analysis (if rts == VRS and model == te)
    // =========================================================================
    if "`rts'" == "vrs" & "`model'" == "te" & "`fdh'" == "" {
        // Use the engine that produced the main scores.  With engine(auto),
        // this remains simplex unless the auto-fallback above switched to IPM.
        local _vrs_eng "`engine'"
        _dea_vrsfrontier, inputs(`inputs') outputs(`outputs') ///
            dmu(`dmu') orientation(`orientation') ///
            engine(`_vrs_eng') touse(`touse') tolerance(`tolerance') ///
            `trace' `debug'
        // Preserve the internal r-class results before mkmat or bootstrap
        // replaces them in the top-level command.
        matrix `_r_vrsfrontier_mat' = r(vrsfrontier)
        local _r_vrs_n_crs = r(n_crs)
        local _r_vrs_n_drs = r(n_drs)
        local _r_vrs_n_irs = r(n_irs)
        local _r_has_vrsfrontier 1
    }

    // =========================================================================
    // 1.6 Simar-Wilson Bootstrap (if requested)
    // =========================================================================
    if "`bootstrap'" != "" {
        _dea_bootstrap, inputs(`inputs') outputs(`outputs') ///
            dmu(`dmu') rts(`rts') orientation(`orientation') ///
            generate(`generate') touse(`touse') engine(`engine') ///
            envars(`envars') reps1(`reps1') reps2(`reps2') level(`level') ///
            `debug'
    }

    // =========================================================================
    // 1.6a Developer-only canonical shadow comparison
    // =========================================================================
    local _contract_checked 0
    local _contract_match .
    local _contract_max_gap .
    local _contract_n_failed .
    local _contract_stage2_failed .
    local _contract_fallbacks .
    local _contract_stage2_fallbacks .
    local _contract_scaling ""
    tempname _r_contract_eff_mat
    if "`contractcheck'" != "" {
        tempvar _contract_eff _contract_status _contract_gap
        local _contract_tol = max(`tolerance', 1e-6)
        local _contract_scaling = cond(`stage'==1, "row", "none")
        quietly _dea_tecontract if `touse', inputs(`inputs') outputs(`outputs') ///
            rts(`rts') orientation(`orientation') stage(`stage') ///
            engine(auto) scaling(`_contract_scaling') tolerance(`_contract_tol') ///
            generate(`_contract_eff') status(`_contract_status')
        local _contract_n_failed = r(n_failed)
        local _contract_stage2_failed = r(stage2_failures)
        local _contract_fallbacks = r(fallbacks)
        local _contract_stage2_fallbacks = r(stage2_fallbacks)
        matrix `_r_contract_eff_mat' = r(efficiency)
        quietly generate double `_contract_gap' = ///
            abs(`generate'-`_contract_eff') if `touse'
        quietly summarize `_contract_gap' if `touse', meanonly
        local _contract_max_gap = r(max)
        local _contract_checked 1
        local _contract_match = (`_contract_n_failed'==0 & ///
            `_contract_max_gap'<=`_contract_tol')
        di as text "Canonical shadow check: " ///
            as result cond(`_contract_match', "matched", "review required") ///
            as text "; max efficiency gap = " ///
            as result %10.3e `_contract_max_gap'
    }

    // =========================================================================
    // 1.7 Store return values
    // =========================================================================
    // NOTE: return matrix MUST come BEFORE return add in Stata 16.
    // Build efficiency matrix from generate variable
    tempname _eff_mat
    qui mkmat `generate' if `touse', matrix(`_eff_mat')
    return matrix efficiency = `_eff_mat'

    if `_r_has_lambda' {
        return matrix lambda = `_r_lambda_mat'
    }
    if `_r_has_slack_in' {
        return matrix slack_in = `_r_slack_in_mat'
    }
    if `_r_has_slack_out' {
        return matrix slack_out = `_r_slack_out_mat'
    }
    if `_r_has_vrsfrontier' {
        return matrix vrsfrontier = `_r_vrsfrontier_mat'
    }
    if `_r_has_multiplier_u0' {
        return matrix multiplier_u0 = `_r_multiplier_u0_mat'
    }
    if `_contract_checked' {
        return matrix contract_efficiency = `_r_contract_eff_mat'
    }

    return add
    return local cmd "dea"
    return local cmdline `"`cmdline'"'
    return local engine "`engine'"
    return local model "`model'"
    return local rts "`rts'"
    return local orientation "`orientation'"
    return local inputs "`inputs'"
    return local outputs "`outputs'"
    return local generate "`generate'"
    return scalar defaults_applied = `_no_option_call'
    return scalar generated_dmu = `_generated_dmu'
    return scalar N          = `n_dmus'
    return scalar N_eff      = `_r_N_eff'
    return scalar mean_eff   = `_r_mean_eff'
    return scalar min_eff    = `_r_min_eff'
    return scalar max_eff    = `_r_max_eff'
    return scalar n_failed   = `_r_n_failed'
    return scalar stage      = `stage'
    return scalar tolerance  = `tolerance'
    return scalar auto_scale_ratio = `_auto_scale_ratio'
    return scalar auto_scale_threshold = `_auto_scale_threshold'
    return scalar auto_preemptive_ipm = `_auto_preemptive_ipm'
    return scalar contract_checked = `_contract_checked'
    return scalar contract_match = `_contract_match'
    return scalar contract_max_gap = `_contract_max_gap'
    return scalar contract_n_failed = `_contract_n_failed'
    return scalar contract_stage2_failures = `_contract_stage2_failed'
    return scalar contract_fallbacks = `_contract_fallbacks'
    return scalar contract_stage2_fallbacks = `_contract_stage2_fallbacks'
    return local contract_scaling "`_contract_scaling'"
    if `_r_has_vrsfrontier' {
        return scalar n_crs = `_r_vrs_n_crs'
        return scalar n_drs = `_r_vrs_n_drs'
        return scalar n_irs = `_r_vrs_n_irs'
    }

    di as text "{hline `c(linesize)'}"
end


********************************************************************************
* 2. SYNTAX PARSERS
********************************************************************************

// ============================================================================
// 2.1 Form A Parser: dea x1 x2 = y1 y2 [if] [in], [options]
// ============================================================================
program define _dea_parse_equal, rclass
    version 16.0

    // ---------------------------------------------------------
    // Split at "=" to get inputs (LHS) and outputs+options (RHS)
    // ---------------------------------------------------------
    local raw0 `"`0'"'

    // Extract input variables before "="
    local invarnames ""
    gettoken word 0 : 0, parse("=,")
    while "`word'" != "=" & "`word'" != "" {
        if "`word'" == "," {
            di as error "Syntax error: comma before '=' sign"
            error 198
        }
        local invarnames `invarnames' `word'
        gettoken word 0 : 0, parse("=,")
    }
    if "`invarnames'" == "" {
        di as error "No input variables specified before '='"
        error 198
    }
    unab invarnames : `invarnames'

    // Now `0' contains: outputs [if] [in], [options]
    #del ;
    syntax varlist(min=1 numeric) [if] [in] [using/]
    [,
        DMU(varname)
        RTS(string)
        ORT(string)
        ORIentation(string)
        STAGE(integer 2)
        FDH
        MODEL(string)
        METHod(string)
        GENerate(name)
        TOLerance(real 1e-7)
        TOL1(real 1e-14)
        TOL2(real 1e-8)
        TRACE
        SAVing(string)
        REPLACE
        DEBUG
        ENGine(string)
        CONTRACTCHECK

        REFerence
        SAVERef
        SAVEWeights
        W(varlist numeric)
        P(varlist numeric)
        GX(varlist numeric)
        GY(varlist numeric)
        DDFBETA(name)
        BAD(varlist numeric)
        DIRection(string)

        Bounded(varlist)
        ORDinal(varlist)
        XLO(varlist numeric)
        XHI(varlist numeric)
        YLO(varlist numeric)
        YHI(varlist numeric)
        XORD(varlist numeric)
        YORD(varlist numeric)

        BOOTstrap
        ENVars(varlist numeric)
        REPS1(integer 100)
        REPS2(integer 2000)
        Level(integer 95)
    ];
    #del cr

    local outputs `varlist'

    // Reconcile orientation: ort() (legacy) vs orientation() (new)
    if "`ort'" != "" & "`orientation'" != "" {
        di as error "Specify only one of: ort() or orientation()"
        exit 198
    }
    if "`ort'" != "" local orientation "`ort'"

    // Reconcile tolerance: tol1/tol2 (legacy) vs tolerance() (new)
    // If user explicitly set tolerance (not default), use it
    // Otherwise fall back to tol1 as primary tolerance
    if `tolerance' == 1e-7 & `tol1' != 1e-14 {
        local tolerance = `tol1'
    }

    // Delegate to common normalizer
    _dea_normalize_options,                           ///
        inputs(`invarnames') outputs(`outputs')       ///
        dmu(`dmu') rts(`rts') orientation(`orientation') ///
        model(`model') stage(`stage') method(`method')   ///
        fdh(`fdh') tolerance(`tolerance')             ///
        generate(`generate') trace(`trace')           ///
        saving(`saving') replace(`replace') debug(`debug') ///
        engine(`engine') contractcheck(`contractcheck') ///
        reference(`reference') saveref(`saveref')     ///
        saveweights(`saveweights')                    ///
        w(`w') p(`p') gx(`gx') gy(`gy') ddfbeta(`ddfbeta') ///
        bad(`bad') direction(`direction')             ///
        bounded(`bounded') ordinal(`ordinal')         ///
        xlo(`xlo') xhi(`xhi') ylo(`ylo') yhi(`yhi')  ///
        xord(`xord') yord(`yord')                     ///
        bootstrap(`bootstrap') envars(`envars')       ///
        reps1(`reps1') reps2(`reps2') level(`level')  ///
        ifcond(`if') incond(`in') using(`using')

    return add
end

// ============================================================================
// 2.2 Form B Parser: dea y1 y2, inputs(x1 x2) [options]
// ============================================================================
program define _dea_parse_standard, rclass
    version 16.0

    #del ;
    syntax varlist(min=1 numeric) [if] [in]
    [,
        Inputs(varlist numeric min=1)
        DMU(varname)
        RTS(string)
        ORT(string)
        ORIentation(string)
        STAGE(integer 2)
        FDH
        MODEL(string)
        METHod(string)
        GENerate(name)
        TOLerance(real 1e-7)
        TRACE
        SAVing(string)
        REPLACE
        DEBUG
        ENGine(string)
        CONTRACTCHECK

        REFerence
        SAVERef
        SAVEWeights
        W(varlist numeric)
        P(varlist numeric)
        GX(varlist numeric)
        GY(varlist numeric)
        DDFBETA(name)
        BAD(varlist numeric)
        DIRection(string)

        Bounded(varlist)
        ORDinal(varlist)
        XLO(varlist numeric)
        XHI(varlist numeric)
        YLO(varlist numeric)
        YHI(varlist numeric)
        XORD(varlist numeric)
        YORD(varlist numeric)

        BOOTstrap
        ENVars(varlist numeric)
        REPS1(integer 100)
        REPS2(integer 2000)
        Level(integer 95)
    ];
    #del cr

    local outputs `varlist'

    if "`inputs'" == "" {
        di as error "inputs() required when using standard syntax"
        di as error "  Syntax: dea outputs, inputs(x1 x2) [options]"
        di as error "  Or use: dea x1 x2 = y1 y2, [options]"
        exit 198
    }

    // Reconcile ort() vs orientation()
    if "`ort'" != "" & "`orientation'" != "" {
        di as error "Specify only one of: ort() or orientation()"
        exit 198
    }
    if "`ort'" != "" local orientation "`ort'"

    // Delegate to common normalizer
    _dea_normalize_options,                           ///
        inputs(`inputs') outputs(`outputs')           ///
        dmu(`dmu') rts(`rts') orientation(`orientation') ///
        model(`model') stage(`stage') method(`method')   ///
        fdh(`fdh') tolerance(`tolerance')             ///
        generate(`generate') trace(`trace')           ///
        saving(`saving') replace(`replace') debug(`debug') ///
        engine(`engine') contractcheck(`contractcheck') ///
        reference(`reference') saveref(`saveref')     ///
        saveweights(`saveweights')                    ///
        w(`w') p(`p') gx(`gx') gy(`gy') ddfbeta(`ddfbeta') ///
        bad(`bad') direction(`direction')             ///
        bounded(`bounded') ordinal(`ordinal')         ///
        xlo(`xlo') xhi(`xhi') ylo(`ylo') yhi(`yhi')  ///
        xord(`xord') yord(`yord')                     ///
        bootstrap(`bootstrap') envars(`envars')       ///
        reps1(`reps1') reps2(`reps2') level(`level')  ///
        ifcond(`if') incond(`in')

    return add
end


********************************************************************************
* 3. OPTION NORMALIZER + VALIDATOR + ENGINE SELECTOR
********************************************************************************

program define _dea_normalize_options, rclass
    version 16.0

    syntax , Inputs(string) Outputs(string) ///
        [DMU(varname)]                      ///
        [RTS(string)]                       ///
        [ORIentation(string)]               ///
        [MODEL(string)]                     ///
        [STAGE(integer 2)]                  ///
        [METHod(string)]                    ///
        [FDH(string)]                       ///
        [TOLerance(real 1e-7)]              ///
        [GENerate(name)]                    ///
        [TRACE(string)]                     ///
        [SAVing(string)]                    ///
        [REPLACE(string)]                   ///
        [DEBUG(string)]                     ///
        [ENGine(string)]                    ///
        [CONTRACTCHECK(string)]             ///
        [REFerence(string)]                 ///
        [SAVERef(string)]                   ///
        [SAVEWeights(string)]               ///
        [W(varlist numeric)]                ///
        [P(varlist numeric)]                ///
        [GX(varlist numeric)]               ///
        [GY(varlist numeric)]               ///
        [DDFBETA(name)]                     ///
        [Bounded(varlist)]                  ///
        [ORDinal(varlist)]                  ///
        [XLO(varlist numeric)]              ///
        [XHI(varlist numeric)]              ///
        [YLO(varlist numeric)]              ///
        [YHI(varlist numeric)]              ///
        [XORD(varlist numeric)]             ///
        [YORD(varlist numeric)]             ///
        [BOOTstrap(string)]                 ///
        [ENVars(varlist numeric)]           ///
        [REPS1(integer 100)]                ///
        [REPS2(integer 2000)]               ///
        [Level(integer 95)]                 ///
        [IFCOND(string)]                    ///
        [INCOND(string)]                    ///
        [USING(string)]                     ///
        [BAD(varlist numeric)]              ///
        [DIRection(string)]

    // -----------------------------------------------------------------
    // 3.1 Normalize RTS
    // -----------------------------------------------------------------
    if "`rts'" == "" {
        if "`fdh'" != "" local rts "vrs"
        else local rts "crs"
    }
    else {
        local rts = lower("`rts'")
        if "`rts'" == "ccr" local rts "crs"
        if "`rts'" == "bcc" local rts "vrs"
        // Theory-facing aliases.  Retain drs/irs for backward compatibility,
        // but normalize the technology restrictions to the established
        // internal names used by the solver backends.
        if "`rts'" == "nirs" local rts "drs"
        if "`rts'" == "ndrs" local rts "irs"
        if !inlist("`rts'", "crs", "vrs", "drs", "irs") {
            di as error "rts() must be one of: crs (ccr), vrs (bcc), nirs (drs), ndrs (irs)"
            exit 198
        }
    }

    // -----------------------------------------------------------------
    // 3.2 Normalize Orientation
    // -----------------------------------------------------------------
    local user_set_ort = ("`orientation'" != "")   // 1=explicit, 0=default
    if "`orientation'" == "" {
        local orientation "input"
    }
    else {
        local orientation = lower("`orientation'")
        // Accept legacy dea.ado abbreviations
        if inlist("`orientation'", "i", "in", "inp", "inpu", "input") {
            local orientation "input"
        }
        else if inlist("`orientation'", "o", "out", "outp", "outpu", "output") {
            local orientation "output"
        }
        else {
            di as error "orientation() must be 'input' or 'output'"
            di as error "  (aliases: i/in/o/out)"
            exit 198
        }
    }

    // FDH does not allow rts()
    if "`fdh'" != "" & "`rts'" != "vrs" {
        if "`rts'" != "crs" {
            // User explicitly set non-VRS RTS with FDH
            di as error "FDH model does not support rts() option"
            exit 198
        }
        // Default: silently set to VRS for FDH
        local rts "vrs"
    }

    // -----------------------------------------------------------------
    // 3.3 Normalize Model
    // -----------------------------------------------------------------
    if "`model'" == "" local model "te"
    local model = lower("`model'")
    if !inlist("`model'", "te", "cost", "revenue", "profit", "ddf", "sbm", "additive", "supereff", "multiplier") & ///
       !inlist("`model'", "allocative", "supersbm", "fdh", "imprecise", "undesirable") & ///
       !inlist("`model'", "congestion", "sbmundesirable") {
        di as error "model() must be one of: te cost revenue profit ddf sbm additive" _n ///
                    "  supereff multiplier allocative supersbm fdh imprecise undesirable" _n ///
                    "  congestion sbmundesirable"
        exit 198
    }
    // model(fdh) implies fdh flag
    if "`model'" == "fdh" local fdh "fdh"

    // -----------------------------------------------------------------
    // 3.3a Validate model(undesirable) options
    // -----------------------------------------------------------------
    if "`model'" == "undesirable" {
        if "`bad'" == "" {
            di as error "model(undesirable) requires bad(varlist) for undesirable output variables"
            exit 198
        }
        // Validate direction
        if "`direction'" == "" local direction "yb"
        local direction = lower("`direction'")
        if !inlist("`direction'", "yb", "xyb", "x") {
            di as error "direction() must be yb (default), xyb, or x for model(undesirable)"
            exit 198
        }
    }

    // -----------------------------------------------------------------
    // 3.3b Validate model(sbmundesirable) options
    // -----------------------------------------------------------------
    if "`model'" == "sbmundesirable" {
        if "`bad'" == "" {
            di as error "model(sbmundesirable) requires bad(varlist) for undesirable output variables"
            exit 198
        }
    }

    // -----------------------------------------------------------------
    // 3.3c Validate model(congestion) options
    // -----------------------------------------------------------------
    if "`model'" == "congestion" {
        // Congestion only meaningful under VRS; silently enforce
        if "`rts'" != "" & "`rts'" != "vrs" {
            di as text "Note: model(congestion) requires VRS. rts() option ignored."
        }
        local rts "vrs"
    }

    // -----------------------------------------------------------------
    // 3.4 Validate Stage
    // -----------------------------------------------------------------
    if !inlist(`stage', 1, 2) {
        di as error "stage() must be 1 or 2"
        exit 198
    }
    if "`contractcheck'" != "" & ("`model'" != "te" | "`fdh'" != "" | ///
        "`bounded'" != "" | "`ordinal'" != "" | "`bootstrap'" != "") {
        di as error "contractcheck is an internal verification option for radial model(te) only"
        exit 198
    }

    // -----------------------------------------------------------------
    // 3.5 Normalize generate()
    // -----------------------------------------------------------------
    if "`ddfbeta'" != "" {
        if "`model'" != "ddf" {
            di as error "ddfbeta() is allowed only with model(ddf)"
            exit 198
        }
        if "`generate'" != "" {
            di as error "Specify only one of generate() or ddfbeta()"
            exit 198
        }
        local generate "`ddfbeta'"
    }
    if "`generate'" == "" local generate "dea_efficiency"

    // -----------------------------------------------------------------
    // 3.6 Engine Selection (version detection)
    // -----------------------------------------------------------------
    // Note: Both IPM and Simplex engines require Stata 16+
    // IPM uses LinearProgram() class; Simplex uses package-supplied Mata LP code.
    // Specialized dispatchers (FDH, imprecise, congestion, SBM-undesirable)
    // may use their own Mata routines after option normalization.
    
    if `c(stata_version)' < 16 {
        di as error "DEA v2.0 requires Stata 16 or later."
        di as error "Current version: `c(stata_version)'"
        di as error "For Stata 13-15, use DEA v1.2.2"
        exit 198
    }
    
    local engine_choice "simplex"
    local engine_label "Two-Phase Revised Simplex (Mata)"
    local auto_fallback ""

    // User override via engine() or method() (engine takes precedence)
    local user_engine = cond("`engine'" != "", "`engine'", "`method'")
    if "`user_engine'" != "" {
        local user_engine = lower("`user_engine'")
        if "`user_engine'" == "ipm" {
            local engine_choice "ipm"
            local engine_label "Interior Point Method (LinearProgram)"
        }
        else if "`user_engine'" == "simplex" {
            local engine_choice "simplex"
            local engine_label "Two-Phase Revised Simplex (Mata)"
        }
        else if "`user_engine'" == "auto" {
            // Auto: try simplex first, fallback to IPM on failure
            local engine_choice "simplex"
            local engine_label "Two-Phase Revised Simplex (Mata) [auto]"
            local auto_fallback "auto_fallback"
        }
        else {
            di as error "engine() must be 'ipm', 'simplex', or 'auto'"
            exit 198
        }
    }
    
    // Assign to engine variable
    local engine "`engine_choice'"

    // Both engines support all models (including undesirable) — no override needed.

    // Check: DRS/IRS with IPM need special handling (new in v2.0)
    // Simplex engine already supports DRS/IRS natively

    // -----------------------------------------------------------------
    // 3.7 Validate model-specific options
    // -----------------------------------------------------------------
    local m_in : word count `inputs'
    local m_out : word count `outputs'

    _dea_validate_model_opts, model(`model') ///
        m_in(`m_in') m_out(`m_out') ///
        w(`w') p(`p') gx(`gx') gy(`gy') ///
        xlo(`xlo') xhi(`xhi') ylo(`ylo') yhi(`yhi') ///
        xord(`xord') yord(`yord') ///
        bounded(`bounded') ordinal(`ordinal') ///
        reference(`reference') saveref(`saveref') saveweights(`saveweights')

    // -----------------------------------------------------------------
    // 3.8 Return if/in conditions for main program to handle touse
    // -----------------------------------------------------------------
    // Note: touse must be created in main program because tempvars
    // don't persist across program boundaries
    
    // Handle DMU variable placeholder
    if "`dmu'" == "" {
        local dmu ""  // Will be created in main program
    }

    // Count inputs/outputs
    local m_in : word count `inputs'
    local m_out : word count `outputs'

    // -----------------------------------------------------------------
    // 3.9 Return all parsed values
    // -----------------------------------------------------------------
    return local inputs      "`inputs'"
    return local outputs     "`outputs'"
    return local dmu         "`dmu'"
    return local rts         "`rts'"
    return local orientation "`orientation'"
    return local model       "`model'"
    return scalar stage      = `stage'
    return local method      "`method'"
    return scalar tolerance  = `tolerance'
    return local generate    "`generate'"
    return local trace       "`trace'"
    return local saving      "`saving'"
    return local replace     "`replace'"
    return local debug       "`debug'"
    return local contractcheck "`contractcheck'"
    return local fdh         "`fdh'"
    return local reference   "`reference'"
    return local saveref     "`saveref'"
    return local saveweights "`saveweights'"
    return local w           "`w'"
    return local p           "`p'"
    return local gx          "`gx'"
    return local gy          "`gy'"
    return local ddfbeta     "`ddfbeta'"
    return local bad         "`bad'"
    return local direction   "`direction'"
    return local xlo         "`xlo'"
    return local xhi         "`xhi'"
    return local ylo         "`ylo'"
    return local yhi         "`yhi'"
    return local xord        "`xord'"
    return local yord        "`yord'"
    return local bounded     "`bounded'"
    return local ordinal     "`ordinal'"
    return local bootstrap   "`bootstrap'"
    return local envars      "`envars'"
    return scalar reps1      = `reps1'
    return scalar reps2      = `reps2'
    return scalar level      = `level'
    return local ifcond      "`ifcond'"
    return local incond      "`incond'"
    return scalar m_in       = `m_in'
    return scalar m_out      = `m_out'
    return local engine      "`engine'"
    return local engine_label "`engine_label'"
    return local auto_fallback "`auto_fallback'"
    return scalar user_set_ort = `user_set_ort'
end


********************************************************************************
* 4. MODEL-SPECIFIC OPTION VALIDATOR
********************************************************************************

program define _dea_validate_model_opts
    syntax , MODEL(string) m_in(integer) m_out(integer) ///
        [W(varlist numeric)] [P(varlist numeric)]       ///
        [GX(varlist numeric)] [GY(varlist numeric)]     ///
        [XLO(varlist numeric)] [XHI(varlist numeric)]   ///
        [YLO(varlist numeric)] [YHI(varlist numeric)]   ///
        [XORD(varlist numeric)] [YORD(varlist numeric)] ///
        [Bounded(varlist)] [ORDinal(varlist)] ///
        [REFerence(string)] [SAVERef(string)] [SAVEWeights(string)]

    // Cost model requires input prices w()
    if "`model'" == "cost" {
        if "`w'" == "" {
            di as error "model(cost) requires w() for input prices"
            exit 198
        }
        local nw : word count `w'
        if `nw' != `m_in' {
            di as error "w() must have `m_in' variables (same as inputs)"
            exit 198
        }
    }

    // Revenue model requires output prices p()
    if "`model'" == "revenue" {
        if "`p'" == "" {
            di as error "model(revenue) requires p() for output prices"
            exit 198
        }
        local np : word count `p'
        if `np' != `m_out' {
            di as error "p() must have `m_out' variables (same as outputs)"
            exit 198
        }
    }

    // Profit model requires both input prices w() and output prices p()
    if "`model'" == "profit" {
        if "`w'" == "" {
            di as error "model(profit) requires w() for input prices"
            exit 198
        }
        if "`p'" == "" {
            di as error "model(profit) requires p() for output prices"
            exit 198
        }
        local nw : word count `w'
        if `nw' != `m_in' {
            di as error "w() must have `m_in' variables (same as inputs)"
            exit 198
        }
        local np : word count `p'
        if `np' != `m_out' {
            di as error "p() must have `m_out' variables (same as outputs)"
            exit 198
        }
    }

    // DDF: direction vectors gx() and gy() are optional
    // If not specified, observed values are used as default direction
    if "`model'" == "ddf" {
        if "`gx'" != "" {
            local ngx : word count `gx'
            if `ngx' != `m_in' {
                di as error "gx() must have `m_in' variables (same as inputs)"
                exit 198
            }
        }
        if "`gy'" != "" {
            local ngy : word count `gy'
            if `ngy' != `m_out' {
                di as error "gy() must have `m_out' variables (same as outputs)"
                exit 198
            }
        }
    }

    // Allocative requires w() or p()
    if "`model'" == "allocative" {
        if "`w'" == "" & "`p'" == "" {
            di as error "model(allocative) requires w() for cost or p() for revenue"
            exit 198
        }
    }

    // Imprecise DEA validation
    if "`model'" == "imprecise" {
        if "`xlo'" == "" & "`xhi'" == "" & "`ylo'" == "" & "`yhi'" == "" & ///
           "`xord'" == "" & "`yord'" == "" & ///
           "`bounded'" == "" & "`ordinal'" == "" {
            di as error "model(imprecise) requires bounds or ordinal specifications"
            exit 198
        }
    }

    // saveref requires reference
    if "`saveref'" != "" & "`reference'" == "" {
        di as error "saveref requires the reference option"
        exit 198
    }

    // saveweights only for multiplier model
    if "`saveweights'" != "" & "`model'" != "multiplier" {
        di as error "saveweights only available with model(multiplier)"
        exit 198
    }
end


********************************************************************************
* 5. MODEL DISPATCHER
********************************************************************************

program define _dea_model_dispatch, rclass
    syntax , Inputs(string) Outputs(string) DMU(varname) ///
        RTS(string) ORIentation(string) MODEL(string) ///
        STAGE(integer) ENGINE(string) TOUSE(varname) ///
        GENerate(name) TOLerance(real) ///
        [METHod(string)] ///
        [W(varlist numeric)] [P(varlist numeric)] ///
        [GX(varlist numeric)] [GY(varlist numeric)] ///
        [DDFBETA(name)] ///
        [BAD(varlist numeric)] [DIRection(string)] ///
        [REFerence] [SAVERef] [SAVEWeights] ///
        [TRACE] [DEBUG] ///
        [SAVing(string)] [REPLACE] ///
        [SBMORT(integer 0)]

    // =================================================================
    // Pre-create lambda variables in this scope (if reference requested)
    // so they survive through _dea_display.
    // tempvars created in sub-programs are destroyed when those programs
    // return, causing "not found" errors during reference set display.
    // =================================================================
    local lambda_vars ""
    if "`reference'" != "" {
        qui count if `touse'
        local _n_for_lam = r(N)
        forvalues j = 1/`_n_for_lam' {
            tempvar _lam`j'
            qui gen double `_lam`j'' = . if `touse'
            local lambda_vars "`lambda_vars' `_lam`j''"
        }
    }

    // =================================================================
    // Pre-create slack vars for TE Stage 2 (so they survive engine return)
    // =================================================================
    local slack_in_vars ""
    local slack_out_vars ""
    if "`model'" == "te" & `stage' == 2 & "`engine'" != "ipm" {
        local _m_in : word count `inputs'
        local _m_out : word count `outputs'
        forvalues j = 1/`_m_in' {
            tempvar _si`j'
            qui gen double `_si`j'' = 0 if `touse'
            local slack_in_vars "`slack_in_vars' `_si`j''"
        }
        forvalues j = 1/`_m_out' {
            tempvar _so`j'
            qui gen double `_so`j'' = 0 if `touse'
            local slack_out_vars "`slack_out_vars' `_so`j''"
        }
    }

    // =================================================================
    // Pre-create multiplier weight variables in this scope.
    // Weight variables are meaningful for the stage-2 multiplier optimum.
    // For stage(1), saveweights is intentionally ignored and only the
    // efficiency score is reported.
    // =================================================================
    local weight_in_vars ""
    local weight_out_vars ""
    local weight_u0_var ""
    if "`model'" == "multiplier" & "`saveweights'" != "" & `stage' == 2 {
        local _m_in : word count `inputs'
        local _m_out : word count `outputs'
        forvalues j = 1/`_m_in' {
            local _wv "_dea_v`j'"
            capture drop `_wv'
            qui gen double `_wv' = . if `touse'
            local _xname : word `j' of `inputs'
            label var `_wv' "Multiplier input weight for `_xname'"
            local weight_in_vars "`weight_in_vars' `_wv'"
        }
        forvalues j = 1/`_m_out' {
            local _wv "_dea_u`j'"
            capture drop `_wv'
            qui gen double `_wv' = . if `touse'
            local _yname : word `j' of `outputs'
            label var `_wv' "Multiplier output weight for `_yname'"
            local weight_out_vars "`weight_out_vars' `_wv'"
        }
        local weight_u0_var "_dea_u0"
        capture drop `weight_u0_var'
        qui gen double `weight_u0_var' = . if `touse'
        label var `weight_u0_var' "Multiplier intercept: u'y - v'x + u0 <= 0"
    }

    // =================================================================
    // Route to appropriate engine
    // =================================================================
    if "`engine'" == "ipm" {
        // IPM engine: _dea_ipm.ado handles all models
        _dea_ipm_solve, inputs(`inputs') outputs(`outputs') ///
            dmu(`dmu') rts(`rts') orientation(`orientation') ///
            model(`model') stage(`stage') touse(`touse') ///
            generate(`generate') tolerance(`tolerance') ///
            w(`w') p(`p') gx(`gx') gy(`gy') ddfbeta(`ddfbeta') ///
            bad(`bad') direction(`direction') ///
            `reference' `saveref' `saveweights' ///
            lambda_vars(`lambda_vars') ///
            weight_in_vars(`weight_in_vars') ///
            weight_out_vars(`weight_out_vars') ///
            weight_u0_var(`weight_u0_var') ///
            sbmort(`sbmort') ///
            `trace' `debug'
    }
    else {
        // Simplex engine: _dea_simplex.ado handles all models
        _dea_simplex_solve, inputs(`inputs') outputs(`outputs') ///
            dmu(`dmu') rts(`rts') orientation(`orientation') ///
            model(`model') stage(`stage') touse(`touse') ///
            generate(`generate') tolerance(`tolerance') ///
            w(`w') p(`p') gx(`gx') gy(`gy') ddfbeta(`ddfbeta') ///
            bad(`bad') direction(`direction') ///
            `reference' `saveref' `saveweights' ///
            lambda_vars(`lambda_vars') ///
            slack_in_vars(`slack_in_vars') ///
            slack_out_vars(`slack_out_vars') ///
            weight_in_vars(`weight_in_vars') ///
            weight_out_vars(`weight_out_vars') ///
            weight_u0_var(`weight_u0_var') ///
            sbmort(`sbmort') ///
            `trace' `debug'
    }

    // =================================================================
    // Display results
    // =================================================================
    // lambda_vars and slack_in/out_vars already live in this scope
    // (lambda: pre-created above; slacks: pre-created or captured from IPM)

    _dea_display, model(`model') rts(`rts') orientation(`orientation') ///
        generate(`generate') touse(`touse') ///
        tolerance(`tolerance') stage(`stage') ///
        dmu(`dmu') `reference' `debug' ///
        sbmort(`sbmort') ///
        lambda_vars(`lambda_vars') ///
        slack_in_vars(`slack_in_vars') ///
        slack_out_vars(`slack_out_vars') ///
        inputs(`inputs') outputs(`outputs')

    // =================================================================
    // Capture display statistics BEFORE any mkmat call may clear r()
    // (mkmat is rclass and resets r(); must save stats first)
    // =================================================================
    local _n_eff_save    = r(n_efficient)
    local _mean_eff_save = r(mean_efficiency)
    local _min_eff_save  = r(min_eff)
    local _max_eff_save  = r(max_eff)
    local _n_failed_save = r(n_failed)

    // =================================================================
    // Save results if requested
    // =================================================================
    if "`saving'" != "" {
        _dea_save, saving(`saving') `replace' ///
            inputs(`inputs') outputs(`outputs') ///
            dmu(`dmu') generate(`generate') touse(`touse') ///
            weight_in_vars(`weight_in_vars') ///
            weight_out_vars(`weight_out_vars') ///
            weight_u0_var(`weight_u0_var') ///
            `saveref'
    }

    // =================================================================
    // Return matrices (lambda, slack_in, slack_out)
    // NOTE: return matrix must come BEFORE return add (Stata 16)
    // Tempvars in _dea_model_dispatch scope survive until after return add.
    // =================================================================
    if "`lambda_vars'" != "" {
        tempname _lam_mat
        capture qui mkmat `lambda_vars' if `touse', matrix(`_lam_mat')
        if _rc == 0 {
            return matrix lambda = `_lam_mat'
        }
    }
    if "`slack_in_vars'" != "" {
        tempname _sin_mat
        capture qui mkmat `slack_in_vars' if `touse', matrix(`_sin_mat')
        if _rc == 0 {
            return matrix slack_in = `_sin_mat'
        }
    }
    if "`slack_out_vars'" != "" {
        tempname _sout_mat
        capture qui mkmat `slack_out_vars' if `touse', matrix(`_sout_mat')
        if _rc == 0 {
            return matrix slack_out = `_sout_mat'
        }
    }
    if "`weight_u0_var'" != "" {
        tempname _u0_mat
        capture qui mkmat `weight_u0_var' if `touse', matrix(`_u0_mat')
        if _rc == 0 {
            return matrix multiplier_u0 = `_u0_mat'
        }
    }

    // Explicitly return display statistics (mkmat may have cleared r())
    return scalar n_efficient     = `_n_eff_save'
    return scalar mean_efficiency = `_mean_eff_save'
    return scalar min_eff         = `_min_eff_save'
    return scalar max_eff         = `_max_eff_save'
    return scalar n_failed        = `_n_failed_save'
end


********************************************************************************
* 6. FDH DISPATCHER (Branch & Bound, engine-independent)
********************************************************************************

program define _dea_fdh_dispatch, rclass
    syntax , Inputs(string) Outputs(string) DMU(varname) ///
        ORIentation(string) STAGE(integer) ///
        TOUSE(varname) GENerate(name) TOLerance(real) ///
        [TRACE] [DEBUG] [REFerence] [SAVERef] ///
        [SAVing(string)] [REPLACE]

    // FDH uses dominance comparison; no LP solver is needed.
    // The Mata implementation is loaded from _dea_ipm.ado.
    
    // Ensure Mata functions are loaded
    capture mata: mata which _dea_fdh()
    if _rc != 0 {
        capture findfile _dea_ipm.ado
        if _rc == 0 {
            qui run "`r(fn)'"
        }
        else {
            di as error "_dea_ipm.ado not found. Please install it in your ado path."
            exit 601
        }
    }

    di as text _n "FDH (Free Disposal Hull) — Dominance Comparison"
    
    // Create efficiency variable
    tempvar efficiency
    qui gen double `efficiency' = . if `touse'
    
    // Call Mata FDH solver
    quietly mata: _dea_fdh("`inputs'", "`outputs'", "`efficiency'", "`touse'", "`orientation'")
    
    // Copy to final variable
    qui gen double `generate' = `efficiency' if `touse'
    
    // Count results
    qui count if `touse'
    local n_dmus = r(N)
    
    qui count if abs(`generate' - 1) < `tolerance' & `touse'
    local n_efficient = r(N)
    
    // Display results
    di as text "{hline 50}"
    di as text "Model: FDH `orientation'"
    di as text "DMUs evaluated:  " as result `n_dmus'
    di as text "Efficient DMUs:  " as result `n_efficient' ///
        as text " (" %5.1f 100*`n_efficient'/`n_dmus' "%)"
    
    qui sum `generate' if `touse', detail
    di as text "Mean efficiency:  " as result %9.6f r(mean)
    di as text "Min efficiency:   " as result %9.6f r(min)
    di as text "Max efficiency:   " as result %9.6f r(max)
    di as text "{hline 50}"
    
    // Display DMU results
    di as text _n "{bf:DMU Efficiency Scores}"
    di as text "{hline 50}"
    
    tempvar eff_rounded rank_var
    qui gen double `eff_rounded' = round(`generate', `tolerance') if `touse'
    qui egen `rank_var' = rank(-`eff_rounded') if `touse', track
    
    preserve
    qui keep if `touse'
    local nobs = _N
    di as text %12s "DMU" " {c |}" %12s "Theta" %8s "Rank"
    di as text "{hline 13}{c +}{hline 20}"
    
    forvalues i = 1/`nobs' {
        local dmu_val = `dmu'[`i']
        local eff_val = `generate'[`i']
        local rank_val = `rank_var'[`i']
        if "`eff_val'" != "." {
            di as text %12s "`dmu_val'" " {c |}" ///
                as result %12.6f `eff_val' %8.0f `rank_val'
        }
        else {
            di as text %12s "`dmu_val'" " {c |}" ///
                as result %12s "FAILED" %8s "."
        }
    }
    di as text "{hline 34}"
    restore
    
    // Return values
    return scalar N = `n_dmus'
    return scalar n_efficient = `n_efficient'
    return local model "fdh"
    return local orientation "`orientation'"
    return add
end


********************************************************************************
* 7. IMPRECISE DEA DISPATCHER
********************************************************************************

program define _dea_imprecise_dispatch, rclass
    syntax , Inputs(string) Outputs(string) DMU(varname) ///
        RTS(string) ORIentation(string) STAGE(integer) ///
        ENGINE(string) TOUSE(varname) ///
        GENerate(name) TOLerance(real) ///
        [Bounded(varlist)] [ORDinal(varlist)] ///
        [XLO(varlist numeric)] [XHI(varlist numeric)] ///
        [YLO(varlist numeric)] [YHI(varlist numeric)] ///
        [XORD(varlist numeric)] [YORD(varlist numeric)] ///
        [TRACE] [DEBUG] ///
        [SAVing(string)] [REPLACE]

    di as text _n "IDEA (Imprecise DEA)"
    di as text "{hline 50}"
    
    // Ensure Mata functions are loaded
    capture mata: mata which _dea_idea()
    if _rc != 0 {
        capture findfile _dea_ipm.ado
        if _rc == 0 {
            qui run "`r(fn)'"
        }
        else {
            di as error "_dea_ipm.ado not found."
            exit 601
        }
    }
    
    // Current IDEA solver uses interval bounds directly. Other imprecise-data
    // option names are parsed for interface compatibility and future extension.
    // Validate interval data specification.
    local m_in : word count `inputs'
    local m_out : word count `outputs'
    
    if "`xlo'" != "" {
        local n_xlo : word count `xlo'
        if `n_xlo' != `m_in' {
            di as error "xlo() must have `m_in' variables (same as inputs)"
            exit 198
        }
    }
    if "`xhi'" != "" {
        local n_xhi : word count `xhi'
        if `n_xhi' != `m_in' {
            di as error "xhi() must have `m_in' variables (same as inputs)"
            exit 198
        }
    }
    if "`ylo'" != "" {
        local n_ylo : word count `ylo'
        if `n_ylo' != `m_out' {
            di as error "ylo() must have `m_out' variables (same as outputs)"
            exit 198
        }
    }
    if "`yhi'" != "" {
        local n_yhi : word count `yhi'
        if `n_yhi' != `m_out' {
            di as error "yhi() must have `m_out' variables (same as outputs)"
            exit 198
        }
    }
    
    // Create efficiency variables
    tempvar eff_pess eff_opt
    qui gen double `eff_pess' = . if `touse'
    qui gen double `eff_opt' = . if `touse'
    
    // Display settings
    di as text "Model: IDEA `rts' `orientation'"
    di as text "Inputs:  `inputs'"
    di as text "Outputs: `outputs'"

    // Legacy IDEA path from dea_original.ado: bounded and ordinal variables
    // are transformed DMU-by-DMU, then the ordinary DEA model is solved on the
    // transformed data.  The integrated interval solver below handles xlo/xhi
    // and ylo/yhi; this path restores the older bounded/ordinal behavior.
    if ("`bounded'" != "" | "`ordinal'" != "" | ///
        "`xord'" != "" | "`yord'" != "") & ///
       ("`xlo'" == "" & "`xhi'" == "" & "`ylo'" == "" & "`yhi'" == "") {
        _dea_imprecise_legacy, inputs(`inputs') outputs(`outputs') ///
            dmu(`dmu') rts(`rts') orientation(`orientation') ///
            stage(`stage') engine(`engine') touse(`touse') ///
            generate(`generate') tolerance(`tolerance') ///
            bounded(`bounded') ordinal(`ordinal' `xord' `yord') ///
            `trace' `debug'
        return add
        exit
    }

    if "`xlo'" != "" | "`xhi'" != "" {
        di as text "Input bounds:  xlo(`xlo') xhi(`xhi')"
    }
    if "`ylo'" != "" | "`yhi'" != "" {
        di as text "Output bounds: ylo(`ylo') yhi(`yhi')"
    }
    
    // Call Mata IDEA solver
    quietly mata: _dea_idea("`inputs'", "`outputs'", ///
        "`xlo'", "`xhi'", "`ylo'", "`yhi'", ///
        "`eff_pess'", "`eff_opt'", ///
        "`touse'", "`orientation'", "`rts'", `tolerance')
    
    // Create output variables
    // Main generate variable = pessimistic (conservative) efficiency
    qui gen double `generate' = `eff_pess' if `touse'
    
    // Also create optimistic variable
    local opt_var = "`generate'_opt"
    capture drop `opt_var'
    qui gen double `opt_var' = `eff_opt' if `touse'
    
    // Count results
    qui count if `touse'
    local n_dmus = r(N)
    
    qui count if abs(`generate' - 1) < `tolerance' & `touse'
    local n_eff_pess = r(N)
    
    qui count if abs(`eff_opt' - 1) < `tolerance' & `touse'
    local n_eff_opt = r(N)
    
    // Summary statistics
    qui sum `generate' if `touse', detail
    local mean_pess = r(mean)
    local min_pess = r(min)
    local max_pess = r(max)
    
    qui sum `eff_opt' if `touse', detail
    local mean_opt = r(mean)
    local min_opt = r(min)
    local max_opt = r(max)
    
    // Display results
    di as text "{hline 50}"
    di as text "DMUs evaluated:  " as result `n_dmus'
    di as text ""
    di as text "{bf:Pessimistic Efficiency} (stored in `generate')"
    di as text "  Efficient DMUs:  " as result `n_eff_pess' ///
        as text " (" %5.1f 100*`n_eff_pess'/`n_dmus' "%)"
    di as text "  Mean:  " as result %9.6f `mean_pess'
    di as text "  Range: " as result %9.6f `min_pess' " - " %9.6f `max_pess'
    di as text ""
    di as text "{bf:Optimistic Efficiency} (stored in `opt_var')"
    di as text "  Efficient DMUs:  " as result `n_eff_opt' ///
        as text " (" %5.1f 100*`n_eff_opt'/`n_dmus' "%)"
    di as text "  Mean:  " as result %9.6f `mean_opt'
    di as text "  Range: " as result %9.6f `min_opt' " - " %9.6f `max_opt'
    di as text "{hline 50}"
    
    // Display DMU results
    di as text _n "{bf:DMU Efficiency Scores}"
    di as text "{hline 60}"
    di as text %12s "DMU" " {c |}" %14s "Pessimistic" %14s "Optimistic" %10s "Interval"
    di as text "{hline 13}{c +}{hline 46}"
    
    preserve
    qui keep if `touse'
    local nobs = _N
    
    forvalues i = 1/`nobs' {
        local dmu_val = `dmu'[`i']
        local pess_val = `generate'[`i']
        local opt_val = `eff_opt'[`i']
        if "`pess_val'" != "." & "`opt_val'" != "." {
            local interval = `opt_val' - `pess_val'
            di as text %12s "`dmu_val'" " {c |}" ///
                as result %14.6f `pess_val' %14.6f `opt_val' %10.6f `interval'
        }
    }
    di as text "{hline 60}"
    restore
    
    // Return values
    return scalar N = `n_dmus'
    return scalar n_efficient_pess = `n_eff_pess'
    return scalar n_efficient_opt = `n_eff_opt'
    return scalar mean_pess = `mean_pess'
    return scalar mean_opt = `mean_opt'
    return local model "idea"
    return local orientation "`orientation'"
    return local rts "`rts'"
    return add
end

capture program drop _dea_imprecise_legacy
program define _dea_imprecise_legacy, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) DMU(varname) ///
        RTS(string) ORIentation(string) STAGE(integer) ///
        ENGINE(string) TOUSE(varname) GENerate(name) TOLerance(real) ///
        [Bounded(varlist)] [ORDinal(varlist)] ///
        [TRACE] [DEBUG]

    tempvar __idea_row
    qui gen long `__idea_row' = _n

    tempfile __idea_base __idea_results __one

    preserve
        qui keep if `touse'
        qui keep `__idea_row' `dmu' `inputs' `outputs' `bounded' `ordinal'
        qui save `__idea_base', replace

        clear
        qui set obs 0
        qui gen long `__idea_row' = .
        qui gen double __idea_eff = .
        qui save `__idea_results', replace

        use `__idea_base', clear
        local __n = _N

        forvalues dmui = 1/`__n' {
            use `__idea_base', clear

            if "`bounded'" != "" {
                foreach boundedVar of varlist `bounded' {
                    capture drop `boundedVar'_2
                    capture confirm string variable `boundedVar'
                    if !_rc {
                        tempvar __comma __lo __hi
                        qui gen int `__comma' = strpos(`boundedVar', ",")
                        qui gen double `__lo' = real(strtrim(subinstr(substr(`boundedVar', 2, `__comma' - 2), "[", "", .)))
                        qui gen double `__hi' = real(strtrim(subinstr(subinstr(substr(`boundedVar', `__comma' + 1, .), "]", "", .), "[", "", .)))
                        qui drop `boundedVar'
                        qui gen double `boundedVar' = `__lo'
                        qui gen double `boundedVar'_2 = `__hi'
                    }
                    else {
                        qui gen double `boundedVar'_2 = `boundedVar'
                    }

                    if strpos(" `inputs' ", " `boundedVar' ") {
                        // For bounded inputs, the evaluated DMU uses the lower
                        // value and reference DMUs use the upper value.
                        qui replace `boundedVar' = `boundedVar'_2 if _n != `dmui'
                    }
                    else if strpos(" `outputs' ", " `boundedVar' ") {
                        // For bounded outputs, the evaluated DMU uses the upper
                        // value and reference DMUs use the lower value.
                        qui replace `boundedVar' = `boundedVar'_2 if _n == `dmui'
                    }
                    qui drop `boundedVar'_2
                }
            }

            if "`ordinal'" != "" {
                foreach ordVar of varlist `ordinal' {
                    if strpos(" `inputs' ", " `ordVar' ") {
                        qui replace `ordVar' = . if `ordVar' < `ordVar'[`dmui']
                        qui replace `ordVar' = 1 if `ordVar' < .
                        qui replace `ordVar' = _N if `ordVar' == .
                    }
                    else if strpos(" `outputs' ", " `ordVar' ") {
                        qui replace `ordVar' = . if `ordVar' <= `ordVar'[`dmui']
                        qui replace `ordVar' = 0 if `ordVar' < .
                        qui replace `ordVar' = 1 if `ordVar' == .
                    }
                }
            }

            capture drop __idea_score
            quietly dea `inputs' = `outputs', rts(`rts') ort(`orientation') ///
                stage(`stage') gen(__idea_score) engine(`engine')

            qui keep in `dmui'
            qui keep `__idea_row' __idea_score
            qui rename __idea_score __idea_eff
            qui save `__one', replace

            use `__idea_results', clear
            qui append using `__one'
            qui save `__idea_results', replace
        }
    restore

    qui merge 1:1 `__idea_row' using `__idea_results', nogen

    capture drop `generate'
    qui gen double `generate' = __idea_eff if `touse'

    local opt_var "`generate'_opt"
    local mid_var "`generate'_mid"
    capture drop `opt_var'
    capture drop `mid_var'
    qui gen double `opt_var' = __idea_eff if `touse'
    qui gen double `mid_var' = __idea_eff if `touse'

    qui count if `touse'
    local n_dmus = r(N)
    qui count if abs(`generate' - 1) < `tolerance' & `touse'
    local n_eff = r(N)
    qui sum `generate' if `touse', detail
    local mean_eff = r(mean)
    local min_eff = r(min)
    local max_eff = r(max)

    di as text "Legacy bounded/ordinal IDEA transformation"
    if "`bounded'" != "" di as text "Bounded variables: `bounded'"
    if "`ordinal'" != "" di as text "Ordinal variables: `ordinal'"
    di as text "{hline 50}"
    di as text "DMUs evaluated:  " as result `n_dmus'
    di as text "Efficient DMUs:  " as result `n_eff' ///
        as text " (" %5.1f 100*`n_eff'/`n_dmus' "%)"
    di as text "Mean:  " as result %9.6f `mean_eff'
    di as text "Range: " as result %9.6f `min_eff' " - " %9.6f `max_eff'
    di as text "{hline 50}"

    return scalar N = `n_dmus'
    return scalar n_efficient = `n_eff'
    return scalar mean_efficiency = `mean_eff'
    return scalar min_eff = `min_eff'
    return scalar max_eff = `max_eff'
    return local model "idea_legacy"
    return local orientation "`orientation'"
    return local rts "`rts'"
    return add
end


********************************************************************************
* 8. VRS FRONTIER ANALYSIS (Scale Efficiency + RTS Determination)
********************************************************************************

program define _dea_vrsfrontier, rclass
    syntax , Inputs(string) Outputs(string) DMU(varname) ///
        ORIentation(string) ENGINE(string) TOUSE(varname) ///
        TOLerance(real) [TRACE] [DEBUG]

    // VRS Frontier: Run CRS → VRS → DRS sequentially
    // Then compute:
    //   SCALE = CRS_TE / VRS_TE
    //   RTS   = Färe-Grosskopf-Lovell (1993) method:
    //           if CRS_TE == VRS_TE → CRS
    //           else if DRS_TE == VRS_TE → DRS
    //           else → IRS
    //   Also: sum(lambda) method from CRS solution

    // Create result variables
    tempvar crs_te vrs_te drs_te scale_eff rts_type lambda_sum
    qui gen double `crs_te'    = . if `touse'
    qui gen double `vrs_te'    = . if `touse'
    qui gen double `drs_te'    = . if `touse'
    qui gen double `scale_eff' = . if `touse'
    qui gen double `rts_type'  = . if `touse'
    qui gen double `lambda_sum' = . if `touse'

    // Call Mata VRS Frontier solver
    if "`engine'" == "simplex" {
        quietly mata: _dea_simplex_envelopment("`inputs'", "`outputs'", ///
            "`crs_te'", "`touse'", "`orientation'", "crs", ///
            `tolerance', "")
        quietly mata: _dea_simplex_envelopment("`inputs'", "`outputs'", ///
            "`vrs_te'", "`touse'", "`orientation'", "vrs", ///
            `tolerance', "")
        quietly mata: _dea_simplex_envelopment("`inputs'", "`outputs'", ///
            "`drs_te'", "`touse'", "`orientation'", "drs", ///
            `tolerance', "")

        qui replace `scale_eff' = `crs_te' / `vrs_te' ///
            if `touse' & `crs_te' < . & `vrs_te' < . & `vrs_te' > 0
        qui replace `rts_type' = 0 if `touse' & ///
            `crs_te' < . & `vrs_te' < . & abs(`crs_te' - `vrs_te') < `tolerance'
        qui replace `rts_type' = -1 if `touse' & missing(`rts_type') & ///
            `drs_te' < . & `vrs_te' < . & abs(`drs_te' - `vrs_te') < `tolerance'
        qui replace `rts_type' = 1 if `touse' & missing(`rts_type') & ///
            `vrs_te' < .
    }
    else {
        quietly mata: _dea_ipm_vrsfrontier("`inputs'", "`outputs'", ///
            "`crs_te'", "`vrs_te'", "`drs_te'", ///
            "`scale_eff'", "`rts_type'", "`lambda_sum'", ///
            "`touse'", "`orientation'", `tolerance')
    }

    // Count DMUs
    qui count if `touse'
    local n_dmus = r(N)

    // Build result matrix
    tempname vrsfrontier
    matrix `vrsfrontier' = J(`n_dmus', 5, 0)

    local i = 0
    qui levelsof `dmu' if `touse', local(dmu_list)
    local _dmu_is_string = substr("`: type `dmu''", 1, 3) == "str"
    foreach d of local dmu_list {
        local ++i
        // levelsof returns text in a local even when the source variable is
        // numeric.  Quote only string identifiers when forming the qualifier.
        if `_dmu_is_string' {
            local _dmu_match `"`dmu' == "`d'""'
        }
        else {
            local _dmu_match "`dmu' == `d'"
        }
        qui summ `crs_te'   if `_dmu_match' & `touse', meanonly
        matrix `vrsfrontier'[`i', 1] = r(mean)
        qui summ `vrs_te'   if `_dmu_match' & `touse', meanonly
        matrix `vrsfrontier'[`i', 2] = r(mean)
        qui summ `drs_te'   if `_dmu_match' & `touse', meanonly
        matrix `vrsfrontier'[`i', 3] = r(mean)
        qui summ `scale_eff' if `_dmu_match' & `touse', meanonly
        matrix `vrsfrontier'[`i', 4] = r(mean)
        qui summ `rts_type'  if `_dmu_match' & `touse', meanonly
        matrix `vrsfrontier'[`i', 5] = r(mean)
    }
    matrix colnames `vrsfrontier' = CRS_TE VRS_TE DRS_TE SCALE RTS

    // Attempt to set row names from DMU values
    capture {
        local rownames ""
        foreach d of local dmu_list {
            local rownames `rownames' `d'
        }
        matrix rownames `vrsfrontier' = `rownames'
    }

    // Display
    di as text _n "VRS Frontier Analysis"
    di as text "{hline 65}"
    di as text "  RTS codes: -1 = DRS, 0 = CRS, 1 = IRS"
    di as text "{hline 65}"
    matrix list `vrsfrontier', noblank nohalf noheader format(%9.6f)
    di as text "{hline 65}"

    // Summary
    qui count if `rts_type' == 0 & `touse'
    local n_crs = r(N)
    qui count if `rts_type' == -1 & `touse'
    local n_drs = r(N)
    qui count if `rts_type' == 1 & `touse'
    local n_irs = r(N)
    di as text "  CRS: `n_crs', DRS: `n_drs', IRS: `n_irs'"

    // Generate permanent variables
    capture drop CRS_TE VRS_TE SCALE RTS_NUM RTS
    qui gen double CRS_TE = `crs_te' if `touse'
    qui gen double VRS_TE = `vrs_te' if `touse'
    qui gen double SCALE  = `scale_eff' if `touse'
    qui gen double RTS_NUM = `rts_type' if `touse'
    qui gen str3 RTS = ""
    qui replace RTS = "drs" if RTS_NUM == -1
    qui replace RTS = "-"   if RTS_NUM == 0
    qui replace RTS = "irs" if RTS_NUM == 1
    format CRS_TE VRS_TE SCALE %9.6f
    label var CRS_TE "CRS Technical Efficiency"
    label var VRS_TE "VRS Technical Efficiency"
    label var SCALE  "Scale Efficiency (CRS_TE/VRS_TE)"
    label var RTS    "Returns to Scale (drs/crs/irs)"

    // Return
    return matrix vrsfrontier = `vrsfrontier'
    return scalar n_crs = `n_crs'
    return scalar n_drs = `n_drs'
    return scalar n_irs = `n_irs'
    return add
end


********************************************************************************
* 9. RESULT DISPLAY
********************************************************************************

program define _dea_display, rclass
    syntax , MODEL(string) RTS(string) ORIentation(string) ///
        GENerate(name) TOUSE(varname) TOLerance(real) ///
        STAGE(integer) DMU(varname) ///
        [REFerence] [DEBUG] ///
        [lambda_vars(string)] [slack_in_vars(string)] [slack_out_vars(string)] ///
        [Inputs(string)] [Outputs(string)] [SBMORT(integer 0)]

    // Summary statistics
    qui sum `generate' if `touse', detail
    local mean_eff = r(mean)
    local min_eff  = r(min)
    local max_eff  = r(max)
    local sd_eff   = r(sd)

    qui count if `touse'
    local n_dmus = r(N)

    // Count efficient DMUs (model-dependent)
    if inlist("`model'", "additive") {
        qui count if abs(`generate') < `tolerance' & `touse'
    }
    else if inlist("`model'", "supereff", "supersbm") {
        if "`orientation'" == "input" {
            qui count if `generate' > 1 + `tolerance' & `touse' & !missing(`generate')
        }
        else {
            qui count if `generate' < 1 - `tolerance' & `touse' & !missing(`generate')
        }
    }
    else {
        qui count if abs(`generate' - 1) < `tolerance' & `touse'
    }
    local n_efficient = r(N)

    // Failed DMUs
    qui count if missing(`generate') & `touse'
    local n_failed = r(N)

    // Display
    di as text _n "DEA Results"
    di as text "{hline 50}"
    if "`model'" == "sbm" & `sbmort' == 1 {
        di as text "Model: `rts' input-oriented (sbm)"
    }
    else if "`model'" == "sbm" & `sbmort' == 2 {
        di as text "Model: `rts' output-oriented (sbm)"
    }
    else if inlist("`model'", "sbm", "supersbm", "profit") {
        di as text "Model: `rts' non-oriented (`model')"
    }
    else {
        di as text "Model: `rts' `orientation' (`model')"
    }
    di as text "Stage: `stage'"
    di as text "DMUs evaluated:  " as result `n_dmus'
    di as text "Efficient DMUs:  " as result `n_efficient' ///
        as text " (" %5.1f 100*`n_efficient'/`n_dmus' "%)"
    if "`model'" == "additive" {
        di as text "Mean inefficiency:" as result %9.6f `mean_eff'
        di as text "Min inefficiency: " as result %9.6f `min_eff'
        di as text "Max inefficiency: " as result %9.6f `max_eff'
    }
    else {
        di as text "Mean efficiency: " as result %9.6f `mean_eff'
        di as text "Min efficiency:  " as result %9.6f `min_eff'
        di as text "Max efficiency:  " as result %9.6f `max_eff'
    }

    if `n_failed' > 0 {
        di as error _n "Warning: `n_failed' DMU(s) failed to converge"
    }
    di as text "{hline 50}"

    // ---- Tabular DMU-level results ----
    di as text _n "{bf:DMU Efficiency Scores}"
    di as text "{hline 50}"

    // Build a matrix-style result display
    tempname dearslt
    qui count if `touse'
    local nobs = r(N)

    // Rank DMUs by efficiency (descending: Rank 1 = most efficient)
    // First round to tolerance to handle numerical precision
    // Use track option: ties get same rank, next rank skips (1, 1, 3 style)
    tempvar eff_rounded rank_var
    qui gen double `eff_rounded' = round(`generate', `tolerance') if `touse'
    if "`model'" == "additive" {
        qui egen `rank_var' = rank(`eff_rounded') if `touse', track
    }
    else {
        qui egen `rank_var' = rank(-`eff_rounded') if `touse', track
    }

    // Show compact list
    preserve
    qui keep if `touse'
    di as text %12s "DMU" " {c |}" %12s "Score" %8s "Rank"
    di as text "{hline 13}{c +}{hline 20}"

    forvalues i = 1/`nobs' {
        local dmu_val = `dmu'[`i']
        local eff_val = `generate'[`i']
        local rank_val = `rank_var'[`i']
        if "`eff_val'" != "." {
            di as text %12s "`dmu_val'" " {c |}" ///
                as result %12.6f `eff_val' %8.0f `rank_val'
        }
        else {
            di as text %12s "`dmu_val'" " {c |}" ///
                as result %12s "FAILED" %8s "."
        }
    }
    di as text "{hline 34}"
    restore

    // ---- Reference set display (if requested) ----
    if "`reference'" != "" & "`lambda_vars'" != "" {
        di as text _n "{bf:Reference Sets} (lambda > `tolerance')"
        di as text "{hline 70}"

        preserve
        qui keep if `touse'
        local n_obs = _N

        forvalues k = 1/`n_obs' {
            local ref_list ""
            local ref_count = 0

            qui scalar _dmu_k = `dmu'[`k']
            local dmu_k = _dmu_k

            forvalues j = 1/`n_obs' {
                local lambda_var : word `j' of `lambda_vars'
                if "`lambda_var'" != "" {
                    qui scalar _lam_val = `lambda_var'[`k']
                    local lam_val = _lam_val
                    if `lam_val' > `tolerance' & `lam_val' < . {
                        qui scalar _dmu_j = `dmu'[`j']
                        local dmu_j = _dmu_j
                        local ref_count = `ref_count' + 1
                        if `ref_count' <= 6 {
                            local lam_str : display %6.4f `lam_val'
                            local ref_list "`ref_list' `dmu_j'(`lam_str')"
                        }
                    }
                }
            }
            if "`ref_list'" != "" {
                local suffix ""
                if `ref_count' > 6 {
                    local more = `ref_count' - 6
                    local suffix " ... (+`more' more)"
                }
                di as text "`dmu_k': " as result "`ref_list'" as text "`suffix'"
            }
        }
        restore
        di as text "{hline 70}"
    }

    // Return
    return scalar N               = `n_dmus'
    return scalar mean_efficiency = `mean_eff'
    return scalar min_eff         = `min_eff'
    return scalar max_eff         = `max_eff'
    return scalar n_efficient     = `n_efficient'
    return scalar n_failed        = `n_failed'
end


********************************************************************************
* 10. RESULT SAVING
********************************************************************************

program define _dea_save
    syntax , SAVing(string) Inputs(string) Outputs(string) ///
        DMU(varname) GENerate(name) TOUSE(varname) ///
        [REPLACE] [SAVERef] ///
        [weight_in_vars(string)] [weight_out_vars(string)] [weight_u0_var(string)]

    preserve
    qui keep if `touse'
    if "`saveref'" != "" {
        qui save "`saving'", `replace'
    }
    else {
        qui keep `dmu' `inputs' `outputs' `generate' `weight_in_vars' `weight_out_vars' `weight_u0_var'
        qui save "`saving'", `replace'
    }
    restore
    di as text _n "Results saved to: `saving'.dta"
end


********************************************************************************
* 11. SIMAR-WILSON BOOTSTRAP DISPATCHER
********************************************************************************

program define _dea_bootstrap, rclass
    syntax , Inputs(string) Outputs(string) DMU(varname) ///
        RTS(string) ORIentation(string) GENerate(name) ///
        TOUSE(varname) ENGINE(string) ///
        [ENVars(varlist numeric)] ///
        [REPS1(integer 100)] [REPS2(integer 2000)] ///
        [Level(integer 95)] [DEBUG]

    di as text _n "{hline 70}"
    di as text "Simar-Wilson Bootstrap"
    di as text "{hline 70}"
    
    // Bootstrap uses the efficiency scores already stored in `generate'.
    // The initial scores may be produced by either simplex or IPM.
    
    // Determine algorithm
    local algorithm = 1
    if "`envars'" != "" {
        local algorithm = 2
    }
    
    // Count observations
    qui count if `touse'
    local n_dmus = r(N)
    
    di as text "Algorithm:   " as result `algorithm' as text cond(`algorithm'==1, " (Bias correction only)", " (Two-stage with env. variables)")
    di as text "DMUs:        " as result `n_dmus'
    di as text "Inputs:      " as result "`inputs'"
    di as text "Outputs:     " as result "`outputs'"
    if "`envars'" != "" {
        di as text "Env. vars:   " as result "`envars'"
    }
    di as text "Score engine:" as result " `engine'"
    di as text "Bootstrap 1: " as result `reps1' as text " replications (bias correction)"
    if `algorithm' == 2 {
        di as text "Bootstrap 2: " as result `reps2' as text " replications (confidence intervals)"
    }
    di as text "Conf. level: " as result `level' "%"
    di as text "{hline 70}"
    
    // =====================================================================
    // Step 1: Initial DEA efficiency scores (already computed in `generate')
    // =====================================================================
    di as text _n "Step 1: Initial DEA efficiency scores"
    
    qui sum `generate' if `touse', detail
    local mean_eff = r(mean)
    local min_eff = r(min)
    local max_eff = r(max)
    di as text "  Mean: " as result %9.6f `mean_eff'
    di as text "  Min:  " as result %9.6f `min_eff'
    di as text "  Max:  " as result %9.6f `max_eff'
    
    // Store original efficiency scores
    tempvar eff_orig
    qui gen double `eff_orig' = `generate' if `touse'
    
    // Transform efficiency for bootstrap
    // Input-oriented: theta in (0,1], use 1/theta >= 1
    // Output-oriented: phi >= 1 already
    tempvar eff_trans
    if "`orientation'" == "input" {
        qui gen double `eff_trans' = 1 / `generate' if `touse' & `generate' > 0
    }
    else {
        qui gen double `eff_trans' = `generate' if `touse'
    }
    
    // =====================================================================
    // Algorithm branching
    // =====================================================================
    
    if `algorithm' == 1 {
        // =================================================================
        // ALGORITHM 1: Bias correction only (no environmental variables)
        // Uses smoothed bootstrap with reflection
        // =================================================================
        di as text _n "Step 2: Smoothed bootstrap bias correction (Algorithm 1)"
        
        // Compute bandwidth using Silverman's rule
        qui sum `eff_trans' if `touse'
        local sd_eff = r(sd)
        local h = 0.9 * `sd_eff' * `n_dmus'^(-0.2)
        di as text "  Bandwidth (Silverman): " as result %9.6f `h'
        
        // Bootstrap accumulator
        tempvar eff_sum eff_count
        qui gen double `eff_sum' = 0 if `touse'
        qui gen double `eff_count' = 0 if `touse'
        
        di as text "  Progress: " _continue
        forvalues b = 1/`reps1' {
            if mod(`b', 20) == 0 {
                di as text "." _continue
            }
            
            // Smoothed bootstrap: sample with replacement + kernel noise
            tempvar idx eff_star noise
            qui gen long `idx' = ceil(runiform() * `n_dmus') if `touse'
            qui gen double `noise' = rnormal(0, `h') if `touse'
            
            // Get resampled efficiency
            qui gen double `eff_star' = .
            forvalues i = 1/`n_dmus' {
                qui sum `idx' if _n == `i', meanonly
                local j = r(mean)
                if `j' >= 1 & `j' <= `n_dmus' {
                    qui sum `eff_trans' if _n == `j', meanonly
                    qui replace `eff_star' = r(mean) + `noise' if _n == `i' & `touse'
                }
            }
            
            // Reflection: if eff* < 1, reflect around 1
            qui replace `eff_star' = 2 - `eff_star' if `eff_star' < 1 & `touse'
            qui replace `eff_star' = max(`eff_star', 1.001) if `touse'
            
            // Convert back to original scale
            if "`orientation'" == "input" {
                qui replace `eff_star' = 1 / `eff_star' if `touse'
            }
            
            // Accumulate
            qui replace `eff_sum' = `eff_sum' + `eff_star' if `touse'
            qui replace `eff_count' = `eff_count' + 1 if `touse'
            
            drop `idx' `eff_star' `noise'
        }
        di as text " done"
        
        // Compute bias
        tempvar eff_boot_mean bias eff_bc
        qui gen double `eff_boot_mean' = `eff_sum' / `eff_count' if `touse'
        qui gen double `bias' = `eff_boot_mean' - `eff_orig' if `touse'
        qui gen double `eff_bc' = `eff_orig' - `bias' if `touse'
        
        // Bounds
        if "`orientation'" == "input" {
            qui replace `eff_bc' = min(`eff_bc', 1) if `touse'
            qui replace `eff_bc' = max(`eff_bc', 0.01) if `touse'
        }
        else {
            qui replace `eff_bc' = max(`eff_bc', 1) if `touse'
        }
    }
    else {
        // =================================================================
        // ALGORITHM 2: Two-stage with environmental variables
        // =================================================================
        di as text _n "Step 2: Truncated regression (Algorithm 2)"
        
        // Count environmental variables
        local n_envars : word count `envars'
        local min_obs = (`n_envars' + 1) * 3  // Rule of thumb: 3x parameters
        
        // Truncated regression: efficiency >= 1 (left-truncated at 1)
        local use_truncreg = 0
        capture which truncreg
        if _rc == 0 & `n_dmus' >= `min_obs' {
            local use_truncreg = 1
        }
        
        if `use_truncreg' {
            di as text "  Using truncated regression (lower limit = 1)"
            qui truncreg `eff_trans' `envars' if `touse', ll(1)
        }
        else {
            if `n_dmus' < `min_obs' {
                di as text "  Note: Too few observations (`n_dmus' < `min_obs') for truncreg, using OLS"
            }
            else {
                di as text "  Note: truncreg not available, using OLS approximation"
            }
            qui regress `eff_trans' `envars' if `touse'
        }
        
        // Store coefficients
        tempname beta
        matrix `beta' = e(b)
        local sigma_hat = e(rmse)
        if missing(`sigma_hat') local sigma_hat = e(sigma)
        if missing(`sigma_hat') {
            qui sum `eff_trans' if `touse'
            local sigma_hat = r(sd)
        }
        
        di as text "  Sigma (residual SD): " as result %9.6f `sigma_hat'
        di as text "  Coefficients:"
        matrix list `beta', noheader
        
        // Get predicted values and residuals
        tempvar xb resid
        qui predict double `xb' if `touse', xb
        qui gen double `resid' = `eff_trans' - `xb' if `touse'
        
        // Bootstrap bias correction
        di as text _n "Step 3: Bootstrap bias correction (B1 = `reps1' reps)"
        
        tempvar eff_sum eff_count
        qui gen double `eff_sum' = 0 if `touse'
        qui gen double `eff_count' = 0 if `touse'
        
        di as text "  Progress: " _continue
        forvalues b = 1/`reps1' {
            if mod(`b', 20) == 0 {
                di as text "." _continue
            }
            
            tempvar eps_star eff_star
            
            // Resample residuals
            qui gen double `eps_star' = .
            forvalues i = 1/`n_dmus' {
                local j = ceil(runiform() * `n_dmus')
                qui sum `resid' if _n == `j', meanonly
                qui replace `eps_star' = r(mean) if _n == `i' & `touse'
            }
            
            // Generate pseudo efficiencies
            qui gen double `eff_star' = `xb' + `eps_star' if `touse'
            qui replace `eff_star' = max(`eff_star', 1.001) if `touse'
            
            if "`orientation'" == "input" {
                qui replace `eff_star' = 1 / `eff_star' if `touse'
            }
            
            qui replace `eff_sum' = `eff_sum' + `eff_star' if `touse'
            qui replace `eff_count' = `eff_count' + 1 if `touse'
            
            drop `eps_star' `eff_star'
        }
        di as text " done"
        
        // Compute bias
        tempvar eff_boot_mean bias eff_bc
        qui gen double `eff_boot_mean' = `eff_sum' / `eff_count' if `touse'
        qui gen double `bias' = `eff_boot_mean' - `eff_orig' if `touse'
        qui gen double `eff_bc' = `eff_orig' - `bias' if `touse'
        
        // Bounds
        if "`orientation'" == "input" {
            qui replace `eff_bc' = min(`eff_bc', 1) if `touse'
            qui replace `eff_bc' = max(`eff_bc', 0.01) if `touse'
        }
        else {
            qui replace `eff_bc' = max(`eff_bc', 1) if `touse'
        }
        
        return matrix beta = `beta'
        return scalar sigma = `sigma_hat'
    }
    
    // =====================================================================
    // Store results
    // =====================================================================
    local bc_var = "`generate'_bc"
    local bias_var = "`generate'_bias"
    
    capture drop `bc_var'
    capture drop `bias_var'
    
    qui gen double `bc_var' = `eff_bc' if `touse'
    qui gen double `bias_var' = `bias' if `touse'
    
    // Summary
    di as text _n "{hline 70}"
    di as text "Bootstrap Results"
    di as text "{hline 70}"
    
    qui sum `generate' if `touse'
    di as text "Original efficiency:     Mean = " as result %9.6f r(mean) as text ", SD = " as result %9.6f r(sd)
    
    qui sum `bc_var' if `touse'
    di as text "Bias-corrected:          Mean = " as result %9.6f r(mean) as text ", SD = " as result %9.6f r(sd)
    
    qui sum `bias_var' if `touse'
    di as text "Bias:                    Mean = " as result %9.6f r(mean) as text ", SD = " as result %9.6f r(sd)
    
    di as text "{hline 70}"
    di as text "Variables created:"
    di as text "  `bc_var'    - Bias-corrected efficiency"
    di as text "  `bias_var'  - Bootstrap bias estimate"
    di as text "{hline 70}"
    
    // Return values
    return scalar N = `n_dmus'
    return scalar reps1 = `reps1'
    return scalar reps2 = `reps2'
    return scalar level = `level'
    return scalar algorithm = `algorithm'
    return local bc_var "`bc_var'"
    return local bias_var "`bias_var'"
    return add
end


********************************************************************************
* 12. IPM ENGINE WRAPPER
********************************************************************************

// Dispatches to Mata routines loaded from _dea_ipm.ado.

capture program drop _dea_ipm_solve
program define _dea_ipm_solve, rclass
    version 16.0

    syntax , Inputs(string) Outputs(string) DMU(varname) ///
        RTS(string) ORIentation(string) MODEL(string) ///
        STAGE(integer) TOUSE(varname) ///
        GENerate(name) TOLerance(real) ///
        [W(varlist numeric)] [P(varlist numeric)] ///
        [GX(varlist numeric)] [GY(varlist numeric)] ///
        [DDFBETA(name)] ///
        [BAD(varlist numeric)] [DIRection(string)] ///
        [REFerence] [SAVERef] [SAVEWeights] ///
        [lambda_vars(string)] ///
        [weight_in_vars(string)] [weight_out_vars(string)] [weight_u0_var(string)] ///
        [TRACE] [DEBUG] [SBMORT(integer 0)]

    // Prepare efficiency variable
    capture drop `generate'

    // Count observations
    qui count if `touse'
    local n_dmus = r(N)
    local m : word count `inputs'
    local s : word count `outputs'

    if "`debug'" != "" {
        di as text "[IPM] Solving: `model' `rts' `orientation'"
        di as text "[IPM] N=`n_dmus', m=`m', s=`s'"
    }

    // =====================================================================
    // Stage 1: Efficiency score computation
    // =====================================================================
    tempvar efficiency_score
    qui gen double `efficiency_score' = . if `touse'

    if "`model'" == "te" {
        // Technical Efficiency (Envelopment form with Dual Fallback)
        _dea_ipm_te, inputs(`inputs') outputs(`outputs') ///
            efficiency(`efficiency_score') touse(`touse') ///
            rts(`rts') orientation(`orientation') ///
            stage(`stage') tolerance(`tolerance') ///
            lambda_vars(`lambda_vars') ///
            `reference' `debug'
    }
    else if "`model'" == "cost" {
        _dea_ipm_cost, inputs(`inputs') outputs(`outputs') ///
            efficiency(`efficiency_score') touse(`touse') ///
            rts(`rts') w(`w') `debug'
    }
    else if "`model'" == "revenue" {
        _dea_ipm_revenue, inputs(`inputs') outputs(`outputs') ///
            efficiency(`efficiency_score') touse(`touse') ///
            rts(`rts') p(`p') `debug'
    }
    else if "`model'" == "profit" {
        _dea_ipm_profit, inputs(`inputs') outputs(`outputs') ///
            efficiency(`efficiency_score') touse(`touse') ///
            rts(`rts') w(`w') p(`p') `debug'
    }
    else if "`model'" == "ddf" {
        _dea_ipm_ddf, inputs(`inputs') outputs(`outputs') ///
            efficiency(`efficiency_score') touse(`touse') ///
            rts(`rts') gx(`gx') gy(`gy') `debug'
    }
    else if "`model'" == "undesirable" {
        _dea_ipm_ddf_undesirable, inputs(`inputs') outputs(`outputs') ///
            bad(`bad') efficiency(`efficiency_score') touse(`touse') ///
            rts(`rts') direction(`direction') `debug'
    }
    else if "`model'" == "sbm" {
        if `sbmort' == 0 {
            _dea_ipm_sbm, inputs(`inputs') outputs(`outputs') ///
                efficiency(`efficiency_score') touse(`touse') ///
                rts(`rts') `debug'
        }
        else {
            local _sbm_ort_str = cond(`sbmort'==1, "input", "output")
            mata: _dea_ipm_sbm_oriented("`inputs'", "`outputs'", ///
                "`efficiency_score'", "`touse'", "`rts'", "`_sbm_ort_str'")
        }
    }
    else if "`model'" == "additive" {
        _dea_ipm_additive, inputs(`inputs') outputs(`outputs') ///
            efficiency(`efficiency_score') touse(`touse') ///
            rts(`rts') `debug'
    }
    else if "`model'" == "supereff" {
        _dea_ipm_supereff, inputs(`inputs') outputs(`outputs') ///
            efficiency(`efficiency_score') touse(`touse') ///
            rts(`rts') orientation(`orientation') ///
            lambda_vars(`lambda_vars') ///
            `reference' `debug'
    }
    else if "`model'" == "multiplier" {
        _dea_ipm_multiplier, inputs(`inputs') outputs(`outputs') ///
            efficiency(`efficiency_score') touse(`touse') ///
            rts(`rts') orientation(`orientation') ///
            `saveweights' m(`m') s(`s') ///
            weight_in_vars(`weight_in_vars') ///
            weight_out_vars(`weight_out_vars') ///
            weight_u0_var(`weight_u0_var') ///
            `debug'
    }
    else if "`model'" == "allocative" {
        _dea_ipm_allocative, inputs(`inputs') outputs(`outputs') ///
            efficiency(`efficiency_score') touse(`touse') ///
            rts(`rts') orientation(`orientation') ///
            w(`w') p(`p') `debug'
    }
    else if "`model'" == "supersbm" {
        _dea_ipm_supersbm, inputs(`inputs') outputs(`outputs') ///
            efficiency(`efficiency_score') touse(`touse') ///
            rts(`rts') `debug'
    }
    else {
        di as error "model(`model') not yet implemented in IPM engine"
        exit 199
    }

    // =====================================================================
    // Stage 2: Slack Maximization (if stage == 2 and model == te)
    // =====================================================================
    // Stage 2 is handled inside _dea_ipm_te when stage==2

    // =====================================================================
    // Generate result variable
    // =====================================================================
    qui gen double `generate' = `efficiency_score' if `touse'
    label var `generate' "DEA efficiency (`model' `rts' `orientation')"

    return add
    return scalar n_dmus = `n_dmus'
end


********************************************************************************
* 13. IPM ENGINE — INDIVIDUAL MODEL WRAPPERS
********************************************************************************
* These wrappers keep the ado-level interface stable while delegating numerical
* work to Mata routines loaded from _dea_ipm.ado.

// ---------------------------------------------------------------------------
// 13.1 Technical Efficiency (Envelopment + Dual Fallback)
// ---------------------------------------------------------------------------
capture program drop _dea_ipm_te
program define _dea_ipm_te, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) ///
        EFFICIENCY(varname) TOUSE(varname) ///
        RTS(string) ORIentation(string) ///
        STAGE(integer) TOLerance(real) ///
        [lambda_vars(string)] ///
        [REFerence] [DEBUG]

    // ---------------------------------------------------------------
    // Uses unified envelopment solver from _dea_ipm.ado
    //   Supports: CRS, VRS, DRS, IRS
    //   Includes: Primal → Dual Fallback
    //   Stage 2:  Slack maximization (optional)
    // ---------------------------------------------------------------

    // Lambda variables are pre-created by _dea_model_dispatch
    // to ensure they persist through the display phase.
    qui count if `touse'
    local n_dmus = r(N)
    local m : word count `inputs'
    local s : word count `outputs'

    // Prepare slack variables if Stage 2
    local slack_in_vars ""
    local slack_out_vars ""
    local do_stage2 = 0
    if `stage' == 2 {
        local do_stage2 = 1
        local i = 0
        foreach v of local inputs {
            local ++i
            tempvar si`i'
            qui gen double `si`i'' = 0 if `touse'
            local slack_in_vars "`slack_in_vars' `si`i''"
        }
        local i = 0
        foreach v of local outputs {
            local ++i
            tempvar so`i'
            qui gen double `so`i'' = 0 if `touse'
            local slack_out_vars "`slack_out_vars' `so`i''"
        }
    }

    // Call unified envelopment solver
    // Use quietly to suppress LinearProgram warnings during Stage 2
    quietly mata: _dea_ipm_envelopment("`inputs'", "`outputs'", ///
        "`efficiency'", "`lambda_vars'", ///
        "`slack_in_vars'", "`slack_out_vars'", ///
        "`touse'", "`orientation'", "`rts'", ///
        `do_stage2', `tolerance')

    // Store slack variable names for downstream use
    if `do_stage2' {
        return local slack_in_vars "`slack_in_vars'"
        return local slack_out_vars "`slack_out_vars'"
    }
    return add
end

// ---------------------------------------------------------------------------
// 13.2-13.10 Other model wrappers
// ---------------------------------------------------------------------------

capture program drop _dea_ipm_cost
program define _dea_ipm_cost, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) ///
        EFFICIENCY(varname) TOUSE(varname) RTS(string) ///
        W(varlist numeric) [DEBUG]
    
    // Cost Efficiency: min w'x s.t. production constraints
    // Cost_Eff = min_cost / observed_cost
    
    if "`debug'" != "" {
        di as text "[Cost] Computing cost efficiency..."
        di as text "[Cost] RTS=`rts'"
        di as text "[Cost] Input prices: `w'"
    }
    
    mata: _dea_ipm_cost("`inputs'", "`outputs'", ///
        "`efficiency'", "`w'", "`touse'", "`rts'")
    
    return add
end

capture program drop _dea_ipm_revenue
program define _dea_ipm_revenue, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) ///
        EFFICIENCY(varname) TOUSE(varname) RTS(string) ///
        P(varlist numeric) [DEBUG]
    
    // Revenue Efficiency: max p'y s.t. production constraints
    // Rev_Eff = observed_revenue / max_revenue
    
    if "`debug'" != "" {
        di as text "[Revenue] Computing revenue efficiency..."
        di as text "[Revenue] RTS=`rts'"
        di as text "[Revenue] Output prices: `p'"
    }
    
    mata: _dea_ipm_revenue("`inputs'", "`outputs'", ///
        "`efficiency'", "`p'", "`touse'", "`rts'")
    
    return add
end

capture program drop _dea_ipm_profit
program define _dea_ipm_profit, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) ///
        EFFICIENCY(varname) TOUSE(varname) RTS(string) ///
        W(varlist numeric) P(varlist numeric) [DEBUG]
    
    // Profit Efficiency: max p'y - w'x s.t. production constraints
    // Profit_Eff = observed_profit / max_profit
    
    if "`debug'" != "" {
        di as text "[Profit] Computing profit efficiency..."
        di as text "[Profit] RTS=`rts'"
        di as text "[Profit] Input prices: `w'"
        di as text "[Profit] Output prices: `p'"
    }
    
    mata: _dea_ipm_profit("`inputs'", "`outputs'", ///
        "`efficiency'", "`w'", "`p'", "`touse'", "`rts'")
    
    return add
end

capture program drop _dea_ipm_ddf
program define _dea_ipm_ddf, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) ///
        EFFICIENCY(varname) TOUSE(varname) RTS(string) ///
        [GX(varlist numeric)] [GY(varlist numeric)] [DEBUG]
    
    // Directional Distance Function
    // max beta s.t. X*lambda <= x0 - beta*gx, Y*lambda >= y0 + beta*gy
    // Default direction: (gx, gy) = (x0, y0)
    
    if "`debug'" != "" {
        di as text "[DDF] Computing directional distance function..."
        di as text "[DDF] RTS=`rts'"
        if "`gx'" != "" di as text "[DDF] Direction gx: `gx'"
        if "`gy'" != "" di as text "[DDF] Direction gy: `gy'"
    }
    
    mata: _dea_ipm_ddf("`inputs'", "`outputs'", ///
        "`efficiency'", "`gx'", "`gy'", "`touse'", "`rts'")

    return add
end

capture program drop _dea_ipm_ddf_undesirable
program define _dea_ipm_ddf_undesirable, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) BAD(string) ///
        EFFICIENCY(varname) TOUSE(varname) RTS(string) ///
        [DIRection(string) DEBUG]

    // DDF with undesirable (bad) outputs — Wang et al. (2022)
    // max β s.t. X'λ+β·gx ≤ x₀, Y'λ-β·gy ≥ y₀, B'λ+β·gb = b₀ (equality!)
    // Direction: yb (default): gx=0, gy=y₀, gb=b₀
    //            xyb:          gx=x₀, gy=y₀, gb=b₀
    //            x:            gx=x₀, gy=0,  gb=0

    if "`direction'" == "" local direction "yb"

    if "`debug'" != "" {
        di as text "[Undesirable DDF] RTS=`rts', direction=`direction'"
        di as text "[Undesirable DDF] Bad outputs: `bad'"
    }

    mata: _dea_ipm_ddf_undesirable("`inputs'", "`outputs'", "`bad'", ///
        "`efficiency'", "`direction'", "`touse'", "`rts'")

    return add
end

capture program drop _dea_ipm_sbm
program define _dea_ipm_sbm, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) ///
        EFFICIENCY(varname) TOUSE(varname) RTS(string) [DEBUG]
    
    // Slacks-Based Measure (Tone, 2001)
    // Uses Charnes-Cooper transformation
    // SBM = [1 - (1/m)*sum(s-/x)] / [1 + (1/s)*sum(s+/y)]
    
    if "`debug'" != "" {
        di as text "[SBM] Computing slacks-based measure..."
        di as text "[SBM] RTS=`rts'"
    }
    
    // Create slack variables
    local m : word count `inputs'
    local s : word count `outputs'
    
    local slack_in_vars ""
    local slack_out_vars ""
    
    forvalues i = 1/`m' {
        tempvar si`i'
        qui gen double `si`i'' = . if `touse'
        local slack_in_vars "`slack_in_vars' `si`i''"
    }
    forvalues r = 1/`s' {
        tempvar so`r'
        qui gen double `so`r'' = . if `touse'
        local slack_out_vars "`slack_out_vars' `so`r''"
    }
    
    mata: _dea_ipm_sbm("`inputs'", "`outputs'", ///
        "`efficiency'", "`slack_in_vars'", "`slack_out_vars'", ///
        "`touse'", "`rts'")
    
    return local slack_in_vars "`slack_in_vars'"
    return local slack_out_vars "`slack_out_vars'"
    return add
end

capture program drop _dea_ipm_additive
program define _dea_ipm_additive, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) ///
        EFFICIENCY(varname) TOUSE(varname) RTS(string) [DEBUG]
    
    // Additive Model: max sum(s-) + sum(s+)
    // Efficiency = 0 if efficient, > 0 if inefficient
    
    if "`debug'" != "" {
        di as text "[Additive] Computing additive model..."
        di as text "[Additive] RTS=`rts'"
    }
    
    // Create slack variables
    local m : word count `inputs'
    local s : word count `outputs'
    
    local slack_in_vars ""
    local slack_out_vars ""
    
    forvalues i = 1/`m' {
        tempvar si`i'
        qui gen double `si`i'' = . if `touse'
        local slack_in_vars "`slack_in_vars' `si`i''"
    }
    forvalues r = 1/`s' {
        tempvar so`r'
        qui gen double `so`r'' = . if `touse'
        local slack_out_vars "`slack_out_vars' `so`r''"
    }
    
    mata: _dea_ipm_additive("`inputs'", "`outputs'", ///
        "`efficiency'", "`slack_in_vars'", "`slack_out_vars'", ///
        "`touse'", "`rts'")
    
    return local slack_in_vars "`slack_in_vars'"
    return local slack_out_vars "`slack_out_vars'"
    return add
end

capture program drop _dea_ipm_supereff
program define _dea_ipm_supereff, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) ///
        EFFICIENCY(varname) TOUSE(varname) RTS(string) ///
        ORIentation(string) [lambda_vars(string)] [REFerence] [DEBUG]
    
    // Super-efficiency: exclude evaluated DMU from reference set
    
    if "`debug'" != "" {
        di as text "[SuperEff] Computing super-efficiency..."
        di as text "[SuperEff] RTS=`rts', Orientation=`orientation'"
    }
    
    // Lambda variables are pre-created by _dea_model_dispatch
    qui count if `touse'
    local n_dmus = r(N)
    
    mata: _dea_ipm_supereff("`inputs'", "`outputs'", ///
        "`efficiency'", "`lambda_vars'", ///
        "`touse'", "`orientation'", "`rts'")
    
    return add
end

capture program drop _dea_ipm_multiplier
program define _dea_ipm_multiplier, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) ///
        EFFICIENCY(varname) TOUSE(varname) RTS(string) ///
        ORIentation(string) M(integer) S(integer) ///
        [SAVEWeights] [weight_in_vars(string)] [weight_out_vars(string)] ///
        [weight_u0_var(string)] [DEBUG]
    
    // Multiplier (Dual) Model
    // Input-oriented: max u'y / v'x
    // Output-oriented: min v'x / u'y
    
    if "`debug'" != "" {
        di as text "[Multiplier] Computing multiplier form..."
        di as text "[Multiplier] RTS=`rts', Orientation=`orientation'"
    }
    
    // Create weight variables if requested
    local v_vars ""
    local u_vars ""
    
    if "`saveweights'" != "" {
        if "`weight_in_vars'" != "" & "`weight_out_vars'" != "" {
            local v_vars "`weight_in_vars'"
            local u_vars "`weight_out_vars'"
        }
        else {
            forvalues i = 1/`m' {
                tempvar v`i'
                qui gen double `v`i'' = . if `touse'
                local v_vars "`v_vars' `v`i''"
            }
            forvalues r = 1/`s' {
                tempvar u`r'
                qui gen double `u`r'' = . if `touse'
                local u_vars "`u_vars' `u`r''"
            }
        }
    }
    
    if "`saveweights'" != "" {
        mata: _dea_ipm_multiplier("`inputs'", "`outputs'", ///
            "`efficiency'", "`v_vars'", "`u_vars'", "`weight_u0_var'", ///
            "`touse'", "`orientation'", "`rts'")
    }

    // By strong duality, multiplier efficiency equals the radial
    // envelopment score. Use the envelopment solver for the reported score;
    // multiplier weights, when requested, are retained from the dual LP above.
    mata: _dea_ipm_envelopment("`inputs'", "`outputs'", ///
        "`efficiency'", "", "", "", "`touse'", "`orientation'", "`rts'", 0, 1e-8)
    
    if "`saveweights'" != "" {
        return local v_vars "`v_vars'"
        return local u_vars "`u_vars'"
    }
    return add
end

capture program drop _dea_ipm_allocative
program define _dea_ipm_allocative, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) ///
        EFFICIENCY(varname) TOUSE(varname) RTS(string) ///
        ORIentation(string) ///
        [W(varlist numeric)] [P(varlist numeric)] [DEBUG]
    
    // Allocative Efficiency = Economic Efficiency / Technical Efficiency
    // For cost (input-oriented): AE = CE / TE
    // For revenue (output-oriented): AE = RE / TE
    
    if "`debug'" != "" {
        di as text "[Allocative] Computing allocative efficiency..."
        di as text "[Allocative] RTS=`rts', Orientation=`orientation'"
    }
    
    // First compute technical efficiency
    tempvar te_score
    qui gen double `te_score' = . if `touse'
    
    mata: _dea_ipm_envelopment("`inputs'", "`outputs'", ///
        "`te_score'", "", "", "", "`touse'", "`orientation'", "`rts'", 0, 1e-8)
    
    // Then compute economic efficiency
    tempvar econ_eff
    qui gen double `econ_eff' = . if `touse'
    
    if "`orientation'" == "input" & "`w'" != "" {
        // Cost efficiency
        mata: _dea_ipm_cost("`inputs'", "`outputs'", ///
            "`econ_eff'", "`w'", "`touse'", "`rts'")
    }
    else if "`orientation'" == "output" & "`p'" != "" {
        // Revenue efficiency
        mata: _dea_ipm_revenue("`inputs'", "`outputs'", ///
            "`econ_eff'", "`p'", "`touse'", "`rts'")
    }
    else {
        di as error "allocative model requires w() for input orientation or p() for output orientation"
        exit 198
    }
    
    // Allocative efficiency = Economic / Technical
    qui replace `efficiency' = `econ_eff' / `te_score' if `touse' & `te_score' > 0
    
    return add
end

capture program drop _dea_ipm_supersbm
program define _dea_ipm_supersbm, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) ///
        EFFICIENCY(varname) TOUSE(varname) RTS(string) [DEBUG]
    
    // Super-efficiency SBM (Tone, 2002)
    // Excludes evaluated DMU from reference set
    
    if "`debug'" != "" {
        di as text "[SuperSBM] Computing super-efficiency SBM..."
        di as text "[SuperSBM] RTS=`rts'"
    }
    
    mata: _dea_ipm_supersbm("`inputs'", "`outputs'", ///
        "`efficiency'", "`touse'", "`rts'")
    
    return add
end


********************************************************************************
* 14. SIMPLEX ENGINE DISPATCHER (routes to _dea_simplex.ado Mata functions)
********************************************************************************

capture program drop _dea_simplex_solve
program define _dea_simplex_solve, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) DMU(varname) ///
        RTS(string) ORIentation(string) MODEL(string) ///
        STAGE(integer) TOUSE(varname) GENerate(name) TOLerance(real) ///
        [METHod(string)] ///
        [W(varlist numeric)] [P(varlist numeric)] ///
        [GX(varlist numeric)] [GY(varlist numeric)] ///
        [DDFBETA(name)] ///
        [BAD(varlist numeric)] [DIRection(string)] ///
        [REFerence] [SAVERef] [SAVEWeights] ///
        [lambda_vars(string)] ///
        [slack_in_vars(string)] [slack_out_vars(string)] ///
        [weight_in_vars(string)] [weight_out_vars(string)] [weight_u0_var(string)] ///
        [SBMORT(integer 0)] ///
        [TRACE] [DEBUG]

    // =================================================================
    // Load Simplex Mata functions
    // =================================================================
    // Check both core and undesirable functions; reload if either is missing
    local _simplex_needs_load 0
    capture mata: mata which _dea_simplex_envelopment()
    if _rc != 0 local _simplex_needs_load 1
    capture mata: mata which _dea_simplex_ddf_undesirable()
    if _rc != 0 local _simplex_needs_load 1
    if `_simplex_needs_load' {
        capture findfile _dea_simplex.ado
        if _rc == 0 {
            qui run "`r(fn)'"
        }
        else {
            di as error "_dea_simplex.ado not found. Please install it in your ado path."
            exit 601
        }
    }
    
    // =================================================================
    // Create efficiency variable
    // =================================================================
    capture drop `generate'
    qui gen double `generate' = . if `touse'
    
    if "`debug'" != "" {
        di as text "[Simplex] Model: `model'  RTS: `rts'  ORT: `orientation'"
    }
    
    // =================================================================
    // Dispatch to appropriate Simplex Mata function
    // =================================================================
    
    if "`model'" == "te" {
        // ---------------------------------------------------------
        // Technical efficiency (envelopment)
        // slack_in_vars / slack_out_vars pre-created in _dea_model_dispatch
        // ---------------------------------------------------------
        quietly mata: _dea_simplex_envelopment("`inputs'", "`outputs'", ///
            "`generate'", "`touse'", "`orientation'", "`rts'", ///
            `tolerance', "`lambda_vars'", "`slack_in_vars'", "`slack_out_vars'")
    }
    else if "`model'" == "multiplier" {
        // ---------------------------------------------------------
        // Multiplier (dual) form. Reported efficiency equals the
        // radial envelopment score by strong duality.  Stage-2
        // saveweights stores one optimal multiplier solution.
        // ---------------------------------------------------------
        if "`saveweights'" != "" & `stage' == 2 & ///
            "`weight_in_vars'" != "" & "`weight_out_vars'" != "" {
            quietly mata: _dea_simplex_multiplier("`inputs'", "`outputs'", ///
                "`generate'", "`weight_in_vars'", "`weight_out_vars'", "`weight_u0_var'", ///
                "`touse'", "`orientation'", "`rts'", `tolerance', "`trace'")
        }
        else {
            quietly mata: _dea_simplex_envelopment("`inputs'", "`outputs'", ///
                "`generate'", "`touse'", "`orientation'", "`rts'", `tolerance', "", "", "")
        }
    }
    else if "`model'" == "cost" {
        // ---------------------------------------------------------
        // Cost efficiency
        // ---------------------------------------------------------
        quietly mata: _dea_simplex_cost("`inputs'", "`outputs'", ///
            "`generate'", "`w'", "`touse'", "`rts'")
    }
    else if "`model'" == "revenue" {
        // ---------------------------------------------------------
        // Revenue efficiency
        // ---------------------------------------------------------
        quietly mata: _dea_simplex_revenue("`inputs'", "`outputs'", ///
            "`generate'", "`p'", "`touse'", "`rts'")
    }
    else if "`model'" == "profit" {
        // ---------------------------------------------------------
        // Profit efficiency
        // ---------------------------------------------------------
        quietly mata: _dea_simplex_profit("`inputs'", "`outputs'", ///
            "`generate'", "`w'", "`p'", "`touse'", "`rts'")
    }
    else if "`model'" == "ddf" {
        // ---------------------------------------------------------
        // Directional distance function
        // ---------------------------------------------------------
        local gx_arg = cond("`gx'" != "", "`gx'", "")
        local gy_arg = cond("`gy'" != "", "`gy'", "")
        quietly mata: _dea_simplex_ddf("`inputs'", "`outputs'", ///
            "`generate'", "`gx_arg'", "`gy_arg'", "`touse'", "`rts'")
    }
    else if "`model'" == "undesirable" {
        // ---------------------------------------------------------
        // DDF with undesirable (bad) outputs — Wang et al. (2022)
        // ---------------------------------------------------------
        quietly mata: _dea_simplex_ddf_undesirable("`inputs'", "`outputs'", ///
            "`bad'", "`generate'", "`direction'", "`touse'", "`rts'")
    }
    else if "`model'" == "sbm" {
        // ---------------------------------------------------------
        // Slacks-based measure
        // sbm_ort: 0=non-oriented, 1=input-oriented, 2=output-oriented
        // ---------------------------------------------------------
        local _sbm_ort_str = cond(`sbmort'==1, "input", cond(`sbmort'==2, "output", ""))
        quietly mata: _dea_simplex_sbm("`inputs'", "`outputs'", ///
            "`generate'", "", "", "`touse'", "`rts'", "`_sbm_ort_str'")
    }
    else if "`model'" == "additive" {
        // ---------------------------------------------------------
        // Additive model
        // ---------------------------------------------------------
        quietly mata: _dea_simplex_additive("`inputs'", "`outputs'", ///
            "`generate'", "", "", "`touse'", "`rts'")
    }
    else if "`model'" == "supereff" {
        // ---------------------------------------------------------
        // Radial super-efficiency (with optional lambda extraction)
        // Lambda variables are pre-created by _dea_model_dispatch
        // ---------------------------------------------------------
        quietly mata: _dea_simplex_supereff("`inputs'", "`outputs'", ///
            "`generate'", "`lambda_vars'", ///
            "`touse'", "`orientation'", "`rts'")
    }
    else if "`model'" == "supersbm" {
        // ---------------------------------------------------------
        // Super-efficiency SBM
        // ---------------------------------------------------------
        quietly mata: _dea_simplex_supersbm("`inputs'", "`outputs'", ///
            "`generate'", "`touse'", "`rts'")
    }
    else if "`model'" == "allocative" {
        // ---------------------------------------------------------
        // Allocative efficiency (composite: TE + Cost/Revenue)
        // ---------------------------------------------------------
        local w_arg = cond("`w'" != "", "`w'", "")
        local p_arg = cond("`p'" != "", "`p'", "")
        quietly mata: _dea_simplex_allocative("`inputs'", "`outputs'", ///
            "`generate'", "`touse'", "`orientation'", "`rts'", ///
            "`w_arg'", "`p_arg'")
    }
    else {
        di as error "model(`model') is not supported by the Simplex engine"
        exit 198
    }

    // =================================================================
    // Count results
    // =================================================================
    qui count if `touse'
    local n_dmus = r(N)
    
    qui count if `generate' < . & `touse'
    local n_solved = r(N)
    
    if "`debug'" != "" {
        di as text "[Simplex] Solved `n_solved' of `n_dmus' DMUs"
    }
    
    return scalar N = `n_dmus'
    return scalar n_solved = `n_solved'
    return add
end


********************************************************************************
* 15. CONGESTION MODEL DISPATCHER
*     Cooper, W.W., Deng, H., Seiford, L.M., and Zhu, J. (2004).
*     "Congestion: Its identification and management with DEA."
*     In Handbook on Data Envelopment Analysis, pp. 177-201.
*     Boston: Kluwer Academic Publishers.
********************************************************************************
//
// Two-step procedure:
//   Step 1: Standard VRS DEA → θ_VRS (strong disposability)
//   Step 2: Congestion model → θ_C  (output equality constraint)
//   Congestion measure = θ_VRS / θ_C  (< 1 indicates congestion)
//

capture program drop _dea_congestion_dispatch
program define _dea_congestion_dispatch, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) DMU(varname) ///
        ENGINE(string) TOUSE(varname) ///
        GENerate(name) TOLerance(real) ///
        [DEBUG]

    // The congestion procedure uses IPM/LinearProgram-based routines internally.
    // The engine() option is accepted by the unified interface but is not used
    // to switch this specialized dispatcher to the Simplex engine.

    // ---------------------------------------------------------------
    // Display header
    // ---------------------------------------------------------------
    di as text "{hline 65}"
    di as text "Congestion DEA — Cooper, Deng, Seiford & Zhu (2004)"
    di as text "{hline 65}"

    // ---------------------------------------------------------------
    // Step 1: Standard VRS and CRS efficiencies
    // ---------------------------------------------------------------
    tempvar _eff_vrs _eff_crs
    quietly gen double `_eff_vrs' = . if `touse'
    quietly gen double `_eff_crs' = . if `touse'

    // VRS efficiency (θ_VRS)
    _dea_ipm_solve, inputs(`inputs') outputs(`outputs') ///
        dmu(`dmu') rts(vrs) orientation(input) ///
        model(te) stage(1) touse(`touse') ///
        generate(`_eff_vrs') tolerance(`tolerance') `debug'

    // CRS efficiency (for scale efficiency)
    _dea_ipm_solve, inputs(`inputs') outputs(`outputs') ///
        dmu(`dmu') rts(crs) orientation(input) ///
        model(te) stage(1) touse(`touse') ///
        generate(`_eff_crs') tolerance(`tolerance') `debug'

    // ---------------------------------------------------------------
    // Step 2: Congestion LP (weak input disposability → θ_C)
    //   min θ  s.t.  Xλ = θ·x₀,  Yλ = y₀,  Σλ = 1,  λ ≥ 0
    // ---------------------------------------------------------------
    tempvar _eff_cong
    quietly gen double `_eff_cong' = . if `touse'

    mata: _dea_congestion_lp("`inputs'", "`outputs'", "`touse'", ///
        `tolerance', "`_eff_cong'")

    // ---------------------------------------------------------------
    // Create output variables
    // ---------------------------------------------------------------
    local v_crs    "`generate'_crs"
    local v_vrs    "`generate'_vrs"
    local v_ceff   "`generate'_cong_eff"
    local v_cong   "`generate'_congestion"
    local v_se     "`generate'_se"

    foreach v in `v_crs' `v_vrs' `v_ceff' `v_cong' `v_se' {
        capture drop `v'
        quietly gen double `v' = . if `touse'
    }

    quietly replace `v_crs'  = `_eff_crs'  if `touse'
    quietly replace `v_vrs'  = `_eff_vrs'  if `touse'
    quietly replace `v_ceff' = `_eff_cong' if `touse'

    // Congestion measure = θ_VRS / θ_C  (1 = no congestion)
    quietly replace `v_cong' = `v_vrs' / `v_ceff' ///
        if `touse' & `v_ceff' > `tolerance'
    quietly replace `v_cong' = 1 if `touse' & `v_ceff' <= `tolerance'

    // Scale efficiency = θ_CRS / θ_VRS
    quietly replace `v_se' = `v_crs' / `v_vrs' ///
        if `touse' & `v_vrs' > `tolerance'

    // Primary generate = congestion measure (1 = efficient, <1 = congested)
    capture drop `generate'
    quietly gen double `generate' = `v_cong' if `touse'
    label var `generate' "Congestion measure (θ_VRS/θ_C; <1=congested)"

    // ---------------------------------------------------------------
    // Summary statistics
    // ---------------------------------------------------------------
    qui count if `touse'
    local N = r(N)
    qui summ `v_vrs' if `touse'
    local mean_vrs = r(mean)
    qui count if `v_cong' < (1 - `tolerance') & `touse'
    local n_cong = r(N)

    di as text "DMUs evaluated:  " as result `N'
    di as text "DMUs congested:  " as result `n_cong' ///
        as text " (congestion measure < 1)"
    di as text "Mean VRS eff.:   " as result %8.6f `mean_vrs'
    di as text "{hline 65}"
    di as text "Output variables (prefix: `generate'):"
    di as text "  `v_crs'    CRS technical efficiency"
    di as text "  `v_vrs'    VRS technical efficiency"
    di as text "  `v_ceff'  Congestion-model efficiency (θ_C)"
    di as text "  `v_cong'  Congestion measure (θ_VRS/θ_C)"
    di as text "  `v_se'     Scale efficiency (θ_CRS/θ_VRS)"

    // Return scalars
    return scalar N           = `N'
    return scalar n_congested = `n_cong'
    return scalar mean_vrs    = `mean_vrs'
    // Standard dea returns
    qui summ `v_cong' if `touse'
    return scalar mean_efficiency = r(mean)
    return scalar n_efficient     = `N' - `n_cong'
    return scalar min_eff         = r(min)
    return scalar max_eff         = r(max)
    return scalar n_failed        = 0
end


********************************************************************************
* 16. SBM-UNDESIRABLE DISPATCHER
*     Tone, K. (2004). Dealing with undesirable outputs in DEA: A Slacks-Based
*     Measure (SBM) approach. GRIPS Discussion Paper I-2003-0005.
*     doi:10.24545/00000955
********************************************************************************
//
// SBM with undesirable (bad) outputs via Charnes-Cooper linearization.
//
// Fractional program:
//   min  (1 - (1/m)*Σ(s⁻ᵢ/xᵢ₀)) / (1 + (1/(s+b))*(Σ(s⁺ᵣ/yᵣ₀)+Σ(sᵇₖ/bₖ₀)))
//   s.t. Xλ + s⁻ = x₀
//        Yλ - s⁺ = y₀  (good outputs: free disposal, ≥)
//        Bλ + sᵇ = b₀  (bad outputs: weak disposability, =)
//        Σλ = 1 (VRS) or unconstrained (CRS)
//        λ, s⁻, s⁺, sᵇ ≥ 0
//

capture program drop _dea_sbmundesirable_dispatch
program define _dea_sbmundesirable_dispatch, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) BAD(string) ///
        DMU(varname) RTS(string) TOUSE(varname) ///
        GENerate(name) TOLerance(real) ///
        [DEBUG]

    // ---------------------------------------------------------------
    // Display header
    // ---------------------------------------------------------------
    local nbad : word count `bad'
    di as text "{hline 65}"
    di as text "SBM-Undesirable DEA — Tone (2004)"
    di as text "{hline 65}"
    di as text "RTS: `rts'    Bad outputs: `nbad'"

    // Create efficiency variable
    capture drop `generate'
    quietly gen double `generate' = . if `touse'

    // Call Mata LP
    mata: _dea_sbmu_solve("`inputs'", "`outputs'", "`bad'", ///
        "`touse'", "`rts'", `tolerance', "`generate'")

    label var `generate' "SBM-Undesirable efficiency (Tone 2004)"

    // Summary
    qui count if `touse'
    local N = r(N)
    qui summ `generate' if `touse'
    local mean_eff = r(mean)
    qui count if `generate' >= (1 - `tolerance') & `touse'
    local n_eff = r(N)

    di as text "DMUs evaluated:  " as result `N'
    di as text "Efficient DMUs:  " as result `n_eff'
    di as text "Mean efficiency: " as result %8.6f `mean_eff'
    di as text "{hline 65}"

    // Standard dea returns
    return scalar mean_efficiency = `mean_eff'
    return scalar n_efficient     = `n_eff'
    return scalar N               = `N'
    qui summ `generate' if `touse'
    return scalar min_eff  = r(min)
    return scalar max_eff  = r(max)
    return scalar n_failed = 0
end


********************************************************************************
* 17. BACKWARD COMPATIBILITY ALIAS
********************************************************************************

// dea_lp → dea (deprecated alias)
// Users who type `dea_lp' get redirected to `dea' with a note

capture program drop dea_lp
program define dea_lp, rclass
    version 16.0
    di as text "Note: dea_lp is now integrated into dea v2.0"
    di as text "      Redirecting to: dea `0'"
    dea `0'
    return add
end


********************************************************************************
* MATA: CONGESTION LP SOLVER
*   min θ  s.t.  Xλ = θ·x₀,  Yλ = y₀,  Σλ = 1,  λ ≥ 0
********************************************************************************
mata:
void _dea_congestion_lp(string scalar inp_names,
                         string scalar out_names,
                         string scalar touse,
                         real scalar   tol,
                         string scalar outvar)
{
    real matrix X, Y
    real scalar n, m, s, k
    real colvector scores
    real rowvector x0, y0, c, lb, ub
    real matrix A_eq
    real colvector b_eq
    class LinearProgram scalar S
    real colvector sol
    real colvector obs_idx

    X = st_data(., tokens(inp_names), touse)
    Y = st_data(., tokens(out_names), touse)
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    obs_idx = selectindex(st_data(., touse))

    scores = J(n, 1, .)

    for (k = 1; k <= n; k++) {
        x0 = X[k, .]
        y0 = Y[k, .]

        // Variables: [λ₁...λₙ, θ]  — n+1 variables
        // Objective: min θ
        c  = (J(1, n, 0), 1)

        // Weak input disposability: Xλ = θ·x₀
        A_eq   = (X', -x0')
        b_eq   = J(m, 1, 0)

        // Output EQUALITY: Yλ = y₀  (congestion constraint)
        A_eq   = A_eq \ (Y', J(s, 1, 0))
        b_eq   = b_eq \ y0'

        // VRS convexity: Σλ = 1
        A_eq   = A_eq \ (J(1, n, 1), 0)
        b_eq   = b_eq \ 1

        // Bounds: λ ≥ 0, θ ∈ [0,1]
        lb = (J(1, n, 0), 0)
        ub = (J(1, n, .), 1)

        S = LinearProgram()
        S.setMaxOrMin("min")
        S.setCoefficients(c)
        S.setEquality(A_eq, b_eq)
        S.setBounds(lb, ub)
        S.setTrace("off")
        (void) S.optimize()

        if (S.converged() & S.errorcode() == 0) {
            sol = S.parameters()
            scores[k] = sol[n + 1]
        }
        else {
            scores[k] = 1   // infeasible → no congestion
        }
    }

    // Store in Stata variable
    real scalar vi
    vi = st_varindex(outvar)
    for (k = 1; k <= n; k++) {
        st_store(obs_idx[k], vi, scores[k])
    }
}
end


********************************************************************************
* MATA: SBM-UNDESIRABLE LP SOLVER
*   Tone (2004) — Charnes-Cooper linearized form
********************************************************************************
mata:
void _dea_sbmu_solve(string scalar inp_names,
                      string scalar out_names,
                      string scalar bad_names,
                      string scalar touse,
                      string scalar rts,
                      real scalar   tol,
                      string scalar outvar)
{
    real matrix X, Y, B
    real scalar n, m, s, b, k
    real colvector scores, obs_idx
    real rowvector x0, y0, b0
    real rowvector c, lb, ub
    real matrix A_ineq, A_eq
    real colvector b_ineq_v, b_eq_v
    class LinearProgram scalar S
    real colvector sol
    real scalar theta, t, vi
    real scalar nv, i, r
    real matrix Aeq2, Aeq3
    real rowvector Aeq4, Aeq5

    X = st_data(., tokens(inp_names), touse)
    Y = st_data(., tokens(out_names), touse)
    B = st_data(., tokens(bad_names), touse)
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    b = cols(B)
    obs_idx = selectindex(st_data(., touse))

    scores = J(n, 1, .)

    // Charnes-Cooper transformation variables (order):
    //   [z₁...zₙ, S⁻₁...S⁻ₘ, S⁺₁...S⁺ₛ, Sᵇ₁...Sᵇᵦ, t]
    //   total: n + m + s + b + 1 variables

    nv = n + m + s + b + 1

    for (k = 1; k <= n; k++) {
        x0 = X[k, .]
        y0 = Y[k, .]
        b0 = B[k, .]

        // Objective: min t - (1/m)*Σ(S⁻ᵢ/xᵢ₀)
        c = J(1, nv, 0)
        c[nv] = 1   // t coefficient
        for (i = 1; i <= m; i++) {
            if (x0[i] > tol) c[n + i] = -(1/m) / x0[i]
        }

        // Equality 1: Xz + S⁻ = t·x₀  →  X·z + I·S⁻ + 0·S⁺ + 0·Sᵇ - x₀·t = 0
        //   X[,k] is n×1 column; z is n×1
        A_eq   = J(m, nv, 0)
        b_eq_v = J(m, 1, 0)
        A_eq[1..m, 1..n]     = X'
        A_eq[1..m, n+1..n+m] = I(m)
        A_eq[1..m, nv]       = -x0'

        // Equality 2: Yz - S⁺ = t·y₀
        Aeq2 = J(s, nv, 0)
        Aeq2[1..s, 1..n]       = Y'
        Aeq2[1..s, n+m+1..n+m+s] = -I(s)
        Aeq2[1..s, nv]         = -y0'
        A_eq   = A_eq \ Aeq2
        b_eq_v = b_eq_v \ J(s, 1, 0)

        // Equality 3: Bz + Sᵇ = t·b₀  (bad output: weak disposability)
        Aeq3 = J(b, nv, 0)
        Aeq3[1..b, 1..n]             = B'
        Aeq3[1..b, n+m+s+1..n+m+s+b] = I(b)
        Aeq3[1..b, nv]               = -b0'
        A_eq   = A_eq \ Aeq3
        b_eq_v = b_eq_v \ J(b, 1, 0)

        // Equality 4: normalizer
        //   t + (1/(s+b))*(Σ(S⁺ᵣ/yᵣ₀) + Σ(Sᵇₖ/bₖ₀)) = 1
        Aeq4 = J(1, nv, 0)
        Aeq4[nv] = 1
        for (r = 1; r <= s; r++) {
            if (y0[r] > tol) Aeq4[n + m + r] = (1/(s+b)) / y0[r]
        }
        for (r = 1; r <= b; r++) {
            if (b0[r] > tol) Aeq4[n + m + s + r] = (1/(s+b)) / b0[r]
        }
        A_eq   = A_eq \ Aeq4
        b_eq_v = b_eq_v \ 1

        // VRS: Σzⱼ = t
        if (rts == "vrs") {
            Aeq5 = J(1, nv, 0)
            Aeq5[1..n] = J(1, n, 1)
            Aeq5[nv]   = -1
            A_eq   = A_eq \ Aeq5
            b_eq_v = b_eq_v \ 0
        }

        // No inequality constraints needed (all equality)
        A_ineq  = J(0, nv, 0)
        b_ineq_v = J(0, 1, 0)

        // Bounds: all ≥ 0, t ≥ 0
        lb = J(1, nv, 0)
        ub = J(1, nv, .)

        S = LinearProgram()
        S.setMaxOrMin("min")
        S.setCoefficients(c)
        if (rows(A_ineq) > 0) {
            S.setInequality(A_ineq, b_ineq_v)
        }
        S.setEquality(A_eq, b_eq_v)
        S.setBounds(lb, ub)
        S.setTrace("off")
        (void) S.optimize()

        if (S.converged() & S.errorcode() == 0) {
            sol = S.parameters()
            theta = S.value()
            scores[k] = theta
        }
        else {
            scores[k] = .
        }
    }

    // Store in Stata variable
    vi = st_varindex(outvar)
    for (k = 1; k <= n; k++) {
        st_store(obs_idx[k], vi, scores[k])
    }
}
end
