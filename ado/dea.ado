*! version 2.0.73 06nov2026
*! Check successful simplex bases; retain certified results and refine verified defects.
*! Validate ordinary TE quantities before replacing results or saving.
*! Optional recoverylog: LP context, model/final availability and command lifecycle.
*! R11-F: specialized model domains, SBM-bad RTS, scaled equality index and explicit conventions.
*! Correct SBM zero terms and scaling; connect super-SBM orientation and validate data.
*! Validate allocative price dimensions and orientation before changing results.
*! Correct interval IDEA self references in both pessimistic and optimistic LPs.
*! Forward the requested tolerance to every legacy bounded/ordinal DEA call.
*! Forward non-TE tolerance and validate original LP constraints before returning results.
*! Return TE numerical diagnostics and preserve small positive radial scores.
*! Preserve specialized-model failures and count unavailable components.
*! Uniform TE stage status, solution provenance, and save-only diagnostics.
*! Return TE slack columns with the expanded input/output variable names and order.
*! Display DMU labels as data, including quotes, macro characters, and literal SMCL.
*! Explicit radialtolerance() recovery preserves scores and reports projection targets.
*! 2026-11-01: connect bounded native retries and per-DMU evidence to TE/VRS.
*! 2026-11-01: validate primal residuals, deliver TE tolerance, and correct output dual recovery.
*! 2026-11-01: return IPM stage-2 projections and separate radial/projection status.
*! 2026-11-01: move dea_lp compatibility entry point to its standalone file.
*! 2026-11-01: require a positive, finite effective tolerance before sample or output changes.
*! 2026-11-01: keep reference display values local; preserve user scalars and matrices.
*! 2026-11-01: exclude missing bad/direction values; preflight numeric IDEA endpoints.
*! 2026-11-01: pass stage to simplex; use radial-only auxiliary VRS fits.
*! 2026-11-01: connect numeric bounds plus legacy options to the shared DMU-wise transformation.
*! 2026-11-01: freeze ordinal comparison values before transforming each DMU dataset.
*! 2026-11-01: restore bounded string inputs/outputs in both syntaxes; preflight endpoints.
*! 2026-11-01: unify IDEA activation; validate option combinations and ordinal variable roles.
*! 2026-11-01: align FDH/IDEA statistics, effective model, saved results and documentation.
*! 2026-11-01: preserve observation rows in VRS tables; return source-row mapping.
*! 2026-11-01: fix undesirable DDF expansion; validate its sample and observed directions.
*! 2026-11-01: classify and rank stored scores by model; separate strict super-efficiency.
*! 2026-11-01: share radial bootstrap LPs, publish SE/CI, connect environmental bootstrap to sequential SW2007.
*! Validate all output names; explicit replacevars; preserve model data and automatic IDs
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
    local _c53_hint ""
    local _c53_guard 0
    capture mata: assert(_c53_context_ready()==53)
    if _rc {
        capture findfile _dea_cooper_context.ado
        if !_rc quietly run "`r(fn)'"
    }
    capture mata: assert(_c53_context_ready()==53)
    if !_rc {
        tempname entry_snapshot
        mata: `entry_snapshot'=_c53_context_save()
        local _c53_guard 1
    }
    capture mata: _dea_rec_abort()
    capture noisily _dea_run public `0'
    local rc=_rc
    if `rc' {
        capture mata: _dea_rec_abort()
        if `_c53_guard' {
            if "`_c53_hint'"=="1" mata: _c53_context_restore(`entry_snapshot')
            mata: mata drop `entry_snapshot'
        }
        return clear
        exit `rc'
    }
    if "`r(recovery_scope)'"=="package_simplex_calls" capture mata: _dea_rec_close()
    if `_c53_guard' {
        if "`_c53_hint'"=="1" mata: _c53_context_restore(`entry_snapshot')
        mata: mata drop `entry_snapshot'
    }
    return add
end

program define _dea_run, rclass
    version 16.0
    // Private entry context is not part of the public command syntax.
    gettoken _dea_call_context 0 : 0

    // =========================================================================
    // 1.0 Load IPM engine Mata functions (if Stata 16+ and not already loaded)
    // =========================================================================
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
        capture noisily _dea_parse_equal `raw0'
        local _c53_parse_rc=_rc
    }
    else {
        capture noisily _dea_parse_standard `raw0'
        local _c53_parse_rc=_rc
    }

    c_local _c53_hint "`_c53_hint'"
    if `_c53_parse_rc' exit `_c53_parse_rc'

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
    local radialtolerance = r(radial_tolerance)
    local generate  "`r(generate)'"
    local trace     "`r(trace)'"
    local recoverylog "`r(recoverylog)'"
    local saving    "`r(saving)'"
    local replace   "`r(replace)'"
    local replacevars "`r(replacevars)'"
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
    local disposability "`r(disposability)'"
    local weakmatrix "`r(weakmatrix)'"
    local xlo       "`r(xlo)'"
    local xhi       "`r(xhi)'"
    local ylo       "`r(ylo)'"
    local yhi       "`r(yhi)'"
    local xord      "`r(xord)'"
    local yord      "`r(yord)'"
    local bounded   "`r(bounded)'"
    local ordinal   "`r(ordinal)'"
    local _idea_method "`r(idea_method)'"
    local _idea_numeric "`r(idea_numeric)'"
    local _idea_lowers "`r(idea_lowers)'"
    local _idea_uppers "`r(idea_uppers)'"
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
    local congmethod "`r(congmethod)'"
    local congweight "`r(congweight)'"
    local congmatrix "`r(congmatrix)'"
    local congranges "`r(congranges)'"
    local congpolicy "`r(congpolicy)'"
    local conganchor "`r(conganchor)'"
    local congreference "`r(congreference)'"
    local congselectlimit "`r(congselectlimit)'"
    local user_set_ort = r(user_set_ort)
    if "`user_set_ort'" == "" | "`user_set_ort'" == "." local user_set_ort 0
    // sbm_ort: 0=non-oriented, 1=input-oriented, 2=output-oriented
    // Integer avoids Stata keyword conflicts (e.g. "input" is a Stata command)
    local sbm_ort_num = 0
    if `user_set_ort' {
        if "`orientation'" == "input"  local sbm_ort_num 1
        if "`orientation'" == "output" local sbm_ort_num 2
    }

    if !("`model'"=="congestion" & "`congmethod'"=="cooper" & "`engine'"=="simplex") {
    if `c(stata_version)' >= 16 {
        // Check core function AND undesirable model function (added in v2.0)
        local _ipm_needs_load 0
        capture mata: assert(_dea_ipm_audit_ready()==4)
        if _rc != 0 local _ipm_needs_load 1
        capture mata: assert(findexternal("_dea_ipm_envelopment()")!=NULL)
        if _rc != 0 local _ipm_needs_load 1
        capture mata: assert(findexternal("_dea_ipm_ddf_undesirable()")!=NULL)
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

    }

    if "`recoverylog'"!="" {
        if !(inlist("`model'","te","sbm","supersbm") | ("`model'"=="undesirable" & "`disposability'"=="weak")) | ///
           "`fdh'`_idea_method'`bootstrap'`contractcheck'`bounded'`ordinal'`xlo'`xhi'`ylo'`yhi'`xord'`yord'"!="" {
            di as error "recoverylog supports ordinary TE, SBM, super-SBM and joint weak-disposal DDF; not bootstrap, IDEA, FDH or contractcheck"
            exit 198
        }
        capture mata: assert(_dea_recovery_ready()==1)
        if _rc {
            quietly findfile _dea_recovery.ado
            quietly run "`r(fn)'"
        }
        mata: _dea_rec_begin(`c(matsize)')
    }

    if "`model'"=="congestion" & "`congmethod'"=="cooper" {
        local _dmuopt ""
        if "`dmu'"!="" local _dmuopt "dmu(`dmu')"
        _dea_congestion_cooper `ifcond' `incond', inputs(`inputs') outputs(`outputs') ///
            generate(`generate') tolerance(`tolerance') engine(`engine') ///
            congweight(`congweight') `congmatrix' `congranges' `_dmuopt' ///
            congpolicy(`congpolicy') conganchor(`conganchor') congreference(`congreference') congselectlimit(`congselectlimit') ///
            saving(`saving') `replace' `replacevars' `debug'
        return add
        return local cmd "dea"
        exit
    }

    if "`model'"=="undesirable" & "`disposability'"=="weak" {
        local _dmuopt ""
        if "`dmu'"!="" local _dmuopt "dmu(`dmu')"
        local _requested_engine "`engine'"
        if "`auto_fallback'"!="" local _requested_engine "auto"
        _dea_weakddf `ifcond' `incond', inputs(`inputs') outputs(`outputs') bad(`bad') ///
            generate(`generate') tolerance(`tolerance') engine(`_requested_engine') ///
            direction(`direction') rts(`rts') `weakmatrix' `_dmuopt' ///
            saving(`saving') `replace' `replacevars' `debug'
        if "`recoverylog'"!="" {
        tempname _rec_log _rec_stats
        mata: _dea_rec_export("`_rec_log'","`_rec_stats'")
        if `_rec_stats'[1,8]>0 return matrix recovery_log=`_rec_log'
        return matrix recovery_stats=`_rec_stats'
        return local recovery_scope "package_simplex_calls"
        }
        return add
        return local cmd "dea"
        return local cmdline `"`cmdline'"'
        exit
    }

    // =========================================================================
    // 1.3 Build estimation sample (must be done in main program)
    // =========================================================================
    tempvar touse
    mark `touse' `ifcond' `incond'
    // Resolve the endpoints actually used by both numeric IDEA paths. Keep
    // the legacy transformation plan and its public return macro unchanged.
    local _bound_vars `_idea_numeric'
    local _bound_lowers `_idea_lowers'
    local _bound_uppers `_idea_uppers'
    local _interval_points ""
    if "`_idea_method'" == "interval" {
        foreach _side in x y {
            local _dims `inputs'
            if "`_side'" == "y" local _dims `outputs'
            local _lo_option "`_side'lo"
            local _hi_option "`_side'hi"
            if "``_lo_option''``_hi_option''" == "" {
                local _interval_points `_interval_points' `_dims'
            }
            else {
                local _n : word count `_dims'
                forvalues _j = 1/`_n' {
                    local _v : word `_j' of `_dims'
                    local _lo : word `_j' of ``_lo_option''
                    local _hi : word `_j' of ``_hi_option''
                    if "`_lo'" == "" local _lo "`_v'"
                    if "`_hi'" == "" local _hi "`_v'"
                    local _bound_vars `_bound_vars' `_v'
                    local _bound_lowers `_bound_lowers' `_lo'
                    local _bound_uppers `_bound_uppers' `_hi'
                }
            }
        }
    }
    local _markout_vars `inputs' `outputs' `envars'
    if "`w'" != "" local _markout_vars `_markout_vars' `w'
    if "`p'" != "" local _markout_vars `_markout_vars' `p'
    if inlist("`model'", "undesirable", "sbmundesirable") local _markout_vars `_markout_vars' `bad'
    if "`model'" == "ddf" local _markout_vars `_markout_vars' `gx' `gy'
    if "`bounded'" != "" {
        local _markout_vars : list _markout_vars - bounded
    }
    // Validate numeric endpoints instead of dropping unused nominal values.
    // If a nominal column also has an unbounded model role, it remains required.
    local _markout_vars : list _markout_vars - _bound_vars
    local _markout_vars `_markout_vars' `_interval_points'
    markout `touse' `_markout_vars'

    // Validate selected bounded values before creating IDs, replacing results,
    // or writing files. Ordinary incomplete model rows were marked out above.
    foreach _bounded_var of local bounded {
        tempvar _bounded_lo _bounded_hi
        _dea_bounded_parse `_bounded_var', touse(`touse') ///
            lower(`_bounded_lo') upper(`_bounded_hi')
    }

    local _n_numeric : word count `_bound_vars'
    if `_n_numeric' > 0 {
        forvalues _j = 1/`_n_numeric' {
            local _v : word `_j' of `_bound_vars'
            local _lo : word `_j' of `_bound_lowers'
            local _hi : word `_j' of `_bound_uppers'
            quietly count if `touse' & (missing(`_lo', `_hi') | `_lo' > `_hi')
            if r(N) {
                di as error "numeric bounds for `_v' must be finite with lower <= upper in every selected observation"
                exit 459
            }
        }
    }

    // Standard SBM admits nonnegative data under its structural-zero rule.
    // Super-SBM ratios require strictly positive quantities; do not silently
    // perturb zeros or drop projected quantities from the super objective.
    if inlist("`model'", "sbm", "supersbm") & "`_idea_method'"=="" & "`fdh'"=="" {
        foreach _side in inputs outputs {
            tempvar _sbm_positive
            quietly generate byte `_sbm_positive'=0 if `touse'
            foreach _v of local `_side' {
                quietly count if `touse' & (`_v'<0 | ("`model'"=="supersbm" & `_v'==0))
                if r(N) {
                    if "`model'"=="supersbm" di as error "model(supersbm) requires strictly positive inputs and outputs"
                    else di as error "model(sbm) requires nonnegative inputs and outputs"
                    exit 459
                }
                quietly replace `_sbm_positive'=1 if `touse' & `_v'>0
            }
            quietly count if `touse' & `_sbm_positive'==0
            if r(N) {
                di as error "model(sbm) requires a positive input and output in every selected DMU"
                exit 459
            }
        }
    }

    // Ordinary radial TE uses the complete selected quantity sample. Check
    // original units before creating IDs/results or writing files. Zero cells
    // and columns are allowed; every DMU must have a positive value per side.
    // Flags: 1 negative quantity, 2 all-zero input row, 4 all-zero output row.
    // min()/rowmax() avoid overflow and do not use a numerical tolerance.
    // An empty sample leaves the existing r(2000) check in charge below.
    if "`_dea_call_context'"=="public" & "`model'"=="te" & "`_idea_method'`fdh'"=="" {
        mata: st_local("_te_domain",strofreal( ///
            (min(st_data(.,tokens("`inputs' `outputs'"),"`touse'"))<0) + ///
            2*(min(rowmax(st_data(.,tokens("`inputs'"),"`touse'")))==0) + ///
            4*(min(rowmax(st_data(.,tokens("`outputs'"),"`touse'")))==0)))
        if mod(`_te_domain',2) {
            di as error "model(te) requires nonnegative inputs and outputs in the selected sample"
            exit 459
        }
        if `_te_domain' {
            di as error "model(te) requires at least one positive input and one positive output in every selected DMU"
            exit 459
        }
    }

    // Specialized models validate their mathematical domains before publication.
    if "`model'"=="sbmundesirable" & "`_idea_method'"=="" & "`fdh'"=="" {
        foreach _v of varlist `inputs' `outputs' `bad' {
            quietly count if `touse' & `_v'<=0
            if r(N) {
                di as error "model(sbmundesirable) requires strictly positive inputs, good and bad outputs (slack denominators)"
                exit 459
            }
        }
    }
    if "`model'"=="congestion" & "`_idea_method'"=="" & "`fdh'"=="" {
        foreach _side in inputs outputs {
            tempvar _positive_cong
            quietly gen byte `_positive_cong'=0 if `touse'
            foreach _v of local `_side' {
                quietly count if `touse' & `_v'<0
                if r(N) {
                    di as error "model(congestion) requires nonnegative quantities"
                    exit 459
                }
                quietly replace `_positive_cong'=1 if `touse' & `_v'>0
            }
            quietly count if `touse' & `_positive_cong'==0
            if r(N) {
                di as error "model(congestion) requires a positive input and output in each selected DMU"
                exit 459
            }
        }
    }

    // Observed DDF directions must be nonnegative and nonzero for each DMU.
    // Do this before replacing any user result column or creating automatic IDs.
    if "`model'" == "undesirable" {
        foreach _v of varlist `inputs' `outputs' `bad' {
            quietly count if `touse' & `_v'<0
            if r(N) {
                di as error "model(undesirable) requires nonnegative inputs, good outputs, and bad outputs"
                exit 459
            }
        }
        local _direction_vars `inputs'
        if "`direction'" == "yb" local _direction_vars `outputs' `bad'
        if "`direction'" == "xyb" local _direction_vars `inputs' `outputs' `bad'
        tempvar _nonzero_direction
        quietly generate byte `_nonzero_direction'=0 if `touse'
        foreach _v of varlist `_direction_vars' {
            quietly replace `_nonzero_direction'=1 if `touse' & `_v'>0
        }
        quietly count if `touse' & `_nonzero_direction'==0
        if r(N) {
            di as error "model(undesirable) requires a nonzero observed direction for every selected DMU"
            exit 459
        }
    }

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
    
    // Use a temporary automatic ID until validation and computation succeed.
    // Never replace an existing user ID implicitly.
    local _generated_dmu 0
    if "`dmu'" == "" {
        capture confirm new variable _dmu
        if _rc {
            di as error "_dmu already exists; specify dmu(_dmu) to use it or supply another dmu() variable"
            exit 110
        }
        if "`generate'" == "_dmu" {
            di as error "generate(_dmu) conflicts with automatic DMU identifiers"
            exit 110
        }
        tempvar _automatic_dmu
        qui gen str32 `_automatic_dmu' = "dmu" + string(_n) if `touse'
        local dmu "`_automatic_dmu'"
        local _generated_dmu 1
    }
    
    // Count DMUs
    qui count if `touse'
    if r(N) == 0 {
        di as error "No observations found"
        exit 2000
    }
    local n_dmus = r(N)

    // Plan persistent names before computation, without modifying user data.
    // A final-engine matrix is copied only after the command has finished solving.
    local _ref_names ""
    if "`saveref'" != "" {
        if !inlist("`model'", "te", "supereff") | "`fdh'" != "" | ///
           "`_idea_method'" != "" {
            di as error "saveref requires radial model(te) or model(supereff), without FDH or imprecise-data options"
            exit 198
        }
        mata: _dea_saveref_plan("`dmu'", "`touse'")
        foreach _rv of local _ref_names {
            if "`_rv'" == "`generate'" {
                di as error "generate() conflicts with planned reference variable `_rv'"
                exit 110
            }
        }
    }

    // Validate the complete output plan before any result variable is changed.
    // Keep this plan aligned with the actual dispatch paths (including IDEA).
    local _planned `generate' `_ref_names'
    local _projection_vars ""
    if `radialtolerance'>0 {
        local _projection_vars _dea_radial _dea_radial2 _dea_radial_adj
        local _planned `_planned' `_projection_vars'
    }
    // The shared normalizer validates and chooses the IDEA path once.
    local _is_idea = ("`_idea_method'" != "")
    local _model_effective "`model'"
    if "`fdh'" != "" local _model_effective "fdh"
    local _idea_auxvars ""
    if `_is_idea' {
        local _model_effective "imprecise"
        local _idea_auxvars `generate'_opt
        local _planned `_planned' `generate'_opt
        if "`_idea_method'" == "legacy" {
            local _idea_auxvars `_idea_auxvars' `generate'_mid
            local _planned `_planned' `generate'_mid
        }
    }
    if "`saving'"!="" & "`_model_effective'"=="te" {
        local _save_columns `dmu' `inputs' `outputs' `generate' `_ref_names'
        foreach _sv in _dea_radial_status _dea_stage2_status _dea_score_source _dea_ref_stage {
            local _collision : list _sv in _save_columns
            if `_collision' {
                di as error "saving() diagnostic name `_sv' conflicts with a result column"
                exit 110
            }
        }
    }
    if "`rts'" == "vrs" & "`_model_effective'" == "te" {
        local _planned `_planned' CRS_TE VRS_TE SCALE RTS_NUM RTS
    }
    if "`model'" == "congestion" & !`_is_idea' & "`fdh'" == "" {
        local _planned `_planned' `generate'_crs `generate'_vrs ///
            `generate'_cong_eff `generate'_congestion `generate'_se
    }
    if "`model'" == "multiplier" & "`saveweights'" != "" & `stage' == 2 & !`_is_idea' & "`fdh'" == "" {
        forvalues j = 1/`m_in' {
            local _planned `_planned' _dea_v`j'
        }
        forvalues j = 1/`m_out' {
            local _planned `_planned' _dea_u`j'
        }
        local _planned `_planned' _dea_u0
    }
    if "`bootstrap'" != "" {
        local _planned `_planned' `generate'_bc
        if "`envars'" == "" local _planned `_planned' `generate'_bias `generate'_se `generate'_lo `generate'_hi
    }
    _dea_output_check, names(`_planned') ///
        protect(`inputs' `outputs' `dmu' `w' `p' `gx' `gy' `bad' ///
            `xlo' `xhi' `ylo' `yhi' `xord' `yord' `bounded' `ordinal' `envars') ///
        `replacevars'
    local _replace_outputs "`r(existing)'"
    foreach _ov of local _replace_outputs {
        quietly drop `_ov'
    }

    // Under engine(auto), avoid sending a radial DEA problem with strongly
    // heterogeneous dimensional scales through the unscaled public simplex
    // path.  Explicit engine(simplex) remains an unconditional user choice.
    local _auto_scale_ratio 1
    local _auto_preemptive_ipm 0
    local _auto_scale_threshold 1e8
    if "`auto_fallback'" != "" & "`_model_effective'" == "te" {
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
    di as text "dea v2.0.72 — Data Envelopment Analysis"
    if "`_model_effective'" == "fdh" local engine_label "FDH dominance comparison"
    if "`_idea_method'" == "interval" local engine_label "Mata LinearProgram (interval IDEA)"
    if inlist("`_model_effective'", "undesirable", "congestion", "sbmundesirable") local engine_label "Stata LinearProgram (specialized LP; engine() does not switch this backend)"
    di as text "Engine: `engine_label'"
    local _reported_orientation "`orientation'"
    if inlist("`_model_effective'", "sbm", "supersbm") & `sbm_ort_num'==0 local _reported_orientation "non-oriented"
    if "`_model_effective'"=="congestion" local _reported_orientation "input"
    if "`_model_effective'"=="sbmundesirable" local _reported_orientation "non-oriented"
    if "`_model_effective'"=="undesirable" local _reported_orientation "directional"
    di as text "Specification: model(`_model_effective') rts(`rts') orientation(`_reported_orientation') stage(`stage')"
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
        di as text "Model: `rts' `orientation' (`_model_effective')"
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
            `trace' `debug' `reference' `saveref'
    }
    // IDEA (Imprecise DEA) — with bounded/ordinal/interval data
    else if `_is_idea' {
        _dea_imprecise_dispatch, inputs(`inputs') outputs(`outputs') ///
            dmu(`dmu') rts(`rts') orientation(`orientation') ///
            stage(`stage') engine(`engine') touse(`touse') ///
            generate(`generate') tolerance(`tolerance') ///
            bounded(`bounded') ordinal(`ordinal') ///
            xlo(`xlo') xhi(`xhi') ylo(`ylo') yhi(`yhi') ///
            ideamethod(`_idea_method') ///
            numvars(`_idea_numeric') lowervars(`_idea_lowers') uppervars(`_idea_uppers') ///
            `trace' `debug'
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
            touse(`touse') generate(`generate') tolerance(`tolerance') radialtolerance(`radialtolerance') ///
            method(`method') ///
            w(`w') p(`p') gx(`gx') gy(`gy') ddfbeta(`ddfbeta') ///
            bad(`bad') direction(`direction') ///
            `reference' `saveref' `saveweights' ///
            sbmort(`sbm_ort_num') ///
            `trace' `debug'
    }

    // =========================================================================
    // 1.5a Auto fallback: simplex failed → retry all with IPM
    // =========================================================================
    // Save r() from _dea_model_dispatch before subsequent calls clear it
    local _r_N_eff    = r(n_efficient)
    local _r_N_super  = r(n_super)
    local _r_mean_eff = r(mean_efficiency)
    local _r_min_eff  = r(min_eff)
    local _r_max_eff  = r(max_eff)
    local _r_n_failed = r(n_failed)
    local _r_n_radial_failed = r(n_radial_failed)
    local _r_n_stage2_failed = r(n_stage2_failed)
    if "`_model_effective'"=="congestion" {
        foreach _cs in n_congested mean_vrs n_failed_primary n_failed_crs n_failed_vrs n_failed_cong n_failed_se {
            local _r_cong_`_cs' = r(`_cs')
        }
    }
    tempname _r_solution_source_mat
    local _r_has_solution_source 0
    capture quietly matrix list r(solution_source)
    if _rc==0 {
        matrix `_r_solution_source_mat'=r(solution_source)
        local _r_has_solution_source 1
    }
    tempname _r_stage_status_mat
    local _r_has_stage_status 0
    capture quietly matrix list r(stage_status)
    if _rc == 0 {
        matrix `_r_stage_status_mat' = r(stage_status)
        local _r_has_stage_status 1
    }
    tempname _r_ipm_trace_mat _r_vrs_ipm_trace_mat
    local _r_has_ipm_trace 0
    local _r_has_vrs_ipm_trace 0
    capture quietly matrix list r(ipm_trace)
    if _rc == 0 {
        matrix `_r_ipm_trace_mat' = r(ipm_trace)
        local _r_has_ipm_trace 1
    }
    tempname _r_projection_mat _r_projection_trace_mat _r_numerical_diagnostics_mat _r_stage_linkage_mat
    foreach _audit in projection projection_trace numerical_diagnostics stage_linkage {
        local _r_has_`_audit' 0
        capture quietly matrix list r(`_audit')
        if _rc==0 {
            matrix `_r_`_audit'_mat'=r(`_audit')
            local _r_has_`_audit' 1
        }
    }
    tempname _r_lambda_mat _r_slack_in_mat _r_slack_out_mat _r_vrsfrontier_mat _r_vrsobs_mat _r_multiplier_u0_mat
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
       "`_idea_method'" != "interval" & ///
       !inlist("`_model_effective'", "congestion", "sbmundesirable") {
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
            if `_is_idea' {
                _dea_imprecise_dispatch, inputs(`inputs') outputs(`outputs') ///
                    dmu(`dmu') rts(`rts') orientation(`orientation') ///
                    stage(`stage') engine(`engine') touse(`touse') ///
                    generate(`generate') tolerance(`tolerance') ///
                    bounded(`bounded') ordinal(`ordinal') ///
                    xlo(`xlo') xhi(`xhi') ylo(`ylo') yhi(`yhi') ///
                    ideamethod(`_idea_method') ///
                    numvars(`_idea_numeric') lowervars(`_idea_lowers') uppervars(`_idea_uppers') ///
                    `trace' `debug'
            }
            else {
                _dea_model_dispatch, inputs(`inputs') outputs(`outputs') ///
                    dmu(`dmu') rts(`rts') orientation(`orientation') ///
                    model(`model') stage(`stage') engine(`engine') ///
                    touse(`touse') generate(`generate') tolerance(`tolerance') radialtolerance(`radialtolerance') ///
                    method(`method') ///
                    w(`w') p(`p') gx(`gx') gy(`gy') ddfbeta(`ddfbeta') ///
                    bad(`bad') direction(`direction') ///
                    `reference' `saveref' `saveweights' ///
                    sbmort(`sbm_ort_num') ///
                    `trace' `debug'
            }
            // Refresh snapshots from the engine that actually supplied the final
            // scores. The initial engine's lambda/slack matrices may be stale.
            local _r_N_eff    = r(n_efficient)
            local _r_N_super  = r(n_super)
            local _r_mean_eff = r(mean_efficiency)
            local _r_min_eff  = r(min_eff)
            local _r_max_eff  = r(max_eff)
            local _r_n_failed = r(n_failed)
            local _r_n_radial_failed = r(n_radial_failed)
            local _r_n_stage2_failed = r(n_stage2_failed)
            foreach _what in lambda slack_in slack_out multiplier_u0 stage_status solution_source ipm_trace projection projection_trace numerical_diagnostics stage_linkage {
                local _r_has_`_what' 0
                capture quietly matrix list r(`_what')
                if _rc == 0 {
                    matrix `_r_`_what'_mat' = r(`_what')
                    local _r_has_`_what' 1
                }
            }
        }
    }

    // =========================================================================
    // 1.5 VRS Frontier analysis (if rts == VRS and model == te)
    // =========================================================================
    if "`rts'" == "vrs" & "`_model_effective'" == "te" {
        // Use the engine that produced the main scores.  With engine(auto),
        // this remains simplex unless the auto-fallback above switched to IPM.
        local _vrs_eng "`engine'"
        if "`recoverylog'"!="" mata: _dea_rec_scope(2)
        _dea_vrsfrontier, inputs(`inputs') outputs(`outputs') ///
            dmu(`dmu') orientation(`orientation') ///
            engine(`_vrs_eng') touse(`touse') tolerance(`tolerance') ///
            `trace' `debug'
        if "`recoverylog'"!="" mata: _dea_rec_scope(1)
        // Preserve the internal r-class results before mkmat or bootstrap
        // replaces them in the top-level command.
        matrix `_r_vrsfrontier_mat' = r(vrsfrontier)
        matrix `_r_vrsobs_mat' = r(vrsobs)
        local _r_vrs_n_crs = r(n_crs)
        local _r_vrs_n_drs = r(n_drs)
        local _r_vrs_n_irs = r(n_irs)
        local _r_has_vrsfrontier 1
        capture quietly matrix list r(vrs_ipm_trace)
        if _rc == 0 {
            matrix `_r_vrs_ipm_trace_mat' = r(vrs_ipm_trace)
            local _r_has_vrs_ipm_trace 1
        }
    }

    // =========================================================================
    // 1.6 Simar-Wilson Bootstrap (if requested)
    // =========================================================================
    local _bootvars ""
    tempname _r_bootstrap_mat
    if "`bootstrap'" != "" & "`envars'" != "" {
        tempname _swprefix
        dea_sw2007 if `touse', inputs(`inputs') outputs(`outputs') ///
            scores(`generate') envars(`envars') dmu(`dmu') ///
            rts(`rts') orientation(`orientation') generate(`_swprefix') ///
            bcreps(`reps1') reps(`reps2') level(`level') tolerance(`tolerance')
        rename `_swprefix'_eff_bc `generate'_bc
        drop `_swprefix'_eff
        local _bootvars `generate'_bc
        quietly mkmat `_bootvars' if `touse', matrix(`_r_bootstrap_mat')
        matrix colnames `_r_bootstrap_mat' = bias_corrected
    }
    else if "`bootstrap'" != "" {
        _dea_bootstrap, inputs(`inputs') outputs(`outputs') ///
            dmu(`dmu') rts(`rts') orientation(`orientation') ///
            generate(`generate') touse(`touse') ///
            reps1(`reps1') level(`level') tolerance(`tolerance') `debug'
        local _bootvars `generate'_bc `generate'_bias `generate'_se `generate'_lo `generate'_hi
        quietly mkmat `_bootvars' if `touse', matrix(`_r_bootstrap_mat')
        matrix colnames `_r_bootstrap_mat' = bias_corrected bias se ci_lo ci_hi
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
    if `_generated_dmu' {
        rename `dmu' _dmu
        local dmu "_dmu"
    }
    if `radialtolerance'>0 {
        if !`_r_has_projection' {
            di as error "projection recovery metadata was not returned"
            exit 498
        }
        foreach _v of local _projection_vars {
            quietly generate double `_v'=.
        }
        mata: st_store(selectindex(st_data(.,"`touse'")),tokens("`_projection_vars'"),st_matrix("`_r_projection_mat'")[.,1..3])
        label variable _dea_radial "Original radial factor (theta input; phi output)"
        label variable _dea_radial2 "Radial factor used by available Stage-2 projection"
        label variable _dea_radial_adj "Stage-2 radial factor minus original radial factor"
    }
    if "`saveref'" != "" {
        if !`_r_has_lambda' {
            di as error "no reference-weight matrix was returned by the selected model"
            exit 498
        }
        mata: _dea_saveref_store("`dmu'", "`touse'", "`_ref_names'", "`_r_lambda_mat'")
        matrix colnames `_r_lambda_mat' = `_ref_names'
    }

    // Compute auxiliary IDEA summaries from the final stored columns, including
    // a possible legacy-engine retry. Primary summaries came from the dispatcher.
    if `_is_idea' {
        quietly count if `touse' & abs(`generate'_opt - 1) <= `tolerance'
        local _r_N_eff_opt = r(N)
        quietly summarize `generate'_opt if `touse', meanonly
        local _r_mean_opt = r(mean)
        local _r_min_opt = r(min)
        local _r_max_opt = r(max)
        quietly count if `touse' & missing(`generate')
        local _r_n_failed_primary = r(N)
        quietly count if `touse' & missing(`generate'_opt)
        local _r_n_failed_opt = r(N)
    }

    // Save exactly once, after the final engine and persistent results are known.
    // Multiplier weights already have permanent names in the dispatcher.
    if "`saving'" != "" {
        local _saved_weights ""
        if "`_model_effective'" == "multiplier" & "`saveweights'" != "" & `stage' == 2 {
            forvalues j = 1/`m_in' {
                local _saved_weights `_saved_weights' _dea_v`j'
            }
            forvalues j = 1/`m_out' {
                local _saved_weights `_saved_weights' _dea_u`j'
            }
            local _saved_weights `_saved_weights' _dea_u0
        }
        local _save_status ""
        if `_r_has_stage_status' & `_r_has_solution_source' {
            local _save_status stage_status_mat(`_r_stage_status_mat') solution_source_mat(`_r_solution_source_mat')
        }
        local _saved_modelvars ""
        if "`_model_effective'"=="congestion" {
            foreach _suffix in crs vrs cong_eff congestion se {
                local _saved_modelvars `_saved_modelvars' `generate'_`_suffix'
            }
        }
        if "`_model_effective'"=="sbmundesirable" local _saved_modelvars `bad'
        _dea_save, modelvars(`_saved_modelvars') saving(`saving') `replace' `_save_status' ///
            inputs(`inputs') outputs(`outputs') dmu(`dmu') ///
            generate(`generate') touse(`touse') ///
            refvars(`_ref_names') weightvars(`_saved_weights') bootvars(`_bootvars') ///
            ideavars(`_idea_auxvars') projectionvars(`_projection_vars')
    }

    // NOTE: return matrix MUST come BEFORE return add in Stata 16.
    // Build efficiency matrix from generate variable
    tempname _eff_mat
    qui mkmat `generate' if `touse', matrix(`_eff_mat')
    if "`recoverylog'"!="" {
        mata: _dea_rec_finalize_data("`generate'",("`engine'"=="ipm" ? 2 : 1))
        tempname _rec_log _rec_stats
        mata: _dea_rec_export("`_rec_log'","`_rec_stats'")
        if `_rec_stats'[1,8]>0 return matrix recovery_log=`_rec_log'
        return matrix recovery_stats=`_rec_stats'
        return local recovery_scope "package_simplex_calls"
    }
    return matrix efficiency = `_eff_mat'
    if `_r_has_solution_source' {
        return matrix solution_source = `_r_solution_source_mat'
    }
    if `_r_has_stage_status' {
        return matrix stage_status = `_r_stage_status_mat'
    }
    if `_r_has_stage_linkage' {
        tempname _link_counts
        mata: st_matrix("`_link_counts'",(sum(st_matrix("`_r_stage_linkage_mat'")[.,15]:==1),sum(st_matrix("`_r_stage_linkage_mat'")[.,10]:<.)))
        return scalar n_stage_link_adjusted = `_link_counts'[1,1]
        return scalar n_stage_link_retries = `_link_counts'[1,2]
        return matrix stage_linkage = `_r_stage_linkage_mat'
    }
    if `_r_has_numerical_diagnostics' {
        tempname _diag_counts
        mata: st_matrix("`_diag_counts'",(sum(st_matrix("`_r_numerical_diagnostics_mat'")[.,5]), ///
            sum(st_matrix("`_r_numerical_diagnostics_mat'")[.,6]:==0), ///
            sum((st_matrix("`_r_numerical_diagnostics_mat'")[.,4]:==0):&(st_matrix("`_r_numerical_diagnostics_mat'")[.,5]:==0))))
        return scalar n_projection_checked = `_diag_counts'[1,1]
        return scalar n_projection_outside_tol = `_diag_counts'[1,2]
        return scalar n_projection_unchecked = `_diag_counts'[1,3]
        return matrix numerical_diagnostics = `_r_numerical_diagnostics_mat'
    }
    if `_r_has_projection' {
        tempname _projection_counts
        mata: st_matrix("`_projection_counts'",(sum(st_matrix("`_r_projection_mat'")[.,6]:<.),sum(st_matrix("`_r_projection_mat'")[.,8]),_dea_ipm_trace_counts(st_matrix("`_r_projection_trace_mat'"))[1]))
        return matrix projection = `_r_projection_mat'
        return matrix projection_trace = `_r_projection_trace_mat'
        return scalar n_projection_attempted = `_projection_counts'[1,1]
        return scalar n_projection_adjusted = `_projection_counts'[1,2]
        return scalar n_projection_solves = `_projection_counts'[1,1]+`_projection_counts'[1,3]
        return local projection_vars "`_projection_vars'"
    }
    foreach _audit in ipm_trace vrs_ipm_trace {
        if `_r_has_`_audit'' {
            tempname _trace_counts
            mata: st_matrix("`_trace_counts'",_dea_ipm_trace_counts(st_matrix("`_r_`_audit'_mat'")))
            return matrix `_audit' = `_r_`_audit'_mat'
            local _counter = cond("`_audit'"=="ipm_trace","ipm","vrs_ipm")
            return scalar n_`_counter'_retries = `_trace_counts'[1,1]
            return scalar n_`_counter'_recovered = `_trace_counts'[1,2]
            return scalar n_`_counter'_conflicts = `_trace_counts'[1,3]
        }
    }
    if "`bootstrap'" != "" {
        return matrix bootstrap = `_r_bootstrap_mat'
    }

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
        return matrix vrsobs = `_r_vrsobs_mat'
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
    return local model "`_model_effective'"
    return local model_requested "`model'"
    return local rts "`rts'"
    return local orientation "`_reported_orientation'"
    return local orientation_requested "`orientation'"
    if inlist("`_model_effective'", "undesirable", "congestion", "sbmundesirable") {
        return local lp_backend "stata_linearprogram"
        if "`_model_effective'"=="congestion" return local formulation "legacy_equalities_ratio"
        if "`_model_effective'"=="sbmundesirable" return local formulation "separable_bad_sbm"
        if "`_model_effective'"=="undesirable" {
            return local formulation "observed_direction_bad_equality"
            return local disposability "equality"
            return local direction "`direction'"
        }
    }
    return local inputs "`inputs'"
    return local outputs "`outputs'"
    return local generate "`generate'"
    return local refvars "`_ref_names'"
    if "`bootstrap'" != "" {
        return scalar bootstrap_reps = `reps1'
        if "`envars'" != "" {
            return scalar bootstrap_reps2 = `reps2'
            return local bootstrap_procedure "sw2007_sequential"
        }
        else return local bootstrap_procedure "radial"
        return scalar bootstrap_level = `level'
        return local bootstrap_engine "simplex"
        return local bootstrap_vars "`_bootvars'"
    }
    return scalar defaults_applied = `_no_option_call'
    return scalar generated_dmu = `_generated_dmu'
    return scalar N          = `n_dmus'
    return scalar N_eff      = `_r_N_eff'
    if inlist("`_model_effective'", "supereff", "supersbm") {
        return scalar N_super = `_r_N_super'
    }
    return scalar mean_eff   = `_r_mean_eff'
    return scalar min_eff    = `_r_min_eff'
    return scalar max_eff    = `_r_max_eff'
    return scalar n_failed   = `_r_n_failed'
    if "`_model_effective'"=="congestion" {
        foreach _cs in n_congested mean_vrs n_failed_primary n_failed_crs n_failed_vrs n_failed_cong n_failed_se {
            return scalar `_cs' = `_r_cong_`_cs''
        }
    }
    if `_r_has_stage_status' {
        return scalar n_radial_failed = `_r_n_radial_failed'
        return scalar n_stage2_failed = `_r_n_stage2_failed'
    }
    if `_is_idea' {
        return local idea_method "`_idea_method'"
        return local idea_vars "`generate' `_idea_auxvars'"
        return local idea_numeric_vars "`_idea_numeric'"
        return scalar N_eff_opt = `_r_N_eff_opt'
        return scalar mean_opt = `_r_mean_opt'
        return scalar min_opt = `_r_min_opt'
        return scalar max_opt = `_r_max_opt'
        return scalar n_failed_primary = `_r_n_failed_primary'
        return scalar n_failed_opt = `_r_n_failed_opt'
    }
    return scalar stage      = `stage'
    return scalar tolerance  = `tolerance'
    return scalar radial_tolerance = `radialtolerance'
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
    syntax varlist(min=1) [if] [in] [using/]
    [,
        DMU(varname)
        Rts(string)
        ORT(string)
        ORientation(string)
        STAGE(integer 2)
        FDH
        MODel(string)
        METhod(string)
        GENerate(name)
        TOLerance(string)
        RADIALTolerance(real 0)
        TOL1(string)
        TOL2(string)
        TRACE
        RECOVERYLog
        SAVing(string)
        REPLACE
        REPLACEVARS
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
        DISPosability(string)
        WEAKMATrix
        CONGMethod(string)
        CONGWeight(string)
        CONGMATrix
        CONGRanges
        CONGPolicy(string)
        CONGAnchor(string)
        CONGREFerence(string)
        CONGSelectlimit(string)

        Bounded(varlist)
        ORDinal(varlist)
        XLO(varlist)
        XHI(varlist)
        YLO(varlist)
        YHI(varlist)
        XORD(varlist numeric)
        YORD(varlist numeric)

        BOOtstrap
        ENVars(varlist numeric)
        REPS1(integer 100)
        REPS2(string)
        Level(real 95)
    ];
    #del cr

    local outputs `varlist'

    // Legacy Form A using filename is an alias for the single final save.
    if `"`using'"' != "" {
        if `"`saving'"' != "" {
            di as error "specify only one of using filename or saving()"
            exit 198
        }
        local saving `"`using'"'
    }

    // Reconcile orientation: ort() (legacy) vs orientation() (new)
    if "`ort'" != "" & "`orientation'" != "" {
        di as error "Specify only one of: ort() or orientation()"
        exit 198
    }
    if "`ort'" != "" local orientation "`ort'"

    // The old two-threshold engine is no longer used. Never ignore a request.
    if "`tol2'" != "" {
        di as error "tol2() is obsolete; use tolerance() for the integrated solver"
        exit 198
    }
    if "`tol1'" != "" & "`tolerance'" != "" {
        di as error "specify only one of tolerance() or its legacy alias tol1()"
        exit 198
    }
    if "`tol1'" != "" {
        local tolerance "`tol1'"
        di as text "note: tol1() is a deprecated alias; use tolerance()"
    }
    if "`tolerance'" == "" local tolerance 1e-7

    // Delegate to common normalizer
    capture noisily _dea_normalize_options,                           ///
        inputs(`invarnames') outputs(`outputs')       ///
        dmu(`dmu') rts(`rts') orientation(`orientation') ///
        model(`model') stage(`stage') method(`method')   ///
        fdh(`fdh') tolerance(`tolerance') radialtolerance(`radialtolerance')             ///
        generate(`generate') trace(`trace') recoverylog(`recoverylog')           ///
        saving(`saving') replace(`replace') replacevars(`replacevars') debug(`debug') ///
        engine(`engine') contractcheck(`contractcheck') ///
        reference(`reference') saveref(`saveref')     ///
        saveweights(`saveweights')                    ///
        w(`w') p(`p') gx(`gx') gy(`gy') ddfbeta(`ddfbeta') ///
        bad(`bad') direction(`direction') disposability(`disposability') weakmatrix(`weakmatrix') ///
        congmethod(`congmethod') congweight(`congweight') congmatrix(`congmatrix') congranges(`congranges') ///
        congpolicy(`congpolicy') conganchor(`conganchor') congreference(`congreference') congselectlimit(`congselectlimit') ///
        bounded(`bounded') ordinal(`ordinal')         ///
        xlo(`xlo') xhi(`xhi') ylo(`ylo') yhi(`yhi')  ///
        xord(`xord') yord(`yord')                     ///
        bootstrap(`bootstrap') envars(`envars')       ///
        reps1(`reps1') reps2(`reps2') level(`level')  ///
        ifcond(`if') incond(`in')
    local _c53_parse_rc=_rc
    c_local _c53_hint "`_c53_hint'"
    if `_c53_parse_rc' exit `_c53_parse_rc'
    return add
end

// ============================================================================
// 2.2 Form B Parser: dea y1 y2, inputs(x1 x2) [options]
// ============================================================================
program define _dea_parse_standard, rclass
    version 16.0

    #del ;
    syntax varlist(min=1) [if] [in]
    [,
        Inputs(varlist min=1)
        DMU(varname)
        Rts(string)
        ORT(string)
        ORientation(string)
        STAGE(integer 2)
        FDH
        MODel(string)
        METhod(string)
        GENerate(name)
        TOLerance(real 1e-7)
        RADIALTolerance(real 0)
        TRACE
        RECOVERYLog
        SAVing(string)
        REPLACE
        REPLACEVARS
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
        DISPosability(string)
        WEAKMATrix
        CONGMethod(string)
        CONGWeight(string)
        CONGMATrix
        CONGRanges
        CONGPolicy(string)
        CONGAnchor(string)
        CONGREFerence(string)
        CONGSelectlimit(string)

        Bounded(varlist)
        ORDinal(varlist)
        XLO(varlist)
        XHI(varlist)
        YLO(varlist)
        YHI(varlist)
        XORD(varlist numeric)
        YORD(varlist numeric)

        BOOtstrap
        ENVars(varlist numeric)
        REPS1(integer 100)
        REPS2(string)
        Level(real 95)
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
    capture noisily _dea_normalize_options,                           ///
        inputs(`inputs') outputs(`outputs')           ///
        dmu(`dmu') rts(`rts') orientation(`orientation') ///
        model(`model') stage(`stage') method(`method')   ///
        fdh(`fdh') tolerance(`tolerance') radialtolerance(`radialtolerance')             ///
        generate(`generate') trace(`trace') recoverylog(`recoverylog')           ///
        saving(`saving') replace(`replace') replacevars(`replacevars') debug(`debug') ///
        engine(`engine') contractcheck(`contractcheck') ///
        reference(`reference') saveref(`saveref')     ///
        saveweights(`saveweights')                    ///
        w(`w') p(`p') gx(`gx') gy(`gy') ddfbeta(`ddfbeta') ///
        bad(`bad') direction(`direction') disposability(`disposability') weakmatrix(`weakmatrix') ///
        congmethod(`congmethod') congweight(`congweight') congmatrix(`congmatrix') congranges(`congranges') ///
        congpolicy(`congpolicy') conganchor(`conganchor') congreference(`congreference') congselectlimit(`congselectlimit') ///
        bounded(`bounded') ordinal(`ordinal')         ///
        xlo(`xlo') xhi(`xhi') ylo(`ylo') yhi(`yhi')  ///
        xord(`xord') yord(`yord')                     ///
        bootstrap(`bootstrap') envars(`envars')       ///
        reps1(`reps1') reps2(`reps2') level(`level')  ///
        ifcond(`if') incond(`in')
    local _c53_parse_rc=_rc
    c_local _c53_hint "`_c53_hint'"
    if `_c53_parse_rc' exit `_c53_parse_rc'
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
        [RADIALTolerance(real 0)]           ///
        [GENerate(name)]                    ///
        [TRACE(string)] [RECOVERYLog(string)]                     ///
        [SAVing(string)]                    ///
        [REPLACE(string)]                   ///
        [REPLACEVARS(string)]               ///
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
        [XLO(varlist)]                      ///
        [XHI(varlist)]                      ///
        [YLO(varlist)]                      ///
        [YHI(varlist)]                      ///
        [XORD(varlist numeric)]             ///
        [YORD(varlist numeric)]             ///
        [BOOTstrap(string)]                 ///
        [ENVars(varlist numeric)]           ///
        [REPS1(integer 100)]                ///
        [REPS2(string)]               ///
        [Level(real 95)]                 ///
        [IFCOND(string)]                    ///
        [INCOND(string)]                    ///
        [BAD(varlist numeric)]              ///
        [DIRection(string)]                 ///
        [DISPosability(string)] [WEAKMATrix(string)] ///
        [CONGMethod(string)] [CONGWeight(string)] ///
        [CONGMATrix(string)] [CONGRanges(string)] ///
        [CONGPolicy(string)] [CONGAnchor(string)] [CONGREFerence(string)] [CONGSelectlimit(string)]

    local _c53_cooper = (lower("`model'")=="congestion" & lower("`congmethod'")=="cooper")
    c_local _c53_hint "`_c53_cooper'"
    // Validate the effective tolerance shared by both public syntax forms.
    // Stata missing values compare larger than ordinary numbers, so <=0 alone
    // does not exclude them. Reject before sample marking or output creation.
    if missing(`tolerance') | `tolerance' <= 0 {
        di as error "tolerance() must be positive and finite"
        exit 198
    }

    if missing(`radialtolerance') | `radialtolerance' < 0 | `radialtolerance' >= 1 {
        di as error "radialtolerance() must be finite and in [0,1)"
        exit 198
    }

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
    // Separate the bad-equality technology from convex joint weak disposal.
    local disposability = lower("`disposability'")
    if "`model'"!="undesirable" & "`disposability'`weakmatrix'"!="" {
        di as error "disposability() and weakmatrix require model(undesirable)"
        exit 198
    }
    if "`disposability'"=="" local disposability "equality"
    if !inlist("`disposability'","equality","weak") {
        di as error "disposability() must be equality or weak"
        exit 198
    }
    if "`weakmatrix'"!="" & "`disposability'"!="weak" {
        di as error "weakmatrix requires disposability(weak)"
        exit 198
    }
    if "`disposability'"=="weak" {
        if "`fdh'`w'`p'`gx'`gy'`reference'`saveref'`saveweights'`trace'"!="" {
            di as error "disposability(weak) does not support FDH, prices, gx/gy, reference/weight exports or trace; use weakmatrix/debug"
            exit 198
        }
        if `user_set_ort' {
            di as error "disposability(weak) uses direction(yb|xyb|x); omit ort()/orientation()"
            exit 198
        }
    }

    if "`congpolicy'`conganchor'`congreference'`congselectlimit'"!="" & ("`model'"!="congestion" | lower("`congmethod'")!="cooper") {
        di as error "Cooper selection options require model(congestion) congmethod(cooper)"
        exit 198
    }
    // Congestion conventions are explicit; legacy callers retain their ratio.
    local congmethod = lower("`congmethod'")
    local congweight = lower("`congweight'")
    if "`model'"!="congestion" & "`congmethod'`congweight'`congmatrix'`congranges'"!="" {
        di as error "congmethod(), congweight(), congmatrix and congranges require model(congestion)"
        exit 198
    }
    if "`congmethod'"=="" local congmethod "legacy"
    if !inlist("`congmethod'","legacy","cooper") {
        di as error "congmethod() must be legacy or cooper"
        exit 198
    }
    if "`congmethod'"=="legacy" & "`congweight'`congmatrix'`congranges'"!="" {
        di as error "congweight(), congmatrix and congranges require congmethod(cooper)"
        exit 198
    }
    if "`congmethod'"=="cooper" {
        if "`congweight'"=="" local congweight "raw"
        if !inlist("`congweight'","raw","normalized") {
            di as error "congweight() must be raw or normalized"
            exit 198
        }
        if "`fdh'`bad'`direction'`w'`p'`gx'`gy'`reference'`saveref'`saveweights'`trace'"!="" {
            di as error "Cooper congestion does not support FDH, bad/direction, price/direction vectors or reference/weight exports or trace (use debug)"
            exit 198
        }
        if `stage'!=2 | (`user_set_ort' & "`orientation'"!="output") {
            di as error "Cooper congestion requires the complete two-stage output-oriented procedure"
            exit 198
        }
    }

    // model(fdh) implies fdh flag
    if "`model'" == "fdh" local fdh "fdh"
    // Resolve the IDEA path before sample marking, automatic IDs or output
    // replacement. Both syntax forms use this same validation and dispatch plan.
    local _idea_interval = ("`xlo'`xhi'`ylo'`yhi'" != "")
    local _idea_legacy = ("`bounded'`ordinal'`xord'`yord'" != "")
    local _idea_method ""
    local _idea_numeric ""
    local _idea_lowers ""
    local _idea_uppers ""
    if `_idea_interval' | `_idea_legacy' | "`model'" == "imprecise" {
        if "`fdh'" != "" | !inlist("`model'", "te", "imprecise") {
            di as error "imprecise-data options require model(te) or model(imprecise), without fdh"
            exit 198
        }
        if !`_idea_interval' & !`_idea_legacy' {
            di as error "model(imprecise) requires bounds or ordinal specifications"
            exit 198
        }
        if `_idea_interval' {
            local _idea_method "interval"
            foreach _side in x y {
                local _needed : word count `inputs'
                if "`_side'" == "y" local _needed : word count `outputs'
                foreach _end in lo hi {
                    local _option "`_side'`_end'"
                    local _count : word count ``_option''
                    if `_count' > 0 & `_count' != `_needed' {
                        di as error "`_option'() must have `_needed' variables (same as the corresponding model varlist)"
                        exit 198
                    }
                    if !`_idea_legacy' {
                        foreach _v of local `_option' {
                            capture confirm numeric variable `_v'
                            if _rc {
                                di as error "`_option'() requires numeric bound variables"
                                exit 109
                            }
                        }
                    }
                }
            }
        }
        if `_idea_legacy' {
            local _idea_method "legacy"
            // Stata list subtraction removes occurrences, so normalize repeated
            // names before checking membership as well as before transforming.
            local xord : list uniq xord
            local yord : list uniq yord
            local _outside : list xord - inputs
            if "`_outside'" != "" {
                di as error "xord() variables must be model inputs: `_outside'"
                exit 198
            }
            local _outside : list yord - outputs
            if "`_outside'" != "" {
                di as error "yord() variables must be model outputs: `_outside'"
                exit 198
            }
            local ordinal `ordinal' `xord' `yord'
            local ordinal : list uniq ordinal
            local bounded : list uniq bounded
            local _io `inputs' `outputs'
            foreach _option in bounded ordinal {
                local _outside : list `_option' - _io
                if "`_outside'" != "" {
                    di as error "`_option'() variables must be model inputs or outputs: `_outside'"
                    exit 198
                }
            }
            foreach _v of local ordinal {
                capture confirm numeric variable `_v'
                if _rc {
                    di as error "ordinal variables must be numeric: `_v'"
                    exit 198
                }
            }
            local _overlap : list bounded & ordinal
            if "`_overlap'" != "" {
                di as error "a variable cannot be both bounded and ordinal: `_overlap'"
                exit 198
            }
            local _shared : list inputs & outputs
            local _transformed `bounded' `ordinal'
            local _ambiguous : list _shared & _transformed
            if "`_ambiguous'" != "" {
                di as error "a transformed IDEA variable must have one input/output role: `_ambiguous'"
                exit 198
            }
            // Numeric lists retain full positional alignment. Existing bounded
            // and ordinal dimensions must use their own names as placeholders;
            // a second source of bounds for the same dimension is ambiguous.
            if `_idea_interval' {
                foreach _side in x y {
                    local _dims `inputs'
                    if "`_side'" == "y" local _dims `outputs'
                    local _lo_option "`_side'lo"
                    local _hi_option "`_side'hi"
                    if "``_lo_option''``_hi_option''" != "" {
                        local _n : word count `_dims'
                        forvalues _j = 1/`_n' {
                            local _v : word `_j' of `_dims'
                            local _lo : word `_j' of ``_lo_option''
                            local _hi : word `_j' of ``_hi_option''
                            if "`_lo'" == "" local _lo "`_v'"
                            if "`_hi'" == "" local _hi "`_v'"
                            local _special : list _v in _transformed
                            if `_special' {
                                if "`_lo'" != "`_v'" | "`_hi'" != "`_v'" {
                                    di as error "bounds for bounded/ordinal variable `_v' must use `_v' itself as the positional placeholder"
                                    exit 198
                                }
                            }
                            else {
                                foreach _source in `_lo' `_hi' {
                                    capture confirm numeric variable `_source'
                                    if _rc {
                                        di as error "numeric bounds for `_v' require numeric endpoint variables"
                                        exit 109
                                    }
                                }
                                local _occurrences 0
                                foreach _candidate of local _io {
                                    if "`_candidate'" == "`_v'" local ++_occurrences
                                }
                                if `_occurrences' != 1 {
                                    di as error "numeric bounded variable `_v' must occur once in the input/output lists"
                                    exit 198
                                }
                                local _idea_numeric `_idea_numeric' `_v'
                                local _idea_lowers `_idea_lowers' `_lo'
                                local _idea_uppers `_idea_uppers' `_hi'
                            }
                        }
                    }
                }
            }
        }
    }
    // Both public syntaxes admit string names so legacy bounded outputs and
    // Form B inputs reach the same worker. Every other model variable must
    // still be numeric, before markout or any output replacement can occur.
    local _numeric_io `inputs' `outputs'
    local _numeric_io : list uniq _numeric_io
    if "`_idea_method'" == "legacy" {
        local _numeric_io : list _numeric_io - bounded
    }
    foreach _v of local _numeric_io {
        capture confirm numeric variable `_v'
        if _rc {
            di as error "model variable `_v' must be numeric; string intervals require bounded() in legacy IDEA"
            exit 109
        }
    }

    // FDH is nonconvex and uses no radial RTS constraint. Preserve the public
    // VRS label for both spellings; CRS/default is accepted as before for fdh.
    if "`fdh'" != "" {
        if !inlist("`rts'", "crs", "vrs") {
            di as error "FDH supports only default/CRS or VRS labels; NIRS/NDRS are not implemented"
            exit 198
        }
        local rts "vrs"
    }

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

    // Environmental inference uses the complete sequential SW2007 module.
    if "`envars'" != "" & "`bootstrap'" == "" {
        di as error "envars() requires bootstrap"
        exit 198
    }
    if "`reps2'" != "" & "`envars'" == "" {
        di as error "reps2() requires bootstrap and envars()"
        exit 198
    }
    if "`envars'" != "" {
        if "`reps2'" == "" local reps2 "2000"
        local _b2 = real("`reps2'")
        if missing(`_b2') | `_b2' < 2 | `_b2' != floor(`_b2') {
            di as error "reps2() must be an integer of at least 2"
            exit 198
        }
        local reps2 = `_b2'
    }
    if "`bootstrap'" != "" {
        if "`model'" != "te" | "`fdh'" != "" | "`_idea_method'" != "" {
            di as error "bootstrap requires radial model(te), without FDH or imprecise-data options"
            exit 198
        }
        if `reps1' < 50 {
            di as error "reps1() must be at least 50"
            exit 198
        }
        if missing(`level') | `level' <= 0 | `level' >= 100 {
            di as error "level() must be strictly between 0 and 100"
            exit 198
        }
    }

    // -----------------------------------------------------------------
    // 3.4 Validate Stage
    // -----------------------------------------------------------------
    if !inlist(`stage', 1, 2) {
        di as error "stage() must be 1 or 2"
        exit 198
    }
    if "`contractcheck'" != "" & ("`model'" != "te" | "`fdh'" != "" | ///
        "`_idea_method'" != "" | "`bootstrap'" != "") {
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
    if `radialtolerance'>0 & ("`engine'"!="ipm" | "`auto_fallback'"!="" | ///
        "`model'"!="te" | `stage'!=2 | "`fdh'"!="" | "`_idea_method'"!="") {
        di as error "positive radialtolerance() requires engine(ipm) model(te) stage(2), without FDH or imprecise-data options"
        exit 198
    }

    if "`model'"=="congestion" & "`congmethod'"=="cooper" {
        if "`auto_fallback'"!="" {
            di as error "Cooper engine(auto) is not implemented; choose simplex or ipm explicitly"
            exit 198
        }
        if "`engine'"!="simplex" & "`congpolicy'`conganchor'`congreference'`congselectlimit'"!="" {
            di as error "Cooper selection options require engine(simplex)"
            exit 198
        }
        if "`congpolicy'"=="" local congpolicy "stage1"
        local congpolicy=lower(strtrim("`congpolicy'"))
        if !inlist("`congpolicy'","stage1","stage2") exit 198
        if "`engine'"=="ipm" local congpolicy ""
        if "`congselectlimit'"=="" local congselectlimit "16"
        if !inlist(real("`congselectlimit'"),16,128) exit 198
        if "`congpolicy'"=="stage2" & "`congreference'"=="" exit 198
        if "`congpolicy'"!="stage2" & "`congreference'"!="" exit 198
    }

    // Both engines support all models (including undesirable) — no override needed.

    // Check: DRS/IRS with IPM need special handling (new in v2.0)
    // Simplex engine already supports DRS/IRS natively

    // -----------------------------------------------------------------
    // 3.7 Validate model-specific options
    // -----------------------------------------------------------------
    local m_in : word count `inputs'
    local m_out : word count `outputs'

    _dea_validate_model_opts, model(`model') orientation(`orientation') ///
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
    return local disposability "`disposability'"
    return local weakmatrix "`weakmatrix'"
    return local congmethod "`congmethod'"
    return local congweight "`congweight'"
    return local congmatrix "`congmatrix'"
    return local congranges "`congranges'"
    return local congpolicy "`congpolicy'"
    return local conganchor "`conganchor'"
    return local congreference "`congreference'"
    return local congselectlimit "`congselectlimit'"
    return local inputs      "`inputs'"
    return local outputs     "`outputs'"
    return local dmu         "`dmu'"
    return local rts         "`rts'"
    return local orientation "`orientation'"
    return local model       "`model'"
    return scalar stage      = `stage'
    return scalar radial_tolerance = `radialtolerance'
    return local method      "`method'"
    return scalar tolerance  = `tolerance'
    return local generate    "`generate'"
    return local recoverylog "`recoverylog'"
    return local trace       "`trace'"
    return local saving      "`saving'"
    return local replace     "`replace'"
    return local replacevars "`replacevars'"
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
    return local idea_method "`_idea_method'"
    return local idea_numeric "`_idea_numeric'"
    return local idea_lowers "`_idea_lowers'"
    return local idea_uppers "`_idea_uppers'"
    return local bootstrap   "`bootstrap'"
    return local envars      "`envars'"
    return scalar reps1      = `reps1'
    return local reps2       "`reps2'"
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
    syntax , MODEL(string) ORIentation(string) m_in(integer) m_out(integer) ///
        [W(varlist numeric)] [P(varlist numeric)]       ///
        [GX(varlist numeric)] [GY(varlist numeric)]     ///
        [XLO(varlist)] [XHI(varlist)]   ///
        [YLO(varlist)] [YHI(varlist)]   ///
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

    // Validate the composite model before creating or replacing results.
    // Missing or mismatched lists must not reach the composite LP engines.
    if "`model'" == "allocative" {
        if "`orientation'" == "input" & "`w'" == "" {
            di as error "model(allocative) with input orientation requires w()"
            exit 198
        }
        if "`orientation'" == "output" & "`p'" == "" {
            di as error "model(allocative) with output orientation requires p()"
            exit 198
        }
        if "`w'" != "" {
            local nw : word count `w'
            if `nw' != `m_in' {
                di as error "w() must have `m_in' variables (same as inputs)"
                exit 198
            }
        }
        if "`p'" != "" {
            local np : word count `p'
            if `np' != `m_out' {
                di as error "p() must have `m_out' variables (same as outputs)"
                exit 198
            }
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
        [SBMORT(integer 0)] [RADIALTolerance(real 0)]

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
    if "`model'" == "te" & `stage' == 2 {
        local _m_in : word count `inputs'
        local _m_out : word count `outputs'
        forvalues j = 1/`_m_in' {
            tempvar _si`j'
            qui gen double `_si`j'' = . if `touse'
            local slack_in_vars "`slack_in_vars' `_si`j''"
        }
        forvalues j = 1/`_m_out' {
            tempvar _so`j'
            qui gen double `_so`j'' = . if `touse'
            local slack_out_vars "`slack_out_vars' `_so`j''"
        }
    }

    // Keep TE status and provenance in this scope for both engines.
    local stage1_status_var ""
    local stage2_status_var ""
    local score_source_var ""
    local ref_stage_var ""
    local ipm_trace_mat ""
    local projection_mat ""
    local projection_trace_mat ""
    local diagnostics_mat ""
    local linkage_mat ""
    if "`model'" == "te" & "`engine'" != "ipm" {
        tempname _linkage_mat
        local linkage_mat "`_linkage_mat'"
    }
    if "`model'" == "te" & "`engine'" == "ipm" {
        tempname _ipm_trace_mat
        local ipm_trace_mat "`_ipm_trace_mat'"
        if `radialtolerance'>0 {
            tempname _projection_mat _projection_trace_mat
            local projection_mat "`_projection_mat'"
            local projection_trace_mat "`_projection_trace_mat'"
        }
    }
    if "`model'" == "te" {
        tempname _diagnostics_mat
        local diagnostics_mat "`_diagnostics_mat'"
        tempvar _score_source _ref_stage
        quietly generate byte `_score_source' = . if `touse'
        quietly generate byte `_ref_stage' = . if `touse'
        local score_source_var "`_score_source'"
        local ref_stage_var "`_ref_stage'"
        tempvar _stage1_status _stage2_status
        quietly generate byte `_stage1_status' = . if `touse'
        quietly generate byte `_stage2_status' = . if `touse'
        local stage1_status_var "`_stage1_status'"
        local stage2_status_var "`_stage2_status'"
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
            foreach _drop_output in `_wv' {
                capture confirm variable `_drop_output', exact
                if !_rc quietly drop `_drop_output'
            }
            qui gen double `_wv' = . if `touse'
            local _xname : word `j' of `inputs'
            label var `_wv' "Multiplier input weight for `_xname'"
            local weight_in_vars "`weight_in_vars' `_wv'"
        }
        forvalues j = 1/`_m_out' {
            local _wv "_dea_u`j'"
            foreach _drop_output in `_wv' {
                capture confirm variable `_drop_output', exact
                if !_rc quietly drop `_drop_output'
            }
            qui gen double `_wv' = . if `touse'
            local _yname : word `j' of `outputs'
            label var `_wv' "Multiplier output weight for `_yname'"
            local weight_out_vars "`weight_out_vars' `_wv'"
        }
        local weight_u0_var "_dea_u0"
        foreach _drop_output in `weight_u0_var' {
            capture confirm variable `_drop_output', exact
            if !_rc quietly drop `_drop_output'
        }
        qui gen double `weight_u0_var' = . if `touse'
        label var `weight_u0_var' "Multiplier intercept: u'y - v'x - u0 <= 0"
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
            slack_in_vars(`slack_in_vars') slack_out_vars(`slack_out_vars') ///
            stage1_status_var(`stage1_status_var') stage2_status_var(`stage2_status_var') ///
            score_source_var(`score_source_var') ref_stage_var(`ref_stage_var') diagnostics_mat(`diagnostics_mat') ///
            ipm_trace_mat(`ipm_trace_mat') radialtolerance(`radialtolerance') ///
            projection_mat(`projection_mat') projection_trace_mat(`projection_trace_mat') ///
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
            stage1_status_var(`stage1_status_var') stage2_status_var(`stage2_status_var') ///
            score_source_var(`score_source_var') ref_stage_var(`ref_stage_var') diagnostics_mat(`diagnostics_mat') linkage_mat(`linkage_mat') ///
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
    // All these variables are pre-created above for both engines.

    if `radialtolerance'>0 {
        mata: st_local("_adjusted",strofreal(sum(st_matrix("`projection_mat'")[.,8])))
        di as text "Stage-2 radial recovery: `_adjusted' projection(s) adjusted; original efficiency scores retained."
    }

    if "`diagnostics_mat'"!="" {
        mata: st_local("_outside",strofreal(sum(st_matrix("`diagnostics_mat'")[.,6]:==0)))
        mata: st_local("_unchecked",strofreal(sum((st_matrix("`diagnostics_mat'")[.,4]:==0):&(st_matrix("`diagnostics_mat'")[.,5]:==0))))
        if `_outside'>0 | `_unchecked'>0 {
            di as text "Returned projection audit: `_outside' outside tolerance; `_unchecked' not checkable."
            di as text "Inspect r(numerical_diagnostics); solver status alone is not an exact certificate."
        }
    }

    _dea_display, model(`model') rts(`rts') orientation(`orientation') ///
        generate(`generate') touse(`touse') ///
        tolerance(`tolerance') stage(`stage') ///
        dmu(`dmu') `reference' `debug' ///
        sbmort(`sbmort') ///
        lambda_vars(`lambda_vars') ///
        slack_in_vars(`slack_in_vars') ///
        slack_out_vars(`slack_out_vars') ///
        stage2_status_var(`stage2_status_var') ///
        inputs(`inputs') outputs(`outputs')

    // =================================================================
    // Capture display statistics BEFORE any mkmat call may clear r()
    // (mkmat is rclass and resets r(); must save stats first)
    // =================================================================
    local _n_eff_save    = r(n_efficient)
    local _n_super_save  = r(n_super)
    local _mean_eff_save = r(mean_efficiency)
    local _min_eff_save  = r(min_eff)
    local _max_eff_save  = r(max_eff)
    local _n_failed_save = r(n_failed)
    local _n_radial_failed_save = r(n_radial_failed)
    local _n_stage2_failed_save = r(n_stage2_failed)

    // =================================================================
    // Return matrices (lambda, slack_in, slack_out)
    // NOTE: return matrix must come BEFORE return add (Stata 16)
    // Tempvars in _dea_model_dispatch scope survive until after return add.
    // =================================================================
    if "`stage1_status_var'" != "" {
        tempname _stage_status_mat
        quietly mkmat `stage1_status_var' `stage2_status_var' if `touse', matrix(`_stage_status_mat')
        matrix colnames `_stage_status_mat' = radial stage2
        return matrix stage_status = `_stage_status_mat'
        tempname _source_mat
        quietly mkmat `score_source_var' `ref_stage_var' if `touse', matrix(`_source_mat')
        matrix colnames `_source_mat' = radial reference
        return matrix solution_source = `_source_mat'
    }
    if "`linkage_mat'"!="" {
        return matrix stage_linkage = `linkage_mat'
    }
    if "`diagnostics_mat'"!="" {
        return matrix numerical_diagnostics = `diagnostics_mat'
    }
    if "`ipm_trace_mat'"!="" {
        return matrix ipm_trace = `ipm_trace_mat'
    }
    if `radialtolerance'>0 {
        return matrix projection = `projection_mat'
        return matrix projection_trace = `projection_trace_mat'
    }
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
            matrix colnames `_sin_mat' = `inputs'
            return matrix slack_in = `_sin_mat'
        }
    }
    if "`slack_out_vars'" != "" {
        tempname _sout_mat
        capture qui mkmat `slack_out_vars' if `touse', matrix(`_sout_mat')
        if _rc == 0 {
            matrix colnames `_sout_mat' = `outputs'
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
    return scalar n_super         = `_n_super_save'
    return scalar mean_efficiency = `_mean_eff_save'
    return scalar min_eff         = `_min_eff_save'
    return scalar max_eff         = `_max_eff_save'
    return scalar n_failed        = `_n_failed_save'
    if "`stage1_status_var'" != "" {
        return scalar n_radial_failed = `_n_radial_failed_save'
        return scalar n_stage2_failed = `_n_stage2_failed_save'
    }
end


********************************************************************************
* 6. FDH DISPATCHER (Branch & Bound, engine-independent)
********************************************************************************

program define _dea_fdh_dispatch, rclass
    syntax , Inputs(string) Outputs(string) DMU(varname) ///
        ORIentation(string) STAGE(integer) ///
        TOUSE(varname) GENerate(name) TOLerance(real) ///
        [TRACE] [DEBUG] [REFerence] [SAVERef]

    // FDH uses dominance comparison; no LP solver is needed.
    // The Mata implementation is loaded from _dea_ipm.ado.
    
    // Ensure Mata functions are loaded
    capture mata: assert(findexternal("_dea_fdh()")!=NULL)
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
    
    qui count if abs(`generate' - 1) <= `tolerance' & `touse'
    local n_efficient = r(N)
    
    qui count if `touse' & missing(`generate')
    local n_failed = r(N)

    // Display results
    di as text "{hline 50}"
    di as text "Model: FDH `orientation'"
    di as text "DMUs evaluated:  " as result `n_dmus'
    di as text "Efficient DMUs:  " as result `n_efficient' ///
        as text " (" %5.1f 100*`n_efficient'/`n_dmus' "%)"
    
    qui sum `generate' if `touse', detail
    local mean_eff = r(mean)
    local min_eff = r(min)
    local max_eff = r(max)
    di as text "Mean efficiency:  " as result %9.6f `mean_eff'
    di as text "Min efficiency:   " as result %9.6f `min_eff'
    di as text "Max efficiency:   " as result %9.6f `max_eff'
    di as text "{hline 50}"
    
    // Display DMU results
    di as text _n "{bf:DMU Efficiency Scores}"
    di as text "{hline 50}"
    
    tempvar eff_rounded rank_var
    qui gen double `eff_rounded' = round(`generate', `tolerance') if `touse'
    qui egen `rank_var' = rank(-`eff_rounded') if `touse', track

    // Display data expressions, never expand label contents into Stata syntax.
    local dmu_text `dmu'
    capture confirm string variable `dmu'
    if _rc {
        tempvar dmu_text
        quietly gen strL `dmu_text' = strtrim(string(`dmu', "%18.0g")) if `touse'
    }
    
    preserve
    qui keep if `touse'
    local nobs = _N
    di as text %12s "DMU" " {c |}" %12s "Theta" %8s "Rank"
    di as text "{hline 13}{c +}{hline 20}"
    
    forvalues i = 1/`nobs' {
        local eff_val = `generate'[`i']
        local rank_val = `rank_var'[`i']
        if !missing(`eff_val') {
            di as text %12s _asis `dmu_text'[`i'] in smcl " {c |}" ///
                as result %12.6f `eff_val' %8.0f `rank_val'
        }
        else {
            di as text %12s _asis `dmu_text'[`i'] in smcl " {c |}" ///
                as result %12s "FAILED" %8s "."
        }
    }
    di as text "{hline 34}"
    restore
    
    // Return values
    return scalar N = `n_dmus'
    return scalar n_efficient = `n_efficient'
    return scalar mean_efficiency = `mean_eff'
    return scalar min_eff = `min_eff'
    return scalar max_eff = `max_eff'
    return scalar n_failed = `n_failed'
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
        GENerate(name) TOLerance(real) IDEAMethod(string) ///
        [NUMvars(varlist numeric)] [LOWERvars(varlist numeric)] [UPPERvars(varlist numeric)] ///
        [Bounded(varlist)] [ORDinal(varlist)] ///
        [XLO(varlist)] [XHI(varlist)] ///
        [YLO(varlist)] [YHI(varlist)] ///
        [TRACE] [DEBUG]

    di as text _n "IDEA (Imprecise DEA)"
    di as text "{hline 50}"
    
    // Ensure Mata functions are loaded
    capture mata: assert(findexternal("_dea_idea()")!=NULL)
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
    
    // The public normalizer validated bounds and selected this path.
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
    if "`ideamethod'" == "legacy" {
        _dea_imprecise_legacy, inputs(`inputs') outputs(`outputs') ///
            dmu(`dmu') rts(`rts') orientation(`orientation') ///
            stage(`stage') engine(`engine') touse(`touse') ///
            generate(`generate') tolerance(`tolerance') ///
            bounded(`bounded') ordinal(`ordinal') ///
            numvars(`numvars') lowervars(`lowervars') uppervars(`uppervars') ///
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
    foreach _drop_output in `opt_var' {
        capture confirm variable `_drop_output', exact
        if !_rc quietly drop `_drop_output'
    }
    qui gen double `opt_var' = `eff_opt' if `touse'
    
    // Count results
    qui count if `touse'
    local n_dmus = r(N)
    
    qui count if abs(`generate' - 1) <= `tolerance' & `touse'
    local n_eff_pess = r(N)
    
    qui count if abs(`eff_opt' - 1) <= `tolerance' & `touse'
    local n_eff_opt = r(N)
    
    qui count if `touse' & missing(`generate', `opt_var')
    local n_failed = r(N)

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

    // Numeric IDs retain the prior general-format representation.
    local dmu_text `dmu'
    capture confirm string variable `dmu'
    if _rc {
        tempvar dmu_text
        quietly gen strL `dmu_text' = strtrim(string(`dmu', "%18.0g")) if `touse'
    }
    
    preserve
    qui keep if `touse'
    local nobs = _N
    
    forvalues i = 1/`nobs' {
        local pess_val = `generate'[`i']
        local opt_val = `eff_opt'[`i']
        if !missing(`pess_val', `opt_val') {
            local interval = `opt_val' - `pess_val'
            di as text %12s _asis `dmu_text'[`i'] in smcl " {c |}" ///
                as result %14.6f `pess_val' %14.6f `opt_val' %10.6f `interval'
        }
        else {
            di as text %12s _asis `dmu_text'[`i'] in smcl " {c |}" ///
                as result %14.6f `pess_val' %14.6f `opt_val' %10s "FAILED"
        }
    }
    di as text "{hline 60}"
    restore
    
    // Return values
    return scalar N = `n_dmus'
    return scalar n_efficient = `n_eff_pess'
    return scalar mean_efficiency = `mean_pess'
    return scalar min_eff = `min_pess'
    return scalar max_eff = `max_pess'
    return scalar n_failed = `n_failed'
    return scalar n_efficient_pess = `n_eff_pess'
    return scalar n_efficient_opt = `n_eff_opt'
    return scalar mean_pess = `mean_pess'
    return scalar mean_opt = `mean_opt'
    return local model "idea"
    return local orientation "`orientation'"
    return local rts "`rts'"
    return add
end

// Parse a bounded model variable without modifying it. The same parser is
// used by preflight and the legacy worker, so whitespace and numeric notation
// cannot be accepted at entry and interpreted differently during computation.
capture program drop _dea_bounded_parse
program define _dea_bounded_parse
    version 16.0
    syntax varname, LOWER(name) UPPER(name) [TOUSE(varname)]

    if "`touse'" == "" {
        tempvar _all
        quietly generate byte `_all' = 1
        local touse "`_all'"
    }
    capture confirm numeric variable `varlist'
    if !_rc {
        quietly generate double `lower' = `varlist' if `touse'
        quietly generate double `upper' = `varlist' if `touse'
    }
    else {
        tempvar _text _comma _valid
        local _number "[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?"
        local _space "[[:space:]]*"
        local _pattern "^\[`_space'`_number'`_space',`_space'`_number'`_space'\]$"
        quietly generate strL `_text' = ustrtrim(`varlist') if `touse'
        quietly generate byte `_valid' = ustrregexm(`_text', "`_pattern'") if `touse'
        quietly generate long `_comma' = strpos(`_text', ",") if `touse'
        quietly generate double `lower' = real(ustrtrim(substr(`_text', 2, `_comma'-2))) ///
            if `touse' & `_valid'
        quietly generate double `upper' = real(ustrtrim(substr(`_text', `_comma'+1, ///
            strlen(`_text')-`_comma'-1))) if `touse' & `_valid'
    }
    quietly count if `touse' & (missing(`lower', `upper') | `lower' > `upper')
    if r(N) {
        local _invalid = r(N)
        di as error "bounded() variable `varlist' has invalid values in `_invalid' selected observation(s)"
        di as error "use [lower,upper] with finite lower <= upper, or nonmissing numeric point bounds"
        exit 459
    }
end

capture program drop _dea_imprecise_legacy
program define _dea_imprecise_legacy, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) DMU(varname) ///
        RTS(string) ORIentation(string) STAGE(integer) ///
        ENGINE(string) TOUSE(varname) GENerate(name) TOLerance(real) ///
        [Bounded(varlist)] [ORDinal(varlist)] ///
        [NUMvars(varlist numeric)] [LOWERvars(varlist numeric)] [UPPERvars(varlist numeric)] ///
        [TRACE] [DEBUG]

    tempvar __idea_row __idea_eff __idea_score
    tempname __idea_cutoff
    qui gen long `__idea_row' = _n

    tempfile __idea_base __idea_results __one

    preserve
        qui keep if `touse'
        qui keep `__idea_row' `dmu' `inputs' `outputs' `bounded' `ordinal' `lowervars' `uppervars'
        // Freeze all endpoint columns before modifying any model dimension.
        // A bound source may itself be another input/output being transformed.
        local __numeric_n : word count `numvars'
        local __numeric_lo ""
        local __numeric_hi ""
        if `__numeric_n' > 0 {
            forvalues j = 1/`__numeric_n' {
                local __lo_source : word `j' of `lowervars'
                local __hi_source : word `j' of `uppervars'
                tempvar __lo_copy __hi_copy
                qui gen double `__lo_copy' = `__lo_source'
                qui gen double `__hi_copy' = `__hi_source'
                local __numeric_lo `__numeric_lo' `__lo_copy'
                local __numeric_hi `__numeric_hi' `__hi_copy'
            }
        }
        // Raw endpoint columns are no longer needed. Keep only model data and
        // the frozen copies so auxiliary names cannot collide with inner DEA
        // outputs such as CRS_TE during a VRS solve.
        qui keep `__idea_row' `dmu' `inputs' `outputs' `bounded' `ordinal' `__numeric_lo' `__numeric_hi'
        qui save `__idea_base', replace

        clear
        qui set obs 0
        qui gen long `__idea_row' = .
        qui gen double `__idea_eff' = .
        qui save `__idea_results', replace

        use `__idea_base', clear
        local __n = _N

        forvalues dmui = 1/`__n' {
            use `__idea_base', clear

            if `__numeric_n' > 0 {
                forvalues j = 1/`__numeric_n' {
                    local __v : word `j' of `numvars'
                    local __lo : word `j' of `__numeric_lo'
                    local __hi : word `j' of `__numeric_hi'
                    // Nominal model variables may have narrower storage types.
                    qui recast double `__v'
                    qui replace `__v' = `__lo'
                    if strpos(" `inputs' ", " `__v' ") {
                        qui replace `__v' = `__hi' if _n != `dmui'
                    }
                    else {
                        qui replace `__v' = `__hi' if _n == `dmui'
                    }
                }
            }

            if "`bounded'" != "" {
                foreach boundedVar of varlist `bounded' {
                    tempvar __idea_upper
                    capture confirm string variable `boundedVar'
                    if !_rc {
                        tempvar __lo
                        _dea_bounded_parse `boundedVar', lower(`__lo') upper(`__idea_upper')
                        qui drop `boundedVar'
                        qui gen double `boundedVar' = `__lo'
                    }
                    else {
                        qui gen double `__idea_upper' = `boundedVar'
                    }

                    if strpos(" `inputs' ", " `boundedVar' ") {
                        // For bounded inputs, the evaluated DMU uses the lower
                        // value and reference DMUs use the upper value.
                        qui replace `boundedVar' = `__idea_upper' if _n != `dmui'
                    }
                    else if strpos(" `outputs' ", " `boundedVar' ") {
                        // For bounded outputs, the evaluated DMU uses the upper
                        // value and reference DMUs use the lower value.
                        qui replace `boundedVar' = `__idea_upper' if _n == `dmui'
                    }
                    qui drop `__idea_upper'
                }
            }

            if "`ordinal'" != "" {
                foreach ordVar of varlist `ordinal' {
                    // Keep the evaluated DMU's original value fixed. In-place
                    // replacement of that row must not change later comparisons.
                    scalar `__idea_cutoff' = `ordVar'[`dmui']
                    if strpos(" `inputs' ", " `ordVar' ") {
                        qui replace `ordVar' = . if `ordVar' < scalar(`__idea_cutoff')
                        qui replace `ordVar' = 1 if `ordVar' < .
                        qui replace `ordVar' = _N if `ordVar' == .
                    }
                    else if strpos(" `outputs' ", " `ordVar' ") {
                        qui replace `ordVar' = . if `ordVar' <= scalar(`__idea_cutoff')
                        qui replace `ordVar' = 0 if `ordVar' < .
                        qui replace `ordVar' = 1 if `ordVar' == .
                    }
                }
            }

            foreach _drop_output in `__idea_score' {
                capture confirm variable `_drop_output', exact
                if !_rc quietly drop `_drop_output'
            }
            // Legacy ordinal transformations may create zero-output reference rows.
            // Keep their historical LP domain separate from public ordinary TE.
            quietly _dea_run legacy_idea `inputs' = `outputs', rts(`rts') ort(`orientation') ///
                stage(`stage') gen(`__idea_score') engine(`engine') dmu(`dmu') ///
                tolerance(`tolerance')

            qui keep in `dmui'
            qui keep `__idea_row' `__idea_score'
            qui rename `__idea_score' `__idea_eff'
            qui save `__one', replace

            use `__idea_results', clear
            qui append using `__one'
            qui save `__idea_results', replace
        }
    restore

    qui merge 1:1 `__idea_row' using `__idea_results', nogen

    foreach _drop_output in `generate' {
        capture confirm variable `_drop_output', exact
        if !_rc quietly drop `_drop_output'
    }
    qui gen double `generate' = `__idea_eff' if `touse'

    local opt_var "`generate'_opt"
    local mid_var "`generate'_mid"
    foreach _drop_output in `opt_var' {
        capture confirm variable `_drop_output', exact
        if !_rc quietly drop `_drop_output'
    }
    foreach _drop_output in `mid_var' {
        capture confirm variable `_drop_output', exact
        if !_rc quietly drop `_drop_output'
    }
    qui gen double `opt_var' = `__idea_eff' if `touse'
    qui gen double `mid_var' = `__idea_eff' if `touse'

    qui count if `touse'
    local n_dmus = r(N)
    qui count if abs(`generate' - 1) <= `tolerance' & `touse'
    local n_eff = r(N)
    qui count if `touse' & missing(`generate')
    local n_failed = r(N)
    qui sum `generate' if `touse', detail
    local mean_eff = r(mean)
    local min_eff = r(min)
    local max_eff = r(max)

    di as text "Legacy bounded/ordinal IDEA transformation"
    if "`numvars'" != "" {
        di as text "Numeric bounded variables: `numvars'"
        di as text "Lower endpoints: `lowervars'"
        di as text "Upper endpoints: `uppervars'"
    }
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
    return scalar n_failed = `n_failed'
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

    tempname vrs_trace
    // VRS decomposition uses only radial scores; no slack-maximization LP is needed.
    if "`engine'" == "simplex" {
        quietly mata: _dea_simplex_envelopment("`inputs'", "`outputs'", ///
            "`crs_te'", "`touse'", "`orientation'", "crs", ///
            `tolerance', "", "", "", 1)
        quietly mata: _dea_simplex_envelopment("`inputs'", "`outputs'", ///
            "`vrs_te'", "`touse'", "`orientation'", "vrs", ///
            `tolerance', "", "", "", 1)
        quietly mata: _dea_simplex_envelopment("`inputs'", "`outputs'", ///
            "`drs_te'", "`touse'", "`orientation'", "drs", ///
            `tolerance', "", "", "", 1)

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
            "`touse'", "`orientation'", `tolerance', "`vrs_trace'")
    }

    // Count DMUs
    qui count if `touse'
    local n_dmus = r(N)

    // IDs are display labels, not aggregation keys. Preserve every selected
    // observation, its missing values, and the order used by r(efficiency).
    tempname vrsfrontier vrsobs
    tempvar source_obs
    quietly generate double `source_obs' = _n if `touse'
    quietly mkmat `crs_te' `vrs_te' `drs_te' `scale_eff' `rts_type' ///
        if `touse', matrix(`vrsfrontier')
    quietly mkmat `source_obs' if `touse', matrix(`vrsobs')
    matrix colnames `vrsfrontier' = CRS_TE VRS_TE DRS_TE SCALE RTS
    matrix colnames `vrsobs' = obsno

    // Display
    di as text _n "VRS Frontier Analysis"
    di as text "{hline 65}"
    di as text "  RTS codes: -1 = DRS, 0 = CRS, 1 = IRS"
    di as text "  Rows follow the selected observations; r(vrsobs) maps source rows."
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
    foreach _drop_output in CRS_TE VRS_TE SCALE RTS_NUM RTS {
        capture confirm variable `_drop_output', exact
        if !_rc quietly drop `_drop_output'
    }
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
    return matrix vrsobs = `vrsobs'
    if "`engine'"!="simplex" return matrix vrs_ipm_trace = `vrs_trace'
    return scalar n_crs = `n_crs'
    return scalar n_drs = `n_drs'
    return scalar n_irs = `n_irs'
    return add
end


********************************************************************************
* 9. RESULT DISPLAY
********************************************************************************

// One score convention for the common model display: DDF/additive scores
// measure inefficiency; both super-efficiency orientations store larger-is-better
// scores (output radial super-efficiency has already been inverted to 1/phi).
program define _dea_score_summary, rclass
    version 16.0
    syntax , MODEL(string) SCORE(varname numeric) TOUSE(varname) ///
        TOLerance(real) RANKvar(name)

    local zero_frontier = inlist("`model'", "ddf", "undesirable", "additive")
    local super_model = inlist("`model'", "supereff", "supersbm")
    local frontier = 1-`zero_frontier'
    tempvar finite boundary rounded
    quietly gen byte `finite' = `touse' & !missing(`score')
    quietly gen byte `boundary' = abs(`score'-`frontier')<=`tolerance' if `finite'
    if `super_model' {
        quietly count if `finite' & `score'>=1-`tolerance'
    }
    else quietly count if `finite' & `boundary'==1
    local n_efficient=r(N)

    local n_super=.
    if `super_model' {
        quietly count if `finite' & `score'>1+`tolerance'
        local n_super=r(N)
    }
    quietly count if `touse' & missing(`score')
    local n_failed=r(N)

    // Snap only the temporary ranking key to the frontier. Stored scores stay
    // unchanged. Elsewhere retain the existing tolerance-grid competition ranks.
    quietly gen double `rounded'=round(`score',`tolerance') if `finite'
    quietly replace `rounded'=`frontier' if `finite' & `boundary'==1
    if `zero_frontier' quietly egen double `rankvar'=rank(`rounded') if `finite', track
    else quietly egen double `rankvar'=rank(-`rounded') if `finite', track

    return scalar n_efficient=`n_efficient'
    return scalar n_super=`n_super'
    return scalar n_failed=`n_failed'
    return scalar zero_frontier=`zero_frontier'
end

program define _dea_display, rclass
    syntax , MODEL(string) RTS(string) ORIentation(string) ///
        GENerate(name) TOUSE(varname) TOLerance(real) ///
        STAGE(integer) DMU(varname) ///
        [REFerence] [DEBUG] ///
        [lambda_vars(string)] [slack_in_vars(string)] [slack_out_vars(string)] ///
        [stage2_status_var(string)] ///
        [Inputs(string)] [Outputs(string)] [SBMORT(integer 0)]

    // Summary statistics
    qui sum `generate' if `touse', detail
    local mean_eff = r(mean)
    local min_eff  = r(min)
    local max_eff  = r(max)
    local sd_eff   = r(sd)

    qui count if `touse'
    local n_dmus = r(N)

    tempvar rank_var
    _dea_score_summary, model(`model') score(`generate') touse(`touse') ///
        tolerance(`tolerance') rankvar(`rank_var')
    local n_efficient=r(n_efficient)
    local n_super=r(n_super)
    local n_failed=r(n_failed)
    local zero_frontier=r(zero_frontier)
    local n_radial_failed=`n_failed'
    local n_stage2_failed=0
    if "`stage2_status_var'" != "" & `stage' == 2 {
        quietly count if `touse' & `stage2_status_var'!=0
        local n_stage2_failed=r(N)
        quietly count if `touse' & (missing(`generate') | `stage2_status_var'!=0)
        local n_failed=r(N)
    }

    // Display
    di as text _n "DEA Results"
    di as text "{hline 50}"
    if inlist("`model'", "sbm", "supersbm") & `sbmort' == 1 {
        di as text "Model: `rts' input-oriented (`model')"
    }
    else if inlist("`model'", "sbm", "supersbm") & `sbmort' == 2 {
        di as text "Model: `rts' output-oriented (`model')"
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
    if inlist("`model'", "supereff", "supersbm") {
        di as text "Strictly super-efficient: " as result `n_super' ///
            as text " (score > 1 + tolerance)"
    }
    if `zero_frontier' {
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
        if "`stage2_status_var'" != "" & `stage' == 2 {
            di as error _n "Warning: `n_failed' DMU(s) have incomplete requested results"
            di as text "Radial scores missing: `n_radial_failed'; Stage 2 unavailable: `n_stage2_failed'"
            di as text "Available radial scores are retained; unavailable projections have missing references and slacks."
        }
        else di as error _n "Warning: `n_failed' DMU(s) failed to converge"
    }
    di as text "{hline 50}"

    // ---- Tabular DMU-level results ----
    if `zero_frontier' di as text _n "{bf:DMU Inefficiency Scores}"
    else di as text _n "{bf:DMU Efficiency Scores}"
    di as text "{hline 50}"

    // Build a matrix-style result display
    tempname dearslt
    qui count if `touse'
    local nobs = r(N)

    // Keep labels in data: quotes/macros/SMCL in an ID are literal display text.
    local dmu_text `dmu'
    capture confirm string variable `dmu'
    if _rc {
        tempvar dmu_text
        quietly gen strL `dmu_text' = strtrim(string(`dmu', "%18.0g")) if `touse'
    }

    // Show compact list
    preserve
    qui keep if `touse'
    di as text %12s "DMU" " {c |}" %12s "Score" %8s "Rank"
    di as text "{hline 13}{c +}{hline 20}"

    forvalues i = 1/`nobs' {
        local eff_val = `generate'[`i']
        local rank_val = `rank_var'[`i']
        if !missing(`generate'[`i']) {
            di as text %12s _asis `dmu_text'[`i'] in smcl " {c |}" ///
                as result %12.6f `eff_val' %8.0f `rank_val'
        }
        else {
            di as text %12s _asis `dmu_text'[`i'] in smcl " {c |}" ///
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
            local ref_count = 0

            forvalues j = 1/`n_obs' {
                local lambda_var : word `j' of `lambda_vars'
                if "`lambda_var'" != "" {
                    local lam_val = `lambda_var'[`k']
                    if `lam_val' > `tolerance' & `lam_val' < . {
                        local ref_count = `ref_count' + 1
                        if `ref_count' == 1 {
                            di as text _asis `dmu_text'[`k'] in smcl ": " _continue
                        }
                        if `ref_count' <= 6 {
                            di as result " " _asis `dmu_text'[`j'] in smcl ///
                                "(" %6.4f `lam_val' ")" _continue
                        }
                    }
                }
            }
            if `ref_count' > 0 {
                if `ref_count' > 6 {
                    local more = `ref_count' - 6
                    di as text " ... (+`more' more)" _continue
                }
                di as text ""
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
    return scalar n_super         = `n_super'
    return scalar n_failed        = `n_failed'
    return scalar n_radial_failed = `n_radial_failed'
    return scalar n_stage2_failed = `n_stage2_failed'
end


********************************************************************************
* 10. RESULT SAVING
********************************************************************************

program define _dea_save
    syntax , SAVing(string) Inputs(string) Outputs(string) ///
        DMU(varname) GENerate(name) TOUSE(varname) ///
        [REPLACE] [REFVARS(varlist numeric)] [WEIGHTVARS(varlist numeric)] ///
        [BOOTVARS(varlist numeric)] [IDEAVARS(varlist numeric)] [PROJECTIONVARS(varlist numeric)] ///
        [stage_status_mat(name)] [solution_source_mat(name)] [MODELVARS(varlist numeric)]

    preserve
    qui keep if `touse'
    qui keep `dmu' `inputs' `outputs' `generate' `refvars' `weightvars' `bootvars' `ideavars' `projectionvars' `modelvars'
    if "`stage_status_mat'"!="" & "`solution_source_mat'"!="" {
        foreach _sv in _dea_radial_status _dea_stage2_status _dea_score_source _dea_ref_stage {
            quietly generate byte `_sv' = .
        }
        mata: st_store(.,tokens("_dea_radial_status _dea_stage2_status _dea_score_source _dea_ref_stage"), ///
            (st_matrix("`stage_status_mat'"),st_matrix("`solution_source_mat'")))
        label variable _dea_radial_status "Radial score: 0 available, 1 missing"
        label variable _dea_stage2_status "Stage 2: 0 success, positive failure, missing not attempted"
        label variable _dea_score_source "Radial source: 0 unavailable, 1 primal, 2 dual"
        label variable _dea_ref_stage "Reference source: 0 unavailable, 1 stage 1, 2 stage 2"
    }
    qui save "`saving'", `replace'
    restore
    di as text _n "Results saved to: `saving'"
end


********************************************************************************
* 11. SIMAR-WILSON BOOTSTRAP DISPATCHER
********************************************************************************

program define _dea_bootstrap, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) DMU(varname) ///
        RTS(string) ORIentation(string) GENerate(name) TOUSE(varname) ///
        REPS1(integer) Level(real) TOLerance(real) [DEBUG]

    // Output names and sample were validated by the main command before solving.
    // The shared worker uses the main engine's point scores as its initial values;
    // all pseudo-frontier LPs use simplex, as in the dedicated command.
    quietly count if `touse' & (missing(`generate') | `generate'<=0)
    if r(N) {
        di as error "bootstrap requires positive, nonmissing radial scores for every selected DMU"
        exit 498
    }
    foreach suffix in bc bias se lo hi {
        quietly generate double `generate'_`suffix' = .
    }
    _dea_bootstrap_core, inputs(`inputs') outputs(`outputs') dmu(`dmu') ///
        touse(`touse') eff(`generate') initial(`generate') ///
        bc(`generate'_bc) se(`generate'_se) lo(`generate'_lo) hi(`generate'_hi) ///
        orientation(`orientation') rts(`rts') reps(`reps1') level(`level') ///
        tolerance(`tolerance') `debug'
    quietly replace `generate'_bias = `generate' - `generate'_bc if `touse'
    label variable `generate'_bc "Bootstrap bias-corrected radial efficiency"
    label variable `generate'_bias "Bootstrap estimated bias (original - corrected)"
    label variable `generate'_se "Bootstrap standard error"
    label variable `generate'_lo "Bootstrap `level'% pivot CI: lower endpoint"
    label variable `generate'_hi "Bootstrap `level'% pivot CI: upper endpoint"
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
        [slack_in_vars(string)] [slack_out_vars(string)] ///
        [stage1_status_var(string)] [stage2_status_var(string)] ///
        [score_source_var(string)] [ref_stage_var(string)] [diagnostics_mat(string)] ///
        [ipm_trace_mat(string)] ///
        [RADIALTolerance(real 0)] [projection_mat(string)] [projection_trace_mat(string)] ///
        [weight_in_vars(string)] [weight_out_vars(string)] [weight_u0_var(string)] ///
        [TRACE] [DEBUG] [SBMORT(integer 0)]

    // Prepare efficiency variable
    foreach _drop_output in `generate' {
        capture confirm variable `_drop_output', exact
        if !_rc quietly drop `_drop_output'
    }

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
            slack_in_vars(`slack_in_vars') slack_out_vars(`slack_out_vars') ///
            stage1_status_var(`stage1_status_var') stage2_status_var(`stage2_status_var') ///
            score_source_var(`score_source_var') ref_stage_var(`ref_stage_var') diagnostics_mat(`diagnostics_mat') ///
            ipm_trace_mat(`ipm_trace_mat') radialtolerance(`radialtolerance') ///
            projection_mat(`projection_mat') projection_trace_mat(`projection_trace_mat') ///
            `reference' `debug'
    }
    else if "`model'" == "cost" {
        _dea_ipm_cost, inputs(`inputs') outputs(`outputs') ///
            efficiency(`efficiency_score') touse(`touse') ///
            rts(`rts') tolerance(`tolerance') w(`w') `debug'
    }
    else if "`model'" == "revenue" {
        _dea_ipm_revenue, inputs(`inputs') outputs(`outputs') ///
            efficiency(`efficiency_score') touse(`touse') ///
            rts(`rts') tolerance(`tolerance') p(`p') `debug'
    }
    else if "`model'" == "profit" {
        _dea_ipm_profit, inputs(`inputs') outputs(`outputs') ///
            efficiency(`efficiency_score') touse(`touse') ///
            rts(`rts') tolerance(`tolerance') w(`w') p(`p') `debug'
    }
    else if "`model'" == "ddf" {
        _dea_ipm_ddf, inputs(`inputs') outputs(`outputs') ///
            efficiency(`efficiency_score') touse(`touse') ///
            rts(`rts') tolerance(`tolerance') gx(`gx') gy(`gy') `debug'
    }
    else if "`model'" == "undesirable" {
        _dea_ipm_ddf_undesirable, inputs(`inputs') outputs(`outputs') ///
            bad(`bad') efficiency(`efficiency_score') touse(`touse') ///
            rts(`rts') tolerance(`tolerance') direction(`direction') `debug'
    }
    else if "`model'" == "sbm" {
        if `sbmort' == 0 {
            _dea_ipm_sbm, inputs(`inputs') outputs(`outputs') ///
                efficiency(`efficiency_score') touse(`touse') ///
                rts(`rts') tolerance(`tolerance') `debug'
        }
        else {
            local _sbm_ort_str = cond(`sbmort'==1, "input", "output")
            mata: _dea_ipm_sbm_oriented("`inputs'", "`outputs'", ///
                "`efficiency_score'", "`touse'", "`rts'", "`_sbm_ort_str'", `tolerance')
        }
    }
    else if "`model'" == "additive" {
        _dea_ipm_additive, inputs(`inputs') outputs(`outputs') ///
            efficiency(`efficiency_score') touse(`touse') ///
            rts(`rts') tolerance(`tolerance') `debug'
    }
    else if "`model'" == "supereff" {
        _dea_ipm_supereff, inputs(`inputs') outputs(`outputs') ///
            efficiency(`efficiency_score') touse(`touse') ///
            rts(`rts') tolerance(`tolerance') orientation(`orientation') ///
            lambda_vars(`lambda_vars') ///
            `reference' `debug'
    }
    else if "`model'" == "multiplier" {
        _dea_ipm_multiplier, inputs(`inputs') outputs(`outputs') ///
            efficiency(`efficiency_score') touse(`touse') ///
            rts(`rts') tolerance(`tolerance') orientation(`orientation') ///
            `saveweights' m(`m') s(`s') ///
            weight_in_vars(`weight_in_vars') ///
            weight_out_vars(`weight_out_vars') ///
            weight_u0_var(`weight_u0_var') ///
            `debug'
    }
    else if "`model'" == "allocative" {
        _dea_ipm_allocative, inputs(`inputs') outputs(`outputs') ///
            efficiency(`efficiency_score') touse(`touse') ///
            rts(`rts') tolerance(`tolerance') orientation(`orientation') ///
            w(`w') p(`p') `debug'
    }
    else if "`model'" == "supersbm" {
        local _sbm_ort_str = cond(`sbmort'==1, "input", cond(`sbmort'==2, "output", ""))
        _dea_ipm_supersbm, inputs(`inputs') outputs(`outputs') ///
            efficiency(`efficiency_score') touse(`touse') ///
            rts(`rts') tolerance(`tolerance') orientation("`_sbm_ort_str'") `debug'
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
        [slack_in_vars(string)] [slack_out_vars(string)] ///
        [stage1_status_var(string)] [stage2_status_var(string)] ///
        [score_source_var(string)] [ref_stage_var(string)] [diagnostics_mat(string)] ///
        [ipm_trace_mat(string)] ///
        [RADIALTolerance(real 0)] [projection_mat(string)] [projection_trace_mat(string)] ///
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

    // Caller-owned slack/status variables survive this wrapper's return.
    local do_stage2 = (`stage' == 2)

    // Call unified envelopment solver
    // Use quietly to suppress LinearProgram warnings during Stage 2
    quietly mata: _dea_ipm_envelopment("`inputs'", "`outputs'", ///
        "`efficiency'", "`lambda_vars'", ///
        "`slack_in_vars'", "`slack_out_vars'", ///
        "`touse'", "`orientation'", "`rts'", ///
        `do_stage2', `tolerance', "`stage1_status_var'", "`stage2_status_var'", "`ipm_trace_mat'", ///
        `radialtolerance', "`projection_mat'", "`projection_trace_mat'", ///
        "`score_source_var'", "`ref_stage_var'", "`diagnostics_mat'")

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
        EFFICIENCY(varname) TOUSE(varname) RTS(string) TOLerance(real) ///
        W(varlist numeric) [DEBUG]
    
    // Cost Efficiency: min w'x s.t. production constraints
    // Cost_Eff = min_cost / observed_cost
    
    if "`debug'" != "" {
        di as text "[Cost] Computing cost efficiency..."
        di as text "[Cost] RTS=`rts'"
        di as text "[Cost] Input prices: `w'"
    }
    
    mata: _dea_ipm_cost("`inputs'", "`outputs'", ///
        "`efficiency'", "`w'", "`touse'", "`rts'", `tolerance')
    
    return add
end

capture program drop _dea_ipm_revenue
program define _dea_ipm_revenue, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) ///
        EFFICIENCY(varname) TOUSE(varname) RTS(string) TOLerance(real) ///
        P(varlist numeric) [DEBUG]
    
    // Revenue Efficiency: max p'y s.t. production constraints
    // Rev_Eff = observed_revenue / max_revenue
    
    if "`debug'" != "" {
        di as text "[Revenue] Computing revenue efficiency..."
        di as text "[Revenue] RTS=`rts'"
        di as text "[Revenue] Output prices: `p'"
    }
    
    mata: _dea_ipm_revenue("`inputs'", "`outputs'", ///
        "`efficiency'", "`p'", "`touse'", "`rts'", `tolerance')
    
    return add
end

capture program drop _dea_ipm_profit
program define _dea_ipm_profit, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) ///
        EFFICIENCY(varname) TOUSE(varname) RTS(string) TOLerance(real) ///
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
        "`efficiency'", "`w'", "`p'", "`touse'", "`rts'", `tolerance')
    
    return add
end

capture program drop _dea_ipm_ddf
program define _dea_ipm_ddf, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) ///
        EFFICIENCY(varname) TOUSE(varname) RTS(string) TOLerance(real) ///
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
        "`efficiency'", "`gx'", "`gy'", "`touse'", "`rts'", `tolerance')

    return add
end

capture program drop _dea_ipm_ddf_undesirable
program define _dea_ipm_ddf_undesirable, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) BAD(string) ///
        EFFICIENCY(varname) TOUSE(varname) RTS(string) TOLerance(real) ///
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
        "`efficiency'", "`direction'", "`touse'", "`rts'", `tolerance')

    return add
end

capture program drop _dea_ipm_sbm
program define _dea_ipm_sbm, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) ///
        EFFICIENCY(varname) TOUSE(varname) RTS(string) TOLerance(real) [DEBUG]
    
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
        "`touse'", "`rts'", `tolerance')
    
    return local slack_in_vars "`slack_in_vars'"
    return local slack_out_vars "`slack_out_vars'"
    return add
end

capture program drop _dea_ipm_additive
program define _dea_ipm_additive, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) ///
        EFFICIENCY(varname) TOUSE(varname) RTS(string) TOLerance(real) [DEBUG]
    
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
        "`touse'", "`rts'", `tolerance')
    
    return local slack_in_vars "`slack_in_vars'"
    return local slack_out_vars "`slack_out_vars'"
    return add
end

capture program drop _dea_ipm_supereff
program define _dea_ipm_supereff, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) ///
        EFFICIENCY(varname) TOUSE(varname) RTS(string) TOLerance(real) ///
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
        "`touse'", "`orientation'", "`rts'", `tolerance')
    
    return add
end

capture program drop _dea_ipm_multiplier
program define _dea_ipm_multiplier, rclass
    version 16.0
    syntax , Inputs(string) Outputs(string) ///
        EFFICIENCY(varname) TOUSE(varname) RTS(string) TOLerance(real) ///
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
            "`touse'", "`orientation'", "`rts'", `tolerance')
    }

    // By strong duality, multiplier efficiency equals the radial
    // envelopment score. Use the envelopment solver for the reported score;
    // multiplier weights, when requested, are retained from the dual LP above.
    mata: _dea_ipm_envelopment("`inputs'", "`outputs'", ///
        "`efficiency'", "", "", "", "`touse'", "`orientation'", "`rts'", 0, `tolerance')
    
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
        EFFICIENCY(varname) TOUSE(varname) RTS(string) TOLerance(real) ///
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
        "`te_score'", "", "", "", "`touse'", "`orientation'", "`rts'", 0, `tolerance')
    
    // Then compute economic efficiency
    tempvar econ_eff
    qui gen double `econ_eff' = . if `touse'
    
    if "`orientation'" == "input" & "`w'" != "" {
        // Cost efficiency
        mata: _dea_ipm_cost("`inputs'", "`outputs'", ///
            "`econ_eff'", "`w'", "`touse'", "`rts'", `tolerance')
    }
    else if "`orientation'" == "output" & "`p'" != "" {
        // Revenue efficiency
        mata: _dea_ipm_revenue("`inputs'", "`outputs'", ///
            "`econ_eff'", "`p'", "`touse'", "`rts'", `tolerance')
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
        EFFICIENCY(varname) TOUSE(varname) RTS(string) TOLerance(real) [ORIentation(string) DEBUG]
    
    // Super-efficiency SBM (Tone, 2002)
    // Excludes evaluated DMU from reference set
    
    if "`debug'" != "" {
        di as text "[SuperSBM] Computing super-efficiency SBM..."
        di as text "[SuperSBM] RTS=`rts'"
    }
    
    mata: _dea_ipm_supersbm("`inputs'", "`outputs'", ///
        "`efficiency'", "`touse'", "`rts'", `tolerance', "`orientation'")
    
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
        [stage1_status_var(string)] [stage2_status_var(string)] ///
        [score_source_var(string)] [ref_stage_var(string)] [diagnostics_mat(string)] [linkage_mat(string)] ///
        [weight_in_vars(string)] [weight_out_vars(string)] [weight_u0_var(string)] ///
        [SBMORT(integer 0)] ///
        [TRACE] [DEBUG]

    // =================================================================
    // Load Simplex Mata functions
    // =================================================================
    // Check both core and undesirable functions; reload if either is missing
    local _simplex_needs_load 0
    capture mata: assert(_dea_simplex_audit_ready()==16 & findexternal("_rsm_duplicate_ratio()")!=NULL & findexternal("_rsm_duplicate_map()")!=NULL & findexternal("_rsm_duplicate_problem()")!=NULL & findexternal("_rsm_duplicate_accept()")!=NULL & findexternal("_rsm_duplicate_retry()")!=NULL & findexternal("_rsm_two_phase_core()")!=NULL & findexternal("_rsm_solve_lp()")!=NULL & findexternal("_rsm_quality_fast()")!=NULL & findexternal("_rsm_affine_quality()")!=NULL & findexternal("_rsm_affine_refine()")!=NULL & findexternal("_rsm_success_fastread()")!=NULL & findexternal("_rsm_success_reset()")!=NULL & findexternal("_rsm_success_read()")!=NULL & findexternal("_rsm_success_counts_read()")!=NULL & findexternal("_rsm_success_check()")!=NULL & findexternal("_rsm_success_apply()")!=NULL & findexternal("_rsm_xgrow()")!=NULL & findexternal("_rsm_xadd()")!=NULL & findexternal("_rsm_xmul()")!=NULL & findexternal("_rsm_xnormalize()")!=NULL & findexternal("_rsm_rowmap()")!=NULL & findexternal("_rsm_row_empty()")!=NULL & findexternal("_rsm_row_accept()")!=NULL & findexternal("_rsm_row_retry()")!=NULL & findexternal("_dea_link_reset()")!=NULL & findexternal("_dea_link_read()")!=NULL & findexternal("_dea_link_store()")!=NULL & findexternal("_dea_link_small()")!=NULL & findexternal("_dea_link_residual()")!=NULL & findexternal("_dea_link_basis()")!=NULL & findexternal("_dea_link_apply()")!=NULL & findexternal("_dea_link_stage2()")!=NULL)
    if _rc != 0 local _simplex_needs_load 1
    capture mata: assert(findexternal("_dea_simplex_envelopment()")!=NULL)
    if _rc != 0 local _simplex_needs_load 1
    capture mata: assert(findexternal("_dea_simplex_ddf_undesirable()")!=NULL)
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
    
    capture mata: assert(_dea_simplex_audit_ready()==16 & findexternal("_rsm_duplicate_ratio()")!=NULL & findexternal("_rsm_duplicate_map()")!=NULL & findexternal("_rsm_duplicate_problem()")!=NULL & findexternal("_rsm_duplicate_accept()")!=NULL & findexternal("_rsm_duplicate_retry()")!=NULL & findexternal("_rsm_two_phase_core()")!=NULL & findexternal("_rsm_solve_lp()")!=NULL & findexternal("_rsm_quality_fast()")!=NULL & findexternal("_rsm_affine_quality()")!=NULL & findexternal("_rsm_affine_refine()")!=NULL & findexternal("_rsm_success_fastread()")!=NULL & findexternal("_rsm_success_reset()")!=NULL & findexternal("_rsm_success_read()")!=NULL & findexternal("_rsm_success_counts_read()")!=NULL & findexternal("_rsm_success_check()")!=NULL & findexternal("_rsm_success_apply()")!=NULL & findexternal("_rsm_xgrow()")!=NULL & findexternal("_rsm_xadd()")!=NULL & findexternal("_rsm_xmul()")!=NULL & findexternal("_rsm_xnormalize()")!=NULL & findexternal("_rsm_rowmap()")!=NULL & findexternal("_rsm_row_empty()")!=NULL & findexternal("_rsm_row_accept()")!=NULL & findexternal("_rsm_row_retry()")!=NULL & findexternal("_dea_link_reset()")!=NULL & findexternal("_dea_link_read()")!=NULL & findexternal("_dea_link_store()")!=NULL & findexternal("_dea_link_small()")!=NULL & findexternal("_dea_link_residual()")!=NULL & findexternal("_dea_link_basis()")!=NULL & findexternal("_dea_link_apply()")!=NULL & findexternal("_dea_link_stage2()")!=NULL)
    if _rc {
        di as error "Incompatible simplex library; reinstall the complete matched DEA package."
        exit 498
    }

    // =================================================================
    // Create efficiency variable
    // =================================================================
    foreach _drop_output in `generate' {
        capture confirm variable `_drop_output', exact
        if !_rc quietly drop `_drop_output'
    }
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
            `tolerance', "`lambda_vars'", "`slack_in_vars'", "`slack_out_vars'", `stage', ///
            "`stage1_status_var'", "`stage2_status_var'", "`score_source_var'", "`ref_stage_var'", "`diagnostics_mat'", "`linkage_mat'")
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
                "`generate'", "`touse'", "`orientation'", "`rts'", `tolerance', "", "", "", `stage')
        }
    }
    else if "`model'" == "cost" {
        // ---------------------------------------------------------
        // Cost efficiency
        // ---------------------------------------------------------
        quietly mata: _dea_simplex_cost("`inputs'", "`outputs'", ///
            "`generate'", "`w'", "`touse'", "`rts'", `tolerance')
    }
    else if "`model'" == "revenue" {
        // ---------------------------------------------------------
        // Revenue efficiency
        // ---------------------------------------------------------
        quietly mata: _dea_simplex_revenue("`inputs'", "`outputs'", ///
            "`generate'", "`p'", "`touse'", "`rts'", `tolerance')
    }
    else if "`model'" == "profit" {
        // ---------------------------------------------------------
        // Profit efficiency
        // ---------------------------------------------------------
        quietly mata: _dea_simplex_profit("`inputs'", "`outputs'", ///
            "`generate'", "`w'", "`p'", "`touse'", "`rts'", `tolerance')
    }
    else if "`model'" == "ddf" {
        // ---------------------------------------------------------
        // Directional distance function
        // ---------------------------------------------------------
        local gx_arg = cond("`gx'" != "", "`gx'", "")
        local gy_arg = cond("`gy'" != "", "`gy'", "")
        quietly mata: _dea_simplex_ddf("`inputs'", "`outputs'", ///
            "`generate'", "`gx_arg'", "`gy_arg'", "`touse'", "`rts'", `tolerance')
    }
    else if "`model'" == "undesirable" {
        // ---------------------------------------------------------
        // DDF with undesirable (bad) outputs — Wang et al. (2022)
        // ---------------------------------------------------------
        quietly mata: _dea_simplex_ddf_undesirable("`inputs'", "`outputs'", ///
            "`bad'", "`generate'", "`direction'", "`touse'", "`rts'", `tolerance')
    }
    else if "`model'" == "sbm" {
        // ---------------------------------------------------------
        // Slacks-based measure
        // sbm_ort: 0=non-oriented, 1=input-oriented, 2=output-oriented
        // ---------------------------------------------------------
        local _sbm_ort_str = cond(`sbmort'==1, "input", cond(`sbmort'==2, "output", ""))
        quietly mata: _dea_simplex_sbm("`inputs'", "`outputs'", ///
            "`generate'", "", "", "`touse'", "`rts'", "`_sbm_ort_str'", `tolerance')
    }
    else if "`model'" == "additive" {
        // ---------------------------------------------------------
        // Additive model
        // ---------------------------------------------------------
        quietly mata: _dea_simplex_additive("`inputs'", "`outputs'", ///
            "`generate'", "", "", "`touse'", "`rts'", `tolerance')
    }
    else if "`model'" == "supereff" {
        // ---------------------------------------------------------
        // Radial super-efficiency (with optional lambda extraction)
        // Lambda variables are pre-created by _dea_model_dispatch
        // ---------------------------------------------------------
        quietly mata: _dea_simplex_supereff("`inputs'", "`outputs'", ///
            "`generate'", "`lambda_vars'", ///
            "`touse'", "`orientation'", "`rts'", `tolerance')
    }
    else if "`model'" == "supersbm" {
        // ---------------------------------------------------------
        // Super-efficiency SBM
        // ---------------------------------------------------------
        local _sbm_ort_str = cond(`sbmort'==1, "input", cond(`sbmort'==2, "output", ""))
        quietly mata: _dea_simplex_supersbm("`inputs'", "`outputs'", ///
            "`generate'", "`touse'", "`rts'", `tolerance', "`_sbm_ort_str'")
    }
    else if "`model'" == "allocative" {
        // ---------------------------------------------------------
        // Allocative efficiency (composite: TE + Cost/Revenue)
        // ---------------------------------------------------------
        local w_arg = cond("`w'" != "", "`w'", "")
        local p_arg = cond("`p'" != "", "`p'", "")
        quietly mata: _dea_simplex_allocative("`inputs'", "`outputs'", ///
            "`generate'", "`touse'", "`orientation'", "`rts'", ///
            "`w_arg'", "`p_arg'", `tolerance')
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
// Historical diagnostic, not the Cooper slack-decomposition algorithm:
//   standard input VRS theta divided by the input/output-equality theta.
//   A gap alone neither establishes nor excludes input congestion.
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
    di as text "Congestion diagnostic — legacy equality-restricted radial ratio"
    di as text "A ratio below 1 is not sufficient to identify input congestion."
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
    // Step 2: historical input/output-equality LP → theta_C
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
        foreach _drop_output in `v' {
            capture confirm variable `_drop_output', exact
            if !_rc quietly drop `_drop_output'
        }
        quietly gen double `v' = . if `touse'
    }

    quietly replace `v_crs'  = `_eff_crs'  if `touse'
    quietly replace `v_vrs'  = `_eff_vrs'  if `touse'
    quietly replace `v_ceff' = `_eff_cong' if `touse'

    // Legacy ratio; one indicates no radial equality-restriction gap.
    quietly replace `v_cong' = `v_vrs' / `v_ceff' ///
        if `touse' & !missing(`v_vrs', `v_ceff') & `v_ceff' > `tolerance'

    // Scale efficiency = θ_CRS / θ_VRS
    quietly replace `v_se' = `v_crs' / `v_vrs' ///
        if `touse' & !missing(`v_crs', `v_vrs') & `v_vrs' > `tolerance'

    // Primary generate retains the historical ratio, not a congestion certificate.
    foreach _drop_output in `generate' {
        capture confirm variable `_drop_output', exact
        if !_rc quietly drop `_drop_output'
    }
    quietly gen double `generate' = `v_cong' if `touse'
    label var `generate' "Legacy equality-restricted ratio (theta_VRS/theta_C)"

    // ---------------------------------------------------------------
    // Summary statistics
    // ---------------------------------------------------------------
    qui count if `touse'
    local N = r(N)
    qui summ `v_vrs' if `touse'
    local mean_vrs = r(mean)
    qui count if `v_cong' < (1 - `tolerance') & `touse'
    local n_cong = r(N)

    // Availability of the main ratio and all requested component results.
    // A CRS-only failure does not erase a usable congestion measure.
    quietly count if `touse' & !missing(`generate') & abs(`generate'-1)<=`tolerance'
    local n_eff = r(N)
    quietly count if `touse' & missing(`generate')
    local n_failed_primary = r(N)
    foreach _component in crs vrs ceff se {
        quietly count if `touse' & missing(`v_`_component'')
        local n_failed_`_component' = r(N)
    }
    quietly count if `touse' & missing(`generate', `v_crs', `v_vrs', `v_ceff', `v_se')
    local n_failed = r(N)
    quietly summarize `generate' if `touse'
    local mean_eff = r(mean)
    local min_eff = r(min)
    local max_eff = r(max)

    di as text "DMUs evaluated:  " as result `N'
    di as text "DMUs with ratio below 1: " as result `n_cong' ///
        as text " (legacy diagnostic; not a congestion certificate)"
    di as text "Mean VRS eff.:   " as result %8.6f `mean_vrs'
    if `n_failed' {
        di as text "Warning: `n_failed' DMU(s) have unavailable requested results."
        di as text "Missing congestion measures: `n_failed_primary'; missing scale efficiencies: `n_failed_se'."
        di as text "Missing CRS/VRS/congestion LP scores: `n_failed_crs'/`n_failed_vrs'/`n_failed_ceff'."
        di as text "Available components are retained; missing results do not imply no congestion."
    }
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
    // Standard dea returns use finite main scores, independently of auxiliary availability.
    return scalar mean_efficiency = `mean_eff'
    return scalar n_efficient     = `n_eff'
    return scalar min_eff         = `min_eff'
    return scalar max_eff         = `max_eff'
    return scalar n_failed        = `n_failed'
    return scalar n_failed_primary = `n_failed_primary'
    return scalar n_failed_crs    = `n_failed_crs'
    return scalar n_failed_vrs    = `n_failed_vrs'
    return scalar n_failed_cong   = `n_failed_ceff'
    return scalar n_failed_se     = `n_failed_se'
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
//        Bλ + sᵇ = b₀  (bad outputs: slack balance, B lambda <= b0)
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
    foreach _drop_output in `generate' {
        capture confirm variable `_drop_output', exact
        if !_rc quietly drop `_drop_output'
    }
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
    qui count if `touse' & !missing(`generate') & abs(`generate'-1)<=`tolerance'
    local n_eff = r(N)
    quietly count if `touse' & missing(`generate')
    local n_failed = r(N)

    di as text "DMUs evaluated:  " as result `N'
    di as text "Efficient DMUs:  " as result `n_eff'
    di as text "Mean efficiency: " as result %8.6f `mean_eff'
    if `n_failed' {
        di as text "Warning: `n_failed' DMU(s) have missing SBM-undesirable scores."
        di as text "Missing scores are excluded from efficiency counts and score summaries."
    }
    di as text "{hline 65}"

    // Standard dea returns
    return scalar mean_efficiency = `mean_eff'
    return scalar n_efficient     = `n_eff'
    return scalar N               = `N'
    qui summ `generate' if `touse'
    return scalar min_eff  = r(min)
    return scalar max_eff  = r(max)
    return scalar n_failed = `n_failed'
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
    real scalar rc, n, m, s, k
    real colvector scores
    real rowvector x0, y0, c, lb, ub, xs, ys
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

    xs=colmax(abs(X));xs=xs:+(xs:==0)
    ys=colmax(abs(Y));ys=ys:+(ys:==0)
    X=X:/xs;Y=Y:/ys
    scores = J(n, 1, .)

    for (k = 1; k <= n; k++) {
        x0 = X[k, .]
        y0 = Y[k, .]

        // Variables: [λ₁...λₙ, θ]  — n+1 variables
        // Objective: min θ
        c  = (J(1, n, 0), 1)

        // Radial input equality: Xλ = θ·x₀
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
        rc = _dea_nonte_optimize(S,tol)

        if (rc == 0) {
            sol = S.parameters()
            scores[k] = sol[n + 1]
        }
        else {
            scores[k] = .   // An unavailable LP result is not evidence of no congestion.
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
    real scalar rc, n, m, s, b, k
    real colvector scores, obs_idx
    real rowvector x0, y0, b0, xscale, yscale, bscale
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

    // Normalize columns and the corresponding slack units together. Ratios,
    // t and peer weights are unchanged; tiny positive values remain positive.
    xscale=colmax(abs(X)); xscale=xscale:+(xscale:==0)
    yscale=colmax(abs(Y)); yscale=yscale:+(yscale:==0)
    bscale=colmax(abs(B)); bscale=bscale:+(bscale:==0)
    X=X:/xscale; Y=Y:/yscale; B=B:/bscale

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
            if (x0[i] > 0) c[n + i] = -(1/m) / x0[i]
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

        // Equality 3: Bz + Sᵇ = t·b₀  (bad output: separable reduction slack)
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
            if (y0[r] > 0) Aeq4[n + m + r] = (1/(s+b)) / y0[r]
        }
        for (r = 1; r <= b; r++) {
            if (b0[r] > 0) Aeq4[n + m + s + r] = (1/(s+b)) / b0[r]
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

        // Transformed nonconstant RTS: z=t*lambda; t is not fixed to one.
        A_ineq  = J(0, nv, 0)
        b_ineq_v = J(0, 1, 0)
        if (rts=="drs" | rts=="irs") {
            Aeq5=J(1,nv,0);Aeq5[1..n]=J(1,n,1);Aeq5[nv]=-1
            A_ineq=(rts=="drs" ? Aeq5 : -Aeq5)
            b_ineq_v=0
        }

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
        rc = _dea_nonte_optimize(S,tol)

        if (rc == 0) {
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


// Persistent reference weights: safe names and sample-aligned storage.
mata:
void _dea_saveref_plan(string scalar dmu, string scalar touse)
{
    real colvector obs
    string colvector ids
    string rowvector names
    string scalar base, name, suffix
    real scalar j, serial
    obs = selectindex(st_data(., touse))
    if (st_isstrvar(dmu)) ids = st_sdata(obs, dmu)
    else ids = strtrim(strofreal(st_data(obs, dmu), "%21.0g"))
    names = J(1, 0, "")
    for (j=1; j<=rows(obs); j++) {
        // Stata 16 strtoname() can retain a backtick, which is not a valid
        // variable-name character and is unsafe when passed through a macro.
        base = strtoname("_ref_" + subinstr(ids[j], char(96), "_"))
        name = base
        serial = 1
        while (any(names :== name)) {
            serial++
            suffix = "_" + strofreal(serial)
            name = substr(base, 1, 32-strlen(suffix)) + suffix
        }
        names = names, name
    }
    st_local("_ref_names", invtokens(names))
}

void _dea_saveref_store(string scalar dmu, string scalar touse,
                       string scalar refvars, string scalar matname)
{
    real colvector obs
    real rowvector vi
    real matrix L
    string colvector ids
    string rowvector names
    real scalar j
    obs = selectindex(st_data(., touse))
    names = tokens(refvars)
    L = st_matrix(matname)
    if (rows(L)!=rows(obs) | cols(L)!=cols(names)) _error(503)
    if (st_isstrvar(dmu)) ids = st_sdata(obs, dmu)
    else ids = strtrim(strofreal(st_data(obs, dmu), "%21.0g"))
    vi = J(1, cols(names), .)
    for (j=1; j<=cols(names); j++) {
        vi[j] = st_addvar("double", names[j])
        st_varlabel(vi[j], usubstr("Reference weight: peer " + strofreal(j) +
                    ", DMU " + ids[j], 1, 80))
        st_global(names[j]+"[dea_peer_id]", ids[j])
        st_global(names[j]+"[dea_peer_index]", strofreal(j))
        st_global(names[j]+"[dea_peer_dmu]", dmu)
    }
    st_store(obs, vi, L)
}
end
