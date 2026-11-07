*! version 2.3.0 06nov2026
*! Return delegated failure counts and documented summary aliases.
*! 2026-11-01: use shared Form A/B parser; preserve if/in expressions.
*! 2026-11-01: preflight output names; explicit replacevars; protect model data.
*! dea_congestion.ado — Explicit Cooper decomposition / legacy radial diagnostic
*!   Wrapper: delegates to dea model(congestion)
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
*! Cooper, W. W., Deng, H., Seiford, L. M., & Zhu, J. (2004).
*!   Congestion: Its identification and management with DEA.
*!   In W. W. Cooper, L. M. Seiford, and J. Zhu (Eds.),
*!   Handbook on Data Envelopment Analysis, pp. 177-201.
*!   Boston: Kluwer Academic Publishers.
*!
*! Cooper, W. W., Seiford, L. M., & Tone, K. (2007). Data Envelopment
*!   Analysis: A Comprehensive Text with Models, Applications, References
*!   and DEA-Solver Software (2nd ed.). New York: Springer.

/**
 * CONGESTION DEA
 * ==============
 *
 * congmethod(cooper) selects the new hybrid VRS worker; see the help equations.
 * The default retains the historical equality-restricted radial ratio. This is not the
 * Cooper slack decomposition and does not identify congestion by itself.
 * Congestion occurs when reducing an input can actually increase output
 * (perverse input-output relationship).
 *
 * Method:
 *   Step 1: Standard VRS DEA (strong disposability) → θ_VRS
 *   Step 2: Congestion model (weak input disposability; input and output
 *           equality constraints) → θ_C
 *
 *   Congestion measure = θ_VRS / θ_C
 *     < 1 : equality restriction gap (θ_C > θ_VRS)
 *     = 1 : no radial gap; does not rule out nonradial congestion
 *
 *   Scale efficiency = θ_CRS / θ_VRS
 *
 * Output variables (stub prefix, default "_cong"):
 *   _cong_crs        CRS technical efficiency
 *   _cong_vrs        VRS technical efficiency
 *   _cong_cong_eff   Congestion-model efficiency (θ_C)
 *   _cong_congestion Congestion measure (θ_VRS / θ_C)
 *   _cong_se         Scale efficiency (θ_CRS / θ_VRS)
 *
 * Note: This command is a convenience wrapper around:
 *   dea inputs = outputs, model(congestion) generate(prefix) ...
 */

capture program drop dea_congestion
program define dea_congestion, rclass
    version 16.0

    // =========================================================================
    // Unified Syntax Parser
    //   Form A: dea_congestion x1 x2 = y1 y2, [options]
    //   Form B: dea_congestion [if] [in], inputs(x1 x2) outputs(y1 y2) [opts]
    // =========================================================================
    _dea_parse_io `0'
    local 0 `"`r(args)'"'

    syntax [if] [in], ///
        Inputs(varlist numeric) ///
        Outputs(varlist numeric) ///
        [ ///
        Generate(name) ///
        REPLACEVARS    ///
        dmu(varname) ///
        engine(string) ///
        TOLerance(real 1e-7) ///
        CONGMethod(string) CONGWeight(string) CONGMATrix CONGRanges ///
        CONGPolicy(string) CONGAnchor(string) CONGREFerence(string) CONGSelectlimit(string) ///
        SAVing(string) REPLACE DEBUG ///
        ]

    // Defaults
    if "`generate'" == "" local generate "_cong"
    if "`engine'"   == "" local engine "ipm"

    // Build options for dea
    local dmu_opt ""
    if "`dmu'" != "" local dmu_opt "dmu(`dmu')"

    if lower("`congmethod'")=="cooper" {
        dea `inputs' = `outputs' `if' `in', model(congestion) ///
            congmethod(cooper) congweight(`congweight') `congmatrix' `congranges' ///
            congpolicy(`congpolicy') conganchor(`conganchor') congreference(`congreference') congselectlimit(`congselectlimit') ///
            generate(`generate') engine(`engine') tolerance(`tolerance') ///
            `dmu_opt' `replacevars' saving(`saving') `replace' `debug'
        return add
        return local cmd "dea_congestion"
        exit
    }

    if "`congpolicy'`conganchor'`congreference'`congselectlimit'"!="" exit 198

    // Mark sample
    marksample touse
    markout `touse' `inputs' `outputs'
    if "`dmu'" != "" markout `touse' `dmu', strok

    qui count if `touse'
    if r(N) < 2 {
        di as error "dea_congestion: at least 2 observations required"
        error 2001
    }

    // =========================================================================
    // Delegate to dea model(congestion)
    // =========================================================================
    dea `inputs' = `outputs' if `touse', ///
        model(congestion) generate(`generate') ///
        engine(`engine') tolerance(`tolerance') `dmu_opt' `replacevars' ///
        congmethod(`congmethod') congweight(`congweight') `congmatrix' `congranges' ///
        saving(`saving') `replace' `debug'

    // Preserve the documented aliases after the core's rclass work.
    local _N_eff = r(N_eff)
    local _mean_eff = r(mean_eff)
    return add
    return scalar n_efficient = `_N_eff'
    return scalar mean_efficiency = `_mean_eff'
end
