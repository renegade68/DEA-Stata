*! version 1.1.0  03nov2026
*! 2026-11-01: accept the minimum option abbreviations advertised in help.
*! 2026-11-01: use shared Form A/B parser; preserve if/in expressions.
*! 2026-11-01: preflight output names; explicit replacevars; protect model data.
*! dea_sbmundesirable.ado — SBM with Undesirable Outputs (Tone, 2004)
*!   Wrapper: delegates to dea model(sbmundesirable)
*!
*! AUTHORS:
*! --------
*!   Choonjoo Lee (corresponding author)
*!     Korea National Defense University
*!     Seoul, Korea
*!     bloom.rampike@gmail.com
*!
*! Changelog:
*!   1.0.1  02may2026  Add wrapper-level return aliases for help compatibility
*!
*! REFERENCES:
*! -----------
*!   Tone, K. (2004). Dealing with undesirable outputs in DEA:
*!     A Slacks-Based Measure (SBM) approach. GRIPS Discussion Paper
*!     I-2003-0005. doi:10.24545/00000955
*!
*!   Tone, K. (2001). A slacks-based measure of efficiency in data
*!     envelopment analysis. European Journal of Operational Research,
*!     130(3), 498-509. doi:10.1016/S0377-2217(99)00407-5
*!
*!   Chambers, R. G., Chung, Y., & Fare, R. (1998). Profit, directional
*!     distance functions, and Nerlovian efficiency. Journal of Optimization
*!     Theory and Applications, 98(2), 351-364.
*!     doi:10.1023/A:1022637501082
*!
*!   Cooper, W. W., Seiford, L. M., & Tone, K. (2007). Data Envelopment
*!     Analysis: A Comprehensive Text with Models, Applications, References
*!     and DEA-Solver Software (2nd ed.). New York: Springer.
*!     doi:10.1007/978-0-387-45283-8

/**
 * SBM-UNDESIRABLE DEA
 * ===================
 *
 * Computes efficiency for production processes that generate both
 * desirable (good) and undesirable (bad) outputs, using the
 * Slacks-Based Measure (SBM) framework of Tone (2004).
 *
 * Model (fractional form):
 *   min  (1 - (1/m)*Σ(s⁻ᵢ/xᵢ₀)) / (1 + (1/(s+b))*(Σ(s⁺ᵣ/yᵣ₀) + Σ(sᵇₖ/bₖ₀)))
 *   s.t. Xλ + s⁻ = x₀      (input slack)
 *        Yλ - s⁺ = y₀      (good output slack, ≥ constraint)
 *        Bλ + sᵇ = b₀      (separable bad-output reduction slack)
 *        Σλ = 1 (VRS), <=1 (NIRS/DRS), >=1 (NDRS/IRS), or unrestricted (CRS)
 *        λ, s⁻, s⁺, sᵇ ≥ 0
 *
 * Solved via Charnes-Cooper linearization.
 *
 * Note: This command is a convenience wrapper around:
 *   dea inputs = outputs, model(sbmundesirable) bad(badvars) ...
 */

capture program drop dea_sbmundesirable
program define dea_sbmundesirable, rclass
    version 16.0

    // =========================================================================
    // Unified Syntax Parser
    //   Form A: dea_sbmundesirable x1 x2 = y1 y2, bad(b1) [options]
    //   Form B: dea_sbmundesirable [if] [in], inputs(x1 x2) outputs(y1 y2)
    //             bad(b1) [options]
    // =========================================================================
    _dea_parse_io `0'
    local 0 `"`r(args)'"'

    syntax [if] [in], ///
        Inputs(varlist numeric) ///
        Outputs(varlist numeric) ///
        BAD(varlist numeric) ///
        [ ///
        Generate(name) ///
        REPLACEVARS    ///
        dmu(varname) ///
        Rts(string) ///
        TOLerance(real 1e-7) ///
        ]

    // Defaults
    if "`generate'" == "" local generate "sbmu_eff"
    if "`rts'"      == "" local rts "crs"

    // Build options for dea
    local dmu_opt ""
    if "`dmu'" != "" local dmu_opt "dmu(`dmu')"

    // Mark sample
    marksample touse
    markout `touse' `inputs' `outputs' `bad'
    if "`dmu'" != "" markout `touse' `dmu', strok

    qui count if `touse'
    if r(N) < 2 {
        di as error "dea_sbmundesirable: at least 2 observations required"
        error 2001
    }

    // =========================================================================
    // Delegate to dea model(sbmundesirable)
    // =========================================================================
    dea `inputs' = `outputs' if `touse', ///
        model(sbmundesirable) bad(`bad') ///
        rts(`rts') generate(`generate') ///
        tolerance(`tolerance') `dmu_opt' `replacevars'

    local _N        = r(N)
    local _N_eff    = r(N_eff)
    local _mean_eff = r(mean_eff)
    local _min_eff  = r(min_eff)
    local _max_eff  = r(max_eff)

    // Return results
    return add
    return scalar n_efficient     = `_N_eff'
    return scalar mean_efficiency = `_mean_eff'
    return scalar N               = `_N'
    return scalar min_eff         = `_min_eff'
    return scalar max_eff         = `_max_eff'
end
