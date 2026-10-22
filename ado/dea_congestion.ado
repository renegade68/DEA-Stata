*! version 2.0.0  26apr2026
*! dea_congestion.ado — Congestion DEA (Cooper et al. 2004)
*!   Wrapper: delegates to dea model(congestion)
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
 * Detects input congestion using the Cooper et al. (2004) two-step procedure.
 * Congestion occurs when reducing an input can actually increase output
 * (perverse input-output relationship).
 *
 * Method:
 *   Step 1: Standard VRS DEA (strong disposability) → θ_VRS
 *   Step 2: Congestion model (weak input disposability; input and output
 *           equality constraints) → θ_C
 *
 *   Congestion measure = θ_VRS / θ_C
 *     < 1 : congestion present (θ_C > θ_VRS)
 *     = 1 : no congestion
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
            local _io_outputs = trim(substr(`"`0'"', 1, `_cpos' - 1))
            local _post = substr(`"`0'"', `_cpos' + 1, .)
            local 0 `", inputs(`_io_inputs') outputs(`_io_outputs') `_post'"'
        }
        else {
            local _io_outputs = trim(`"`0'"')
            local 0 `", inputs(`_io_inputs') outputs(`_io_outputs')"'
        }
    }

    syntax [if] [in], ///
        Inputs(varlist numeric) ///
        Outputs(varlist numeric) ///
        [ ///
        Generate(name) ///
        dmu(varname) ///
        engine(string) ///
        TOLerance(real 1e-7) ///
        ]

    // Defaults
    if "`generate'" == "" local generate "_cong"
    if "`engine'"   == "" local engine "ipm"

    // Build options for dea
    local dmu_opt ""
    if "`dmu'" != "" local dmu_opt "dmu(`dmu')"

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
        engine(`engine') tolerance(`tolerance') `dmu_opt'

    // Return results
    return add
end
