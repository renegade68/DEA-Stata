*! version 1.0.0  01nov2026
*! Compatibility entry point for the former IPM DEA command (not generic LP)
program define dea_lp, rclass
    version 16.0
    syntax anything(name=spec equalok) [if] [in] [using/] ///
        [, Rts(string) ENGine(string) METhod(string) *]

    // Preserve the historical default technology and solver. Explicit engine()
    // or method() uses the current dea precedence rules.
    if "`rts'" == "" local rts "vrs"
    if "`engine'`method'" == "" local engine "ipm"
    local solver
    if "`engine'" != "" local solver "engine(`engine')"
    if "`method'" != "" local solver "`solver' method(`method')"
    local source
    if `"`using'"' != "" local source `"using "`using'""'
    dea `spec' `if' `in' `source', rts(`rts') `solver' `options'
    return add
    return local cmd "dea_lp"
end
