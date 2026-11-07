*! version 2.0.5  05nov2026
*! Dataset-table LP/MILP interface using the shared continuous LP contract
program define lp, rclass
    version 16.0
    syntax anything(name=variables) [if] [in] [using/] ///
        [, REL(name) RHS(name) MIN MAX INTVars(string) ROWNAMES(name) ///
        ENGine(string) SCALING(string) TOLerance(string) TOL1(string) TOL2(string) ///
        LB(name) UB(name) MAXNodes(integer 10000) INTTol(real 1e-7) ///
        SAVing(string) REPLACE TRACE]

    if "`min'`max'" != "min" & "`min'`max'" != "max" {
        di as error "specify exactly one of min or max"
        exit 198
    }
    if "`engine'" == "" local engine "simplex"
    local engine = lower("`engine'")
    if !inlist("`engine'", "simplex", "ipm", "auto") {
        di as error "engine() must be simplex, ipm, or auto"
        exit 198
    }
    if "`scaling'" == "" local scaling "none"
    local scaling = lower("`scaling'")
    if !inlist("`scaling'", "none", "row") {
        di as error "scaling() must be none or row"
        exit 198
    }
    // The modern solvers have one LP tolerance, not the old pair of pivot and
    // matrix-inversion tolerances. Accept either old spelling without silently
    // choosing between conflicting explicit requests.
    local tol 1e-8
    local seen 0
    foreach opt in tolerance tol1 tol2 {
        if "``opt''" != "" {
            capture confirm number ``opt''
            if _rc {
                di as error "`opt'() must be positive and finite"
                exit 198
            }
            if missing(``opt'') | ``opt'' <= 0 {
                di as error "`opt'() must be positive and finite"
                exit 198
            }
            if `seen' & ``opt'' != `tol' {
                di as error "conflicting tolerances; specify one tolerance() for the modern solver"
                exit 198
            }
            local tol = ``opt''
            local seen 1
        }
    }
    if missing(`inttol') | `inttol' <= 0 | missing(`maxnodes') | `maxnodes' < 1 {
        di as error "inttol() must be positive and finite; maxnodes() must be a positive integer"
        exit 198
    }
    if `"`saving'"' != "" {
        if lower(substr(`"`saving'"', -4, 4)) != ".dta" local saving `"`saving'.dta"'
        if "`replace'" == "" confirm new file `"`saving'"'
    }
    if "`rel'" == "" local rel "rel"
    if "`rhs'" == "" local rhs "rhs"

    // Automatic restoration on both normal and error exits includes using-data.
    preserve
    if `"`using'"' != "" quietly use `"`using'"', clear
    unab variables : `variables'
    confirm numeric variable `variables' `rhs'
    confirm string variable `rel'
    local n : word count `variables'
    local duplicate : list posof "`rhs'" in variables
    if `duplicate' {
        di as error "rhs() must not be a decision-variable column"
        exit 198
    }
    local integers
    if "`intvars'" != "" {
        unab intvars : `intvars'
        local intvars : list uniq intvars
        foreach v of local intvars {
            local j : list posof "`v'" in variables
            if !`j' {
                di as error "intvars() must be a subset of the decision-variable columns"
                exit 198
            }
            local integers `integers' `j'
        }
    }
    marksample touse, novarlist
    quietly keep if `touse'
    if _N == 0 error 2000
    foreach v of local variables {
        quietly count if missing(`v')
        if r(N) {
            di as error "selected objective and constraint coefficients must be finite"
            exit 459
        }
    }
    quietly count if missing(`rhs') & _n > 1
    if r(N) {
        di as error "selected constraint right-hand sides must be finite"
        exit 459
    }
    // The objective-row RHS and relation are placeholders, as in legacy lp.
    quietly replace `rhs' = 0 in 1
    tempname problem C AE BE AU BU row value sol report slack
    local rowopt
    if "`rownames'" != "" local rowopt "rownames(`rownames')"
    quietly mkmat `variables' `rhs', matrix(`problem') `rowopt'
    matrix `C' = `problem'[1,1..`n']
    local m = _N - 1
    local nslacks 0
    local slacknames
    if `m' > 0 {
        forvalues i = 2/`=_N' {
            local op = strtrim(`rel'[`i'])
            if !inlist("`op'", "=", "<=", "<", ">=", ">") {
                di as error "constraint row `i' has an invalid relation"
                exit 198
            }
            matrix `row' = `problem'[`i',1..`n']
            matrix `value' = `problem'[`i',`=`n'+1']
            if "`op'" == "=" {
                matrix `AE' = nullmat(`AE') \ `row'
                matrix `BE' = nullmat(`BE') \ `value'
            }
            else {
                if inlist("`op'", ">", ">=") {
                    matrix `row' = -`row'
                    matrix `value' = -`value'
                }
                matrix `AU' = nullmat(`AU') \ `row'
                matrix `BU' = nullmat(`BU') \ `value'
                local ++nslacks
                local slacknames `slacknames' slack`nslacks'
            }
        }
    }
    local constraints
    capture confirm matrix `AE'
    if !_rc local constraints "aeq(`AE') beq(`BE')"
    capture confirm matrix `AU'
    if !_rc local constraints "`constraints' aub(`AU') bub(`BU')"
    local bounds
    if "`lb'" != "" local bounds "lb(`lb')"
    if "`ub'" != "" local bounds "`bounds' ub(`ub')"
    if "`integers'" == "" {
        // Restore an adapter partially cleared after an earlier public call.
        capture mata: assert(_dea_lpcontract_ready()==9 & findexternal("_rsm_public_success()")!=NULL & findexternal("_dea_lpcontract_simplex()")!=NULL & findexternal("_dea_lpcontract_dispatch()")!=NULL)
        if _rc {
            quietly findfile _dea_lpcontract.ado
            quietly run "`r(fn)'"
        }
        capture mata: assert(_dea_lpcontract_ready()==9 & findexternal("_rsm_public_success()")!=NULL & findexternal("_dea_lpcontract_simplex()")!=NULL & findexternal("_dea_lpcontract_dispatch()")!=NULL)
        if _rc {
            di as error "Incompatible LP adapter; reinstall the complete matched DEA package."
            exit 498
        }
        _dea_lpcontract, c(`C') `constraints' `bounds' sense(`min'`max') ///
            engine(`engine') scaling(`scaling') tolerance(`tol') `trace'
    }
    else {
        milp_bb, c(`C') `constraints' `bounds' sense(`min'`max') ///
            integer(`integers') lpengine(`engine') lpscaling(`scaling') ///
            lptol(`tol') inttol(`inttol') maxnodes(`maxnodes') `trace'
    }
    local status = r(status)
    local objective = r(objective)
    local status_text "`r(status_text)'"
    matrix `sol' = r(x)
    return add
    matrix colnames `sol' = `variables'
    matrix `report' = (`objective'), `sol'
    if `nslacks' {
        matrix `slack' = (`BU' - `AU'*`sol'')'
        matrix `report' = `report', `slack'
    }
    matrix colnames `report' = z `variables' `slacknames'
    matrix rownames `report' = opt_val
    di as text "LP/MILP: `status_text'"
    matrix list `report', noheader
    if `"`saving'"' != "" {
        quietly clear
        quietly set obs `n'
        quietly generate str32 parameter = ""
        quietly generate double value = .
        forvalues j = 1/`n' {
            local v : word `j' of `variables'
            quietly replace parameter = "`v'" in `j'
            quietly replace value = `sol'[1,`j'] in `j'
        }
        quietly generate double objective = `objective'
        quietly generate int status = `status'
        quietly save `"`saving'"', `replace'
    }
    return matrix x = `sol'
    return matrix lprslt = `report'
    return matrix problem = `problem'
    return scalar nvars = `n'
    return scalar nconstraints = `m'
    return scalar nslacks = `nslacks'
    return scalar tolerance = `tol'
    return local cmd "lp"
    return local variables "`variables'"
    return local intvars "`intvars'"
    return local problem_type = cond("`integers'"=="", "LP", "MILP")
end
