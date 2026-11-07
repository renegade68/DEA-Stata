*! version 1.0.0  01nov2026
*! Validate the complete public output plan before any existing data are changed.
program define _dea_output_check, rclass
    version 16.0
    syntax , NAMES(string) [PROTECT(varlist) REPLACEVARS]

    local seen ""
    local existing ""
    foreach v of local names {
        confirm name `v'
        if "`v'" == "_all" {
            di as error "_all cannot be used as an output variable name"
            exit 7
        }
        if `: list v in seen' {
            di as error "output variable `v' is requested more than once; choose distinct output names"
            exit 110
        }
        if `: list v in protect' {
            di as error "output variable `v' is also a model or identifier variable; choose another output name"
            exit 110
        }
        // Exact lookup: a new name must not match an abbreviation of another variable.
        capture confirm variable `v', exact
        if !_rc {
            if "`replacevars'" == "" {
                di as error "variable `v' already exists; specify replacevars or choose another output name"
                exit 110
            }
            local existing `existing' `v'
        }
        else {
            confirm new variable `v'
        }
        local seen `seen' `v'
    }
    return local existing "`existing'"
end
