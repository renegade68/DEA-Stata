*! version 1.0.0  01nov2026
*! Normalize Form A without treating equality in if expressions as an IO separator.
program define _dea_parse_io, rclass
    version 16.0
    local raw0 `"`0'"'
    local scan `"`0'"'
    local lhs ""
    local form_a 0
    while `"`scan'"' != "" {
        gettoken word scan : scan, parse(" =,")
        if inlist(`"`word'"', "if", "in", ",") | substr(`"`word'"', 1, 3) == "if(" continue, break
        if `"`word'"' == "=" {
            local form_a 1
            continue, break
        }
        local lhs `lhs' `word'
    }
    if !`form_a' {
        return local args `"`raw0'"'
        exit
    }

    // Let Stata parse qualifiers (including functions containing commas) and
    // expand RHS varlists; never split the RHS at its first comma by hand.
    unab lhs : `lhs'
    local 0 `"`scan'"'
    syntax varlist(numeric min=1) [if] [in] [, *]
    return local args `"`if' `in', inputs(`lhs') outputs(`varlist') `options'"'
end
