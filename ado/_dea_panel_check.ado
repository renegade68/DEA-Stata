*! version 1.0.0  01nov2026
*! Validate selected DMU-period keys without sorting or modifying the dataset.
program define _dea_panel_check, rclass
    version 16.0
    syntax , DMU(varname) TIME(varname numeric) TOUSE(varname numeric)

    tempname duplicates
    mata: st_numscalar("`duplicates'", _dea_panel_duplicates("`dmu'", "`time'", "`touse'"))
    if `duplicates' > 0 {
        di as error "dmu() and time() must uniquely identify observations in the analysis sample"
        di as error "Repeated DMU-period keys: " `duplicates' " extra observation(s)"
        di as error "Resolve these repeated keys before estimation; observations are not automatically aggregated."
        exit 459
    }
    return scalar duplicates = 0
end

capture mata: mata drop _dea_panel_duplicates()
mata:
real scalar _dea_panel_duplicates(string scalar dmu, string scalar time,
                                  string scalar touse)
{
    real colvector ids, periods, index
    string colvector labels
    real matrix keys
    real scalar n, i, group

    periods = st_data(., time, touse)
    n = rows(periods)
    if (n < 2) return(0)
    if (substr(st_vartype(dmu), 1, 3) == "str") {
        // Group exact strings, including strL, without truncation or trimming.
        labels = st_sdata(., dmu, touse)
        index = order(labels, 1)
        ids = J(n, 1, .)
        group = 1
        ids[index[1]] = group
        for (i = 2; i <= n; i++) {
            if (labels[index[i]] != labels[index[i-1]]) group++
            ids[index[i]] = group
        }
    }
    else ids = st_data(., dmu, touse)

    // Retain double precision for numeric IDs and periods; formatted labels
    // must never determine identity. Count excess rows, not all group members.
    keys = sort((ids, periods), (1, 2))
    return(sum((keys[2..n, 1] :== keys[1..n-1, 1]) :&
               (keys[2..n, 2] :== keys[1..n-1, 2])))
}
end
