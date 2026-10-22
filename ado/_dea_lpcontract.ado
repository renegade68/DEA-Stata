*! version 0.2.0  20jul2026
*! Loader for the internal canonical LP contract
version 16.0

capture mata: mata which _rsm_solve_lp()
if _rc {
    capture findfile _dea_simplex.ado
    if _rc {
        di as error "_dea_lpcontract: _dea_simplex.ado not found"
        exit 601
    }
    quietly run "`r(fn)'"
}

capture findfile _dea_lpcontract_impl.do
if _rc {
    di as error "_dea_lpcontract: implementation file not found"
    exit 601
}
quietly run "`r(fn)'"
