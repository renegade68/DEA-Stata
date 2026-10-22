*! version 1.2.0  20jul2026
*! Loader for milp_bb and the shared continuous-LP adapter
version 16.0

capture mata: mata which _dea_lpcontract_ipm()
if _rc {
    capture findfile _dea_lpcontract.ado
    if _rc {
        di as error "milp_bb: _dea_lpcontract.ado not found"
        exit 601
    }
    quietly run "`r(fn)'"
}

capture findfile milp_bb_impl.do
if _rc {
    di as error "milp_bb: implementation file not found"
    exit 601
}
quietly run "`r(fn)'"
