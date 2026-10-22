*! version 1.14.0  20jul2026
*! Loader for optimal-set DEA reference-weight diagnostics
version 16.0

capture findfile _dea_tecontract.ado
if _rc {
    di as error "dea_reference: _dea_tecontract.ado not found"
    exit 601
}
quietly run "`r(fn)'"

capture findfile milp_bb_impl.do
if _rc {
    di as error "dea_reference: milp_bb_impl.do not found"
    exit 601
}
capture program drop milp_bb
quietly run "`r(fn)'"

capture findfile dea_reference_impl.do
if _rc {
    di as error "dea_reference: implementation file not found"
    exit 601
}
quietly run "`r(fn)'"
