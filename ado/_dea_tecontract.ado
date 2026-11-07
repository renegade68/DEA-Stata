*! version 0.7.0  01nov2026
*! Loader for the internal radial DEA canonical translator
version 16.0

capture mata: assert(findexternal("_dea_lpcontract_dispatch()")!=NULL)
if _rc {
    capture findfile _dea_lpcontract.ado
    if _rc {
        di as error "_dea_tecontract: _dea_lpcontract.ado not found"
        exit 601
    }
    quietly run "`r(fn)'"
}

capture findfile _dea_tecontract_impl.do
if _rc {
    di as error "_dea_tecontract: implementation file not found"
    exit 601
}
quietly run "`r(fn)'"
