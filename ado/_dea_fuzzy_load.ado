*! version 1.0.0 09nov2026
program define _dea_fuzzy_load
 version 16.0
 quietly _dea_crossgame_load
 capture mata: assert(_t120_ready()==12001 & findexternal("_t120_eval()")!=NULL & findexternal("_t120_relerr()")!=NULL & findexternal("_t120_fit_reason()")!=NULL)
 if _rc {
  quietly findfile _dea_twostage_exact.ado
  quietly run "`r(fn)'"
 }
 capture mata: assert(_f135_ready()==13501)
 if _rc {
  quietly findfile _dea_fuzzy_exact.ado
  quietly run "`r(fn)'"
 }
 mata: assert(_f135_ready()==13501)
end
