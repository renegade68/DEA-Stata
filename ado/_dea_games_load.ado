*! version 1.0.0 09nov2026
*! Versioned existing-helper load chain; call before replacing user variables.
program define _dea_games_load
 version 16.0
 quietly _dea_crossgame_load
 capture mata: assert(_t122_ready()==12201)
 if _rc {
  quietly findfile _dea_twostage_exact.ado
  quietly run "`r(fn)'"
  quietly findfile _dea_independent_exact.ado
  quietly run "`r(fn)'"
 }
 capture mata: assert(_dea_games_ready()==12901 & _t128_ready()==12801 & _g129_ready()==13001)
 if _rc {
  quietly findfile _dea_games_exact.ado
  quietly run "`r(fn)'"
  quietly findfile _dea_games_quality.ado
  quietly run "`r(fn)'"
  quietly findfile _dea_games.ado
  quietly run "`r(fn)'"
 }
 mata: assert(_dea_games_ready()==12901 & _t128_ready()==12801 & _g129_ready()==13001)
end
