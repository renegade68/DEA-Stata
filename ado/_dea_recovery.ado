*! version 0.1.0 03nov2026
*! Public simplex recovery context; research candidate, default off.
version 16.0
capture mata: mata drop _dea_recovery_ready()
capture mata: mata drop _dea_rec_abort()
capture mata: mata drop _dea_rec_begin()
capture mata: mata drop _dea_rec_close()
capture mata: mata drop _dea_rec_scope()
capture mata: mata drop _dea_rec_rts()
capture mata: mata drop _dea_rec_context()
capture mata: mata drop _dea_rec_phase()
capture mata: mata drop _dea_rec_emit()
capture mata: mata drop _dea_rec_status()
capture mata: mata drop _dea_rec_group()
capture mata: mata drop _dea_rec_final()
capture mata: mata drop _dea_rec_finalize_data()
capture mata: mata drop _dea_rec_export()
mata:
mata set matastrict on
real scalar _dea_recovery_ready() return(1)
void _dea_rec_abort()
{
    external real scalar _dea_rec_on,_dea_rec_limit,_dea_rec_scope_id
    external real matrix _dea_rec_rows
    external real rowvector _dea_rec_counts,_dea_rec_ctx
    _dea_rec_on=0;_dea_rec_limit=0;_dea_rec_scope_id=1
    _dea_rec_rows=J(0,25,.);_dea_rec_counts=J(1,8,0);_dea_rec_ctx=J(1,6,.)
}
void _dea_rec_begin(real scalar limit)
{
    external real scalar _dea_rec_on,_dea_rec_limit
    external real matrix _dea_rec_rows
    _dea_rec_abort()
    _dea_rec_limit=min((10000,max((1,floor(limit)))))
    _dea_rec_rows=J(_dea_rec_limit,25,.)
    _dea_rec_on=1
}
void _dea_rec_close()
{
    external real scalar _dea_rec_on
    _dea_rec_on=0
}
void _dea_rec_scope(real scalar scope)
{
    external real scalar _dea_rec_scope_id
    _dea_rec_scope_id=scope
}
real scalar _dea_rec_rts(string scalar rts)
{
    return(rts=="crs" ? 1 : (rts=="vrs" ? 2 : (rts=="drs" ? 3 : 4)))
}
void _dea_rec_context(real scalar obs,real scalar model,string scalar rts)
{
    external real scalar _dea_rec_on,_dea_rec_scope_id
    external real rowvector _dea_rec_ctx
    if(!_dea_rec_on) return
    _dea_rec_ctx=(obs,_dea_rec_scope_id,model,_dea_rec_rts(rts),1,1)
}
void _dea_rec_phase(real scalar stage,real scalar form)
{
    external real scalar _dea_rec_on
    external real rowvector _dea_rec_ctx
    if(!_dea_rec_on) return
    _dea_rec_ctx[5..6]=(stage,form)
}
void _dea_rec_emit(real rowvector row,real scalar m,real scalar n)
{
    external real scalar _dea_rec_on,_dea_rec_limit
    external real matrix _dea_rec_rows
    external real rowvector _dea_rec_counts,_dea_rec_ctx
    real scalar k
    if(!_dea_rec_on) return
    _dea_rec_counts[1]=_dea_rec_counts[1]+1
    _dea_rec_counts[2..6]=_dea_rec_counts[2..6]+row[(1,2,5,6,7)]
    k=_dea_rec_counts[1]
    if(k<=_dea_rec_limit) {
        _dea_rec_rows[k,.]=(_dea_rec_ctx,k,row,m,n,J(1,6,.))
        _dea_rec_counts[8]=k
    }
    else _dea_rec_counts[7]=_dea_rec_counts[7]+1
}
void _dea_rec_status(real scalar status)
{
    external real scalar _dea_rec_on,_dea_rec_limit
    external real matrix _dea_rec_rows
    external real rowvector _dea_rec_counts
    real scalar k
    if(!_dea_rec_on) return
    k=_dea_rec_counts[1]
    if(k>0 & k<=_dea_rec_limit) _dea_rec_rows[k,20]=status
}
void _dea_rec_group(real scalar obs,real scalar model,string scalar rts,
    real scalar score_available,real scalar reference_available)
{
    external real scalar _dea_rec_on,_dea_rec_scope_id
    external real matrix _dea_rec_rows
    external real rowvector _dea_rec_counts
    real colvector ix
    real scalar k
    if(!_dea_rec_on) return
    k=_dea_rec_counts[8];if(!k) return
    ix=selectindex((_dea_rec_rows[1..k,1]:==obs):&
        (_dea_rec_rows[1..k,2]:==_dea_rec_scope_id):&
        (_dea_rec_rows[1..k,3]:==model):&(_dea_rec_rows[1..k,4]:==_dea_rec_rts(rts)))
    if(rows(ix)) {
        _dea_rec_rows[ix,21]=J(rows(ix),1,score_available)
        _dea_rec_rows[ix,22]=J(rows(ix),1,reference_available)
    }
}
void _dea_rec_final(real scalar obs,real scalar score,real scalar source)
{
    external real scalar _dea_rec_on
    external real matrix _dea_rec_rows
    external real rowvector _dea_rec_counts
    real colvector ix
    real scalar k
    if(!_dea_rec_on) return
    k=_dea_rec_counts[8];if(!k) return
    ix=selectindex(_dea_rec_rows[1..k,1]:==obs)
    if(rows(ix)) {
        _dea_rec_rows[ix,23]=J(rows(ix),1,score<.)
        _dea_rec_rows[ix,24]=J(rows(ix),1,score)
        _dea_rec_rows[ix,25]=J(rows(ix),1,(score<. ? source : 0))
    }
}
void _dea_rec_finalize_data(string scalar scorevar,real scalar source)
{
    external real scalar _dea_rec_on
    external real matrix _dea_rec_rows
    external real rowvector _dea_rec_counts
    real scalar k,j,obs,score
    if(!_dea_rec_on) return
    k=_dea_rec_counts[8]
    for(j=1;j<=k;j++) {
        obs=_dea_rec_rows[j,1]
        if(obs>=1 & obs<=st_nobs() & obs==floor(obs)) {
            score=st_data(obs,scorevar)
            _dea_rec_rows[j,23..25]=(score<.,score,(score<. ? source : 0))
        }
    }
}
void _dea_rec_export(string scalar logname,string scalar statsname)
{
    external real matrix _dea_rec_rows
    external real rowvector _dea_rec_counts
    real scalar k
    string rowvector names
    st_matrix(statsname,_dea_rec_counts)
    st_matrixcolstripe(statsname,(J(8,1,""),tokens("lp_calls attempted recovered eta_iters factor_iters total_iters omitted stored")'))
    k=_dea_rec_counts[8]
    if(k) {
        st_matrix(logname,_dea_rec_rows[1..k,.])
        names=tokens("obsno scope model rts stage form lp_seq attempted recovered eta_status factor_status eta_iters factor_iters total_iters core_status eta_reason factor_reason n_rows n_cols caller_status model_score model_ref final_available final_score final_source")
        st_matrixcolstripe(logname,(J(25,1,""),names'))
    }
}
_dea_rec_abort()

end
