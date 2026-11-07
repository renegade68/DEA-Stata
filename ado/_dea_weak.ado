*! version 1.0.11  05nov2026
*! Convex joint weak-disposal DDF; active z and idle w, total intensity z+w.
version 16.0
capture mata: assert(_dea_nonte_lp_ready()==1)
if _rc {
    quietly findfile _dea_nonte_lp.ado
    quietly run "`r(fn)'"
}
capture mata: assert(_dea_native_retry_ready()==1)
if _rc {
    quietly findfile _dea_native_retry.ado
    quietly run "`r(fn)'"
}
capture mata: assert(_dea_simplex_audit_ready()==16 & findexternal("_rsm_duplicate_ratio()")!=NULL & findexternal("_rsm_duplicate_map()")!=NULL & findexternal("_rsm_duplicate_problem()")!=NULL & findexternal("_rsm_duplicate_accept()")!=NULL & findexternal("_rsm_duplicate_retry()")!=NULL & findexternal("_rsm_two_phase_core()")!=NULL & findexternal("_rsm_solve_lp()")!=NULL & findexternal("_rsm_quality_fast()")!=NULL & findexternal("_rsm_affine_quality()")!=NULL & findexternal("_rsm_affine_refine()")!=NULL & findexternal("_rsm_success_fastread()")!=NULL & findexternal("_rsm_success_reset()")!=NULL & findexternal("_rsm_success_read()")!=NULL & findexternal("_rsm_success_counts_read()")!=NULL & findexternal("_rsm_success_check()")!=NULL & findexternal("_rsm_success_apply()")!=NULL & findexternal("_rsm_xgrow()")!=NULL & findexternal("_rsm_xadd()")!=NULL & findexternal("_rsm_xmul()")!=NULL & findexternal("_rsm_xnormalize()")!=NULL & findexternal("_rsm_rowmap()")!=NULL & findexternal("_rsm_row_empty()")!=NULL & findexternal("_rsm_row_accept()")!=NULL & findexternal("_rsm_row_retry()")!=NULL & findexternal("_dea_link_reset()")!=NULL & findexternal("_dea_link_read()")!=NULL & findexternal("_dea_link_store()")!=NULL & findexternal("_dea_link_small()")!=NULL & findexternal("_dea_link_residual()")!=NULL & findexternal("_dea_link_basis()")!=NULL & findexternal("_dea_link_apply()")!=NULL & findexternal("_dea_link_stage2()")!=NULL)
if _rc {
    quietly findfile _dea_simplex.ado
    quietly run "`r(fn)'"
}
capture mata: assert(_dea_simplex_audit_ready()==16 & findexternal("_rsm_duplicate_ratio()")!=NULL & findexternal("_rsm_duplicate_map()")!=NULL & findexternal("_rsm_duplicate_problem()")!=NULL & findexternal("_rsm_duplicate_accept()")!=NULL & findexternal("_rsm_duplicate_retry()")!=NULL & findexternal("_rsm_two_phase_core()")!=NULL & findexternal("_rsm_solve_lp()")!=NULL & findexternal("_rsm_quality_fast()")!=NULL & findexternal("_rsm_affine_quality()")!=NULL & findexternal("_rsm_affine_refine()")!=NULL & findexternal("_rsm_success_fastread()")!=NULL & findexternal("_rsm_success_reset()")!=NULL & findexternal("_rsm_success_read()")!=NULL & findexternal("_rsm_success_counts_read()")!=NULL & findexternal("_rsm_success_check()")!=NULL & findexternal("_rsm_success_apply()")!=NULL & findexternal("_rsm_xgrow()")!=NULL & findexternal("_rsm_xadd()")!=NULL & findexternal("_rsm_xmul()")!=NULL & findexternal("_rsm_xnormalize()")!=NULL & findexternal("_rsm_rowmap()")!=NULL & findexternal("_rsm_row_empty()")!=NULL & findexternal("_rsm_row_accept()")!=NULL & findexternal("_rsm_row_retry()")!=NULL & findexternal("_dea_link_reset()")!=NULL & findexternal("_dea_link_read()")!=NULL & findexternal("_dea_link_store()")!=NULL & findexternal("_dea_link_small()")!=NULL & findexternal("_dea_link_residual()")!=NULL & findexternal("_dea_link_basis()")!=NULL & findexternal("_dea_link_apply()")!=NULL & findexternal("_dea_link_stage2()")!=NULL)
if _rc {
    di as error "Incompatible simplex library; reinstall the complete matched DEA package."
    exit 498
}

capture mata: assert(_dea_wb_ready()==3)
if _rc {
    quietly findfile _dea_weak_bound.ado
    quietly run "`r(fn)'"
}
capture mata: mata drop _dea_weak_ready()
capture mata: mata drop _dea_weak_compute()
mata:
mata set matastrict on
real scalar _dea_weak_ready() return(6)

void _dea_weak_compute(string scalar inputs,string scalar outputs,string scalar bad,
    string scalar touse,string scalar resultvars,string scalar direction,
    string scalar rts,string scalar engine,real scalar tol,real scalar matrices,real scalar debug,
    string scalar destinations)
{
    external real scalar _dea_rec_on
    real matrix X,Y,B,E,U,R,D,Z,W,TX,TY,TB,SI,SO,meta,EA,UA,F
    real colvector obs,e,u,gx,gy,gb,eb,ub
    real rowvector sx,sy,sb,c,lo,hi,sol,checks,active,idle
    real scalar n,m,s,q,k,v,rc,beta,eps,retry,sretry,backend,failed,eff
    string rowvector dest
    class LinearProgram scalar S
    struct _DEA_NativeRetry scalar recovered
    struct RSM_Result scalar res
    struct _DEA_WeakBound scalar bound
    X=st_data(.,tokens(inputs),touse);Y=st_data(.,tokens(outputs),touse)
    B=st_data(.,tokens(bad),touse);obs=selectindex(st_data(.,touse))
    n=rows(X);m=cols(X);s=cols(Y);q=cols(B);v=2*n+1
    sx=colmax(X);sx=sx:+(sx:==0);sy=colmax(Y);sy=sy:+(sy:==0)
    sb=colmax(B);sb=sb:+(sb:==0);X=X:/sx;Y=Y:/sy;B=B:/sb
    eps=min((tol,1e-9))
    R=J(n,2,.);R[,2]=J(n,1,1);D=J(n,10,.);D[,1]=J(n,1,1);D[,2..4]=J(n,3,0)
    F=J(n,7,.);F[,1]=J(n,1,0)
    Z=W=J(matrices ? n : 0,matrices ? n : 0,.)
    TX=SI=J(matrices ? n : 0,m,.)
    TY=SO=J(matrices ? n : 0,s,.);TB=J(matrices ? n : 0,q,.)
    c=(J(1,2*n,0),1);lo=J(1,v,0);hi=J(1,v,.)
    failed=eff=0
    for(k=1;k<=n;k++) {
        if (_dea_rec_on) _dea_rec_context(obs[k],4,rts)
        gx=(direction=="yb" ? J(m,1,0) : X[k,.]')
        gy=(direction=="x" ? J(s,1,0) : Y[k,.]')
        gb=(direction=="x" ? J(q,1,0) : B[k,.]')
        // X'(z+w)+beta gx<=x0; Y'z-beta gy>=y0; B'z+beta gb=b0.
        U=(X',X',gx)\(-Y',J(s,n,0),gy);u=X[k,.]'\(-Y[k,.]')
        E=(B',J(q,n,0),gb);e=B[k,.]'
        if(rts=="vrs") {
            E=E\(J(1,2*n,1),0);e=e\1
        }
        if(rts=="drs") {
            U=U\(J(1,2*n,1),0);u=u\1
        }
        if(rts=="irs") {
            U=U\(J(1,2*n,-1),0);u=u\(-1)
        }
        rc=1;retry=0;sretry=0;backend=1
        if(engine!="ipm") {
            res=_dea_nonte_simplex(U,u,E,e,-c,eps)
            rc=res.status
            if(rc==0) sol=res.x'
            if(rc!=0) {
                sretry=1
                // Alternate coordinates: total lambda, active mu, beta.
                UA=(X',J(m,n,0),gx)\(J(s,n,0),-Y',gy)\(-I(n),I(n),J(n,1,0))
                ub=X[k,.]'\(-Y[k,.]')\J(n,1,0)
                EA=(J(q,n,0),B',gb);eb=B[k,.]'
                if(rts=="vrs") {
                    EA=EA\(J(1,n,1),J(1,n+1,0));eb=eb\1
                }
                if(rts=="drs") {
                    UA=UA\(J(1,n,1),J(1,n+1,0));ub=ub\1
                }
                if(rts=="irs") {
                    UA=UA\(J(1,n,-1),J(1,n+1,0));ub=ub\(-1)
                }
                if (_dea_rec_on) _dea_rec_phase(1,2)
                res=_dea_nonte_simplex(UA,ub,EA,eb,-c,eps)
                rc=res.status
                if(rc==0) {
                    sol=res.x'
                    sol=(sol[n+1..2*n],sol[1..n]-sol[n+1..2*n],sol[v]);sretry=2
                }
            }
            if(rc==0) {
                checks=_dea_lpcheck(sol,E,e,U,u,lo,hi)
                if(!_dea_lpcheck_ok(checks,eps)) {
                    rc=1
                    if(sretry) sretry=1
                }
            }
        }
        // Research recovery: an exactly feasible self point, certified by a
        // conservative dual upper bound. Restricted to VRS/input direction.
        if(engine!="ipm" & rc!=0 & rts=="vrs" & direction=="x") {
            bound=_dea_wb_solve(X,Y,B,k,eps)
            F[k,.]=(1,bound.upper,bound.solver_status,bound.margin,bound.denominator,bound.corrections,bound.attempts)
            if(bound.status==0) {
                F[k,1]=3
                if(bound.upper<=eps) {
                    sol=J(1,v,0);sol[k]=1
                    checks=_dea_lpcheck(sol,E,e,U,u,lo,hi)
                    if(_dea_lpcheck_ok(checks,eps)) {
                        rc=0;F[k,1]=2
                    }
                }
            }
        }
        if(engine=="ipm" | (engine=="auto" & rc!=0)) {
            backend=2
            S=LinearProgram();S.setMaxOrMin("max");S.setCoefficients(c)
            S.setEquality(E,e);S.setInequality(U,u);S.setBounds(lo,hi);S.setTrace("off")
            rc=_dea_nonte_optimize(S,eps,1e-12)
            if(rc==0) sol=S.parameters()
            else {
                retry=1
                recovered=_dea_native_retry(c,"max",E,e,U,u,lo,hi,min((eps,1e-12)),"")
                rc=recovered.status
                if(rc==0) {
                    sol=recovered.x;retry=2
                }
            }
        }
        D[k,2..4]=(backend,sretry,retry)
        if(rc==0) {
            // Keep the accepted point, including tolerance-size roundoff.
            // No clipping, RHS relaxation or replacement of missing by zero.
            checks=_dea_lpcheck(sol,E,e,U,u,lo,hi)
            if(!_dea_lpcheck_ok(checks,eps)) rc=1
        }
        if(rc!=0) {
            if (_dea_rec_on) _dea_rec_group(obs[k],4,rts,0,.)
            if (_dea_rec_on) _dea_rec_final(obs[k],.,0)
            failed++;continue
        }
        beta=sol[v];
        if (_dea_rec_on) _dea_rec_group(obs[k],4,rts,1,.)
        if (_dea_rec_on) _dea_rec_final(obs[k],beta,(F[k,1]==2 ? 3 : backend))
        R[k,.]=(beta,0);D[k,.]=(0,backend,sretry,retry,checks)
        if(abs(beta)<=tol) eff++
        if(matrices) {
            active=sol[1..n];idle=sol[n+1..2*n];Z[k,.]=active;W[k,.]=idle
            TX[k,.]=(X[k,.]-beta*gx'):*sx
            TY[k,.]=(Y[k,.]+beta*gy'):*sy
            TB[k,.]=(B[k,.]-beta*gb'):*sb
            SI[k,.]=(X[k,.]-beta*gx'-(active+idle)*X):*sx
            SO[k,.]=(active*Y-Y[k,.]-beta*gy'):*sy
        }
    }
    st_store(obs,tokens(resultvars),R)
    dest=tokens(destinations)
    meta=(n,failed,eff,sum(D[,4]:>0),sum(D[,4]:==2),
        sum((D[,1]:==0):&(D[,2]:==1)),sum((D[,1]:==0):&(D[,2]:==2)),
        sum(D[,3]:>0),sum(D[,3]:==2),sum(D[,2]:==2),sum(F[,1]:>0),sum(F[,1]:==2))
    st_matrix(dest[1],meta)
    if(matrices) {
        st_matrix(dest[2],Z);st_matrix(dest[3],W);st_matrix(dest[4],Z+W)
        st_matrix(dest[5],TX);st_matrix(dest[6],TY);st_matrix(dest[7],TB)
        st_matrix(dest[8],SI);st_matrix(dest[9],SO);st_matrix(dest[10],D)
        st_matrix(dest[11],obs)
        st_matrix(dest[12],F)
    }
    if(debug) printf("Weak-disposal DDF: %g/%g scores; accepted simplex/native %g/%g\n",n-failed,n,meta[6],meta[7])
}
end
