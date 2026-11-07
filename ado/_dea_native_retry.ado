*! version 1.0.0  01nov2026
*! Bounded native-only retries; each candidate is checked in original coordinates.
version 16.0
capture mata: assert(_dea_lpcheck_ok(J(1,6,0),1e-8))
if _rc {
    quietly findfile _dea_lpcheck.ado
    quietly run "`r(fn)'"
}
capture mata: mata drop _DEA_NativeRetry()
capture mata: mata drop _dea_native_retry_ready()
capture mata: mata drop _dea_native_retry()
capture mata: mata drop _dea_native_record()
mata:
mata set matastrict on

struct _DEA_NativeRetry {
    real scalar status
    real rowvector x
    real matrix history
}

real scalar _dea_native_retry_ready()
{
    return(1)
}

// History columns: native code, checked status, convergence tolerance,
// iterations, maximum scaled primal residual, original objective.
// Native success rejected by the residual check is coded 90.
real rowvector _dea_native_record(real scalar native, real scalar checked,
    real scalar tol, real scalar iterations, real rowvector x,
    real rowvector c, real matrix AE, real colvector BE,
    real matrix AU, real colvector BU, real rowvector lb, real rowvector ub)
{
    real rowvector d
    real scalar residual, objective
    residual=.;objective=.
    if (cols(x)==cols(c) & !hasmissing(x)) {
        d=_dea_lpcheck(x,AE,BE,AU,BU,lb,ub)
        if (!hasmissing(d)) residual=max(d[4..6])
        objective=c*x'
    }
    return((native,checked,tol,iterations,residual,objective))
}

struct _DEA_NativeRetry scalar _dea_native_retry(
    real rowvector c, string scalar sense,
    real matrix AE, real colvector BE, real matrix AU, real colvector BU,
    real rowvector lb, real rowvector ub, real scalar tol, string scalar trace)
{
    struct _DEA_NativeRetry scalar out
    class LinearProgram scalar S
    real matrix EA,UA
    real colvector EB,UB
    real rowvector cc,x,d
    real scalar v,i,den,t,rc,status,changed

    out.status=3;out.x=J(1,cols(c),.);out.history=J(3,6,.)
    for (v=1;v<=3;v++) {
        EA=AE;EB=BE;UA=AU;UB=BU;cc=c;t=tol;changed=0
        if (v==1 & rows(AE)==0) continue
        if (v<=2 & rows(AE)>0) {
            UA=AU\AE\(-AE);UB=BU\BE\(-BE)
            EA=J(0,cols(c),0);EB=J(0,1,0)
        }
        if (v==2) {
            den=max(abs(c))
            if (den>0 & den!=1) { cc=c/den;changed=1; }
            for (i=1;i<=rows(UA);i++) {
                den=max((abs(UA[i,.]),abs(UB[i])))
                if (den>0 & den!=1) {
                    UA[i,.]=UA[i,.]/den;UB[i]=UB[i]/den;changed=1
                }
            }
            // Do not repeat an identical formulation with identical settings.
            if (!changed) continue
        }
        if (v==3) {
            t=min((tol,1e-12))
            if (t==tol) continue
        }
        S=LinearProgram();S.setMaxOrMin(sense);S.setCoefficients(cc)
        if (rows(EA)>0) S.setEquality(EA,EB)
        if (rows(UA)>0) S.setInequality(UA,UB)
        S.setBounds(lb,ub);S.setTol(t)
        S.setTrace(trace=="" ? "off" : "on")
        rc=S.optimize();x=J(1,cols(c),.)
        status=(S.errorcode()==0 ? 3 : S.errorcode())
        if (S.converged() & S.errorcode()==0) {
            x=S.parameters();d=_dea_lpcheck(x,AE,BE,AU,BU,lb,ub)
            status=(_dea_lpcheck_ok(d,tol) & c*x'<. ? 0 : 90)
        }
        out.history[v,.]=_dea_native_record(S.errorcode(),status,t,S.iterations(),
            x,c,AE,BE,AU,BU,lb,ub)
        if (status==0) { out.status=0;out.x=x;return(out); }
    }
    return(out)
}
end
