*! version 0.2.0  20jul2026
*! Internal canonical LP contract for DEA solver adapters
program define _dea_lpcontract, rclass
    version 16.0
    syntax, C(name) [SENSE(string) Aeq(name) beq(name) Aub(name) bub(name) ///
        LB(name) UB(name) ENGINE(string) TOLerance(real 1e-8) ///
        SCALING(string) FORCESIMPLEXFAIL TRACE]

    if "`sense'" == "" local sense "min"
    local sense = lower("`sense'")
    if !inlist("`sense'", "min", "max") {
        di as error "_dea_lpcontract: sense() must be min or max"
        exit 198
    }

    if "`engine'" == "" local engine "auto"
    local engine = lower("`engine'")
    if !inlist("`engine'", "simplex", "ipm", "auto") {
        di as error "_dea_lpcontract: engine() must be simplex, ipm, or auto"
        exit 198
    }
    if `tolerance' <= 0 {
        di as error "_dea_lpcontract: tolerance() must be positive"
        exit 198
    }
    if "`scaling'" == "" local scaling "none"
    local scaling = lower("`scaling'")
    if !inlist("`scaling'", "none", "row") {
        di as error "_dea_lpcontract: scaling() must be none or row"
        exit 198
    }
    if "`forcesimplexfail'" != "" & "`engine'" != "auto" {
        di as error "_dea_lpcontract: forcesimplexfail is a test hook for engine(auto)"
        exit 198
    }

    foreach M in c aeq beq aub bub lb ub {
        if "``M''" != "" {
            capture confirm matrix ``M''
            if _rc {
                di as error "_dea_lpcontract: matrix ``M'' not found"
                exit 111
            }
        }
    }
    if ("`aeq'" == "") != ("`beq'" == "") {
        di as error "_dea_lpcontract: aeq() and beq() must be specified together"
        exit 198
    }
    if ("`aub'" == "") != ("`bub'" == "") {
        di as error "_dea_lpcontract: aub() and bub() must be specified together"
        exit 198
    }

    if inlist("`engine'", "simplex", "auto") {
        capture mata: mata which _rsm_solve_lp()
        if _rc {
            capture findfile _dea_simplex.ado
            if _rc {
                di as error "_dea_lpcontract: _dea_simplex.ado not found"
                exit 601
            }
            quietly run "`r(fn)'"
        }
    }

    mata: _dea_lpcontract_solve("`c'", "`sense'", "`aeq'", "`beq'", ///
        "`aub'", "`bub'", "`lb'", "`ub'", "`engine'", ///
        `tolerance', "`scaling'", "`forcesimplexfail'", "`trace'")

    local status = r(status)
    local objective = r(objective)
    local iterations = r(iterations)
    local fallback = r(fallback_used)
    local eqvio = r(max_eq_violation)
    local ineqvio = r(max_ineq_violation)
    local boundvio = r(max_bound_violation)
    local objscale = r(objective_scale)
    local rowscalemin = r(row_scale_min)
    local rowscalemax = r(row_scale_max)
    local status_text "`r(status_text)'"
    local engine_used "`r(engine)'"
    tempname solution
    matrix `solution' = r(x)

    return scalar status = `status'
    return scalar objective = `objective'
    return scalar iterations = `iterations'
    return scalar fallback_used = `fallback'
    return scalar max_eq_violation = `eqvio'
    return scalar max_ineq_violation = `ineqvio'
    return scalar max_bound_violation = `boundvio'
    return scalar objective_scale = `objscale'
    return scalar row_scale_min = `rowscalemin'
    return scalar row_scale_max = `rowscalemax'
    return matrix x = `solution'
    return local status_text "`status_text'"
    return local engine "`engine_used'"
    return local engine_requested "`engine'"
    return local scaling "`scaling'"
    return local sense "`sense'"
    return local cmd "_dea_lpcontract"
end

capture mata: mata drop _DEA_LPResult()
capture mata: mata drop _dea_lpcontract_ipm()
capture mata: mata drop _dea_lpcontract_simplex()
capture mata: mata drop _dea_lpcontract_diagnostics()
capture mata: mata drop _dea_lpcontract_dispatch()
capture mata: mata drop _dea_lpcontract_publish()
capture mata: mata drop _dea_lpcontract_solve()

mata:
mata set matastrict on

struct _DEA_LPResult {
    real scalar status
    real scalar objective
    real rowvector x
    real scalar iterations
    string scalar engine
    string scalar status_text
}

struct _DEA_LPResult scalar _dea_lpcontract_ipm(
    real rowvector c, string scalar sense,
    real matrix Aeq, real colvector beq,
    real matrix Aub, real colvector bub,
    real rowvector lb, real rowvector ub,
    real scalar tol, string scalar trace)
{
    struct _DEA_LPResult scalar out
    class LinearProgram scalar S
    real scalar rc, err
    real colvector pars

    out.status=3; out.objective=.; out.x=J(1,cols(c),.)
    out.iterations=.; out.engine="ipm"; out.status_text="failed"
    S=LinearProgram(); S.setMaxOrMin(sense); S.setCoefficients(c)
    if (rows(Aeq)>0) S.setEquality(Aeq,beq)
    if (rows(Aub)>0) S.setInequality(Aub,bub)
    S.setBounds(lb,ub); S.setTol(tol)
    S.setTrace(trace=="" ? "off" : "on")
    rc=S.optimize(); err=S.errorcode()
    if (err==1) {
        out.status=1; out.status_text="infeasible"
    }
    else if (err==2) {
        out.status=2; out.status_text="unbounded"
    }
    else if (S.converged() & err==0) {
        out.status=0; out.status_text="optimal"; out.objective=S.value()
        pars=S.parameters(); out.x=(rows(pars)>1 ? pars' : pars)
    }
    return(out)
}

struct _DEA_LPResult scalar _dea_lpcontract_simplex(
    real rowvector c, string scalar sense,
    real matrix Aeq, real colvector beq,
    real matrix Aub, real colvector bub,
    real rowvector lb, real rowvector ub,
    real scalar tol)
{
    struct _DEA_LPResult scalar out
    struct RSM_Result scalar raw
    real scalar n, j
    real rowvector cs
    real matrix As
    real colvector bs, beqs

    n=cols(c); out.status=3; out.objective=.; out.x=J(1,n,.)
    out.iterations=.; out.engine="simplex"; out.status_text="failed"

    /* The revised-simplex core uses y>=0.  Shift z=lb+y and express
       finite upper bounds as additional <= rows. */
    cs=(sense=="max" ? -c : c)
    As=Aub; bs=(rows(Aub)>0 ? bub-Aub*lb' : J(0,1,0))
    beqs=beq
    if (rows(Aeq)>0) beqs=beq-Aeq*lb'
    for (j=1; j<=n; j++) {
        if (ub[j]<.) {
            As=As\(I(n)[j,.])
            bs=bs\(ub[j]-lb[j])
        }
    }

    if (rows(As)==0 & rows(Aeq)==0) {
        if (any(cs:<-tol)) {
            out.status=2; out.status_text="unbounded"
        }
        else {
            out.status=0; out.status_text="optimal"; out.x=lb
            out.objective=c*out.x'; out.iterations=0
        }
        return(out)
    }

    raw=_rsm_solve_lp(As,bs,Aeq,beqs,cs,tol)
    out.iterations=raw.iters
    if (raw.status==0) {
        out.status=0; out.status_text="optimal"; out.x=lb+raw.x'
        out.objective=c*out.x'
    }
    else if (raw.status==1) {
        out.status=1; out.status_text="infeasible"
    }
    else if (raw.status==2) {
        out.status=2; out.status_text="unbounded"
    }
    return(out)
}

struct _DEA_LPResult scalar _dea_lpcontract_dispatch(
    real rowvector c, string scalar sense,
    real matrix Aeq, real colvector beq,
    real matrix Aub, real colvector bub,
    real rowvector lb, real rowvector ub,
    string scalar engine, real scalar tol, string scalar scaling,
    real scalar forcefail, string scalar trace,
    real scalar fallback, real scalar objscale,
    real scalar rowscalemin, real scalar rowscalemax)
{
    struct _DEA_LPResult scalar out
    real rowvector cs
    real matrix Aeqs, Aubs
    real colvector beqs, bubs
    real scalar i, den, mult

    cs=c; Aeqs=Aeq; beqs=beq; Aubs=Aub; bubs=bub
    fallback=0; objscale=1; rowscalemin=1; rowscalemax=1
    if (scaling=="row") {
        objscale=max(abs(c))
        if (objscale<=0) objscale=1
        cs=c/objscale
        rowscalemin=.; rowscalemax=0
        for (i=1; i<=rows(Aeqs); i++) {
            den=max((abs(Aeqs[i,.]),abs(beqs[i])))
            mult=(den>0 ? 1/den : 1)
            Aeqs[i,.]=Aeqs[i,.]*mult; beqs[i]=beqs[i]*mult
            rowscalemin=(rowscalemin==. ? mult : min((rowscalemin,mult)))
            rowscalemax=max((rowscalemax,mult))
        }
        for (i=1; i<=rows(Aubs); i++) {
            den=max((abs(Aubs[i,.]),abs(bubs[i])))
            mult=(den>0 ? 1/den : 1)
            Aubs[i,.]=Aubs[i,.]*mult; bubs[i]=bubs[i]*mult
            rowscalemin=(rowscalemin==. ? mult : min((rowscalemin,mult)))
            rowscalemax=max((rowscalemax,mult))
        }
        if (rowscalemin==.) rowscalemin=1
        if (rowscalemax==0) rowscalemax=1
    }

    if (engine=="ipm") out=_dea_lpcontract_ipm(cs,sense,Aeqs,beqs,Aubs,bubs,lb,ub,tol,trace)
    else {
        out=_dea_lpcontract_simplex(cs,sense,Aeqs,beqs,Aubs,bubs,lb,ub,tol)
        if (forcefail) {
            out.status=3; out.status_text="forced simplex failure"
            out.objective=.; out.x=J(1,cols(c),.)
        }
        if (engine=="auto" & out.status==3) {
            fallback=1
            out=_dea_lpcontract_ipm(cs,sense,Aeqs,beqs,Aubs,bubs,lb,ub,tol,trace)
        }
    }
    if (out.status==0) out.objective=c*out.x'
    return(out)
}

real rowvector _dea_lpcontract_diagnostics(
    real rowvector x, real matrix Aeq, real colvector beq,
    real matrix Aub, real colvector bub,
    real rowvector lb, real rowvector ub)
{
    real scalar eqv, ubv, lbv, bnduv, j
    eqv=(rows(Aeq)>0 ? max(abs(Aeq*x'-beq)) : 0)
    ubv=(rows(Aub)>0 ? max((Aub*x'-bub)\0) : 0)
    lbv=max((lb-x,J(1,cols(x),0)))
    bnduv=0
    for (j=1; j<=cols(ub); j++) {
        if (ub[j]<.) bnduv=max((bnduv,x[j]-ub[j]))
    }
    return((eqv,ubv,max((lbv,bnduv))))
}

void _dea_lpcontract_publish(struct _DEA_LPResult scalar out,
    real matrix Aeq, real colvector beq,
    real matrix Aub, real colvector bub,
    real rowvector lb, real rowvector ub,
    real scalar fallback, real scalar objscale,
    real scalar rowscalemin, real scalar rowscalemax)
{
    real rowvector diag
    diag=(out.status==0 ? _dea_lpcontract_diagnostics(out.x,Aeq,beq,Aub,bub,lb,ub) : J(1,3,.))
    st_numscalar("r(status)",out.status)
    st_numscalar("r(objective)",out.objective)
    st_numscalar("r(iterations)",out.iterations)
    st_numscalar("r(fallback_used)",fallback)
    st_numscalar("r(objective_scale)",objscale)
    st_numscalar("r(row_scale_min)",rowscalemin)
    st_numscalar("r(row_scale_max)",rowscalemax)
    st_numscalar("r(max_eq_violation)",diag[1])
    st_numscalar("r(max_ineq_violation)",diag[2])
    st_numscalar("r(max_bound_violation)",diag[3])
    st_matrix("r(x)",out.x)
    st_global("r(status_text)",out.status_text)
    st_global("r(engine)",out.engine)
}

void _dea_lpcontract_solve(string scalar cn, string scalar sense,
    string scalar aeqn, string scalar beqn,
    string scalar aubn, string scalar bubn,
    string scalar lbn, string scalar ubn,
    string scalar engine, real scalar tol, string scalar scaling,
    string scalar forcefail_s, string scalar trace)
{
    struct _DEA_LPResult scalar out
    real rowvector c, lb, ub
    real colvector beq, bub
    real matrix Aeq, Aub
    real scalar n, bad, fallback, objscale, rowscalemin, rowscalemax

    c=st_matrix(cn); if (rows(c)>1) c=c'
    n=cols(c)
    Aeq=(aeqn=="" ? J(0,n,0) : st_matrix(aeqn))
    beq=(beqn=="" ? J(0,1,0) : st_matrix(beqn))
    Aub=(aubn=="" ? J(0,n,0) : st_matrix(aubn))
    bub=(bubn=="" ? J(0,1,0) : st_matrix(bubn))
    lb=(lbn=="" ? J(1,n,0) : st_matrix(lbn)); if (rows(lb)>1) lb=lb'
    ub=(ubn=="" ? J(1,n,.) : st_matrix(ubn)); if (rows(ub)>1) ub=ub'
    if (rows(beq)==1 & cols(beq)>1) beq=beq'
    if (rows(bub)==1 & cols(bub)>1) bub=bub'

    bad=(n<1 | rows(c)!=1 | cols(lb)!=n | cols(ub)!=n |
        cols(Aeq)!=n | rows(Aeq)!=rows(beq) |
        cols(Aub)!=n | rows(Aub)!=rows(bub) |
        any(c:>=.) | any(Aeq:>=.) | any(beq:>=.) |
        any(Aub:>=.) | any(bub:>=.) | any(lb:>=.) | any(lb:>ub))
    if (bad) {
        out.status=3; out.objective=.; out.x=J(1,max((n,1)),.)
        out.iterations=.; out.engine="none"; out.status_text="invalid input"
        _dea_lpcontract_publish(out,Aeq,beq,Aub,bub,lb,ub,0,1,1,1)
        st_global("r(cmd)","_dea_lpcontract")
        st_global("r(sense)",sense)
        st_global("r(engine_requested)",engine)
        return
    }

    out=_dea_lpcontract_dispatch(c,sense,Aeq,beq,Aub,bub,lb,ub,
        engine,tol,scaling,forcefail_s!="",trace,
        fallback,objscale,rowscalemin,rowscalemax)
    _dea_lpcontract_publish(out,Aeq,beq,Aub,bub,lb,ub,fallback,
        objscale,rowscalemin,rowscalemax)
    st_global("r(cmd)","_dea_lpcontract")
    st_global("r(sense)",sense)
    st_global("r(engine_requested)",engine)
    st_global("r(scaling)",scaling)
}
end
