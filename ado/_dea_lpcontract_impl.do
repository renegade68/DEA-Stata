*! version 0.4.9  05nov2026
*! Internal canonical LP contract for DEA solver adapters
capture program drop _dea_lpcontract
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
    if missing(`tolerance') | `tolerance' <= 0 {
        di as error "_dea_lpcontract: tolerance() must be positive and finite"
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
        capture mata: assert(_dea_simplex_audit_ready()==16 & findexternal("_rsm_duplicate_ratio()")!=NULL & findexternal("_rsm_duplicate_map()")!=NULL & findexternal("_rsm_duplicate_problem()")!=NULL & findexternal("_rsm_duplicate_accept()")!=NULL & findexternal("_rsm_duplicate_retry()")!=NULL & findexternal("_rsm_two_phase_core()")!=NULL & findexternal("_rsm_solve_lp()")!=NULL & findexternal("_rsm_quality_fast()")!=NULL & findexternal("_rsm_affine_quality()")!=NULL & findexternal("_rsm_affine_refine()")!=NULL & findexternal("_rsm_success_fastread()")!=NULL & findexternal("_rsm_success_reset()")!=NULL & findexternal("_rsm_success_read()")!=NULL & findexternal("_rsm_success_counts_read()")!=NULL & findexternal("_rsm_success_check()")!=NULL & findexternal("_rsm_success_apply()")!=NULL & findexternal("_rsm_xgrow()")!=NULL & findexternal("_rsm_xadd()")!=NULL & findexternal("_rsm_xmul()")!=NULL & findexternal("_rsm_xnormalize()")!=NULL & findexternal("_rsm_rowmap()")!=NULL & findexternal("_rsm_row_empty()")!=NULL & findexternal("_rsm_row_accept()")!=NULL & findexternal("_rsm_row_retry()")!=NULL & findexternal("_dea_link_reset()")!=NULL & findexternal("_dea_link_read()")!=NULL & findexternal("_dea_link_store()")!=NULL & findexternal("_dea_link_small()")!=NULL & findexternal("_dea_link_residual()")!=NULL & findexternal("_dea_link_basis()")!=NULL & findexternal("_dea_link_apply()")!=NULL & findexternal("_dea_link_stage2()")!=NULL)
        if _rc {
            capture findfile _dea_simplex.ado
            if _rc {
                di as error "_dea_lpcontract: _dea_simplex.ado not found"
                exit 601
            }
            quietly run "`r(fn)'"
        }
    }

    if inlist("`engine'", "simplex", "auto") {
        capture mata: assert(_dea_simplex_audit_ready()==16 & findexternal("_rsm_duplicate_ratio()")!=NULL & findexternal("_rsm_duplicate_map()")!=NULL & findexternal("_rsm_duplicate_problem()")!=NULL & findexternal("_rsm_duplicate_accept()")!=NULL & findexternal("_rsm_duplicate_retry()")!=NULL & findexternal("_rsm_two_phase_core()")!=NULL & findexternal("_rsm_solve_lp()")!=NULL & findexternal("_rsm_quality_fast()")!=NULL & findexternal("_rsm_affine_quality()")!=NULL & findexternal("_rsm_affine_refine()")!=NULL & findexternal("_rsm_success_fastread()")!=NULL & findexternal("_rsm_success_reset()")!=NULL & findexternal("_rsm_success_read()")!=NULL & findexternal("_rsm_success_counts_read()")!=NULL & findexternal("_rsm_success_check()")!=NULL & findexternal("_rsm_success_apply()")!=NULL & findexternal("_rsm_xgrow()")!=NULL & findexternal("_rsm_xadd()")!=NULL & findexternal("_rsm_xmul()")!=NULL & findexternal("_rsm_xnormalize()")!=NULL & findexternal("_rsm_rowmap()")!=NULL & findexternal("_rsm_row_empty()")!=NULL & findexternal("_rsm_row_accept()")!=NULL & findexternal("_rsm_row_retry()")!=NULL & findexternal("_dea_link_reset()")!=NULL & findexternal("_dea_link_read()")!=NULL & findexternal("_dea_link_store()")!=NULL & findexternal("_dea_link_small()")!=NULL & findexternal("_dea_link_residual()")!=NULL & findexternal("_dea_link_basis()")!=NULL & findexternal("_dea_link_apply()")!=NULL & findexternal("_dea_link_stage2()")!=NULL)
        if _rc {
            di as error "Incompatible simplex library; reinstall the complete matched DEA package."
            exit 498
        }

    }

    // Direct internal callers must not enter a partially cleared adapter.
    capture mata: assert(_dea_lpcontract_ready()==9 & findexternal("_rsm_public_success()")!=NULL & findexternal("_dea_lpcontract_simplex()")!=NULL & findexternal("_dea_lpcontract_dispatch()")!=NULL)
    if _rc {
        di as error "Incomplete LP adapter; save data and clear all, then reload the matched DEA package."
        exit 498
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
    local seqvio = r(max_scaled_eq_violation)
    local siqvio = r(max_scaled_ineq_violation)
    local sbvio = r(max_scaled_bound_violation)
    local reported_status = r(reported_status)
    local native_retries = r(native_retries)
    local conflict = r(status_conflict)
    local initial_status = r(initial_status)
    local certificate "`r(certificate)'"
    tempname attempts
    matrix `attempts' = r(attempts)
    matrix colnames `attempts' = engine variant native checked tolerance iterations residual objective
    local objscale = r(objective_scale)
    local rowscalemin = r(row_scale_min)
    local rowscalemax = r(row_scale_max)
    local status_text "`r(status_text)'"
    local engine_used "`r(engine)'"
    tempname solution
    matrix `solution' = r(x)
    tempname candidate
    matrix `candidate' = r(candidate_x)

    return scalar status = `status'
    return scalar objective = `objective'
    return scalar iterations = `iterations'
    return scalar fallback_used = `fallback'
    return scalar max_eq_violation = `eqvio'
    return scalar max_ineq_violation = `ineqvio'
    return scalar max_bound_violation = `boundvio'
    return scalar max_scaled_eq_violation = `seqvio'
    return scalar max_scaled_ineq_violation = `siqvio'
    return scalar max_scaled_bound_violation = `sbvio'
    return scalar reported_status = `reported_status'
    return scalar initial_status = `initial_status'
    return scalar native_retries = `native_retries'
    return scalar status_conflict = `conflict'
    return matrix attempts = `attempts'
    return local certificate "`certificate'"
    return scalar feasibility_tolerance = `tolerance'
    return matrix candidate_x = `candidate'
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

capture mata: mata drop _dea_lpcontract_ready()
capture mata: mata drop _DEA_LPResult()
capture mata: mata drop _dea_lpcontract_ipm()
capture mata: mata drop _rsm_public_success()
capture mata: mata drop _dea_lpcontract_simplex()
capture mata: mata drop _dea_lpcontract_diagnostics()
capture mata: mata drop _dea_lpcontract_dispatch()
capture mata: mata drop _dea_lpcontract_pow2()
capture mata: mata drop _dea_lpcontract_verify()
capture mata: mata drop _dea_lpcontract_publish()
capture mata: mata drop _dea_lpcontract_solve()
capture mata: mata drop _dea_lpcontract_retry()
capture mata: mata drop _dea_lpcontract_certificate()

mata:
mata set matastrict on
real scalar _dea_lpcontract_ready() return(9)

struct _DEA_LPResult {
    real scalar status
    real scalar objective
    real rowvector x
    real rowvector candidate
    real rowvector residual
    real scalar reported_status
    real scalar status_conflict
    real matrix attempts
    string scalar certificate
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
    out.candidate=out.x; out.residual=J(1,6,.); out.reported_status=3
    out.status_conflict=0;out.certificate="none";out.attempts=J(0,8,.)
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
    out.reported_status=out.status;out.iterations=S.iterations()
    out.attempts=(2,0,S.errorcode(),out.status,tol,out.iterations,.,out.objective)
    return(out)
}

// Check actual user coordinates before publishing. Basis column order is
// unchanged by positive row/objective scaling and lower-bound translation.
void _rsm_public_success(struct RSM_Result scalar raw,
 real rowvector c,string scalar sense,real matrix E,real colvector b,
 real matrix U,real colvector u,real rowvector lb,real rowvector ub,
 real scalar tol,real rowvector actual)
{
 external real scalar _rsm_success_policy
 struct RSM_Prob scalar P
 struct RSM_Result scalar full
 real matrix A
 real colvector B,L,slack
 real rowvector cost
 real scalar n,j,ni
 actual=lb+raw.x';n=cols(c)
 if(_rsm_success_policy==0) return
 cost=(sense=="max" ? -c : c);A=U;B=u
 for(j=1;j<=n;j++) if(ub[j]<.) {A=A\I(n)[j,.];B=B\ub[j];}
 ni=rows(A)
 _rsm_std_form(P,A,B,E,b,cost)
 slack=(ni>0 ? rowmax((J(ni,1,0),B-A*actual')) : J(0,1,0))
 L=lb'\J(ni,1,0);full=raw
 full.x=actual'\slack;full.objval=cost*actual'
 _rsm_success_apply(P,full,raw.basiscols,tol,L)
 raw.status=full.status;raw.msg=full.msg;raw.recovery=full.recovery
 if(full.status==0) actual=full.x[1..n]'
 else {actual=J(1,n,.);raw.x=J(n,1,.);raw.objval=.;}
}

struct _DEA_LPResult scalar _dea_lpcontract_simplex(
    real rowvector c, string scalar sense,
    real matrix Aeq, real colvector beq,
    real matrix Aub, real colvector bub,
    real rowvector lb, real rowvector ub,
    real scalar tol,| real rowvector rawc,real matrix rawE,real colvector rawb,real matrix rawU,real colvector rawu)
{
    struct _DEA_LPResult scalar out
    struct RSM_Result scalar raw
    real scalar n, j
    real rowvector cs,actual
    real matrix As
    real colvector bs, beqs

    n=cols(c); out.status=3; out.objective=.; out.x=J(1,n,.)
    out.candidate=out.x; out.residual=J(1,6,.); out.reported_status=3
    out.status_conflict=0;out.certificate="none";out.attempts=J(0,8,.)
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
        if (any(cs:<0)) {
            out.status=2; out.status_text="unbounded"
        }
        else {
            out.status=0; out.status_text="optimal"; out.x=lb
            out.objective=c*out.x'; out.iterations=0
        }
        out.reported_status=out.status
        out.attempts=(1,0,out.reported_status,out.status,tol,out.iterations,.,out.objective)
        return(out)
    }

    raw=_rsm_solve_lp(As,bs,Aeq,beqs,cs,tol,0)
    out.iterations=raw.iters
    out.reported_status=raw.status
    if(raw.status==0) {
        if(args()==9) _rsm_public_success(raw,c,sense,Aeq,beq,Aub,bub,lb,ub,tol,actual)
        else _rsm_public_success(raw,rawc,sense,rawE,rawb,rawU,rawu,lb,ub,tol,actual)
    }
    if (raw.status==3 & raw.msg=="Residual check failed" & rows(raw.x)>=n) {
        out.candidate=lb+raw.x[1..n]'
        out.status_text="residual check failed"
    }
    if (raw.status==0) {
        out.status=0; out.status_text="optimal"; out.x=actual
        out.objective=c*out.x'
    }
    else if (raw.status==1) {
        out.status=1; out.status_text="infeasible"
    }
    else if (raw.status==2) {
        out.status=2; out.status_text="unbounded"
    }
    out.attempts=(1,0,out.reported_status,out.status,tol,out.iterations,.,out.objective)
    return(out)
}

void _dea_lpcontract_verify(struct _DEA_LPResult scalar out,
    real rowvector c, real matrix AE, real colvector BE,
    real matrix AU, real colvector BU, real rowvector lb, real rowvector ub,
    real scalar tol)
{
    real scalar a
    if (out.status==0) out.candidate=out.x
    if (cols(out.candidate)==cols(c) & !hasmissing(out.candidate)) {
        out.residual=_dea_lpcheck(out.candidate,AE,BE,AU,BU,lb,ub)
    }
    if (out.status==0) {
        out.objective=c*out.x'
        if (!_dea_lpcheck_ok(out.residual,tol) | out.objective>=.) {
            out.status=3; out.status_text="residual check failed"
            out.x=J(1,cols(c),.); out.objective=.
        }
    }
    a=rows(out.attempts)
    if (a>0) {
        out.attempts[a,4]=(out.status_text=="residual check failed" ? 90 : out.status)
        if (!hasmissing(out.residual)) out.attempts[a,7]=max(out.residual[4..6])
        if (!hasmissing(out.candidate)) out.attempts[a,8]=c*out.candidate'
    }
}

void _dea_lpcontract_retry(struct _DEA_LPResult scalar out,
    real rowvector c, string scalar sense, real matrix AE, real colvector BE,
    real matrix AU, real colvector BU, real rowvector lb, real rowvector ub,
    real scalar tol, string scalar trace)
{
    struct _DEA_NativeRetry scalar retry
    real scalar i
    if (out.status==0 | out.certificate!="none") return
    retry=_dea_native_retry(c,sense,AE,BE,AU,BU,lb,ub,tol,trace)
    for (i=1;i<=3;i++) {
        if (retry.history[i,1]<.) out.attempts=out.attempts\(2,i,retry.history[i,.])
    }
    if (retry.status==0 & !any(out.attempts[.,3]:==2)) {
        out.status=0;out.x=retry.x;out.objective=c*out.x'
        out.candidate=out.x;out.reported_status=0;out.status_text="optimal"
        out.iterations=out.attempts[rows(out.attempts),6]
        _dea_lpcontract_verify(out,c,AE,BE,AU,BU,lb,ub,tol)
    }
    else if (any(out.attempts[.,3]:==2) &
        (any(out.attempts[.,3]:==1) | retry.status==0)) {
        out.status=3;out.status_text="conflicting native status";out.status_conflict=1
        out.x=J(1,cols(c),.);out.objective=.
    }
    else if (any((out.attempts[.,4]:!=1):&(out.attempts[.,4]:!=2))) {
        // Numerical failures or rejected candidates do not confirm an earlier
        // infeasible/unbounded report. Keep the uncertainty visible.
        out.status=3
        if (out.status_text!="residual check failed") out.status_text="native retries exhausted"
        out.x=J(1,cols(c),.);out.objective=.
    }
}

// These are structural proofs, using signs/zeros only: no residual tolerance
// can turn a small nonzero coefficient into a recession direction.
void _dea_lpcontract_certificate(struct _DEA_LPResult scalar out,
    real rowvector c, string scalar sense, real matrix AE, real colvector BE,
    real matrix AU, real colvector BU, real rowvector lb, real rowvector ub)
{
    real scalar i,j,nonnegative
    nonnegative=all(lb:>=0)
    for (i=1;i<=rows(AE);i++) {
        if (all(AE[i,.]:==0) & BE[i]!=0) out.certificate="zero_row"
        else if (nonnegative & ((all(AE[i,.]:>=0) & BE[i]<0) |
            (all(AE[i,.]:<=0) & BE[i]>0))) out.certificate="nonnegative_row"
    }
    for (i=1;i<=rows(AU);i++) {
        if (all(AU[i,.]:==0) & BU[i]<0) out.certificate="zero_row"
        else if (nonnegative & all(AU[i,.]:>=0) & BU[i]<0) out.certificate="nonnegative_row"
    }
    if (out.certificate!="none") {
        out.status=1;out.status_text="infeasible";out.x=J(1,cols(c),.);out.objective=.
        return
    }
    // The origin is an exact feasible point; a coordinate ray is exact too.
    if (!all(lb:<=0) | !all(ub:>=0) | (rows(BE)>0 & !all(BE:==0)) | (rows(BU)>0 & !all(BU:>=0))) return
    for (j=1;j<=cols(c);j++) {
        if (ub[j]==. & (rows(AE)==0 | all(AE[.,j]:==0)) & (rows(AU)==0 | all(AU[.,j]:<=0)) &
            ((sense=="max" & c[j]>0) | (sense=="min" & c[j]<0))) {
            out.status=2;out.status_text="unbounded";out.certificate="zero_origin_ray"
            out.x=J(1,cols(c),.);out.objective=.
            return
        }
    }
}


// Power-of-two initial scaling with exact coefficient roundtrip.
// If the requested binary normalization loses a stored coefficient, use 1.
real scalar _dea_lpcontract_pow2(real rowvector values)
{
    real scalar den,mult,next,k
    real rowvector scaled
    if(hasmissing(values) | cols(values)==0) return(1)
    den=max(abs(values));if(den==0) return(1)
    mult=1
    for(k=0;k<2048;k++) {
        if(den>=2) {next=mult/2;if(next==0) return(1);mult=next;den=den/2;}
        else if(den<1) {next=mult*2;if(next>=.) return(1);mult=next;den=den*2;}
        else break
    }
    if(k==2048 | mult<=0 | mult>=. | 1/mult>=.) return(1)
    scaled=values*mult
    if(hasmissing(scaled) | any(scaled/mult:!=values)) return(1)
    return(mult)
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
    struct _DEA_LPResult scalar out, first
    real rowvector cs
    real matrix Aeqs, Aubs
    real colvector beqs, bubs
    real scalar i, den, mult

    cs=c; Aeqs=Aeq; beqs=beq; Aubs=Aub; bubs=bub
    fallback=0; objscale=1; rowscalemin=1; rowscalemax=1
    if (scaling=="row") {
        mult=_dea_lpcontract_pow2(c)
        objscale=1/mult
        cs=c*mult
        rowscalemin=.; rowscalemax=0
        for (i=1; i<=rows(Aeqs); i++) {
            mult=_dea_lpcontract_pow2((Aeqs[i,.],beqs[i]))
            Aeqs[i,.]=Aeqs[i,.]*mult; beqs[i]=beqs[i]*mult
            rowscalemin=(rowscalemin==. ? mult : min((rowscalemin,mult)))
            rowscalemax=max((rowscalemax,mult))
        }
        for (i=1; i<=rows(Aubs); i++) {
            mult=_dea_lpcontract_pow2((Aubs[i,.],bubs[i]))
            Aubs[i,.]=Aubs[i,.]*mult; bubs[i]=bubs[i]*mult
            rowscalemin=(rowscalemin==. ? mult : min((rowscalemin,mult)))
            rowscalemax=max((rowscalemax,mult))
        }
        if (rowscalemin==.) rowscalemin=1
        if (rowscalemax==0) rowscalemax=1
    }

    if (engine=="ipm") {
        out=_dea_lpcontract_ipm(cs,sense,Aeqs,beqs,Aubs,bubs,lb,ub,tol,trace)
        _dea_lpcontract_verify(out,c,Aeq,beq,Aub,bub,lb,ub,tol)
        _dea_lpcontract_certificate(out,c,sense,Aeq,beq,Aub,bub,lb,ub)
        _dea_lpcontract_retry(out,c,sense,Aeq,beq,Aub,bub,lb,ub,tol,trace)
    }
    else {
        out=_dea_lpcontract_simplex(cs,sense,Aeqs,beqs,Aubs,bubs,lb,ub,tol,c,Aeq,beq,Aub,bub)
        _dea_lpcontract_verify(out,c,Aeq,beq,Aub,bub,lb,ub,tol)
        _dea_lpcontract_certificate(out,c,sense,Aeq,beq,Aub,bub,lb,ub)
        if (forcefail) {
            out.status=3;out.status_text="forced simplex failure";out.certificate="none"
            out.objective=.;out.x=J(1,cols(c),.)
            out.attempts[rows(out.attempts),4]=3
        }
        if (engine=="auto" & out.status!=0 & out.certificate=="none") {
            first=out;fallback=1
            out=_dea_lpcontract_ipm(cs,sense,Aeqs,beqs,Aubs,bubs,lb,ub,tol,trace)
            _dea_lpcontract_verify(out,c,Aeq,beq,Aub,bub,lb,ub,tol)
            _dea_lpcontract_certificate(out,c,sense,Aeq,beq,Aub,bub,lb,ub)
            _dea_lpcontract_retry(out,c,sense,Aeq,beq,Aub,bub,lb,ub,tol,trace)
            if (out.certificate=="none" & (first.status==1 | first.status==2)) {
                if ((first.status==1 & out.status==2) |
                    (first.status==2 & (out.status==0 | out.status==1))) {
                    out.status=3;out.status_text="conflicting solver status";out.status_conflict=1
                    out.x=J(1,cols(c),.);out.objective=.
                }
                else if (out.status==3) out.status_text="unverified solver status"
            }
            out.attempts=first.attempts\out.attempts
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
    diag=out.residual
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
    st_numscalar("r(max_scaled_eq_violation)",diag[4])
    st_numscalar("r(max_scaled_ineq_violation)",diag[5])
    st_numscalar("r(max_scaled_bound_violation)",diag[6])
    st_numscalar("r(reported_status)",out.reported_status)
    st_numscalar("r(initial_status)",(rows(out.attempts)>0 ? out.attempts[1,3] : .))
    st_numscalar("r(native_retries)",sum((out.attempts[.,1]:==2):&(out.attempts[.,2]:>0)))
    st_numscalar("r(status_conflict)",out.status_conflict)
    st_matrix("r(attempts)",(rows(out.attempts)>0 ? out.attempts : J(1,8,.)))
    st_global("r(certificate)",out.certificate)
    st_matrix("r(candidate_x)",out.candidate)
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

    bad=(n<1 | rows(c)!=1 | rows(lb)!=1 | rows(ub)!=1 |
        cols(lb)!=n | cols(ub)!=n | cols(Aeq)!=n | cols(Aub)!=n |
        cols(beq)!=1 | cols(bub)!=1 | rows(Aeq)!=rows(beq) |
        rows(Aub)!=rows(bub))
    if (!bad) bad=(any(c:>=.) | any(Aeq:>=.) | any(beq:>=.) |
        any(Aub:>=.) | any(bub:>=.) | any(lb:>=.) | any(lb:>ub) | any(ub:>.))
    if (bad) {
        out.status=3; out.objective=.; out.x=J(1,max((n,1)),.)
        out.status_conflict=0;out.certificate="none";out.attempts=J(0,8,.)
        out.iterations=.; out.engine="none"; out.status_text="invalid input"
        out.candidate=out.x; out.residual=J(1,6,.); out.reported_status=3
        // Do not pass malformed vectors through the typed diagnostic API.
        // Status is nonoptimal, so the published residuals remain missing.
        _dea_lpcontract_publish(out,J(0,n,0),J(0,1,0),J(0,n,0),
            J(0,1,0),J(1,n,0),J(1,n,.),0,1,1,1)
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
