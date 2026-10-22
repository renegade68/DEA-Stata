*! version 1.0.0  20jul2026
*! Dynamic slacks-based DEA with interperiod carry-over constraints

capture program drop dea_dynamic
program define dea_dynamic, rclass
    version 16.0

    syntax [if] [in], INPUTS(varlist numeric) OUTPUTS(varlist numeric) ///
        DMU(varname) TIME(varname numeric) CARRYOVER(varlist numeric) ///
        [CARRYTYPE(string) RTS(string) ORT(string) ORIENTATION(string) ///
         GENERATE(name) TOLERANCE(real 1e-7) DETAIL DEBUG]

    if "`ort'" != "" & "`orientation'" != "" {
        di as error "specify only one of ort() or orientation()"
        exit 198
    }
    if "`ort'" != "" local orientation "`ort'"
    if "`orientation'" == "" local orientation "nonoriented"
    if "`rts'" == "" local rts "vrs"
    if "`generate'" == "" local generate "_dyn"

    local orientation = lower("`orientation'")
    local rts = lower("`rts'")
    if inlist("`orientation'", "i", "in") local orientation "input"
    if inlist("`orientation'", "o", "out") local orientation "output"
    if inlist("`orientation'", "non", "non-oriented", "non_oriented") local orientation "nonoriented"
    if !inlist("`orientation'", "input", "output", "nonoriented") {
        di as error "orientation() must be input, output, or nonoriented"
        exit 198
    }
    if !inlist("`rts'", "crs", "vrs") {
        di as error "rts() must be crs or vrs"
        exit 198
    }
    if `tolerance' <= 0 {
        di as error "tolerance() must be positive"
        exit 198
    }

    local q : word count `carryover'
    if `"`carrytype'"' == "" {
        forvalues h = 1/`q' {
            local carrytype "`carrytype' desirable"
        }
    }
    local nqtype : word count `carrytype'
    if `nqtype' != `q' {
        di as error "carrytype() must contain one type for each carryover variable"
        exit 198
    }
    local normalized_types
    foreach typ of local carrytype {
        local typ = lower("`typ'")
        if "`typ'" == "good" local typ "desirable"
        if "`typ'" == "bad" local typ "undesirable"
        if "`typ'" == "free" local typ "discretionary"
        if "`typ'" == "fixed" local typ "nondiscretionary"
        if !inlist("`typ'", "desirable", "undesirable", "discretionary", "nondiscretionary") {
            di as error "carry types are desirable, undesirable, discretionary, or nondiscretionary"
            exit 198
        }
        local normalized_types "`normalized_types' `typ'"
    }
    local carrytype : list retok normalized_types

    confirm name `generate'_overall
    confirm name `generate'_period
    foreach v in `generate'_overall `generate'_period {
        capture confirm new variable `v'
        if _rc {
            di as error "variable `v' already exists"
            exit 110
        }
    }

    marksample touse, novarlist
    markout `touse' `inputs' `outputs' `carryover' `dmu' `time', strok
    quietly count if `touse'
    if r(N) == 0 error 2000
    local N = r(N)

    foreach v of varlist `inputs' `outputs' `carryover' {
        quietly count if `touse' & `v' <= 0
        if r(N) {
            di as error "inputs, outputs, and carry-overs must be strictly positive; check `v'"
            exit 459
        }
    }

    tempvar gid order tagd tagt dup
    quietly egen long `gid' = group(`dmu') if `touse'
    quietly egen byte `tagd' = tag(`gid') if `touse'
    quietly count if `tagd'
    local n_dmus = r(N)
    quietly egen byte `tagt' = tag(`time') if `touse'
    quietly count if `tagt'
    local T = r(N)
    if `T' < 2 {
        di as error "dynamic DEA requires at least two time periods"
        exit 459
    }
    if `N' != `n_dmus' * `T' {
        di as error "a strictly balanced DMU-by-time panel is required"
        exit 459
    }
    quietly bysort `gid' `time': generate byte `dup' = _N if `touse'
    quietly count if `touse' & `dup' != 1
    if r(N) {
        di as error "dmu() and time() must uniquely identify observations"
        exit 459
    }

    quietly generate long `order' = _n
    quietly generate double `generate'_overall = .
    quietly generate double `generate'_period = .
    sort `gid' `time'

    mata: _dea_dynamic_compute("`inputs'", "`outputs'", "`carryover'", ///
        "`carrytype'", "`touse'", "`generate'_overall", ///
        "`generate'_period", "`orientation'", "`rts'", `tolerance', ///
        "`debug'" != "")

    sort `order'
    quietly summarize `generate'_overall if `touse', meanonly
    local mean_eff = r(mean)
    quietly count if `touse' & missing(`generate'_overall)
    local n_failed = r(N) / `T'

    label variable `generate'_overall "Overall dynamic SBM efficiency"
    label variable `generate'_period "Period dynamic SBM efficiency"

    di as text ""
    di as text "Dynamic DEA (slacks-based measure)"
    di as text "{hline 70}"
    di as text "DMUs:        " as result `n_dmus'
    di as text "Periods:     " as result `T'
    di as text "Orientation: " as result "`orientation'"
    di as text "RTS:         " as result upper("`rts'")
    di as text "Carry types: " as result "`carrytype'"
    di as text "Mean overall efficiency: " as result %9.6f `mean_eff'
    di as text "Failed DMUs: " as result `n_failed'
    di as text "{hline 70}"
    if "`detail'" != "" {
        sort `gid' `time'
        list `dmu' `time' `generate'_period `generate'_overall if `touse', ///
            noobs sepby(`dmu') abbreviate(20)
        sort `order'
    }

    return scalar n_dmus = `n_dmus'
    return scalar T = `T'
    return scalar mean_eff = `mean_eff'
    return scalar n_failed = `n_failed'
    return local rts "`rts'"
    return local orientation "`orientation'"
    return local carrytype "`carrytype'"
    return local cmd "dea_dynamic"
end

capture mata: mata drop _dea_dynamic_compute()
capture mata: mata drop _dea_dynamic_one()

mata:
mata set matastrict on

real rowvector _dea_dynamic_one(real matrix X, real matrix Y, real matrix Z,
    string rowvector ctype, real scalar o, real scalar T, string scalar orientation,
    string scalar rts, real scalar tol, real scalar debug)
{
    real scalar n, m, s, q, nlam, nsi, nso, nsz, nk, nv
    real scalar t, j, i, h, row, col, rc, kappa, obj, at, bt, signz
    real scalar rho, newrho, iter
    real rowvector c, acoef, bcoef, lb, ub, sol, result
    real matrix Aeq
    real colvector beq
    class LinearProgram scalar LP

    n = rows(X)/T; m = cols(X); s = cols(Y); q = cols(Z)
    nlam=n*T; nsi=m*T; nso=s*T; nsz=q*T
    nk = 0
    nv=nlam+nsi+nso+nsz+nk
    c=J(1,nv,0); acoef=J(1,nv,0); bcoef=J(1,nv,0)
    lb=J(1,nv,0); ub=J(1,nv,.)

    /* Free-link deviations are unrestricted; fixed-link deviations are zero. */
    for (t=1; t<=T; t++) for (h=1; h<=q; h++) {
        col=nlam+nsi+nso+(t-1)*q+h
        if (ctype[h]=="discretionary") ub[col]=0
        if (ctype[h]=="nondiscretionary") ub[col]=0
    }

    /* m+s+q technology rows per term, T-1 continuity rows, and VRS rows. */
    Aeq=J(T*(m+s+q-sum(ctype:=="discretionary"))+(T-1)*q+(rts=="vrs")*T+nk,nv,0)
    beq=J(rows(Aeq),1,0); row=0
    for (t=1; t<=T; t++) {
        for (i=1; i<=m; i++) {
            row++; for (j=1; j<=n; j++) Aeq[row,(t-1)*n+j]=X[(j-1)*T+t,i]
            Aeq[row,nlam+(t-1)*m+i]=1
            if (nk) Aeq[row,nv]=-X[(o-1)*T+t,i]
            else beq[row]=X[(o-1)*T+t,i]
        }
        for (i=1; i<=s; i++) {
            row++; for (j=1; j<=n; j++) Aeq[row,(t-1)*n+j]=Y[(j-1)*T+t,i]
            Aeq[row,nlam+nsi+(t-1)*s+i]=-1
            if (nk) Aeq[row,nv]=-Y[(o-1)*T+t,i]
            else beq[row]=Y[(o-1)*T+t,i]
        }
        for (h=1; h<=q; h++) {
            if (ctype[h]=="discretionary") continue
            row++; for (j=1; j<=n; j++) Aeq[row,(t-1)*n+j]=Z[(j-1)*T+t,h]
            signz = (ctype[h]=="desirable" ? -1 : 1)
            if (ctype[h]=="nondiscretionary") signz=0
            Aeq[row,nlam+nsi+nso+(t-1)*q+h]=signz
            if (nk) Aeq[row,nv]=-Z[(o-1)*T+t,h]
            else beq[row]=Z[(o-1)*T+t,h]
        }
    }
    for (t=1; t<T; t++) for (h=1; h<=q; h++) {
        row++
        for (j=1; j<=n; j++) {
            Aeq[row,(t-1)*n+j]=Z[(j-1)*T+t,h]
            Aeq[row,t*n+j]=-Z[(j-1)*T+t,h]
        }
    }
    if (rts=="vrs") for (t=1; t<=T; t++) {
        row++; Aeq[row,(t-1)*n+1..t*n]=J(1,n,1)
        if (nk) Aeq[row,nv]=-1
        else beq[row]=1
    }

    /* Equal period and item weights, as in Tone and Tsutsui (2010). */
    for (t=1; t<=T; t++) {
        for (i=1; i<=m; i++) {
            col=nlam+(t-1)*m+i
            if (orientation=="input") c[col]=1/(T*(m+sum(ctype:=="undesirable"))*X[(o-1)*T+t,i])
            if (orientation=="nonoriented") acoef[col]=1/(T*(m+sum(ctype:=="undesirable"))*X[(o-1)*T+t,i])
        }
        for (i=1; i<=s; i++) {
            col=nlam+nsi+(t-1)*s+i
            if (orientation=="output") c[col]=1/(T*(s+sum(ctype:=="desirable"))*Y[(o-1)*T+t,i])
            if (orientation=="nonoriented") bcoef[col]=1/(T*(s+sum(ctype:=="desirable"))*Y[(o-1)*T+t,i])
        }
        for (h=1; h<=q; h++) {
            col=nlam+nsi+nso+(t-1)*q+h
            if (ctype[h]=="undesirable" & orientation=="input")
                c[col]=1/(T*(m+sum(ctype:=="undesirable"))*Z[(o-1)*T+t,h])
            if (ctype[h]=="undesirable" & orientation=="nonoriented")
                acoef[col]=1/(T*(m+sum(ctype:=="undesirable"))*Z[(o-1)*T+t,h])
            if (ctype[h]=="desirable" & orientation=="output")
                c[col]=1/(T*(s+sum(ctype:=="desirable"))*Z[(o-1)*T+t,h])
            if (ctype[h]=="desirable" & orientation=="nonoriented")
                bcoef[col]=1/(T*(s+sum(ctype:=="desirable"))*Z[(o-1)*T+t,h])
        }
    }

    result=J(1,T+2,.); rho=1
    if (orientation=="nonoriented") {
        for (iter=1; iter<=100; iter++) {
            c=acoef+rho*bcoef
            LP=LinearProgram(); LP.setMaxOrMin("max")
            LP.setCoefficients(c); LP.setEquality(Aeq,beq); LP.setBounds(lb,ub); LP.setTol(tol)
            (void) LP.optimize(); rc=LP.errorcode()
            if (rc!=0 | !LP.converged()) break
            sol=LP.parameters()
            newrho=(1-acoef*sol')/(1+bcoef*sol')
            if (abs(newrho-rho)<=tol) {
                rho=newrho
                break
            }
            rho=newrho
        }
        obj=rho
    }
    else {
        LP=LinearProgram(); LP.setMaxOrMin("max")
        LP.setCoefficients(c); LP.setEquality(Aeq,beq); LP.setBounds(lb,ub); LP.setTol(tol)
        (void) LP.optimize(); rc=LP.errorcode()
        if (rc==0 & LP.converged()) {
            sol=LP.parameters()
            obj=LP.value()
        }
    }
    if (rc!=0 | !LP.converged()) {
        if (debug) printf("DMU %g: LP return code %g\n",o,rc)
        return(result)
    }
    kappa=1
    if (kappa<=tol) return(result)

    if (orientation=="input") result[1]=1-obj
    else if (orientation=="output") result[1]=1/(1+obj)
    else result[1]=obj
    for (t=1; t<=T; t++) {
        at=0; bt=0
        for (i=1; i<=m; i++) at=at+sol[nlam+(t-1)*m+i]/kappa/(m+sum(ctype:=="undesirable"))/X[(o-1)*T+t,i]
        for (i=1; i<=s; i++) bt=bt+sol[nlam+nsi+(t-1)*s+i]/kappa/(s+sum(ctype:=="desirable"))/Y[(o-1)*T+t,i]
        for (h=1; h<=q; h++) {
            col=nlam+nsi+nso+(t-1)*q+h
            if (ctype[h]=="undesirable") at=at+sol[col]/kappa/(m+sum(ctype:=="undesirable"))/Z[(o-1)*T+t,h]
            if (ctype[h]=="desirable") bt=bt+sol[col]/kappa/(s+sum(ctype:=="desirable"))/Z[(o-1)*T+t,h]
        }
        if (orientation=="input") result[t+1]=1-at
        else if (orientation=="output") result[t+1]=1/(1+bt)
        else result[t+1]=(1-at)/(1+bt)
    }
    result[T+2]=rc
    return(result)
}

void _dea_dynamic_compute(string scalar inputs, string scalar outputs,
    string scalar carries, string scalar carrytypes, string scalar touse,
    string scalar overallvar, string scalar periodvar, string scalar orientation,
    string scalar rts, real scalar tol, real scalar debug)
{
    real matrix X,Y,Z,R
    real colvector use, ov, pv
    real scalar N,T,n,o,t
    string rowvector typ
    X=st_data(.,tokens(inputs),touse); Y=st_data(.,tokens(outputs),touse)
    Z=st_data(.,tokens(carries),touse); typ=tokens(carrytypes)
    use=selectindex(st_data(.,touse)); N=rows(X)
    /* Sorted input has one consecutive block per DMU; infer T from the
       number of repeated observations represented by the first carry panel. */
    T=strtoreal(st_local("T")); n=strtoreal(st_local("n_dmus"))
    ov=J(N,1,.); pv=J(N,1,.)
    for (o=1;o<=n;o++) {
        R=_dea_dynamic_one(X,Y,Z,typ,o,T,orientation,rts,tol,debug)
        for (t=1;t<=T;t++) {
            ov[(o-1)*T+t]=R[1]
            pv[(o-1)*T+t]=R[t+1]
        }
    }
    st_store(use,st_varindex(overallvar),ov)
    st_store(use,st_varindex(periodvar),pv)
}
end
