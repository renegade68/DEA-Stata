*! version 0.5.0  20jul2026
*! Internal stage-1 radial technical-efficiency translator
program define _dea_tecontract, rclass
    version 16.0
    syntax [if] [in], INPUTS(varlist numeric min=1) ///
        OUTPUTS(varlist numeric min=1) GENerate(name) ///
        [STATUS(name) RTS(string) ORIentation(string) STAGE(integer 1) ///
         ENGINE(string) SCALING(string) SLACKNormalization(string) ///
         RADIALTolerance(real 0) TOLerance(real 1e-8) TRACE]

    if "`rts'" == "" local rts "vrs"
    local rts = lower("`rts'")
    if !inlist("`rts'", "crs", "vrs", "drs", "irs") {
        di as error "_dea_tecontract: rts() must be crs, vrs, drs, or irs"
        exit 198
    }
    if "`orientation'" == "" local orientation "input"
    local orientation = lower("`orientation'")
    if inlist("`orientation'", "i", "in", "inp", "inpu") local orientation "input"
    if inlist("`orientation'", "o", "out", "outp", "outpu") local orientation "output"
    if !inlist("`orientation'", "input", "output") {
        di as error "_dea_tecontract: orientation() must be input or output"
        exit 198
    }
    if !inlist(`stage', 1, 2) {
        di as error "_dea_tecontract: stage() must be 1 or 2"
        exit 198
    }
    if "`engine'" == "" local engine "auto"
    local engine = lower("`engine'")
    if !inlist("`engine'", "simplex", "ipm", "auto") {
        di as error "_dea_tecontract: engine() must be simplex, ipm, or auto"
        exit 198
    }
    if "`scaling'" == "" local scaling "none"
    local scaling = lower("`scaling'")
    if !inlist("`scaling'", "none", "row") {
        di as error "_dea_tecontract: scaling() must be none or row"
        exit 198
    }
    if `tolerance' <= 0 {
        di as error "_dea_tecontract: tolerance() must be positive"
        exit 198
    }
    if "`slacknormalization'" == "" local slacknormalization "raw"
    local slacknormalization = lower("`slacknormalization'")
    if !inlist("`slacknormalization'", "raw", "range") {
        di as error "_dea_tecontract: slacknormalization() must be raw or range"
        exit 198
    }
    if `radialtolerance' < 0 {
        di as error "_dea_tecontract: radialtolerance() must be nonnegative"
        exit 198
    }

    confirm new variable `generate'
    if "`status'" != "" confirm new variable `status'
    marksample touse, novarlist
    markout `touse' `inputs' `outputs'
    quietly count if `touse'
    if r(N)==0 {
        di as error "_dea_tecontract: no complete observations"
        exit 2000
    }
    quietly generate double `generate' = .
    if "`status'" != "" quietly generate byte `status' = .

    mata: _dea_tecontract_solve("`inputs'", "`outputs'", "`touse'", ///
        "`generate'", "`status'", "`rts'", "`orientation'", `stage', "`engine'", ///
        "`scaling'", "`slacknormalization'", `radialtolerance', ///
        `tolerance', "`trace'")

    local N = r(N)
    local failed = r(n_failed)
    local fallbacks = r(fallbacks)
    local mean_eff = r(mean_efficiency)
    local maxeq = r(max_eq_violation)
    local maxineq = r(max_ineq_violation)
    local maxbound = r(max_bound_violation)
    local stage2failed = r(stage2_failures)
    local stage2fallbacks = r(stage2_fallbacks)
    tempname E S G R R2 RA L SI SO PX PY TS S2 SW
    matrix `E' = r(efficiency)
    matrix `S' = r(status)
    matrix `G' = r(engine_code)
    matrix `R' = r(radial)
    matrix `R2' = r(radial_stage2)
    matrix `RA' = r(radial_adjustment)
    matrix `L' = r(lambda)
    matrix `SI' = r(slack_in)
    matrix `SO' = r(slack_out)
    matrix `PX' = r(projected_inputs)
    matrix `PY' = r(projected_outputs)
    matrix `TS' = r(total_slack)
    matrix `S2' = r(stage2_status)
    matrix `SW' = r(slack_weights)

    return scalar N = `N'
    return scalar n_failed = `failed'
    return scalar fallbacks = `fallbacks'
    return scalar mean_efficiency = `mean_eff'
    return scalar max_eq_violation = `maxeq'
    return scalar max_ineq_violation = `maxineq'
    return scalar max_bound_violation = `maxbound'
    return scalar stage2_failures = `stage2failed'
    return scalar stage2_fallbacks = `stage2fallbacks'
    return matrix efficiency = `E'
    return matrix status = `S'
    return matrix engine_code = `G'
    return matrix radial = `R'
    return matrix radial_stage2 = `R2'
    return matrix radial_adjustment = `RA'
    return matrix lambda = `L'
    return matrix slack_in = `SI'
    return matrix slack_out = `SO'
    return matrix projected_inputs = `PX'
    return matrix projected_outputs = `PY'
    return matrix total_slack = `TS'
    return matrix stage2_status = `S2'
    return matrix slack_weights = `SW'
    return local cmd "_dea_tecontract"
    return local inputs "`inputs'"
    return local outputs "`outputs'"
    return local rts "`rts'"
    return local orientation "`orientation'"
    return local model "te"
    return scalar stage = `stage'
    return local engine_requested "`engine'"
    return local scaling "`scaling'"
    return local slack_normalization "`slacknormalization'"
    return scalar radial_tolerance = `radialtolerance'
    return local generate "`generate'"
    return local statusvar "`status'"
end

capture mata: mata drop _dea_tecontract_solve()
mata:
mata set matastrict on

void _dea_tecontract_solve(string scalar inputs, string scalar outputs,
    string scalar touse, string scalar effvar, string scalar statusvar,
    string scalar rts, string scalar orientation, real scalar stage,
    string scalar engine, string scalar scaling, string scalar slacknorm,
    real scalar radialtol, real scalar tol, string scalar trace)
{
    real matrix X, Y, Aub, Aeq, Lambda, Sin, Sout, Xproj, Yproj, Aeq2, Aub2
    real colvector obs, x0, y0, bub, beq, eff, statuses, engines, radials, lam
    real colvector beq2, bub2, total_slack, stage2statuses, radial2, radialadj
    real rowvector c, lb, ub, diag, c2, lb2, ub2
    real scalar n, m, s, k, nv, nv2, fallback, fallback2, objscale, radial, radialtry
    real rowvector slackweights, inputscale, outputscale
    real scalar rowscalemin, rowscalemax, maxeq, maxineq, maxbound
    real scalar fallbacks, failed, stage2failed, stage2fallbacks
    struct _DEA_LPResult scalar out, out2

    X=st_data(.,tokens(inputs),touse)
    Y=st_data(.,tokens(outputs),touse)
    obs=selectindex(st_data(.,touse))
    n=rows(X); m=cols(X); s=cols(Y); nv=n+1
    eff=J(n,1,.); radials=J(n,1,.); statuses=J(n,1,3); engines=J(n,1,0)
    Lambda=J(n,n,.); Sin=J(n,m,.); Sout=J(n,s,.)
    Xproj=J(n,m,.); Yproj=J(n,s,.)
    total_slack=J(n,1,.); stage2statuses=J(n,1,.)
    radial2=J(n,1,.); radialadj=J(n,1,.)
    if (slacknorm=="range") {
        inputscale=colmax(abs(X)); inputscale=inputscale+(inputscale:==0)
        outputscale=colmax(abs(Y)); outputscale=outputscale+(outputscale:==0)
        slackweights=(1:/inputscale,1:/outputscale)
    }
    else slackweights=J(1,m+s,1)
    fallbacks=0; failed=0; stage2failed=0; stage2fallbacks=0
    maxeq=0; maxineq=0; maxbound=0

    for (k=1; k<=n; k++) {
        x0=X[k,.]'; y0=Y[k,.]'
        if (orientation=="input") {
            c=(1,J(1,n,0)); lb=J(1,nv,0); ub=J(1,nv,.)
            Aub=(-x0,X') \ (J(s,1,0),-Y')
            bub=J(m,1,0) \ (-y0)
        }
        else {
            c=(-1,J(1,n,0)); lb=(1,J(1,n,0)); ub=J(1,nv,.)
            Aub=(J(m,1,0),X') \ (y0,-Y')
            bub=x0 \ J(s,1,0)
        }
        if (rts=="vrs") {
            Aeq=(0,J(1,n,1)); beq=1
        }
        else {
            Aeq=J(0,nv,0); beq=J(0,1,0)
        }
        if (rts=="drs") {
            Aub=Aub \ (0,J(1,n,1)); bub=bub \ 1
        }
        else if (rts=="irs") {
            Aub=Aub \ (0,J(1,n,-1)); bub=bub \ (-1)
        }
        out=_dea_lpcontract_dispatch(c,"min",Aeq,beq,Aub,bub,lb,ub,
            engine,tol,scaling,0,trace,fallback,objscale,
            rowscalemin,rowscalemax)
        /* A radial DEA technology always contains the evaluated DMU, so the
           stage-1 primal has a feasible reference solution.  Under auto,
           verify a simplex infeasible/unbounded classification with IPM
           before treating it as terminal.  Keep this DEA-specific invariant
           out of the generic LP adapter, where either status may be genuine. */
        if (engine=="auto" & (out.status==1 | out.status==2)) {
            out=_dea_lpcontract_dispatch(c,"min",Aeq,beq,Aub,bub,lb,ub,
                "ipm",tol,scaling,0,trace,fallback2,objscale,
                rowscalemin,rowscalemax)
            fallback=1
        }
        statuses[k]=out.status
        engines[k]=(out.engine=="simplex" ? 1 : (out.engine=="ipm" ? 2 : 0))
        fallbacks=fallbacks+fallback
        if (out.status==0) {
            radial=out.x[1]
            lam=out.x[2..nv]'
            radials[k]=radial
            if (stage==2) {
                nv2=n+m+s
                c2=(J(1,n,0),-slackweights)
                lb2=J(1,nv2,0); ub2=J(1,nv2,.)
                Aeq2=(X',I(m),J(m,s,0)) \
                     (Y',J(s,m,0),-I(s))
                beq2=(orientation=="input" ? radial*x0 : x0) \
                     (orientation=="input" ? y0 : radial*y0)
                Aub2=J(0,nv2,0); bub2=J(0,1,0)
                if (rts=="vrs") {
                    Aeq2=Aeq2 \ (J(1,n,1),J(1,m+s,0)); beq2=beq2 \ 1
                }
                else if (rts=="drs") {
                    Aub2=(J(1,n,1),J(1,m+s,0)); bub2=1
                }
                else if (rts=="irs") {
                    Aub2=(J(1,n,-1),J(1,m+s,0)); bub2=-1
                }
                out2=_dea_lpcontract_dispatch(c2,"min",Aeq2,beq2,Aub2,bub2,
                    lb2,ub2,engine,tol,scaling,0,trace,fallback2,objscale,
                    rowscalemin,rowscalemax)
                /* The fixed-radial stage-2 program is feasible whenever the
                   accepted stage-1 solution is feasible.  Verify suspect
                   simplex terminal classifications under auto as well. */
                if (engine=="auto" & (out2.status==1 | out2.status==2)) {
                    out2=_dea_lpcontract_dispatch(c2,"min",Aeq2,beq2,Aub2,bub2,
                        lb2,ub2,"ipm",tol,scaling,0,trace,fallback,objscale,
                        rowscalemin,rowscalemax)
                    fallback2=1
                }
                /* A barrier solution may report the radial optimum a few
                   digits beyond the feasible optimal face.  Retry only after
                   failure, at the conservative endpoint of the declared
                   band, and report the adjustment explicitly. */
                if (out2.status!=0 & radialtol>0) {
                    radialtry=(orientation=="input" ?
                        radial+radialtol*(1+abs(radial)) :
                        max((1,radial-radialtol*(1+abs(radial)))))
                    beq2=(orientation=="input" ? radialtry*x0 : x0) \
                         (orientation=="input" ? y0 : radialtry*y0)
                    if (rts=="vrs") beq2=beq2\1
                    out2=_dea_lpcontract_dispatch(c2,"min",Aeq2,beq2,Aub2,bub2,
                        lb2,ub2,engine,tol,scaling,0,trace,fallback2,objscale,
                        rowscalemin,rowscalemax)
                    if (out2.status==0) radial=radialtry
                }
                stage2statuses[k]=out2.status
                fallbacks=fallbacks+fallback2
                stage2fallbacks=stage2fallbacks+fallback2
                if (out2.status==0) {
                    lam=out2.x[1..n]'
                    Sin[k,.]=out2.x[(n+1)..(n+m)]
                    Sout[k,.]=out2.x[(n+m+1)..nv2]
                    total_slack[k]=-out2.objective
                    radial2[k]=radial
                    radialadj[k]=radial-radials[k]
                    engines[k]=(out2.engine=="simplex" ? 1 :
                        (out2.engine=="ipm" ? 2 : 0))
                }
                else {
                    statuses[k]=out2.status
                    stage2failed++; failed++
                    continue
                }
            }
            else {
                stage2statuses[k]=.
                radial2[k]=radial
                radialadj[k]=0
            }
            Lambda[k,.]=lam'
            Xproj[k,.]=lam'*X
            Yproj[k,.]=lam'*Y
            if (stage==1 & orientation=="input") {
                Sin[k,.]=(radial*x0-X'*lam)'
                Sout[k,.]=(Y'*lam-y0)'
            }
            else if (stage==1) {
                Sin[k,.]=(x0-X'*lam)'
                Sout[k,.]=(Y'*lam-radial*y0)'
            }
            Sin[k,.]=Sin[k,.]:*(abs(Sin[k,.]):>10*tol)
            Sout[k,.]=Sout[k,.]:*(abs(Sout[k,.]):>10*tol)
            eff[k]=(orientation=="input" ? radial : 1/radial)
            diag=_dea_lpcontract_diagnostics(out.x,Aeq,beq,Aub,bub,lb,ub)
            maxeq=max((maxeq,diag[1])); maxineq=max((maxineq,diag[2]))
            maxbound=max((maxbound,diag[3]))
        }
        else failed++
    }

    st_store(obs,st_varindex(effvar),eff)
    if (statusvar!="") st_store(obs,st_varindex(statusvar),statuses)
    st_numscalar("r(N)",n); st_numscalar("r(n_failed)",failed)
    st_numscalar("r(fallbacks)",fallbacks)
    st_numscalar("r(mean_efficiency)",(failed<n ? mean(select(eff,eff:<.)) : .))
    st_numscalar("r(max_eq_violation)",maxeq)
    st_numscalar("r(max_ineq_violation)",maxineq)
    st_numscalar("r(max_bound_violation)",maxbound)
    st_numscalar("r(stage2_failures)",stage2failed)
    st_numscalar("r(stage2_fallbacks)",stage2fallbacks)
    st_matrix("r(efficiency)",eff)
    st_matrix("r(status)",statuses)
    st_matrix("r(engine_code)",engines)
    st_matrix("r(radial)",radials)
    st_matrix("r(radial_stage2)",radial2)
    st_matrix("r(radial_adjustment)",radialadj)
    st_matrix("r(lambda)",Lambda)
    st_matrix("r(slack_in)",Sin)
    st_matrix("r(slack_out)",Sout)
    st_matrix("r(projected_inputs)",Xproj)
    st_matrix("r(projected_outputs)",Yproj)
    st_matrix("r(total_slack)",total_slack)
    st_matrix("r(stage2_status)",stage2statuses)
    st_matrix("r(slack_weights)",slackweights)
}
end
