*! version 2.3.0  02nov2026
*! Validated scale-free LPs; complete Nash iterations; matrix/detail and diagnostics.
*! 2026-11-01: honor if/in for evaluation and reference samples; reject empty samples.
*! 2026-11-01: preflight output names; explicit replacevars; protect model data.
*! dea_crossgame.ado — Cross-Efficiency Game DEA
*! Includes: Aggressive, Benevolent, Nash Equilibrium strategies
*! Requires: Stata 16+ (LinearProgram)
*!
*! AUTHORS:
*! --------
*!   Choonjoo Lee (corresponding author)
*!     Korea National Defense University
*!     Seoul, Korea
*!     bloom.rampike@gmail.com
*!
*! REFERENCES:
*! -----------
*!   Doyle, J., & Green, R. (1994). Efficiency and cross-efficiency in
*!     DEA: Derivations, meanings and uses. Journal of the Operational
*!     Research Society, 45(5), 567-578. doi:10.1057/jors.1994.84
*!
*!   Liang, L., Wu, J., Cook, W. D., & Zhu, J. (2008). The DEA game
*!     cross-efficiency model and its Nash equilibrium. Operations
*!     Research, 56(5), 1278-1288. doi:10.1287/opre.1070.0487
*!
*!   Wu, J., Liang, L., & Yang, F. (2009). Achievement and benchmarking
*!     of countries at the Summer Olympics using cross efficiency
*!     evaluation method. European Journal of Operational Research,
*!     197(2), 722-730. doi:10.1016/j.ejor.2008.06.030
*!
*! Changelog:
*!   2.2.0  03may2026  Implement Liang et al. (2008) model (5) for strategy(nash)
*!   2.1.0  28apr2026  Add Form A syntax (x1 x2 = y1 y2)

/**
 * CROSS-EFFICIENCY GAME DEA
 * =========================
 *
 * Models DMU competition/cooperation through cross-efficiency evaluation.
 *
 * Strategies:
 *   - Aggressive: Minimize the weighted-output secondary surrogate
 *   - Benevolent: Maximize the weighted-output secondary surrogate
 *   - Nash: Iterate model (5), with explicit numerical stopping status
 *
 * References:
 *   - Doyle & Green (1994): Cross-efficiency derivations, meanings, and uses
 *   - Liang et al. (2008): DEA game cross-efficiency model and Nash equilibrium
 *   - Wu et al. (2009): Cross-efficiency benchmarking application
 */

capture program drop dea_crossgame
program define dea_crossgame, rclass
    version 16.0

    // =========================================================================
    // Unified Syntax Parser (Form A / Form B)
    //   Form A: dea_crossgame x1 x2 = y1 y2 [if] [in], [options]
    //   Form B: dea_crossgame [if] [in], inputs(x1 x2) outputs(y1 y2) [options]
    // =========================================================================
    _dea_parse_io `0'
    local 0 `"`r(args)'"'

    syntax [if] [in], ///
        Inputs(varlist numeric) ///
        Outputs(varlist numeric) ///
        [DMU(varname) ///
         STRAtegy(string) ///
         MAXiter(integer 100) ///
         GENerate(string) ///
         REPLACEVARS      ///
         TOLerance(real 1e-6) ///
         MATrix ///
         DETail ///
         DEBUG]

    if missing(`tolerance') | `tolerance'<=0 {
        di as error "tolerance() must be positive and finite"
        exit 198
    }
    if missing(`maxiter') | `maxiter'<1 {
        di as error "maxiter() must be a positive integer"
        exit 198
    }
    capture mata: assert(_dea_nonte_lp_ready()==1)
    if _rc {
        quietly findfile _dea_nonte_lp.ado
        quietly run "`r(fn)'"
    }

    // Defaults
    if "`strategy'" == "" local strategy "all"
    if "`generate'" == "" local generate "_cg"

    // Validate strategy
    if !inlist("`strategy'", "aggressive", "benevolent", "nash", "all") {
        di as error "strategy() must be aggressive, benevolent, nash, or all"
        exit 198
    }

    local _planned ""
    if inlist("`strategy'", "all", "aggressive") local _planned `_planned' `generate'_agg
    if inlist("`strategy'", "all", "benevolent") local _planned `_planned' `generate'_ben
    if inlist("`strategy'", "all", "nash") local _planned `_planned' `generate'_nash
    // Validate every output before any existing variable can be replaced.
    _dea_output_check, names(`_planned') ///
        protect(`inputs' `outputs' `dmu') `replacevars'
    local _replace_outputs "`r(existing)'"

    // Mark sample — do NOT include dmu in markout (string DMU causes r(109))
    marksample touse, novarlist
    markout `touse' `inputs' `outputs'
    if "`dmu'"!="" quietly replace `touse'=0 if missing(`dmu')

    // Count
    qui count if `touse'
    local n = r(N)
    if `n' == 0 {
        di as error "No observations in estimation sample"
        exit 2000
    }
    // Aggressive/benevolent averages exclude self; all requests both of them.
    // Check the marked sample before IDs, output replacement, or any LP solve.
    if `n' < 2 & inlist("`strategy'", "aggressive", "benevolent", "all") {
        di as error "strategy(`strategy') requires at least 2 DMUs after if/in and missing-value exclusions"
        exit 2001
    }

    foreach _side in inputs outputs {
        tempvar _positive
        quietly generate byte `_positive'=0 if `touse'
        foreach _v of local `_side' {
            quietly count if `touse' & `_v'<0
            if r(N) {
                di as error "crossgame quantities must be nonnegative"
                exit 459
            }
            quietly replace `_positive'=1 if `touse' & `_v'>0
        }
        quietly count if `touse' & `_positive'==0
        if r(N) {
            di as error "each selected observation requires a positive input and output"
            exit 459
        }
    }

    // Display header
    di as text ""
    di as text "{hline 70}"
    di as text "Cross-Efficiency Game DEA"
    di as text "{hline 70}"
    di as text "Strategy:      " as result "`strategy'"
    di as text "DMUs:          " as result "`n'"
    di as text "Max iterations:" as result "`maxiter'" as text " (for Nash)"
    di as text "{hline 70}"

    // Create DMU identifier if not provided
    if "`dmu'" == "" {
        tempvar dmu_temp
        qui gen `dmu_temp' = _n if `touse'
        local dmu "`dmu_temp'"
    }

    // The sample may refer to an existing output; mark it before replacement.
    foreach _ov of local _replace_outputs {
        quietly drop `_ov'
    }

    // Create output variables
    if "`strategy'" == "all" | "`strategy'" == "aggressive" {
        qui gen double `generate'_agg = . if `touse'
    }
    if "`strategy'" == "all" | "`strategy'" == "benevolent" {
        qui gen double `generate'_ben = . if `touse'
    }
    if "`strategy'" == "all" | "`strategy'" == "nash" {
        qui gen double `generate'_nash = . if `touse'
    }

    // Call Mata function; tempname isolates all optional matrix results.
    tempname cg
    mata: _dea_crossgame_compute("`inputs'", "`outputs'", ///
        "`dmu'", "`touse'", "`generate'", ///
        "`strategy'", `maxiter', `tolerance', ///
        "`matrix'" != "", "`detail'" != "", "`debug'" != "", "`cg'")

    di as text ""
    di as text "{hline 70}"

    // Return values
    return local strategy "`strategy'"
    return scalar n = `n'
    return scalar tolerance=`tolerance'
    return scalar maxiter=`maxiter'
    return scalar n_failed_primary=`cg'_meta[1,1]
    if inlist("`strategy'","all","aggressive") {
        return scalar n_failed_agg=`cg'_meta[1,2]
        return scalar n_missing_agg=`cg'_meta[1,10]
    }
    if inlist("`strategy'","all","benevolent") {
        return scalar n_failed_ben=`cg'_meta[1,3]
        return scalar n_missing_ben=`cg'_meta[1,11]
    }
    if inlist("`strategy'","all","nash") {
        return scalar nash_status=`cg'_meta[1,4]
        return scalar nash_converged=(`cg'_meta[1,4]==0)
        return scalar nash_iterations=`cg'_meta[1,5]
        return scalar nash_attempts=`cg'_meta[1,6]
        return scalar nash_failed_pairs=`cg'_meta[1,7]
        return scalar nash_last_change=`cg'_meta[1,8]
        return scalar nash_initial_missing=`cg'_meta[1,9]
    }
    if "`matrix'"!="" {
        matrix colnames `cg'_v=`inputs'
        matrix colnames `cg'_u=`outputs'
        return matrix self_eff=`cg'_self
        return matrix weights_input=`cg'_v
        return matrix weights_output=`cg'_u
        return matrix ce_initial=`cg'_initial
        if inlist("`strategy'","all","aggressive") {
            matrix colnames `cg'_va=`inputs'
            matrix colnames `cg'_ua=`outputs'
            return matrix ce_agg=`cg'_agg
            return matrix weights_agg_input=`cg'_va
            return matrix weights_agg_output=`cg'_ua
            return matrix status_agg=`cg'_sa
            return matrix count_agg=`cg'_pa
        }
        if inlist("`strategy'","all","benevolent") {
            matrix colnames `cg'_vb=`inputs'
            matrix colnames `cg'_ub=`outputs'
            return matrix ce_ben=`cg'_ben
            return matrix weights_ben_input=`cg'_vb
            return matrix weights_ben_output=`cg'_ub
            return matrix status_ben=`cg'_sb
            return matrix count_ben=`cg'_pb
        }
        if inlist("`strategy'","all","nash") {
            matrix colnames `cg'_trace=iteration max_change missing_pairs complete
            return matrix ce_nash=`cg'_nash
            return matrix nash_history=`cg'_history
            return matrix nash_trace=`cg'_trace
        }
    }
    // Mata-created suffixes are not individually registered Stata tempnames.
    // Return copies above, then remove every private workspace matrix.
    foreach _m in meta self v u initial agg va ua sa pa ben vb ub sb pb nash history trace {
        capture matrix drop `cg'_`_m'
    }
end


// Reloads must replace old helper definitions in the same Stata session.
capture mata: mata drop _dea_cg_dmu_read()
capture mata: mata drop _dea_crossgame_compute()
capture mata: mata drop _dea_cg_multiplier()
capture mata: mata drop _dea_cg_aggressive()
capture mata: mata drop _dea_cg_benevolent()
capture mata: mata drop _dea_cg_nash()
capture mata: mata drop _dea_cg_lp()
capture mata: mata drop _dea_cg_peers()
capture mata: mata drop _dea_cg_average()
capture mata: mata drop _dea_cg_show()

mata:
mata set matastrict on
// mode: 0 primary, 1 aggressive numerator surrogate, 2 benevolent surrogate,
    // 3 Liang et al. (2008) model (5), with protected peer d and score theta.
real scalar _dea_cg_lp(real matrix X0, real matrix Y0, real scalar j,
    real scalar mode, real scalar own, real scalar d, real scalar theta,
    real scalar tol, real rowvector v, real rowvector u, real scalar value)
{
    real matrix X,Y,E,A
    real rowvector sx,sy,c,sol
    real colvector b,rs
    real scalar k,m,s,rc
    class LinearProgram scalar S
    m=cols(X0);s=cols(Y0);sx=X0[j,.];sy=Y0[j,.]
    v=J(1,m,.);u=J(1,s,.);value=.
    for(k=1;k<=m;k++) if(sx[k]<=0) sx[k]=1
    for(k=1;k<=s;k++) if(sy[k]<=0) sy[k]=1
    X=X0:/sx;Y=Y0:/sy
    if(hasmissing(X) | hasmissing(Y)) return(1)
    E=(X[j,.],J(1,s,0));b=1
    A=(-X,Y)
    if(mode==1 | mode==2) {
        if(own>=.) return(1)
        E=E\(J(1,m,0),Y[j,.]);b=b\own
        c=(J(1,m,0),colsum(Y)-Y[j,.])
        if(mode==2) c=-c
    }
    else c=(J(1,m,0),-Y[j,.])
    if(mode==3) {
        if(theta>=. | theta<0) return(1)
        A=A\(theta*X[d,.],-Y[d,.])
    }
    rs=rowmax(abs(A))
    for(k=1;k<=rows(rs);k++) if(rs[k]==0) rs[k]=1
    A=A:/rs
    S=LinearProgram();S.setMaxOrMin("min");S.setCoefficients(c)
    S.setEquality(E,b);S.setInequality(A,J(rows(A),1,0))
    S.setBounds(J(1,m+s,0),J(1,m+s,.));S.setTrace("off")
    rc=_dea_nonte_optimize(S,tol,1e-12)
    if(rc) return(1)
    sol=S.parameters();value=S.value()
    if(mode==0 | mode==3) {
        value=-value
        if(value<0 | value>1+tol | value>=.) {
            value=.;return(1)
        }
    }
    v=sol[1..m]:/sx;u=sol[(m+1)..(m+s)]:/sy
    return(0)
}

real matrix _dea_cg_peers(real matrix X,real matrix Y,real matrix V,
    real matrix U,real scalar tol)
{
    real matrix CE
    real scalar i,j,num,den,q,n
    n=rows(X);CE=J(n,n,.)
    for(i=1;i<=n;i++) for(j=1;j<=n;j++) {
        num=U[i,.]*Y[j,.]';den=V[i,.]*X[j,.]'
        if(num>=0 & num<. & den>0 & den<.) {
            q=num/den
            if(q<. & q<=1+tol) CE[i,j]=q
        }
    }
    return(CE)
}

real colvector _dea_cg_average(real matrix CE,real scalar exclude,
    real colvector counts)
{
    real colvector scores
    real scalar n,i,j,total
    n=rows(CE);counts=J(n,1,0);scores=J(n,1,.)
    for(j=1;j<=n;j++) {
        total=0
        for(i=1;i<=n;i++) if((!exclude | i!=j) & CE[i,j]<.) {
            total=total+CE[i,j];counts[j]=counts[j]+1
        }
        if(counts[j]) scores[j]=total/counts[j]
    }
    return(scores)
}

void _dea_cg_show(real matrix CE,string colvector labels,string scalar title)
{
    real scalar i,j
    printf("\n%s (row d, column j; selected observation order)\n",title)
    printf("%12s", "d / j")
    for(j=1;j<=rows(CE);j++) printf(" %9.0f",j)
    printf("\n")
    for(i=1;i<=rows(CE);i++) {
        printf("%12s",labels[i])
        for(j=1;j<=cols(CE);j++) printf(" %9.5f",CE[i,j])
        printf("\n")
    }
}

void _dea_crossgame_compute(string scalar inputs,string scalar outputs,
    string scalar dmu_var,string scalar touse,string scalar gen_prefix,
    string scalar strategy,real scalar maxiter,real scalar tol,
    real scalar show_matrix,real scalar show_detail,real scalar show_debug,
    string scalar prefix)
{
    real matrix X,Y,V,U,Va,Ua,Vb,Ub,C0,Ca,Cb,Cn,Cnew,hist,trace
    real rowvector sx,sy,v,u,meta
    real colvector obs_indices,self,sa,sb,pa,pb,dummy,a,b,theta,next
    string colvector labels
    real scalar n,m,s,i,j,d,k,rc,value,ns,iterations,attempts,failed,delta,initmissing
    real scalar wanta,wantb,wantn
    X = st_data(., tokens(inputs), touse)
    Y=st_data(.,tokens(outputs),touse);obs_indices=selectindex(st_data(.,touse))
    if(substr(st_vartype(dmu_var),1,3)=="str") labels=st_sdata(.,dmu_var,touse)
    else labels=strtrim(strofreal(st_data(.,dmu_var,touse),"%21.17g"))
    n=rows(X);m=cols(X);s=cols(Y)
    sx=colmax(X);sy=colmax(Y)
    for(k=1;k<=m;k++) if(sx[k]==0) sx[k]=1
    for(k=1;k<=s;k++) if(sy[k]==0) sy[k]=1
    X=X:/sx;Y=Y:/sy
    wanta=(strategy=="all" | strategy=="aggressive")
    wantb=(strategy=="all" | strategy=="benevolent")
    wantn=(strategy=="all" | strategy=="nash")
    self=J(n,1,.);V=J(n,m,.);U=J(n,s,.)
    for(i=1;i<=n;i++) {
        rc=_dea_cg_lp(X,Y,i,0,.,.,.,tol,v,u,value)
        if(!rc) {
            self[i]=value;V[i,.]=v;U[i,.]=u
        }
    }
    C0=_dea_cg_peers(X,Y,V,U,tol)
    Va=V;Ua=U;Vb=V;Ub=U;sa=J(n,1,0);sb=J(n,1,0)
    a=J(n,1,.);b=a;pa=J(n,1,0);pb=pa
    for(i=1;i<=n;i++) if(self[i]<.) {
        if(wanta) {
            rc=_dea_cg_lp(X,Y,i,1,self[i],.,.,tol,v,u,value)
            sa[i]=(rc ? 3 : 2)
            if(!rc) {
                Va[i,.]=v;Ua[i,.]=u
            }
        }
        if(wantb) {
            rc=_dea_cg_lp(X,Y,i,2,self[i],.,.,tol,v,u,value)
            sb[i]=(rc ? 3 : 2)
            if(!rc) {
                Vb[i,.]=v;Ub[i,.]=u
            }
        }
    }
    if(wanta) {
        Ca=_dea_cg_peers(X,Y,Va,Ua,tol);a=_dea_cg_average(Ca,1,pa)
    }
    if(wantb) {
        Cb=_dea_cg_peers(X,Y,Vb,Ub,tol);b=_dea_cg_average(Cb,1,pb)
    }
    ns=.;iterations=0;attempts=0;failed=0;delta=.;initmissing=sum(C0:>=.)
    theta=J(n,1,.);Cn=J(n,n,.)
    if(wantn) {
        // Never run a different game by silently averaging an incomplete matrix.
        ns=2
        if(!hasmissing(C0)) {
            Cn=C0;theta=colsum(C0)'/n;ns=1
        }
        hist=theta';trace=(0,.,initmissing,ns==1)
        if(ns==1) for(k=1;k<=maxiter;k++) {
            attempts=k;Cnew=J(n,n,.)
            for(d=1;d<=n;d++) for(j=1;j<=n;j++) {
                rc=_dea_cg_lp(X,Y,j,3,.,d,theta[d],tol,v,u,value)
                if(!rc) Cnew[d,j]=value
            }
            failed=sum(Cnew:>=.)
            if(failed) {
                ns=3;trace=trace\(k,.,failed,0);break
            }
            next=colsum(Cnew)'/n;delta=max(abs(next-theta))
            iterations=k;Cn=Cnew;theta=next;hist=hist\theta'
            trace=trace\(k,delta,0,1)
            if(show_debug) printf("  Nash iteration %g: max change = %12.8g\n",k,delta)
            if(delta<tol) {
                ns=0;break
            }
        }
    }
    printf("\nCross-Efficiency Game Results (CRS, input orientation)\n")
    printf("%20s | %10s", "DMU", "Self")
    if(wanta) printf(" %10s","Aggressive")
    if(wantb) printf(" %10s","Benevolent")
    if(wantn) printf(" %10s","Nash state")
    printf("\n")
    for(i=1;i<=n;i++) {
        printf("%20s | %10.5f",labels[i],self[i])
        if(wanta) {
            printf(" %10.5f",a[i]);st_store(obs_indices[i],gen_prefix+"_agg",a[i])
        }
        if(wantb) {
            printf(" %10.5f",b[i]);st_store(obs_indices[i],gen_prefix+"_ben",b[i])
        }
        if(wantn) {
            printf(" %10.5f",theta[i]);st_store(obs_indices[i],gen_prefix+"_nash",theta[i])
        }
        printf("\n")
    }
    if(sum(self:>=.)) printf("Note: %g primary weight solve(s) unavailable.\n",sum(self:>=.))
    if(wanta & sum(sa:==3)) printf("Note: %g aggressive secondary solve(s) unavailable; primary weights retained.\n",sum(sa:==3))
    if(wantb & sum(sb:==3)) printf("Note: %g benevolent secondary solve(s) unavailable; primary weights retained.\n",sum(sb:==3))
    if(wanta & sum(pa)<n*(n-1)) printf("Note: aggressive means use available peers; %g pairs unavailable.\n",n*(n-1)-sum(pa))
    if(wantb & sum(pb)<n*(n-1)) printf("Note: benevolent means use available peers; %g pairs unavailable.\n",n*(n-1)-sum(pb))
    if(wantn) {
        printf("Nash status %g: %g complete iteration(s), %g attempted.\n",ns,iterations,attempts)
        if(ns==0) printf("Successive complete states satisfy tolerance; this is a numerical stopping test.\n")
        if(ns==1) printf("Iteration limit reached; returning the last complete state, not confirmed convergence.\n")
        if(ns==2) printf("Initialization incomplete (%g missing cells); Nash scores unavailable.\n",initmissing)
        if(ns==3) printf("Iteration failed (%g pair solves); returning the last complete state, not an equilibrium.\n",failed)
    }
    if(show_detail) {
        printf("\nWeight stages: 0=primary unavailable, 2=secondary accepted, 3=primary fallback\n")
        for(i=1;i<=n;i++) {
            printf("%20s",labels[i])
            if(wanta) printf(" | aggressive stage %g, peers %g",sa[i],pa[i])
            if(wantb) printf(" | benevolent stage %g, peers %g",sb[i],pb[i])
            printf("\n")
        }
        if(wantn) {
            printf("Nash trace: iteration, max change, missing pairs, complete state\n")
            for(i=1;i<=rows(trace);i++) printf("%6.0f %14.8g %8.0f %3.0f\n",trace[i,1],trace[i,2],trace[i,3],trace[i,4])
        }
    }
    meta=(sum(self:>=.),(wanta ? sum(sa:==3) : .),(wantb ? sum(sb:==3) : .),
    ns,iterations,attempts,failed,delta,initmissing,
    (wanta ? n*(n-1)-sum(pa) : .),(wantb ? n*(n-1)-sum(pb) : .))
    st_matrix(prefix+"_meta",meta)
    if(show_matrix) {
        st_matrix(prefix+"_self",self);st_matrix(prefix+"_v",V:/sx);st_matrix(prefix+"_u",U:/sy)
        st_matrix(prefix+"_initial",C0)
        if(wanta) {
            _dea_cg_show(Ca,labels,"Aggressive cross-efficiency matrix")
            st_matrix(prefix+"_agg",Ca);st_matrix(prefix+"_va",Va:/sx);st_matrix(prefix+"_ua",Ua:/sy)
            st_matrix(prefix+"_sa",sa);st_matrix(prefix+"_pa",pa)
        }
        if(wantb) {
            _dea_cg_show(Cb,labels,"Benevolent cross-efficiency matrix")
            st_matrix(prefix+"_ben",Cb);st_matrix(prefix+"_vb",Vb:/sx);st_matrix(prefix+"_ub",Ub:/sy)
            st_matrix(prefix+"_sb",sb);st_matrix(prefix+"_pb",pb)
        }
        if(wantn) {
            _dea_cg_show(Cn,labels,"Nash last complete matrix (protected peer d, target j)")
            st_matrix(prefix+"_nash",Cn);st_matrix(prefix+"_history",hist);st_matrix(prefix+"_trace",trace)
        }
    }
}

end
