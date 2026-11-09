*! version 1.1.0  09nov2026
*! Shared input DEA and explicit additive / historical allocation diagnostics.
version 16.0
capture mata: assert(_dea_nonte_lp_ready()==1)
if _rc {
    quietly findfile _dea_nonte_lp.ado
    quietly run "`r(fn)'"
}
capture mata: mata drop _dea_games_ready()
capture mata: mata drop _dea_games_scores()
capture mata: mata drop _dea_games_rank()
capture mata: mata drop _dea_games_coalition()
capture mata: mata drop _dea_games_allocate()
capture mata: mata drop _dea_games_compute()
mata:
mata set matastrict on
real scalar _dea_games_ready() return(12901)

// raw score, bounded score, accepted, endpoint adjustment. One LP per DMU.
real matrix _dea_games_scores(real matrix X, real matrix Y,
    string scalar rts, real scalar tol, real scalar debug, | real matrix N)
{
    real matrix xs,ys,D
    real rowvector scale,c,ub
    real scalar n,m,s,i,j,v,eps,nrc
    class LinearProgram scalar S
    n=rows(X);m=cols(X);s=cols(Y);D=J(n,4,.)
    D[,3..4]=J(n,2,0);N=J(n,4,.)
    scale=colmax(X);for(j=1;j<=m;j++) if(scale[j]==0) scale[j]=1
    xs=X:/scale
    scale=colmax(Y);for(j=1;j<=s;j++) if(scale[j]==0) scale[j]=1
    ys=Y:/scale
    c=(J(1,n,0),1);ub=(J(1,n,.),1);eps=min((tol,1e-8))
    for(i=1;i<=n;i++) {
        S=LinearProgram();S.setMaxOrMin("min");S.setCoefficients(c)
        S.setInequality((xs',-xs[i,.]')\(-ys',J(s,1,0)),J(m,1,0)\(-ys[i,.]'))
        if(rts=="vrs") S.setEquality((J(1,n,1),0),1)
        // theta<=1 is implied by the included self reference, including VRS.
        S.setBounds(J(1,n+1,0),ub);S.setTrace("off")
        // Match residual acceptance to the endpoint envelope; a loose user
        // tolerance must not accept a candidate that the score guard rejects.
        nrc=_dea_nonte_optimize(S,eps)
        N[i,.]=((nrc ? (S.converged() & S.errorcode()==0 ? 2 : 1) : 0),S.converged(),S.errorcode(),S.value())
        if(nrc==0) {
            v=S.value();D[i,1]=v
            if(v>=-eps & v<=1+eps) {
                D[i,2]=min((1,max((0,v))));D[i,3]=1
                D[i,4]=(D[i,2]!=v)
            }
        }
        if(nrc==0 & D[i,3]==0) N[i,1]=3
        if(debug) printf("  DMU row %g: accepted=%g, endpoint_adjusted=%g\n",i,D[i,3],D[i,4])
    }
    return(D)
}

real colvector _dea_games_rank(real colvector v)
{
    real colvector r
    real scalar i
    r=J(rows(v),1,.)
    if(hasmissing(v)) return(r)
    for(i=1;i<=rows(v);i++) r[i]=1+sum(v:>v[i])
    return(r)
}
real scalar _dea_games_coalition(real colvector e,real scalar k)
{
    real scalar i,v
    v=0
    for(i=1;i<=rows(e);i++) if(mod(floor(k/2^(i-1)),2)) {
        if(e[i]>=.) return(.)
        v=v+e[i]
    }
    return(v)
}

// Historical rule, not a Nash optimizer. Preserve its 0.2 transform, .01 floor,
// and upper cap; report the gain discarded by that cap instead of hiding it.
// summary: total, allocated, unallocated, mean, product, log_product,
//          zero_gain_count, capped_count, product_underflow, complete.
real matrix _dea_games_allocate(real colvector d,real scalar tol,real rowvector meta)
{
    real matrix R
    real colvector w,proposed
    real scalar n,i,G,lp
    n=rows(d);R=(d,J(n,3,.));meta=J(1,10,.)
    meta[10]=0
    for(i=1;i<=n;i++) if(d[i]<. & d[i]>=0 & d[i]<=1) {
        R[i,2]=(d[i]>=1-tol ? 1 : d[i]/(1-0.2*(1-d[i])))
    }
    if(n==0 | hasmissing(R[,1..2])) return(R)
    G=sum(R[,2]-d)
    w=rowmax((J(n,1,.01),1:-d))
    proposed=G:*w:/sum(w)
    for(i=1;i<=n;i++) {
        R[i,3]=min((1,d[i]+proposed[i]))
        R[i,4]=R[i,3]-d[i]
    }
    meta[1]=G;meta[2]=sum(R[,4]);meta[3]=G-meta[2]
    meta[4]=mean(R[,4]);meta[7]=sum(R[,4]:==0)
    meta[8]=sum((d+proposed):>1);meta[9]=0;meta[10]=1
    if(meta[7]) meta[5]=0
    else {
        lp=sum(log(R[,4]));meta[6]=lp;meta[5]=exp(lp)
        meta[9]=(meta[5]==0)
    }
    return(R)
}

void _dea_games_compute(string scalar inputs,string scalar outputs,
    string scalar dmu,string scalar touse,string scalar gen,string scalar mode,
    string scalar rts,real scalar tol,real scalar matrices,real scalar detail,
    real scalar debug,string scalar dest)
{
    real matrix X,Y,D,R,N,RQ,RR,NR,A,NA
    struct _T128_Scores scalar G
    struct _G129_Derived scalar B
    struct _C58_Work scalar W,WC
    real colvector obs,e
    string colvector labels
    string rowvector suffix
    real rowvector meta,alloc,oldalloc,out,coalitions
    real scalar n,i,k,complete,grand,shown,withheld,allocationq
    X=st_data(.,tokens(inputs),touse);Y=st_data(.,tokens(outputs),touse)
    obs=selectindex(st_data(.,touse));n=rows(X)
    if(substr(st_vartype(dmu),1,3)=="str") labels=st_sdata(.,dmu,touse)
    else labels=strtrim(strofreal(st_data(.,dmu,touse),"%21.0g"))
    D=_dea_games_scores(X,Y,rts,tol,debug,N)
    G=_t128_certify(X,Y,rts,tol,D,N);D=G.D;e=D[,2]
    complete=!hasmissing(e);grand=.;shown=0;coalitions=J(1,4,0);allocationq=0
    if(mode=="cooperative") {
        _t128_rank(G);NR=(e,e,e-e,_dea_games_rank(e));R=(e,e,e-e,G.ranks)
        RQ=(G.Q[.,7],G.Q[.,7],G.Q[.,7],J(n,1,G.rankdiag[1]))
        RR=(G.Q[.,8],G.Q[.,8],G.Q[.,8],J(n,1,G.rankdiag[2]))
        suffix=("_indiv","_shapley","_contrib","_rank")
        NA=(complete ? sum(e) : .);W=_c84_work();A=_g129_coalition(G,-1,tol,W);grand=A[3]
        if(detail & n<=15) {
            shown=1;WC=_c84_work()
            printf("\nCoalition values (mask; selected observation order; quality)\n")
            for(k=0;k<=2^n-1;k++) {
                out=_g129_coalition(G,k,tol,WC)
                coalitions[out[1]+1]=coalitions[out[1]+1]+1
                printf("%10.0f | %21.15g | %g\n",k,out[3],out[1])
            }
        }
        else if(detail) printf("Coalition listing skipped for n>15; individual scores remain available.\n")
        printf("\n%21s | %12s %12s %12s %8s\n","DMU","Individual","Shapley","Contribution","Rank")
        if(detail) printf("Certified ranks compare exact efficiencies; rounded displayed scores may coincide.\n")
    }
    else {
        NR=_dea_games_allocate(e,tol,oldalloc);NA=oldalloc[1..6]'
        B=_g129_nash(G,tol,NR,oldalloc);R=B.R;RQ=B.Q;RR=B.why;alloc=B.meta;A=B.A
        suffix=("_disagree","_coop","_bargain","_gain")
        printf("\n%21s | %12s %12s %12s %12s\n","DMU","Disagree","Scenario","Allocation","Gain")
        if(detail) printf("Historical rule uses exact decimal constants1/5 and1/100. Gains are exported before rounding: subtracting displayed allocation and disagreement need not reproduce a small gain.\n")
        if(alloc[10]) allocationq=(all(RQ[.,2..4]:==2) & all(A[1..4,1]:==2) ? 2 : 1)
        if(any(RQ[.,2..4]:==3) | any(A[1..4,1]:==3)) allocationq=3
        st_matrix(dest+"_alloc",alloc)
    }
    for(i=1;i<=n;i++) printf("%21s | %12.7g %12.7g %12.7g %12.7g\n",labels[i],R[i,1],R[i,2],R[i,3],R[i,4])
    if(any(G.Q[.,7]:==1)) printf("Unverified numerical DEA scores: %g. Related numerical results remain unverified.\n",sum(G.Q[.,7]:==1))
    withheld=sum(RQ:==3)+sum(A[.,1]:==3)+coalitions[4]
    if(withheld) printf("Warning: %g output(s) withheld after exact arithmetic or export-quality checks; other valid results are retained. Use matrix for quality/reason diagnostics.\n",withheld)
    if(!complete) printf("Unavailable DEA scores: %g. Full-sample totals, allocation or ranks may be missing.\n",sum(D[,3]:==0))
    if(mode=="cooperative") printf("Grand value: %12.7g. Additive game: Shapley=individual, contribution=0 where known.\n",grand)
    else printf("Scenario gain=%12.7g; allocated=%12.7g; unallocated after cap=%12.7g\n",alloc[1],alloc[2],alloc[3])
    for(i=1;i<=4;i++) st_store(obs,st_varindex(gen+suffix[i]),R[,i])
    meta=(sum(D[,3]:==0),sum(D[,4]),complete,grand,shown)
    st_matrix(dest+"_meta",meta);st_matrix(dest+"_aggregate",A)
    st_matrix(dest+"_summary",(sum(G.Q[.,7]:==2),sum(G.Q[.,7]:==1),sum(G.Q[.,7]:==0),sum(G.native[.,3]:==0),sum(G.native[.,4]),G.rankdiag[1..2],sum(RQ:==3),sum(A[.,1]:==3),coalitions,allocationq))
    if(matrices) {
        st_matrix(dest+"_results",R);st_matrix(dest+"_lp",D);st_matrix(dest+"_scorediag",G.Q)
        st_matrix(dest+"_native",G.native);st_matrix(dest+"_quality",RQ);st_matrix(dest+"_reason",RR)
        st_matrix(dest+"_numerical",NR);st_matrix(dest+"_numaggregate",NA)
        if(any(G.Q[.,7]:==2)) {st_matrix(dest+"_lambda",G.lambda);st_matrix(dest+"_weights",G.weights);}
    }
}

end
