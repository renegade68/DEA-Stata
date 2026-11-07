*! version 1.0.0  03nov2026
*! Cooper hybrid VRS congestion: lexicographic BCC, then slack decomposition.
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
capture mata: mata drop _dea_cooper_ready()
capture mata: mata drop _dea_cooper_lp()
capture mata: mata drop _dea_cooper_range()
capture mata: mata drop _dea_cooper_compute()
mata:
mata set matastrict on
real scalar _dea_cooper_ready() return(1)

// Maximize an unchanged LP. Project only tolerance-size bound roundoff, then
// recheck the original rows AND objective before passing the point downstream.
// adj counts changed coordinates; no RHS relaxation or objective tradeoff.
real scalar _dea_cooper_lp(real rowvector c,real matrix E,real colvector b,
    real rowvector lo,real rowvector hi,real scalar tol,
    real rowvector z,real scalar adj, | real matrix U,real colvector ub,real scalar retry)
{
    class LinearProgram scalar S
    real rowvector raw
    struct _DEA_NativeRetry scalar recovered
    real scalar j,eps
    eps=min((tol,1e-9));adj=0;retry=0
    if(args()<9) {
        U=J(0,cols(c),0);ub=J(0,1,0)
    }
    if(hasmissing(c) | hasmissing(E) | hasmissing(b) | hasmissing(U) | hasmissing(ub)) return(1)
    S=LinearProgram();S.setMaxOrMin("max");S.setCoefficients(c)
    S.setEquality(E,b);if(rows(U)) S.setInequality(U,ub)
    S.setBounds(lo,hi);S.setTrace("off")
    if (_dea_nonte_optimize(S,eps,1e-12)) {
        // Bounded, algebraically equivalent native retries. Original rows and
        // objective remain the acceptance contract; no solver/model fallback.
        retry=1
        recovered=_dea_native_retry(c,"max",E,b,U,ub,lo,hi,min((eps,1e-12)),"")
        if(recovered.status!=0) return(1)
        raw=recovered.x;retry=2
    }
    else raw=S.parameters()
    z=raw
    for(j=1;j<=cols(c);j++) {
        if(lo[j]<. & z[j]<lo[j]) z[j]=lo[j]
        if(hi[j]<. & z[j]>hi[j]) z[j]=hi[j]
    }
    adj=sum(z:!=raw)
    if(!_dea_lpcheck_ok(_dea_lpcheck(z,E,b,U,ub,lo,hi),eps)) return(1)
    if(abs(c*(z-raw)')>eps*(1+abs(c*raw'))) return(1)
    return(0)
}

// Range of linear functionals over a fixed optimal face. Each row of Q is one
// functional. Missing endpoints are explicit; these LPs never replace a score.
real matrix _dea_cooper_range(real matrix Q,real matrix E,real colvector b,
    real rowvector lo,real rowvector hi,real scalar tol, | real matrix U,real colvector ub)
{
    real matrix R
    real rowvector z
    real scalar j,adj
    if(args()<7) {
        U=J(0,cols(Q),0);ub=J(0,1,0)
    }
    R=J(2,rows(Q),.)
    for(j=1;j<=rows(Q);j++) {
        if(!_dea_cooper_lp(-Q[j,.],E,b,lo,hi,tol,z,adj,U,ub)) R[1,j]=Q[j,.]*z'
        if(!_dea_cooper_lp(Q[j,.],E,b,lo,hi,tol,z,adj,U,ub)) R[2,j]=Q[j,.]*z'
    }
    return(R)
}

void _dea_cooper_compute(string scalar inputs,string scalar outputs,
    string scalar touse,string scalar resultvars,string scalar weight,
    real scalar tol,real scalar matrices,real scalar ranges,real scalar debug,
    string scalar destinations)
{
    real matrix X,Y,xs,ys,R,SI,TI,CI,SO,TX,TY,L1,L2,ST,PM,PMX,CM,CMX,RS
    real matrix E,F,Q,rr,U,meta
    real rowvector sx,sy,c,lo,hi,z,w,sm,sp,delta,tx,ty,c2,lo2,hi2,z2
    real colvector obs,b,bb
    real scalar n,m,s,v,i,phi,adj,eps,multi1,multi2,failedrange,eff,change,retry,rc
    real rowvector rawz
    string rowvector dest
    X=st_data(.,tokens(inputs),touse);Y=st_data(.,tokens(outputs),touse)
    obs=selectindex(st_data(.,touse));n=rows(X);m=cols(X);s=cols(Y);v=n+m+s+1
    sx=colmax(X);sy=colmax(Y);xs=X:/sx;ys=Y:/sy;eps=min((tol,1e-9))
    R=J(n,6,.);ST=J(n,7,.);ST[,4..7]=J(n,4,0)
    SI=TI=CI=TX=J(n,m,.);SO=TY=J(n,s,.)
    L1=L2=J(matrices ? n : 0,matrices ? n : 0,.)
    PM=PMX=J(ranges ? n : 0,m+s,.);CM=CMX=J(ranges ? n : 0,m,.)
    RS=J(ranges ? n : 0,2,.)
    eff=multi1=multi2=failedrange=0
    for(i=1;i<=n;i++) {
        // Phase 1a: (lambda,sminus,splus,phi). Self observation included.
        E=(xs',I(m),J(m,s,0),J(m,1,0))\
          (ys',J(s,m,0),-I(s),-ys[i,.]')\
          (J(1,n,1),J(1,m+s+1,0))
        b=xs[i,.]'\J(s,1,0)\1
        lo=(J(1,n+m+s,0),1);hi=(J(1,n,1),xs[i,.],J(1,s+1,.))
        c=(J(1,v-1,0),1);R[i,6]=1;ST[i,1]=1
        rc=_dea_cooper_lp(c,E,b,lo,hi,tol,z,adj,J(0,v,0),J(0,1,0),retry)
        ST[i,6]=ST[i,6]+(retry>0);ST[i,7]=ST[i,7]+(retry==2 & !rc)
        if(rc) continue
        ST[i,1]=0;ST[i,4]=ST[i,4]+adj;phi=z[v];R[i,2]=phi;R[i,3]=1/phi
        // Phase 1b: fix phi, maximize weighted slacks (no epsilon blending).
        w=(sx,sy)
        if(weight=="normalized") w=1:/(xs[i,.],ys[i,.])
        w=w:/max(w);c=(J(1,n,0),w,0);lo[v]=hi[v]=phi
        R[i,6]=2;ST[i,2]=1
        rc=_dea_cooper_lp(c,E,b,lo,hi,tol,z,adj,J(0,v,0),J(0,1,0),retry)
        ST[i,6]=ST[i,6]+(retry>0);ST[i,7]=ST[i,7]+(retry==2 & !rc)
        if(rc) continue
        ST[i,4]=ST[i,4]+adj
        // Reconstruct a coherent peer projection from a simplex-normalized
        // lambda. Recompute the slacks; reject if the original phase-1 LP or
        // its objective changes beyond the acceptance envelope. Tiny negative
        // slacks are clipped only after this full-candidate check below.
        rawz=z;z[1..n]=z[1..n]:/sum(z[1..n])
        tx=z[1..n]*xs;ty=z[1..n]*ys
        sm=xs[i,.]-tx;sp=ty-phi*ys[i,.]
        sm=rowmax((sm' , J(m,1,0)))';sp=rowmax((sp',J(s,1,0)))'
        z[n+1..n+m]=sm;z[n+m+1..n+m+s]=sp
        change=max(abs(z-rawz));ST[i,5]=change
        if(!_dea_lpcheck_ok(_dea_lpcheck(z,E,b,J(0,v,0),J(0,1,0),lo,hi),eps)) continue
        if(abs(c*(z-rawz)')>eps*(1+abs(c*rawz'))) continue
        ST[i,2]=0;ST[i,4]=ST[i,4]+sum(z:!=rawz)
        tx=xs[i,.]-sm
        // Y lambda is the chosen output projection; its difference from
        // phi*y0+splus has just been checked against the original equality.
        SI[i,.]=sm:*sx;SO[i,.]=sp:*sy;TX[i,.]=tx:*sx;TY[i,.]=ty:*sy
        R[i,4]=mean((sm:/xs[i,.])')
        if(matrices) L1[i,.]=z[1..n]
        if(abs(phi-1)<=eps & max((sm:/xs[i,.],sp:/ys[i,.]))<=eps) eff++
        if(ranges) {
            Q=(xs',J(m,m+s+1,0))\(ys',J(s,m+s+1,0))
            rr=_dea_cooper_range(Q,E\c,b\(c*z'),lo,hi,tol)
            PM[i,.]=rr[1,.]:*(sx,sy);PMX[i,.]=rr[2,.]:*(sx,sy)
            RS[i,1]=hasmissing(rr)
            if(!RS[i,1]) multi1=multi1+(max((rr[2,.]-rr[1,.]):/(xs[i,.],ys[i,.]))>10*eps)
        }
        // Phase 2: maximize technical slack at the SELECTED phase-1 projection.
        // Eliminate delta=X'mu-tx. The same LP now has only n variables:
        // tx<=X'mu<=x0, Y'mu=ty, sum(mu)=1. This avoids redundant tiny
        // slack bounds at efficient points; no inequality/RHS is relaxed.
        F=ys'\J(1,n,1);bb=ty'\1;U=xs'\(-xs')
        lo2=J(1,n,0);hi2=J(1,n,1)
        w=sx;if(weight=="normalized") w=1:/xs[i,.]
        w=w:/max(w);c2=(xs*w')'
        R[i,6]=3;ST[i,3]=1
        rc=_dea_cooper_lp(c2,F,bb,lo2,hi2,tol,z2,adj,U,xs[i,.]'\(-tx'),retry)
        ST[i,6]=ST[i,6]+(retry>0);ST[i,7]=ST[i,7]+(retry==2 & !rc)
        if(rc) continue
        delta=z2*xs-tx
        // Recovered technical slacks must also satisfy their original bounds.
        rawz=delta;delta=rowmax((delta',J(m,1,0)))'
        delta=rowmin((delta',sm'))'
        if(max(abs(delta-rawz))>eps) continue
        if(abs(w*(delta-rawz)')>eps*(1+abs(w*rawz'))) continue
        if(max(abs(z2*xs-delta-tx))>eps) continue
        ST[i,3]=0;ST[i,4]=ST[i,4]+adj+sum(delta:!=rawz)
        TI[i,.]=delta:*sx;CI[i,.]=(sm-delta):*sx
        R[i,1]=mean(((sm-delta):/xs[i,.])');R[i,5]=mean((delta:/xs[i,.])');R[i,6]=0
        if(matrices) L2[i,.]=z2[1..n]
        if(ranges) {
            Q=xs'
            rr=_dea_cooper_range(Q,F\c2,bb\(c2*z2'),lo2,hi2,tol,U,xs[i,.]'\(-tx'))
            CM[i,.]=(xs[i,.]-rr[2,.]):*sx;CMX[i,.]=(xs[i,.]-rr[1,.]):*sx
            RS[i,2]=hasmissing(rr)
            if(!RS[i,2]) multi2=multi2+(max((rr[2,.]-rr[1,.]):/xs[i,.])>10*eps)
        }
    }
    if(ranges) failedrange=sum(rowsum(RS:==1):>0)
    if(debug) for(i=1;i<=n;i++) printf("  DMU observation %g: status=%g, adjusted coordinates=%g\n",obs[i],R[i,6],ST[i,4])
    st_store(obs,tokens(resultvars),R)
    dest=tokens(destinations)
    meta=(n,sum(R[,6]:!=0),eff,sum((R[,1]:>eps):&(R[,1]:<.)),sum(ST[,4]),multi1,multi2,failedrange,
        (ranges ? sum(rowsum(RS:==0):==2) : .))\
        (sum(ST[,6]),sum(ST[,7]),J(1,7,.))
    st_matrix(dest[1],meta)
    if(matrices) {
        st_matrix(dest[2],SI);st_matrix(dest[3],TI);st_matrix(dest[4],CI)
        st_matrix(dest[5],SO);st_matrix(dest[6],TX);st_matrix(dest[7],TY)
        st_matrix(dest[8],L1);st_matrix(dest[9],L2);st_matrix(dest[10],ST);st_matrix(dest[11],obs)
    }
    if(ranges) {
        st_matrix(dest[12],PM);st_matrix(dest[13],PMX);st_matrix(dest[14],CM)
        st_matrix(dest[15],CMX);st_matrix(dest[16],RS)
    }
}
end
