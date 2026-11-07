*! Research prototype 0.1.0 03nov2026; VRS/x only, no shared solver changes.
version 16.0
capture mata: mata drop _dea_wb_pad()
capture mata: mata drop _dea_wb_dot()
capture mata: mata drop _dea_wb_bound()
capture mata: mata drop _dea_wb_solve()
capture mata: mata drop _dea_wb_ready()
capture mata: mata drop _DEA_WeakBound()
mata:
mata set matastrict on
struct _DEA_WeakBound {
    real scalar status, solver_status, upper, margin, denominator, corrections, attempts
    real rowvector prices
}
real scalar _dea_wb_ready() return(3)

// Enclose one rounded binary operation. 2^-50=8u and 2^-1022 is the
// binary64 minimum normal. This deliberately dominates rounding/underflow;
// finite interval endpoints are required. No compiler-vectorized dot product.
real rowvector _dea_wb_pad(real scalar z)
{
    real scalar radius
    real rowvector a
    if (missing(z)) return(J(1,2,.))
    radius=abs(z)*2^(-50)+2^(-1022)
    a=(z-radius,z+radius)
    if (hasmissing(a) | missing(radius)) return(J(1,2,.))
    return(a)
}
real rowvector _dea_wb_dot(real rowvector a,real rowvector b)
{
    real scalar j
    real rowvector v,p,t
    if(cols(a)!=cols(b) | hasmissing(a) | hasmissing(b)) return(J(1,2,.))
    v=(0,0)
    for(j=1;j<=cols(a);j++) {
        p=_dea_wb_pad(a[j]*b[j])
        if(hasmissing(p)) return(J(1,2,.))
        t=_dea_wb_pad(v[1]+p[1]);v[1]=t[1]
        t=_dea_wb_pad(v[2]+p[2]);v[2]=t[2]
        if(hasmissing(v)) return(J(1,2,.))
    }
    return(v)
}
// prices=(nonnegative input v, nonnegative good u, free bad h, free kappa).
// The exported stored doubles define the candidate. Free-price subtraction
// is performed before entry; feasibility is recomputed from this candidate.
struct _DEA_WeakBound scalar _dea_wb_bound(real matrix X,real matrix Y,
    real matrix B,real scalar k,real rowvector prices)
{
    struct _DEA_WeakBound scalar R
    real scalar n,m,s,q,j,iter,minimum
    real rowvector t,den,obj
    R.status=1;R.solver_status=.;R.upper=.;R.margin=.;R.denominator=.;R.corrections=0;R.attempts=0
    R.prices=prices;n=rows(X);m=cols(X);s=cols(Y);q=cols(B)
    if(n==0 | m==0 | s==0 | q==0 | rows(Y)!=n | rows(B)!=n | k<1 | k>n | k!=floor(k)) return(R)
    if(cols(prices)!=m+s+q+1 | hasmissing(X) | hasmissing(Y) | hasmissing(B) | hasmissing(prices)) return(R)
    if(any(X:<0) | any(Y:<0) | any(B:<0) | any(prices[1..m+s]:<0)) return(R)
    den=_dea_wb_dot(X[k,.],prices[1..m])
    if(hasmissing(den) | den[1]<=0) return(R)
    R.denominator=den[1]
    for(iter=0;iter<=3;iter++) {
        minimum=.
        for(j=1;j<=n;j++) {
            t=_dea_wb_dot((X[j,.],-Y[j,.],B[j,.],1),prices)
            if(hasmissing(t)) return(R)
            minimum=min((minimum,t[1]))
            t=_dea_wb_dot((X[j,.],J(1,s+q,0),1),prices)
            if(hasmissing(t)) return(R)
            minimum=min((minimum,t[1]))
        }
        if(minimum>=0 & minimum<.) break
        if(iter==3) return(R)
        // kappa is free only under VRS. Increase it, never alter any LP row.
        t=_dea_wb_pad(prices[m+s+q+1]-minimum)
        if(hasmissing(t)) return(R)
        prices[m+s+q+1]=t[2];R.corrections++
    }
    R.margin=minimum;R.prices=prices
    obj=_dea_wb_dot((X[k,.],-Y[k,.],B[k,.],1),prices)
    if(hasmissing(obj) | obj[2]<0) return(R)
    t=_dea_wb_pad(obj[2]/den[1])
    if(hasmissing(t) | t[2]<0) return(R)
    R.upper=t[2];R.status=0
    return(R)
}
struct _DEA_WeakBound scalar _dea_wb_solve(real matrix X,real matrix Y,
    real matrix B,real scalar k,real scalar tol)
{
    external real scalar _dea_rec_on
    struct _DEA_WeakBound scalar R
    struct RSM_Result scalar fit
    real scalar n,m,s,q,v,attempts
    real matrix U
    real colvector u
    real rowvector c,p
    R.status=1;R.solver_status=.;R.upper=.;R.margin=.;R.denominator=.;R.corrections=0;R.attempts=0
    R.prices=J(1,0,.)
    n=rows(X);m=cols(X);s=cols(Y);q=cols(B);v=m+s+2*q+2
    c=(X[k,.],-Y[k,.],B[k,.],-B[k,.],1,-1)
    U=(-X,Y,-B,B,J(n,1,-1),J(n,1,1))\(-X,J(n,s+2*q,0),J(n,1,-1),J(n,1,1))\(-X[k,.],J(1,s+2*q+2,0))
    u=J(2*n,1,0)\(-1)
    if (_dea_rec_on) _dea_rec_phase(1,3)
    fit=_dea_nonte_simplex(U,u,J(0,v,0),J(0,1,0),c,tol)
    attempts=1
    if(fit.status!=0 & n>1) {
        // A bounded equivalent retry changes peer row ordering only.
        U=U[(n..1,2*n..n+1,2*n+1),.]
        if (_dea_rec_on) _dea_rec_phase(1,4)
        fit=_dea_nonte_simplex(U,u,J(0,v,0),J(0,1,0),c,tol)
        attempts=2
    }
    R.solver_status=fit.status;R.attempts=attempts
    if(fit.status!=0) return(R)
    p=fit.x'
    p=(p[1..m+s],p[m+s+1..m+s+q]-p[m+s+q+1..m+s+2*q],p[v-1]-p[v])
    R=_dea_wb_bound(X,Y,B,k,p);R.solver_status=fit.status;R.attempts=attempts
    return(R)
}
end
