*! version 1.0.0  01nov2026
*! Shared primal residual check; does not certify optimality or infeasibility.
version 16.0
capture mata: mata drop _dea_lpcheck()
capture mata: mata drop _dea_lpcheck_ok()
mata:
mata set matastrict on

// Returns absolute equality/inequality/bound violations, then scaled ones.
// Row denominator: max(|b|, sum_j |a_j| max(1,|x_j|)). A zero row with
// zero RHS has zero residual. Positive row rescaling leaves the ratio invariant.
// Each finite bound uses max(1,|bound|,|x_j|). Missing bounds are unbounded.
real rowvector _dea_lpcheck(real rowvector x,
    real matrix AE, real colvector BE, real matrix AU, real colvector BU,
    real rowvector lb, real rowvector ub)
{
    real rowvector d, a, w
    real matrix A
    real colvector b
    real scalar n, kind, i, j, scale, rhs, err, den, v, absolute
    n=cols(x); d=J(1,6,0)
    if (n<1 | cols(AE)!=n | cols(AU)!=n | rows(AE)!=rows(BE) |
        rows(AU)!=rows(BU) | cols(lb)!=n | cols(ub)!=n) return(J(1,6,.))
    if (hasmissing(x) | hasmissing(AE) | hasmissing(BE) |
        hasmissing(AU) | hasmissing(BU)) return(J(1,6,.))
    w=rowmax((J(n,1,1),abs(x)'))'
    for (kind=1;kind<=2;kind++) {
        A=(kind==1 ? AE : AU); b=(kind==1 ? BE : BU)
        for (i=1;i<=rows(A);i++) {
            scale=max((abs(A[i,.]),abs(b[i])))
            if (scale==0) continue
            a=A[i,.]/scale; rhs=b[i]/scale
            err=a*x'-rhs
            den=sum(abs(a):*w)
            if (err>=. | den>=.) return(J(1,6,.))
            den=max((abs(rhs),den))
            v=(kind==1 ? abs(err) : max((0,err)))
            absolute=v*scale
            if (absolute>=.) return(J(1,6,.))
            d[kind]=max((d[kind],absolute))
            d[kind+3]=max((d[kind+3],v/den))
        }
    }
    for (j=1;j<=n;j++) {
        if (lb[j]<.) {
            err=lb[j]-x[j]
            if (err>=.) return(J(1,6,.))
            v=max((0,err)); den=max((1,abs(lb[j]),abs(x[j])))
            d[3]=max((d[3],v));d[6]=max((d[6],v/den))
        }
        if (ub[j]<.) {
            err=x[j]-ub[j]
            if (err>=.) return(J(1,6,.))
            v=max((0,err)); den=max((1,abs(ub[j]),abs(x[j])))
            d[3]=max((d[3],v));d[6]=max((d[6],v/den))
        }
    }
    return(d)
}

real scalar _dea_lpcheck_ok(real rowvector d, real scalar tol)
{
    if (cols(d)!=6 | tol<=0 | tol>=. | hasmissing(d)) return(0)
    return(max(d[4..6])<=tol)
}
end
