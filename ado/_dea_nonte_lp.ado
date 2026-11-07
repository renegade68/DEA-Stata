*! version 1.0.0  02nov2026
*! Native non-TE LP acceptance against the unmodified model and requested tolerance.
version 16.0
capture mata: assert(_dea_lpcheck_ok(J(1,6,0),1e-8))
if _rc {
    quietly findfile _dea_lpcheck.ado
    quietly run "`r(fn)'"
}
capture mata: mata drop _dea_nonte_lp_ready()
capture mata: mata drop _dea_nonte_optimize()
mata:
mata set matastrict on
real scalar _dea_nonte_lp_ready() return(1)

// 0: finite, converged, primal-accepted; 1: unavailable/rejected.
// This is not an infeasibility or optimality certificate. Snapshot the original
// constraints before optimize(), using the documented LinearProgram getters.
// Retain a model's stricter convergence cap; never loosen it to the user value.
real scalar _dea_nonte_optimize(class LinearProgram scalar S,
    real scalar tol, | real scalar cap)
{
    real matrix E, U, B
    real rowvector c, x
    real scalar n, attempt, use_tol, maxiter, accepted
    string scalar sense
    if (args()<3) cap=1e-8
    if (tol<=0 | tol>=. | cap<=0 | cap>=.) return(1)
    c=S.getCoefficients(); n=cols(c)
    E=S.getEquality(); U=S.getInequality(); B=S.getBounds()
    if (n<1 | hasmissing(c)) return(1)
    if (rows(E)==0) E=J(0,n+1,0)
    if (rows(U)==0) U=J(0,n+1,0)
    if (rows(B)==0) B=J(2,n,.)
    if (cols(E)!=n+1 | cols(U)!=n+1 | rows(B)!=2 | cols(B)!=n) return(1)
    sense=S.getMaxOrMin(); maxiter=S.getMaxiter()
    use_tol=min((tol,cap))
    for (attempt=1;attempt<=2;attempt++) {
        S.setTol(use_tol)
        (void) S.optimize()
        if (!S.converged() | S.errorcode()!=0) return(1)
        x=S.parameters()
        accepted=0
        if (cols(x)==n & !hasmissing(x) & S.value()<.) {
            if (c*x'<.) accepted=_dea_lpcheck_ok(_dea_lpcheck(x,E[,1..n],E[,n+1],
                U[,1..n],U[,n+1],B[1,.],B[2,.]),tol)
        }
        if (accepted) return(0)
        // Only a converged but rejected candidate gets one stricter restart.
        // Rebuild exactly the same LP. Do not relax rows, bounds or user tol.
        if (attempt==2 | use_tol<=min((tol,1e-12))) return(1)
        use_tol=min((tol,1e-12))
        S=LinearProgram()
        S.setMaxOrMin(sense); S.setCoefficients(c)
        if (rows(E)) S.setEquality(E[,1..n],E[,n+1])
        if (rows(U)) S.setInequality(U[,1..n],U[,n+1])
        S.setBounds(B[1,.],B[2,.]); S.setMaxiter(maxiter); S.setTrace("off")
    }
    return(1)
}
end
