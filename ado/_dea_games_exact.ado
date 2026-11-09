*! version 0.1.0 09nov2026
*! T128 isolated score/rank candidate. Derived aggregation is still numerical.
version 16.0
capture mata: mata drop _T128_Scores()
capture mata: mata drop _t128_ready()
capture mata: mata drop _t128_certify()
capture mata: mata drop _t128_rank()
mata:
mata set matastrict on
struct _T128_Scores {
 real matrix D,native,Q,lambda,weights
 real colvector ranks
 real rowvector rankdiag
 transmorphic scalar exact
}
real scalar _t128_ready() {
 string rowvector names
 real scalar j
 names=("_t128_certify()","_t128_rank()","_t122_stage()","_t122_dualdot()","_t122_stop()","_t120_eval()","_t120_relerr()","_t120_fit_reason()","_c84_work()","_c84_gt()","_c84_value()","_dea_games_scores()","_dea_games_rank()","_dea_games_allocate()","_dea_games_coalition()","_dea_games_compute()")
 for(j=1;j<=cols(names);j++) if(findexternal(names[j])==NULL) return(0)
 return(12801)
}
// Q: row,native_site,converged,error,value,score,quality,reason,score,27 helper diagnostics.
// quality:0 unavailable,1 native numerical,2 exact original-data score+export certified.
struct _T128_Scores scalar _t128_certify(real matrix X,real matrix Y,string scalar rts,
 real scalar tol,real matrix D,real matrix N) {
 struct _T128_Scores scalar G
 struct _T122_Result scalar R
 struct _C58_Work scalar W
 struct _C84_Q scalar q
 real scalar i,n,p,reason
 n=rows(X);p=cols(X)+cols(Y)+2*(rts=="vrs")
 G.D=D;G.native=D;G.Q=J(n,36,.);G.exact=asarray_create("real",1)
 G.lambda=J((n<=32 & p<=16 ? n : 0),(n<=32 & p<=16 ? n : 0),.)
 G.weights=J((n<=32 & p<=16 ? n : 0),(n<=32 & p<=16 ? p : 0),.)
 G.ranks=J(n,1,.);G.rankdiag=(0,0,0,0)
 for(i=1;i<=n;i++) {
  G.Q[i,1..9]=(i,N[i,.],D[i,2],D[i,3],.,D[i,2])
  R=_t122_stage(X,Y,i,rts,"input",tol)
  G.Q[i,10..36]=R.diag;reason=R.diag[21]
  if(R.ok) {
   W=_c84_work();q=_t120_eval(R.c*R.T,R.primary,W);q.N=_c58b_neg(q.N)
   if(W.error) reason=1200+W.error
   else if(q.N[1]<=0 | _c84_value(q)!=R.score) reason=1301
   else {
    asarray(G.exact,i,q);G.D[i,.]=(R.score,R.score,1,0)
    G.Q[i,7]=2;G.Q[i,9]=R.score
    G.lambda[i,.]=R.lambda;G.weights[i,.]=R.weights
   }
  }
  G.Q[i,8]=reason
 }
 return(G)
}
// Compare certified rationals, not their rounded doubles. A failed exact
// comparison withholds ranks and preserves the valid individual scores.
void _t128_rank(struct _T128_Scores scalar G) {
 struct _C58_Work scalar W
 struct _C84_Q scalar a,b
 real scalar i,j,n
 n=rows(G.D);G.ranks=_dea_games_rank(G.D[.,2])
 if(hasmissing(G.ranks)) return
 if(any(G.Q[.,7]:!=2)) {G.rankdiag=(1,1001,0,0);return;}
 G.ranks=J(n,1,1);W=_c84_work()
 for(i=1;i<=n;i++) for(j=1;j<=n;j++) if(i!=j) {
  a=asarray(G.exact,i);b=asarray(G.exact,j)
  if(_c84_gt(b,a,W)) G.ranks[i]=G.ranks[i]+1
  if(W.error) {G.ranks=J(n,1,.);G.rankdiag=(3,1200+W.error,W.ops,W.maxlimbs);return;}
 }
 G.rankdiag=(2,0,W.ops,W.maxlimbs)
}
end
