*! version 1.0.0 09nov2026
*! Internal exact two-stage support; use the public DEA commands.
version 16.0
capture mata: mata drop _T122_Result()
capture mata: mata drop _t122_ready()
capture mata: mata drop _t122_stop()
capture mata: mata drop _t122_dualdot()
capture mata: mata drop _t122_stage()
mata:
mata set matastrict on
struct _T122_Result {
 real scalar ok,score,radial
 real rowvector diag,weights,lambda,c
 real matrix E,U,T,G
 real colvector b,h,gh
 struct _C58_Cert scalar primary
}
real scalar _t122_ready() {
 string rowvector names
 real scalar j
 names=("_t122_stage()","_t122_dualdot()","_t122_stop()","_t120_eval()","_t120_relerr()","_t120_fit_reason()")
 for(j=1;j<=cols(names);j++) if(findexternal(names[j])==NULL) return(0)
 return(12201)
}
struct _T122_Result scalar _t122_stop(struct _T122_Result scalar R,real scalar reason, |real scalar error) {if(args()>=3 & error!=0 & error<.) {R.diag[18]=error;reason=800+error;} R.diag[21]=reason;return(R);}
// Dot a row with the full exact dual. Share the basis denominator and align
// dyadic row/cost exponents; do not accumulate redundant rational denominators.
struct _C84_Q scalar _t122_dualdot(real rowvector c,struct _C58_Cert scalar fit,struct _C58_Work scalar W) {
 struct _C84_Q scalar q
 transmorphic scalar terms
 real scalar i,e,emin,p
 real rowvector a,N,ex
 p=rows(fit.basis);ex=J(1,p,0);emin=.;terms=asarray_create("real",1)
 for(i=1;i<=p;i++) {
  a=_c70_decode(c[fit.basis[i]],e,W);asarray(terms,i,a)
  ex[i]=e+fit.row_exponent[p+1]-fit.row_exponent[i]
  if(a[1]) emin=min((emin,ex[i]))
 }
 if(emin>=.) emin=0
 N=(0,0)
 for(i=1;i<=p;i++) {
  a=asarray(terms,i)
  if(a[1]) N=_c58b_add(N,_c58b_shift(_c58b_mul(a,_c58b_norm(fit.dual_numerator[i,.]),W),ex[i]-emin,W),W)
 }
 q.N=_c58b_shift(N,max((0,emin)),W);q.D=_c58b_shift(fit.basis_denominator,max((0,-emin)),W)
 return(q)
}
struct _T122_Result scalar _t122_stage(real matrix X,real matrix Y,real scalar k,string scalar rts,string scalar ort,real scalar tol) {
 struct _T122_Result scalar R
 struct _C58_Work scalar W
 struct _C84_Q scalar radial,score,one,zero,q,v,w,cmp,objective,lambdaq
 real scalar n,m,s,p,j,mr,vrs,inp,quality_ok,lambdaerr,scoreerr
 real rowvector ex,ey,intercept,coeff,guard,eg,ec
 real matrix EE,EU
 real colvector Eb,Eh
 R.ok=0;R.score=R.radial=.;R.diag=J(1,27,.);R.diag[1..2]=(0,1)
 n=rows(X);m=cols(X);s=cols(Y);vrs=(rts=="vrs");inp=(ort=="input");p=m+s+2*vrs
 R.weights=J(1,p,.);R.lambda=J(1,n,.)
 if(n<1 | rows(Y)!=n | min((m,s))<1 | k<1 | k>n | k!=floor(k) | k>=. | (rts!="crs" & rts!="vrs") | (ort!="input" & ort!="output")) return(_t122_stop(R,101))
 if(hasmissing((X,Y)) | min((X,Y))<0 | tol<=0 | tol>=.) return(_t122_stop(R,102))
 if(min(rowmax(X))<=0 | min(rowmax(Y))<=0) return(_t122_stop(R,103))
 R.diag[23]=sum(colmax((X,Y)):==0);R.diag[26..27]=(vrs,inp)
 if(n>32) return(_t122_stop(R,201))
 if(p>16) return(_t122_stop(R,202))
 ex=(X[k,.],J(1,s+2*vrs,0));ey=(J(1,m,0),Y[k,.],J(1,2*vrs,0));intercept=J(1,p,0)
 R.U=(-X,Y)
 if(vrs) {intercept[p-1..p]=(1,-1);R.U=(R.U,J(n,1,1),J(n,1,-1));}
 R.E=(inp ? ex : ey);R.b=1;R.h=J(n,1,0);R.c=(inp ? -ey-intercept : ex-intercept)
 R.diag[2]=2;R.T=_c58_map(R.E,R.U,R.c,0,mr);R.diag[22]=mr
 if(!rows(R.T)) return(_t122_stop(R,300+mr))
 R.diag[2]=3
 _c58_fit(R.c*R.T,R.E*R.T,R.b,R.U*R.T,R.h,tol,R.primary,R.G,R.gh,1)
 R.diag[3..4]=(R.primary.status,R.primary.reason);R.diag[7]=R.primary.ops;R.diag[9]=R.primary.maxlimbs;R.diag[24..25]=(R.primary.repair_stop,R.primary.repair_steps)
 if(R.primary.status!=1) return(_t122_stop(R,_t120_fit_reason(R.primary,400)))
 W=_c84_work();one=_c84_number(1,W);zero=_c84_number(0,W);R.diag[2]=5
 objective=_t120_eval(R.c*R.T,R.primary,W);radial=objective;if(inp) radial.N=_c58b_neg(radial.N)
 if(radial.N[1]<=0) return(_t122_stop(R,602,W.error))
 score=(inp ? radial : _c84_div(one,radial,W))
 if(_c84_gt(score,one,W)) return(_t122_stop(R,603,W.error))
 // The technology-row duals are the exact envelopment lambda vector.
 // Verify the original envelopment constraints independently of score range.
 coeff=J(1,rows(R.G),0);coeff[2..n+1]=J(1,n,1);q=_t122_dualdot(coeff,R.primary,W)
 if(vrs) {cmp=_c84_sub(q,one,W);if(cmp.N[1]!=0) return(_t122_stop(R,601,W.error));}
 for(j=1;j<=m;j++) {
  coeff[2..n+1]=X[.,j]';q=_t122_dualdot(coeff,R.primary,W);v=_c84_number(X[k,j],W)
  if(inp) v=_c84_mul(v,radial,W)
  if(_c84_gt(q,v,W)) return(_t122_stop(R,601,W.error))
 }
 for(j=1;j<=s;j++) {
  coeff[2..n+1]=Y[.,j]';q=_t122_dualdot(coeff,R.primary,W);v=_c84_number(Y[k,j],W)
  if(!inp) v=_c84_mul(v,radial,W)
  if(_c84_gt(v,q,W)) return(_t122_stop(R,601,W.error))
 }
 R.diag[2]=6;R.radial=_c84_value(radial);R.score=_c84_value(score);R.weights=_c58_point(R.T,R.primary,1)
 if(R.radial<=0 | R.radial>=. | R.score<=0 | R.score>=.) return(_t122_stop(R,701,W.error))
 if(hasmissing(R.weights) | min(R.weights)<0) return(_t122_stop(R,702,W.error))
 for(j=1;j<=p;j++) if(R.primary.numerator[j,1]!=0 & R.weights[j]==0) return(_t122_stop(R,703,W.error))
 R.lambda=J(1,n,0);quality_ok=1;lambdaerr=0
 for(j=1;j<=n;j++) {
  coeff=J(1,rows(R.G),0);coeff[j+1]=1;lambdaq=_t122_dualdot(coeff,R.primary,W)
  if(lambdaq.N[1]<0) return(_t122_stop(R,601,W.error))
  R.lambda[j]=_c84_value(lambdaq)
  if(R.lambda[j]<0 | R.lambda[j]>=. | (lambdaq.N[1]!=0 & R.lambda[j]==0)) return(_t122_stop(R,703,W.error))
  if(lambdaq.N[1]!=0) lambdaerr=max((lambdaerr,_t120_relerr(_c84_number(R.lambda[j],W),lambdaq,tol,W,quality_ok)))
 }
 if(W.error) {R.diag[18]=W.error;return(_t122_stop(R,800+W.error));}
 R.diag[2]=7
 guard=_g114_guard(R.E,R.b,R.U,R.h,R.c,R.weights,objective,tol)
 R.diag[11..13]=guard[2..4]
 if(!guard[1]) {R.diag[18]=guard[6];return(_t122_stop(R,guard[6] ? 800+guard[6] : 901));}
 q=_c84_dot(R.c,R.weights,0,W);R.diag[10]=_t120_relerr(q,objective,tol,W,quality_ok)
 EE=J(0,n+1,0);Eb=J(0,1,0)
 if(vrs) {EE=(J(1,n,1),0);Eb=1;}
 if(inp) {EU=(X',-X[k,.]')\(-Y',J(s,1,0));Eh=J(m,1,0)\(-Y[k,.]');}
 else {EU=(X',J(m,1,0))\(-Y',Y[k,.]')\(J(1,n,0),-1);Eh=X[k,.]'\J(s,1,0)\(-1);}
 ec=(J(1,n,0),1);eg=_g114_guard(EE,Eb,EU,Eh,ec,(R.lambda,R.radial),radial,tol)
 R.diag[14..15]=eg[2..3];R.diag[16]=lambdaerr
 if(!eg[1]) {R.diag[18]=eg[6];return(_t122_stop(R,eg[6] ? 800+eg[6] : 902));}
 scoreerr=max((_t120_relerr(_c84_number(R.score,W),score,tol,W,quality_ok),_t120_relerr(_c84_number(R.radial,W),radial,tol,W,quality_ok)))
 R.diag[17]=scoreerr;R.diag[18]=W.error;R.diag[20]=W.ops
 if(W.error) return(_t122_stop(R,800+W.error))
 if(!quality_ok | hasmissing((R.diag[10],R.diag[16..17]))) return(_t122_stop(R,903))
 R.ok=1;R.diag[1..2]=(1,0);R.diag[21]=0
 return(R)
}
end
