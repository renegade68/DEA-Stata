*! version 1.0.0 09nov2026
*! Bounded original-center fuzzy certificate; alpha=1 only.
version 16.0
capture mata: mata drop _F135_Center()
capture mata: mata drop _f135_ready()
capture mata: mata drop _f135_stop()
capture mata: mata drop _f135_dualdot()
capture mata: mata drop _f135_stage()
capture mata: mata drop _f135_apply()
mata:
mata set matastrict on
struct _F135_Center {
 real scalar ok,score,radial
 real rowvector diag,weights,lambda,c
 real matrix E,U,T,G
 real colvector b,h,gh
 struct _C58_Cert scalar primary
}
real scalar _f135_ready() {
 string rowvector names
 real scalar j
 names=("_f135_apply()","_f135_stage()","_f135_dualdot()","_f135_stop()","_t120_eval()","_t120_relerr()","_t120_fit_reason()")
 for(j=1;j<=cols(names);j++) if(findexternal(names[j])==NULL) return(0)
 return(13501)
}
struct _F135_Center scalar _f135_stop(struct _F135_Center scalar R,real scalar reason, |real scalar error) {if(args()>=3 & error!=0 & error<.) {R.diag[18]=error;reason=800+error;} R.diag[21]=reason;return(R);}
// Dot a row with the full exact dual. Share the basis denominator and align
// dyadic row/cost exponents; do not accumulate redundant rational denominators.
struct _C84_Q scalar _f135_dualdot(real rowvector c,struct _C58_Cert scalar fit,struct _C58_Work scalar W) {
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
struct _F135_Center scalar _f135_stage(real matrix X,real matrix Y,real scalar k,string scalar rts,string scalar ort,real scalar tol) {
 struct _F135_Center scalar R
 struct _C58_Work scalar W
 struct _C84_Q scalar radial,score,one,zero,q,v,w,cmp,objective,lambdaq
 real scalar n,m,s,p,j,mr,vrs,drs,irs,nw,inp,quality_ok,lambdaerr,scoreerr
 real rowvector ex,ey,intercept,coeff,guard,eg,ec
 real matrix EE,EU
 real colvector Eb,Eh
 R.ok=0;R.score=R.radial=.;R.diag=J(1,27,.);R.diag[1..2]=(0,1)
 n=rows(X);m=cols(X);s=cols(Y);vrs=(rts=="vrs");drs=(rts=="drs");irs=(rts=="irs");nw=2*vrs+drs+irs;inp=(ort=="input");p=m+s+nw
 R.weights=J(1,p,.);R.lambda=J(1,n,.)
 if(n<1 | rows(Y)!=n | min((m,s))<1 | k<1 | k>n | k!=floor(k) | k>=. | (rts!="crs" & !vrs & !drs & !irs) | (ort!="input" & ort!="output")) return(_f135_stop(R,101))
 if(hasmissing((X,Y)) | min((X,Y))<0 | tol<=0 | tol>=.) return(_f135_stop(R,102))
 if(min(rowmax(X))<=0 | min(rowmax(Y))<=0) return(_f135_stop(R,103))
 R.diag[23]=sum(colmax((X,Y)):==0);R.diag[26..27]=(vrs+2*drs+3*irs,inp)
 if(n>32) return(_f135_stop(R,201))
 if(p>16) return(_f135_stop(R,202))
 ex=(X[k,.],J(1,s+nw,0));ey=(J(1,m,0),Y[k,.],J(1,nw,0));intercept=J(1,p,0)
 R.U=(-X,Y)
 if(vrs) {intercept[p-1..p]=(1,-1);R.U=(R.U,J(n,1,1),J(n,1,-1));}
 if(drs | irs) {intercept[p]=(irs ? 1 : -1);R.U=(R.U,J(n,1,intercept[p]));}
 R.E=(inp ? ex : ey);R.b=1;R.h=J(n,1,0);R.c=(inp ? -ey-intercept : ex-intercept)
 R.diag[2]=2;R.T=_c58_map(R.E,R.U,R.c,0,mr);R.diag[22]=mr
 if(!rows(R.T)) return(_f135_stop(R,300+mr))
 R.diag[2]=3
 _c58_fit(R.c*R.T,R.E*R.T,R.b,R.U*R.T,R.h,tol,R.primary,R.G,R.gh,1)
 R.diag[3..4]=(R.primary.status,R.primary.reason);R.diag[7]=R.primary.ops;R.diag[9]=R.primary.maxlimbs;R.diag[24..25]=(R.primary.repair_stop,R.primary.repair_steps)
 if(R.primary.status!=1) return(_f135_stop(R,_t120_fit_reason(R.primary,400)))
 W=_c84_work();one=_c84_number(1,W);zero=_c84_number(0,W);R.diag[2]=5
 objective=_t120_eval(R.c*R.T,R.primary,W);radial=objective;if(inp) radial.N=_c58b_neg(radial.N)
 if(radial.N[1]<=0) return(_f135_stop(R,602,W.error))
 score=(inp ? radial : _c84_div(one,radial,W))
 if(_c84_gt(score,one,W)) return(_f135_stop(R,603,W.error))
 // The technology-row duals are the exact envelopment lambda vector.
 // Verify the original envelopment constraints independently of score range.
 coeff=J(1,rows(R.G),0);coeff[2..n+1]=J(1,n,1);q=_f135_dualdot(coeff,R.primary,W)
 if(vrs) {cmp=_c84_sub(q,one,W);if(cmp.N[1]!=0) return(_f135_stop(R,601,W.error));}
 if(drs & _c84_gt(q,one,W) | irs & _c84_gt(one,q,W)) return(_f135_stop(R,601,W.error));
 for(j=1;j<=m;j++) {
  coeff[2..n+1]=X[.,j]';q=_f135_dualdot(coeff,R.primary,W);v=_c84_number(X[k,j],W)
  if(inp) v=_c84_mul(v,radial,W)
  if(_c84_gt(q,v,W)) return(_f135_stop(R,601,W.error))
 }
 for(j=1;j<=s;j++) {
  coeff[2..n+1]=Y[.,j]';q=_f135_dualdot(coeff,R.primary,W);v=_c84_number(Y[k,j],W)
  if(!inp) v=_c84_mul(v,radial,W)
  if(_c84_gt(v,q,W)) return(_f135_stop(R,601,W.error))
 }
 R.diag[2]=6;R.radial=_c84_value(radial);R.score=_c84_value(score);R.weights=_c58_point(R.T,R.primary,1)
 if(R.radial<=0 | R.radial>=. | R.score<=0 | R.score>=.) return(_f135_stop(R,701,W.error))
 if(hasmissing(R.weights) | min(R.weights)<0) return(_f135_stop(R,702,W.error))
 for(j=1;j<=p;j++) if(R.primary.numerator[j,1]!=0 & R.weights[j]==0) return(_f135_stop(R,703,W.error))
 R.lambda=J(1,n,0);quality_ok=1;lambdaerr=0
 for(j=1;j<=n;j++) {
  coeff=J(1,rows(R.G),0);coeff[j+1]=1;lambdaq=_f135_dualdot(coeff,R.primary,W)
  if(lambdaq.N[1]<0) return(_f135_stop(R,601,W.error))
  R.lambda[j]=_c84_value(lambdaq)
  if(R.lambda[j]<0 | R.lambda[j]>=. | (lambdaq.N[1]!=0 & R.lambda[j]==0)) return(_f135_stop(R,703,W.error))
  if(lambdaq.N[1]!=0) lambdaerr=max((lambdaerr,_t120_relerr(_c84_number(R.lambda[j],W),lambdaq,tol,W,quality_ok)))
 }
 if(W.error) {R.diag[18]=W.error;return(_f135_stop(R,800+W.error));}
 R.diag[2]=7
 guard=_g114_guard(R.E,R.b,R.U,R.h,R.c,R.weights,objective,tol)
 R.diag[11..13]=guard[2..4]
 if(!guard[1]) {R.diag[18]=guard[6];return(_f135_stop(R,guard[6] ? 800+guard[6] : 901));}
 q=_c84_dot(R.c,R.weights,0,W);R.diag[10]=_t120_relerr(q,objective,tol,W,quality_ok)
 EE=J(0,n+1,0);Eb=J(0,1,0)
 if(vrs) {EE=(J(1,n,1),0);Eb=1;}
 if(inp) {EU=(X',-X[k,.]')\(-Y',J(s,1,0));Eh=J(m,1,0)\(-Y[k,.]');}
 else {EU=(X',J(m,1,0))\(-Y',Y[k,.]')\(J(1,n,0),-1);Eh=X[k,.]'\J(s,1,0)\(-1);}
 if(drs | irs) {EU=EU\(J(1,n,(drs ? 1 : -1)),0);Eh=Eh\(drs ? 1 : -1);}
 ec=(J(1,n,0),1);eg=_g114_guard(EE,Eb,EU,Eh,ec,(R.lambda,R.radial),radial,tol)
 R.diag[14..15]=eg[2..3];R.diag[16]=lambdaerr
 if(!eg[1]) {R.diag[18]=eg[6];return(_f135_stop(R,eg[6] ? 800+eg[6] : 902));}
 scoreerr=max((_t120_relerr(_c84_number(R.score,W),score,tol,W,quality_ok),_t120_relerr(_c84_number(R.radial,W),radial,tol,W,quality_ok)))
 R.diag[17]=scoreerr;R.diag[18]=W.error;R.diag[20]=W.ops
 if(W.error) return(_f135_stop(R,800+W.error))
 if(!quality_ok | hasmissing((R.diag[10],R.diag[16..17]))) return(_f135_stop(R,903))
 R.ok=1;R.diag[1..2]=(1,0);R.diag[21]=0
 return(R)
}
// N records the unchanged native endpoint fits (six fields per endpoint).
// Certification failure never becomes an infeasibility certificate or zero.
void _f135_apply(string scalar xs,string scalar ys,string scalar use,
 string scalar pess,string scalar opt,string scalar rts,string scalar ort,
 real scalar tol,string scalar native_name,string scalar stats_name,
 string scalar diag_name,string scalar counts_name) {
 real matrix X,Y,N,D
 real colvector P,O,obs
 real rowvector stats,counts
 real scalar n,k,quality,recovered,unresolved
 struct _F135_Center scalar R
 X=st_data(.,tokens(xs),use);Y=st_data(.,tokens(ys),use);obs=selectindex(st_data(.,use))
 N=st_matrix(native_name);n=rows(X);assert(rows(N)==n & cols(N)==12)
 P=N[.,1];O=N[.,7];D=J(n,43,.);counts=J(1,4,0);unresolved=0
 for(k=1;k<=n;k++) {
  R=_f135_stage(X,Y,k,rts,ort,tol)
  quality=(P[k]<. | O[k]<. ? 1 : 0);recovered=0
  if(R.ok) {quality=2;recovered=(P[k]>=. | O[k]>=.);P[k]=O[k]=R.score;}
  else unresolved=unresolved+N[k,5]+N[k,11]
  D[k,.]=(obs[k],quality,R.diag[21],recovered,N[k,.],R.diag)
  counts[quality+1]=counts[quality+1]+1;counts[4]=counts[4]+recovered
 }
 stats=st_matrix(stats_name);stats[5]=unresolved
 st_matrix(stats_name,stats);st_matrix(diag_name,D);st_matrix(counts_name,counts)
 st_store(obs,st_varindex(pess),P);st_store(obs,st_varindex(opt),O)
}

end
