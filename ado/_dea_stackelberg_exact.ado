*! version 1.0.0 09nov2026
*! Internal exact two-stage support; use the public DEA commands.
version 16.0
capture mata: mata drop _t121_ready()
capture mata: mata drop _t121_error()
capture mata: mata drop _t121_nc()
mata:
mata set matastrict on
real scalar _t121_ready() {
 string rowvector names
 real scalar j
 names=("_t121_nc()","_t121_error()","_t120_rel()","_t120_eval()","_t120_relerr()","_t120_stop()","_t120_fit_reason()")
 for(j=1;j<=cols(names);j++) if(findexternal(names[j])==NULL) return(0)
 return(12101)
}
// A true zero is verified by exact equality, not a relative error or epsilon.
real scalar _t121_error(struct _C84_Q scalar a,struct _C84_Q scalar truth,real scalar tol,struct _C58_Work scalar W,real scalar ok) {
 if(truth.N[1]==0) {
  if(a.N[1]!=0 | W.error) {ok=0;return(1);}
  return(0)
 }
 return(_t120_relerr(a,truth,tol,W,ok))
}
struct _T120_Result scalar _t121_nc(real matrix X,real matrix Z,real matrix Y,real scalar k,real scalar lead,real scalar tol) {
 struct _T120_Result scalar R
 struct _C58_Work scalar W
 struct _C84_Q scalar primary,den,mid,num,one,cmp,s1,s2,overall,ed,em,en,r1,r2,ro,target
 real scalar n,m,d,s,p,mr,j,scoreerr,quality_ok
 real rowvector ex,ez,ey,guard,obj
 real matrix FR
 R.ok=0;R.diag=J(1,27,.);R.diag[1]=0;R.diag[2]=1
 R.scores=J(1,3,.);R.weights=J(1,cols(X)+cols(Z)+cols(Y),.)
 n=rows(X);m=cols(X);d=cols(Z);s=cols(Y);p=m+d+s
 if(n<1 | rows(Z)!=n | rows(Y)!=n | min((m,d,s))<1 | k<1 | k>n | k!=floor(k) | k>=. | (lead!=1 & lead!=2)) return(_t120_stop(R,101))
 if(hasmissing((X,Z,Y)) | min((X,Z,Y))<0 | tol<=0 | tol>=.) return(_t120_stop(R,102))
 if(min(rowmax(X))<=0 | min(rowmax(Z))<=0 | min(rowmax(Y))<=0) return(_t120_stop(R,103))
 R.diag[23]=sum(colmax((X,Z,Y)):==0)
 if(n>32) return(_t120_stop(R,201))
 if(p>16) return(_t120_stop(R,202))
 ex=(X[k,.],J(1,d+s,0));ez=(J(1,m,0),Z[k,.],J(1,s,0));ey=(J(1,m+d,0),Y[k,.])
 R.E=ez;R.b=1;R.U=(-X,Z,J(n,s,0))\(J(n,m,0),-Z,Y);R.h=J(2*n,1,0)
 // All weights are shared; the nonleader side can be extended freely in the primary.
 R.c=(lead==1 ? ex : -ey);R.c2=(lead==1 ? -ey : -ez)
 R.diag[2]=2;R.T=_c58_map(R.E,R.U,R.c,0,mr);R.diag[22]=mr
 if(!rows(R.T)) return(_t120_stop(R,300+mr))
 for(j=1;j<=p;j++) if((R.c2[j]*R.T[j,j])/R.T[j,j]!=R.c2[j]) return(_t120_stop(R,304))
 R.diag[2]=3
 _c58_fit(R.c*R.T,R.E*R.T,R.b,R.U*R.T,R.h,tol,R.primary,R.G,R.gh,1)
 R.diag[3..4]=(R.primary.status,R.primary.reason);R.diag[7]=R.primary.ops;R.diag[9]=R.primary.maxlimbs
 R.diag[24..25]=(R.primary.repair_stop,R.primary.repair_steps)
 if(R.primary.status!=1) return(_t120_stop(R,_t120_fit_reason(R.primary,400)))
 // Only homogeneous inequalities/bounds with positive optimal duals define the cone.
 // Stage1 retains wZk=1. Stage2 replaces it by vXk=1, not a rounded score equality.
 R.F=(lead==1 ? ez : ex)*R.T;R.fb=1
 if(rows(R.primary.face)) {
  if(any(R.gh[R.primary.face]:!=0)) return(_t120_stop(R,604))
  R.F=R.F\R.G[R.primary.face,.];R.fb=R.fb\R.gh[R.primary.face]
 }
 R.diag[10]=rows(R.primary.face);R.diag[2]=4
 _c58_fit(R.c2*R.T,R.F,R.fb,R.U*R.T,R.h,tol,R.secondary,R.G2,R.gh2,1)
 R.diag[5..6]=(R.secondary.status,R.secondary.reason);R.diag[8]=R.secondary.ops;R.diag[9]=max((R.diag[9],R.secondary.maxlimbs))
 R.diag[26..27]=(R.secondary.repair_stop,R.secondary.repair_steps)
 if(R.secondary.status!=1) return(_t120_stop(R,_t120_fit_reason(R.secondary,500)))
 W=_c84_work();one=_c84_number(1,W)
 primary=_t120_eval((lead==1 ? ex : ey)*R.T,R.primary,W)
 den=_t120_eval(ex*R.T,R.secondary,W);mid=_t120_eval(ez*R.T,R.secondary,W);num=_t120_eval(ey*R.T,R.secondary,W)
 R.diag[2]=5
 cmp=(lead==1 ? _c84_sub(den,primary,W) : _c84_sub(num,_c84_mul(primary,mid,W),W))
 if(W.error) {R.diag[18]=W.error;return(_t120_stop(R,800+W.error));}
 if(cmp.N[1]!=0) return(_t120_stop(R,601))
 cmp=_c84_sub((lead==1 ? mid : den),one,W)
 if(W.error) {R.diag[18]=W.error;return(_t120_stop(R,800+W.error));}
 if(cmp.N[1]!=0 | den.N[1]<=0 | mid.N[1]<=0 | num.N[1]<0 | primary.N[1]<=0) return(_t120_stop(R,602))
 s1=_c84_div(mid,den,W);s2=_c84_div(num,mid,W);overall=_c84_div(num,den,W)
 quality_ok=(!_c84_gt(s1,one,W) & !_c84_gt(s2,one,W) & !_c84_gt(overall,one,W))
 if(W.error) {R.diag[18]=W.error;return(_t120_stop(R,800+W.error));}
 if(!quality_ok) return(_t120_stop(R,603))
 R.scores=(_c84_value(overall),_c84_value(s1),_c84_value(s2));R.diag[2]=6
 if(hasmissing(R.scores) | min(R.scores)<0) return(_t120_stop(R,701))
 if((overall.N[1]!=0 & R.scores[1]==0) | R.scores[2]==0 | (s2.N[1]!=0 & R.scores[3]==0)) return(_t120_stop(R,701))
 R.weights=_c58_point(R.T,R.secondary,1)
 if(hasmissing(R.weights) | min(R.weights)<0) return(_t120_stop(R,702))
 for(j=1;j<=p;j++) if(R.secondary.numerator[j,1]!=0 & R.weights[j]==0) return(_t120_stop(R,703))
 FR=R.F
 for(j=1;j<=p;j++) FR[.,j]=FR[.,j]/R.T[j,j]
 if(any(FR*R.T:!=R.F)) return(_t120_stop(R,704))
 obj=(lead==1 ? ey : ez);target=(lead==1 ? num : mid)
 guard=_g114_guard(FR,R.fb,R.U,R.h,obj,R.weights,target,tol)
 R.diag[11..13]=guard[2..4];R.diag[19]=guard[7];R.diag[2]=7
 if(!guard[1]) {R.diag[18]=guard[6];return(_t120_stop(R,guard[6] ? 800+guard[6] : 901));}
 ed=_c84_dot(ex,R.weights,0,W);em=_c84_dot(ez,R.weights,0,W);en=_c84_dot(ey,R.weights,0,W)
 if(W.error) {R.diag[18]=W.error;return(_t120_stop(R,800+W.error));}
 if(ed.N[1]<=0 | em.N[1]<=0 | en.N[1]<0) return(_t120_stop(R,902))
 ro=_c84_div(en,ed,W);r1=_c84_div(em,ed,W);r2=_c84_div(en,em,W);quality_ok=1
 R.diag[14..16]=(_t121_error(ro,overall,tol,W,quality_ok),_t121_error(r1,s1,tol,W,quality_ok),_t121_error(r2,s2,tol,W,quality_ok))
 scoreerr=max((_t121_error(_c84_number(R.scores[1],W),overall,tol,W,quality_ok),_t121_error(_c84_number(R.scores[2],W),s1,tol,W,quality_ok),_t121_error(_c84_number(R.scores[3],W),s2,tol,W,quality_ok)))
 R.diag[17]=scoreerr;R.diag[18]=W.error;R.diag[20]=W.ops
 if(W.error) return(_t120_stop(R,800+W.error))
 if(!quality_ok | hasmissing(R.diag[14..17])) return(_t120_stop(R,903))
 R.ok=1;R.diag[1..2]=(1,0);R.diag[21]=0
 return(R)
}
end
