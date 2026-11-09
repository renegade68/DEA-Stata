*! version 1.0.0 09nov2026
*! Internal exact two-stage support; use the public DEA commands.
version 16.0
capture mata: mata drop _T120_Result()
capture mata: mata drop _t120_ready()
capture mata: mata drop _t120_eval()
capture mata: mata drop _t120_relerr()
capture mata: mata drop _t120_rel()
capture mata: mata drop _t120_stop()
capture mata: mata drop _t120_fit_reason()
mata:
mata set matastrict on
struct _T120_Result {
 real scalar ok
 real rowvector scores,weights,diag,c,c2
 real matrix E,U,T,G,G2,F
 real colvector b,h,gh,gh2,fb
 struct _C58_Cert scalar primary,secondary
}
// Detailed private reasons never assert mathematical infeasibility.
struct _T120_Result scalar _t120_stop(struct _T120_Result scalar R,real scalar detail) {
 R.diag[21]=detail
 return(R)
}
real scalar _t120_fit_reason(struct _C58_Cert scalar f,real scalar base) {
 if(f.reason>=11 & f.reason<.) return(base+f.reason)
 if(f.ops>=80000000 & f.ops<. | f.repair_stop==7) return(base+21)
 if(f.repair_steps>=16 & f.repair_steps<. | f.repair_stop==6) return(base+22)
 return(base)
}
real scalar _t120_ready() return(12001)
// Exact dot against the certificate's common rational denominator. Do not
// multiply a fresh denominator at every summand or round the primary target.
struct _C84_Q scalar _t120_eval(real rowvector c,struct _C58_Cert scalar fit,struct _C58_Work scalar W)
{
 struct _C84_Q scalar q
 transmorphic scalar terms
 real rowvector ex,N,a
 real scalar i,e,emin
 terms=asarray_create("real",1);ex=J(1,cols(c),0);emin=.
 for(i=1;i<=cols(c);i++) {
  a=_c70_decode(c[i],e,W);asarray(terms,i,a);ex[i]=e
  if(a[1]) emin=min((emin,e))
 }
 if(emin>=.) emin=0
 N=(0,0)
 for(i=1;i<=cols(c);i++) {
  a=asarray(terms,i)
  if(a[1]) N=_c58b_add(N,_c58b_shift(_c58b_mul(a,_c58b_norm(fit.numerator[i,.]),W),ex[i]-emin,W),W)
 }
 q.N=_c58b_shift(N,max((0,emin)),W);q.D=_c58b_shift(fit.denominator,max((0,-emin)),W)
 return(q)
}
real scalar _t120_relerr(struct _C84_Q scalar a,struct _C84_Q scalar truth,real scalar tol,struct _C58_Work scalar W,real scalar ok)
{
 struct _C84_Q scalar q
 if(truth.N[1]==0) {ok=0;return(.);}
 q=_c84_div(_c84_abs(_c84_sub(a,truth,W)),_c84_abs(truth),W)
 if(_c84_gt(q,_c84_number(tol,W),W) | W.error) ok=0
 return(W.error ? . : _c84_value(q))
}
struct _T120_Result scalar _t120_rel(real matrix X,real matrix Z,real matrix Y,real scalar k,real scalar tol)
{
 struct _T120_Result scalar R
 struct _C58_Work scalar W
 struct _C84_Q scalar primary,den,mid,num,one,zero,cmp,qe,qt,s1,s2,overall,ed,em,en,r1,r2,ro
 real scalar n,m,d,s,p,mr,j,e,scoreerr,quality_ok
 real rowvector ex,ez,ey,guard
 real matrix FR
 R.ok=0;R.diag=J(1,27,.);R.diag[1]=0;R.diag[2]=1
 R.scores=J(1,3,.);R.weights=J(1,cols(X)+cols(Z)+cols(Y),.)
 n=rows(X);m=cols(X);d=cols(Z);s=cols(Y);p=m+d+s
 if(n<1 | rows(Z)!=n | rows(Y)!=n | min((m,d,s))<1 | k<1 | k>n | k!=floor(k) | k>=.) return(_t120_stop(R,101))
 if(hasmissing((X,Z,Y)) | min((X,Z,Y))<0 | tol<=0 | tol>=.) return(_t120_stop(R,102))
 if(min(rowmax(X))<=0 | min(rowmax(Z))<=0 | min(rowmax(Y))<=0) return(_t120_stop(R,103))
 R.diag[23]=sum(colmax((X,Z,Y)):==0)
 if(n>32) return(_t120_stop(R,201))
 if(p>16) return(_t120_stop(R,202))
 ex=(X[k,.],J(1,d+s,0));ez=(J(1,m,0),Z[k,.],J(1,s,0));ey=(J(1,m+d,0),Y[k,.])
 R.E=ex;R.b=1;R.U=(-X,Z,J(n,s,0))\(J(n,m,0),-Z,Y);R.h=J(2*n,1,0);R.c=-ey;R.c2=-ez
 R.diag[2]=2;R.T=_c58_map(R.E,R.U,R.c,0,mr)
 R.diag[22]=mr
 if(!rows(R.T)) return(_t120_stop(R,300+mr))
 // The secondary cost must also round-trip under the same exact dyadic map.
 for(j=1;j<=p;j++) if((R.c2[j]*R.T[j,j])/R.T[j,j]!=R.c2[j]) return(_t120_stop(R,304))
 R.diag[2]=3
 _c58_fit(R.c*R.T,R.E*R.T,R.b,R.U*R.T,R.h,tol,R.primary,R.G,R.gh,1)
 R.diag[3..4]=(R.primary.status,R.primary.reason);R.diag[7]=R.primary.ops;R.diag[9]=R.primary.maxlimbs
 R.diag[24..25]=(R.primary.repair_stop,R.primary.repair_steps)
 if(R.primary.status!=1) return(_t120_stop(R,_t120_fit_reason(R.primary,400)))
 R.F=R.E*R.T;R.fb=R.b
 if(rows(R.primary.face)) {R.F=R.F\R.G[R.primary.face,.];R.fb=R.fb\R.gh[R.primary.face];}
 R.diag[10]=rows(R.primary.face);R.diag[2]=4
 _c58_fit(R.c2*R.T,R.F,R.fb,R.U*R.T,R.h,tol,R.secondary,R.G2,R.gh2,1)
 R.diag[5..6]=(R.secondary.status,R.secondary.reason);R.diag[8]=R.secondary.ops;R.diag[9]=max((R.diag[9],R.secondary.maxlimbs))
 R.diag[26..27]=(R.secondary.repair_stop,R.secondary.repair_steps)
 if(R.secondary.status!=1) return(_t120_stop(R,_t120_fit_reason(R.secondary,500)))
 W=_c84_work();one=_c84_number(1,W);zero=_c84_number(0,W);qt=_c84_number(tol,W)
 primary=_t120_eval(ey*R.T,R.primary,W)
 den=_t120_eval(ex*R.T,R.secondary,W);mid=_t120_eval(ez*R.T,R.secondary,W);num=_t120_eval(ey*R.T,R.secondary,W)
 R.diag[2]=5
 cmp=_c84_sub(num,primary,W)
 if(W.error) {R.diag[18]=W.error;return(_t120_stop(R,800+W.error));}
 if(cmp.N[1]!=0) return(_t120_stop(R,601))
 cmp=_c84_sub(den,one,W)
 if(W.error) {R.diag[18]=W.error;return(_t120_stop(R,800+W.error));}
 if(cmp.N[1]!=0 | mid.N[1]<=0 | num.N[1]<=0) return(_t120_stop(R,602))
 s1=_c84_div(mid,den,W);s2=_c84_div(num,mid,W);overall=_c84_div(num,den,W)
 quality_ok=(!_c84_gt(s1,one,W) & !_c84_gt(s2,one,W) & !_c84_gt(overall,one,W))
 if(W.error) {R.diag[18]=W.error;return(_t120_stop(R,800+W.error));}
 if(!quality_ok) return(_t120_stop(R,603))
 R.scores=(_c84_value(overall),_c84_value(s1),_c84_value(s2));R.diag[2]=6
 if(hasmissing(R.scores) | min(R.scores)<=0) return(_t120_stop(R,701))
 R.weights=_c58_point(R.T,R.secondary,1)
 if(hasmissing(R.weights) | min(R.weights)<0) return(_t120_stop(R,702))
 for(j=1;j<=p;j++) if(R.secondary.numerator[j,1]!=0 & R.weights[j]==0) return(_t120_stop(R,703))
 FR=R.F
 for(j=1;j<=p;j++) FR[.,j]=FR[.,j]/R.T[j,j]
 if(any(FR*R.T:!=R.F)) return(_t120_stop(R,704))
 guard=_g114_guard(FR,R.fb,R.U,R.h,ey,R.weights,primary,tol)
 R.diag[11..13]=guard[2..4];R.diag[19]=guard[7];R.diag[2]=7
 if(!guard[1]) {R.diag[18]=guard[6];return(_t120_stop(R,guard[6] ? 800+guard[6] : 901));}
 // Operand-scaled feasibility plus relative quality of ALL final ratios.
 ed=_c84_dot(ex,R.weights,0,W);em=_c84_dot(ez,R.weights,0,W);en=_c84_dot(ey,R.weights,0,W)
 if(W.error) {R.diag[18]=W.error;return(_t120_stop(R,800+W.error));}
 if(ed.N[1]<=0 | em.N[1]<=0 | en.N[1]<=0) return(_t120_stop(R,902))
 ro=_c84_div(en,ed,W);r1=_c84_div(em,ed,W);r2=_c84_div(en,em,W)
 quality_ok=1
 R.diag[14..16]=(_t120_relerr(ro,overall,tol,W,quality_ok),_t120_relerr(r1,s1,tol,W,quality_ok),_t120_relerr(r2,s2,tol,W,quality_ok))
 scoreerr=max((_t120_relerr(_c84_number(R.scores[1],W),overall,tol,W,quality_ok),_t120_relerr(_c84_number(R.scores[2],W),s1,tol,W,quality_ok),_t120_relerr(_c84_number(R.scores[3],W),s2,tol,W,quality_ok)))
 R.diag[17]=scoreerr;R.diag[18]=W.error;R.diag[20]=W.ops
 if(W.error) return(_t120_stop(R,800+W.error))
 if(!quality_ok | hasmissing(R.diag[14..17])) return(_t120_stop(R,903))
 R.ok=1;R.diag[1..2]=(1,0);R.diag[21]=0
 return(R)
}
end
