*! version 1.0.0 08nov2026
*! Bounded failure-only recovery. No disk I/O or command-persistent state.
version 16.0
capture mata: mata drop _g114_ready()
capture mata: mata drop _g114_opt()
capture mata: mata drop _g114_optc()
capture mata: mata drop _g114_guard()
capture mata: mata drop _g114_recover()
capture mata: mata drop _g114_dispatch()
capture mata: mata drop _g114_report()
mata:
mata set matastrict on
real scalar _g114_ready() return(11401)
struct _C84_Q scalar _g114_opt(real rowvector c,struct _C58_Cert scalar f,struct _C58_Work scalar W)
{
 struct _C84_Q scalar q,a,b
 real scalar i
 q=_c84_number(0,W)
 for(i=1;i<=cols(c);i++) {
  a=_c84_number(c[i],W);b.N=_c58b_norm(f.numerator[i,.]);b.D=f.denominator
  q=_c84_sub(q,_c84_mul(a,b,W),W)
 }
 q.N[1]=-q.N[1];return(q)
}
struct _C84_Q scalar _g114_optc(transmorphic scalar C,real scalar ce,struct _C58_Cert scalar f,struct _C58_Work scalar W)
{
 struct _C84_Q scalar q
 real scalar i
 q.N=(0,0);q.D=f.denominator
 for(i=1;i<=rows(f.numerator);i++) q.N=_c58b_add(q.N,_c58b_mul(asarray(C,(1,i)),_c58b_norm(f.numerator[i,.]),W),W)
 q.N=_c58b_shift(_c58b_neg(q.N),max((0,ce)),W);q.D=_c58b_shift(q.D,max((0,-ce)),W)
 return(q)
}
// Exact binary export residuals. Equality and inequality scales use original
// operand magnitudes. Objective compares against the unrounded rational optimum.
real rowvector _g114_guard(real matrix E,real colvector b,real matrix U,real colvector h,real matrix O,real rowvector x,struct _C84_Q scalar opt,real scalar tol)
{
 struct _C58_Work scalar W
 struct _C84_Q scalar zero,one,qt,q,r,den,diff,rhs,scale,total,neg
 real scalar i,eq,iq,gap,ok,bounds
 W=_c84_work();zero=_c84_number(0,W);one=_c84_number(1,W);qt=_c84_number(tol,W)
 if(hasmissing(x)) return((0,.,.,.,0,4,0,0))
 eq=iq=0;ok=1;bounds=all(x:>=0)
 for(i=1;i<=rows(E);i++) {
  q=_c84_dot(E[i,.],x,0,W);rhs=_c84_number(b[i],W);diff=_c84_abs(_c84_sub(q,rhs,W))
  den=_c84_dot(E[i,.],x,1,W);scale=_c84_abs(rhs);if(_c84_gt(scale,den,W)) den=scale
  if(den.N[1]==0) den=one
  r=_c84_div(diff,den,W);eq=max((eq,_c84_value(r)));if(_c84_gt(r,qt,W)) ok=0
 }
 for(i=1;i<=rows(U);i++) {
  q=_c84_dot(U[i,.],x,0,W);rhs=_c84_number(h[i],W);diff=_c84_sub(q,rhs,W)
  if(diff.N[1]<=0) continue
  den=_c84_dot(U[i,.],x,1,W);scale=_c84_abs(rhs);if(_c84_gt(scale,den,W)) den=scale
  if(den.N[1]==0) den=one
  r=_c84_div(diff,den,W);iq=max((iq,_c84_value(r)));if(_c84_gt(r,qt,W)) ok=0
 }
 total=zero
 for(i=1;i<=rows(O);i++) {q=_c84_dot(O[i,.],x,0,W);neg=q;neg.N[1]=-neg.N[1];total=_c84_sub(total,neg,W);}
 scale=_c84_abs(opt);if(!_c84_gt(scale,one,W)) scale=one
 diff=_c84_abs(_c84_sub(total,opt,W));r=_c84_div(diff,scale,W);gap=_c84_value(r)
 if(_c84_gt(r,qt,W)) ok=0
 return((ok & bounds & !W.error,eq,iq,gap,bounds,W.error,W.ops,W.maxlimbs))
}

// diag is a command-local 32-column record; its caller fills native provenance.
real scalar _g114_recover(real matrix X,real matrix Y,real rowvector sx,real rowvector sy,real scalar j,real scalar mode,real scalar own,real scalar d,real scalar theta,real scalar tol,real rowvector v,real rowvector u,real scalar value,real rowvector diag)
{
 real scalar n,m,s,map_reason,ce,ne,primary_value,rc
 real matrix E,U,T,G,G2,F,O,FR
 real colvector b,h,h2,fb
 real rowvector c,c2,od,x,exported,guard
 transmorphic scalar C
 struct _C58_Cert scalar p,f
 struct _C58_Work scalar W
 struct _C84_Q scalar opt
 n=rows(X);m=cols(X);s=cols(Y);primary_value=.;rc=1
 v=J(1,m,.);u=J(1,s,.);value=.;diag[15]=0
 if(n>32 | m+s+(mode==3)>16 | hasmissing(X) | hasmissing(Y) | min(X)<0 | min(Y)<0) {diag[16]=1;return(1);}
 E=(X[j,.],J(1,s,0));b=1;U=(-X,Y);h=J(n,1,0);c=(J(1,m,0),-Y[j,.])
 if(mode==3) {
  E=(E,0)\(-X[d,.],J(1,s,0),1);b=b\0
  U=(U,J(n,1,0))\(J(1,m,0),-Y[d,.],theta);h=h\0;c=(c,0)
 }
 T=_c58_map(E,U,c,0,map_reason)
 if(!rows(T)) {diag[16]=2;diag[18]=map_reason;return(1);}
 _c58_fit(c*T,E*T,b,U*T,h,tol,p,G,h2,1)
 diag[17]=p.status;diag[18]=p.reason;diag[21]=p.ops;diag[23]=p.maxlimbs;diag[24]=p.repair_steps
 if(p.status!=1) {diag[16]=3;return(1);}
 primary_value=-_c58_value(c*T,p,1)
 F=E*T;fb=b;f=p;O=c
 if(mode==1 | mode==2) {
  diag[12]=primary_value;diag[13]=primary_value-own
  if(rows(p.face)) {F=F\G[p.face,.];fb=fb\h2[p.face];}
  _c70_build(J(n,m,0),Y,j,0,mode==1,T,c2,C,ce,od)
  if(od[2]!=1) {diag[16]=4;diag[20]=od[2];return(1);}
  _c70_fit(c2,F,fb,U*T,h,tol,C,ce,od,f,G2,h2,1)
  diag[19]=f.status;diag[20]=f.reason;diag[22]=f.ops;diag[23]=max((p.maxlimbs,f.maxlimbs));diag[25]=f.repair_steps
  if(f.status!=1) {diag[16]=5;return(1);}
  O=(J(n,m,0),Y);O[j,.]=J(1,m+s,0);if(mode==2) O=-O
 }
 x=_c58_point(T,f,1)
 if(hasmissing(x)) {diag[16]=6;return(1);}
 v=x[1..m]:*sx;u=x[(m+1)..(m+s)]:*sy
 exported=(v:/sx,u:/sy)
 if(mode==3) exported=(exported,x[m+s+1])
 W=_c84_work()
 if(mode==3) opt=_g114_opt(c*T,f,W)
 else opt=_g114_optc(C,ce,f,W)
 FR=F
 for(ne=1;ne<=cols(T);ne++) FR[.,ne]=F[.,ne]/T[ne,ne]
 guard=_g114_guard(FR,fb,U,h,O,exported,opt,tol)
 diag[26..30]=guard[2..6]
 diag[16]=6
 if(guard[1] & !W.error & !hasmissing(v) & !hasmissing(u)) {
  value=(mode==3 ? -_c84_value(opt) : _c84_value(opt));diag[16]=7
  if(value<. & (mode!=3 | (value>=0 & value<=1+tol))) {rc=0;diag[16]=0;}
 }
 diag[15]=rc==0;diag[31]=value
 return(rc)
}
real scalar _g114_dispatch(real matrix X,real matrix Y,real matrix XR,real matrix YR,real rowvector sx,real rowvector sy,real scalar j,real scalar mode,real scalar own,real scalar d,real scalar theta,real scalar tol,real rowvector v,real rowvector u,real scalar value,real scalar iteration,real colvector obs,real matrix D,real rowvector counts,real scalar source)
{
 real scalar rc
 real rowvector nd,diag
 source=0
 rc=_dea_cg_lp(X,Y,j,mode,own,d,theta,tol,v,u,value,nd)
 if(!rc) {source=1;return(0);}
 if(mode==0) return(1)
 counts[4+mode]=counts[4+mode]+1
 counts[3]=counts[3]+1
 if(counts[1]>=256) {counts[4]=counts[4]+1;return(1);}
 counts[1]=counts[1]+1
 diag=J(1,32,.)
 diag[1..4]=(iteration,j,obs[j],mode)
 if(mode==3) {diag[5..6]=(d,obs[d]);diag[32]=0;}
 diag[7..10]=nd;diag[11]=own;diag[14]=theta
 rc=_g114_recover(XR,YR,sx,sy,j,mode,own,d,theta,tol,v,u,value,diag)
 D=D\diag
 if(!rc) {
  source=2;counts[2]=counts[2]+1;counts[3]=counts[3]-1
  counts[7+mode]=counts[7+mode]+1
  if(mode!=3 & diag[13]!=0 & diag[13]<.) {counts[11]=counts[11]+1;counts[12]=max((counts[12],abs(diag[13])));}
  counts[13]=max((counts[13],max(diag[26..28])))
 }
 return(rc)
}
void _g114_report(real matrix D,real rowvector C,string scalar prefix,real scalar detail)
{
 real scalar i,n,cap
 string scalar who
 if(C[2]) printf("Note: exact-certified recovery accepted %g native-unavailable LP(s): aggressive %g, benevolent %g, game %g.\n",C[2],C[8],C[9],C[10])
 if(C[11]) printf("Note: %g recovered secondary solve(s) use a certified primary face; original self scores retained. Maximum displayed target change = %12.8g.\n",C[11],C[12])
 if(C[3]) printf("Warning: %g native-unavailable LP(s) remain unresolved (%g beyond the recovery attempt limit); fallback or last complete Nash state retained. This is not a proof of infeasibility.\n",C[3],C[4])
 if(sum(C[5..7])) printf("Note: native optimizer warnings describe that solver attempt; consult recovery provenance and the complete-iteration Nash status for the final result.\n")
 if(detail & rows(D)) {
  printf("\nRecovery: selected evaluator / observation, mode (1 aggressive, 2 benevolent, 3 game), iteration, peer, accepted, reason, native error, target change\n")
  for(i=1;i<=rows(D);i++) printf("%5.0f / %5.0f %3.0f %5.0f %5.0f %3.0f %3.0f %6.0f %14.8g\n",D[i,2],D[i,3],D[i,4],D[i,1],D[i,5],D[i,15],D[i,16],D[i,8],D[i,13])
 }
 cap=st_numscalar("c(matsize)");n=(cap>=32 ? min((rows(D),cap)) : 0)
 if(n) st_matrix(prefix+"_recovery",D[1..n,.])
 if(sum(C[5..7])>n) printf("Note: %g recovery record(s) omitted from the stored matrix because of attempt or Stata matrix-size limits; aggregate counts include them.\n",sum(C[5..7])-n)
 st_matrix(prefix+"_recovery_counts",(C,n,sum(C[5..7])-n))
}
end
