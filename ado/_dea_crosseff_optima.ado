*! version 1.1.0 07nov2026
*! Opt-in bounded optimal-face weight and denominator diagnostics; no selection policy.
version 16.0
capture mata: mata drop _c77_optima_ready()
capture mata: mata drop _c73_ready()
capture mata: mata drop _c73_limits()
capture mata: mata drop _c73_allow()
capture mata: mata drop _c73_charge()
capture mata: mata drop _c73_solve()
capture mata: mata drop _c73_second()
capture mata: mata drop _c73_fraction()
capture mata: mata drop _c73_range()
capture mata: mata drop _c73_link()
capture mata: mata drop _c73_face()
capture mata: mata drop _c73_compute()
capture mata: mata drop _C73_Budget()
mata:
mata set matastrict on
// Independent, opt-in original-data diagnostics. No reported weights are changed.
struct _C73_Budget {
 real scalar ops,rowops,calls,rowcalls,reason
 real rowvector limits
}
real scalar _c73_ready() return(7301)
real scalar _c77_optima_ready() return(7701)
real rowvector _c73_limits() return((1000000000,200000000,4096,224))
real scalar _c73_allow(struct _C73_Budget scalar B,real scalar reserve)
{
 if(B.ops+reserve>B.limits[1] | B.rowops+reserve>B.limits[2] |
 B.calls>=B.limits[3] | B.rowcalls>=B.limits[4]) {B.reason=4;return(0);}
 return(1)
}
void _c73_charge(struct _C73_Budget scalar B,real scalar ops)
{
 B.ops=B.ops+ops;B.rowops=B.rowops+ops
}
real scalar _c73_solve(real rowvector c,real matrix E,real colvector b,
 real matrix U,real colvector ub,real scalar tol,struct _C58_Cert scalar fit,
 real matrix G,real colvector h,struct _C73_Budget scalar B)
{
 if(!_c73_allow(B,81000000)) return(-1)
 B.calls=B.calls+1;B.rowcalls=B.rowcalls+1
 _c58_fit(c,E,b,U,ub,tol,fit,G,h)
 _c73_charge(B,fit.ops)
 // T73 linear certificate observation point.
 if(fit.status==1) return(1)
 if(fit.trace62[16]==1) return(2)
 B.reason=(fit.reason>=11 ? 5 : 3)
 return(-1)
}
real scalar _c73_second(real matrix X,real matrix Y,real scalar index,
 real scalar out,real scalar aggressive,real matrix T,real matrix E,real colvector b,
 real matrix U,real colvector ub,real scalar tol,struct _C58_Cert scalar fit,
 real matrix G,real colvector h,struct _C73_Budget scalar B)
{
 real rowvector c,od
 real scalar ce
 transmorphic scalar C
 if(!_c73_allow(B,81000000)) return(-1)
 B.calls=B.calls+1;B.rowcalls=B.rowcalls+1
 _c70_build(X,Y,index,out,aggressive,T,c,C,ce,od)
 if(od[2]!=1) {_c73_charge(B,od[5]);B.reason=5;return(-1);}
 _c70_fit(c,E,b,U,ub,tol,C,ce,od,fit,G,h)
 _c73_charge(B,fit.ops)
 // T73 secondary certificate observation point.
 if(fit.status==1) return(1)
 if(fit.trace62[16]==1) return(2)
 B.reason=(fit.reason>=11 ? 5 : 3)
 return(-1)
}
void _c73_fraction(real rowvector c,struct _C58_Cert scalar fit,
 real rowvector N,real rowvector D,struct _C58_Work scalar W)
{
 transmorphic scalar A
 real rowvector ex
 real scalar re,j
 A=_c58i_rows(c,J(1,1,0),ex,re,W);N=(0,0);D=fit.denominator
 if(W.error) return
 for(j=1;j<=cols(c);j++) N=_c58b_add(N,_c58b_mul(asarray(A,(1,j)),_c58b_norm(fit.numerator[j,.]),W),W)
 N=_c58b_shift(N,max((0,ex[1])),W);D=_c58b_shift(D,max((0,-ex[1])),W)
}
// kind: 1=coordinate,2=peer denominator. States use exact signs, not displayed doubles.
real rowvector _c73_range(real rowvector c,real matrix E,real colvector b,
 real matrix U,real colvector ub,real scalar tol,real scalar cutoff,
 real scalar kind,struct _C73_Budget scalar B)
{
 struct _C58_Cert scalar lo,hi
 struct _C58_Work scalar W
 real matrix G
 real colvector h
 real scalar ls,hs,state,denstate,cutstate,e,low,high,cmp,lcmp,hcmp
 real rowvector N1,D1,N2,D2,v,a,td,z
 ls=_c73_solve(c,E,b,U,ub,tol,lo,G,h,B)
 hs=_c73_solve(-c,E,b,U,ub,tol,hi,G,h,B)
 state=denstate=cutstate=-1;low=high=.;cmp=lcmp=hcmp=.
 if(!_c73_allow(B,2000000)) return((ls,hs,state,denstate,cutstate,low,high))
 W.ops=0;W.limit=1000000;W.limbs=512;W.maxlimbs=0;W.error=0
 if(ls==1) {_c73_fraction(c,lo,N1,D1,W);z=_c58d_double(N1,D1);if(!W.error & !z[2]) low=z[1];}
 if(hs==1) {_c73_fraction(c,hi,N2,D2,W);z=_c58d_double(N2,D2);if(!W.error & !z[2]) high=z[1];}
 if(ls==1 & hs==1) {
  v=_c58b_add(_c58b_mul(N1,D2,W),_c58b_neg(_c58b_mul(N2,D1,W)),W);cmp=v[1]
  if(!W.error & cmp>0) W.error=5
  if(!W.error) state=(cmp==0 ? 1 : 2)
 }
 if(ls==2 | hs==2) state=2
 if(kind==2 & ls==1 & !W.error) {
  if(N1[1]>0) denstate=1
  else if(N1[1]==0 & hs==2) denstate=3
  else if(N1[1]==0 & hs==1) denstate=(N2[1]>0 ? 3 : (N2[1]==0 ? 2 : -1))
  else W.error=5
  a=_c70_decode(cutoff,e,W);td=_c58b_shift((1,1),max((0,-e)),W);a=_c58b_shift(a,max((0,e)),W)
  v=_c58b_add(_c58b_mul(N1,td,W),_c58b_neg(_c58b_mul(a,D1,W)),W);lcmp=v[1]
  if(hs==1) {v=_c58b_add(_c58b_mul(N2,td,W),_c58b_neg(_c58b_mul(a,D2,W)),W);hcmp=v[1];}
  if(lcmp>0) cutstate=1
  else if(hs==1) cutstate=(hcmp<=0 ? 2 : 3)
  else if(hs==2) cutstate=3
 }
 _c73_charge(B,W.ops)
 if(W.error) {state=denstate=cutstate=-1;B.reason=5;}
 // T73 range observation point.
 return((ls,hs,state,denstate,cutstate,low,high))
}
// Check stored binary64 ORIGINAL weights against the exact face. This tests
// weights only, not the separately rounded self score or public CE division.
real scalar _c73_link(real rowvector w,real matrix T,real matrix E,
 real colvector b,real matrix U,real colvector ub,struct _C73_Budget scalar B)
{
 struct _C58_Work scalar W
 transmorphic scalar N,A
 real scalar n,d,i,j,k,e,te,emin,re,value,ok
 real rowvector ex,a,es,D,v,nn
 real matrix G
 real colvector h
 if(hasmissing(w)) return(3)
 if(!_c73_allow(B,2000000)) return(-1)
 W.ops=0;W.limit=1000000;W.limbs=512;W.maxlimbs=0;W.error=0
 n=cols(T);d=rows(T);N=asarray_create("real",1);es=J(1,n,0);emin=0
 for(j=1;j<=n;j++) {
  k=(j<=d ? j : d);value=w[k]
  if(n>d & k==d) value=(j==d ? max((0,value)) : max((0,-value)))
  a=_c70_decode(value,e,W);v=_c70_decode(abs(T[k,j]),te,W)
  if(cols(v)!=2 | v[1]!=1 | v[2]!=1) W.error=5
  es[j]=e-te;emin=min((emin,es[j]));asarray(N,j,a)
  if(W.error) break
 }
 if(W.error) {_c73_charge(B,W.ops);B.reason=5;return(-1);}
 for(j=1;j<=n;j++) asarray(N,j,_c58b_shift(asarray(N,j),es[j]-emin,W))
 D=_c58b_shift((1,1),-emin,W);G=E\U\(-I(n));h=b\ub\J(n,1,0)
 if(W.error) {_c73_charge(B,W.ops);B.reason=5;return(-1);}
 A=_c58i_rows(G,h,ex,re,W);ok=1
 if(W.error) {_c73_charge(B,W.ops);B.reason=5;return(-1);}
 for(i=1;i<=rows(G);i++) {
  v=_c58b_neg(_c58b_mul(asarray(A,(i,n+1)),D,W))
  for(j=1;j<=n;j++) {
   nn=_c58b_shift(asarray(N,j),-re,W)
   v=_c58b_add(v,_c58b_mul(asarray(A,(i,j)),nn,W),W)
  }
  if((i<=rows(E) & v[1]!=0) | (i>rows(E) & v[1]>0)) ok=0
  if(W.error) break
 }
 _c73_charge(B,W.ops)
 if(W.error) {B.reason=5;return(-1);}
 return(ok ? 1 : 2)
}
void _c73_face(real scalar obs,real scalar f,real matrix T,real matrix E,
 real colvector b,real matrix U,real colvector ub,real matrix X,real colvector peers,real scalar tol,
 real scalar cutoff,real rowvector weights,real scalar stage,
 real matrix FD,real matrix CD,real matrix DD,struct _C73_Budget scalar B)
{
 real scalar j,d,n,state,done,ddone,multi,link,startops,startcalls
 real rowvector q,c
 d=rows(T);n=rows(X);done=ddone=multi=0;startops=B.ops;startcalls=B.calls
 for(j=1;j<=d;j++) {
  q=_c73_range(T[j,.],E,b,U,ub,tol,cutoff,1,B)
  CD=CD\(obs,f,j,q[1..3],q[6..7])
  if(q[3]>0) done=done+1
  if(q[3]==2) multi=1
 }
 state=(multi ? 2 : (done==d ? 1 : -1))
 for(j=1;j<=n;j++) {
  c=(X[j,.],J(1,d-cols(X),0))*T
  q=_c73_range(c,E,b,U,ub,tol,cutoff,2,B)
  DD=DD\(obs,f,peers[j],q[1..2],q[4..5],q[6..7])
  if(q[4]>0 & q[5]>0) ddone=ddone+1
 }
 link=_c73_link(weights,T,E,b,U,ub,B)
 FD=FD\(obs,f,1,state,done,ddone,link,stage,(done==d & ddone==n & link>=0 ? 0 : B.reason),B.rowcalls,B.rowops)
}
void _c73_compute(string scalar inputs,string scalar outputs,string scalar touse,
 string scalar rts,string scalar orientation,string scalar method,string scalar vrscross,
 real scalar tol,real scalar supported,string scalar vn,string scalar un,string scalar wn,
 string scalar sn,string scalar fn,string scalar cn,string scalar dn,string scalar sumname,
 | real scalar ce_requested,string scalar cen,string scalar ces)
{
 struct _C73_Budget scalar B
 struct _C58_Cert scalar p,q
 transmorphic scalar CB
 pointer scalar ce_start,ce_begin,ce_face,ce_finish
 real matrix CE
 real scalar ce_supported
 real matrix X,Y,T,E,U,G,F1,F2,FD,CD,DD,V,W
 real colvector b,ub,h,b1,b2,obs,stage
 real rowvector c,weights,S
 real scalar n,m,s,d,free,out,nf,i,k,rc,sc,reason,cutoff,requested
 X=st_data(.,tokens(inputs),touse);Y=st_data(.,tokens(outputs),touse);obs=selectindex(st_data(.,touse))
 n=rows(X);m=cols(X);s=cols(Y);free=rts=="vrs";out=orientation=="output";d=m+s+free
 nf=(method=="arbitrary" ? 1 : 2);requested=n*nf
 FD=J(0,11,.);CD=J(0,8,.);DD=J(0,9,.)
 if(args()<18) ce_requested=0
 ce_supported=supported & vrscross=="legacy" & m+s+2*free<=15
 if(ce_requested) {
  ce_start=findexternal("_c77_start()");ce_begin=findexternal("_c77_begin()")
  ce_face=findexternal("_c77_face()");ce_finish=findexternal("_c77_finish()")
  CB=(*ce_start)();CE=J(0,20,.)
  if(!ce_supported) CE=(.,0,.,-2,-2,-1,-1,.,.,-1,-1,0,0,0,0,0,.,2,0,0)
 }
 B.ops=B.rowops=B.calls=B.rowcalls=B.reason=0;B.limits=_c73_limits()
 if(!supported | vrscross!="legacy") {
  FD=(.,0,-2,-1,0,0,0,.,1,0,0)
  st_matrix(fn,FD);st_matrix(cn,J(1,8,.));st_matrix(dn,J(1,9,.))
  if(ce_requested) (*ce_finish)(CE,n*requested,CB,cen,ces)
  st_matrix(sumname,(requested,0,0,0,requested,requested,0,0,0,0,0,0,0,0,0,0,0,0));return
 }
 V=(st_matrix(vn),st_matrix(un));W=st_matrix(wn);stage=st_matrix(sn)
 cutoff=(free ? tol : 0)
 for(i=1;i<=n;i++) {
  if(ce_requested & ce_supported) (*ce_begin)(CB,n)
  B.rowops=B.rowcalls=B.reason=0
  c=(out ? (X[i,.],J(1,s,0)) : (J(1,m,0),-Y[i,.]))
  E=(out ? (J(1,m,0),Y[i,.]) : (X[i,.],J(1,s,0)));U=(-X,Y)
  if(free) {c=(c,(out ? 1 : -1));E=(E,0);U=(U,J(n,1,(out ? -1 : 1)));}
  b=1;ub=J(n,1,0);weights=V[i,.];if(free) weights=(weights,W[i])
  T=_c58_map(E,U,c,free,reason)
  if(rows(T)==0) {
   for(k=1;k<=nf;k++) {
    FD=FD\(obs[i],k,-1,-1,0,0,0,stage[i],2,0,0)
    if(ce_requested & ce_supported) (*ce_face)(obs[i],k,-1,J(0,0,.),J(0,0,.),J(0,1,.),J(0,0,.),J(0,1,.),X,Y,obs,tol,cutoff,0,stage[i],CE,CB)
   }
   continue
  }
  rc=_c73_solve(c*T,E*T,b,U*T,ub,tol,p,G,h,B)
  if(rc!=1) {
   FD=FD\(obs[i],1,rc,(rc==2 ? 0 : -1),0,0,0,stage[i],B.reason,B.rowcalls,B.rowops)
   if(ce_requested & ce_supported) (*ce_face)(obs[i],1,rc,J(0,0,.),J(0,0,.),J(0,1,.),J(0,0,.),J(0,1,.),X,Y,obs,tol,cutoff,0,stage[i],CE,CB)
   if(nf==2) {
    FD=FD\(obs[i],2,-1,-1,0,0,0,stage[i],6,0,0)
    if(ce_requested & ce_supported) (*ce_face)(obs[i],2,-1,J(0,0,.),J(0,0,.),J(0,1,.),J(0,0,.),J(0,1,.),X,Y,obs,tol,cutoff,0,stage[i],CE,CB)
   }
   continue
  }
  F1=E*T;b1=b
  if(rows(p.face)) {F1=F1\G[p.face,.];b1=b1\h[p.face];}
  // T73 F1 face observation point.
  _c73_face(obs[i],1,T,F1,b1,U*T,ub,X,obs,tol,cutoff,weights,stage[i],FD,CD,DD,B)
  if(ce_requested & ce_supported) (*ce_face)(obs[i],1,1,T,F1,b1,U*T,ub,X,Y,obs,tol,cutoff,FD[rows(FD),7],stage[i],CE,CB)
  if(nf==2) {
   sc=_c73_second(X,Y,i,out,method=="aggressive",T,F1,b1,U*T,ub,tol,q,G,h,B)
   if(sc!=1) {FD=FD\(obs[i],2,sc,(sc==2 ? 0 : -1),0,0,0,stage[i],(sc==2 ? 0 : B.reason),B.rowcalls,B.rowops);
    if(ce_requested & ce_supported) (*ce_face)(obs[i],2,sc,J(0,0,.),J(0,0,.),J(0,1,.),J(0,0,.),J(0,1,.),X,Y,obs,tol,cutoff,0,stage[i],CE,CB)
    continue;
   }
   F2=F1;b2=b1
   if(rows(q.face)) {F2=F2\G[q.face,.];b2=b2\h[q.face];}
   // T73 F2 face observation point.
   _c73_face(obs[i],2,T,F2,b2,U*T,ub,X,obs,tol,cutoff,weights,stage[i],FD,CD,DD,B)
  if(ce_requested & ce_supported) (*ce_face)(obs[i],2,1,T,F2,b2,U*T,ub,X,Y,obs,tol,cutoff,FD[rows(FD),7],stage[i],CE,CB)
  }
 }
 S=(requested,sum(FD[.,3]:==1),sum(FD[.,4]:==1),sum(FD[.,4]:==2),sum(FD[.,4]:<0),0,
 sum(FD[.,3]:==2),sum(FD[.,7]:==1),sum(FD[.,7]:==2),sum(FD[.,7]:==3),
 sum(DD[.,6]:==1),sum(DD[.,6]:==2),sum(DD[.,6]:==3),sum(DD[.,6]:<0),sum(DD[.,7]:==3),B.calls,B.ops,sum(FD[.,7]:==-1))
 if(rows(CD)==0) CD=J(1,8,.)
 if(rows(DD)==0) DD=J(1,9,.)
 if(ce_requested) (*ce_finish)(CE,n*requested,CB,cen,ces)
 st_matrix(fn,FD);st_matrix(cn,CD);st_matrix(dn,DD);st_matrix(sumname,S)
}

end
