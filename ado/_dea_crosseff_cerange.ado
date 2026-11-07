*! version 1.0.0 07nov2026
*! Bounded exact CE-range core and public diagnostic adapter; no weight selection.
version 16.0
capture mata: mata drop _c76_ready()
capture mata: mata drop _c76_budget()
capture mata: mata drop _c76_begin()
capture mata: mata drop _c76_allow()
capture mata: mata drop _c76_charge()
capture mata: mata drop _c76_work()
capture mata: mata drop _c76_empty()
capture mata: mata drop _c76_endempty()
capture mata: mata drop _c76_result()
capture mata: mata drop _c76_certpoint()
capture mata: mata drop _c76_pack()
capture mata: mata drop _c76_dot()
capture mata: mata drop _c76_cmp()
capture mata: mata drop _c76_check()
capture mata: mata drop _c76_positive()
capture mata: mata drop _c76_homogenize()
capture mata: mata drop _c76_dehomogenize()
capture mata: mata drop _c76_display()
capture mata: mata drop _c76_solve()
capture mata: mata drop _c76_endpoint()
capture mata: mata drop _c76_cell()
capture mata: mata drop _c76_coeff()
capture mata: mata drop _c76_original()
capture mata: mata drop _c76_peer()
capture mata: mata drop _C76_Budget()
capture mata: mata drop _C76_Result()
capture mata: mata drop _C76_End()
capture mata: mata drop _C76_Point()
capture mata: mata drop _c77_ready()
capture mata: mata drop _c77_start()
capture mata: mata drop _c77_begin()
capture mata: mata drop _c77_face()
capture mata: mata drop _c77_finish()
mata:
mata set matastrict on
// T76 internal research helper. Minimize costs; caller preserves solver context.
struct _C76_Point {
 real matrix N
 real rowvector D
}
struct _C76_End {
 real scalar state,value,display
 real rowvector N,D
 struct _C76_Point scalar point,ray
}
struct _C76_Result {
 real scalar state,denstate,cutstate,reason,rawlower_state,rawupper_state
 struct _C76_Point scalar anchor
 struct _C76_End scalar lower,upper
}
struct _C76_Budget {
 real scalar ops,rowops,calls,rowcalls,reason,maxlimbs,evaluators,peers
 real scalar mathlimit,mathlimbs
 real rowvector limits
}
real scalar _c76_ready() return(7601)
struct _C76_Budget scalar _c76_budget()
{
 struct _C76_Budget scalar B
 B.ops=B.rowops=B.calls=B.rowcalls=B.reason=B.maxlimbs=0
 B.evaluators=0;B.peers=1
 B.limits=(1000000000,200000000,4096,512)
 B.mathlimit=1000000;B.mathlimbs=512
 return(B)
}
real scalar _c76_begin(struct _C76_Budget scalar B,real scalar peers)
{
 if(peers<1 | peers>32 | peers!=floor(peers) | B.evaluators>=32) {B.reason=2;return(0);}
 B.evaluators=B.evaluators+1;B.peers=peers;B.rowops=B.rowcalls=0
 return(1)
}
real scalar _c76_allow(struct _C76_Budget scalar B,real scalar reserve,real scalar lp)
{
 if(B.ops+reserve>B.limits[1] | B.rowops+reserve>B.limits[2] |
 (lp & (B.calls>=B.limits[3] | B.rowcalls>=B.limits[4]))) {B.reason=4;return(0);}
 return(1)
}
void _c76_charge(struct _C76_Budget scalar B,struct _C58_Work scalar W)
{
 B.ops=B.ops+W.ops;B.rowops=B.rowops+W.ops;B.maxlimbs=max((B.maxlimbs,W.maxlimbs))
 if(W.error) B.reason=5
}
struct _C58_Work scalar _c76_work(struct _C76_Budget scalar B)
{
 struct _C58_Work scalar W
 W.ops=W.maxlimbs=W.error=0;W.limit=B.mathlimit;W.limbs=B.mathlimbs
 return(W)
}
struct _C76_Point scalar _c76_empty()
{
 struct _C76_Point scalar P
 P.N=J(0,0,.);P.D=(1,1);return(P)
}
struct _C76_End scalar _c76_endempty()
{
 struct _C76_End scalar R
 R.state=-1;R.value=.;R.display=0;R.N=(0,0);R.D=(1,1)
 R.point=R.ray=_c76_empty();return(R)
}
struct _C76_Result scalar _c76_result()
{
 struct _C76_Result scalar R
 R.state=R.denstate=R.cutstate=-1;R.reason=0
 R.rawlower_state=R.rawupper_state=0
 R.anchor=_c76_empty();R.lower=R.upper=_c76_endempty();return(R)
}
struct _C76_Point scalar _c76_certpoint(struct _C58_Cert scalar q,real scalar raybase)
{
 struct _C76_Point scalar P
 if(raybase) {P.N=q.ray_num;P.D=q.ray_den;}
 else {P.N=q.numerator;P.D=q.denominator;}
 return(P)
}
struct _C76_Point scalar _c76_pack(transmorphic scalar A,real scalar n,
 real rowvector D,struct _C58_Work scalar W)
{
 struct _C76_Point scalar P
 real rowvector a
 real scalar j,k
 P.D=D;k=2
 for(j=1;j<=n;j++) {a=asarray(A,j);k=max((k,cols(a)));}
 if(!_c58b_tick(W,n*k,k-1)) return(_c76_empty())
 P.N=J(n,k,0)
 for(j=1;j<=n;j++) {a=asarray(A,j);P.N[j,1..cols(a)]=a;}
 return(P)
}
void _c76_dot(real rowvector c,struct _C76_Point scalar P,
 real rowvector N,real rowvector D,struct _C58_Work scalar W)
{
 transmorphic scalar A
 real rowvector es,a
 real scalar e,emin,j
 N=(0,0);D=P.D;emin=0;es=J(1,cols(c),0);A=asarray_create("real",1)
 for(j=1;j<=cols(c);j++) {a=_c70_decode(c[j],e,W);es[j]=e;emin=min((emin,e));asarray(A,j,a);}
 for(j=1;j<=cols(c);j++) {
  a=_c58b_shift(asarray(A,j),es[j]-emin,W)
  N=_c58b_add(N,_c58b_mul(a,_c58b_norm(P.N[j,.]),W),W)
 }
 D=_c58b_shift(D,-emin,W)
}
real scalar _c76_cmp(real rowvector N,real rowvector D,real rowvector M,
 real rowvector H,struct _C58_Work scalar W)
{
 real rowvector a
 a=_c58b_add(_c58b_mul(N,H,W),_c58b_neg(_c58b_mul(M,D,W)),W)
 return(a[1])
}
// Exact original rows and sign bounds; malformed evidence is an internal error.
void _c76_check(struct _C76_Point scalar P,real matrix E,real colvector b,
 real matrix U,real colvector ub,struct _C58_Work scalar W)
{
 real rowvector N,D,H,K,a
 real scalar j,e,sgn
 if(W.error) return
 if(rows(P.N)!=cols(E) | P.D[1]!=1) _error(3498)
 for(j=1;j<=rows(P.N);j++) if(P.N[j,1]<0) _error(3498)
 for(j=1;j<=rows(E)+rows(U);j++) {
  _c76_dot((j<=rows(E) ? E[j,.] : U[j-rows(E),.]),P,N,D,W)
  H=_c70_decode((j<=rows(E) ? b[j] : ub[j-rows(E)]),e,W)
  H=_c58b_shift(H,max((0,e)),W);K=_c58b_shift((1,1),max((0,-e)),W)
  sgn=_c76_cmp(N,D,H,K,W)
  if(!W.error & ((j<=rows(E) & sgn!=0) | (j>rows(E) & sgn>0))) _error(3498)
  if(W.error) return
 }
}
struct _C76_Point scalar _c76_positive(struct _C58_Cert scalar q,real scalar status,
 real rowvector c,struct _C58_Work scalar W)
{
 struct _C76_Point scalar P
 transmorphic scalar A
 real rowvector N,D,a
 real scalar j
 P=_c76_certpoint(q,status==2);_c76_dot(c,P,N,D,W)
 if(W.error | N[1]>0) return(P)
 if(status!=2) _error(3498)
 A=asarray_create("real",1)
 for(j=1;j<=rows(P.N);j++) {
  a=_c58b_add(_c58b_norm(P.N[j,.]),_c58b_mul(_c58b_norm(q.ray[j,.]),P.D,W),W)
  asarray(A,j,a)
 }
 P=_c76_pack(A,rows(P.N),P.D,W)
 if(W.error) return(P)
 _c76_dot(c,P,N,D,W);if(!W.error & N[1]<=0) _error(3498)
 return(P)
}
struct _C76_Point scalar _c76_homogenize(struct _C76_Point scalar P,
 real rowvector den,struct _C58_Work scalar W)
{
 transmorphic scalar A
 real rowvector N,D,H
 real scalar j,n
 _c76_dot(den,P,N,D,W);if(W.error) return(_c76_empty())
 if(N[1]<=0) _error(3498)
 n=rows(P.N);A=asarray_create("real",1)
 for(j=1;j<=n;j++) asarray(A,j,_c58b_mul(_c58b_norm(P.N[j,.]),D,W))
 asarray(A,n+1,_c58b_mul(P.D,D,W));H=_c58b_mul(P.D,N,W)
 return(_c76_pack(A,n+1,H,W))
}
struct _C76_Point scalar _c76_dehomogenize(struct _C76_Point scalar P)
{
 struct _C76_Point scalar Q
 real scalar n
 n=rows(P.N);Q.N=P.N[1..n-1,.];Q.D=_c58b_norm(P.N[n,.])
 if(Q.D[1]<=0) _error(3498)
 return(Q)
}
// Conversion work is conservatively reserved/charged at its full 1m ceiling.
void _c76_display(struct _C76_End scalar R,struct _C76_Budget scalar B)
{
 real rowvector z
 if(!_c76_allow(B,1000000,0)) {R.display=5;return;}
 z=_c58d_double(R.N,R.D)
 B.ops=B.ops+1000000;B.rowops=B.rowops+1000000
 R.value=z[1]
 if(z[2]==3) R.display=4
 else if(z[2]) {R.display=5;B.reason=5;}
 else R.display=(z[1]==0 & R.N[1]!=0 ? 3 : (z[3] ? 2 : 1))
}
real scalar _c76_solve(real rowvector c,real matrix E,real colvector b,
 real matrix U,real colvector ub,real scalar tol,struct _C58_Cert scalar q,
 real matrix G,real colvector h,struct _C76_Budget scalar B)
{
 if(cols(c)>16 | rows(E)+rows(U)>64 | rows(E)+rows(U)+cols(c)>256) {B.reason=2;return(-2);}
 if(!_c76_allow(B,81000000,1)) return(-1)
 B.calls=B.calls+1;B.rowcalls=B.rowcalls+1
 _c58_fit(c,E,b,U,ub,tol,q,G,h)
 B.ops=B.ops+q.ops;B.rowops=B.rowops+q.ops;B.maxlimbs=max((B.maxlimbs,q.maxlimbs))
 // T76 solve observation point.
 if(q.status==1) return(1)
 if(q.trace62[16]==1) return(2)
 B.reason=(q.reason>=11 ? 5 : 3);return(-1)
}
void _c76_endpoint(real rowvector num,real rowvector den,real matrix E,
 real colvector b,real matrix U,real colvector ub,real matrix CE,real colvector cb,
 real matrix CU,real colvector cub,real scalar sense,real scalar tol,
 struct _C76_Point scalar anchor,struct _C76_End scalar R,struct _C76_Budget scalar B)
{
 struct _C58_Cert scalar q,tq
 struct _C58_Work scalar W
 struct _C76_Point scalar P,Q
 real matrix G,OE
 real colvector h,ob
 real rowvector c,t,N,D,M,H
 real scalar status,ts,n,sgn
 n=cols(CE);c=(num,0);R=_c76_endempty()
 status=_c76_solve(sense*c,CE,cb,CU,cub,tol,q,G,h,B)
 if(status<0) return
 if(!_c76_allow(B,2000000,0)) return
 W=_c76_work(B)
 if(status==2) {
  R.ray.N=q.ray;R.ray.D=(1,1);R.point=anchor
  _c76_check(R.ray,CE,J(rows(CE),1,0),CU,J(rows(CU),1,0),W)
  _c76_dot(c,R.ray,N,D,W)
  if(!W.error & sense*N[1]>=0) _error(3498)
  if(!W.error) R.state=3
  _c76_charge(B,W);return
 }
 P=_c76_certpoint(q,0);_c76_dot(c,P,R.N,R.D,W)
 _c76_charge(B,W);if(W.error) return
 // Dual-active rows retain the entire exact endpoint face, with dyadic RHS.
 if(P.N[n,1]==0) {
  OE=CE;ob=cb
  if(rows(q.face)) {OE=OE\G[q.face,.];ob=ob\h[q.face];}
  t=J(1,n,0);t[n]=1
  ts=_c76_solve(-t,OE,ob,CU,cub,tol,tq,G,h,B)
  if(ts<0) return
  if(!_c76_allow(B,2000000,0)) return
  W=_c76_work(B)
  if(ts==2) P=_c76_positive(tq,ts,t,W)
  else if(tq.numerator[n,1]>0) P=_c76_certpoint(tq,0)
  _c76_charge(B,W);if(W.error) return
 }
 if(!_c76_allow(B,2000000,0)) return
 W=_c76_work(B)
 _c76_check(P,CE,cb,CU,cub,W);_c76_dot(c,P,N,D,W)
 sgn=_c76_cmp(N,D,R.N,R.D,W);if(!W.error & sgn!=0) _error(3498)
 if(!W.error & P.N[n,1]>0) {
  Q=_c76_dehomogenize(P);_c76_check(Q,E,b,U,ub,W)
  _c76_dot(den,Q,N,D,W);if(!W.error & N[1]<=0) _error(3498)
  R.point=Q;if(!W.error) R.state=1
 }
 else if(!W.error) {R.point=P;R.state=2;}
 _c76_charge(B,W)
 if(W.error) {R.state=-1;return;}
 _c76_display(R,B)
}
// face_state:1 certified finite face,2 exact unbounded secondary,other unavailable.
void _c76_cell(real rowvector num,real rowvector den,real matrix E,
 real colvector b,real matrix U,real colvector ub,real scalar cutoff,real scalar tol,
 real scalar face_state,struct _C76_Budget scalar B,struct _C76_Result scalar R)
{
 struct _C58_Cert scalar lo,hi
 struct _C58_Work scalar W
 struct _C76_Point scalar P
 real matrix G,CE,CU
 real colvector h,cb,cub
 real rowvector N,D,M,H,K,L
 real scalar n,ls,hs,e,lcmp,hcmp,cmp
 R=_c76_result();B.reason=0;n=cols(num)
 if(face_state!=1) {R.state=-3;R.reason=8;return;}
 if(n<1 | n+1>16 | rows(E)+rows(U)+1>64 | cols(den)!=n | cols(E)!=n |
 cols(U)!=n | rows(b)!=rows(E) | rows(ub)!=rows(U) | cols(b)!=1 | cols(ub)!=1 |
 cutoff<0 | missing(cutoff) | missing(tol) | tol<=0 | hasmissing((num,den)) |
 hasmissing(E) | hasmissing(U) | hasmissing(b) | hasmissing(ub)) {R.state=-2;R.reason=2;return;}
 CE=(E,-b)\(den,0);cb=J(rows(E),1,0)\1;CU=(U,-ub);cub=J(rows(U),1,0)
 ls=_c76_solve(den,E,b,U,ub,tol,lo,G,h,B)
 hs=_c76_solve(-den,E,b,U,ub,tol,hi,G,h,B)
 if(ls<0 | hs<0) {R.reason=B.reason;return;}
 if(ls==2) {R.state=-2;R.reason=7;return;}
 if(!_c76_allow(B,2000000,0)) {R.reason=B.reason;return;}
 W=_c76_work(B);P=_c76_certpoint(lo,0);_c76_dot(den,P,N,D,W)
 if(!W.error & N[1]<0) {R.state=-2;R.reason=7;_c76_charge(B,W);return;}
 if(hs==1) {P=_c76_certpoint(hi,0);_c76_dot(den,P,M,H,W);}
 else {M=(1,1);H=(1,1);}
 if(!W.error) R.denstate=(N[1]>0 ? 1 : (M[1]==0 ? 2 : 3))
 K=_c70_decode(cutoff,e,W);K=_c58b_shift(K,max((0,e)),W);L=_c58b_shift((1,1),max((0,-e)),W)
 lcmp=_c76_cmp(N,D,K,L,W);hcmp=(hs==2 ? 1 : _c76_cmp(M,H,K,L,W))
 if(!W.error) R.cutstate=(lcmp>0 ? 1 : (hcmp<=0 ? 2 : 3))
 _c76_charge(B,W)
 if(W.error) {R.denstate=R.cutstate=-1;R.reason=B.reason;return;}
 if(R.denstate==2) {R.state=0;return;}
 if(!_c76_allow(B,2000000,0)) {R.reason=B.reason;return;}
 W=_c76_work(B);P=_c76_positive(hi,hs,den,W)
 _c76_check(P,E,b,U,ub,W)
 if(!W.error) R.anchor=_c76_homogenize(P,den,W)
 if(!W.error) _c76_check(R.anchor,CE,cb,CU,cub,W)
 _c76_charge(B,W);if(W.error) {R.reason=B.reason;return;}
 _c76_endpoint(num,den,E,b,U,ub,CE,cb,CU,cub,1,tol,R.anchor,R.lower,B)
 _c76_endpoint(num,den,E,b,U,ub,CE,cb,CU,cub,-1,tol,R.anchor,R.upper,B)
 if(R.lower.state<0 | R.upper.state<0) {R.reason=B.reason;return;}
 if(R.lower.state==3 | R.upper.state==3) {R.state=2;return;}
 if(!_c76_allow(B,2000000,0)) {R.reason=B.reason;return;}
 W=_c76_work(B);cmp=_c76_cmp(R.lower.N,R.lower.D,R.upper.N,R.upper.D,W)
 if(!W.error & cmp>0) _error(3498)
 if(!W.error) R.state=(cmp==0 ? 1 : 2)
 if(R.state==1 & (R.lower.state!=1 | R.upper.state!=1)) _error(3498)
 _c76_charge(B,W);if(W.error) R.state=-1
 R.reason=B.reason
}
// T73 maps have one signed power-of-two entry per split/scaled column.
real scalar _c76_coeff(real rowvector raw,real matrix T,real rowvector c,
 struct _C76_Budget scalar B)
{
 struct _C58_Work scalar W
 real rowvector a
 real scalar j,k,i,e
 c=J(1,cols(T),.)
 if(rows(T)!=cols(raw) | cols(T)>15 | rows(T)>16 | hasmissing(raw) | hasmissing(T)) {B.reason=6;return(0);}
 if(!_c76_allow(B,2000000,0)) return(0)
 W=_c76_work(B)
 for(j=1;j<=cols(T);j++) {
  k=0
  for(i=1;i<=rows(T);i++) if(T[i,j]!=0) {if(k) {W.error=4;break;}k=i;}
  if(!k | W.error) {W.error=4;break;}
  a=_c70_decode(T[k,j],e,W)
  if(W.error) break
  if(cols(a)!=2 | abs(a[2])!=1) {W.error=4;break;}
  c[j]=raw[k]*T[k,j]
  if(missing(c[j]) | c[j]/T[k,j]!=raw[k] | (raw[k]!=0 & c[j]==0)) {W.error=4;break;}
 }
 _c76_charge(B,W);if(W.error) {B.reason=(W.error==4 ? 6 : 5);return(0);}
 return(1)
}
struct _C76_Point scalar _c76_original(real matrix T,struct _C76_Point scalar P,
 struct _C58_Work scalar W)
{
 transmorphic scalar A,V
 real matrix ex
 real rowvector a,N,D
 real scalar i,j,e,emin
 A=asarray_create("real",2);V=asarray_create("real",1);emin=0;ex=J(rows(T),cols(T),0)
 for(i=1;i<=rows(T);i++) for(j=1;j<=cols(T);j++) {
  a=_c70_decode(T[i,j],e,W);asarray(A,(i,j),a);ex[i,j]=e;emin=min((emin,e))
 }
 for(i=1;i<=rows(T);i++) {
  N=(0,0)
  for(j=1;j<=cols(T);j++) {
   a=_c58b_shift(asarray(A,(i,j)),ex[i,j]-emin,W)
   N=_c58b_add(N,_c58b_mul(a,_c58b_norm(P.N[j,.]),W),W)
  }
  asarray(V,i,N)
 }
 D=_c58b_shift(P.D,-emin,W)
 return(_c76_pack(V,rows(T),D,W))
}
void _c76_peer(real rowvector numerator,real rowvector denominator,real matrix T,
 real matrix E,real colvector b,real matrix U,real colvector ub,real scalar cutoff,
 real scalar tol,real scalar face_state,struct _C76_Budget scalar B,
 struct _C76_Result scalar R,struct _C76_Point scalar rawlower,struct _C76_Point scalar rawupper)
{
 struct _C58_Work scalar W
 real rowvector num,den
 R=_c76_result();rawlower=rawupper=_c76_empty();B.reason=0
 if(face_state!=1) {R.state=-3;R.reason=8;return;}
 if(!_c76_coeff(numerator,T,num,B) | !_c76_coeff(denominator,T,den,B)) {
  R.state=(B.reason==6 ? -2 : -1);R.reason=B.reason;return
 }
 _c76_cell(num,den,E,b,U,ub,cutoff,tol,face_state,B,R)
 R.rawlower_state=(R.lower.state==1 ? -1 : 0);R.rawupper_state=(R.upper.state==1 ? -1 : 0)
 if(R.lower.state!=1 & R.upper.state!=1) return
 if(!_c76_allow(B,2000000,0)) {R.reason=B.reason;return;}
 W=_c76_work(B)
 if(R.lower.state==1) rawlower=_c76_original(T,R.lower.point,W)
 if(R.upper.state==1 & !W.error) rawupper=_c76_original(T,R.upper.point,W)
 _c76_charge(B,W)
 if(W.error) {rawlower=rawupper=_c76_empty();R.reason=B.reason;}
 else {R.rawlower_state=(R.lower.state==1);R.rawupper_state=(R.upper.state==1);}
}

real scalar _c77_ready() return(7701)
struct _C76_Budget scalar _c77_start() return(_c76_budget())
void _c77_begin(struct _C76_Budget scalar B,real scalar peers)
{
 if(!_c76_begin(B,peers)) _error(3498)
}
void _c77_face(real scalar obs,real scalar face,real scalar face_state,
 real matrix T,real matrix E,real colvector b,real matrix U,real colvector ub,
 real matrix X,real matrix Y,real colvector peers,real scalar tol,real scalar cutoff,
 real scalar stored_link,real scalar stage,real matrix CE,struct _C76_Budget scalar B)
{
 struct _C76_Result scalar R
 struct _C76_Point scalar rawlo,rawhi
 real rowvector num,den
 real scalar j,m,s,extra
 m=cols(X);s=cols(Y);extra=rows(T)-m-s
 for(j=1;j<=rows(X);j++) {
  if(face_state==1) {
   num=(J(1,m,0),Y[j,.],J(1,extra,0));den=(X[j,.],J(1,s+extra,0))
   _c76_peer(num,den,T,E,b,U,ub,cutoff,tol,1,B,R,rawlo,rawhi)
  }
  else {R=_c76_result();R.state=-3;R.reason=8;}
  // T77 public range observation point.
  CE=CE\(obs,face,peers[j],face_state,R.state,R.denstate,R.cutstate,
   R.lower.value,R.upper.value,R.lower.state,R.upper.state,R.lower.display,R.upper.display,
   R.rawlower_state,R.rawupper_state,stored_link,stage,R.reason,B.rowcalls,B.rowops)
 }
}
void _c77_finish(real matrix CE,real scalar requested,struct _C76_Budget scalar B,
 string scalar matrixname,string scalar sumname)
{
 real rowvector S
 S=J(1,20,0);S[1]=requested;S[19..20]=(B.calls,B.ops)
 if(rows(CE)==1 & CE[1,2]==0) S[7]=requested
 else {
  S[2]=sum((CE[.,5]:==1):|(CE[.,5]:==2));S[3]=sum(CE[.,5]:==1);S[4]=sum(CE[.,5]:==2)
  S[5]=sum(CE[.,5]:==0);S[6]=sum(CE[.,5]:==-1);S[7]=sum(CE[.,5]:==-2);S[8]=sum(CE[.,5]:==-3)
  S[9]=sum((CE[.,10]:==2):|(CE[.,11]:==2));S[10]=sum((CE[.,10]:==3):|(CE[.,11]:==3))
  S[11]=sum(CE[.,6]:==3);S[12]=sum(CE[.,7]:==3)
  S[13]=sum((CE[.,12]:==3):|(CE[.,12]:==4):|(CE[.,13]:==3):|(CE[.,13]:==4))
  S[14]=sum((CE[.,12]:==5):|(CE[.,13]:==5));S[15]=sum(CE[.,16]:==2)
  S[16]=sum(CE[.,16]:==3);S[17]=sum(CE[.,16]:==-1)
  S[18]=sum((CE[.,14]:==-1):|(CE[.,15]:==-1))
 }
 st_matrix(matrixname,CE);st_matrix(sumname,S)
}

end
