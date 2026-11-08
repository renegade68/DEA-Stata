*! version 1.1.0 08nov2026
*! Failure-only certified original-unit weights and bounded exact peer arithmetic.
version 16.0
capture mata: mata drop _c84_work()
capture mata: mata drop _c84_log()
capture mata: mata drop _c84_scale()
capture mata: mata drop _c84_coefficient()
capture mata: mata drop _c84_build()
capture mata: mata drop _c84_dyad()
capture mata: mata drop _c84_number()
capture mata: mata drop _c84_abs()
capture mata: mata drop _c84_sub()
capture mata: mata drop _c84_div()
capture mata: mata drop _c84_gt()
capture mata: mata drop _c84_value()
capture mata: mata drop _c84_dot()
capture mata: mata drop _c84_objective()
capture mata: mata drop _c84_point()
capture mata: mata drop _c84_guard()
capture mata: mata drop _c84_ready()
capture mata: mata drop _c84_mul()
capture mata: mata drop _c84_peer()
capture mata: mata drop _c84_multiplier()
capture mata: mata drop _C84_Model()
capture mata: mata drop _C84_Q()
mata:
mata set matastrict on
struct _C84_Model {
 real rowvector c,qc,ve,re,lo,hi
 real matrix E,U,qE,qU,S
 real colvector b,ub,qb,qub
 real scalar oe,ok,ops
}
struct _C84_Q {
 real rowvector N,D
}
struct _C58_Work scalar _c84_work()
{
 struct _C58_Work scalar W
 W.ops=0;W.limit=10000000;W.limbs=512;W.maxlimbs=0;W.error=0
 return(W)
}
real scalar _c84_log(real scalar v0)
{
 real scalar e,v
 e=0;v=abs(v0)
 if(v==0 | missing(v)) return(.)
 while(v>=2) {v=v/2;e++;}
 while(v<1) {v=v*2;e--;}
 return(e)
}
real scalar _c84_scale(real scalar v0,real scalar e,real scalar ops)
{
 real scalar i,v
 v=v0
 for(i=1;i<=abs(e);i++) {
  v=v*(e>0 ? 2 : .5);ops++
  if(missing(v) | ops>2000000) return(.)
 }
 return(v)
}
real scalar _c84_coefficient(real scalar v,real scalar e,real scalar ops)
{
 real scalar q,back
 q=_c84_scale(v,e,ops)
 if(missing(q) | (v!=0 & q==0)) return(.)
 back=_c84_scale(q,-e,ops)
 return(back==v ? q : .)
}
struct _C84_Model scalar _c84_build(real matrix X,real matrix Y,real scalar index,
 real scalar free,real scalar out,real scalar variant)
{
 struct _C84_Model scalar A
 real matrix R,B
 real rowvector ex
 real colvector h,col
 real scalar i,j,d,e,anchor
 d=cols(X)+cols(Y)+free
 A.ok=0;A.ops=0;A.ve=J(1,d,0);R=(X,Y)
 A.c=(out ? (X[index,.],J(1,cols(Y),0)) : (J(1,cols(X),0),-Y[index,.]))
 A.E=(out ? (J(1,cols(X),0),Y[index,.]) : (X[index,.],J(1,cols(Y),0)))
 A.U=(-X,Y)
 if(free) {A.c=(A.c,(out ? 1 : -1));A.E=(A.E,0);A.U=(A.U,J(rows(X),1,out ? -1 : 1));}
 A.b=1;A.ub=J(rows(X),1,0);A.lo=J(1,d,0);A.hi=J(1,d,.);if(free) A.lo[d]=.
 A.S=I(d);if(free) A.S=(A.S,-I(d)[.,d])
 if(rows(X)>32 | d+free>16 | hasmissing(R) | min(R)<0) return(A)
 for(j=1;j<=cols(R);j++) {
  col=select(R[.,j],R[.,j]:>0)
  if(!rows(col)) continue
  if(variant==1) {anchor=R[index,j];if(anchor==0) anchor=max(col);A.ve[j]=-_c84_log(anchor);}
  else {
   ex=J(1,rows(col),0)
   for(i=1;i<=rows(col);i++) ex[i]=_c84_log(col[i])
   A.ve[j]=-floor((min(ex)+max(ex))/2)
  }
 }
 ex=J(1,0,.)
 for(j=1;j<=d;j++) if(A.c[j]!=0) ex=(ex,_c84_log(A.c[j])+A.ve[j])
 A.oe=(cols(ex) ? -max(ex) : 0)
 A.qc=J(1,d,.)
 for(j=1;j<=d;j++) A.qc[j]=_c84_coefficient(A.c[j],A.ve[j]+A.oe,A.ops)
 B=A.E\A.U;h=A.b\A.ub;A.re=J(1,rows(B),0)
 for(i=1;i<=rows(B);i++) {
  ex=J(1,0,.)
  for(j=1;j<=d;j++) if(B[i,j]!=0) ex=(ex,_c84_log(B[i,j])+A.ve[j])
  if(h[i]!=0) ex=(ex,_c84_log(h[i]))
  e=(cols(ex) ? -floor((max(ex)+min(ex))/2) : 0);A.re[i]=e
  for(j=1;j<=d;j++) B[i,j]=_c84_coefficient(B[i,j],A.ve[j]+e,A.ops)
  h[i]=_c84_coefficient(h[i],e,A.ops)
 }
 A.qE=B[1,.];A.qU=B[2..rows(B),.];A.qb=h[1];A.qub=h[2..rows(B)]
 A.ok=!(hasmissing(B) | hasmissing(h) | hasmissing(A.qc))
 return(A)
}
struct _C84_Q scalar _c84_dyad(real rowvector N,real scalar e,struct _C58_Work scalar W)
{
 struct _C84_Q scalar q
 q.N=_c58b_shift(N,max((0,e)),W);q.D=_c58b_shift((1,1),max((0,-e)),W)
 return(q)
}
struct _C84_Q scalar _c84_number(real scalar v,struct _C58_Work scalar W)
{
 real rowvector N
 real scalar e
 N=_c70_decode(v,e,W)
 return(_c84_dyad(N,e,W))
}
struct _C84_Q scalar _c84_abs(struct _C84_Q scalar q)
{
 struct _C84_Q scalar z
 z=q;z.N[1]=abs(z.N[1]);return(z)
}
struct _C84_Q scalar _c84_sub(struct _C84_Q scalar a,struct _C84_Q scalar b,struct _C58_Work scalar W)
{
 struct _C84_Q scalar q
 q.N=_c58b_add(_c58b_mul(a.N,b.D,W),_c58b_neg(_c58b_mul(b.N,a.D,W)),W)
 q.D=_c58b_mul(a.D,b.D,W);return(q)
}
struct _C84_Q scalar _c84_div(struct _C84_Q scalar a,struct _C84_Q scalar b,struct _C58_Work scalar W)
{
 struct _C84_Q scalar q
 q.N=_c58b_mul(a.N,b.D,W);q.D=_c58b_mul(a.D,b.N,W)
 if(q.D[1]<0) {q.D[1]=-q.D[1];q.N[1]=-q.N[1];}
 if(q.D[1]==0) W.error=4
 return(q)
}
real scalar _c84_gt(struct _C84_Q scalar a,struct _C84_Q scalar b,struct _C58_Work scalar W)
{
 struct _C84_Q scalar z
 z=_c84_sub(a,b,W);return(z.N[1]>0)
}
real scalar _c84_value(struct _C84_Q scalar q)
{
 real rowvector z
 z=_c58d_double(q.N,q.D);return(z[2] ? . : z[1])
}
struct _C84_Q scalar _c84_dot(real rowvector a,real rowvector b,real scalar absolute,struct _C58_Work scalar W)
{
 transmorphic scalar terms
 real rowvector ex,N,aa,bb
 real scalar i,e,f,emin
 terms=asarray_create("real",1);ex=J(1,cols(a),0);emin=.
 for(i=1;i<=cols(a);i++) {
  N=(0,0)
  if(a[i]!=0 & b[i]!=0) {
   aa=_c70_decode(a[i],e,W);bb=_c70_decode(b[i],f,W)
   N=_c58b_mul(aa,bb,W);if(absolute) N[1]=abs(N[1])
   ex[i]=e+f;emin=min((emin,e+f))
  }
  asarray(terms,i,N)
 }
 if(emin>=.) emin=0
 N=(0,0)
 for(i=1;i<=cols(a);i++) {
  aa=asarray(terms,i)
  if(aa[1]!=0) N=_c58b_add(N,_c58b_shift(aa,ex[i]-emin,W),W)
 }
 return(_c84_dyad(N,emin,W))
}
// The original primary objective is derived from the new exact certificate.
// No T82 certificate, known score, or Python oracle enters a decision.
struct _C84_Q scalar _c84_objective(struct _C84_Model scalar A,struct _C58_Cert scalar fit,struct _C58_Work scalar W)
{
 real rowvector c,ex,a,N
 real scalar i,e,emin
 transmorphic scalar terms
 struct _C84_Q scalar q
 c=A.qc*A.S;ex=J(1,cols(c),0);emin=.
 terms=asarray_create("real",1)
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
 e=emin-A.oe;q.N=_c58b_shift(N,max((0,e)),W)
 q.D=_c58b_shift(fit.denominator,max((0,-e)),W);return(q)
}
real rowvector _c84_point(real rowvector ve,real scalar free,struct _C58_Cert scalar fit,real rowvector diag)
{
 struct _C58_Work scalar W
 real rowvector x,N,D,z
 real scalar j,d
 d=cols(ve);x=J(1,d,.);diag=J(1,4,0)
 for(j=1;j<=d;j++) {
  W=_c84_work();N=_c58b_norm(fit.numerator[j,.]);D=fit.denominator
  if(free & j==d) N=_c58b_add(N,_c58b_neg(_c58b_norm(fit.numerator[d+1,.])),W)
  N=_c58b_shift(N,max((0,ve[j])),W);D=_c58b_shift(D,max((0,-ve[j])),W)
  z=_c58d_double(N,D);diag[3]=diag[3]+W.ops;diag[4]=max((diag[4],W.maxlimbs))
  if(W.error | z[2] | (N[1]!=0 & z[1]==0)) {diag[2]=(W.error ? 10+W.error : z[2] ? z[2] : 4);return(J(1,d,.));}
  x[j]=z[1]
 }
 diag[1]=1;return(x)
}
real rowvector _c84_guard(struct _C84_Model scalar A,real rowvector x,
 struct _C84_Q scalar opt,real scalar tol)
{
 struct _C58_Work scalar W
 struct _C84_Q scalar one,zero,qt,q,r,den,diff,scale
 real scalar i,ok,eq,iq,gap,bounds
 W=_c84_work();one=_c84_number(1,W);zero=_c84_number(0,W);qt=_c84_number(tol,W)
 if(hasmissing(x)) return((0,.,.,.,0,4,0,0))
 q=_c84_dot(A.E,x,0,W);diff=_c84_abs(_c84_sub(q,one,W))
 eq=_c84_value(diff);ok=!_c84_gt(diff,qt,W);iq=0
 for(i=1;i<=rows(A.U);i++) {
  q=_c84_dot(A.U[i,.],x,0,W)
  if(q.N[1]<=0) continue
  den=_c84_dot(A.U[i,.],x,1,W);r=_c84_div(q,den,W)
  iq=max((iq,_c84_value(r)));if(_c84_gt(r,qt,W)) ok=0
 }
 bounds=all(select(x,A.lo:<.):>=0)
 scale=_c84_abs(opt);if(!_c84_gt(scale,one,W)) scale=one
 q=_c84_dot(A.c,x,0,W);diff=_c84_abs(_c84_sub(q,opt,W));r=_c84_div(diff,scale,W)
 gap=_c84_value(r);if(_c84_gt(r,qt,W)) ok=0
 return((ok & bounds & !W.error,eq,iq,gap,bounds,W.error,W.ops,W.maxlimbs))
}
real scalar _c84_ready() return(8403)
struct _C84_Q scalar _c84_mul(struct _C84_Q scalar a,struct _C84_Q scalar b,struct _C58_Work scalar W)
{
 struct _C84_Q scalar q
 q.N=_c58b_mul(a.N,b.N,W);q.D=_c58b_mul(a.D,b.D,W);return(q)
}
// Peer status:1finite,2denominator policy,3export range/loss,
// 4arithmetic limit,5invalid operands,6Lim-Zhu appraisal outside [0,1+tol].
real scalar _c84_peer(real rowvector v,real rowvector x,real rowvector u,real rowvector y,
 real scalar w,real scalar free,real scalar lz,real scalar tol,real rowvector diag)
{
 struct _C58_Work scalar W
 struct _C84_Q scalar vx,uy,qw,den,cut,q,zero,scale
 real scalar z,status
 diag=(5,0,0)
 if(hasmissing((v,x,u,y,w,tol)) | min((v,x,u,y))<0 | tol<=0) return(.)
 W=_c84_work();vx=_c84_dot(v,x,0,W);uy=_c84_dot(u,y,0,W)
 den=vx;zero=_c84_number(0,W);cut=_c84_number(free ? tol : 0,W)
 if(lz) {
  qw=_c84_number(w,W);den=_c84_sub(vx,qw,W)
  scale=_c84_abs(qw)
  if(_c84_gt(vx,scale,W)) scale=vx
  if(_c84_gt(uy,scale,W)) scale=uy
  cut=_c84_mul(_c84_number(tol,W),scale,W)
 }
 status=2;z=.
 if(_c84_gt(den,zero,W) & _c84_gt(den,cut,W)) {
  q=_c84_div(uy,den,W)
  if(lz & _c84_gt(q,_c84_number(1+tol,W),W)) status=6
  else {
   z=_c84_value(q);status=(missing(z) | (q.N[1]!=0 & z==0) ? 3 : 1)
   if(status!=1) z=.
  }
 }
 if(W.error) {status=4;z=.;}
 diag=(status,W.ops,W.maxlimbs);return(z)
}
// Failure-only raw primary recovery. Caller preserves accepted evaluators.
// One exact primary fit, at most two numerical proposals, no oracle input.
void _c84_multiplier(real matrix X,real matrix Y,real scalar index,real scalar free,
 real scalar out,string scalar method,real scalar tol,real rowvector x,
 real scalar score,real scalar stage,real rowvector rd,real rowvector scaling,
 real rowvector pd,real rowvector sd,real rowvector pt,real rowvector st,
 real rowvector link,real rowvector od)
{
 struct _C84_Model scalar A,last_valid
 struct _C58_Cert scalar p,s
 struct _C58_Work scalar W
 struct _C84_Q scalar opt
 transmorphic scalar exactC
 real rowvector sol,diag,g,ed,c2,trial,check
 real matrix G,G2,F,Traw
 real colvector h,h2,fb
 real scalar attempt,rc,obj,ops,j,ce,d,decoder_extended,chosen,last_rep
 real colvector decoder_values
 d=cols(X)+cols(Y)+free;x=J(1,d,.);score=.;stage=0
 rd=(1,0,0,.,.,0,.,.,0,0,6,.,.,.,0,0,0)
 scaling=J(1,0,.)
 pd=sd=(0,.,.,.,.,.,.,.,.);pt=st=_c62_empty_trace()
 link=(0,.,0,0,.,.,.,.,.,.);od=J(1,6,0)
 if(rows(X)>32 | d+free>16) {rd[11]=10;return;}
 scaling=J(1,d+rows(X)+2,.);chosen=0;last_rep=0
 for(attempt=1;attempt<=2;attempt++) {
  A=_c84_build(X,Y,index,free,out,attempt);rd[2]=rd[2]+1;rd[15]=rd[15]+A.ops
  if(!A.ok) {rd[11]=5;continue;}
  last_valid=A;last_rep=attempt
  rc=_c57_solve(A.qc,A.qE,A.qb,A.qU,A.qub,A.lo,A.hi,tol,sol,obj,diag)
  if(rc!=0 | hasmissing(sol)) {rd[11]=6;continue;}
  trial=sol;ops=0
  for(j=1;j<=d;j++) {
   trial[j]=_c84_scale(sol[j],A.ve[j],ops)
   if(sol[j]!=0 & trial[j]==0) trial[j]=.
  }
  if(hasmissing(trial)) {rd[11]=8;rd[17]=1;continue;}
  W=_c84_work();opt=_c84_dot(A.c,trial,0,W);g=_c84_guard(A,trial,opt,tol)
  if(!g[1] | W.error) {rd[11]=9;continue;}
  chosen=attempt;break
 }
 // Existing successful proposals keep their representation and selected path.
 // After both proposals fail, certify the last valid exact representation once.
 // No numerical point, known score, or data-dependent special case is required.
 if(!chosen) {
  if(!last_rep) return
  A=last_valid;chosen=last_rep
 }
 rd[3]=chosen;scaling=(A.oe,A.ve,A.re)
  _c58_fit(A.qc*A.S,A.qE*A.S,A.qb,A.qU*A.S,A.qub,tol,p,G,h,1)
  rd[4..5]=(p.status,p.reason);pd=_c60_certdiag(p,1);pt=p.trace62
  link[1]=p.status;link[5]=p.ops;link[7]=p.solver_status;link[9]=p.reason
  if(p.status!=1) {rd[11]=7;return;}
  // Broader secondary decoding only when this primary model was outside the
  // old coefficient domain. Existing recovered primary/secondary rows keep it.
  decoder_values=abs(vec(G)\h\(A.qc*A.S)')
  decoder_extended=any((decoder_values:>0):&(decoder_values:<1.1830521861667747e-271)) | any(decoder_values:>8.452712498170644e270)

  trial=_c84_point(A.ve,free,p,ed);rd[16]=rd[16]+ed[3]
  if(!ed[1]) {rd[11]=8;rd[17]=1;return;}
  W=_c84_work();opt=_c84_objective(A,p,W)
  // Convert the reciprocal directly for output orientation.
  obj=_c84_value(out ? _c84_div(_c84_number(1,W),opt,W) : opt)
  score=(out ? obj : -obj);g=_c84_guard(A,trial,opt,tol);rd[12..14]=g[2..4]
  if(!g[1] | W.error | missing(score) | score<=0) {rd[11]=9;score=.;return;}
  x=trial;rd[6]=1;stage=1;rd[10]=1;rd[11]=1;link[4]=4
  if(method!="arbitrary") {
   stage=3;rd[10]=3;rd[11]=3;F=A.qE*A.S;fb=A.qb
   if(rows(p.face)) {F=F\G[p.face,.];fb=fb\h[p.face];}
   link[10]=rows(p.face);Traw=A.S;ops=0
   for(j=1;j<=rows(Traw);j++) Traw[j,.]=Traw[j,.]*_c84_scale(1,A.ve[j],ops)
   if(hasmissing(Traw)) {rd[11]=12;return;}
   _c70_build(X,Y,index,out,method=="aggressive",Traw,c2,exactC,ce,od)
   if(od[2]!=1) {rd[11]=12;return;}
   _c70_fit(c2,F,fb,A.qU*A.S,A.qub,tol,exactC,ce,od,s,G2,h2,decoder_extended)
   rd[7..8]=(s.status,s.reason);sd=_c60_certdiag(s,rows(F));st=s.trace62
   link[2]=s.status;link[6]=s.ops;link[8]=s.solver_status
   if(s.trace62[16]==1) rd[11]=4
   if(s.status==1) {
    trial=_c84_point(A.ve,free,s,ed);rd[16]=rd[16]+ed[3]
    g=_c84_guard(A,trial,opt,tol)
    if(ed[1] & g[1]) {
     x=trial;rd[9]=1;stage=2;rd[10]=2;rd[11]=2;rd[12..14]=g[2..4]
     link[3]=1;link[4]=1;od[4]=1
    }
    else {rd[11]=11;link[4]=5;}
   }
  }
}

end
