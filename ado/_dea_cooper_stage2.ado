*! version 1.0.0 06nov2026
*! Internal optional nearest-reference selector on the original technical face.
version 16.0
mata:
mata set matastrict on
// A verified decomposition vertex and its exact largest reference discrepancy
// define a feasible lifted vertex. No oracle basis or rounded RHS is used.
struct _C44_Cert scalar _c50_lift(string scalar tag,real rowvector c,
 real matrix SG,real colvector sh,real scalar newne,real scalar oldne,
 real scalar oldnu,real scalar n,real rowvector reference,
 struct _C44_Cert scalar base)
{
 struct _C44_Work scalar W
 struct _C44_Cert scalar out
 real rowvector a,b,q,best,left,right
 real colvector ids
 real scalar v,j,e,pow,bestpow,bestrow,row,cmp,oldid,newnu
 v=cols(c)-1;newnu=oldnu+2*n
 W.ops=0;W.limit=1000000;W.limbs=128;W.maxlimbs=0;W.error=0
 best=(0,0);bestpow=0;bestrow=0
 for(j=1;j<=n;j++) {
  a=_c44_decode(reference[j],e,W);pow=max((0,-e))
  b=_c44_shift(_c44_mul(base.denominator,a,W),max((0,e)),W)
  q=_c44_add(_c44_shift(_c44_norm(base.numerator[v-n+j,.]),pow,W),_c44_neg(b),W)
  row=(q[1]<0 ? n+j : j);if(q[1]<0) q=_c44_neg(q)
  left=_c44_shift(q,bestpow,W);right=_c44_shift(best,pow,W);cmp=_c44_cmp(left,right)
  if(bestrow==0 | cmp>0 | (cmp==0 & row<bestrow)) {best=q;bestpow=pow;bestrow=row;}
 }
 ids=(newne ? (1::newne) : J(0,1,.))
 for(j=1;j<=rows(base.basis);j++) {
  oldid=base.basis[j]
  if(oldid<=oldne) row=oldid
  else if(oldid<=oldne+oldnu) row=newne+oldid-oldne
  else if(oldid<=oldne+oldnu+v) row=newne+newnu+oldid-oldne-oldnu
  else row=newne+newnu+v+1+oldid-oldne-oldnu-v
  ids=ids\row
 }
 ids=ids\(newne+oldnu+bestrow)
 savebin(tag+"_lift_plan.bin",(bestrow,W.ops,W.maxlimbs,W.error))
 if(W.error) ids=J(0,1,.)
 out=_c44_certify(c,SG,sh,newne,ids,10000000)
 out.solver_status=.
 return(out)
}

void _c50_select(string scalar tag,real rowvector c,real matrix E,real colvector b,
 real matrix U,real colvector ub,real rowvector lo,real rowvector hi,real scalar n,
 real rowvector reference,struct _C44_Cert scalar base,real matrix G,real colvector h,
 struct _C44_Cert scalar chosen,real scalar maxpivots)
{
 real matrix F,A,Q,V,SG
 real colvector fb,vb,sh
 real rowvector objective,lower,upper
 real scalar v,status
 chosen=base;v=cols(c)
 if(base.status!=1 | n<1 | n>v | cols(reference)!=n |
  (maxpivots!=16 & maxpivots!=128)) {savebin(tag+"_state.bin",(-1,0));return;}
 if(hasmissing(reference)) {savebin(tag+"_state.bin",(2,0));return;}
 _c44_face(E,b,G,h,base,F,fb)
 Q=(J(n,v-n,0),I(n));A=(F,J(rows(F),1,0))
 V=(U,J(rows(U),1,0))\(Q,J(n,1,-1))\(-Q,J(n,1,-1))
 vb=ub\reference'\(-reference');objective=(J(1,v,0),-1)
 lower=(lo,0);upper=(hi,.)
 _c48_propose(tag,objective,A,fb,V,vb,lower,upper,chosen,SG,sh)
 if(chosen.status==-1 & rows(chosen.basis)==0) chosen=_c50_lift(tag,objective,SG,sh,rows(A),rows(E),rows(U),n,reference,base)
 if(chosen.status==0) _c45_repair(tag,objective,SG,sh,rows(A),chosen,maxpivots,80000000)
 _c44_emit(tag,objective,A,fb,V,vb,lower,upper,chosen)
 savebin(tag+"_state.bin",(chosen.status==1,maxpivots))
}
void _c50_pipeline(string scalar source,real scalar n,real scalar m,real scalar s,
 real rowvector previous_lambda,real rowvector previous_mu,real scalar maxpivots)
{
 struct _C44_Cert scalar fit,chosen
 real matrix E,U,G,Fr,Fs,Fa,Ua,B,xs,ys,Ed,Ud,Fd,Q
 real colvector b,ub,h,br,bs,ba,uba,hB,bd,ubd,bf
 real rowvector c,lo,hi,cs,ct,la,ha,ca,cd,ld,hd,state
 real scalar v,v1,j,signum,count,ok
 string scalar tag
 c=loadbin(source+"/radial/c.bin");E=loadbin(source+"/radial/E.bin");b=loadbin(source+"/radial/b.bin")
 U=loadbin(source+"/radial/U.bin");ub=loadbin(source+"/radial/ub.bin")
 lo=loadbin(source+"/radial/lo.bin");hi=loadbin(source+"/radial/hi.bin")
 cs=loadbin(source+"/slack/c.bin");ct=loadbin(source+"/decomposition/c.bin");ct=ct[cols(ct)-n+1..cols(ct)]
 v=cols(c);xs=E[1..m,1..n]';ys=E[m+1..m+s,1..n]';state=J(1,6,0)
 ok=_c44_step("radial",c,E,b,U,ub,lo,hi,fit,G,h);state[1]=ok
 if(!ok) {savebin("pipeline.bin",state);return;}
 _c44_face(E,b,G,h,fit,Fr,br)
 ok=_c44_step("slack",cs,Fr,br,U,ub,lo,hi,fit,G,h);state[2]=ok
 if(!ok) {savebin("pipeline.bin",state);return;}
 _c44_face(Fr,br,G,h,fit,Fs,bs)
 la=lo;ha=hi
 if(!hasmissing(previous_lambda)) {
  Fa=(Fs,J(rows(Fs),1,0));ba=bs;Q=(I(n),J(n,v-n,0))
  Ua=(Q,J(n,1,-1))\(-Q,J(n,1,-1));uba=previous_lambda'\(-previous_lambda')
  la=(lo,0);ha=(hi,.);ca=(J(1,v,0),-1)
  ok=_c44_step("anchor",ca,Fa,ba,Ua,uba,la,ha,fit,G,h);state[3]=ok
  if(!ok) {savebin("pipeline.bin",state);return;}
 }
 else state[3]=2
 // The exact inverse identity and original primal proof establish this unique
 // selected point. Other stage-1 constraints are implied, not numerically relaxed.
 B=G[fit.basis,.];hB=h[fit.basis];v1=cols(B)
 Ed=(B,J(v1,n,0))\(-ys',J(s,v1-n,0),ys')\(J(1,v1,0),J(1,n,1))
 bd=hB\J(s,1,0)\1
 Ud=(J(m,v1,0),xs')\(xs',J(m,v1-n,0),-xs')
 ubd=b[1..m]\J(m,1,0);ld=(la,J(1,n,0));hd=(ha,J(1,n,1));cd=(J(1,v1,0),ct)
 ok=_c44_step("decomposition",cd,Ed,bd,Ud,ubd,ld,hd,fit,G,h);state[4]=ok
 if(!ok) {savebin("pipeline.bin",state);return;}
 _c44_face(Ed,bd,G,h,fit,Fd,bf)
 _c50_select("stage2_anchor",cd,Ed,bd,Ud,ubd,ld,hd,n,previous_mu,fit,G,h,chosen,maxpivots)
 Q=(xs',J(m,v-n,0))\(ys',J(s,v-n,0));count=0
 for(j=1;j<=rows(Q);j++) for(signum=-1;signum<=1;signum=signum+2) {
  tag="projection_"+strofreal(j)+"_"+strofreal(signum)
  count=count+_c44_step(tag,signum*Q[j,.],Fs,bs,U,ub,lo,hi,fit,G,h)
 }
 state[5]=count;Q=(J(m,v1,0),xs');count=0
 for(j=1;j<=m;j++) for(signum=-1;signum<=1;signum=signum+2) {
  tag="congestion_"+strofreal(j)+"_"+strofreal(signum)
  count=count+_c44_step(tag,signum*Q[j,.],Fd,bf,Ud,ubd,ld,hd,fit,G,h)
 }
 state[6]=count;savebin("pipeline.bin",state)
}
end
