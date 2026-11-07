*! version 1.0.0 06nov2026
version 16.0
quietly findfile _dea_cooper_live_core.ado
quietly run "`r(fn)'"
capture mata: mata drop _c53_data_ready()
capture mata: mata drop _c53_compute()
capture mata: mata drop _c53_pipeline()
capture mata: mata drop _c53_aux()
capture mata: mata drop _c53_model()
capture mata: mata drop _C53_Live()
mata:
mata set matastrict on
struct _C53_Live {
 struct _C51_Result scalar point
 real rowvector pmin,pmax,cmin,cmax,state,so,tx,ty
 real scalar complete,conversion,rounded
 string scalar range_reason
}
// Construct the same binary64 original LP as the public Cooper equations.
// No computed optimum or rounded projection enters a subsequent model RHS.
void _c53_model(real matrix X,real matrix Y,real scalar index,string scalar weight,
 real rowvector c,real matrix E,real colvector b,real matrix U,real colvector ub,
 real rowvector lo,real rowvector hi,real rowvector cs,real rowvector ct,
 real rowvector sx,real rowvector sy)
{
 real matrix xs,ys
 real rowvector w
 real scalar n,m,s,v
 n=rows(X);m=cols(X);s=cols(Y);v=n+m+s+1
 sx=colmax(X);sy=colmax(Y);xs=X:/sx;ys=Y:/sy
 E=(xs',I(m),J(m,s,0),J(m,1,0))\(ys',J(s,m,0),-I(s),-ys[index,.]')\(J(1,n,1),J(1,m+s+1,0))
 b=xs[index,.]'\J(s,1,0)\1;U=J(0,v,0);ub=J(0,1,0)
 lo=(J(1,n+m+s,0),1);hi=(J(1,n,1),xs[index,.],J(1,s+1,.));c=(J(1,v-1,0),1)
 w=(sx,sy);if(weight=="normalized") w=1:/(xs[index,.],ys[index,.])
 w=w:/max(w);cs=(J(1,n,0),w,0)
 w=sx;if(weight=="normalized") w=1:/xs[index,.]
 w=w:/max(w);ct=(xs*w')'
}
void _c53_aux(struct _C44_Cert scalar base,struct _C44_Cert scalar chosen,
 real matrix xs,real matrix ys,real rowvector sx,real rowvector sy,
 real scalar n,real scalar m,real scalar s,struct _C53_Live scalar L)
{
 struct _C44_Cert scalar point
 struct _C44_Work scalar W
 real scalar j,v
 real rowvector N,D,q
 point=base;v=rows(base.numerator)
 if(L.point.applied) {point.numerator=chosen.numerator[1..v,.];point.denominator=chosen.denominator;}
 D=point.denominator;W.ops=0;W.limit=10000000;W.limbs=512;W.maxlimbs=0;W.error=0
 for(j=1;j<=m;j++) {
  q=(xs[,j]',J(1,v-n,0))
  L.tx[j]=_c52_endpoint(q,point,1,0,0,sx[j],L.point)
 }
 for(j=1;j<=s;j++) {
  L.so[j]=_c51_scaled(_c44_norm(point.numerator[n+m+j,.]),D,sy[j],W,L.point)
  q=(ys[,j]',J(1,v-n,0))
  L.ty[j]=_c52_endpoint(q,point,1,0,0,sy[j],L.point)
 }
 if(W.error) L.point.conversion=2
 if(L.point.conversion) {
  L.point.lambda=L.point.mu=J(1,n,.);L.point.scores=J(1,5,.)
  L.point.total=L.point.technical=L.point.congestion=J(1,m,.)
  L.point.usable=L.point.applied=L.point.rounded=0;L.point.distance=.
  L.point.source="unavailable";L.point.reason="conversion_unavailable"
  L.tx=J(1,m,.);L.so=L.ty=J(1,s,.)
 }
}
struct _C53_Live scalar _c53_pipeline(real rowvector c,real matrix E,real colvector b,real matrix U,real colvector ub,real rowvector lo,real rowvector hi,real rowvector cs,real rowvector ct,real scalar n,real scalar m,real scalar s,
 real scalar index,real rowvector previous_lambda,real rowvector previous_mu,
 string scalar policy,real scalar maxpivots,real rowvector xu,real rowvector yu,real scalar ranges)
{
 struct _C53_Live scalar L
 struct _C51_Result scalar RR
 struct _C44_Cert scalar fit,base,chosen
 real matrix G,Fr,Fs,Fa,Ua,B,xs,ys,Ed,Ud,Fd,Q,BG
 real colvector h,br,bs,ba,uba,hB,bd,ubd,bf,Bh
 real rowvector la,ha,ca,cd,ld,hd,units
 real scalar v,v1,j,signum,count,ok,value
 string scalar tag
 L.point=_c52_empty(n,m);L.state=J(1,6,0);L.pmin=L.pmax=J(1,m+s,.);L.cmin=L.cmax=J(1,m,.)
 L.tx=J(1,m,.);L.so=L.ty=J(1,s,.)
 L.complete=L.conversion=L.rounded=0;L.range_reason=(ranges ? "pipeline_unavailable" : "not_requested")
 v=n+m+s+1
 if(rows(c)!=1 | cols(c)!=v | rows(E)!=m+s+1 | cols(E)!=v | rows(b)!=rows(E) | cols(b)!=1 | cols(U)!=v | rows(ub)!=rows(U) | cols(ub)!=1 | rows(lo)!=1 | cols(lo)!=v | rows(hi)!=1 | cols(hi)!=v | rows(cs)!=1 | cols(cs)!=v | rows(ct)!=1 | cols(ct)<n) {L.point.reason="input_model_dimensions";return(L);}
 xs=E[1..m,1..n]';ys=E[m+1..m+s,1..n]'
 if(hasmissing(c) | hasmissing(E) | hasmissing(b) | hasmissing(U) | hasmissing(ub) | hasmissing(lo) | hasmissing(cs) | hasmissing(ct) | any(vec(xs):<=0) | any(vec(ys):<=0)) {L.point.reason="input_model_values";return(L);}
 ct=ct[cols(ct)-n+1..cols(ct)]
 ok=_c52_step("radial",c,E,b,U,ub,lo,hi,fit,G,h);L.state[1]=ok
 if(!ok) {_c52_save("pipeline.bin",L.state);return(L);}
 _c44_face(E,b,G,h,fit,Fr,br)
 ok=_c52_step("slack",cs,Fr,br,U,ub,lo,hi,fit,G,h);L.state[2]=ok
 if(!ok) {_c52_save("pipeline.bin",L.state);return(L);}
 _c44_face(Fr,br,G,h,fit,Fs,bs);la=lo;ha=hi
 if(!hasmissing(previous_lambda)) {
  Fa=(Fs,J(rows(Fs),1,0));ba=bs;Q=(I(n),J(n,v-n,0))
  Ua=(Q,J(n,1,-1))\(-Q,J(n,1,-1));uba=previous_lambda'\(-previous_lambda')
  la=(lo,0);ha=(hi,.);ca=(J(1,v,0),-1)
  ok=_c52_step("anchor",ca,Fa,ba,Ua,uba,la,ha,fit,G,h);L.state[3]=ok
  if(!ok) {_c52_save("pipeline.bin",L.state);return(L);}
 }
 else L.state[3]=2
 B=G[fit.basis,.];hB=h[fit.basis];v1=cols(B)
 Ed=(B,J(v1,n,0))\(-ys',J(s,v1-n,0),ys')\(J(1,v1,0),J(1,n,1));bd=hB\J(s,1,0)\1
 Ud=(J(m,v1,0),xs')\(xs',J(m,v1-n,0),-xs');ubd=b[1..m]\J(m,1,0)
 ld=(la,J(1,n,0));hd=(ha,J(1,n,1));cd=(J(1,v1,0),ct)
 ok=_c52_step("decomposition",cd,Ed,bd,Ud,ubd,ld,hd,fit,G,h);L.state[4]=ok
 if(!ok) {_c52_save("pipeline.bin",L.state);return(L);}
 base=fit;BG=G;Bh=h
 _c44_face(Ed,bd,G,h,fit,Fd,bf)
 chosen=base
 if(policy=="stage2") _c52_select("stage2_anchor",cd,Ed,bd,Ud,ubd,ld,hd,n,previous_mu,base,G,h,chosen,maxpivots)
 L.point=_c52_point(cd,BG,Bh,rows(Ed),base,chosen,xs,index,n,m,s,xu,policy,previous_mu,maxpivots)
 if(L.point.usable) _c53_aux(base,chosen,xs,ys,xu,yu,n,m,s,L)
 if(!ranges) {_c52_save("pipeline.bin",L.state);return(L);}
 RR=_c52_empty(n,m);units=(xu,yu)
 Q=(xs',J(m,v-n,0))\(ys',J(s,v-n,0));count=0
 for(j=1;j<=rows(Q);j++) for(signum=-1;signum<=1;signum=signum+2) {
  tag="projection_"+strofreal(j)+"_"+strofreal(signum)
  ok=_c52_step(tag,signum*Q[j,.],Fs,bs,U,ub,lo,hi,fit,G,h);count=count+ok
  if(ok) {
   value=_c52_endpoint(signum*Q[j,.],fit,signum,0,0,units[j],RR)
   if(signum<0) L.pmin[j]=value
   else L.pmax[j]=value
  }
 }
 L.state[5]=count;Q=(J(m,v1,0),xs');count=0
 for(j=1;j<=m;j++) for(signum=-1;signum<=1;signum=signum+2) {
  tag="congestion_"+strofreal(j)+"_"+strofreal(signum)
  ok=_c52_step(tag,signum*Q[j,.],Fd,bf,Ud,ubd,ld,hd,fit,G,h);count=count+ok
  if(ok) {
   value=_c52_endpoint(signum*Q[j,.],fit,signum,1,xs[index,j],xu[j],RR)
   if(signum>0) L.cmin[j]=value
   else L.cmax[j]=value
  }
 }
 L.state[6]=count;L.conversion=RR.conversion;L.range_reason="range_certificate_unavailable"
 if(RR.conversion) L.range_reason="range_conversion_unavailable"
 if(L.state[5]==2*(m+s) & count==2*m & !RR.conversion) {L.complete=1;L.rounded=RR.rounded;L.range_reason="certified";}
 else {L.pmin=L.pmax=J(1,m+s,.);L.cmin=L.cmax=J(1,m,.);}
 _c52_save("pipeline.bin",L.state);return(L)
}
void _c53_compute(string scalar inputs,string scalar outputs,string scalar touse,
 string scalar resultvars,string scalar weight,real scalar tol,real scalar matrices,
 real scalar ranges,real scalar debug,string scalar destinations,string scalar policy,
 string scalar anchor,string scalar reference,real scalar cap,string scalar auditroot)
{
 struct _C53_Live scalar L
 real matrix X,Y,E,U,R,SI,TI,CI,SO,TX,TY,L1,L2,ST,CS,PM,PMX,CM,CMX,RS,A,M,meta
 real colvector obs,b,ub
 real rowvector c,lo,hi,cs,ct,sx,sy,l1,l2,stage
 real scalar n,m,s,i,fail,eps,eff,multi1,multi2,reason
 string rowvector dest,obsnames
 string scalar audit
 X=st_data(.,tokens(inputs),touse);Y=st_data(.,tokens(outputs),touse);obs=selectindex(st_data(.,touse))
 n=rows(X);m=cols(X);s=cols(Y);eps=min((tol,1e-9));dest=tokens(destinations)
 A=M=J(n,n,.);if(anchor!="") A=st_matrix(anchor);if(reference!="") M=st_matrix(reference)
 R=J(n,6,.);SI=TI=CI=TX=J(n,m,.);SO=TY=J(n,s,.);L1=L2=J(n,n,.)
 ST=J(n,7,0);CS=J(n,18,.);PM=PMX=J(n,m+s,.);CM=CMX=J(n,m,.);RS=J(n,2,.)
 fail=eff=multi1=multi2=0
 for(i=1;i<=n;i++) {
  audit=""
  if(auditroot!="") audit=auditroot+"/dmu_"+strofreal(i)
  st_local("audit",audit)
  _c53_model(X,Y,i,weight,c,E,b,U,ub,lo,hi,cs,ct,sx,sy)
  L=_c53_pipeline(c,E,b,U,ub,lo,hi,cs,ct,n,m,s,i,A[i,.],M[i,.],policy,cap,sx,sy,ranges)
  stage=L.state;reason=0
  if(!stage[1]) reason=1
  else if(!stage[2] | !stage[3]) reason=2
  else if(!stage[4]) reason=3
  else if(!L.point.usable) reason=4
  R[i,.]=(L.point.scores,reason);fail=fail+(reason!=0)
  ST[i,1..3]=(stage[1]!=1,stage[2]!=1 | stage[3]==0,stage[4]!=1)
  CS[i,.]=(L.point.usable,L.point.base_cert,L.point.policy_status,L.point.applied,L.point.selected_cert,L.point.conversion,L.point.rounded,L.point.cap,L.point.distance,L.complete,L.conversion,L.rounded,stage)
  if(L.point.usable) {
   SI[i,.]=L.point.total;TI[i,.]=L.point.technical;CI[i,.]=L.point.congestion
   SO[i,.]=L.so;TX[i,.]=L.tx;TY[i,.]=L.ty;L1[i,.]=L.point.lambda;L2[i,.]=L.point.mu
   if(abs(R[i,2]-1)<=eps & max((SI[i,.]:/(X[i,.]),SO[i,.]:/(Y[i,.])))<=eps) eff++
  }
  if(ranges) {
   PM[i,.]=L.pmin;PMX[i,.]=L.pmax;CM[i,.]=L.cmin;CMX[i,.]=L.cmax;RS[i,.]=J(1,2,!L.complete)
   if(L.complete) {
    multi1=multi1+(max((L.pmax-L.pmin):/(X[i,.],Y[i,.]))>10*eps)
    multi2=multi2+(max((L.cmax-L.cmin):/X[i,.])>10*eps)
   }
  }
  if(debug) printf("  DMU observation %g: result=%g; selection=%g; conversion=%g; ranges=%g\n",obs[i],reason,L.point.policy_status,L.point.conversion,L.complete)
 }
 st_local("audit","");st_store(obs,tokens(resultvars),R)
 meta=(n,fail,eff,sum((R[,1]:>eps):&(R[,1]:<.)),0,multi1,multi2,(ranges ? sum(CS[,10]:!=1) : .),(ranges ? sum(CS[,10]:==1) : .))\J(1,9,0)
 meta=meta\(sum(CS[,3]:==0),sum(CS[,3]:==1),sum(CS[,3]:==2),sum(CS[,3]:==3),sum(CS[,3]:==4),sum(CS[,6]:>0),sum(CS[,15]:==2),0,0)
 st_matrix(dest[1],meta);st_matrix(dest[17],CS)
 if(matrices) {
  st_matrix(dest[2],SI);st_matrix(dest[3],TI);st_matrix(dest[4],CI);st_matrix(dest[5],SO);st_matrix(dest[6],TX);st_matrix(dest[7],TY)
  st_matrix(dest[8],L1);st_matrix(dest[9],L2);st_matrix(dest[10],ST);st_matrix(dest[11],obs);st_matrix(dest[17],CS)
  obsnames="obs":+strofreal(obs')
  for(i=2;i<=17;i++) if(i<12 | i==17) st_matrixrowstripe(dest[i],(J(n,1,""),obsnames'))
  st_matrixcolstripe(dest[8],(J(n,1,""),obsnames'));st_matrixcolstripe(dest[9],(J(n,1,""),obsnames'))
 }
 if(ranges) {
  st_matrix(dest[12],PM);st_matrix(dest[13],PMX);st_matrix(dest[14],CM);st_matrix(dest[15],CMX);st_matrix(dest[16],RS)
  for(i=12;i<=16;i++) st_matrixrowstripe(dest[i],(J(n,1,""),obsnames'))
  st_matrixcolstripe(dest[12],((J(m,1,"input")\J(s,1,"output")),(tokens(inputs),tokens(outputs))'))
  st_matrixcolstripe(dest[13],((J(m,1,"input")\J(s,1,"output")),(tokens(inputs),tokens(outputs))'))
 }
}
real scalar _c53_data_ready() return(5401)

end
