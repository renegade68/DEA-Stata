*! version 1.1.0 09nov2026
*! Bounded exact game aggregates; decimal rule constants 1/5 and 1/100.
version 16.0
capture mata: mata drop _G129_Derived()
capture mata: mata drop _g129_ready()
capture mata: mata drop _g129_gcd()
capture mata: mata drop _g129_idiv()
capture mata: mata drop _g129_zero()
capture mata: mata drop _g129_norm()
capture mata: mata drop _g129_add()
capture mata: mata drop _g129_sub()
capture mata: mata drop _g129_mul()
capture mata: mata drop _g129_div()
capture mata: mata drop _g129_out()
capture mata: mata drop _g129_bad()
capture mata: mata drop _g129_log()
capture mata: mata drop _g129_coalition()
capture mata: mata drop _g129_nash()
capture mata: mata drop _g130_sum()
capture mata: mata drop _g130_meta()
mata:
mata set matastrict on
struct _G129_Derived {
 real matrix R,Q,why,A
 real rowvector meta
}
real scalar _g129_ready() {
 string rowvector names
 real scalar j
 names=("_g130_sum()","_g130_meta()","_g129_gcd()","_g129_idiv()","_g129_zero()","_g129_norm()","_g129_add()","_g129_sub()","_g129_mul()","_g129_div()","_g129_out()","_g129_bad()","_g129_log()","_g129_coalition()","_g129_nash()","_c107_gcd()","_c107_div()","_c58d_bits()","_c58b_div()")
 for(j=1;j<=cols(names);j++) if(findexternal(names[j])==NULL) return(0)
 return(13001)
}
// Unit gcd/division shortcuts retain input validation and charge the budget.
real rowvector _g129_gcd(real rowvector a,real rowvector b,struct _C58_Work scalar W) {
 if(!_c107_valid(a,W) | !_c107_valid(b,W)) return((0,0))
 if((cols(a)==2 & abs(a[1])==1 & a[2]==1) | (cols(b)==2 & abs(b[1])==1 & b[2]==1)) return((1,1))
 return(_c107_gcd(a,b,W))
}
real rowvector _g129_idiv(real rowvector a,real rowvector b,struct _C58_Work scalar W) {
 real rowvector q,check
 if(!_c107_valid(a,W) | !_c107_valid(b,W)) return((0,0))
 if(cols(b)==2 & b[1]==1 & b[2]==1) return(a)
 // Existing base-2^26 quotient routine: exact remainder and reconstruction.
 // Reserve a scratch limb; at full width keep aligned binary division.
 if(max((cols(a),cols(b)))<=W.limbs) {
  q=_c58b_div(a,b,W)
  if(!W.error) {
   check=_c58b_mul(q,b,W)
   if(!W.error & _c58b_cmp(check,a)!=0) W.error=3
  }
  if(W.error) return((0,0))
  return(q)
 }
 return(_c107_div(a,b,W))
}
struct _C84_Q scalar _g129_zero() {
 struct _C84_Q scalar q
 q.N=(0,0);q.D=(1,1);return(q)
}
struct _C84_Q scalar _g129_norm(struct _C84_Q scalar a,struct _C58_Work scalar W) {
 real rowvector g
 struct _C84_Q scalar q
 q=a
 if(W.error) return(_g129_zero())
 if(q.D[1]<=0) {W.error=4;return(_g129_zero());}
 if(q.N[1]==0) return(_g129_zero())
 g=_g129_gcd(q.N,q.D,W);q.N=_g129_idiv(q.N,g,W);q.D=_g129_idiv(q.D,g,W)
 return(q)
}
struct _C84_Q scalar _g129_add(struct _C84_Q scalar a,struct _C84_Q scalar b,struct _C58_Work scalar W) {
 struct _C84_Q scalar q
 real rowvector g,ad,bd,h
 if(W.error) return(_g129_zero())
 if(a.N[1]==0) return(b)
 if(b.N[1]==0) return(a)
 g=_g129_gcd(a.D,b.D,W);ad=_g129_idiv(a.D,g,W);bd=_g129_idiv(b.D,g,W)
 q.N=_c58b_add(_c58b_mul(a.N,bd,W),_c58b_mul(b.N,ad,W),W)
 if(q.N[1]==0) return(_g129_zero())
 // Reduced operands imply any remaining common factor divides g.
 h=_g129_gcd(q.N,g,W);q.N=_g129_idiv(q.N,h,W)
 q.D=_c58b_mul(ad,_g129_idiv(b.D,h,W),W)
 return(q)
}
struct _C84_Q scalar _g129_sub(struct _C84_Q scalar a,struct _C84_Q scalar b,struct _C58_Work scalar W) {
 struct _C84_Q scalar neg
 neg=b;neg.N=_c58b_neg(neg.N);return(_g129_add(a,neg,W))
}
struct _C84_Q scalar _g129_mul(struct _C84_Q scalar a,struct _C84_Q scalar b,struct _C58_Work scalar W) {
 struct _C84_Q scalar q
 real rowvector g,h
 if(W.error | a.N[1]==0 | b.N[1]==0) return(_g129_zero())
 g=_g129_gcd(a.N,b.D,W);h=_g129_gcd(b.N,a.D,W)
 q.N=_c58b_mul(_g129_idiv(a.N,g,W),_g129_idiv(b.N,h,W),W)
 q.D=_c58b_mul(_g129_idiv(a.D,h,W),_g129_idiv(b.D,g,W),W)
 return(q)
}
struct _C84_Q scalar _g129_div(struct _C84_Q scalar a,struct _C84_Q scalar b,struct _C58_Work scalar W) {
 struct _C84_Q scalar inv
 if(W.error) return(_g129_zero())
 if(b.N[1]==0) {W.error=4;return(_g129_zero());}
 inv.N=b.D;inv.D=b.N
 if(inv.D[1]<0) {inv.D[1]=1;inv.N[1]=-inv.N[1];}
 return(_g129_mul(a,inv,W))
}
// quality,reason,value,relative export error,export operations,max limbs.
// 0 unavailable,1 numerical,2 exact expression + export checked,3 withheld.
real rowvector _g129_bad(real scalar reason,struct _C58_Work scalar W) return((3,reason,.,.,W.ops,W.maxlimbs))
real rowvector _g129_out(struct _C84_Q scalar q,real scalar tol) {
 struct _C58_Work scalar W
 struct _C84_Q scalar error
 real scalar value,reason,relative
 W=_c84_work();reason=0;relative=.;value=.
 if(q.D[1]<=0 | q.N[1]<0 | tol<=0 | tol>=.) return(_g129_bad(1102,W))
 value=_c84_value(q)
 if(value<0 | value>=.) reason=1102
 else if(q.N[1]!=0 & value==0) reason=1101
 else if(q.N[1]==0) relative=0
 else {
  error=_g129_div(_c84_abs(_g129_sub(_c84_number(value,W),q,W)),q,W)
  relative=_c84_value(error)
  if(_c84_gt(error,_c84_number(tol,W),W) | relative>=.) reason=1103
 }
 if(W.error) reason=1200+W.error
 return(((reason ? 3 : 2),reason,(reason ? . : value),relative,W.ops,W.maxlimbs))
}
// Logarithms are numerical, not rationally certified. Normalize the fraction
// before taking log so a representable log survives a nonrepresentable gain.
real rowvector _g129_log(struct _C84_Q scalar source) {
 struct _C84_Q scalar q
 struct _C58_Work scalar W
 real scalar e,m,v
 W=_c84_work();q=source
 if(q.N[1]==0) return((0,1004,.,.,0,0))
 if(q.N[1]<0 | q.D[1]<=0) return(_g129_bad(1102,W))
 e=_c58d_bits(q.N)-_c58d_bits(q.D)
 if(e>=0) q.D=_c58f_restore(q.D,e,W)
 else q.N=_c58f_restore(q.N,-e,W)
 m=_c84_value(q);v=ln(m)+e*ln(2)
 if(W.error) return(_g129_bad(1200+W.error,W))
 if(m<=0 | m>=. | v>=.) return(_g129_bad(1102,W))
 return((1,1005,v,.,W.ops,W.maxlimbs))
}
// mask<0 selects the entire sample, otherwise selected binary-mask members.
// Caller shares W across a coalition listing; a limit is not silently reset.
real rowvector _g129_coalition(struct _T128_Scores scalar G,real scalar mask,real scalar tol,struct _C58_Work scalar W) {
 struct _C84_Q scalar q,a
 real scalar i,n,v,exact,missing
 n=rows(G.D);v=0;exact=1;missing=0;q=_g129_zero()
 for(i=1;i<=n;i++) if(mask<0 | mod(floor(mask/2^(i-1)),2)) {
  if(G.D[i,2]>=.) missing=1
  else v=v+G.D[i,2]
  if(G.Q[i,7]!=2) exact=0
 }
 if(missing) return((0,1001,.,.,0,0))
 if(!exact) return((1,1001,v,.,0,0))
 if(mask==0) return((2,0,0,0,0,0))
 for(i=1;i<=n;i++) if(mask<0 | mod(floor(mask/2^(i-1)),2)) {
  a=asarray(G.exact,i);a=_g129_norm(a,W);q=_g129_add(q,a,W)
  if(W.error) return(_g129_bad(1200+W.error,W))
 }
 return(_g129_out(q,tol))
}
// Pairwise reduction bounds depth to five for n<=32; source arrays are immutable.
struct _C84_Q scalar _g130_sum(transmorphic scalar source,real scalar count,struct _C58_Work scalar W) {
 transmorphic scalar current,next
 struct _C84_Q scalar a,b
 real scalar i,k,n
 n=count
 if(W.error) return(_g129_zero())
 if(n<0 | n>32 | n!=floor(n)) {W.error=4;return(_g129_zero());}
 if(n==0) return(_g129_zero())
 current=source
 while(n>1 & !W.error) {
  next=asarray_create("real",1);k=0
  for(i=1;i<=n;i=i+2) {
   a=asarray(current,i)
   if(i<n) {b=asarray(current,i+1);a=_g129_add(a,b,W);}
   k=k+1;asarray(next,k,a)
   if(W.error) break
  }
  current=next;n=k
 }
 if(W.error) return(_g129_zero())
 return(asarray(current,1))
}
void _g130_meta(struct _G129_Derived scalar B,real scalar zeros,real scalar caps) {
 B.meta=(B.A[.,3]',zeros,caps,(B.A[5,2]==1101),(!hasmissing(B.R) & !hasmissing(B.A[1..4,3])))
}
struct _G129_Derived scalar _g129_nash(struct _T128_Scores scalar G,real scalar tol,real matrix numerical,real rowvector oldmeta) {
 struct _G129_Derived scalar B
 struct _C58_Work scalar W,WP,failedwork
 struct _C84_Q scalar one,coef,floorq,cut,d,c,w,delta,total,wsum,unit,share,b,g,loss,allocated,unallocated,meanq,product,inc
 transmorphic scalar ds,increases,weights,slacks,gains,losses
 real rowvector out,logq
 real scalar n,i,j,allcert,failed,first,zeros,caps,logsum,logfailed,rowfailed,reason,capped
 n=rows(G.D);B.R=numerical;B.meta=oldmeta;B.Q=(B.R:<.);B.why=J(n,4,1001);B.Q[.,1]=G.Q[.,7];B.why[.,1]=G.Q[.,8]
 B.A=J(6,6,.);for(j=1;j<=6;j++) B.A[j,.]=(oldmeta[j]<.,1001,oldmeta[j],.,0,0)
 if(n<1 | n>32) return(B)
 ds=asarray_create("real",1);increases=asarray_create("real",1);weights=asarray_create("real",1)
 slacks=asarray_create("real",1);gains=asarray_create("real",1);losses=asarray_create("real",1)
 allcert=all(G.Q[.,7]:==2);failed=0;first=1
 floorq.N=(1,1);floorq.D=(1,100)
 // One bounded scenario/term context per row. Export failure does not erase
 // its exact expression, and partial certification retains known scenarios.
 for(i=1;i<=n;i++) if(G.Q[i,7]==2) {
  W=_c84_work();one=_c84_number(1,W);coef=_g129_div(one,_c84_number(5,W),W)
  cut=_g129_sub(one,_c84_number(tol,W),W);d=asarray(G.exact,i);d=_g129_norm(d,W)
  delta=_g129_sub(one,d,W)
  if(!_c84_gt(cut,d,W)) c=one
  else c=_g129_div(d,_g129_sub(one,_g129_mul(coef,delta,W),W),W)
  if(W.error | d.N[1]<0 | delta.N[1]<0) {
   reason=(W.error ? 1200+W.error : 1102);out=_g129_bad(reason,W)
   if(!failed) {failed=reason;failedwork=W;}
  }
  else {
   out=_g129_out(c,tol)
   if(allcert) {
    w=floorq;if(_c84_gt(delta,floorq,W)) w=delta
    inc=_g129_sub(c,d,W)
    if(W.error) {if(!failed) {failed=1200+W.error;failedwork=W;}}
    else {asarray(ds,i,d);asarray(slacks,i,delta);asarray(weights,i,w);asarray(increases,i,inc);}
   }
  }
  B.R[i,2]=out[3];B.Q[i,2]=out[1];B.why[i,2]=out[2]
 }
 if(!allcert) return(B)
 // Each aggregate has one bounded context. No failed context is retried.
 if(!failed) {
  W=_c84_work();total=_g130_sum(increases,n,W)
  if(W.error) {failed=1200+W.error;failedwork=W;}
  else {B.A[1,.]=_g129_out(total,tol);first=2;}
 }
 if(!failed) {
  W=_c84_work();wsum=_g130_sum(weights,n,W)
  if(W.error) {failed=1200+W.error;failedwork=W;}
 }
 if(!failed) {
  W=_c84_work();unit=_g129_div(total,wsum,W)
  if(W.error) {failed=1200+W.error;failedwork=W;}
 }
 if(failed) {
  B.R[.,3..4]=J(n,2,.);B.Q[.,3..4]=J(n,2,3);B.why[.,3..4]=J(n,2,failed)
  for(j=first;j<=6;j++) B.A[j,.]=_g129_bad(failed,failedwork)
  _g130_meta(B,.,.);return(B)
 }
 zeros=0;caps=0;rowfailed=0;logsum=0;logfailed=0
 // q_i=unit*w_i; g_i=min(q_i,1-d_i); loss_i=max(0,q_i-(1-d_i)).
 // Thus sum(g_i)=total-sum(loss_i), without repeatedly adding large gains.
 for(i=1;i<=n;i++) {
  W=_c84_work();d=asarray(ds,i);w=asarray(weights,i);delta=asarray(slacks,i)
  share=_g129_mul(unit,w,W);capped=_c84_gt(share,delta,W)
  if(capped) {b=one;g=delta;loss=_g129_sub(share,delta,W);}
  else {g=share;b=_g129_add(d,share,W);loss=_g129_zero();}
  if(W.error | g.N[1]<0 | loss.N[1]<0) {
   reason=(W.error ? 1200+W.error : 1102)
   if(!rowfailed) {rowfailed=reason;failedwork=W;}
   B.R[i,3..4]=(.,.);B.Q[i,3..4]=(3,3);B.why[i,3..4]=(reason,reason)
   continue
  }
  asarray(gains,i,g);asarray(losses,i,loss);caps=caps+capped;zeros=zeros+(g.N[1]==0)
  out=_g129_out(b,tol);B.R[i,3]=out[3];B.Q[i,3]=out[1];B.why[i,3]=out[2]
  out=_g129_out(g,tol);B.R[i,4]=out[3];B.Q[i,4]=out[1];B.why[i,4]=out[2]
  if(g.N[1]!=0) {logq=_g129_log(g);if(logq[1]==3) logfailed=logq[2];else logsum=logsum+logq[3];}
 }
 if(rowfailed) {
  for(j=2;j<=6;j++) B.A[j,.]=_g129_bad(rowfailed,failedwork)
  _g130_meta(B,.,.);return(B)
 }
 W=_c84_work();unallocated=_g130_sum(losses,n,W)
 if(W.error) for(j=2;j<=4;j++) B.A[j,.]=_g129_bad(1200+W.error,W)
 else {
  B.A[3,.]=_g129_out(unallocated,tol)
  W=_c84_work();allocated=_g129_sub(total,unallocated,W)
  if(W.error | allocated.N[1]<0) {
   reason=(W.error ? 1200+W.error : 1102);B.A[2,.]=B.A[4,.]=_g129_bad(reason,W)
  }
  else {
   B.A[2,.]=_g129_out(allocated,tol);meanq=_g129_div(allocated,_c84_number(n,W),W)
   B.A[4,.]=(W.error ? _g129_bad(1200+W.error,W) : _g129_out(meanq,tol))
  }
 }
 WP=_c84_work();product=_c84_number(1,WP)
 if(zeros) {product=_g129_zero();B.A[6,.]=(0,1004,.,.,0,0);}
 else {
  for(i=1;i<=n;i++) {g=asarray(gains,i);product=_g129_mul(product,g,WP);if(WP.error) break;}
  B.A[6,.]=((logfailed ? 3 : 1),(logfailed ? logfailed : 1005),(logfailed ? . : logsum),.,0,0)
 }
 B.A[5,.]=(WP.error ? _g129_bad(1200+WP.error,WP) : _g129_out(product,tol))
 _g130_meta(B,zeros,caps);return(B)
}

end
