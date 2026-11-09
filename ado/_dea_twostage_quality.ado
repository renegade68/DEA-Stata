*! version 1.0.0 09nov2026
*! Internal exact two-stage support; use the public DEA commands.
version 16.0
capture mata: mata drop _t123_ready()
capture mata: mata drop _t123_stageq()
capture mata: mata drop _t123_aggregate()
mata:
mata set matastrict on
real scalar _t123_ready() {
 string rowvector names
 real scalar j
 names=("_t123_aggregate()","_t123_stageq()","_t120_rel()","_t121_nc()","_t122_stage()","_dea_ts_compute()","_dea_ts_stage()","_dea_ts_joint()","_dea_ts_link()","_c84_work()","_c84_number()","_c84_sub()","_c84_mul()","_c84_div()","_c84_value()","_c84_gt()","_t120_eval()","_t120_relerr()")
 for(j=1;j<=cols(names);j++) if(findexternal(names[j])==NULL) return(0)
 return(12301)
}
struct _C84_Q scalar _t123_stageq(struct _T122_Result scalar R,string scalar ort,struct _C58_Work scalar W) {
 struct _C84_Q scalar q
 q=_t120_eval(R.c*R.T,R.primary,W)
 if(ort=="input") q.N=_c58b_neg(q.N)
 else q=_c84_div(_c84_number(1,W),q,W)
 return(q)
}
// quality, reason, exported value, relative error, operations, limbs,
// arithmetic error, legacy rounded-stage aggregate. No epsilon substitution.
real rowvector _t123_aggregate(struct _T122_Result scalar A,struct _T122_Result scalar B,
 string scalar mode,string scalar ort,real scalar tol,real scalar legacy) {
 struct _C58_Work scalar W
 struct _C84_Q scalar a,b,q,zero,one
 real rowvector out
 real scalar good,value,reason
 out=(legacy<.,1001,legacy,.,0,0,0,legacy)
 if(!A.ok | !B.ok) return(out)
 out[1]=3;out[3]=.;reason=0;good=1
 if((mode!="product" & mode!="mean" & mode!="independent") | (ort!="input" & ort!="output") | tol<=0 | tol>=.) {out[2]=1002;return(out);}
 W=_c84_work();zero=_c84_number(0,W);one=_c84_number(1,W)
 a=_t123_stageq(A,ort,W);b=_t123_stageq(B,ort,W)
 if(mode=="mean") q=_c84_div(_c84_sub(a,_c84_sub(zero,b,W),W),_c84_number(2,W),W)
 else q=_c84_mul(a,b,W)
 if(!W.error) {
  value=_c84_value(q)
  if(q.N[1]<=0 | _c84_gt(q,one,W) | value<0 | value>=.) reason=1102
  else if(value==0) reason=1101
  else {
   out[4]=_t120_relerr(_c84_number(value,W),q,tol,W,good)
   if(!good | out[4]>=.) reason=1103
  }
 }
 out[5..7]=(W.ops,W.maxlimbs,W.error)
 if(W.error) reason=1200+W.error
 out[2]=reason
 if(reason==0) {out[1]=2;out[3]=value;}
 return(out)
}
end
