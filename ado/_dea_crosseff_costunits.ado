*! version 1.0.0 08nov2026
*! Bounded failure-only secondary cost recovery and explicit objective units.
version 16.0
capture mata: mata drop _C107_Result()
capture mata: mata drop _C107_Cost()
capture mata: mata drop _c108_needed()
capture mata: mata drop _c108_join()
capture mata: mata drop _c107_ready()
capture mata: mata drop _c107_valid()
capture mata: mata drop _c107_gcd()
capture mata: mata drop _c107_div()
capture mata: mata drop _c107_cost()
capture mata: mata drop _c107_attempt()
mata:
mata set matastrict on
// A certificate's dual numerator is in working-objective units. Recover the
// original dual by multiplying by positive factor; primal/ray units do not change.
// Keep this context outside _C58_Cert: existing primary certificates keep their ABI.
struct _C107_Cost {
 transmorphic scalar original,working
 real rowvector factor
 real scalar exponent,n,valid,error,ops,maxlimbs
}
struct _C107_Result {
 struct _C107_Cost scalar cost
 struct _C58_Cert scalar fit
 real scalar attempted,complete,ops,maxlimbs,original_value
 real scalar prior_ops,prior_pivots,new_pivots,total_pivots
}
real scalar _c107_ready() return(10702)
// Canonical signed base-2^26 inputs only. Refuse malformed or over-width inputs
// before arithmetic; do not silently round fractional digits or truncate limbs.
real scalar _c107_valid(real rowvector a,struct _C58_Work scalar W)
{
 real scalar n
 if(W.error) return(0)
 n=cols(a)
 if(n<2) {W.error=4;return(0);}
 if(!_c58b_tick(W,4*n,n-1)) return(0)
 if(hasmissing(a) | any(a:!=floor(a)) | abs(a[1])>1 |
    any(a[2..n]:<0) | any(a[2..n]:>=67108864)) {W.error=4;return(0);}
 if((a[1]==0 & (n!=2 | a[2]!=0)) | (a[1]!=0 & a[n]==0)) {W.error=4;return(0);}
 return(1)
}
// Binary gcd: subtraction and exact power-of-two removal only. Every loop
// consumes the caller's existing work budget; no unbounded integer conversion.
real rowvector _c107_gcd(real rowvector a0,real rowvector b0,struct _C58_Work scalar W)
{
 real rowvector a,b,t
 real scalar ea,eb,e,order
 if(!_c107_valid(a0,W) | !_c107_valid(b0,W)) return((0,0))
 a=a0;b=b0
 if(a[1]) a[1]=1
 if(b[1]) b[1]=1
 if(a[1]==0) return(b)
 if(b[1]==0) return(a)
 a=_c58f_split(a,ea,W);b=_c58f_split(b,eb,W);e=min((ea,eb))
 while(b[1]!=0 & !W.error) {
  if(!_c58b_tick(W,cols(a)+cols(b),max((cols(a),cols(b)))-1)) break
  order=_c58b_cmp(a,b)
  if(order>0) {t=a;a=b;b=t;}
  b=_c58b_add(b,_c58b_neg(a),W)
  b=_c58f_split(b,eb,W)
 }
 if(W.error) return((0,0))
 return(_c58f_restore(a,e,W))
}
// Exact signed division by aligned binary subtraction. Alignment never exceeds
// the dividend. This also handles a 128-limb divisor without a 129-limb trial
// product. A nonzero remainder is an error, never a rounded quotient.
real rowvector _c107_div(real rowvector a0,real rowvector b0,struct _C58_Work scalar W)
{
 real rowvector a,b,p,q,bit
 real scalar s,k,order
 if(!_c107_valid(a0,W) | !_c107_valid(b0,W)) return((0,0))
 if(b0[1]==0) {W.error=3;return((0,0));}
 if(a0[1]==0) return((0,0))
 a=a0;b=b0;s=a[1]*b[1];a[1]=b[1]=1;q=(0,0)
 while(a[1]!=0 & !W.error) {
  if(!_c58b_tick(W,cols(a)+cols(b)+52,max((cols(a),cols(b)))-1)) break
  order=_c58b_cmp(a,b)
  if(order<0) {W.error=3;break;}
  k=_c58d_bits(a)-_c58d_bits(b);p=_c58f_restore(b,k,W)
  if(W.error) break
  if(!_c58b_tick(W,cols(a)+cols(p),max((cols(a),cols(p)))-1)) break
  if(_c58b_cmp(p,a)>0) {k=k-1;p=_c58f_restore(b,k,W);}
  a=_c58b_add(a,_c58b_neg(p),W)
  bit=_c58f_restore((1,1),k,W);q=_c58b_add(q,bit,W)
 }
 if(W.error) return((0,0))
 q[1]=s;return(q)
}
struct _C107_Cost scalar _c107_cost(transmorphic scalar C,real scalar n,
 real scalar ce,struct _C58_Work scalar W)
{
 struct _C107_Cost scalar out
 real rowvector a,g,q
 real scalar j
 out.original=asarray_create("real",2);out.working=asarray_create("real",2)
 out.factor=(0,0);out.exponent=ce;out.n=n;out.valid=0
 out.error=W.error;out.ops=W.ops;out.maxlimbs=W.maxlimbs
 if(W.error) return(out)
 if(missing(n) | n<1 | n>64 | n!=floor(n) | missing(ce) | ce!=floor(ce)) W.error=4
 // Copy each integer explicitly: the original and working arrays never alias C.
 if(!W.error) for(j=1;j<=n;j++) {
  if(!asarray_contains(C,(1,j))) {W.error=4;break;}
  a=asarray(C,(1,j));if(!_c107_valid(a,W)) break
  asarray(out.original,(1,j),a)
 }
 g=(0,0)
 if(!W.error) for(j=1;j<=n;j++) {
  a=asarray(out.original,(1,j));g=_c107_gcd(g,a,W)
  if(W.error | (cols(g)==2 & g[2]==1)) break
 }
 if(!W.error) {
  if(g[1]==0) g=(1,1) // The all-zero objective uses the positive unit factor.
  for(j=1;j<=n;j++) {
   a=asarray(out.original,(1,j));q=_c107_div(a,g,W)
   if(W.error) break
   // Exact division with zero remainder proves a=g*q; no over-width product.
   asarray(out.working,(1,j),q)
  }
 }
 if(!W.error) {out.factor=g;out.valid=1;}
 else out.working=asarray_create("real",2)
 out.error=W.error;out.ops=W.ops;out.maxlimbs=W.maxlimbs
 return(out)
}
// Research entry point for a failed secondary certificate. The numerical c is
// retained as proposal only; all exact secondary checks share cost.working.
// This does not replace an accepted fit, change the public command, or reset
// work/pivot limits. The caller must retain its old result unless complete==1.
struct _C107_Result scalar _c107_attempt(real rowvector c,real matrix G,
 real colvector h,real scalar ne,real colvector choices,transmorphic scalar C,
 real scalar ce,real scalar spent,real scalar used_pivots,real scalar wide_decode)
{
 struct _C107_Result scalar out
 struct _C58_Work scalar W
 real rowvector report,journal
 real scalar used
 out.prior_ops=spent;out.prior_pivots=used_pivots;out.new_pivots=0;out.total_pivots=used_pivots
 out.attempted=out.complete=0;out.ops=spent;out.maxlimbs=0;out.original_value=.
 // Initialize all legacy certificate fields even when factor discovery fails.
 out.fit=_c70b_certify(c,G,h,ne,J(0,1,.),0,C,ce,wide_decode)
 W.ops=0;W.limit=0;W.limbs=128;W.maxlimbs=0;W.error=4
 out.cost=_c107_cost(C,cols(c),ce,W)
 if(missing(spent) | spent<0 | spent!=floor(spent) | spent>=80000000 |
    missing(used_pivots) | used_pivots<0 | used_pivots!=floor(used_pivots) | used_pivots>=16 |
    cols(c)<1 | cols(c)>16 | cols(G)!=cols(c) | rows(G)>256 |
    rows(h)!=rows(G) | cols(h)!=1 | hasmissing(c) | hasmissing(G) | hasmissing(h) |
    missing(ne) | ne<0 | ne>rows(G) | ne!=floor(ne) |
    (wide_decode!=0 & wide_decode!=1)) return(out)
 W.error=0;W.limit=min((10000000,80000000-spent));out.attempted=1
 out.cost=_c107_cost(C,cols(c),ce,W)
 out.ops=spent+W.ops;out.maxlimbs=W.maxlimbs
 if(!out.cost.valid) return(out)
 out.fit=_c70b_certify(c,G,h,ne,choices,max((0,min((10000000-W.ops,80000000-out.ops)))),out.cost.working,ce,wide_decode)
 out.fit.ops=out.fit.ops+out.ops;out.fit.maxlimbs=max((out.fit.maxlimbs,W.maxlimbs))
 if(out.fit.status==0) _c70r_repair("",c,G,h,ne,out.fit,16-used_pivots,80000000,out.cost.working,ce,wide_decode)
 used=used_pivots+(out.fit.repair_steps<. ? out.fit.repair_steps : 0)
 journal=_c62_initial(out.fit)
 if(out.fit.status!=1) {
  _c70_recover(c,G,h,ne,out.fit,report,journal,used,out.cost.working,ce,wide_decode)
  out.fit.recovery=report
  // Legacy journals count this attempt. Keep inherited pivots in the context.
  journal[13]=journal[13]-used_pivots
 }
 out.new_pivots=journal[13];out.total_pivots=used_pivots+out.new_pivots
 out.fit.trace62=journal;out.ops=out.fit.ops;out.maxlimbs=out.fit.maxlimbs
 out.complete=(out.fit.status==1 | journal[2]==2)
 if(out.fit.status==1) out.original_value=_c70_value(out.cost.original,ce,out.fit)
 return(out)
}
// Only a failed secondary certificate with a recorded integer-width refusal is
// eligible. Successful certificates and certified rays are never retried.
real scalar _c108_needed(struct _C58_Cert scalar prior)
{
 real rowvector t
 t=prior.trace62
 if(prior.status!=-1 | cols(t)!=36 | prior.ops>=80000000 | t[13]>=16 | t[16]==1) return(0)
 return(prior.reason==12 | t[27]==12 | t[33]==12)
}
// Preserve original attempt accounting in the fixed public trace schema.
// The caller MUST retain retry.cost alongside the returned normalized dual.
// Existing raw consumers use only primal values, signs and diagnostics.
struct _C58_Cert scalar _c108_join(struct _C58_Cert scalar prior,
 struct _C107_Result scalar retry)
{
 struct _C58_Cert scalar out
 real rowvector a,b,z
 real scalar k
 out=prior
 if(!retry.attempted) return(out)
 a=prior.trace62;b=retry.fit.trace62
 if(retry.prior_ops!=prior.ops | retry.prior_pivots!=a[13] |
    retry.ops<prior.ops | retry.total_pivots>16) return(out)
 if(retry.complete) out=retry.fit
 // A failed factor search has no certificate journal. Record its charged
 // work without inventing a feasible basis or discarding the prior failure.
 if(!retry.cost.valid) {
  z=a;z[1]=1;z[2]=3;z[3]=5;z[4]=2
  z[21]=z[21]+retry.ops-prior.ops;z[22]=retry.ops
 }
 else {
  z=b;z[1]=1;z[2]=(retry.complete ? (b[16]==1 ? 2 : 1) : 3)
  z[3]=5;z[4]=(retry.complete ? 0 : 2);z[5..6]=a[5..6]
  z[7]=a[7]
  for(k=8;k<=12;k++) z[k]=((k==8 | k==10) ? max((a[k],b[k])) : a[k]+b[k])
  // The factor attempt's first check/repair is recovery of the original
  // secondary objective, not part of the original pre-dispatcher work.
  z[12]=z[12]+b[7]
  z[13]=a[13]+b[13];z[25]=a[25]+b[25]
  z[17]=a[17]
  for(k=18;k<=21;k++) z[k]=a[k]+b[k]
  z[21]=z[21]+b[17]-prior.ops
  z[22]=retry.ops
  // Retain old stage observations when the new attempt did not reach them.
  for(k=26;k<=36;k++) if(z[k]>=. & a[k]<.) z[k]=a[k]
 }
 out.ops=retry.ops;out.maxlimbs=max((prior.maxlimbs,retry.maxlimbs))
 out.solver_status=prior.solver_status;z[23]=out.maxlimbs
 out.trace62=z
 return(out)
}

end
