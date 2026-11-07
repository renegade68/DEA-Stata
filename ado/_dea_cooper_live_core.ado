*! version 1.0.0 06nov2026
*! Internal Cooper LP/certificate/export support loaded as a matched bundle.
version 16.0
quietly findfile _dea_cooper_live.ado
local _c52_dir = substr("`r(fn)'",1,length("`r(fn)'")-length("_dea_cooper_live.ado"))
tempname _c52_f
file open `_c52_f' using "`_c52_dir'_dea_lpcheck.ado", read text
forvalues _c52_i=1/1 {
 file read `_c52_f' _c52_line
}
file close `_c52_f'
if `"`_c52_line'"' != "*! version 1.0.0  01nov2026" {
 di as error "Cooper dependency version mismatch: _dea_lpcheck"
 exit 499
}
tempname _c52_f
file open `_c52_f' using "`_c52_dir'_dea_nonte_lp.ado", read text
forvalues _c52_i=1/1 {
 file read `_c52_f' _c52_line
}
file close `_c52_f'
if `"`_c52_line'"' != "*! version 1.0.0  02nov2026" {
 di as error "Cooper dependency version mismatch: _dea_nonte_lp"
 exit 499
}
tempname _c52_f
file open `_c52_f' using "`_c52_dir'_dea_teaudit.ado", read text
forvalues _c52_i=1/1 {
 file read `_c52_f' _c52_line
}
file close `_c52_f'
if `"`_c52_line'"' != "*! version 1.0.0  02nov2026" {
 di as error "Cooper dependency version mismatch: _dea_teaudit"
 exit 499
}
tempname _c52_f
file open `_c52_f' using "`_c52_dir'_dea_recovery.ado", read text
forvalues _c52_i=1/1 {
 file read `_c52_f' _c52_line
}
file close `_c52_f'
if `"`_c52_line'"' != "*! version 0.1.0 03nov2026" {
 di as error "Cooper dependency version mismatch: _dea_recovery"
 exit 499
}
tempname _c52_f
file open `_c52_f' using "`_c52_dir'_dea_simplex.ado", read text
forvalues _c52_i=1/6 {
 file read `_c52_f' _c52_line
}
file close `_c52_f'
if `"`_c52_line'"' != "*! Version 2.0.25 - 2026-11-05" {
 di as error "Cooper dependency version mismatch: _dea_simplex"
 exit 499
}
tempname _c52_f
file open `_c52_f' using "`_c52_dir'_dea_cooper_cert.ado", read text
forvalues _c52_i=1/1 {
 file read `_c52_f' _c52_line
}
file close `_c52_f'
if `"`_c52_line'"' != "*! version 1.0.0 06nov2026" {
 di as error "Cooper dependency version mismatch: _dea_cooper_cert"
 exit 499
}
tempname _c52_f
file open `_c52_f' using "`_c52_dir'_dea_cooper_repair.ado", read text
forvalues _c52_i=1/1 {
 file read `_c52_f' _c52_line
}
file close `_c52_f'
if `"`_c52_line'"' != "*! version 1.0.0 06nov2026" {
 di as error "Cooper dependency version mismatch: _dea_cooper_repair"
 exit 499
}
tempname _c52_f
file open `_c52_f' using "`_c52_dir'_dea_cooper_linktest.ado", read text
forvalues _c52_i=1/1 {
 file read `_c52_f' _c52_line
}
file close `_c52_f'
if `"`_c52_line'"' != "*! version 1.0.0 06nov2026" {
 di as error "Cooper dependency version mismatch: _dea_cooper_linktest"
 exit 499
}
tempname _c52_f
file open `_c52_f' using "`_c52_dir'_dea_cooper_stage2.ado", read text
forvalues _c52_i=1/1 {
 file read `_c52_f' _c52_line
}
file close `_c52_f'
if `"`_c52_line'"' != "*! version 1.0.0 06nov2026" {
 di as error "Cooper dependency version mismatch: _dea_cooper_stage2"
 exit 499
}
tempname _c52_f
file open `_c52_f' using "`_c52_dir'_dea_cooper_result.ado", read text
forvalues _c52_i=1/1 {
 file read `_c52_f' _c52_line
}
file close `_c52_f'
if `"`_c52_line'"' != "*! version 1.0.0 06nov2026" {
 di as error "Cooper dependency version mismatch: _dea_cooper_result"
 exit 499
}
capture mata: mata drop _c52_ready()
capture mata: mata drop _c52_bridge()
capture mata: mata drop _c52_pipeline()
capture mata: mata drop _c52_endpoint()
capture mata: mata drop _c52_objective()
capture mata: mata drop _c52_point()
capture mata: mata drop _c52_parent()
capture mata: mata drop _c52_empty()
capture mata: mata drop _c52_load()
capture mata: mata drop _c52_select()
capture mata: mata drop _c52_lift()
capture mata: mata drop _c52_step()
capture mata: mata drop _c52_propose()
capture mata: mata drop _c52_emit()
capture mata: mata drop _c52_repair()
capture mata: mata drop _c52_log_cert()
capture mata: mata drop _c52_save()
capture mata: mata drop _c51_bridge()
capture mata: mata drop _c51_result()
capture mata: mata drop _c51_values()
capture mata: mata drop _c51_associated()
capture mata: mata drop _c51_model()
capture mata: mata drop _c51_same()
capture mata: mata drop _c51_has()
capture mata: mata drop _c51_load()
capture mata: mata drop _c51_scaled()
capture mata: mata drop _c51_decode_export()
capture mata: mata drop _c51_take()
capture mata: mata drop _c51_double()
capture mata: mata drop _c51_bits()
capture mata: mata drop _c51_valid()
capture mata: mata drop _c50_pipeline()
capture mata: mata drop _c50_select()
capture mata: mata drop _c50_lift()
capture mata: mata drop _c44_pipeline()
capture mata: mata drop _c44_face()
capture mata: mata drop _c44_step()
capture mata: mata drop _c48_propose()
capture mata: mata drop _c44_emit()
capture mata: mata drop _c44_unshifted()
capture mata: mata drop _c45_repair()
capture mata: mata drop _c45_save()
capture mata: mata drop _c45_plan()
capture mata: mata drop _c45_ratio_compare()
capture mata: mata drop _c45_rowcol()
capture mata: mata drop _c47_update()
capture mata: mata drop _dea_cooper_cert_ready()
capture mata: mata drop _c44_solve()
capture mata: mata drop _c44_certify()
capture mata: mata drop _c44_pack()
capture mata: mata drop _c44_pick()
capture mata: mata drop _c46_rows()
capture mata: mata drop _c44_shift()
capture mata: mata drop _c44_decode()
capture mata: mata drop _c44_div()
capture mata: mata drop _c44_mul()
capture mata: mata drop _c44_add()
capture mata: mata drop _c44_neg()
capture mata: mata drop _c44_cmp()
capture mata: mata drop _c44_int()
capture mata: mata drop _c44_tick()
capture mata: mata drop _c44_norm()
capture mata: mata drop _C52_Live()
capture mata: mata drop _C51_Result()
capture mata: mata drop _C44_Cert()
capture mata: mata drop _C44_Work()
quietly run "`_c52_dir'_dea_lpcheck.ado"
quietly run "`_c52_dir'_dea_nonte_lp.ado"
quietly run "`_c52_dir'_dea_teaudit.ado"
quietly run "`_c52_dir'_dea_recovery.ado"
quietly run "`_c52_dir'_dea_simplex.ado"
quietly run "`_c52_dir'_dea_cooper_cert.ado"
quietly run "`_c52_dir'_dea_cooper_repair.ado"
quietly run "`_c52_dir'_dea_cooper_linktest.ado"
quietly run "`_c52_dir'_dea_cooper_stage2.ado"
quietly run "`_c52_dir'_dea_cooper_result.ado"

mata:
mata set matastrict on
void _c52_save(string scalar file,real matrix A)
{
 string scalar directory
 real matrix B
 real scalar f
 directory=st_local("audit")
 if(directory=="") return
 file=directory+"/"+file
 if(fileexists(file)) _error(602)
 B=bufio();bufbyteorder(B,2);f=fopen(file,"w")
 fbufput(B,f,"%4b",rows(A));fbufput(B,f,"%4b",cols(A));fbufput(B,f,"%8z",A);fclose(f)
}
void _c52_log_cert(string scalar tag,struct _C44_Cert scalar fit)
{
 _c52_save(tag+"_meta.bin",(fit.status,fit.reason,fit.solver_status,fit.ops,fit.maxlimbs))
 _c52_save(tag+"_basis.bin",fit.basis);_c52_save(tag+"_face.bin",fit.face)
 _c52_save(tag+"_num.bin",fit.numerator);_c52_save(tag+"_den.bin",fit.denominator)
 _c52_save(tag+"_basis_den.bin",fit.basis_denominator);_c52_save(tag+"_rhs_exp.bin",fit.rhs_exponent)
 _c52_save(tag+"_dual.bin",fit.dual_numerator);_c52_save(tag+"_exponent.bin",fit.row_exponent)
 _c52_save(tag+"_primal_sign.bin",fit.primal_sign);_c52_save(tag+"_dual_sign.bin",fit.dual_sign)
 _c52_save(tag+"_phase.bin",fit.phase);_c52_save(tag+"_phase_cost.bin",fit.phase_cost)
}
void _c52_repair(string scalar tag,real rowvector c,real matrix G,real colvector h,
 real scalar ne,struct _C44_Cert scalar fit,real scalar maxpivots,real scalar totalbudget)
{
 struct _C44_Work scalar W
 real matrix seen,trace
 real colvector ids
 real rowvector plan,key
 real scalar step,stop,total,j,cycle,maxlimbs,solverstatus
 trace=J(0,11,.);seen=J(0,cols(c),.);step=0;stop=0;total=fit.ops;maxlimbs=fit.maxlimbs
 _c52_log_cert(tag+"_repair0",fit)
 if(rows(fit.basis)==cols(c)) seen=sort(fit.basis,1)'
 while(fit.status==0) {
  if(step>=maxpivots) {stop=6;break;}
  if(total>=totalbudget) {stop=7;break;}
  W.ops=0;W.limit=min((10000000,totalbudget-total));W.limbs=128;W.maxlimbs=0;W.error=0
  plan=_c45_plan(G,h,ne,fit,W);total=total+W.ops;maxlimbs=max((maxlimbs,W.maxlimbs))
  if(W.error) {stop=8;break;}
  if(plan[5]) {stop=plan[5];break;}
  if(total>=totalbudget) {stop=7;break;}
  ids=fit.basis;ids[plan[4]]=plan[2];key=sort(ids,1)';cycle=0
  for(j=1;j<=rows(seen);j++) if(all(seen[j,.]:==key)) cycle=1
  if(cycle) {stop=5;break;}
  seen=seen\key;solverstatus=fit.solver_status
  fit=_c47_update(c,G,h,ne,fit,plan[4],plan[2],min((10000000,totalbudget-total)));fit.solver_status=solverstatus
  total=total+fit.ops;maxlimbs=max((maxlimbs,fit.maxlimbs));step=step+1
  _c52_log_cert(tag+"_repair"+strofreal(step),fit)
  trace=trace\(step,plan[1..4],W.ops,fit.ops,fit.status,fit.reason,fit.phase,total)
 }
 if(fit.status==-1 & !stop) stop=9
 _c52_save(tag+"_repair_trace.bin",trace);_c52_save(tag+"_repair_info.bin",(stop,step,total,maxlimbs))
}
void _c52_emit(string scalar tag,real rowvector c,real matrix E,real colvector b,
 real matrix U,real colvector ub,real rowvector lo,real rowvector hi,
 struct _C44_Cert scalar fit)
{
 _c52_save(tag+"_c.bin",c);_c52_save(tag+"_E.bin",E);_c52_save(tag+"_b.bin",b)
 _c52_save(tag+"_U.bin",U);_c52_save(tag+"_ub.bin",ub);_c52_save(tag+"_lo.bin",lo);_c52_save(tag+"_hi.bin",hi)
 _c52_save(tag+"_meta.bin",(fit.status,fit.reason,fit.solver_status,fit.ops,fit.maxlimbs))
 _c52_save(tag+"_basis.bin",fit.basis);_c52_save(tag+"_face.bin",fit.face)
 _c52_save(tag+"_num.bin",fit.numerator);_c52_save(tag+"_den.bin",fit.denominator)
 _c52_save(tag+"_basis_den.bin",fit.basis_denominator);_c52_save(tag+"_rhs_exp.bin",fit.rhs_exponent)
 _c52_save(tag+"_dual.bin",fit.dual_numerator);_c52_save(tag+"_exponent.bin",fit.row_exponent)
 _c52_save(tag+"_primal_sign.bin",fit.primal_sign);_c52_save(tag+"_dual_sign.bin",fit.dual_sign)
 _c52_save(tag+"_phase.bin",fit.phase);_c52_save(tag+"_phase_cost.bin",fit.phase_cost)
}
// Preserve the latest fully checked, nonsingular basis when a later numerical
// proposal is unavailable. Status 0 is only a repair seed, never an optimum.
// Existing successful proposals and a final status-0 proposal retain priority.
// No new tolerance, numerical solve, oracle basis, or relaxed acceptance is used.
void _c52_propose(string scalar tag,real rowvector c,real matrix E,real colvector b,
 real matrix U,real colvector ub,real rowvector lo,real rowvector hi,
 struct _C44_Cert scalar fit,real matrix G,real colvector h)
{
 struct _C44_Cert scalar retained
 real matrix retainedG
 real colvector retainedh
 real scalar attempt,tol,last,selected,restored,ops,limbs
 last=0;selected=0;restored=0;ops=0;limbs=0
 for(attempt=1;attempt<=3;attempt++) {
  tol=(attempt==1 ? 1e-12 : 1e-9)
  if(attempt<3) fit=_c44_solve(c,E,b,U,ub,lo,hi,tol,10000000,G,h)
  else fit=_c44_unshifted(c,E,b,U,ub,lo,hi,1e-12,10000000,G,h)
  _c52_emit(tag+"_attempt"+strofreal(attempt),c,E,b,U,ub,lo,hi,fit)
  selected=attempt;ops=ops+fit.ops;limbs=max((limbs,fit.maxlimbs))
  if(fit.status==1) break
  if(fit.status==0) {retained=fit;retainedG=G;retainedh=h;last=attempt;}
 }
 if(fit.status==-1 & last>0) {
  fit=retained;G=retainedG;h=retainedh;selected=last;restored=1
 }
 // Attempt work is separate from the existing per-repair 80-million budget.
 // Repair's starting certificate cost is already included in proposal ops.
 _c52_save(tag+"_proposal.bin",(min((attempt,3)),selected,restored,ops,limbs))
}
real scalar _c52_step(string scalar tag,real rowvector c,real matrix E,real colvector b,
 real matrix U,real colvector ub,real rowvector lo,real rowvector hi,
 struct _C44_Cert scalar fit,real matrix G,real colvector h)
{
 _c52_propose(tag,c,E,b,U,ub,lo,hi,fit,G,h)
 if(fit.status==0) _c52_repair(tag,c,G,h,rows(E),fit,16,80000000)
 _c52_emit(tag,c,E,b,U,ub,lo,hi,fit)
 return(fit.status==1)
}

struct _C44_Cert scalar _c52_lift(string scalar tag,real rowvector c,
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
 _c52_save(tag+"_lift_plan.bin",(bestrow,W.ops,W.maxlimbs,W.error))
 if(W.error) ids=J(0,1,.)
 out=_c44_certify(c,SG,sh,newne,ids,10000000)
 out.solver_status=.
 return(out)
}

void _c52_select(string scalar tag,real rowvector c,real matrix E,real colvector b,
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
  (maxpivots!=16 & maxpivots!=128)) {_c52_save(tag+"_state.bin",(-1,0));return;}
 if(hasmissing(reference)) {_c52_save(tag+"_state.bin",(2,0));return;}
 _c44_face(E,b,G,h,base,F,fb)
 Q=(J(n,v-n,0),I(n));A=(F,J(rows(F),1,0))
 V=(U,J(rows(U),1,0))\(Q,J(n,1,-1))\(-Q,J(n,1,-1))
 vb=ub\reference'\(-reference');objective=(J(1,v,0),-1)
 lower=(lo,0);upper=(hi,.)
 _c52_propose(tag,objective,A,fb,V,vb,lower,upper,chosen,SG,sh)
 if(chosen.status==-1 & rows(chosen.basis)==0) chosen=_c52_lift(tag,objective,SG,sh,rows(A),rows(E),rows(U),n,reference,base)
 if(chosen.status==0) _c52_repair(tag,objective,SG,sh,rows(A),chosen,maxpivots,80000000)
 _c52_emit(tag,objective,A,fb,V,vb,lower,upper,chosen)
 _c52_save(tag+"_state.bin",(chosen.status==1,maxpivots))
}

// Original-LP/reference file framing is checked before any allocation or read.
// This reader is separate from the preserved T51 checkpoint reader.
real matrix _c52_load(string scalar file)
{
 real matrix B,A
 real scalar f,n,m,bytes
 f=_fopen(file,"r");if(f<0) return(J(0,0,.))
 fseek(f,0,1);bytes=ftell(f);fseek(f,0,-1)
 if(bytes<8 | bytes>16000008) {fclose(f);return(J(0,0,.));}
 B=bufio();bufbyteorder(B,2);n=fbufget(B,f,"%4b");m=fbufget(B,f,"%4b")
 if(n<0 | m<0 | n*m>2000000 | bytes!=8+8*n*m) {fclose(f);return(J(0,0,.));}
 A=(n*m ? fbufget(B,f,"%8z",n,m) : J(n,m,.));fclose(f);return(A)
}
struct _C52_Live {
 struct _C51_Result scalar point
 real rowvector pmin,pmax,cmin,cmax,state
 real scalar complete,conversion,rounded
 string scalar range_reason
}
struct _C51_Result scalar _c52_empty(real scalar n,real scalar m)
{
 struct _C51_Result scalar R
 R.lambda=R.mu=J(1,n,.);R.scores=J(1,5,.);R.total=R.technical=R.congestion=J(1,m,.)
 R.usable=R.applied=R.conversion=R.rounded=R.base_cert=R.selected_cert=0
 R.cap=R.distance=.;R.policy_status=4;R.reason="base_certificate_unavailable";R.source="unavailable"
 return(R)
}
// Check the selected point against the original model with integers, before
// stripping the selector coordinate. No rounded intermediate RHS is used.
real scalar _c52_parent(real matrix G,real colvector h,real scalar ne,
 real matrix N,real rowvector D)
{
 struct _C44_Work scalar W
 transmorphic scalar A
 real rowvector ex,z,t
 real scalar re,i,j
 W.ops=0;W.limit=10000000;W.limbs=512;W.maxlimbs=0;W.error=0
 A=_c46_rows(G,h,ex,re,W)
 for(i=1;i<=rows(G);i++) {
  z=(0,0)
  for(j=1;j<=cols(G);j++) z=_c44_add(z,_c44_mul(asarray(A,(i,j)),_c44_norm(N[j,.]),W),W)
  z=_c44_shift(z,-re,W);t=_c44_mul(asarray(A,(i,cols(G)+1)),D,W)
  z=_c44_add(z,_c44_neg(t),W)
  if((i<=ne & z[1]!=0) | (i>ne & z[1]>0)) return(0)
 }
 return(!W.error)
}
struct _C51_Result scalar _c52_point(real rowvector c,real matrix G,real colvector h,
 real scalar ne,struct _C44_Cert scalar base,struct _C44_Cert scalar chosen,
 real matrix X,real scalar index,real scalar n,real scalar m,real scalar s,
 real rowvector units,string scalar policy,real rowvector reference,real scalar cap)
{
 struct _C51_Result scalar R
 real matrix N
 real rowvector D
 real scalar v
 R=_c52_empty(n,m);v=cols(c)
 if(base.status!=1) return(R)
 R.base_cert=1;R.source="baseline";R.policy_status=0;R.reason="not_requested"
 N=base.numerator;D=base.denominator
 if(policy=="stage2") {
  R.policy_status=2;R.reason="reference_missing"
  if(!hasmissing(reference)) {
   R.policy_status=3;R.reason="selection_uncertified";R.cap=cap
   if(chosen.status==1 & rows(chosen.numerator)==v+1) {
    R.policy_status=4;R.reason="selection_parent_mismatch"
    if(_c51_associated(c,n,base,chosen)) if(_c52_parent(G,h,ne,chosen.numerator[1..v,.],chosen.denominator)) {
     R.policy_status=R.applied=R.selected_cert=1;R.source="stage2";R.reason="certified"
     N=chosen.numerator[1..v,.];D=chosen.denominator
     R.distance=_c51_take(_c44_norm(chosen.numerator[v+1,.]),D,R)
    }
   }
  }
 }
 _c51_values(N,D,X,index,n,m,s,units,R)
 if(R.conversion) {
  R.lambda=R.mu=J(1,n,.);R.scores=J(1,5,.);R.total=R.technical=R.congestion=J(1,m,.)
  R.distance=.;R.applied=0;R.rounded=0;R.source="unavailable";R.reason="conversion_unavailable";return(R)
 }
 R.usable=1;return(R)
}
void _c52_objective(real rowvector c,struct _C44_Cert scalar fit,
 real rowvector N,real rowvector D,struct _C44_Work scalar W)
{
 transmorphic scalar C
 real rowvector ex
 real scalar re,j
 C=_c46_rows(c,J(1,1,0),ex,re,W);N=(0,0);D=fit.denominator
 for(j=1;j<=cols(c);j++) N=_c44_add(N,_c44_mul(asarray(C,(1,j)),_c44_norm(fit.numerator[j,.]),W),W)
 N=_c44_shift(N,max((0,ex[1])),W);D=_c44_shift(D,max((0,-ex[1])),W)
}
real scalar _c52_endpoint(real rowvector c,struct _C44_Cert scalar fit,
 real scalar signum,real scalar congestion,real scalar xi,real scalar scale,
 struct _C51_Result scalar R)
{
 struct _C44_Work scalar W
 real rowvector N,D,a
 real scalar e,value
 if(fit.status!=1) return(.)
 W.ops=0;W.limit=10000000;W.limbs=512;W.maxlimbs=0;W.error=0
 _c52_objective(c,fit,N,D,W)
 if(signum<0) N=_c44_neg(N)
 if(congestion) {
  a=_c44_decode(xi,e,W)
  a=_c44_mul(_c44_shift(a,max((0,e)),W),D,W)
  a=_c44_add(a,_c44_neg(_c44_shift(N,max((0,-e)),W)),W)
  D=_c44_shift(D,max((0,-e)),W);N=a
 }
 value=_c51_scaled(N,D,scale,W,R)
 if(W.error) R.conversion=2
 return(value)
}
struct _C52_Live scalar _c52_pipeline(string scalar source,real scalar n,real scalar m,real scalar s,
 real scalar index,real rowvector previous_lambda,real rowvector previous_mu,
 string scalar policy,real scalar maxpivots,real rowvector xu,real rowvector yu,real scalar ranges)
{
 struct _C52_Live scalar L
 struct _C51_Result scalar RR
 struct _C44_Cert scalar fit,base,chosen
 real matrix E,U,G,Fr,Fs,Fa,Ua,B,xs,ys,Ed,Ud,Fd,Q,BG
 real colvector b,ub,h,br,bs,ba,uba,hB,bd,ubd,bf,Bh
 real rowvector c,lo,hi,cs,ct,la,ha,ca,cd,ld,hd,units
 real scalar v,v1,j,signum,count,ok,value
 string scalar tag
 L.point=_c52_empty(n,m);L.state=J(1,6,0);L.pmin=L.pmax=J(1,m+s,.);L.cmin=L.cmax=J(1,m,.)
 L.complete=L.conversion=L.rounded=0;L.range_reason=(ranges ? "pipeline_unavailable" : "not_requested")
 c=_c52_load(source+"/radial/c.bin");E=_c52_load(source+"/radial/E.bin");b=_c52_load(source+"/radial/b.bin")
 U=_c52_load(source+"/radial/U.bin");ub=_c52_load(source+"/radial/ub.bin")
 lo=_c52_load(source+"/radial/lo.bin");hi=_c52_load(source+"/radial/hi.bin")
 cs=_c52_load(source+"/slack/c.bin");ct=_c52_load(source+"/decomposition/c.bin")
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
void _c52_bridge(string scalar source,real scalar n,real scalar m,real scalar s,real scalar index,
 string scalar policy,string scalar reference1,string scalar reference2,real rowvector xu,real rowvector yu,
 real scalar ranges,real scalar cap,string rowvector names)
{
 struct _C52_Live scalar L
 real rowvector l1,l2
 real matrix A
 l1=l2=J(1,n,.)
 if(reference1!="") {A=_c52_load(reference1);if(cols(A)!=n | rows(A)<index) _error(198);l1=A[index,.];}
 if(reference2!="") {A=_c52_load(reference2);if(cols(A)!=n | rows(A)<index) _error(198);l2=A[index,.];}
 L=_c52_pipeline(source,n,m,s,index,l1,l2,policy,cap,xu,yu,ranges)
 st_matrix(names[1],L.point.lambda);st_matrix(names[2],L.point.mu);st_matrix(names[3],L.point.scores)
 st_matrix(names[4],L.point.total);st_matrix(names[5],L.point.technical);st_matrix(names[6],L.point.congestion)
 st_matrix(names[7],L.pmin);st_matrix(names[8],L.pmax);st_matrix(names[9],L.cmin);st_matrix(names[10],L.cmax)
 st_matrix(names[11],(L.point.usable,L.point.policy_status,L.point.applied,L.point.conversion,L.point.rounded,L.point.base_cert,L.point.selected_cert,L.point.cap,L.point.distance,L.complete,L.conversion,L.rounded))
 st_matrix(names[12],L.state)
 st_local("_c52_reason",L.point.reason);st_local("_c52_source",L.point.source);st_local("_c52_range_reason",L.range_reason)
}

real scalar _c52_ready() return(52)
end
