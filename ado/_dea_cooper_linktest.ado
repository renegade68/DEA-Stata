*! version 1.0.0 06nov2026
*! Internal original-model faces and proposal support for Cooper.
version 16.0
mata:
mata set matastrict on
// Additional bounded proposal for nonnegative variables. Keep original RHS;
// positive lower bounds become explicit inequalities instead of a rounded shift.
struct _C44_Cert scalar _c44_unshifted(real rowvector c,real matrix E,real colvector b,
 real matrix U,real colvector ub,real rowvector lo,real rowvector hi,
 real scalar tol,real scalar budget,real matrix G,real colvector h)
{
 struct RSM_Result scalar fit
 struct _C44_Cert scalar out
 real matrix As
 real colvector bs,choices,upper,ineqmap
 real scalar n,j,ne,nu
 n=cols(c);ne=rows(E);nu=rows(U);G=E\U\(-I(n));h=b\ub\(-lo')
 upper=vec(selectindex(hi':<.))
 for(j=1;j<=rows(upper);j++) {G=G\I(n)[upper[j],.];h=h\hi[upper[j]];}
 choices=J(0,1,.)
 if(any(lo:<0) | hasmissing(lo)) return(_c44_certify(c,G,h,ne,choices,budget))
 As=U;bs=ub;ineqmap=(nu ? (ne+1::ne+nu) : J(0,1,.))
 for(j=1;j<=n;j++) if(lo[j]>0) {
  As=As\(-I(n)[j,.]);bs=bs\(-lo[j]);ineqmap=ineqmap\(ne+nu+j)
 }
 for(j=1;j<=rows(upper);j++) {
  As=As\I(n)[upper[j],.];bs=bs\hi[upper[j]];ineqmap=ineqmap\(ne+nu+n+j)
 }
 fit=_rsm_solve_lp(As,bs,E,b,-c,tol,0)
 if(fit.status==0 & rows(fit.basiscols)>0) {
  for(j=1;j<=n;j++) if(lo[j]==0 & !any(fit.basiscols:==j)) choices=choices\(ne+nu+j)
  for(j=1;j<=rows(ineqmap);j++) if(!any(fit.basiscols:==n+j)) choices=choices\ineqmap[j]
  if(ne) choices=choices\(1::ne)
 }
 out=_c44_certify(c,G,h,ne,choices,budget);out.solver_status=fit.status
 return(out)
}
void _c44_emit(string scalar tag,real rowvector c,real matrix E,real colvector b,
 real matrix U,real colvector ub,real rowvector lo,real rowvector hi,
 struct _C44_Cert scalar fit)
{
 savebin(tag+"_c.bin",c);savebin(tag+"_E.bin",E);savebin(tag+"_b.bin",b)
 savebin(tag+"_U.bin",U);savebin(tag+"_ub.bin",ub);savebin(tag+"_lo.bin",lo);savebin(tag+"_hi.bin",hi)
 savebin(tag+"_meta.bin",(fit.status,fit.reason,fit.solver_status,fit.ops,fit.maxlimbs))
 savebin(tag+"_basis.bin",fit.basis);savebin(tag+"_face.bin",fit.face)
 savebin(tag+"_num.bin",fit.numerator);savebin(tag+"_den.bin",fit.denominator)
 savebin(tag+"_basis_den.bin",fit.basis_denominator);savebin(tag+"_rhs_exp.bin",fit.rhs_exponent)
 savebin(tag+"_dual.bin",fit.dual_numerator);savebin(tag+"_exponent.bin",fit.row_exponent)
 savebin(tag+"_primal_sign.bin",fit.primal_sign);savebin(tag+"_dual_sign.bin",fit.dual_sign)
 savebin(tag+"_phase.bin",fit.phase);savebin(tag+"_phase_cost.bin",fit.phase_cost)
}
// Preserve the latest fully checked, nonsingular basis when a later numerical
// proposal is unavailable. Status 0 is only a repair seed, never an optimum.
// Existing successful proposals and a final status-0 proposal retain priority.
// No new tolerance, numerical solve, oracle basis, or relaxed acceptance is used.
void _c48_propose(string scalar tag,real rowvector c,real matrix E,real colvector b,
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
  _c44_emit(tag+"_attempt"+strofreal(attempt),c,E,b,U,ub,lo,hi,fit)
  selected=attempt;ops=ops+fit.ops;limbs=max((limbs,fit.maxlimbs))
  if(fit.status==1) break
  if(fit.status==0) {retained=fit;retainedG=G;retainedh=h;last=attempt;}
 }
 if(fit.status==-1 & last>0) {
  fit=retained;G=retainedG;h=retainedh;selected=last;restored=1
 }
 // Attempt work is separate from the existing per-repair 80-million budget.
 // Repair's starting certificate cost is already included in proposal ops.
 savebin(tag+"_proposal.bin",(min((attempt,3)),selected,restored,ops,limbs))
}
real scalar _c44_step(string scalar tag,real rowvector c,real matrix E,real colvector b,
 real matrix U,real colvector ub,real rowvector lo,real rowvector hi,
 struct _C44_Cert scalar fit,real matrix G,real colvector h)
{
 _c48_propose(tag,c,E,b,U,ub,lo,hi,fit,G,h)
 if(fit.status==0) _c45_repair(tag,c,G,h,rows(E),fit,16,80000000)
 _c44_emit(tag,c,E,b,U,ub,lo,hi,fit)
 return(fit.status==1)
}
void _c44_face(real matrix E,real colvector b,real matrix G,real colvector h,
 struct _C44_Cert scalar fit,real matrix F,real colvector fb)
{
 F=E;fb=b
 if(rows(fit.face)) {F=F\G[fit.face,.];fb=fb\h[fit.face];}
}
void _c44_pipeline(string scalar source,real scalar n,real scalar m,real scalar s,
 real rowvector previous_lambda)
{
 struct _C44_Cert scalar fit
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
