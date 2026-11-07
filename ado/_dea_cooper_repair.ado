*! version 1.0.0 06nov2026
*! Internal bounded exact primal/dual basis exchange.
version 16.0
mata:
mata set matastrict on
// Exact one-row inverse update. The previous certificate supplies a verified
// integer inverse J/Db. Every new inverse is checked against the full new B.
// q=a*J, s=q[k], Db_new=abs(s),
// J_new=sign(s)*(s*J-J[,k]*q+Db*J[,k]*e_k')/Db.
struct _C44_Cert scalar _c47_update(real rowvector c,real matrix G,
 real colvector h,real scalar ne,struct _C44_Cert scalar prior,
 real scalar position,real scalar added,real scalar budget)
{
 struct _C44_Work scalar W
 struct _C44_Cert scalar out
 transmorphic scalar A,C,V,Q,N,q
 real rowvector ex,ce,a,b,v,pivot,den,pden,old_den
 real colvector ids
 real scalar n,i,j,k,sgn,ok,re,cre
 W.ops=0;W.limit=budget;W.limbs=128;W.maxlimbs=0;W.error=0
 out.status=-1;out.reason=1;out.solver_status=prior.solver_status;out.ops=0;out.maxlimbs=0
 out.phase=0;out.phase_cost=J(1,6,0);out.inverse_numerator=asarray_create("real",2)
 out.basis=out.face=out.primal_sign=out.dual_sign=J(0,1,.)
 out.numerator=out.dual_numerator=J(0,129,0)
 out.denominator=out.basis_denominator=out.row_exponent=J(1,0,.);out.rhs_exponent=.
 n=cols(c)
 if(n<1 | n>64 | rows(G)>256 | cols(G)!=n | rows(h)!=rows(G) | cols(h)!=1 |
  hasmissing(G) | hasmissing(h) | hasmissing(c) | missing(budget) | budget<0 |
  ne<0 | ne>rows(G) | (prior.status!=0 & prior.status!=1) |
  rows(prior.basis)!=n | cols(prior.row_exponent)!=n+1 |
  position<1 | position>n | position!=floor(position) |
  added<1 | added>rows(G) | added!=floor(added)) return(out)
 for(i=1;i<=n;i++) if(prior.basis[i]<1 | prior.basis[i]>rows(G) | prior.basis[i]!=floor(prior.basis[i])) return(out)
 if(cols(prior.basis_denominator)<2) return(out)
 old_den=prior.basis_denominator
 if(old_den[1]!=1) return(out)
 out.phase=1;A=_c46_rows(G,h,ex,re,W);C=_c46_rows(c,J(1,1,0),ce,cre,W);out.phase_cost[1]=W.ops
 if(!W.error) {
  out.phase=2
  if(re!=prior.rhs_exponent | any(ex[prior.basis]:!=prior.row_exponent[1..n])) W.error=5
  // The old inverse column links the new determinant to the verified old Db.
  for(i=1;i<=n;i++) {
   v=(0,0)
   for(j=1;j<=n;j++) v=_c44_add(v,_c44_mul(asarray(A,(prior.basis[i],j)),asarray(prior.inverse_numerator,(j,position)),W),W)
   v=_c44_add(v,_c44_neg(_c44_mul(old_den,_c44_int(i==position),W)),W)
   if(v[1]!=0 & !W.error) W.error=5
  }
  ids=prior.basis:+J(n,1,0);ids[position]=added;out.basis=ids;out.phase_cost[2]=W.ops
 }
 if(!W.error) {
  out.phase=3;q=asarray_create("real",2);V=asarray_create("real",2)
  for(j=1;j<=n;j++) {
   v=(0,0)
   for(i=1;i<=n;i++) v=_c44_add(v,_c44_mul(asarray(A,(added,i)),asarray(prior.inverse_numerator,(i,j)),W),W)
   asarray(q,(1,j),v)
  }
  pivot=asarray(q,(1,position));sgn=pivot[1]
  if(sgn==0 & !W.error) {out.reason=2;out.ops=W.ops;out.maxlimbs=W.maxlimbs;return(out);}
  den=_c44_mul(pivot,_c44_int(sgn),W)
  for(i=1;i<=n;i++) for(j=1;j<=n;j++) {
   if(j==position) v=_c44_mul(asarray(prior.inverse_numerator,(i,j)),_c44_int(sgn),W)
   else {
    a=_c44_mul(pivot,asarray(prior.inverse_numerator,(i,j)),W)
    b=_c44_mul(asarray(prior.inverse_numerator,(i,position)),asarray(q,(1,j)),W)
    v=_c44_div(_c44_add(a,_c44_neg(b),W),old_den,W)
    v=_c44_mul(v,_c44_int(sgn),W)
   }
   asarray(V,(i,j),v)
  }
  out.phase_cost[3]=W.ops
 }
 if(!W.error) {
  out.phase=4
  for(i=1;i<=n;i++) for(j=1;j<=n;j++) {
   v=(0,0)
   for(k=1;k<=n;k++) v=_c44_add(v,_c44_mul(asarray(A,(ids[i],k)),asarray(V,(k,j)),W),W)
   v=_c44_add(v,_c44_neg(_c44_mul(den,_c44_int(i==j),W)),W)
   if(v[1]!=0 & !W.error) W.error=5
  }
  out.phase_cost[4]=W.ops
 }
 if(!W.error) {
  out.phase=5;N=asarray_create("real",2);Q=asarray_create("real",2)
  pden=_c44_shift(den,-re,W)
  for(i=1;i<=n;i++) {
   v=(0,0)
   for(j=1;j<=n;j++) v=_c44_add(v,_c44_mul(asarray(V,(i,j)),asarray(A,(ids[j],n+1)),W),W)
   asarray(N,(i,1),v)
  }
  out.primal_sign=J(rows(G),1,.);out.dual_sign=J(n,1,.);ok=1
  for(i=1;i<=rows(G);i++) {
   v=_c44_neg(_c44_mul(asarray(A,(i,n+1)),den,W))
   for(j=1;j<=n;j++) v=_c44_add(v,_c44_mul(asarray(A,(i,j)),asarray(N,(j,1)),W),W)
   out.primal_sign[i]=v[1]
   if((i<=ne & v[1]!=0) | (i>ne & v[1]>0)) ok=0
  }
  out.phase_cost[5]=W.ops
 }
 if(!W.error) {
  out.phase=6
  for(i=1;i<=n;i++) {
   v=(0,0)
   for(j=1;j<=n;j++) v=_c44_add(v,_c44_mul(asarray(C,(1,j)),asarray(V,(j,i)),W),W)
   asarray(Q,(i,1),v);out.dual_sign[i]=v[1]
   if(ids[i]>ne & v[1]<0) ok=0
  }
  for(j=1;j<=n;j++) {
   v=_c44_neg(_c44_mul(asarray(C,(1,j)),den,W))
   for(i=1;i<=n;i++) v=_c44_add(v,_c44_mul(asarray(A,(ids[i],j)),asarray(Q,(i,1)),W),W)
   if(v[1]!=0 & !W.error) W.error=5
  }
  out.phase_cost[6]=W.ops
 }
 if(!W.error) {
  out.status=ok;out.reason=(ok ? 0 : 3);out.phase=0
  out.numerator=_c44_pack(N,1,n,128);out.dual_numerator=_c44_pack(Q,1,n,128)
  out.denominator=pden;out.basis_denominator=den;out.rhs_exponent=re;out.row_exponent=(ex[ids],ce)
  out.inverse_numerator=V
  if(ok) out.face=select(ids,(ids:>ne):&(out.dual_sign:>0))
 }
 else out.reason=10+W.error
 out.ops=W.ops;out.maxlimbs=W.maxlimbs
 return(out)
}

real rowvector _c45_rowcol(transmorphic scalar A,real scalar i,
 struct _C44_Cert scalar fit,real scalar k,real scalar n,struct _C44_Work scalar W)
{
 real scalar j
 real rowvector q
 q=(0,0)
 for(j=1;j<=n;j++) q=_c44_add(q,_c44_mul(asarray(A,(i,j)),asarray(fit.inverse_numerator,(j,k)),W),W)
 return(q)
}
real scalar _c45_ratio_compare(real rowvector a,real rowvector b,
 real rowvector c,real rowvector d,struct _C44_Work scalar W)
{
 real rowvector x,y
 x=_c44_mul(a,d,W);y=_c44_mul(c,b,W)
 return(_c44_cmp(x,y))
}
// Plan fields: mode(1 primal,2 dual),added row,removed row,position,stop reason.
// Only exact integer signs and cross-products determine the exchange.
real rowvector _c45_plan(real matrix G,real colvector h,real scalar ne,
 struct _C44_Cert scalar fit,struct _C44_Work scalar W)
{
 transmorphic scalar A
 real rowvector ex,q,a,b,bestnum,bestden,num,den
 real scalar n,i,j,k,pos,add,cmp,primal,dual,re
 n=cols(G);pos=0;add=0
 if(fit.status!=0) return((0,0,0,0,1))
 if(ne) if(any(fit.primal_sign[1..ne]:!=0)) return((0,0,0,0,3))
 primal=!any(fit.primal_sign:>0)
 dual=!any((fit.basis:>ne):&(fit.dual_sign:<0))
 if(!primal & !dual) return((0,0,0,0,2))
 A=_c46_rows(G,h,ex,re,W)
 if(!W.error & re!=fit.rhs_exponent) {W.error=5;return((0,0,0,0,8));}
 bestnum=(0,0);bestden=(1,1)
 if(primal & !dual) {
  for(k=1;k<=n;k++) if(fit.basis[k]>ne & fit.dual_sign[k]<0) {
   if(pos==0) pos=k
   else if(fit.basis[k]<fit.basis[pos]) pos=k
  }
  for(i=1;i<=ne;i++) {
   q=_c45_rowcol(A,i,fit,pos,n,W)
   if(q[1]!=0) return((0,0,0,0,4))
  }
  for(i=ne+1;i<=rows(G);i++) if(!any(fit.basis:==i)) {
   q=_c45_rowcol(A,i,fit,pos,n,W)
   if(q[1]<0) {
    den=_c44_neg(q);num=_c44_mul(asarray(A,(i,n+1)),fit.basis_denominator,W)
    for(j=1;j<=n;j++) num=_c44_add(num,_c44_neg(_c44_mul(asarray(A,(i,j)),_c44_norm(fit.numerator[j,.]),W)),W)
    if(num[1]<0 & !W.error) {W.error=5;return((0,0,0,0,8));}
    cmp=(add ? _c45_ratio_compare(num,den,bestnum,bestden,W) : -1)
    if(cmp<0 | (cmp==0 & i<add)) {add=i;bestnum=num;bestden=den;}
   }
  }
  if(W.error) return((0,0,0,0,8))
  if(!add) return((0,0,0,0,4))
  return((1,add,fit.basis[pos],pos,0))
 }
 if(dual & !primal) {
  for(i=ne+1;i<=rows(G);i++) if(fit.primal_sign[i]>0) {add=i;break;}
  for(k=1;k<=n;k++) if(fit.basis[k]>ne) {
   den=_c45_rowcol(A,add,fit,k,n,W)
   if(den[1]>0) {
    num=_c44_norm(fit.dual_numerator[k,.])
    cmp=(pos ? _c45_ratio_compare(num,den,bestnum,bestden,W) : -1)
    if(pos==0) {pos=k;bestnum=num;bestden=den;}
    else if(cmp<0 | (cmp==0 & fit.basis[k]<fit.basis[pos])) {pos=k;bestnum=num;bestden=den;}
   }
  }
  if(W.error) return((0,0,0,0,8))
  if(!pos) return((0,0,0,0,4))
  return((2,add,fit.basis[pos],pos,0))
 }
 return((0,0,0,0,1))
}
void _c45_save(string scalar tag,struct _C44_Cert scalar fit)
{
 savebin(tag+"_meta.bin",(fit.status,fit.reason,fit.solver_status,fit.ops,fit.maxlimbs))
 savebin(tag+"_basis.bin",fit.basis);savebin(tag+"_face.bin",fit.face)
 savebin(tag+"_num.bin",fit.numerator);savebin(tag+"_den.bin",fit.denominator)
 savebin(tag+"_basis_den.bin",fit.basis_denominator);savebin(tag+"_rhs_exp.bin",fit.rhs_exponent)
 savebin(tag+"_dual.bin",fit.dual_numerator);savebin(tag+"_exponent.bin",fit.row_exponent)
 savebin(tag+"_primal_sign.bin",fit.primal_sign);savebin(tag+"_dual_sign.bin",fit.dual_sign)
 savebin(tag+"_phase.bin",fit.phase);savebin(tag+"_phase_cost.bin",fit.phase_cost)
}
void _c45_repair(string scalar tag,real rowvector c,real matrix G,real colvector h,
 real scalar ne,struct _C44_Cert scalar fit,real scalar maxpivots,real scalar totalbudget)
{
 struct _C44_Work scalar W
 real matrix seen,trace
 real colvector ids
 real rowvector plan,key
 real scalar step,stop,total,j,cycle,maxlimbs,solverstatus
 trace=J(0,11,.);seen=J(0,cols(c),.);step=0;stop=0;total=fit.ops;maxlimbs=fit.maxlimbs
 _c45_save(tag+"_repair0",fit)
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
  _c45_save(tag+"_repair"+strofreal(step),fit)
  trace=trace\(step,plan[1..4],W.ops,fit.ops,fit.status,fit.reason,fit.phase,total)
 }
 if(fit.status==-1 & !stop) stop=9
 savebin(tag+"_repair_trace.bin",trace);savebin(tag+"_repair_info.bin",(stop,step,total,maxlimbs))
}
end
