*! version 1.0.0 06nov2026
*! Exact intended peer objective; bounded certificate path; old primary path unchanged.
version 16.0
capture mata: mata drop _c70_exact_ready()
capture mata: mata drop _c70_decode()
capture mata: mata drop _c70_build()
capture mata: mata drop _c70_value()
capture mata: mata drop _c70_fit()
capture mata: mata drop _c70b_certify()
capture mata: mata drop _c70b_solve()
capture mata: mata drop _c70u_update()
capture mata: mata drop _c70r_repair()
capture mata: mata drop _c70_recover()
capture mata: mata drop _c70_ray()
mata:
mata set matastrict on
real scalar _c70_exact_ready() return(7001)
// Decode any finite Stata binary64 operand by exact binary steps. The old
// coefficient decoder and its supported domain remain unchanged.
real rowvector _c70_decode(real scalar x,real scalar exponent,struct _C58_Work scalar W)
{
 real scalar m,k,i,reconstructed
 exponent=0
 if(x==0) return((0,0))
 if(missing(x)) {W.error=4;return((0,0));}
 m=abs(x);k=0
 while(m>=9007199254740992) {m=m/2;exponent=exponent+1;k=k+1;}
 while(m!=floor(m)) {m=m*2;exponent=exponent-1;k=k+1;}
 while(mod(m,2)==0) {m=m/2;exponent=exponent+1;k=k+1;}
 reconstructed=m
 for(i=1;i<=abs(exponent);i++) reconstructed=reconstructed*(exponent>0 ? 2 : .5)
 if(reconstructed!=abs(x) | !_c58b_tick(W,k+abs(exponent)+1,3)) {if(!W.error) W.error=4;return((0,0));}
 return(_c58b_int(sign(x)*m))
}
// C is the EXACT maximizing objective -c_min after T. All columns share ce.
// The numerical minimizing vector is a proposal; it is never substituted
// into exact stationarity, repair, or terminal recession checks.
void _c70_build(real matrix X,real matrix Y,real scalar index,
 real scalar output,real scalar aggressive,real matrix T,
 real rowvector proposal,transmorphic scalar C,real scalar ce,real rowvector od)
{
 struct _C58_Work scalar W
 transmorphic scalar terms
 real matrix A,exponents
 real rowvector a,b,total,check
 real scalar n,m,s,d,i,j,k,tt,e,te,minimum,maximum,bits,ops,value,part,t,bb,error,rounded,lost,sgn
 n=rows(X);m=cols(X);s=cols(Y);d=cols(T)
 W.ops=0;W.limit=10000000;W.limbs=128;W.maxlimbs=0;W.error=0
 C=asarray_create("real",2);terms=asarray_create("real",2)
 proposal=J(1,d,.);od=(1,-1,4,3,0,0);ce=0
 if(n>32 | d>16 | index<1 | index>n | rows(T)<m+s | hasmissing(X) | hasmissing(Y) | hasmissing(T)) {od[2]=-2;return;}
 A=(X,Y);exponents=J(n,d,0);minimum=.
 for(j=1;j<=d;j++) {
  k=0
  for(tt=1;tt<=rows(T);tt++) if(T[tt,j]!=0) {if(k) {W.error=4;break;}k=tt;}
  for(i=1;i<=n;i++) asarray(terms,(i,j),(0,0))
  if(W.error) break
  if(!k | k>m+s | (output & k>m)) continue
  b=_c70_decode(T[k,j],te,W)
  if(cols(b)!=2 | abs(b[2])!=1) {W.error=4;break;}
  sgn= -1*(k<=m ? 1 : -1)*(aggressive ? -1 : 1)*b[1]
  for(i=1;i<=n;i++) if(i!=index & A[i,k]!=0) {
   a=_c70_decode(A[i,k],e,W);a[1]=a[1]*sgn;asarray(terms,(i,j),a)
   exponents[i,j]=e+te;minimum=min((minimum,e+te))
   if(W.error) break
  }
  if(W.error) break
 }
 if(minimum>=.) minimum=0
 ce=minimum;maximum=0
 if(!W.error) for(j=1;j<=d;j++) {
  total=(0,0)
  for(i=1;i<=n;i++) {
   a=asarray(terms,(i,j))
   if(a[1]) total=_c58b_add(total,_c58b_shift(a,exponents[i,j]-minimum,W),W)
  }
  total=_c58b_norm(total);asarray(C,(1,j),total)
  if(total[1]) maximum=max((maximum,_c58d_bits(total)))
 }
 if(maximum>0) maximum=maximum-1
 rounded=lost=0
 if(!W.error) {
  proposal=J(1,d,0)
  for(j=1;j<=d;j++) {
   a=asarray(C,(1,j));value=0
   for(k=2;k<=cols(a);k++) if(a[k]) {
    ops=W.ops;part=_c70_scale(a[k],26*(k-2)-maximum,ops,W.limit);W.ops=ops
    if(part>=.) {W.error=1;break;}
    t=value+part;value=t
    if(!_c58b_tick(W,1,cols(a)-1)) break
   }
   proposal[j]= -a[1]*value
   if(a[1]) {
    if(value==0) lost=1
    else {
     b=_c70_decode(-proposal[j],e,W)
     if(e+maximum>=0) b=_c58b_shift(b,e+maximum,W)
     else a=_c58b_shift(a,-e-maximum,W)
     check=_c58b_add(a,_c58b_neg(b),W)
     if(check[1]) rounded=1
    }
   }
   if(W.error) break
  }
 }
 if(W.error | hasmissing(proposal)) {proposal=J(1,d,.);od[5..6]=(W.ops,W.maxlimbs);return;}
 od=(1,1,(lost ? 3 : (rounded ? 2 : 1)),3,W.ops,W.maxlimbs)
}
real scalar _c70_value(transmorphic scalar C,real scalar ce,struct _C58_Cert scalar fit)
{
 struct _C58_Work scalar W
 real rowvector N,D,z
 real scalar j
 if(fit.status!=1) return(.)
 W.ops=0;W.limit=10000000;W.limbs=512;W.maxlimbs=0;W.error=0
 N=(0,0);D=fit.denominator
 for(j=1;j<=rows(fit.numerator);j++) N=_c58b_add(N,_c58b_mul(asarray(C,(1,j)),_c58b_norm(fit.numerator[j,.]),W),W)
 N=_c58b_shift(_c58b_neg(N),max((0,ce)),W);D=_c58b_shift(D,max((0,-ce)),W)
 z=_c58d_double(N,D)
 return(W.error | z[2] ? . : z[1])
}
void _c70_fit(real rowvector c,real matrix E,real colvector b,
 real matrix U,real colvector ub,real scalar tol,
 transmorphic scalar exactC,real scalar exact_exponent,real rowvector build,
 struct _C58_Cert scalar fit,real matrix G,real colvector h)
{
 real rowvector recovery,journal
 real scalar used
 fit=_c70b_solve(-c,E,b,U,ub,J(1,cols(c),0),J(1,cols(c),.),tol,max((0,10000000-build[5])),G,h,exactC,exact_exponent)
 fit.ops=fit.ops+build[5];fit.maxlimbs=max((fit.maxlimbs,build[6]))
 fit.phase_cost[1]=fit.phase_cost[1]+build[5]
 if(fit.status==0) _c70r_repair("",-c,G,h,rows(E),fit,16,80000000,exactC,exact_exponent)
 journal=_c62_initial(fit)
 used=(fit.repair_steps<. ? fit.repair_steps : 0)
 if(fit.status!=1) {_c70_recover(-c,G,h,rows(E),fit,recovery,journal,used,exactC,exact_exponent);fit.recovery=recovery;}
 fit.trace62=journal
}

struct _C58_Cert scalar _c70b_certify(real rowvector c,real matrix G,
 real colvector h,real scalar ne,real colvector choices,real scalar budget,transmorphic scalar C0,real scalar ce0)
{
 struct _C58_Work scalar W
 struct _C58_Cert scalar out
 transmorphic scalar A,T,C,N,Q
 real rowvector ex,ce,a,b,v,tmp,pivot,previous,den,pden,cn,qn
 real colvector ids
 real scalar n,i,j,k,p,sgn,ok,re,cre
 W.ops=0;W.limit=budget;W.limbs=128;W.maxlimbs=0;W.error=0
 out.repair_stop=out.repair_steps=.
 out.recovery=J(1,15,0)
 out.trace62=J(1,36,0);out.ray=out.ray_num=J(0,129,0);out.ray_den=J(1,0,.)
 out.status=-1;out.reason=1;out.solver_status=.;out.ops=0;out.maxlimbs=0
 out.phase=0;out.phase_cost=J(1,6,0);out.inverse_numerator=asarray_create("real",2)
 out.basis=out.face=out.primal_sign=out.dual_sign=J(0,1,.)
 out.numerator=out.dual_numerator=J(0,129,0);out.denominator=out.basis_denominator=J(1,0,.);out.rhs_exponent=.;out.row_exponent=J(1,0,.)
 n=cols(c)
 if(n<1 | n>64 | rows(G)>256 | cols(G)!=n | rows(h)!=rows(G) | cols(h)!=1 |
  hasmissing(G) | hasmissing(h) | hasmissing(c) | rows(choices)<n |
  missing(budget) | budget<0 | ne<0 | ne>rows(G)) return(out)
 for(i=1;i<=rows(choices);i++) if(choices[i]<1 | choices[i]>rows(G) | choices[i]!=floor(choices[i])) return(out)
 out.phase=1;A=_c58i_rows(G,h,ex,re,W);C=C0;ce=ce0;out.phase_cost[1]=W.ops
 if(!W.error) {
  out.phase=2
  ids=(rows(choices)==n ? choices : _c58b_pick(A,choices,n,W))
  out.phase_cost[2]=W.ops
  if(rows(ids)!=n) {out.reason=(W.error ? 10+W.error : 2);out.ops=W.ops;out.maxlimbs=W.maxlimbs;return(out);}
  out.basis=ids;T=asarray_create("real",2);previous=(1,1)
  for(i=1;i<=n;i++) {
   for(j=1;j<=n;j++) asarray(T,(i,j),asarray(A,(ids[i],j)))
   for(j=1;j<=n;j++) asarray(T,(i,n+j),_c58b_int(i==j))
  }
  // Fraction-free Gauss-Jordan. Every division must have exact zero remainder.
  out.phase=3
  for(k=1;k<=n;k++) {
   p=0
   for(i=k;i<=n;i++) {a=asarray(T,(i,k));if(a[1]!=0) {p=i;break;}}
   if(p==0) {out.reason=2;out.ops=W.ops;out.maxlimbs=W.maxlimbs;return(out);}
   if(p!=k) for(j=1;j<=2*n;j++) {tmp=asarray(T,(k,j));asarray(T,(k,j),asarray(T,(p,j)));asarray(T,(p,j),tmp);}
   pivot=asarray(T,(k,k))
   for(i=1;i<=n;i++) if(i!=k) {
    a=asarray(T,(i,k))
    for(j=1;j<=2*n;j++) if(j!=k) {
     b=_c58b_mul(pivot,asarray(T,(i,j)),W);v=_c58b_mul(a,asarray(T,(k,j)),W)
     v=_c58b_div(_c58b_add(b,_c58b_neg(v),W),previous,W);asarray(T,(i,j),v)
    }
    asarray(T,(i,k),(0,0))
   }
   previous=pivot
   if(W.error) break
  }
  out.phase_cost[3]=W.ops
  if(!W.error) {
   den=asarray(T,(n,n));sgn=den[1];den[1]=1;N=asarray_create("real",2);Q=asarray_create("real",2)
   // Exact inverse identity proves rank and uniqueness; not a floating rank test.
   out.phase=4
   for(i=1;i<=n;i++) for(j=1;j<=n;j++) {
    v=(0,0)
    for(k=1;k<=n;k++) v=_c58b_add(v,_c58b_mul(asarray(A,(ids[i],k)),asarray(T,(k,n+j)),W),W)
    b=_c58b_int(i==j);a=den;a[1]=sgn
    v=_c58b_add(v,_c58b_neg(_c58b_mul(a,b,W)),W)
    if(v[1]!=0 & !W.error) W.error=5
   }
   out.phase_cost[4]=W.ops;if(!W.error) out.phase=5
   // Solve only after the coefficient inverse is certified. Primal and dual
   // denominators differ: Dx=Db*2^-re, N=Adj*H; Q=cint*Adj uses Db.
   pden=_c58b_shift(den,-re,W)
   for(i=1;i<=n;i++) {
    a=(0,0)
    for(j=1;j<=n;j++) a=_c58b_add(a,_c58b_mul(asarray(T,(i,n+j)),asarray(A,(ids[j],n+1)),W),W)
    a[1]=a[1]*sgn;asarray(N,(i,1),a)
   }
   out.primal_sign=J(rows(G),1,.);out.dual_sign=J(n,1,.);ok=1
   for(i=1;i<=rows(G);i++) {
    v=_c58b_neg(_c58b_mul(asarray(A,(i,n+1)),den,W))
    for(j=1;j<=n;j++) v=_c58b_add(v,_c58b_mul(asarray(A,(i,j)),asarray(N,(j,1)),W),W)
    out.primal_sign[i]=v[1]
    if((i<=ne & v[1]!=0) | (i>ne & v[1]>0)) ok=0
   }
   out.phase_cost[5]=W.ops;if(!W.error) out.phase=6
   for(i=1;i<=n;i++) {
    qn=(0,0)
    for(j=1;j<=n;j++) qn=_c58b_add(qn,_c58b_mul(asarray(C,(1,j)),asarray(T,(j,n+i)),W),W)
    qn[1]=qn[1]*sgn;asarray(Q,(i,1),qn);out.dual_sign[i]=qn[1]
    if(ids[i]>ne & qn[1]<0) ok=0
   }
   // Check stationarity again directly in integer row units.
   for(j=1;j<=n;j++) {
    v=_c58b_neg(_c58b_mul(asarray(C,(1,j)),den,W))
    for(i=1;i<=n;i++) v=_c58b_add(v,_c58b_mul(asarray(A,(ids[i],j)),asarray(Q,(i,1)),W),W)
    if(v[1]!=0 & !W.error) W.error=5
   }
   out.phase_cost[6]=W.ops
   if(!W.error) {
    out.status=ok;out.reason=(ok ? 0 : 3)
    out.numerator=_c58b_pack(N,1,n,128);out.dual_numerator=_c58b_pack(Q,1,n,128)
    out.denominator=pden;out.basis_denominator=den;out.rhs_exponent=re;out.row_exponent=(ex[ids],ce)
    for(i=1;i<=n;i++) for(j=1;j<=n;j++) {
     a=asarray(T,(i,n+j));a[1]=a[1]*sgn;asarray(out.inverse_numerator,(i,j),a)
    }
    out.phase=0
    if(ok) out.face=select(ids,(ids:>ne):&(out.dual_sign:>0))
   }
  }
 }
 if(W.error) {out.status=-1;out.reason=10+W.error;}
 out.ops=W.ops;out.maxlimbs=W.maxlimbs
 return(out)
}

struct _C58_Cert scalar _c70b_solve(real rowvector c,real matrix E,real colvector b,
 real matrix U,real colvector ub,real rowvector lo,real rowvector hi,
 real scalar tol,real scalar budget,real matrix G,real colvector h,transmorphic scalar C0,real scalar ce0)
{
 struct RSM_Result scalar fit
 struct _C58_Cert scalar out
 real matrix As
 real colvector bs,beqs,choices,upper,ineqmap
 real scalar n,j,ne,nu
 n=cols(c);ne=rows(E);nu=rows(U)
 G=E\U\(-I(n));h=b\ub\(-lo')
 upper=vec(selectindex(hi':<.));ineqmap=(nu ? (ne+1::ne+nu) : J(0,1,.))
 for(j=1;j<=rows(upper);j++) {
  G=G\I(n)[upper[j],.];h=h\hi[upper[j]];ineqmap=ineqmap\rows(G)
 }
 As=U;bs=ub-U*lo';beqs=b-E*lo'
 for(j=1;j<=rows(upper);j++) {As=As\I(n)[upper[j],.];bs=bs\(hi[upper[j]]-lo[upper[j]]);}
 fit=_rsm_solve_lp(As,bs,E,beqs,-c,tol,0)
 choices=J(0,1,.)
 if(fit.status==0 & rows(fit.basiscols)>0) {
  for(j=1;j<=n;j++) if(!any(fit.basiscols:==j)) choices=choices\(ne+nu+j)
  for(j=1;j<=rows(ineqmap);j++) if(!any(fit.basiscols:==n+j)) choices=choices\ineqmap[j]
  if(ne) choices=choices\(1::ne)
 }
 out=_c70b_certify(c,G,h,ne,choices,budget,C0,ce0);out.solver_status=fit.status
 return(out)
}

struct _C58_Cert scalar _c70u_update(real rowvector c,real matrix G,
 real colvector h,real scalar ne,struct _C58_Cert scalar prior,
 real scalar position,real scalar added,real scalar budget,transmorphic scalar C0,real scalar ce0)
{
 struct _C58_Work scalar W
 struct _C58_Cert scalar out
 transmorphic scalar A,C,V,Q,N,q
 real rowvector ex,ce,a,b,v,pivot,den,pden,old_den
 real colvector ids
 real scalar n,i,j,k,sgn,ok,re,cre
 W.ops=0;W.limit=budget;W.limbs=128;W.maxlimbs=0;W.error=0
 out.repair_stop=out.repair_steps=.
 out.recovery=J(1,15,0)
 out.trace62=J(1,36,0);out.ray=out.ray_num=J(0,129,0);out.ray_den=J(1,0,.)
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
 out.phase=1;A=_c58i_rows(G,h,ex,re,W);C=C0;ce=ce0;out.phase_cost[1]=W.ops
 if(!W.error) {
  out.phase=2
  if(re!=prior.rhs_exponent | any(ex[prior.basis]:!=prior.row_exponent[1..n])) W.error=5
  // The old inverse column links the new determinant to the verified old Db.
  for(i=1;i<=n;i++) {
   v=(0,0)
   for(j=1;j<=n;j++) v=_c58b_add(v,_c58b_mul(asarray(A,(prior.basis[i],j)),asarray(prior.inverse_numerator,(j,position)),W),W)
   v=_c58b_add(v,_c58b_neg(_c58b_mul(old_den,_c58b_int(i==position),W)),W)
   if(v[1]!=0 & !W.error) W.error=5
  }
  ids=prior.basis:+J(n,1,0);ids[position]=added;out.basis=ids;out.phase_cost[2]=W.ops
 }
 if(!W.error) {
  out.phase=3;q=asarray_create("real",2);V=asarray_create("real",2)
  for(j=1;j<=n;j++) {
   v=(0,0)
   for(i=1;i<=n;i++) v=_c58b_add(v,_c58b_mul(asarray(A,(added,i)),asarray(prior.inverse_numerator,(i,j)),W),W)
   asarray(q,(1,j),v)
  }
  pivot=asarray(q,(1,position));sgn=pivot[1]
  if(sgn==0 & !W.error) {out.reason=2;out.ops=W.ops;out.maxlimbs=W.maxlimbs;return(out);}
  den=_c58b_mul(pivot,_c58b_int(sgn),W)
  for(i=1;i<=n;i++) for(j=1;j<=n;j++) {
   if(j==position) v=_c58b_mul(asarray(prior.inverse_numerator,(i,j)),_c58b_int(sgn),W)
   else {
    a=_c58b_mul(pivot,asarray(prior.inverse_numerator,(i,j)),W)
    b=_c58b_mul(asarray(prior.inverse_numerator,(i,position)),asarray(q,(1,j)),W)
    v=_c58b_div(_c58b_add(a,_c58b_neg(b),W),old_den,W)
    v=_c58b_mul(v,_c58b_int(sgn),W)
   }
   asarray(V,(i,j),v)
  }
  out.phase_cost[3]=W.ops
 }
 if(!W.error) {
  out.phase=4
  for(i=1;i<=n;i++) for(j=1;j<=n;j++) {
   v=(0,0)
   for(k=1;k<=n;k++) v=_c58b_add(v,_c58b_mul(asarray(A,(ids[i],k)),asarray(V,(k,j)),W),W)
   v=_c58b_add(v,_c58b_neg(_c58b_mul(den,_c58b_int(i==j),W)),W)
   if(v[1]!=0 & !W.error) W.error=5
  }
  out.phase_cost[4]=W.ops
 }
 if(!W.error) {
  out.phase=5;N=asarray_create("real",2);Q=asarray_create("real",2)
  pden=_c58b_shift(den,-re,W)
  for(i=1;i<=n;i++) {
   v=(0,0)
   for(j=1;j<=n;j++) v=_c58b_add(v,_c58b_mul(asarray(V,(i,j)),asarray(A,(ids[j],n+1)),W),W)
   asarray(N,(i,1),v)
  }
  out.primal_sign=J(rows(G),1,.);out.dual_sign=J(n,1,.);ok=1
  for(i=1;i<=rows(G);i++) {
   v=_c58b_neg(_c58b_mul(asarray(A,(i,n+1)),den,W))
   for(j=1;j<=n;j++) v=_c58b_add(v,_c58b_mul(asarray(A,(i,j)),asarray(N,(j,1)),W),W)
   out.primal_sign[i]=v[1]
   if((i<=ne & v[1]!=0) | (i>ne & v[1]>0)) ok=0
  }
  out.phase_cost[5]=W.ops
 }
 if(!W.error) {
  out.phase=6
  for(i=1;i<=n;i++) {
   v=(0,0)
   for(j=1;j<=n;j++) v=_c58b_add(v,_c58b_mul(asarray(C,(1,j)),asarray(V,(j,i)),W),W)
   asarray(Q,(i,1),v);out.dual_sign[i]=v[1]
   if(ids[i]>ne & v[1]<0) ok=0
  }
  for(j=1;j<=n;j++) {
   v=_c58b_neg(_c58b_mul(asarray(C,(1,j)),den,W))
   for(i=1;i<=n;i++) v=_c58b_add(v,_c58b_mul(asarray(A,(ids[i],j)),asarray(Q,(i,1)),W),W)
   if(v[1]!=0 & !W.error) W.error=5
  }
  out.phase_cost[6]=W.ops
 }
 if(!W.error) {
  out.status=ok;out.reason=(ok ? 0 : 3);out.phase=0
  out.numerator=_c58b_pack(N,1,n,128);out.dual_numerator=_c58b_pack(Q,1,n,128)
  out.denominator=pden;out.basis_denominator=den;out.rhs_exponent=re;out.row_exponent=(ex[ids],ce)
  out.inverse_numerator=V
  if(ok) out.face=select(ids,(ids:>ne):&(out.dual_sign:>0))
 }
 else out.reason=10+W.error
 out.ops=W.ops;out.maxlimbs=W.maxlimbs
 return(out)
}

void _c70r_repair(string scalar tag,real rowvector c,real matrix G,real colvector h,
 real scalar ne,struct _C58_Cert scalar fit,real scalar maxpivots,real scalar totalbudget,transmorphic scalar C0,real scalar ce0)
{
 struct _C58_Work scalar W
 real matrix seen,trace
 real colvector ids
 real rowvector plan,key
 real scalar step,stop,total,j,cycle,maxlimbs,solverstatus
 trace=J(0,11,.);seen=J(0,cols(c),.);step=0;stop=0;total=fit.ops;maxlimbs=fit.maxlimbs
 if(rows(fit.basis)==cols(c)) seen=sort(fit.basis,1)'
 while(fit.status==0) {
  if(step>=maxpivots) {stop=6;break;}
  if(total>=totalbudget) {stop=7;break;}
  W.ops=0;W.limit=min((10000000,totalbudget-total));W.limbs=128;W.maxlimbs=0;W.error=0
  plan=_c58r_plan(G,h,ne,fit,W);total=total+W.ops;maxlimbs=max((maxlimbs,W.maxlimbs))
  if(W.error) {stop=8;break;}
  if(plan[5]) {stop=plan[5];break;}
  if(total>=totalbudget) {stop=7;break;}
  ids=fit.basis;ids[plan[4]]=plan[2];key=sort(ids,1)';cycle=0
  for(j=1;j<=rows(seen);j++) if(all(seen[j,.]:==key)) cycle=1
  if(cycle) {stop=5;break;}
  seen=seen\key;solverstatus=fit.solver_status
  fit=_c70u_update(c,G,h,ne,fit,plan[4],plan[2],min((10000000,totalbudget-total)),C0,ce0);fit.solver_status=solverstatus
  total=total+fit.ops;maxlimbs=max((maxlimbs,fit.maxlimbs));step=step+1
  trace=trace\(step,plan[1..4],W.ops,fit.ops,fit.status,fit.reason,fit.phase,total)
 }
 if(fit.status==-1 & !stop) stop=9
 fit.ops=total;fit.maxlimbs=maxlimbs
 fit.repair_stop=stop;fit.repair_steps=step
}

void _c70_recover(real rowvector c,real matrix G,real colvector h,
 real scalar ne,struct _C58_Cert scalar fit,real rowvector report,real rowvector journal,real scalar used_pivots,transmorphic scalar C0,real scalar ce0)
{
 struct _C58_Work scalar W
 struct _C58_Cert scalar trial,aux
 real matrix P,terminal_ray
 real colvector ph,ids,preferred
 real rowvector pc,a
 real scalar n,i,j,total,width,pivots,remaining,tau,active,solverstatus,success,why,mode,hastrial,rayok
 journal=_c62_initial(fit);journal[1]=1
 report=(1,0,0,0,0,.,.,.,0,0,0,fit.maxlimbs,.,fit.status,fit.reason)
 n=cols(c);total=fit.ops;width=fit.maxlimbs
 pivots=(args()>=8 ? used_pivots : (fit.repair_steps<. ? fit.repair_steps : 0));journal[13]=pivots;solverstatus=fit.solver_status
 success=0;why=1;mode=0;hastrial=0;rayok=0
 if(n>16 | rows(G)>256 | rows(G)<n+ne | total>=80000000 | pivots>=16 |
    (fit.status!=0 & (fit.status!=-1 | fit.reason>=11))) {report[4]=1;journal[2]=-1;journal[4]=1;return;}
 if(any(G[rows(G)-n+1..rows(G),.]:!=-I(n)) | any(h[rows(G)-n+1..rows(G)]:!=0)) {report[4]=3;journal[2]=-1;journal[4]=3;return;}
 // First retry: equality-first rank selection from the old basis, then
 // original nonnegative bounds and remaining original constraints.
 if(rows(fit.basis)==n) {
  report[5]=1;preferred=fit.basis
  for(i=rows(G)-n+1;i<=rows(G);i++) if(!any(preferred:==i)) preferred=preferred\i
  for(i=ne+1;i<=rows(G);i++) if(!any(preferred:==i)) preferred=preferred\i
  W.ops=0;W.limit=min((10000000,80000000-total));W.limbs=128;W.maxlimbs=0;W.error=0
  ids=_c61_basis(G,h,ne,preferred,W);total=total+W.ops;width=max((width,W.maxlimbs))
  journal[26]=-1;journal[27]=(W.error ? 10+W.error : 2)
  if(!W.error & rows(ids)==n & total<80000000) {
   trial=_c70b_certify(c,G,h,ne,ids,min((10000000,80000000-total)),C0,ce0);hastrial=1
   total=total+trial.ops;width=max((width,trial.maxlimbs));report[6]=trial.status
   if(trial.status==0 & total<80000000 & pivots<16) {
    remaining=trial.ops
    _c70r_repair("",c,G,h,ne,trial,16-pivots,80000000-total+trial.ops,C0,ce0)
    total=total+trial.ops-remaining
    width=max((width,trial.maxlimbs))
    pivots=pivots+trial.repair_steps;report[10]=trial.repair_steps;journal[9]=trial.repair_steps
   }
   journal[26..28]=(trial.status,trial.reason,trial.repair_stop)
   if(trial.status==1) {success=1;mode=1;}
  }
 }
 journal[8]=report[5];journal[19]=total-journal[17]
 // The equality seed can still be neither primal nor dual feasible.
 // Use a fresh, provably primal-feasible auxiliary vertex instead.
 if(!success & total<80000000 & pivots<16) {
  P=J(0,n+1,.);ph=J(0,1,.)
  if(ne) {P=(G[1..ne,.],J(ne,1,-1))\(-G[1..ne,.],J(ne,1,-1));ph=h[1..ne]\(-h[1..ne]);}
  if(rows(G)>ne+n) {P=P\(G[ne+1..rows(G)-n,.],J(rows(G)-ne-n,1,-1));ph=ph\h[ne+1..rows(G)-n];}
  P=P\(-I(n+1));ph=ph\J(n+1,1,0)
  pc=(J(1,n,0),-1);tau=max((0,max(-ph)));active=0
  for(i=1;i<=rows(ph)-n-1;i++) if(-ph[i]==tau) {active=i;break;}
  if(active==0) active=rows(P)
  ids=(rows(P)-n::rows(P)-1)\active
  aux=_c58b_certify(pc,P,ph,0,ids,min((10000000,80000000-total)))
  total=total+aux.ops;width=max((width,aux.maxlimbs))
  if(aux.status==0 & total<80000000 & pivots<16) {
   remaining=aux.ops
   _c58r_repair("",pc,P,ph,0,aux,16-pivots,80000000-total+remaining)
   total=total+aux.ops-remaining
   width=max((width,aux.maxlimbs))
   pivots=pivots+aux.repair_steps;report[9]=aux.repair_steps
  }
  journal[10]=1;journal[11]=report[9];journal[20]=total-journal[17]-journal[19]
  journal[29..31]=(aux.status,aux.reason,aux.repair_stop)
  report[7]=aux.status;report[13]=aux.repair_stop
  if(aux.status==1) {
   a=_c58b_norm(aux.numerator[n+1,.]);report[8]=(a[1]==0)
   if(a[1]==0 & total<80000000) {
    W.ops=0;W.limit=min((10000000,80000000-total));W.limbs=128;W.maxlimbs=0;W.error=0
    ids=_c61_active(G,h,ne,aux,W);total=total+W.ops;width=max((width,W.maxlimbs))
    journal[32]=-1;journal[33]=(W.error ? 10+W.error : 2)
    if(!W.error & rows(ids)==n & total<80000000) {
     trial=_c70b_certify(c,G,h,ne,ids,min((10000000,80000000-total)),C0,ce0);hastrial=1
     total=total+trial.ops;width=max((width,trial.maxlimbs))
     if(trial.status==0 & total<80000000 & pivots<16) {
      remaining=trial.ops
      _c70r_repair("",c,G,h,ne,trial,16-pivots,80000000-total+remaining,C0,ce0)
      total=total+trial.ops-remaining
      width=max((width,trial.maxlimbs))
      pivots=pivots+trial.repair_steps;report[10]=trial.repair_steps;journal[12]=trial.repair_steps
     }
     journal[32..34]=(trial.status,trial.reason,trial.repair_stop)
     if(trial.status==1) {success=1;mode=2;}
    }
   }
  }
 }
 if(!success & hastrial & total<80000000) {
  W.ops=0;W.limit=min((10000000,80000000-total));W.limbs=128;W.maxlimbs=0;W.error=0
  rayok=_c70_ray(c,G,h,ne,trial,W,terminal_ray,C0,ce0)
  total=total+W.ops;width=max((width,W.maxlimbs));journal[18]=W.ops
  if(rayok) {fit.ray=terminal_ray;fit.ray_num=trial.numerator;fit.ray_den=trial.denominator;}
 }
 journal[2]=(success ? 1 : 3);journal[3]=mode;journal[4]=(success ? 0 : 2)
 journal[13]=pivots;journal[21]=total-journal[17]-journal[18]-journal[19]-journal[20]
 if(rayok) {journal[2]=2;journal[3]=4;journal[4]=0;journal[16]=1;}
 journal[22]=total;journal[23]=width;journal[24]=(success ? trial.repair_stop : fit.repair_stop)
 report[2]=success;report[3]=mode;report[4]=(success ? 0 : 2)
 report[11]=total-fit.ops;report[12]=width
 if(success) {fit=trial;fit.solver_status=solverstatus;}
 fit.ops=total;fit.maxlimbs=width
}

real scalar _c70_ray(real rowvector c,real matrix G,real colvector h,
 real scalar ne,struct _C58_Cert scalar fit,struct _C58_Work scalar W,
 real matrix ray,transmorphic scalar C0,real scalar ce0)
{
 transmorphic scalar A,C,R
 real rowvector ex,ce,a,v,d
 real scalar n,i,j,k,re,cre,ok
 ray=J(0,129,0);n=cols(c)
 if(fit.status!=0 | rows(fit.basis)!=n | rows(fit.primal_sign)!=rows(G) |
    rows(fit.dual_sign)!=n | hasmissing(fit.primal_sign) | hasmissing(fit.dual_sign)) return(0)
 if(ne) if(any(fit.primal_sign[1..ne]:!=0)) return(0)
 if(rows(G)>ne) if(any(fit.primal_sign[ne+1..rows(G)]:>0)) return(0)
 A=_c58i_rows(G,h,ex,re,W);C=C0;ce=ce0
 if(W.error) return(0)
 for(k=1;k<=n;k++) if(fit.basis[k]>ne & fit.dual_sign[k]<0) {
  R=asarray_create("real",2)
  for(j=1;j<=n;j++) asarray(R,(j,1),_c58b_neg(asarray(fit.inverse_numerator,(j,k))))
  ok=1
  for(i=1;i<=rows(G);i++) {
   v=(0,0)
   for(j=1;j<=n;j++) v=_c58b_add(v,_c58b_mul(asarray(A,(i,j)),asarray(R,(j,1)),W),W)
   if((i<=ne & v[1]!=0) | (i>ne & v[1]>0)) ok=0
  }
  v=(0,0)
  for(j=1;j<=n;j++) v=_c58b_add(v,_c58b_mul(asarray(C,(1,j)),asarray(R,(j,1)),W),W)
  if(v[1]<=0) ok=0
  if(W.error) return(0)
  if(ok) {ray=_c58b_pack(R,1,n,128);return(1);}
 }
 return(0)
}
end
