*! version 1.0.2 08nov2026
*! Bounded exact cross-efficiency certificates, recovery and terminal ray checks.
version 16.0
capture mata: mata drop _c58f_residual_sign()
capture mata: mata drop _c58f_quotient()
capture mata: mata drop _c58f_update_quotient()
capture mata: mata drop _c58f_restore()
capture mata: mata drop _c58f_rowcol()
capture mata: mata drop _c58f_split()
capture mata: mata drop _c58f_wide_product()
capture mata: mata drop _c58f_ratio()
capture mata: mata drop _c58f_identity_needs()
capture mata: mata drop _c58f_identity_cell()
capture mata: mata drop _c58_value()
capture mata: mata drop _c58_point()
capture mata: mata drop _c58_fit()
capture mata: mata drop _c58_map()
capture mata: mata drop _c60_certdiag()
capture mata: mata drop _c58_cert_ready()
capture mata: mata drop _c58d_double()
capture mata: mata drop _c58d_bits()
capture mata: mata drop _c58d_valid()
capture mata: mata drop _c58r_repair()
capture mata: mata drop _c58r_plan()
capture mata: mata drop _c58r_ratio_compare()
capture mata: mata drop _c58r_rowcol()
capture mata: mata drop _c58u_update()
capture mata: mata drop _c58b_solve()
capture mata: mata drop _c58b_certify()
capture mata: mata drop _c58b_pack()
capture mata: mata drop _c58b_pick()
capture mata: mata drop _c58i_rows()
capture mata: mata drop _c58b_shift()
capture mata: mata drop _c58b_decode()
capture mata: mata drop _c58b_div()
capture mata: mata drop _c58b_mul()
capture mata: mata drop _c58b_add()
capture mata: mata drop _c58b_neg()
capture mata: mata drop _c58b_cmp()
capture mata: mata drop _c58b_int()
capture mata: mata drop _c58b_tick()
capture mata: mata drop _c58b_norm()
capture mata: mata drop _C58_Cert()
capture mata: mata drop _C58_Work()
capture mata: mata drop _c61_pick()
capture mata: mata drop _c61_basis()
capture mata: mata drop _c61_active()
capture mata: mata drop _c61_recover()
capture mata: mata drop _c62_empty_trace()
capture mata: mata drop _c62_initial()
capture mata: mata drop _c62_ray()
mata:
mata set matastrict on
struct _C58_Work {
 real scalar ops,limit,limbs,maxlimbs,error
}
struct _C58_Cert {
 real scalar status,reason,solver_status,ops,maxlimbs
 real scalar repair_stop,repair_steps
 real rowvector recovery,trace62,ray_den
 real matrix ray,ray_num
 real scalar phase
 real rowvector phase_cost
 transmorphic scalar inverse_numerator
 real colvector basis,face,primal_sign,dual_sign
 real matrix numerator,dual_numerator
 real rowvector denominator,basis_denominator,row_exponent
 real scalar rhs_exponent
}
// Signed base-2^26 digits, little endian; first element is sign.
real rowvector _c58b_norm(real rowvector a0)
{
 real rowvector a

 real scalar n
 a=a0
 n=cols(a)
 while(n>2 & a[n]==0) n=n-1
 a=a[1..n]
 if(n==2 & a[2]==0) a[1]=0
 return(a)
}
real scalar _c58b_tick(struct _C58_Work scalar W,real scalar cost,real scalar limbs)
{
 if(W.error) return(0)
 W.ops=W.ops+cost;W.maxlimbs=max((W.maxlimbs,limbs))
 if(W.ops>W.limit) W.error=1
 if(limbs>W.limbs) W.error=2
 return(W.error==0)
}
real rowvector _c58b_int(real scalar x)
{
 real scalar a,b
 if(x==0) return((0,0))
 a=abs(x);b=floor(a/67108864)
 if(b==0) return((sign(x),a))
 if(b<67108864) return((sign(x),a-b*67108864,b))
 return((sign(x),a-b*67108864,b-floor(b/67108864)*67108864,floor(b/67108864)))
}
real scalar _c58b_cmp(real rowvector a,real rowvector b)
{
 real scalar i
 if(cols(a)!=cols(b)) return(sign(cols(a)-cols(b)))
 for(i=cols(a);i>=2;i--) if(a[i]!=b[i]) return(sign(a[i]-b[i]))
 return(0)
}
real rowvector _c58b_neg(real rowvector a0)
{
 real rowvector a

 a=a0
 a[1]=-a[1];return(a)
}
real rowvector _c58b_add(real rowvector a0,real rowvector b0,struct _C58_Work scalar W)
{
 real rowvector a,b

 real scalar n,i,x,y,v,carry,sgn,cmp
 real rowvector out,tmp
 a=a0;b=b0
 if(!_c58b_tick(W,max((cols(a),cols(b))),max((cols(a),cols(b)))-1)) return((0,0))
 if(a[1]==0) return(b)
 if(b[1]==0) return(a)
 n=max((cols(a),cols(b)));out=J(1,n+1,0);carry=0;sgn=a[1]
 if(a[1]!=b[1]) {
  cmp=_c58b_cmp(a,b)
  if(cmp==0) return((0,0))
  if(cmp<0) {tmp=a;a=b;b=tmp;sgn=a[1];}
 }
 out[1]=sgn
 for(i=2;i<=n;i++) {
  x=(i<=cols(a) ? a[i] : 0);y=(i<=cols(b) ? b[i] : 0)
  v=x+(a[1]==b[1] ? y : -y)+carry
  if(a[1]==b[1]) {carry=floor(v/67108864);v=v-carry*67108864;}
  else {carry=(v<0 ? -1 : 0);if(v<0) v=v+67108864;}
  out[i]=v
 }
 if(carry<0) {W.error=3;return((0,0));}
 out[n+1]=carry;out=_c58b_norm(out)
 if(!_c58b_tick(W,0,cols(out)-1)) return((0,0))
 return(out)
}
real rowvector _c58b_mul(real rowvector a0,real rowvector b0,struct _C58_Work scalar W)
{
 real rowvector a,b

 real scalar i,j,k,v,carry,na,nb
 real rowvector out
 a=a0;b=b0
 if(W.error | a[1]==0 | b[1]==0) return((0,0))
 if(cols(a)==2 & a[2]==1) {b[1]=a[1]*b[1];return(b);}
 if(cols(b)==2 & b[2]==1) {a[1]=a[1]*b[1];return(a);}
 na=cols(a)-1;nb=cols(b)-1
 if(!_c58b_tick(W,na*nb,na+nb)) return((0,0))
 out=J(1,na+nb+1,0);out[1]=a[1]*b[1]
 for(i=2;i<=cols(a);i++) {
  carry=0
  for(j=2;j<=cols(b);j++) {
   k=i+j-2;v=out[k]+a[i]*b[j]+carry
   carry=floor(v/67108864);out[k]=v-carry*67108864
  }
  out[i+nb]=carry
 }
 return(_c58b_norm(out))
}
real rowvector _c58b_div(real rowvector a0,real rowvector b0,struct _C58_Work scalar W)
{
 real rowvector a,b

 real scalar i,q,na,nb,ha,hb,tries,sgn
 real rowvector rem,out,p,next,one
 a=a0;b=b0
 if(W.error) return((0,0))
 if(b[1]==0) {W.error=3;return((0,0));}
 if(a[1]==0) return((0,0))
 if(cols(b)==2 & b[2]==1) {a[1]=a[1]*b[1];return(a);}
 sgn=a[1]*b[1];a[1]=b[1]=1;out=J(1,cols(a),0);out[1]=sgn;rem=(0,0)
 for(i=cols(a);i>=2;i--) {
  if(!_c58b_tick(W,1,cols(rem))) return((0,0))
  rem=_c58b_norm((1,a[i],rem[2..cols(rem)]));q=0
  if(_c58b_cmp(rem,b)>=0) {
   na=cols(rem)-1;nb=cols(b)-1
   ha=rem[na+1]*67108864+(na>1 ? rem[na] : 0)
   hb=b[nb+1]*67108864+(nb>1 ? b[nb] : 0)
   if(na-nb<0 | na-nb>1) {W.error=3;return((0,0));}
   q=min((67108863,max((0,floor(ha/hb*(na==nb ? 1 : 67108864))))))
   p=_c58b_mul(b,_c58b_int(q),W);tries=0
   while(_c58b_cmp(p,rem)>0 & !W.error) {
    q=q-1;p=_c58b_add(p,_c58b_neg(b),W);tries=tries+1
    if(tries>8 | q<0) W.error=3
   }
   next=_c58b_add(p,b,W)
   while(_c58b_cmp(next,rem)<=0 & !W.error) {
    q=q+1;p=next;next=_c58b_add(p,b,W);tries=tries+1
    if(tries>8 | q>=67108864) W.error=3
   }
   rem=_c58b_add(rem,_c58b_neg(p),W)
  }
  out[i]=q
  if(W.error) return((0,0))
 }
 if(rem[1]!=0) {W.error=3;return((0,0));}
 return(_c58b_norm(out))
}
// Default coefficient domain is unchanged. Failure-only callers may opt in to
// finite binary64 magnitudes, including subnormals; integer budgets are unchanged.
real rowvector _c58b_decode(real scalar x,real scalar e,struct _C58_Work scalar W, | real scalar wide_decode)
{
 real scalar m,k,reconstructed,i
 e=0
 if(x==0) return((0,0))
 if(missing(x) | missing(abs(x)) | (wide_decode!=1 & (abs(x)<1.1830521861667747e-271 | abs(x)>8.452712498170644e270))) {W.error=4;return((0,0));}
 // Do not form powers via logarithms or general exponentiation: that path can
 // round even nominal powers of two. Each step here is exact binary scaling.
 m=abs(x);k=0
 while(m>=9007199254740992) {m=m/2;e=e+1;k=k+1;}
 while(m!=floor(m)) {m=m*2;e=e-1;k=k+1;}
 while(mod(m,2)==0) {m=m/2;e=e+1;k=k+1;}
 reconstructed=m
 for(i=1;i<=abs(e);i++) reconstructed=reconstructed*(e>0 ? 2 : .5)
 if(reconstructed!=abs(x)) {W.error=4;return((0,0));}
 if(!_c58b_tick(W,k+abs(e)+1,3)) return((0,0))
 return(_c58b_int(sign(x)*m))
}
real rowvector _c58b_shift(real rowvector a0,real scalar e,struct _C58_Work scalar W)
{
 real rowvector a

 real scalar q,k,f,i
 a=a0
 if(a[1]==0) return(a)
 if(e<0 | e!=floor(e)) {W.error=4;return((0,0));}
 q=floor(e/26);k=mod(e,26)
 if(!_c58b_tick(W,q+1,cols(a)-1+q+1)) return((0,0))
 if(q) a=(a[1],J(1,q,0),a[2..cols(a)])
 f=1
 for(i=1;i<=k;i++) f=f*2
 return(_c58b_mul(a,_c58b_int(f),W))
}
// Coefficient-only row scaling. The RHS has a separate shared dyadic scale.
// A_i = G_i * 2^-ex_i, H_i = h_i * 2^(-ex_i-rhs_exponent).
// rhs_exponent <= 0; A and H are exact integers. Never round h.
transmorphic scalar _c58i_rows(real matrix G,real colvector h,
 real rowvector exponents,real scalar rhs_exponent,struct _C58_Work scalar W, | real scalar wide_decode)
{
 transmorphic scalar T
 real scalar i,j,e,emin
 real rowvector a,es,he
 T=asarray_create("real",2);exponents=he=J(1,rows(G),0);rhs_exponent=0
 for(i=1;i<=rows(G);i++) {
  es=J(1,cols(G),0);emin=.
  for(j=1;j<=cols(G);j++) {
   a=_c58b_decode(G[i,j],e,W,wide_decode);asarray(T,(i,j),a);es[j]=e
   if(a[1]!=0) emin=min((emin,e))
  }
  if(emin==.) emin=0
  exponents[i]=emin
  for(j=1;j<=cols(G);j++) {
   a=asarray(T,(i,j));a=_c58b_shift(a,es[j]-emin,W);asarray(T,(i,j),a)
  }
  a=_c58b_decode(h[i],e,W,wide_decode);asarray(T,(i,cols(G)+1),a);he[i]=e-emin
  if(a[1]!=0) rhs_exponent=min((rhs_exponent,he[i]))
  if(W.error) return(T)
 }
 for(i=1;i<=rows(G);i++) {
  a=asarray(T,(i,cols(G)+1));a=_c58b_shift(a,he[i]-rhs_exponent,W)
  asarray(T,(i,cols(G)+1),a)
  if(W.error) return(T)
 }
 return(T)
}
// Exact original-row selection when the numerical basis omitted dependent rows.
real colvector _c58b_pick(transmorphic scalar A,real colvector choices,
 real scalar n,struct _C58_Work scalar W)
{
 transmorphic scalar T
 real scalar i,j,k,p,m,tmpidx
 real rowvector a,b,v,pivot,previous,tmp
 real colvector ids
 m=rows(choices);ids=choices;T=asarray_create("real",2);previous=(1,1)
 for(i=1;i<=m;i++) for(j=1;j<=n;j++) asarray(T,(i,j),asarray(A,(ids[i],j)))
 for(k=1;k<=n;k++) {
  p=0
  for(i=k;i<=m;i++) {a=asarray(T,(i,k));if(a[1]!=0) {p=i;break;}}
  if(p==0) return(J(0,1,.))
  if(p!=k) {
   tmpidx=ids[p];ids[p]=ids[k];ids[k]=tmpidx
   for(j=k;j<=n;j++) {tmp=asarray(T,(k,j));asarray(T,(k,j),asarray(T,(p,j)));asarray(T,(p,j),tmp);}
  }
  pivot=asarray(T,(k,k))
  for(i=k+1;i<=m;i++) {
   a=asarray(T,(i,k))
   for(j=k+1;j<=n;j++) {
    b=_c58b_mul(pivot,asarray(T,(i,j)),W)
    v=_c58b_mul(a,asarray(T,(k,j)),W)
    v=_c58b_div(_c58b_add(b,_c58b_neg(v),W),previous,W);asarray(T,(i,j),v)
   }
   asarray(T,(i,k),(0,0))
  }
  previous=pivot
  if(W.error) return(J(0,1,.))
 }
 return(ids[1..n])
}
real matrix _c58b_pack(transmorphic scalar T,real scalar col,real scalar n,real scalar cap)
{
 real matrix out
 real rowvector a
 real scalar i
 out=J(n,cap+1,0)
 for(i=1;i<=n;i++) {a=asarray(T,(i,col));out[i,1..cols(a)]=a;}
 return(out)
}
struct _C58_Cert scalar _c58b_certify(real rowvector c,real matrix G,
 real colvector h,real scalar ne,real colvector choices,real scalar budget, | real scalar wide_decode)
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
 out.phase=1;A=_c58i_rows(G,h,ex,re,W,wide_decode);C=_c58i_rows(c,J(1,1,0),ce,cre,W,wide_decode);out.phase_cost[1]=W.ops
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
// Actual common simplex basis; its floating lower-bound shift is a PROPOSAL.
// Certification above uses only original coefficients, bounds and RHS values.
struct _C58_Cert scalar _c58b_solve(real rowvector c,real matrix E,real colvector b,
 real matrix U,real colvector ub,real rowvector lo,real rowvector hi,
 real scalar tol,real scalar budget,real matrix G,real colvector h, | real scalar wide_decode)
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
 out=_c58b_certify(c,G,h,ne,choices,budget,wide_decode);out.solver_status=fit.status
 return(out)
}

// Exact one-row inverse update. The previous certificate supplies a verified
// integer inverse J/Db. Every new inverse is checked against the full new B.
// q=a*J, s=q[k], Db_new=abs(s),
// J_new=sign(s)*(s*J-J[,k]*q+Db*J[,k]*e_k')/Db.
// Exact dyadic factoring on otherwise over-width checks.
real rowvector _c58f_split(real rowvector a0,real scalar exponent,struct _C58_Work scalar W) {
 real rowvector a,out
 real scalar q,k,f,i,v,carry,n
 a=_c58b_norm(a0);exponent=0
 if(W.error | a[1]==0) return(a)
 if(!_c58b_tick(W,2*(cols(a)-1)+26,cols(a)-1)) return((0,0))
 q=0
 while(a[q+2]==0) q=q+1
 k=0;v=a[q+2]
 while(mod(v,2)==0) {v=v/2;k=k+1;}
 exponent=26*q+k;a=(a[1],a[q+2..cols(a)]);n=cols(a)
 if(k==0) return(a)
 f=1
 for(i=1;i<=k;i++) f=f*2
 out=J(1,n,0);out[1]=a[1];carry=0
 for(i=n;i>=2;i--) {
  v=a[i]+carry*67108864;out[i]=floor(v/f);carry=v-out[i]*f
 }
 if(carry!=0) {W.error=3;return((0,0));}
 return(_c58b_norm(out))
}
real scalar _c58f_wide_product(real rowvector a,real rowvector b,struct _C58_Work scalar W) {
 if(a[1]==0 | b[1]==0) return(0)
 if((cols(a)==2 & a[2]==1) | (cols(b)==2 & b[2]==1)) return(0)
 return(cols(a)+cols(b)-2>W.limbs)
}
real scalar _c58f_ratio(real rowvector a,real rowvector b,real rowvector c,real rowvector d,struct _C58_Work scalar W) {
 real rowvector aa,bb,cc,dd,x,y
 real scalar ea,eb,ec,ed,ex,ey,bx,by,emin
 aa=_c58f_split(a,ea,W);bb=_c58f_split(b,eb,W);cc=_c58f_split(c,ec,W);dd=_c58f_split(d,ed,W)
 x=_c58b_mul(aa,dd,W);y=_c58b_mul(cc,bb,W)
 if(W.error) return(0)
 if(x[1]==0 | y[1]==0) return((x[1]!=0)-(y[1]!=0))
 ex=ea+ed;ey=ec+eb;bx=_c58d_bits(x)+ex;by=_c58d_bits(y)+ey
 if(!_c58b_tick(W,cols(x)+cols(y),max((cols(x),cols(y)))-1)) return(0)
 if(bx!=by) return(sign(bx-by))
 emin=min((ex,ey));x=_c58b_shift(x,ex-emin,W);y=_c58b_shift(y,ey-emin,W)
 return(W.error ? 0 : _c58b_cmp(x,y))
}
real scalar _c58f_identity_needs(transmorphic scalar A,transmorphic scalar V,real colvector ids,real scalar i,real scalar j,real scalar n,struct _C58_Work scalar W) {
 real scalar k
 for(k=1;k<=n;k++) if(_c58f_wide_product(asarray(A,(ids[i],k)),asarray(V,(k,j)),W)) return(1)
 return(0)
}
real rowvector _c58f_identity_cell(transmorphic scalar A,transmorphic scalar V,real colvector ids,real rowvector den,real scalar i,real scalar j,real scalar n,struct _C58_Work scalar W) {
 transmorphic scalar X,Y
 real rowvector a,b,v,ex
 real scalar k,ea,eb,emin
 X=asarray_create("real",1);Y=asarray_create("real",1);ex=J(1,n+1,0);emin=.
 for(k=1;k<=n+1;k++) {
  if(k<=n) {a=asarray(A,(ids[i],k));b=asarray(V,(k,j));}
  else {a=(i==j ? _c58b_neg(den) : (0,0));b=(1,1);}
  a=_c58f_split(a,ea,W);b=_c58f_split(b,eb,W);ex[k]=ea+eb
  asarray(X,k,a);asarray(Y,k,b)
  if(a[1]!=0 & b[1]!=0) emin=min((emin,ex[k]))
 }
 if(emin>=.) return((0,0))
 v=(0,0)
 for(k=1;k<=n+1;k++) {
  a=asarray(X,k);b=asarray(Y,k)
  if(a[1]!=0 & b[1]!=0) {
   a=_c58b_mul(a,b,W);a=_c58b_shift(a,ex[k]-emin,W);v=_c58b_add(v,a,W)
  }
 }
 return(v)
}
// Preserve the full integer value after exact dot-product factoring.
real rowvector _c58f_restore(real rowvector a0,real scalar exponent,struct _C58_Work scalar W) {
 real rowvector a,out
 real scalar n,q,k,f,j,v,carry
 a=_c58b_norm(a0)
 if(W.error) return((0,0))
 if(exponent<0 | exponent>=. | exponent!=floor(exponent)) {W.error=4;return((0,0));}
 if(a[1]==0) return(a)
 n=ceil((_c58d_bits(a)+exponent)/26)
 if(!_c58b_tick(W,cols(a)+n,n)) return((0,0))
 q=floor(exponent/26);k=mod(exponent,26);f=1
 for(j=1;j<=k;j++) f=f*2
 out=J(1,n+1,0);out[1]=a[1];carry=0
 for(j=2;j<=cols(a);j++) {
  v=a[j]*f+carry;carry=floor(v/67108864);out[j+q]=v-carry*67108864
 }
 if(carry) out[cols(a)+q+1]=carry
 return(_c58b_norm(out))
}
real rowvector _c58f_rowcol(transmorphic scalar A,real scalar i,
 transmorphic scalar V,real scalar k,real scalar n,struct _C58_Work scalar W) {
 transmorphic scalar X,Y
 real rowvector a,b,q,ex
 real scalar j,ea,eb,emin
 if(W.error) return((0,0))
 X=asarray_create("real",1);Y=asarray_create("real",1);ex=J(1,n,0);emin=.
 for(j=1;j<=n;j++) {
  a=_c58f_split(asarray(A,(i,j)),ea,W);b=_c58f_split(asarray(V,(j,k)),eb,W)
  if(W.error) return((0,0))
  asarray(X,j,a);asarray(Y,j,b);ex[j]=ea+eb
  if(a[1]!=0 & b[1]!=0) emin=min((emin,ex[j]))
 }
 if(emin>=.) return((0,0))
 q=(0,0)
 for(j=1;j<=n;j++) {
  a=asarray(X,j);b=asarray(Y,j)
  if(a[1]!=0 & b[1]!=0) {
   a=_c58b_mul(a,b,W);a=_c58f_restore(a,ex[j]-emin,W);q=_c58b_add(q,a,W)
   if(W.error) return((0,0))
  }
 }
 return(_c58f_restore(q,emin,W))
}
// Bounded exact (a*b-c*d)/den without full-width products.
real rowvector _c58f_quotient(real rowvector a0,real rowvector b0,
 real rowvector c0,real rowvector d0,real rowvector den0,struct _C58_Work scalar W) {
 real rowvector a,b,c,d,den,x,y,q
 real scalar ea,eb,ec,ed,ez,ex,ey,emin,eq,shift
 if(W.error) return((0,0))
 den=_c58b_norm(den0)
 if(den[1]==0) {W.error=3;return((0,0));}
 a=_c58f_split(a0,ea,W);b=_c58f_split(b0,eb,W)
 c=_c58f_split(c0,ec,W);d=_c58f_split(d0,ed,W);den=_c58f_split(den,ez,W)
 if(W.error) return((0,0))
 ex=ea+eb;ey=ec+ed;emin=.
 if(a[1]!=0 & b[1]!=0) emin=ex
 if(c[1]!=0 & d[1]!=0) emin=min((emin,ey))
 if(emin>=.) return((0,0))
 x=_c58b_mul(a,b,W);y=_c58b_mul(c,d,W)
 if(W.error) return((0,0))
 if(x[1]!=0) x=_c58f_restore(x,ex-emin,W)
 if(y[1]!=0) y=_c58f_restore(y,ey-emin,W)
 q=_c58b_add(x,_c58b_neg(y),W)
 if(W.error | q[1]==0) return((0,0))
 q=_c58f_split(q,eq,W);shift=emin+eq-ez
 if(W.error) return((0,0))
 // Both reduced numerator and denominator are odd; a negative exponent
 // cannot produce an integer quotient, even if an odd factor cancels.
 if(shift<0) {W.error=3;return((0,0));}
 q=_c58b_div(q,den,W)
 if(W.error) return((0,0))
 return(_c58f_restore(q,shift,W))
}
real rowvector _c58f_update_quotient(real rowvector a,real rowvector b,
 real rowvector c,real rowvector d,real rowvector den,struct _C58_Work scalar W) {
 real rowvector x,y
 if(_c58f_wide_product(a,b,W) | _c58f_wide_product(c,d,W)) return(_c58f_quotient(a,b,c,d,den,W))
 x=_c58b_mul(a,b,W);y=_c58b_mul(c,d,W)
 return(_c58b_div(_c58b_add(x,_c58b_neg(y),W),den,W))
}
// Certification needs the residual sign, so a common positive power of2
// may be removed without reconstructing a potentially over-width residual.
real scalar _c58f_residual_sign(transmorphic scalar A,real scalar i,
 transmorphic scalar N,real rowvector rhs,real rowvector den,real scalar n,struct _C58_Work scalar W) {
 transmorphic scalar X,Y
 real scalar j
 real rowvector v
 X=asarray_create("real",2);Y=asarray_create("real",2)
 for(j=1;j<=n;j++) {asarray(X,(1,j),asarray(A,(i,j)));asarray(Y,(j,1),asarray(N,(j,1)));}
 asarray(X,(1,n+1),rhs);asarray(Y,(n+1,1),_c58b_neg(den))
 v=_c58f_identity_cell(X,Y,(1),(0,0),1,1,n+1,W)
 return(v[1])
}
struct _C58_Cert scalar _c58u_update(real rowvector c,real matrix G,
 real colvector h,real scalar ne,struct _C58_Cert scalar prior,
 real scalar position,real scalar added,real scalar budget, | real scalar wide_decode)
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
 out.phase=1;A=_c58i_rows(G,h,ex,re,W,wide_decode);C=_c58i_rows(c,J(1,1,0),ce,cre,W,wide_decode);out.phase_cost[1]=W.ops
 if(!W.error) {
  out.phase=2
  if(re!=prior.rhs_exponent | any(ex[prior.basis]:!=prior.row_exponent[1..n])) W.error=5
  // The old inverse column links the new determinant to the verified old Db.
  for(i=1;i<=n;i++) {
   if(_c58f_identity_needs(A,prior.inverse_numerator,prior.basis,i,position,n,W)) v=_c58f_identity_cell(A,prior.inverse_numerator,prior.basis,old_den,i,position,n,W)
   else {
   v=(0,0)
   for(j=1;j<=n;j++) v=_c58b_add(v,_c58b_mul(asarray(A,(prior.basis[i],j)),asarray(prior.inverse_numerator,(j,position)),W),W)
   v=_c58b_add(v,_c58b_neg(_c58b_mul(old_den,_c58b_int(i==position),W)),W)
   }
   if(v[1]!=0 & !W.error) W.error=5
  }
  ids=prior.basis:+J(n,1,0);ids[position]=added;out.basis=ids;out.phase_cost[2]=W.ops
 }
 if(!W.error) {
  out.phase=3;q=asarray_create("real",2);V=asarray_create("real",2)
  for(j=1;j<=n;j++) {
   if(_c58f_identity_needs(A,prior.inverse_numerator,(added),1,j,n,W)) v=_c58f_rowcol(A,added,prior.inverse_numerator,j,n,W)
   else {
   v=(0,0)
   for(i=1;i<=n;i++) v=_c58b_add(v,_c58b_mul(asarray(A,(added,i)),asarray(prior.inverse_numerator,(i,j)),W),W)
   }
   asarray(q,(1,j),v)
  }
  pivot=asarray(q,(1,position));sgn=pivot[1]
  if(sgn==0 & !W.error) {out.reason=2;out.ops=W.ops;out.maxlimbs=W.maxlimbs;return(out);}
  den=_c58b_mul(pivot,_c58b_int(sgn),W)
  for(i=1;i<=n;i++) for(j=1;j<=n;j++) {
   if(j==position) v=_c58b_mul(asarray(prior.inverse_numerator,(i,j)),_c58b_int(sgn),W)
   else {
    v=_c58f_update_quotient(pivot,asarray(prior.inverse_numerator,(i,j)),asarray(prior.inverse_numerator,(i,position)),asarray(q,(1,j)),old_den,W)
    v=_c58b_mul(v,_c58b_int(sgn),W)
   }
   asarray(V,(i,j),v)
  }
  out.phase_cost[3]=W.ops
 }
 if(!W.error) {
  out.phase=4
  for(i=1;i<=n;i++) for(j=1;j<=n;j++) {
   if(_c58f_identity_needs(A,V,ids,i,j,n,W)) v=_c58f_identity_cell(A,V,ids,den,i,j,n,W)
   else {
   v=(0,0)
   for(k=1;k<=n;k++) v=_c58b_add(v,_c58b_mul(asarray(A,(ids[i],k)),asarray(V,(k,j)),W),W)
   v=_c58b_add(v,_c58b_neg(_c58b_mul(den,_c58b_int(i==j),W)),W)
   }
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
   if(_c58f_identity_needs(A,N,(i),1,1,n,W) | _c58f_wide_product(asarray(A,(i,n+1)),den,W)) out.primal_sign[i]=_c58f_residual_sign(A,i,N,asarray(A,(i,n+1)),den,n,W)
   else {
   v=_c58b_neg(_c58b_mul(asarray(A,(i,n+1)),den,W))
   for(j=1;j<=n;j++) v=_c58b_add(v,_c58b_mul(asarray(A,(i,j)),asarray(N,(j,1)),W),W)
   out.primal_sign[i]=v[1]
   }
   if((i<=ne & out.primal_sign[i]!=0) | (i>ne & out.primal_sign[i]>0)) ok=0
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

real rowvector _c58r_rowcol(transmorphic scalar A,real scalar i,
 struct _C58_Cert scalar fit,real scalar k,real scalar n,struct _C58_Work scalar W)
{
 real scalar j
 real rowvector q
 for(j=1;j<=n;j++) if(_c58f_wide_product(asarray(A,(i,j)),asarray(fit.inverse_numerator,(j,k)),W)) return(_c58f_rowcol(A,i,fit.inverse_numerator,k,n,W))
 q=(0,0)
 for(j=1;j<=n;j++) q=_c58b_add(q,_c58b_mul(asarray(A,(i,j)),asarray(fit.inverse_numerator,(j,k)),W),W)
 return(q)
}
real scalar _c58r_ratio_compare(real rowvector a,real rowvector b,
 real rowvector c,real rowvector d,struct _C58_Work scalar W)
{
 real rowvector x,y
 if(_c58f_wide_product(a,d,W) | _c58f_wide_product(c,b,W)) return(_c58f_ratio(a,b,c,d,W))
 x=_c58b_mul(a,d,W);y=_c58b_mul(c,b,W)
 return(_c58b_cmp(x,y))
}
// Plan fields: mode(1 primal,2 dual),added row,removed row,position,stop reason.
// Only exact integer signs and cross-products determine the exchange.
real rowvector _c58r_plan(real matrix G,real colvector h,real scalar ne,
 struct _C58_Cert scalar fit,struct _C58_Work scalar W, | real scalar wide_decode)
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
 A=_c58i_rows(G,h,ex,re,W,wide_decode)
 if(!W.error & re!=fit.rhs_exponent) {W.error=5;return((0,0,0,0,8));}
 bestnum=(0,0);bestden=(1,1)
 if(primal & !dual) {
  for(k=1;k<=n;k++) if(fit.basis[k]>ne & fit.dual_sign[k]<0) {
   if(pos==0) pos=k
   else if(fit.basis[k]<fit.basis[pos]) pos=k
  }
  for(i=1;i<=ne;i++) {
   q=_c58r_rowcol(A,i,fit,pos,n,W)
   if(q[1]!=0) return((0,0,0,0,4))
  }
  for(i=ne+1;i<=rows(G);i++) if(!any(fit.basis:==i)) {
   q=_c58r_rowcol(A,i,fit,pos,n,W)
   if(q[1]<0) {
    den=_c58b_neg(q);num=_c58b_mul(asarray(A,(i,n+1)),fit.basis_denominator,W)
    for(j=1;j<=n;j++) num=_c58b_add(num,_c58b_neg(_c58b_mul(asarray(A,(i,j)),_c58b_norm(fit.numerator[j,.]),W)),W)
    if(num[1]<0 & !W.error) {W.error=5;return((0,0,0,0,8));}
    cmp=(add ? _c58r_ratio_compare(num,den,bestnum,bestden,W) : -1)
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
   den=_c58r_rowcol(A,add,fit,k,n,W)
   if(den[1]>0) {
    num=_c58b_norm(fit.dual_numerator[k,.])
    cmp=(pos ? _c58r_ratio_compare(num,den,bestnum,bestden,W) : -1)
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
void _c58r_repair(string scalar tag,real rowvector c,real matrix G,real colvector h,
 real scalar ne,struct _C58_Cert scalar fit,real scalar maxpivots,real scalar totalbudget, | real scalar wide_decode)
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
  plan=_c58r_plan(G,h,ne,fit,W,wide_decode);total=total+W.ops;maxlimbs=max((maxlimbs,W.maxlimbs))
  if(W.error) {stop=8;break;}
  if(plan[5]) {stop=plan[5];break;}
  if(total>=totalbudget) {stop=7;break;}
  ids=fit.basis;ids[plan[4]]=plan[2];key=sort(ids,1)';cycle=0
  for(j=1;j<=rows(seen);j++) if(all(seen[j,.]:==key)) cycle=1
  if(cycle) {stop=5;break;}
  seen=seen\key;solverstatus=fit.solver_status
  fit=_c58u_update(c,G,h,ne,fit,plan[4],plan[2],min((10000000,totalbudget-total)),wide_decode);fit.solver_status=solverstatus
  total=total+fit.ops;maxlimbs=max((maxlimbs,fit.maxlimbs));step=step+1
  trace=trace\(step,plan[1..4],W.ops,fit.ops,fit.status,fit.reason,fit.phase,total)
 }
 if(fit.status==-1 & !stop) stop=9
 fit.ops=total;fit.maxlimbs=maxlimbs
 fit.repair_stop=stop;fit.repair_steps=step
}
real scalar _c58d_valid(real rowvector a)
{
 real scalar i,nonzero
 if(cols(a)<2 | cols(a)>513 | hasmissing(a)) return(0)
 if(!any(a[1]:==(-1,0,1))) return(0)
 nonzero=0
 for(i=2;i<=cols(a);i++) {
  if(a[i]<0 | a[i]>=67108864 | a[i]!=floor(a[i])) return(0)
  nonzero=nonzero | a[i]!=0
 }
 return((a[1]!=0)==nonzero)
}
real scalar _c58d_bits(real rowvector a)
{
 real scalar top,k
 top=a[cols(a)];k=1
 while(top>=2) {top=floor(top/2);k=k+1;}
 return(26*(cols(a)-2)+k)
}
// Correct nearest/even binary64 conversion. Status:0ok,1bad input,2resource,
// 3outside symmetric export range |x|<2^1023. Working limbs are export-only,not solver caps.
real rowvector _c58d_double(real rowvector a0,real rowvector d0)
{
 struct _C58_Work scalar W
 real rowvector a,d,left,right,rem,trial,twice
 real scalar e,spacing,k,q,power,cmp,sgn,value,inexact
 if(!_c58d_valid(a0) | !_c58d_valid(d0)) return((.,1,0))
 a=_c58b_norm(a0);d=_c58b_norm(d0)
 if(d[1]!=1) return((.,1,0))
 if(a[1]==0) return((0,0,0))
 sgn=a[1];a[1]=1
 W.ops=0;W.limit=1000000;W.limbs=1024;W.maxlimbs=0;W.error=0
 e=_c58d_bits(a)-_c58d_bits(d)
 left=_c58b_shift(a,max((0,-e)),W);right=_c58b_shift(d,max((0,e)),W)
 if(_c58b_cmp(left,right)<0) e=e-1
 if(e>=1023) return((.,3,0))
 if(e< -1075) return((sgn*0,0,1))
 spacing=(e>= -1022 ? e-52 : -1074)
 rem=_c58b_shift(a,max((0,-spacing)),W);d=_c58b_shift(d,max((0,spacing)),W)
 q=0;power=4503599627370496
 for(k=52;k>=0;k--) {
  trial=_c58b_shift(d,k,W)
  if(_c58b_cmp(rem,trial)>=0) {rem=_c58b_add(rem,_c58b_neg(trial),W);q=q+power;}
  power=power/2
 }
 inexact=rem[1]!=0;twice=_c58b_shift(rem,1,W);cmp=_c58b_cmp(twice,d)
 if(cmp>0 | (cmp==0 & mod(q,2)==1)) q=q+1
 if(W.error) return((.,2,inexact))
 value=q
 for(k=1;k<=abs(spacing);k++) value=value*(spacing>=0 ? 2 : .5)
 value=sgn*value
 if(missing(value) | abs(value)>=8.98846567431158e307) return((.,3,inexact))
 return((value,0,inexact))
}


real scalar _c58_cert_ready() return(6203)
// Exact dyadic scaling only. Free variables use opposite columns. Each
// coefficient must round-trip under the column scaling or certification stops.
real matrix _c58_map(real matrix E,real matrix U,real rowvector c,real scalar free,|real scalar map_reason)
{
 real matrix T,A,B
 real scalar n,j,v,s,k
 map_reason=0
 n=cols(c);T=I(n);if(free) T=(T,-I(n)[.,n])
 A=(E\U\c)*T
 for(j=1;j<=cols(T);j++) {
  v=max(abs(A[.,j]));s=1;k=0
  if(v>0 & v<.) {
   while(v>1 & k<900) {v=v/2;s=s/2;k=k+1;}
   while(v<.5 & k<900) {v=v*2;s=s*2;k=k+1;}
  }
  if(k>=900 | s==0 | s>=.) {map_reason=(k>=900 ? 1 : 2);return(J(0,0,.));}
  B=A[.,j]*s
  if(hasmissing(B) | any(B/s:!=A[.,j])) {map_reason=3;return(J(0,0,.));}
  T[.,j]=T[.,j]*s
 }
 return(T)
}
void _c58_fit(real rowvector c,real matrix E,real colvector b,
 real matrix U,real colvector ub,real scalar tol,
 struct _C58_Cert scalar fit,real matrix G,real colvector h, | real scalar wide_decode)
{
 real rowvector recovery,journal
 fit=_c58b_solve(-c,E,b,U,ub,J(1,cols(c),0),J(1,cols(c),.),tol,10000000,G,h,wide_decode)
 if(fit.status==0) _c58r_repair("",-c,G,h,rows(E),fit,16,80000000,wide_decode)
 journal=_c62_initial(fit)
 if(fit.status!=1) {_c61_recover(-c,G,h,rows(E),fit,recovery,journal,(fit.repair_steps<. ? fit.repair_steps : 0),wide_decode);fit.recovery=recovery;}
 fit.trace62=journal
}
// Convert the certified internal fraction once, combining split columns
// before binary64 rounding. Exported doubles are checked separately.
real rowvector _c58_point(real matrix T,struct _C58_Cert scalar fit, | real scalar wide_decode)
{
 struct _C58_Work scalar W
 real rowvector out,a,N,D,z
 real scalar i,j,e
 out=J(1,rows(T),.);if(fit.status!=1) return(out)
 for(i=1;i<=rows(T);i++) {
  W.ops=0;W.limit=10000000;W.limbs=512;W.maxlimbs=0;W.error=0
  N=(0,0);D=fit.denominator
  // T entries are signed powers of two. Use a shared dyadic exponent.
  real scalar emin
  emin=0
  for(j=1;j<=cols(T);j++) if(T[i,j]!=0) {
   a=_c58b_decode(T[i,j],e,W,wide_decode);emin=min((emin,e))
  }
  for(j=1;j<=cols(T);j++) if(T[i,j]!=0) {
   a=_c58b_decode(T[i,j],e,W,wide_decode)
   a=_c58b_shift(_c58b_mul(a,_c58b_norm(fit.numerator[j,.]),W),e-emin,W)
   N=_c58b_add(N,a,W)
  }
  D=_c58b_shift(D,-emin,W)
  z=_c58d_double(N,D)
  if(W.error | z[2]) return(J(1,rows(T),.))
  out[i]=z[1]
 }
 return(out)
}
real scalar _c58_value(real rowvector c,struct _C58_Cert scalar fit, | real scalar wide_decode)
{
 struct _C58_Work scalar W
 transmorphic scalar C
 real rowvector ex,N,D,z
 real scalar re,j
 if(fit.status!=1) return(.)
 W.ops=0;W.limit=10000000;W.limbs=512;W.maxlimbs=0;W.error=0
 C=_c58i_rows(c,J(1,1,0),ex,re,W,wide_decode);N=(0,0);D=fit.denominator
 for(j=1;j<=cols(c);j++) N=_c58b_add(N,_c58b_mul(asarray(C,(1,j)),_c58b_norm(fit.numerator[j,.]),W),W)
 N=_c58b_shift(N,max((0,ex[1])),W);D=_c58b_shift(D,max((0,-ex[1])),W)
 z=_c58d_double(N,D)
 return(W.error | z[2] ? . : z[1])
}

// Observations of an already computed certificate; never changes its state.
real rowvector _c60_certdiag(struct _C58_Cert scalar fit,real scalar ne)
{
 real rowvector d
 real scalar n
 d=(1,fit.status,fit.reason,fit.repair_stop,fit.repair_steps,fit.maxlimbs,.,.,.)
 // Partial arrays after an arithmetic/resource failure are not established signs.
 if(fit.status!=0 & fit.status!=1) return(d)
 n=rows(fit.primal_sign)
 if(ne>0 & n>=ne) if(!hasmissing(fit.primal_sign[1..ne])) d[7]=any(fit.primal_sign[1..ne]:!=0)
 if(n>ne) if(!hasmissing(fit.primal_sign[ne+1..n])) d[8]=any(fit.primal_sign[ne+1..n]:>0)
 if(rows(fit.basis)>0 & rows(fit.basis)==rows(fit.dual_sign)) if(!hasmissing(fit.dual_sign)) d[9]=any((fit.basis:>ne):&(fit.dual_sign:<0))
 return(d)
}

// Bounded failure-only recovery. Maximize c*x subject to G*x<=h,
// with the first ne rows exact equalities. All signs/ranks use integers.
// Earlier independent rows are retained, so equality rows precede proposals.
real colvector _c61_pick(transmorphic scalar A,real colvector choices,
 real scalar n,struct _C58_Work scalar W)
{
 transmorphic scalar T
 real colvector selected
 real rowvector unused,pivot,previous,a,b,v
 real scalar i,j,k,p,col,m
 T=asarray_create("real",2);m=rows(choices);selected=J(0,1,.)
 unused=J(1,n,1);previous=(1,1)
 for(i=1;i<=m;i++) for(j=1;j<=n;j++) asarray(T,(i,j),asarray(A,(choices[i],j)))
 for(i=1;i<=m;i++) {
  col=0
  for(j=1;j<=n;j++) if(unused[j]) {a=asarray(T,(i,j));if(a[1]!=0) {col=j;break;}}
  if(!col) continue
  selected=selected\choices[i]
  if(rows(selected)==n) return(selected)
  unused[col]=0;pivot=asarray(T,(i,col))
  for(k=i+1;k<=m;k++) {
   a=asarray(T,(k,col))
   for(j=1;j<=n;j++) if(unused[j]) {
    b=_c58b_mul(pivot,asarray(T,(k,j)),W)
    v=_c58b_mul(a,asarray(T,(i,j)),W)
    v=_c58b_div(_c58b_add(b,_c58b_neg(v),W),previous,W)
    asarray(T,(k,j),v)
   }
   asarray(T,(k,col),(0,0))
  }
  previous=pivot
  if(W.error) return(J(0,1,.))
 }
 return(J(0,1,.))
}
real colvector _c61_basis(real matrix G,real colvector h,real scalar ne,
 real colvector preferred,struct _C58_Work scalar W, | real scalar wide_decode)
{
 transmorphic scalar A
 real colvector choices
 real rowvector ex
 real scalar i,re
 choices=(ne ? (1::ne) : J(0,1,.))
 for(i=1;i<=rows(preferred);i++) if(!any(choices:==preferred[i])) choices=choices\preferred[i]
 A=_c58i_rows(G,h,ex,re,W,wide_decode)
 if(W.error) return(J(0,1,.))
 return(_c61_pick(A,choices,cols(G),W))
}
// Recover the ORIGINAL active rows from an exact auxiliary point. No rounded
// point or approximate active-set tolerance is used to construct the seed.
real colvector _c61_active(real matrix G,real colvector h,real scalar ne,
 struct _C58_Cert scalar aux,struct _C58_Work scalar W, | real scalar wide_decode)
{
 transmorphic scalar A,P
 real colvector choices
 real rowvector ex,v,a,N
 real scalar i,j,re,wide,signum
 choices=J(0,1,.);A=_c58i_rows(G,h,ex,re,W,wide_decode)
 for(i=1;i<=rows(G);i++) {
  P=asarray_create("real",2)
  wide=_c58f_wide_product(asarray(A,(i,cols(G)+1)),aux.denominator,W);v=(0,0)
  if(!wide) v=_c58b_neg(_c58b_mul(asarray(A,(i,cols(G)+1)),aux.denominator,W))
  for(j=1;j<=cols(G);j++) {
   N=_c58b_shift(_c58b_norm(aux.numerator[j,.]),-re,W);asarray(P,(j,1),N)
   if(_c58f_wide_product(asarray(A,(i,j)),N,W)) wide=1
   if(!wide) v=_c58b_add(v,_c58b_mul(asarray(A,(i,j)),N,W),W)
  }
  if(W.error) return(J(0,1,.))
  signum=(wide ? _c58f_residual_sign(A,i,P,asarray(A,(i,cols(G)+1)),aux.denominator,cols(G),W) : v[1])
  if(W.error) return(J(0,1,.))
  if((i<=ne & signum!=0) | (i>ne & signum>0)) {W.error=5;return(J(0,1,.));}
  if(signum==0) choices=choices\i
 }
 if(W.error) return(J(0,1,.))
 return(_c61_pick(A,choices,cols(G),W))
}
// Auxiliary feasibility LP: E*x-b<=t, -E*x+b<=t, U*x-u<=t,
// x>=0 and t>=0. A known primal-feasible vertex is x=0,
// t=max(0,abs(b),-u). Exact t=0 is required before returning to the original LP.
void _c61_recover(real rowvector c,real matrix G,real colvector h,
 real scalar ne,struct _C58_Cert scalar fit,real rowvector report,|real rowvector journal,real scalar used_pivots, real scalar wide_decode)
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
  ids=_c61_basis(G,h,ne,preferred,W,wide_decode);total=total+W.ops;width=max((width,W.maxlimbs))
  journal[26]=-1;journal[27]=(W.error ? 10+W.error : 2)
  if(!W.error & rows(ids)==n & total<80000000) {
   trial=_c58b_certify(c,G,h,ne,ids,min((10000000,80000000-total)),wide_decode);hastrial=1
   total=total+trial.ops;width=max((width,trial.maxlimbs));report[6]=trial.status
   if(trial.status==0 & total<80000000 & pivots<16) {
    remaining=trial.ops
    _c58r_repair("",c,G,h,ne,trial,16-pivots,80000000-total+trial.ops,wide_decode)
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
  aux=_c58b_certify(pc,P,ph,0,ids,min((10000000,80000000-total)),wide_decode)
  total=total+aux.ops;width=max((width,aux.maxlimbs))
  if(aux.status==0 & total<80000000 & pivots<16) {
   remaining=aux.ops
   _c58r_repair("",pc,P,ph,0,aux,16-pivots,80000000-total+remaining,wide_decode)
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
    ids=_c61_active(G,h,ne,aux,W,wide_decode);total=total+W.ops;width=max((width,W.maxlimbs))
    journal[32]=-1;journal[33]=(W.error ? 10+W.error : 2)
    if(!W.error & rows(ids)==n & total<80000000) {
     trial=_c58b_certify(c,G,h,ne,ids,min((10000000,80000000-total)),wide_decode);hastrial=1
     total=total+trial.ops;width=max((width,trial.maxlimbs))
     if(trial.status==0 & total<80000000 & pivots<16) {
      remaining=trial.ops
      _c58r_repair("",c,G,h,ne,trial,16-pivots,80000000-total+remaining,wide_decode)
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
  rayok=_c62_ray(c,G,h,ne,trial,W,terminal_ray,wide_decode)
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

// Post-failure exact unboundedness observation. A feasible ORIGINAL point and
// an improving ORIGINAL recession direction are both required. No tolerance.
real rowvector _c62_empty_trace()
{
 real rowvector z
 z=J(1,36,0);z[26..36]=J(1,11,.);z[4]=z[5]=z[6]=z[15]=z[24]=.
 return(z)
}
real rowvector _c62_initial(struct _C58_Cert scalar fit)
{
 real rowvector z
 z=_c62_empty_trace();z[5..6]=(fit.status,fit.reason)
 z[7]=z[13]=(fit.repair_steps<. ? fit.repair_steps : 0)
 z[17]=z[22]=fit.ops;z[23]=fit.maxlimbs;z[24]=fit.repair_stop
 return(z)
}
real scalar _c62_ray(real rowvector c,real matrix G,real colvector h,
 real scalar ne,struct _C58_Cert scalar fit,struct _C58_Work scalar W,
 real matrix ray, | real scalar wide_decode)
{
 transmorphic scalar A,C,R
 real rowvector ex,ce,a,v,d
 real scalar n,i,j,k,re,cre,ok
 ray=J(0,129,0);n=cols(c)
 if(fit.status!=0 | rows(fit.basis)!=n | rows(fit.primal_sign)!=rows(G) |
    rows(fit.dual_sign)!=n | hasmissing(fit.primal_sign) | hasmissing(fit.dual_sign)) return(0)
 if(ne) if(any(fit.primal_sign[1..ne]:!=0)) return(0)
 if(rows(G)>ne) if(any(fit.primal_sign[ne+1..rows(G)]:>0)) return(0)
 A=_c58i_rows(G,h,ex,re,W,wide_decode);C=_c58i_rows(c,J(1,1,0),ce,cre,W,wide_decode)
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
