*! version 1.0.0 06nov2026
*! Bounded integer primal/dual certificate for the Cooper simplex path.
version 16.0
capture mata: assert(_dea_simplex_audit_ready()==16)
if _rc {
    quietly findfile _dea_simplex.ado
    quietly run "`r(fn)'"
}
mata:
mata set matastrict on
struct _C44_Work {
 real scalar ops,limit,limbs,maxlimbs,error
}
struct _C44_Cert {
 real scalar status,reason,solver_status,ops,maxlimbs
 real scalar phase
 real rowvector phase_cost
 transmorphic scalar inverse_numerator
 real colvector basis,face,primal_sign,dual_sign
 real matrix numerator,dual_numerator
 real rowvector denominator,basis_denominator,row_exponent
 real scalar rhs_exponent
}
// Signed base-2^26 digits, little endian; first element is sign.
real rowvector _c44_norm(real rowvector a0)
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
real scalar _c44_tick(struct _C44_Work scalar W,real scalar cost,real scalar limbs)
{
 if(W.error) return(0)
 W.ops=W.ops+cost;W.maxlimbs=max((W.maxlimbs,limbs))
 if(W.ops>W.limit) W.error=1
 if(limbs>W.limbs) W.error=2
 return(W.error==0)
}
real rowvector _c44_int(real scalar x)
{
 real scalar a,b
 if(x==0) return((0,0))
 a=abs(x);b=floor(a/67108864)
 if(b==0) return((sign(x),a))
 if(b<67108864) return((sign(x),a-b*67108864,b))
 return((sign(x),a-b*67108864,b-floor(b/67108864)*67108864,floor(b/67108864)))
}
real scalar _c44_cmp(real rowvector a,real rowvector b)
{
 real scalar i
 if(cols(a)!=cols(b)) return(sign(cols(a)-cols(b)))
 for(i=cols(a);i>=2;i--) if(a[i]!=b[i]) return(sign(a[i]-b[i]))
 return(0)
}
real rowvector _c44_neg(real rowvector a0)
{
 real rowvector a

 a=a0
 a[1]=-a[1];return(a)
}
real rowvector _c44_add(real rowvector a0,real rowvector b0,struct _C44_Work scalar W)
{
 real rowvector a,b

 real scalar n,i,x,y,v,carry,sgn,cmp
 real rowvector out,tmp
 a=a0;b=b0
 if(!_c44_tick(W,max((cols(a),cols(b))),max((cols(a),cols(b)))-1)) return((0,0))
 if(a[1]==0) return(b)
 if(b[1]==0) return(a)
 n=max((cols(a),cols(b)));out=J(1,n+1,0);carry=0;sgn=a[1]
 if(a[1]!=b[1]) {
  cmp=_c44_cmp(a,b)
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
 out[n+1]=carry;out=_c44_norm(out)
 if(!_c44_tick(W,0,cols(out)-1)) return((0,0))
 return(out)
}
real rowvector _c44_mul(real rowvector a0,real rowvector b0,struct _C44_Work scalar W)
{
 real rowvector a,b

 real scalar i,j,k,v,carry,na,nb
 real rowvector out
 a=a0;b=b0
 if(W.error | a[1]==0 | b[1]==0) return((0,0))
 if(cols(a)==2 & a[2]==1) {b[1]=a[1]*b[1];return(b);}
 if(cols(b)==2 & b[2]==1) {a[1]=a[1]*b[1];return(a);}
 na=cols(a)-1;nb=cols(b)-1
 if(!_c44_tick(W,na*nb,na+nb)) return((0,0))
 out=J(1,na+nb+1,0);out[1]=a[1]*b[1]
 for(i=2;i<=cols(a);i++) {
  carry=0
  for(j=2;j<=cols(b);j++) {
   k=i+j-2;v=out[k]+a[i]*b[j]+carry
   carry=floor(v/67108864);out[k]=v-carry*67108864
  }
  out[i+nb]=carry
 }
 return(_c44_norm(out))
}
real rowvector _c44_div(real rowvector a0,real rowvector b0,struct _C44_Work scalar W)
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
  if(!_c44_tick(W,1,cols(rem))) return((0,0))
  rem=_c44_norm((1,a[i],rem[2..cols(rem)]));q=0
  if(_c44_cmp(rem,b)>=0) {
   na=cols(rem)-1;nb=cols(b)-1
   ha=rem[na+1]*67108864+(na>1 ? rem[na] : 0)
   hb=b[nb+1]*67108864+(nb>1 ? b[nb] : 0)
   if(na-nb<0 | na-nb>1) {W.error=3;return((0,0));}
   q=min((67108863,max((0,floor(ha/hb*(na==nb ? 1 : 67108864))))))
   p=_c44_mul(b,_c44_int(q),W);tries=0
   while(_c44_cmp(p,rem)>0 & !W.error) {
    q=q-1;p=_c44_add(p,_c44_neg(b),W);tries=tries+1
    if(tries>8 | q<0) W.error=3
   }
   next=_c44_add(p,b,W)
   while(_c44_cmp(next,rem)<=0 & !W.error) {
    q=q+1;p=next;next=_c44_add(p,b,W);tries=tries+1
    if(tries>8 | q>=67108864) W.error=3
   }
   rem=_c44_add(rem,_c44_neg(p),W)
  }
  out[i]=q
  if(W.error) return((0,0))
 }
 if(rem[1]!=0) {W.error=3;return((0,0));}
 return(_c44_norm(out))
}
// Exact normal binary64 decomposition x=integer*2^e. Unsupported range=>unknown.
real rowvector _c44_decode(real scalar x,real scalar e,struct _C44_Work scalar W)
{
 real scalar m,k,reconstructed,i
 e=0
 if(x==0) return((0,0))
 if(missing(x) | abs(x)<1.1830521861667747e-271 | abs(x)>8.452712498170644e270) {W.error=4;return((0,0));}
 // Do not form powers via logarithms or general exponentiation: that path can
 // round even nominal powers of two. Each step here is exact binary scaling.
 m=abs(x);k=0
 while(m>=9007199254740992) {m=m/2;e=e+1;k=k+1;}
 while(m!=floor(m)) {m=m*2;e=e-1;k=k+1;}
 while(mod(m,2)==0) {m=m/2;e=e+1;k=k+1;}
 reconstructed=m
 for(i=1;i<=abs(e);i++) reconstructed=reconstructed*(e>0 ? 2 : .5)
 if(reconstructed!=abs(x)) {W.error=4;return((0,0));}
 if(!_c44_tick(W,k+abs(e)+1,3)) return((0,0))
 return(_c44_int(sign(x)*m))
}
real rowvector _c44_shift(real rowvector a0,real scalar e,struct _C44_Work scalar W)
{
 real rowvector a

 real scalar q,k,f,i
 a=a0
 if(a[1]==0) return(a)
 if(e<0 | e!=floor(e)) {W.error=4;return((0,0));}
 q=floor(e/26);k=mod(e,26)
 if(!_c44_tick(W,q+1,cols(a)-1+q+1)) return((0,0))
 if(q) a=(a[1],J(1,q,0),a[2..cols(a)])
 f=1
 for(i=1;i<=k;i++) f=f*2
 return(_c44_mul(a,_c44_int(f),W))
}
// Coefficient-only row scaling. The RHS has a separate shared dyadic scale.
// A_i = G_i * 2^-ex_i, H_i = h_i * 2^(-ex_i-rhs_exponent).
// rhs_exponent <= 0; A and H are exact integers. Never round h.
transmorphic scalar _c46_rows(real matrix G,real colvector h,
 real rowvector exponents,real scalar rhs_exponent,struct _C44_Work scalar W)
{
 transmorphic scalar T
 real scalar i,j,e,emin
 real rowvector a,es,he
 T=asarray_create("real",2);exponents=he=J(1,rows(G),0);rhs_exponent=0
 for(i=1;i<=rows(G);i++) {
  es=J(1,cols(G),0);emin=.
  for(j=1;j<=cols(G);j++) {
   a=_c44_decode(G[i,j],e,W);asarray(T,(i,j),a);es[j]=e
   if(a[1]!=0) emin=min((emin,e))
  }
  if(emin==.) emin=0
  exponents[i]=emin
  for(j=1;j<=cols(G);j++) {
   a=asarray(T,(i,j));a=_c44_shift(a,es[j]-emin,W);asarray(T,(i,j),a)
  }
  a=_c44_decode(h[i],e,W);asarray(T,(i,cols(G)+1),a);he[i]=e-emin
  if(a[1]!=0) rhs_exponent=min((rhs_exponent,he[i]))
  if(W.error) return(T)
 }
 for(i=1;i<=rows(G);i++) {
  a=asarray(T,(i,cols(G)+1));a=_c44_shift(a,he[i]-rhs_exponent,W)
  asarray(T,(i,cols(G)+1),a)
  if(W.error) return(T)
 }
 return(T)
}
// Exact original-row selection when the numerical basis omitted dependent rows.
real colvector _c44_pick(transmorphic scalar A,real colvector choices,
 real scalar n,struct _C44_Work scalar W)
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
    b=_c44_mul(pivot,asarray(T,(i,j)),W)
    v=_c44_mul(a,asarray(T,(k,j)),W)
    v=_c44_div(_c44_add(b,_c44_neg(v),W),previous,W);asarray(T,(i,j),v)
   }
   asarray(T,(i,k),(0,0))
  }
  previous=pivot
  if(W.error) return(J(0,1,.))
 }
 return(ids[1..n])
}
real matrix _c44_pack(transmorphic scalar T,real scalar col,real scalar n,real scalar cap)
{
 real matrix out
 real rowvector a
 real scalar i
 out=J(n,cap+1,0)
 for(i=1;i<=n;i++) {a=asarray(T,(i,col));out[i,1..cols(a)]=a;}
 return(out)
}
struct _C44_Cert scalar _c44_certify(real rowvector c,real matrix G,
 real colvector h,real scalar ne,real colvector choices,real scalar budget)
{
 struct _C44_Work scalar W
 struct _C44_Cert scalar out
 transmorphic scalar A,T,C,N,Q
 real rowvector ex,ce,a,b,v,tmp,pivot,previous,den,pden,cn,qn
 real colvector ids
 real scalar n,i,j,k,p,sgn,ok,re,cre
 W.ops=0;W.limit=budget;W.limbs=128;W.maxlimbs=0;W.error=0
 out.status=-1;out.reason=1;out.solver_status=.;out.ops=0;out.maxlimbs=0
 out.phase=0;out.phase_cost=J(1,6,0);out.inverse_numerator=asarray_create("real",2)
 out.basis=out.face=out.primal_sign=out.dual_sign=J(0,1,.)
 out.numerator=out.dual_numerator=J(0,129,0);out.denominator=out.basis_denominator=J(1,0,.);out.rhs_exponent=.;out.row_exponent=J(1,0,.)
 n=cols(c)
 if(n<1 | n>64 | rows(G)>256 | cols(G)!=n | rows(h)!=rows(G) | cols(h)!=1 |
  hasmissing(G) | hasmissing(h) | hasmissing(c) | rows(choices)<n |
  missing(budget) | budget<0 | ne<0 | ne>rows(G)) return(out)
 for(i=1;i<=rows(choices);i++) if(choices[i]<1 | choices[i]>rows(G) | choices[i]!=floor(choices[i])) return(out)
 out.phase=1;A=_c46_rows(G,h,ex,re,W);C=_c46_rows(c,J(1,1,0),ce,cre,W);out.phase_cost[1]=W.ops
 if(!W.error) {
  out.phase=2
  ids=(rows(choices)==n ? choices : _c44_pick(A,choices,n,W))
  out.phase_cost[2]=W.ops
  if(rows(ids)!=n) {out.reason=(W.error ? 10+W.error : 2);out.ops=W.ops;out.maxlimbs=W.maxlimbs;return(out);}
  out.basis=ids;T=asarray_create("real",2);previous=(1,1)
  for(i=1;i<=n;i++) {
   for(j=1;j<=n;j++) asarray(T,(i,j),asarray(A,(ids[i],j)))
   for(j=1;j<=n;j++) asarray(T,(i,n+j),_c44_int(i==j))
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
     b=_c44_mul(pivot,asarray(T,(i,j)),W);v=_c44_mul(a,asarray(T,(k,j)),W)
     v=_c44_div(_c44_add(b,_c44_neg(v),W),previous,W);asarray(T,(i,j),v)
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
    for(k=1;k<=n;k++) v=_c44_add(v,_c44_mul(asarray(A,(ids[i],k)),asarray(T,(k,n+j)),W),W)
    b=_c44_int(i==j);a=den;a[1]=sgn
    v=_c44_add(v,_c44_neg(_c44_mul(a,b,W)),W)
    if(v[1]!=0 & !W.error) W.error=5
   }
   out.phase_cost[4]=W.ops;if(!W.error) out.phase=5
   // Solve only after the coefficient inverse is certified. Primal and dual
   // denominators differ: Dx=Db*2^-re, N=Adj*H; Q=cint*Adj uses Db.
   pden=_c44_shift(den,-re,W)
   for(i=1;i<=n;i++) {
    a=(0,0)
    for(j=1;j<=n;j++) a=_c44_add(a,_c44_mul(asarray(T,(i,n+j)),asarray(A,(ids[j],n+1)),W),W)
    a[1]=a[1]*sgn;asarray(N,(i,1),a)
   }
   out.primal_sign=J(rows(G),1,.);out.dual_sign=J(n,1,.);ok=1
   for(i=1;i<=rows(G);i++) {
    v=_c44_neg(_c44_mul(asarray(A,(i,n+1)),den,W))
    for(j=1;j<=n;j++) v=_c44_add(v,_c44_mul(asarray(A,(i,j)),asarray(N,(j,1)),W),W)
    out.primal_sign[i]=v[1]
    if((i<=ne & v[1]!=0) | (i>ne & v[1]>0)) ok=0
   }
   out.phase_cost[5]=W.ops;if(!W.error) out.phase=6
   for(i=1;i<=n;i++) {
    qn=(0,0)
    for(j=1;j<=n;j++) qn=_c44_add(qn,_c44_mul(asarray(C,(1,j)),asarray(T,(j,n+i)),W),W)
    qn[1]=qn[1]*sgn;asarray(Q,(i,1),qn);out.dual_sign[i]=qn[1]
    if(ids[i]>ne & qn[1]<0) ok=0
   }
   // Check stationarity again directly in integer row units.
   for(j=1;j<=n;j++) {
    v=_c44_neg(_c44_mul(asarray(C,(1,j)),den,W))
    for(i=1;i<=n;i++) v=_c44_add(v,_c44_mul(asarray(A,(ids[i],j)),asarray(Q,(i,1)),W),W)
    if(v[1]!=0 & !W.error) W.error=5
   }
   out.phase_cost[6]=W.ops
   if(!W.error) {
    out.status=ok;out.reason=(ok ? 0 : 3)
    out.numerator=_c44_pack(N,1,n,128);out.dual_numerator=_c44_pack(Q,1,n,128)
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
struct _C44_Cert scalar _c44_solve(real rowvector c,real matrix E,real colvector b,
 real matrix U,real colvector ub,real rowvector lo,real rowvector hi,
 real scalar tol,real scalar budget,real matrix G,real colvector h)
{
 struct RSM_Result scalar fit
 struct _C44_Cert scalar out
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
 out=_c44_certify(c,G,h,ne,choices,budget);out.solver_status=fit.status
 return(out)
}
real scalar _dea_cooper_cert_ready() return(1)
end
