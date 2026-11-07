*! version 1.0.0 06nov2026
*! Internal exact fraction/result conversion and checkpoint reader support.
version 16.0
capture mata: assert(_dea_cooper_cert_ready()==1)
if _rc {
 quietly findfile _dea_cooper_cert.ado
 quietly run "`r(fn)'"
}
capture program drop _dea_cooper_result
program define _dea_cooper_result, rclass
 version 16.0
 syntax using/, N(integer) Inputs(integer) Outputs(integer) Index(integer) [POLicy(string) SELection(string) REFerence(string) XSCALE(string)]
 if "`policy'"=="" local policy "stage1"
 local policy=lower(strtrim("`policy'"))
 if !inlist("`policy'","stage1","stage2") {
  di as error "policy() must be stage1 or stage2"
  exit 198
 }
 if `n'<1 | `n'>64 | `inputs'<1 | `inputs'>64 | `outputs'<1 | `outputs'>64 | 2*`n'+`inputs'+`outputs'+1>64 | `index'<1 | `index'>`n' {
  di as error "invalid result dimensions or DMU index"
  exit 198
 }
 if "`policy'"=="stage2" & `"`reference'"'=="" {
  di as error "reference() is required for the experimental stage2 adapter"
  exit 198
 }
 if "`xscale'"=="" {
  forvalues j=1/`inputs' {
   local xscale "`xscale' 1"
  }
 }
 local k: word count `xscale'
 if `k'!=`inputs' {
  di as error "xscale() must contain one positive multiplier per input"
  exit 198
 }
 foreach x of local xscale {
  capture confirm number `x'
  if _rc exit 198
  if `x'<=0 | missing(`x') exit 198
 }
 tempname l1 l2 sc it tt co info
 local _c51_names "`l1' `l2' `sc' `it' `tt' `co' `info'"
 mata: _c51_bridge(st_local("using"),strtoreal(st_local("n")),strtoreal(st_local("inputs")),strtoreal(st_local("outputs")),strtoreal(st_local("index")),st_local("policy"),st_local("selection"),st_local("reference"),strtoreal(tokens(st_local("xscale"))),tokens(st_local("_c51_names")))
 matrix colnames `sc'=congestion phi te total technical
 return matrix lambda_stage1=`l1'
 return matrix lambda_stage2=`l2'
 return matrix scores=`sc'
 return matrix input_total=`it'
 return matrix input_technical=`tt'
 return matrix input_congestion=`co'
 return scalar usable=`info'[1,1]
 return scalar selection_status=`info'[1,2]
 return scalar selection_applied=`info'[1,3]
 return scalar conversion_status=`info'[1,4]
 return scalar rounded_values=`info'[1,5]
 return scalar base_certificate=`info'[1,6]
 return scalar selection_certificate=`info'[1,7]
 return scalar selection_maxpivots=`info'[1,8]
 return scalar reference_distance=`info'[1,9]
 return local requested_policy "`policy'"
 return local result_source "`_c51_source'"
 return local selection_reason "`_c51_reason'"
 return local component_units "lp_units_times_xscale"
 return local certificate_scope "internal_exact_point"
 return local range_scope "original_technical_optimal_face"
 return local cmd "_dea_cooper_result"
end
mata:
mata set matastrict on
real scalar _c51_valid(real rowvector a)
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
real scalar _c51_bits(real rowvector a)
{
 real scalar top,k
 top=a[cols(a)];k=1
 while(top>=2) {top=floor(top/2);k=k+1;}
 return(26*(cols(a)-2)+k)
}
// Correct nearest/even binary64 conversion. Status:0ok,1bad input,2resource,
// 3outside symmetric export range |x|<2^1023. Working limbs are export-only,not solver caps.
real rowvector _c51_double(real rowvector a0,real rowvector d0)
{
 struct _C44_Work scalar W
 real rowvector a,d,left,right,rem,trial,twice
 real scalar e,spacing,k,q,power,cmp,sgn,value,inexact
 if(!_c51_valid(a0) | !_c51_valid(d0)) return((.,1,0))
 a=_c44_norm(a0);d=_c44_norm(d0)
 if(d[1]!=1) return((.,1,0))
 if(a[1]==0) return((0,0,0))
 sgn=a[1];a[1]=1
 W.ops=0;W.limit=1000000;W.limbs=1024;W.maxlimbs=0;W.error=0
 e=_c51_bits(a)-_c51_bits(d)
 left=_c44_shift(a,max((0,-e)),W);right=_c44_shift(d,max((0,e)),W)
 if(_c44_cmp(left,right)<0) e=e-1
 if(e>=1023) return((.,3,0))
 if(e< -1075) return((sgn*0,0,1))
 spacing=(e>= -1022 ? e-52 : -1074)
 rem=_c44_shift(a,max((0,-spacing)),W);d=_c44_shift(d,max((0,spacing)),W)
 q=0;power=4503599627370496
 for(k=52;k>=0;k--) {
  trial=_c44_shift(d,k,W)
  if(_c44_cmp(rem,trial)>=0) {rem=_c44_add(rem,_c44_neg(trial),W);q=q+power;}
  power=power/2
 }
 inexact=rem[1]!=0;twice=_c44_shift(rem,1,W);cmp=_c44_cmp(twice,d)
 if(cmp>0 | (cmp==0 & mod(q,2)==1)) q=q+1
 if(W.error) return((.,2,inexact))
 value=q
 for(k=1;k<=abs(spacing);k++) value=value*(spacing>=0 ? 2 : .5)
 value=sgn*value
 if(missing(value) | abs(value)>=8.98846567431158e307) return((.,3,inexact))
 return((value,0,inexact))
}
struct _C51_Result {
 real rowvector lambda,mu,scores,total,technical,congestion
 real scalar usable,policy_status,applied,conversion,rounded,base_cert,selected_cert,cap,distance
 string scalar reason,source
}
real scalar _c51_take(real rowvector N,real rowvector D,struct _C51_Result scalar R)
{
 real rowvector z
 z=_c51_double(N,D)
 if(z[2]) R.conversion=z[2]
 R.rounded=R.rounded+z[3]
 return(z[1])
}
real rowvector _c51_decode_export(real scalar x,real scalar e,struct _C44_Work scalar W)
{
 real scalar m,k,reconstructed,i
 e=0
 if(x==0) return((0,0))
 if(missing(x)) {W.error=4;return((0,0));}
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
real scalar _c51_scaled(real rowvector N,real rowvector D,real scalar scale,
 struct _C44_Work scalar W,struct _C51_Result scalar R)
{
 real rowvector a,nn,dd
 real scalar e
 a=_c51_decode_export(scale,e,W);nn=_c44_shift(_c44_mul(N,a,W),max((0,e)),W)
 dd=_c44_shift(D,max((0,-e)),W)
 return(_c51_take(nn,dd,R))
}
real matrix _c51_load(string scalar file)
{
 real matrix B,A
 real scalar f,n,m
 B=bufio();bufbyteorder(B,2);f=fopen(file,"r")
 n=fbufget(B,f,"%4b");m=fbufget(B,f,"%4b")
 if(n<0 | m<0 | n*m>2000000) {fclose(f);return(J(0,0,.));}
 A=fbufget(B,f,"%8z",n,m);fclose(f);return(A)
}
real scalar _c51_has(string scalar prefix)
{
 string rowvector keys
 real scalar i
 keys=("c","E","b","U","ub","lo","hi","basis")
 for(i=1;i<=cols(keys);i++) if(!fileexists(prefix+"_"+keys[i]+".bin")) return(0)
 return(1)
}
real scalar _c51_same(real matrix A,real matrix B)
{
 if(rows(A)!=rows(B) | cols(A)!=cols(B)) return(0)
 if(rows(A)*cols(A)==0) return(1)
 return(all(vec(A:==B)))
}
void _c51_model(real matrix E,real colvector b,real matrix U,real colvector ub,
 real rowvector lo,real rowvector hi,real matrix G,real colvector h)
{
 real scalar v,j
 real colvector upper
 v=cols(lo);G=E\U\(-I(v));h=b\ub\(-lo');upper=vec(selectindex(hi':<.))
 for(j=1;j<=rows(upper);j++) {G=G\I(v)[upper[j],.];h=h\hi[upper[j]];}
}
real scalar _c51_associated(real rowvector c,real scalar n,
 struct _C44_Cert scalar base,struct _C44_Cert scalar selected)
{
 struct _C44_Work scalar W
 transmorphic scalar C
 real scalar i,v,re
 real rowvector a,b,x,y,ex
 W.ops=0;W.limit=10000000;W.limbs=512;W.maxlimbs=0;W.error=0
 v=cols(c)
 for(i=1;i<=v-n;i++) {
  a=_c44_mul(_c44_norm(base.numerator[i,.]),selected.denominator,W)
  b=_c44_mul(_c44_norm(selected.numerator[i,.]),base.denominator,W)
  if(a[1]!=b[1] | _c44_cmp(a,b)!=0) return(0)
 }
 C=_c46_rows(c,J(1,1,0),ex,re,W);x=y=(0,0)
 for(i=1;i<=v;i++) {
  a=asarray(C,(1,i));x=_c44_add(x,_c44_mul(a,_c44_norm(base.numerator[i,.]),W),W)
  y=_c44_add(y,_c44_mul(a,_c44_norm(selected.numerator[i,.]),W),W)
 }
 x=_c44_mul(x,selected.denominator,W);y=_c44_mul(y,base.denominator,W)
 return(!W.error & x[1]==y[1] & _c44_cmp(x,y)==0)
}
void _c51_values(real matrix N,real rowvector D,real matrix X,real scalar index,
 real scalar n,real scalar m,real scalar s,real rowvector units,struct _C51_Result scalar R)
{
 struct _C44_Work scalar W
 transmorphic scalar A,T,DELTA
 real scalar i,j,v,re,e,emin
 real rowvector ex,ei,a,b,z,delta,ns,nd,den,pt,pu,pd,nt,nu,co,common
 v=rows(N);W.ops=0;W.limit=10000000;W.limbs=512;W.maxlimbs=0;W.error=0
 A=_c46_rows(X',J(m,1,0),ex,re,W);ei=J(1,m,0);T=asarray_create("real",1);DELTA=asarray_create("real",1);emin=0
 for(j=1;j<=m;j++) {
  a=_c44_decode(X[index,j],e,W);ei[j]=e;asarray(T,j,a)
  emin=min((emin,-e,ex[j]-e));delta=(0,0)
  for(i=1;i<=n;i++) {
   b=_c44_add(_c44_norm(N[v-n+i,.]),_c44_neg(_c44_norm(N[i,.])),W)
   delta=_c44_add(delta,_c44_mul(asarray(A,(j,i)),b,W),W)
  }
  asarray(DELTA,j,delta)
 }
 pt=pu=(0,0);pd=(1,1)
 for(j=1;j<=m;j++) {
  a=asarray(T,j);ns=_c44_norm(N[n+j,.]);delta=asarray(DELTA,j)
  nt=_c44_shift(ns,-ei[j]-emin,W);nu=_c44_shift(delta,ex[j]-ei[j]-emin,W)
  pt=_c44_add(_c44_mul(pt,a,W),_c44_mul(nt,pd,W),W)
  pu=_c44_add(_c44_mul(pu,a,W),_c44_mul(nu,pd,W),W);pd=_c44_mul(pd,a,W)
  den=_c44_shift(D,max((0,-ex[j])),W);nd=_c44_shift(delta,max((0,ex[j])),W)
  co=_c44_add(_c44_shift(ns,max((0,-ex[j])),W),_c44_neg(nd),W)
  if(delta[1]<0 | co[1]<0) {R.conversion=1;return;}
  R.total[j]=_c51_scaled(ns,D,units[j],W,R);R.technical[j]=_c51_scaled(nd,den,units[j],W,R);R.congestion[j]=_c51_scaled(co,den,units[j],W,R)
 }
 common=_c44_mul(_c44_shift(_c44_mul(D,pd,W),-emin,W),_c44_int(m),W)
 if(W.error) {R.conversion=2;return;}
 for(i=1;i<=n;i++) {
  R.lambda[i]=_c51_take(_c44_norm(N[i,.]),D,R)
  R.mu[i]=_c51_take(_c44_norm(N[v-n+i,.]),D,R)
 }
 a=_c44_norm(N[n+m+s+1,.]);if(a[1]!=1) {R.conversion=1;return;}
 R.scores[1]=_c51_take(_c44_add(pt,_c44_neg(pu),W),common,R)
 R.scores[2]=_c51_take(a,D,R);R.scores[3]=_c51_take(D,a,R)
 R.scores[4]=_c51_take(pt,common,R);R.scores[5]=_c51_take(pu,common,R)
 if(W.error) R.conversion=2
}
struct _C51_Result scalar _c51_result(string scalar directory,real scalar n,real scalar m,
 real scalar s,real scalar index,string scalar policy,string scalar selection,string scalar reference,real rowvector units)
{
 struct _C51_Result scalar R
 struct _C44_Cert scalar base,chosen
 real matrix E,U,G,radial,X,F,A,V,SG,cn,E2,U2
 real colvector b,ub,h,fb,vb,sh,b2,ub2
 real rowvector c,lo,hi,cc,ll,hh,lo2,hi2,ref,state,meta,D
 real matrix N
 real scalar v,ne,good
 string scalar prefix
 R.lambda=R.mu=J(1,n,.);R.scores=J(1,5,.);R.total=R.technical=R.congestion=J(1,m,.)
 R.usable=R.applied=R.conversion=R.rounded=R.base_cert=R.selected_cert=0
 R.cap=R.distance=.;R.policy_status=4;R.reason="base_certificate_unavailable";R.source="unavailable"
 prefix=directory+"/decomposition"
 if(!_c51_has(prefix) | !fileexists(directory+"/radial_E.bin")) return(R)
 c=_c51_load(prefix+"_c.bin");E=_c51_load(prefix+"_E.bin");b=_c51_load(prefix+"_b.bin");U=_c51_load(prefix+"_U.bin");ub=_c51_load(prefix+"_ub.bin");lo=_c51_load(prefix+"_lo.bin");hi=_c51_load(prefix+"_hi.bin")
 v=cols(c);ne=rows(E);radial=_c51_load(directory+"/radial_E.bin")
 if(v<2*n+m+s+1 | v>2*n+m+s+2 | rows(radial)!=m+s+1 | cols(radial)!=n+m+s+1 | rows(U)!=2*m) return(R)
 if(rows(c)!=1 | cols(E)!=v | cols(U)!=v | rows(b)!=ne | cols(b)!=1 | rows(ub)!=rows(U) | cols(ub)!=1 | rows(lo)!=1 | cols(lo)!=v | rows(hi)!=1 | cols(hi)!=v) return(R)
 X=radial[1..m,1..n]';if(hasmissing(X) | any(vec(X):<=0)) return(R)
 if(!_c51_same(U,(J(m,v-n,0),X')\(X',J(m,v-2*n,0),-X'))) return(R)
 _c51_model(E,b,U,ub,lo,hi,G,h)
 base=_c44_certify(c,G,h,ne,_c51_load(prefix+"_basis.bin"),10000000)
 if(base.status!=1) return(R)
 R.base_cert=1;R.source="baseline";N=base.numerator;D=base.denominator
 R.policy_status=0;R.reason="not_requested"
 if(policy=="stage2") {
  R.policy_status=4;R.reason="reference_unavailable"
  if(fileexists(reference)) {
   cn=_c51_load(reference)
   if(cols(cn)==n & rows(cn)>=index) {
    ref=cn[index,.]
    if(hasmissing(ref)) {R.policy_status=2;R.reason="reference_missing";}
    else {
     R.policy_status=3;R.reason="selection_uncertified"
     if(_c51_has(selection) & fileexists(selection+"_state.bin") & fileexists(selection+"_meta.bin")) {
      state=_c51_load(selection+"_state.bin");meta=_c51_load(selection+"_meta.bin")
      if(cols(state)==2 & rows(state)==1 & cols(meta)==5 & rows(meta)==1) {
       R.cap=((state[2]==16 | state[2]==128) ? state[2] : .)
       if(state[1]==1 & meta[1]==1 & (R.cap==16 | R.cap==128)) {
        R.policy_status=4;R.reason="selection_model_mismatch"
        F=E;fb=b;if(rows(base.face)) {F=F\G[base.face,.];fb=fb\h[base.face];}
        A=(F,J(rows(F),1,0));cc=(J(1,v,0),-1)
        V=(U,J(rows(U),1,0))\(J(n,v-n,0),I(n),J(n,1,-1))\(J(n,v-n,0),-I(n),J(n,1,-1))
        vb=ub\ref'\(-ref');ll=(lo,0);hh=(hi,.)
        good=_c51_same(_c51_load(selection+"_c.bin"),cc) & _c51_same(_c51_load(selection+"_E.bin"),A) & _c51_same(_c51_load(selection+"_b.bin"),fb) & _c51_same(_c51_load(selection+"_U.bin"),V) & _c51_same(_c51_load(selection+"_ub.bin"),vb) & _c51_same(_c51_load(selection+"_lo.bin"),ll) & _c51_same(_c51_load(selection+"_hi.bin"),hh)
        if(good) {
         _c51_model(A,fb,V,vb,ll,hh,SG,sh)
         chosen=_c44_certify(cc,SG,sh,rows(A),_c51_load(selection+"_basis.bin"),10000000)
         R.reason="selection_certificate_rejected"
         if(chosen.status==1) if(_c51_associated(c,n,base,chosen)) {
          R.policy_status=1;R.applied=R.selected_cert=1;R.source="stage2";R.reason="certified"
          N=chosen.numerator[1..v,.];D=chosen.denominator
          R.distance=_c51_take(_c44_norm(chosen.numerator[v+1,.]),D,R)
         }
        }
       }
      }
     }
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
void _c51_bridge(string scalar directory,real scalar n,real scalar m,real scalar s,
 real scalar index,string scalar policy,string scalar selection,string scalar reference,real rowvector units,string rowvector names)
{
 struct _C51_Result scalar R
 R=_c51_result(directory,n,m,s,index,policy,selection,reference,units)
 st_matrix(names[1],R.lambda);st_matrix(names[2],R.mu);st_matrix(names[3],R.scores)
 st_matrix(names[4],R.total);st_matrix(names[5],R.technical);st_matrix(names[6],R.congestion)
 st_matrix(names[7],(R.usable,R.policy_status,R.applied,R.conversion,R.rounded,R.base_cert,R.selected_cert,R.cap,R.distance))
 st_local("_c51_reason",R.reason);st_local("_c51_source",R.source)
}

end
