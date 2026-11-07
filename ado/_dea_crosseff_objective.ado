*! version 1.0.0 06nov2026
*! Peer-only numerical objectives and scoped build diagnostics; no solver load.
version 16.0
capture mata: mata drop _c70_objective_ready()
capture mata: mata drop _c70_scale()
capture mata: mata drop _c70_numeric()
capture mata: mata drop _c70_counts()
capture mata: mata drop _c70_detail()
mata:
mata set matastrict on
real scalar _c70_objective_ready() return(7001)
real scalar _c70_scale(real scalar x,real scalar exponent,real scalar ops,real scalar limit)
{
 real scalar value,k
 value=x
 for(k=1;k<=abs(exponent);k++) {
  ops=ops+1
  if(ops>limit) return(.)
  value=value*(exponent>0 ? 2 : .5)
  if(missing(value)) return(.)
  if(value==0) break
 }
 return(value)
}
void _c70_numeric(real matrix X,real matrix Y,real scalar index,
 real scalar output,real scalar free,real scalar aggressive,
 real rowvector old,real rowvector c,real rowvector od)
{
 real matrix A
 real scalar n,m,s,d,i,j,maximum,exponent,ops,value,total,t,b,error,rounded,lost,compatible,back
 real rowvector sums,prior
 n=rows(X);m=cols(X);s=cols(Y);d=m+s+free
 od=(1,-1,4,2,0,0);c=J(1,d,.)
 if(index<1 | index>n | cols(old)!=d | hasmissing(X) | hasmissing(Y)) return
 A=(X,(output ? J(n,s,0) : Y));maximum=0;ops=0
 for(i=1;i<=n;i++) if(i!=index) maximum=max((maximum,max(A[i,.])))
 exponent=0;value=maximum
 if(value>0) {
  while(value>.25/max((1,n-1))) {value=value/2;exponent=exponent-1;ops=ops+1;}
  while(value<=.125/max((1,n-1))) {value=value*2;exponent=exponent+1;ops=ops+1;}
 }
 sums=J(1,m+s,0);rounded=lost=0
 for(i=1;i<=n;i++) if(i!=index) for(j=1;j<=m+s;j++) if(A[i,j]>0) {
  value=_c70_scale(A[i,j],exponent,ops,10000000)
  if(value>=.) {od[5]=ops;return;}
  if(value==0) lost=1
  total=sums[j];t=total+value;b=t-total
  error=(total-(t-b))+(value-b);ops=ops+6
  if(error!=0) rounded=1
  sums[j]=t
 }
 maximum=max(sums)
 while(maximum>0 & maximum<1) {sums=sums*2;maximum=maximum*2;exponent=exponent+1;ops=ops+m+s+1;}
 c=(sums[1..m],-sums[m+1..m+s]);if(free) c=(c,0)
 if(aggressive) c=-c
 // Preserve the old numerical path only when it is exactly the same objective
 // up to this dyadic scale; never retain total-minus-self cancellation.
 compatible=!hasmissing(old) & !rounded & !lost
 if(compatible) {
  prior=J(1,d,0)
  for(j=1;j<=d;j++) {
   prior[j]=_c70_scale(old[j],exponent,ops,10000000)
   back=_c70_scale(prior[j],-exponent,ops,10000000)
   if(prior[j]>=. | back!=old[j]) compatible=0
  }
  if(compatible) if(all(prior:==c)) c=old
 }
 if(ops>10000000 | hasmissing(c)) {c=J(1,d,.);od[5]=ops;return;}
 od=(1,1,(lost ? 3 : (rounded ? 2 : 1)),2,ops,0)
}
void _c70_counts(real matrix D)
{
 st_numscalar("r(n_objective_attempted)",sum(D[.,1]:==1))
 st_numscalar("r(n_objective_build_failed)",sum((D[.,1]:==1):&(D[.,2]:!=1)))
 st_numscalar("r(n_objective_numeric_rounded)",sum(D[.,3]:==2))
 st_numscalar("r(n_objective_numeric_lost)",sum(D[.,3]:==3))
 st_numscalar("r(n_objective_exact_verified)",sum(D[.,4]:==1))
 st_numscalar("r(n_objective_warning)",sum((D[.,1]:==1):&((D[.,2]:!=1):|((D[.,3]:>1):&(D[.,4]:!=1)))))
}
void _c70_detail(real matrix D,real colvector obs,string colvector labels)
{
 real scalar i
 for(i=1;i<=rows(D);i++) if(D[i,1]) {
  if(D[i,2]!=1 | D[i,3]>1) printf("  Peer objective: observation %g, DMU %s; build %g, numerical representation %g, scope %g; operations %g, limbs %g\n",obs[i],labels[i],D[i,2],D[i,3],D[i,4],D[i,5],D[i,6])
 }
}
end
