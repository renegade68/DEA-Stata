*! version 1.1.0 06nov2026
*! Cross-efficiency quantity and recovery observations; no solver dependency.
version 16.0
capture mata: mata drop _c60_diag_ready()
capture mata: mata drop _c60_loss()
capture mata: mata drop _c60_prepare()
capture mata: mata drop _c60_computation()
capture mata: mata drop _c60_empty_cert()
capture mata: mata drop _c60_loss_detail()
capture mata: mata drop _c60_stop_label()
capture mata: mata drop _c60_cert_detail()
capture mata: mata drop _c62_empty_recovery()
capture mata: mata drop _c62_recovery_detail()
capture mata: mata drop _c66_native_detail()
mata:
mata set matastrict on
real scalar _c60_diag_ready() return(7001)
real matrix _c60_loss(real matrix A)
{
 real rowvector s
 real matrix scaled
 real scalar j
 s=colmax(A)
 for(j=1;j<=cols(s);j++) if(s[j]==0) s[j]=1
 scaled=A:/s
 return((A:>0):&(scaled:==0))
}
void _c60_prepare(string scalar inputs,string scalar outputs,string scalar touse,
 real scalar matrices,string scalar in_name,string scalar out_name,
 string scalar obs_name,string scalar counts_name)
{
 real matrix LI,LO
 real colvector obs,lost
 LI=_c60_loss(st_data(.,tokens(inputs),touse))
 LO=_c60_loss(st_data(.,tokens(outputs),touse))
 lost=rowsum((LI,LO));obs=selectindex(st_data(.,touse))
 st_matrix(counts_name,(sum(LI)+sum(LO),sum(LI),sum(LO),sum(lost:>0)))
 if(matrices) {st_matrix(in_name,LI);st_matrix(out_name,LO);st_matrix(obs_name,obs);}
}
void _c60_computation(string scalar stage_name,string scalar in_name,
 string scalar out_name,string scalar obs_name,string scalar result_name)
{
 real colvector stage
 stage=st_matrix(stage_name)
 st_matrix(result_name,(st_matrix(obs_name),stage:>0,stage,rowsum((st_matrix(in_name),st_matrix(out_name)))))
}
real matrix _c60_empty_cert(real scalar n)
{
 real matrix D
 D=J(n,9,.);D[.,1]=J(n,1,0)
 return(D)
}
void _c60_loss_detail(string scalar inputs,string scalar outputs,
 string scalar touse,string scalar id)
{
 real matrix LI,LO
 real colvector obs
 string rowvector vi,vo
 string scalar label,a,b
 real scalar i,j
 vi=tokens(inputs);vo=tokens(outputs);LI=_c60_loss(st_data(.,vi,touse));LO=_c60_loss(st_data(.,vo,touse))
 obs=selectindex(st_data(.,touse))
 for(i=1;i<=rows(obs);i++) if(sum(LI[i,.])+sum(LO[i,.])>0) {
  a="";b=""
  for(j=1;j<=cols(vi);j++) if(LI[i,j]) a=a+" "+vi[j]
  for(j=1;j<=cols(vo);j++) if(LO[i,j]) b=b+" "+vo[j]
  label=strofreal(obs[i],"%21.17g")
  if(id!="") {
   if(substr(st_vartype(id),1,3)=="str") label=st_sdata(obs[i],id)
   else label=strofreal(st_data(obs[i],id),"%21.17g")
  }
  printf("  Global scaling loss: observation %g, DMU %s; inputs:%s; outputs:%s\n",obs[i],label,a,b)
 }
}
string scalar _c60_stop_label(real scalar code)
{
 if(code>=.) return("not attempted")
 if(code==0) return("certificate obtained")
 if(code==1) return("proposal no longer eligible")
 if(code==2) return("proposed basis has primal and dual violations")
 if(code==3) return("proposed basis violates an original equality")
 if(code==4) return("no admissible exact exchange")
 if(code==5) return("repeated basis")
 if(code==6) return("pivot limit")
 if(code==7) return("repair operation limit")
 if(code==8) return("repair arithmetic/resource limit")
 if(code==9) return("updated certificate unavailable")
 return("unknown code")
}
void _c60_cert_detail(real matrix P,real matrix S,real colvector obs)
{
 real scalar i,k
 real rowvector d
 for(i=1;i<=rows(P);i++) for(k=1;k<=2;k++) {
  d=(k==1 ? P[i,.] : S[i,.])
  if(d[1] & d[2]!=1) printf("  Observation %g, certificate stage %g: status %g, reason %g; repair %s (%g), steps %g\n",obs[i],k,d[2],d[3],_c60_stop_label(d[4]),d[4],d[5])
 }
}

real matrix _c62_empty_recovery(real scalar n)
{
 real matrix z
 z=J(n,36,0);z[.,26..36]=J(n,11,.);z[.,4]=z[.,5]=z[.,6]=z[.,15]=z[.,24]=J(n,1,.)
 return(z[.,(1..13,16..24,26..34)])
}
void _c62_recovery_detail(real matrix P,real matrix S,real colvector obs)
{
 real scalar i,k
 real rowvector d
 for(i=1;i<=rows(P);i++) for(k=1;k<=2;k++) {
  d=(k==1 ? P[i,.] : S[i,.])
  if(d[1]) printf("  Observation %g, recovery stage %g: outcome %g, path %g, stop %g; exact ray %g; total exchanges %g, integer operations %g\n",obs[i],k,d[2],d[3],d[4],d[16],d[13],d[22])
 }
}

void _c66_native_detail(real matrix D,real colvector obs,string colvector labels)
{
 real scalar i,k,state,rc
 string scalar why
 for(i=1;i<=rows(D);i++) for(k=1;k<=2;k++) {
  state=D[i,2*k-1];rc=D[i,2*k]
  if(state==0 | state==1 | state>=.) continue
  if(state==2) why="prepared coefficients outside finite numeric range; optimizer not called"
  else if(state==3) why="numerical solver result not accepted"
  else if(state==4) why="self score outside accepted positive finite range"
  else why="unknown numerical state"
  printf("  Numerical warning: observation %g, DMU %s, stage %g: %s; solver rc %g\n",obs[i],labels[i],k,why,rc)
 }
}
end
