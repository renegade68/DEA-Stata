#!/usr/bin/env python3
"""Cooper hybrid: independent eliminated-slack LPs, conditional optimal faces,
source example, alias/domain/output contracts and isolated failure injection.
Requires licensed Stata 16+ and NumPy/SciPy; optional independent R lpSolve.
"""
from pathlib import Path
import argparse,hashlib,itertools,json,re,shutil,subprocess,time
import numpy as np
from scipy.optimize import linprog
from test_core_retries import WRITER,readbin
from test_output_protection import unchanged
from cooper_certified.exact import audit_pipeline
REPO=Path(__file__).resolve().parents[1]
MATS='input_total input_technical input_congestion output_slack target_inputs target_outputs lambda_stage1 lambda_stage2 status obs_index projection_min projection_max congestion_min congestion_max range_status'.split()
PAPER='input double(id x y)\n1 1 .5\n2 2 2\n3 3 2\n4 5 1\n5 4 1\n6 4 1.2\n7 4.5 1.2\nend\n'

class Oracle:
 def __init__(self):self.count=0;self.retries=0
 def solve(self,c,U=None,b=None,E=None,e=None,bounds=(0,None)):
  self.count+=1
  f=linprog(c,A_ub=U,b_ub=b,A_eq=E,b_eq=e,bounds=bounds,method='highs-ds',options=dict(presolve=False,primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9))
  for method,pre in [('highs-ipm',False),('highs-ds',True),('highs-ipm',True)]:
   valid=f.success and np.isfinite(f.x).all() and (U is None or np.max(U@f.x-b,initial=0)<2e-8) and (E is None or np.max(abs(E@f.x-e),initial=0)<2e-8)
   if valid:break
   self.retries+=1
   f=linprog(c,A_ub=U,b_ub=b,A_eq=E,b_eq=e,bounds=bounds,method=method,options=dict(presolve=pre,primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9))
  assert f.success,f.message
  if U is not None:assert np.max(U@f.x-b,initial=0)<2e-8
  if E is not None:assert np.max(abs(E@f.x-e),initial=0)<2e-8
  return f.fun,f.x

def audit(d,X,Y,weight,ranges=True,exact=False):
 a=readbin(d/'scores.bin');z={n:readbin(d/(n+'.bin')) for n in MATS if (d/(n+'.bin')).exists()};meta=readbin(d/'meta.bin').ravel()
 xs=X.max(0);ys=Y.max(0);x=X/xs;y=Y/ys;n,m=x.shape;s=y.shape[1];o=Oracle();gap=0;rangegap=0
 assert int(meta[1])==sum(a[:,5]!=0);assert np.all((a[:,5]!=0)==~np.isfinite(a[:,0]))
 for k in range(n):
  # Independent multiplier-free formulation: first radial LP, then ONLY lambda.
  U=np.vstack([np.c_[x.T,np.zeros(m)],np.c_[-y.T,y[k]]]);rhs=np.r_[x[k],np.zeros(s)];E=np.r_[np.ones(n),0][None,:]
  f,_=o.solve(np.r_[np.zeros(n),-1],U,rhs,E,np.ones(1),[(0,None)]*n+[(1,None)]);phi=-f
  if np.isfinite(a[k,1]):assert abs(a[k,1]-phi)<2e-7;gap=max(gap,abs(a[k,1]-phi));assert abs(a[k,2]-1/phi)<2e-7
  if not np.isfinite(z['input_total'][k]).all():assert a[k,5] in ([1,2,3,4] if exact else [1,2]);continue
  sm=z['input_total'][k]/xs;sp=z['output_slack'][k]/ys;tx=z['target_inputs'][k]/xs;ty=z['target_outputs'][k]/ys;l=z['lambda_stage1'][k]
  wx=xs if weight=='raw' else 1/x[k];wy=ys if weight=='raw' else 1/y[k];scale=max(max(wx),max(wy));wx=wx/scale;wy=wy/scale
  c=x@wx-y@wy;U=np.vstack([x.T,-y.T]);b=np.r_[x[k],-phi*y[k]];E=np.ones((1,n));e=np.ones(1)
  f,_=o.solve(c,U,b,E,e);actual=c@l;g=abs(f-actual)/(1+abs(f));assert g<2e-7,(d.name,k,'slack objective',g);gap=max(gap,g)
  assert np.max(abs(x.T@l+sm-x[k]))<2e-7 and np.max(abs(y.T@l-sp-a[k,1]*y[k]))<2e-7
  assert abs(l.sum()-1)<2e-7 and min(l)>-2e-7 and min(sm)>-2e-7 and min(sp)>-2e-7
  assert np.max(abs(tx-x[k]+sm))<2e-7 and np.max(abs(ty-a[k,1]*y[k]-sp))<2e-7
  assert abs(a[k,3]-np.mean(sm/x[k]))<2e-7
  if ranges and not exact and z['range_status'][k,0]==0:
   Ef=np.vstack([E,c]);ef=np.r_[e,actual];brange=np.r_[x[k],-a[k,1]*y[k]]
   for j,q in enumerate(np.r_[x.T,y.T]):
    mn,_=o.solve(q,U,brange,Ef,ef);negmx,_=o.solve(-q,U,brange,Ef,ef)
    for typ,want in [('projection_min',mn),('projection_max',-negmx)]:
     rg=abs(z[typ][k,j]/np.r_[xs,ys][j]-want);assert rg<2e-6,(k,j,typ,rg);rangegap=max(rangegap,rg)
  if a[k,5]!=0:assert a[k,5]==3;continue
  delta=z['input_technical'][k]/xs;co=z['input_congestion'][k]/xs;l=z['lambda_stage2'][k]
  wx=xs if weight=='raw' else 1/x[k];wx=wx/max(wx);c=-(x@wx);U=np.r_[x.T,-x.T];b=np.r_[x[k],-tx];E=np.r_[y.T,np.ones((1,n))];e=np.r_[ty,1]
  f,_=o.solve(c,U,b,E,e);actual=c@l;g=abs(f-actual)/(1+abs(f));assert g<2e-7,(d.name,k,'technical objective',g);gap=max(gap,g)
  assert np.max(abs(x.T@l-delta-tx))<2e-7 and np.max(abs(y.T@l-ty))<2e-7
  assert abs(l.sum()-1)<2e-7 and min(l)>-2e-7 and min(delta)>-2e-7 and min(co)>-2e-7
  assert np.max(abs(sm-delta-co))<2e-7 and abs(a[k,0]-np.mean(co/x[k]))<2e-7 and abs(a[k,4]-np.mean(delta/x[k]))<2e-7
  if ranges and not exact and z['range_status'][k,1]==0:
   Ef=np.r_[E,c[None,:]];ef=np.r_[e,actual]
   for j in range(m):
    mn,_=o.solve(x[:,j],U,b,Ef,ef);negmx,_=o.solve(-x[:,j],U,b,Ef,ef)
    for typ,want in [('congestion_min',x[k,j]+negmx),('congestion_max',x[k,j]-mn)]:
     rg=abs(z[typ][k,j]/xs[j]-want);assert rg<2e-6,(k,j,typ,rg);rangegap=max(rangegap,rg)
 return dict(independent_LPs=o.count,oracle_retries=o.retries,failed_DMU=int(meta[1]),range_failed_DMU=int(meta[-1]),projection_multiple=int(meta[-3]),decomposition_multiple=int(meta[-2]),max_objective_gap=gap,max_range_gap=rangegap)

R_CODE=r'''
suppressPackageStartupMessages(library(lpSolve));a<-commandArgs(TRUE);d<-read.csv(a[1]);m<-as.integer(a[2]);s<-as.integer(a[3]);z<-as.matrix(d[,2:(m+s+1)]);z<-z[d$selected==1,,drop=FALSE];sc<-apply(z,2,max);z<-sweep(z,2,sc,'/');x<-z[,1:m,drop=FALSE];y<-z[,(m+1):(m+s),drop=FALSE];n<-nrow(z);ans<-matrix(NA,n,4)
solve<-function(c,A,dir,b){f<-lp('max',c,A,dir,b,transpose.constraints=TRUE);stopifnot(f$status==0);f}
for(k in 1:n){
 f<-solve(c(rep(0,n),1),rbind(cbind(t(x),0),cbind(t(y),-y[k,]),c(rep(1,n),0)),c(rep('<=',m),rep('>=',s),'='),c(x[k,],rep(0,s),1));phi<-f$objval
 w<-if(a[4]=='raw')sc else 1/z[k,];w<-w/max(w);cc<-as.vector(-x%*%w[1:m]+y%*%w[(m+1):(m+s)])
 f<-solve(cc,rbind(t(x),t(y),rep(1,n)),c(rep('<=',m),rep('>=',s),'='),c(x[k,],phi*y[k,],1));tx<-as.vector(t(x)%*%f$solution);ty<-as.vector(t(y)%*%f$solution);sm<-x[k,]-tx
 w<-if(a[4]=='raw')sc[1:m] else 1/x[k,];w<-w/max(w)
 f2<-solve(as.vector(x%*%w),rbind(t(x),t(x),t(y),rep(1,n)),c(rep('<=',m),rep('>=',m),rep('=',s),'='),c(x[k,],tx,ty,1));delta<-as.vector(t(x)%*%f2$solution)-tx
 ans[k,]<-c(phi,mean(sm/x[k,]),mean(delta/x[k,]),mean((sm-delta)/x[k,]))
}
write.table(ans,a[5],row.names=FALSE,col.names=FALSE,sep=',');cat('lpSolve',as.character(packageVersion('lpSolve')),'\n')
'''

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--only');ap.add_argument('--rscript');a=ap.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir();hashes={};results=[]
 for n in sorted('ado/'+f.name for f in (REPO/'ado').iterdir() if f.is_file() and f.suffix in ['.ado','.do','.sthlp']):raw=(REPO/n).read_bytes();(pkg/Path(n).name).write_bytes(raw);hashes[n]=hashlib.sha256(raw).hexdigest()
 (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes),indent=2)+'\n')
 # Observe original LPs/certificates in a private copy. Only the diagnostic
 # directory is injected; equations, proposals, budgets and return logic stay.
 observer=out/'exact_observer';shutil.copytree(pkg,observer)
 f=observer/'_dea_congestion_exact.ado';t=f.read_text();needle='    if "`engine\'"!="simplex" exit 198'
 assert needle in t
 t=t.replace(needle,'    local audit "${COOPER_TEST_AUDIT}"\n'+needle,1);f.write_text(t)

 def run(name,body,check=None,source=pkg):
  if a.only and not re.search(a.only,name):return
  d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset type double\nset processors 1\nset linesize 180\n'
  for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/key).mkdir();head+=f'sysdir set {key} "{d/key}"\n'
  (d/'test.do').write_text(head+f'adopath ++ "{source}"\n'+WRITER+body+'\ndi "TEST_PASS"\n');r=dict(test=name,passed=False);t=time.monotonic()
  try:
   p=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=1200);log=(d/'test.log').read_text(errors='replace');assert p.returncode==0 and 'TEST_PASS' in log.splitlines(),"\n".join(log.splitlines()[-18:])
   if check:r.update(check(d))
   r['passed']=True
  except Exception as e:r['error']=repr(e)
  r['seconds']=round(time.monotonic()-t,3);results.append(r);(out/'results.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(r),flush=True)
 def dump():
  b='matrix meta=(r(N),r(n_failed),r(n_projection_multiple),r(n_decomposition_multiple),r(n_range_failed))\nmata: savebin("meta.bin",st_matrix("meta"))\n'
  for n in MATS:b+=f'mata: savebin("{n}.bin",st_matrix("r({n})"))\n'
  return b+'mata: savebin("scores.bin",st_data(.,tokens("res res_phi res_eff res_total res_technical res_status"),"selected"))\n'
 fixtures=[('paper',np.array([[1,.5],[2,2],[3,2],[5,1],[4,1],[4,1.2],[4.5,1.2]]),1,1),('constant',np.array([[1,1,1],[2,3,1]]),2,1),('nonradial',np.array([[1,1,2],[2,1,1]]),2,1),('outputmix',np.array([[1,1,1,2],[2,1,1,1]]),2,2),('multiple',np.array([[1,3,2],[3,1,2],[4,2,2],[2,4,2],[5,5,1]]),2,1),('singleton',np.array([[3,4,5]]),2,1)]
 industrial=np.array([[389.0, 19.86, 856.02, 90.43, 3.81, 70.47], [412.3, 21.16, 866.85, 94.28, 4.13, 82.07], [423.5, 17.08, 956.04, 104.66, 5.56, 117.78], [417.3, 18.1, 1082.94, 121.24, 9.5, 168.29], [570.0, 12.61, 1273.2, 140.72, 21.44, 273.99], [600.5, 13.45, 1230.72, 129.08, 20.95, 212.89], [641.1, 15.91, 1410.66, 134.83, 30.99, 273.19], [715.3, 23.72, 1728.16, 150.58, 41.29, 407.29], [736.0, 25.97, 2109.57, 157.07, 37.88, 481.02], [745.0, 18.24, 2291.08, 156.53, 41.3, 492.49], [756.0, 14.4, 2533.27, 170.39, 58.93, 704.48], [743.0, 17.5, 2899.16, 184.87, 102.75, 1191.05], [684.0, 25.08, 3520.74, 193.26, 164.27, 1792.0], [691.0, 25.45, 4949.93, 196.88, 198.77, 2183.1], [673.0, 29.35, 4604.0, 195.25, 231.34, 2530.87], [634.0, 23.05, 4722.29, 195.06, 194.9, 2399.09], [596.0, 25.02, 4760.28, 197.81, 203.96, 2668.69]])
 fixtures.extend([('textile',industrial[:,:3],2,1),('automobile',industrial[:,3:],2,1)])
 rng=np.random.default_rng(20261103)
 for m,s in [(1,1),(2,2),(3,2)]:
  z=rng.uniform(.5,7,(8,m+s));z[-1,:m]=z[0,:m]*1.6;z[-1,m:]=z[0,m:]*.7
  for variant in ['base','units','duplicate','near','permuted','sample']:
   zz=z.copy()
   if variant=='units':zz*=10.**np.linspace(-4,4,m+s)
   if variant in ['duplicate','near']:zz[1]=zz[0]*(1 if variant=='duplicate' else 1+1e-8)
   if variant=='permuted':zz=zz[::-1]
   fixtures.append((f'random{m}{s}_{variant}',zz,m,s))
 for key,z,m,s in fixtures:
  n=len(z);sel=np.ones(n);xs=' '.join(f'x{i+1}' for i in range(m));ys=' '.join(f'y{i+1}' for i in range(s))
  if key.endswith('sample'):sel[[1,3]]=0
  f=out/(key+'.csv');np.savetxt(f,np.c_[np.arange(1,n+1),z,sel],delimiter=',',header=','.join(['id']+xs.split()+ys.split()+['selected']),comments='',fmt='%.17g')
  for weight,engine in itertools.product(['raw','normalized'],['simplex','ipm','auto']):
   keyname=f'{key}_{weight}_{engine}';body=f'import delimited using "{f}",clear asdouble\ndea {xs}={ys} if selected,model(congestion) congmethod(cooper) congweight({weight}) congranges dmu(id) gen(res) engine({engine})\n'+dump()
   body+='assert abs(res-res_total+res_technical)<1e-8 if res_status==0\nassert missing(res) if !selected\n'
   if engine=='auto':
    run(keyname,body.split('\ndea ',1)[0]+'\n'+unchanged('dea '+xs+'='+ys+' if selected,model(congestion) congmethod(cooper) engine(auto) gen(res)',198))
    continue
   exact=(engine=='simplex')
   if exact:
    auditroot=out/keyname/'audit'
    body='global COOPER_TEST_AUDIT "'+str(auditroot)+'"\nmkdir "'+str(auditroot)+'"\n'+body

   def check(d,z=z,sel=sel,m=m,weight=weight,key=key,exact=exact):
    r=audit(d,z[sel==1,:m],z[sel==1,m:],weight,exact=exact)
    if exact:r.update(audit_pipeline(d,z[sel==1,:m],z[sel==1,m:],weight))
    if key in ['paper','constant','nonradial','outputmix','singleton','multiple']:assert r['failed_DMU']==0
    scores=readbin(d/'scores.bin')
    if key=='paper':assert np.allclose(readbin(d/'input_congestion.bin').ravel(),[0,0,0,2,1,1,1.5],atol=1e-7)
    if key=='constant':assert abs(scores[-1,0])<1e-7
    if key in ['nonradial','outputmix']:assert abs(scores[-1,0]-.25)<1e-7
    if key=='multiple':assert r['projection_multiple']>0 and r['decomposition_multiple']>0
    return r
   run(keyname,body,check,source=observer if exact else pkg)
  if a.rscript and key in ['paper','constant','nonradial','outputmix','singleton','textile','automobile']:
   (out/'compare.R').write_text(R_CODE)
   for w in ['raw','normalized']:
    def rcheck(d,f=f,m=m,s=s,w=w,key=key):
     p=subprocess.run([a.rscript,str(out/'compare.R'),str(f),str(m),str(s),w,str(d/'r.csv')],capture_output=True,text=True,timeout=120);(d/'console.txt').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr
     rr=np.loadtxt(d/'r.csv',delimiter=',',ndmin=2);ss=readbin(out/f'{key}_{w}_ipm'/'scores.bin')[:,[1,3,4,0]];ok=np.isfinite(ss);assert np.isfinite(rr).all();g=np.max(abs(rr[ok]-ss[ok]),initial=0);assert g<1e-6;return dict(R_values=int(ok.sum()),R_native_missing=int((~ok).sum()),R_max_gap=float(g),R_over_2e7=int((abs(rr[ok]-ss[ok])>2e-7).sum()))
    run(f'R_{key}_{w}','',rcheck)
 # Entire public command spellings, sample/ID/order and replacement contracts.
 for command in ['dea x=y,model(congestion)','dea y,inputs(x) model(congestion)','dea_congestion x=y,','dea_congestion,inputs(x) outputs(y)']:
  for short in [False,True]:
   opts='congm(cooper) congw(normalized) congr' if short else 'congmethod(cooper) congweight(normalized) congranges'
   run('syntax_'+str(len(results)),PAPER+command+' '+opts+' dmu(id) gen(res)\nassert res_status==0\nassert abs(res-.25)<1e-7 in 5\n')
 for invalid in ['congmethod(no)','congmethod(legacy) congranges','congweight(raw)','congmethod(cooper) congweight(no)','congmethod(cooper) stage(1)','congmethod(cooper) ort(input)','congmethod(cooper) fdh','congmethod(cooper) bad(y)','congmethod(cooper) direction(x)','congmethod(cooper) reference','congmethod(cooper) saveref','congmethod(cooper) saveweights','congmethod(cooper) bootstrap','congmethod(cooper) w(x)','congmethod(cooper) radialtolerance(.1)']:
  run('reject_option_'+str(len(results)),PAPER+unchanged('dea x=y,model(congestion) gen(res) '+invalid,198))
 for opt in ['congmethod(cooper)','congweight(raw)','congmatrix','congranges']:run('wrong_model_'+opt.split('(')[0],PAPER+unchanged('dea x=y, '+opt,198))
 for field,value in itertools.product(['x','y'],['0','-1']):run('domain_'+field+value,PAPER+f'replace {field}={value} in 2\ngen str8 res="old"\n'+unchanged('dea_congestion x=y,congm(cooper) gen(res) replacevars',459))
 for opt,rc in [('gen(x) replacevars',110),('gen(id) replacevars dmu(id)',110),('gen(_dmu)',110),('gen(abcdefghijklmnopqrstuvwxyzabc)',7)]:run('collision_'+str(len(results)),PAPER+unchanged('dea_congestion x=y,congm(cooper) '+opt,rc))
 for suffix in ['', '_phi','_eff','_total','_technical','_status']:
  run('protect'+suffix,PAPER+f'gen str8 res{suffix}="old"\n'+unchanged('dea_congestion x=y,congm(cooper) gen(res) dmu(id)',110))
 run('sample_missing',PAPER+'gen str3 label=string(id)\nreplace label="" in 2\nreplace x=. in 3\ngen byte selected=inrange(id,2,6)&!missing(x)&label!=""\ndea_congestion x=y if id>=2 in 2/6,congm(cooper) congmat dmu(label) gen(res)\nassert r(N)==3\nmatrix ob=r(obs_index)\nassert ob[1,1]==4 & ob[3,1]==6\nassert missing(res) if !selected\nassert res_status==0 if selected\n')
 run('empty',PAPER+unchanged('dea x=y if id<0,model(congestion) congm(cooper)',2000))
 run('missing_id',PAPER+'gen _dmu=999\n'+unchanged('dea x=y,model(congestion) congm(cooper)',110))
 run('saving',PAPER+'gen str8 res="old"\ndea_congestion x=y if id>=3,congm(cooper) gen(res) replacevars saving("space = a, b")\nassert _dmu==id if id>=3\nassert missing(_dmu) if id<3\nassert res_status==0 if id>=3\nlocal avg=r(mean_congestion)\npreserve\nuse "space = a, b.dta",clear\nassert _N==5\nassert res_status==0\nsummarize res,meanonly\nassert abs(r(mean)-`avg\')<1e-8\nrestore\n'+unchanged('dea_congestion x=y,congm(cooper) dmu(id) gen(res) replacevars saving("space = a, b")',602)+'dea_congestion x=y,congm(cooper) dmu(id) gen(res) replacevars saving("space = a, b") replace\nassert res_status==0\n')
 run('saving_badpath',PAPER+'gen str8 res="old"\n'+unchanged('dea x=y,model(congestion) congm(cooper) gen(res) replacevars saving("absent/path/result")',603))
 for tol in ['0','-1','.','.a']:run('tolerance_'+tol,PAPER+unchanged('dea_congestion x=y,congm(cooper) tol('+tol+')',198))
 run('small_matsize',PAPER+'set matsize 10\ndea_congestion x=y,congm(cooper) dmu(id) gen(res) congmat\nassert res_status==0\n')
 # Diagnostics do not change the selected primary solution; tolerance, file atomicity.
 run('optional_diagnostics',PAPER+'dea_congestion x=y,congm(cooper) dmu(id) gen(a)\ndea_congestion x=y,congm(cooper) dmu(id) gen(b) congmat\ndea_congestion x=y,congm(cooper) dmu(id) gen(c) congr\nassert a==b & b==c\nassert a_status==b_status & b_status==c_status\n')
 run('saving_write_failure',PAPER+'mkdir target.dta\ngen str8 res="old"\n'+unchanged('dea_congestion x=y,congm(cooper) gen(res) replacevars saving(target.dta) replace',608))
 for tol in ['1e-4','1e-9','1e-12']:
  run('accepted_tolerance_'+tol,PAPER+f'dea_congestion x=y,congm(cooper) dmu(id) gen(res) tol({tol})\nlocal failed=r(n_failed)\ncount if res_status>0\nassert r(N)==`failed\'\nassert abs(res-.25)<1e-7 if id==5 & res_status==0\ncount if res_status==0\nassert r(N)>0\n')
 run('underflow',PAPER+'replace x=1e-200 in 1\nreplace x=1e200 in 2\n'+unchanged('dea_congestion x=y,congm(cooper)',459))
 run('reject_trace',PAPER+unchanged('dea x=y,model(congestion) congm(cooper) trace',198))
 private=out/'fault_initial';shutil.copytree(pkg,private);p=private/'_dea_cooper.ado';t=p.read_text().replace('if (_dea_nonte_optimize(S,eps,1e-12)) {','if (1) {',1);p.write_text(t)
 run('equivalent_retry',PAPER+'dea_congestion x=y,congm(cooper) congmat dmu(id) gen(res)\nassert r(n_native_retry)==21\nassert r(n_native_recovered)==21\nassert res_status==0\nassert abs(res-.25)<1e-7 in 5\n',source=private)
 # Failure at each stage: successful upstream fields remain; no 1 substitution.
 for phase in [1,2,3]:
  private=out/f'fault{phase}';shutil.copytree(pkg,private);p=private/'_dea_cooper.ado';t=p.read_text()
  needle={1:'rc=_dea_cooper_lp(c,E,b,lo,hi,tol,z,adj,J(0,v,0),J(0,1,0),retry)',2:'rc=_dea_cooper_lp(c,E,b,lo,hi,tol,z,adj,J(0,v,0),J(0,1,0),retry)',3:"rc=_dea_cooper_lp(c2,F,bb,lo2,hi2,tol,z2,adj,U,xs[i,.]'\\(-tx'),retry)"}[phase]
  pos=[m.start() for m in re.finditer(re.escape(needle),t)][1 if phase==2 else 0];t=t[:pos]+t[pos:].replace(needle,'if(i==5) continue\n        '+needle,1);p.write_text(t)
  body=PAPER+'dea_congestion x=y,congm(cooper) congranges dmu(id) gen(res)\nassert r(n_failed)==1\nassert r(N_complete)==6\nassert missing(res) in 5\n'+f'assert res_status=={phase} in 5\nassert res_status==0 if id!=5\n'+('assert missing(res_phi) in 5\n' if phase==1 else 'assert abs(res_phi-2)<1e-7 in 5\n')+('assert abs(res_total-.5)<1e-7 in 5\n' if phase==3 else 'assert missing(res_total) in 5\n')
  run(f'injected_phase{phase}',body,source=private)
 private=out/'fault_ranges';shutil.copytree(pkg,private);p=private/'_dea_cooper.ado';t=p.read_text().replace('R=J(2,rows(Q),.)','R=J(2,rows(Q),.)\n    return(R)',1);p.write_text(t)
 run('injected_ranges',PAPER+'dea_congestion x=y,congm(cooper) congranges dmu(id) gen(res)\nassert r(n_failed)==0\nassert r(n_range_failed)==7\nassert res_status==0\nassert abs(res-.25)<1e-7 in 5\n',source=private)
 summary=dict(total=len(results),passed=sum(r['passed'] for r in results),failed=sum(not r['passed'] for r in results));(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary);return bool(summary['failed'])
if __name__=='__main__':raise SystemExit(main())
