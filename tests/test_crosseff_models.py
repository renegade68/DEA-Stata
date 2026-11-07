#!/usr/bin/env python3
"""Cross-efficiency: independent multiplier LPs, optimal faces and public diagnostics.

Alternative optimal weights need not yield identical peer scores. Check submitted
weights, both objectives, and reconstruction instead of requiring one solver's basis.
"""
from _support import package_files
from pathlib import Path
import argparse,csv,hashlib,itertools,json,math,os,re,subprocess,time
import numpy as np
from scipy.optimize import linprog
from test_core_retries import WRITER,readbin
from test_output_protection import unchanged
REPO=Path(__file__).resolve().parents[1]

def fixture(m,s,variant):
 rng=np.random.default_rng(55212+10*m+s);x=rng.uniform(1,8,(6,m));y=rng.uniform(1,8,(6,s));ids=np.arange(1,7);selected=np.ones(6)
 if variant=='near':
  x=np.ones((3,m));x[1]*=1.000001;x[2]*=2;y=np.ones((3,s));ids=np.arange(1,4);selected=np.ones(3)
 if variant=='singleton':x=x[:1];y=y[:1];ids=ids[:1];selected=selected[:1]
 if variant=='zero':
  if m>1:x[0,0]=0
  if s>1:y[1,0]=0
 if variant=='zerocol':
  if m>1:x[:,1]=0
  if s>1:y[:,1]=0
 if variant=='duplicate':x=np.vstack([x,x[0]]);y=np.vstack([y,y[0]]);ids=np.arange(1,8);selected=np.ones(7)
 if variant=='sample':selected[[0,4]]=0;x[0,0]=np.nan;y[4,0]=-1
 if variant=='units':x=x*np.array([1e9,1e-9])[:m];y=y*np.array([1e-10,1e10])[:s]
 if variant=='tiny':x=x*1e-12;y=y*1e-12
 z=np.c_[ids,x,y,selected]
 if variant=='permuted':z=z[rng.permutation(len(z))]
 return z,['id']+[f'x{i+1}' for i in range(m)]+[f'y{i+1}' for i in range(s)]+['selected']

def problem(x,y,k,rts,ori):
 n,m=x.shape;s=y.shape[1];vrs=rts=='vrs';out=ori=='output';d=m+s+int(vrs)
 c=np.r_[x[k],np.zeros(s)] if out else np.r_[np.zeros(m),-y[k]]
 if vrs:c=np.r_[c,1 if out else -1]
 e=np.r_[np.zeros(m),y[k]] if out else np.r_[x[k],np.zeros(s)]
 if vrs:e=np.r_[e,0]
 u=np.c_[-x,y]
 if vrs:u=np.c_[u,np.full(n,-1 if out else 1)]
 bounds=[(0,None)]*(m+s)+([(None,None)] if vrs else [])
 return c,e[None,:],np.array([1.]),u,np.zeros(n),bounds

def fit(p,c=None,lock=None):
 c0,e,b,u,h,bd=p
 if c is None:c=c0
 if lock is not None:e=np.vstack([e,c0]);b=np.r_[b,lock]
 return linprog(c,A_eq=e,b_eq=b,A_ub=u,b_ub=h,bounds=bd,method='highs',options=dict(primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9))

def secondary(x,y,k,rts,ori,method):
 peers=np.arange(len(x))!=k
 a=np.array([math.fsum(x[peers,j]) for j in range(x.shape[1])]);b=np.array([math.fsum(y[peers,j]) for j in range(y.shape[1])])
 q=np.r_[a,np.zeros(y.shape[1]) if ori=='output' else -b]
 if rts=='vrs':q=np.r_[q,0]
 return -q if method=='aggressive' else q

def verify(d,m,s,rts,ori,method,peer_rule='legacy'):
 z=readbin(d/'data.bin');x=z[:,1:1+m];y=z[:,1+m:1+m+s];n=len(z)
 sx=np.maximum(x.max(0),1e-300);sy=np.maximum(y.max(0),1e-300);sx=np.where(x.max(0)>0,sx,1);sy=np.where(y.max(0)>0,sy,1);xx=x/sx;yy=y/sy
 v=readbin(d/'v.bin');u=readbin(d/'u.bin');w=readbin(d/'w.bin').ravel();st=readbin(d/'status.bin').ravel();peers=readbin(d/'peers.bin').ravel();ce=readbin(d/'ce.bin');scores=readbin(d/'scores.bin');meta=readbin(d/'meta.bin').ravel()
 assert np.isin(st,[0,1,2,3]).all();assert ce.shape==(n,n);assert len(w)==n
 oracle=[];max_self=0.;max_sec=0.;unbounded=0;fallback_feasible=0;max_residual=0.;records=[]
 for k in range(n):
  p=problem(xx,yy,k,rts,ori);f=fit(p);assert f.success,(k,f.message)
  expected=-f.fun if ori=='input' else 1/f.fun;oracle.append(expected)
  assert st[k]!=0,('primary unavailable',k,expected)
  actual=scores[k,1];delta=abs(actual-expected);max_self=max(max_self,delta);assert delta<3e-6,(k,actual,expected)
  a=np.r_[v[k]*sx,u[k]*sy,w[k]] if rts=='vrs' else np.r_[v[k]*sx,u[k]*sy]
  c,e,b,con,rhs,bounds=p;err=max(float(np.max(con@a-rhs)),float(np.max(abs(e@a-b))),float(np.max(-a[:m+s])));max_residual=max(max_residual,err);assert err<3e-6,(k,err,a)
  # The reported self score must agree with the returned primary/final weights.
  val=-actual if ori=='input' else 1/actual
  assert abs(c@a-val)<3e-6*max(1,abs(val)),(k,c@a,val)
  if method=='arbitrary':assert st[k]==1
  else:
   q=secondary(xx,yy,k,rts,ori,method);g=fit(p,q,lock=val)
   assert g.status in [0,2,3],g.message
   if st[k]==2:
    assert g.success,(k,'accepted secondary',g.message)
    gap=abs(q@a-g.fun)/max(1,abs(g.fun));max_sec=max(max_sec,gap);assert gap<3e-6,(k,gap,q@a,g.fun)
   else:
    assert st[k]==3
    if g.status==3:unbounded+=1
    if g.success:fallback_feasible+=1
   records.append(dict(row=k+1,native_status=int(st[k]),independent_secondary_status=int(g.status)))
 num=u@y.T;den=v@x.T;expected_ce=np.full((n,n),np.nan)
 if peer_rule=='limzhu':
  # Independently reconstruct the public input VRS convention from original-unit weights.
  scale=np.maximum(np.maximum(abs(den),abs(w[:,None])),num);vx=den.copy();den=den-w[:,None]
  with np.errstate(divide='ignore',invalid='ignore'):ratio=num/den
  valid=np.isfinite(num)&np.isfinite(vx)&np.isfinite(den)&(vx>=0)&(num>=0)&(den>0)&(den>meta[-1]*scale)&np.isfinite(ratio)&(ratio<=1+meta[-1]);expected_ce[valid]=ratio[valid]
 else:
  valid=np.isfinite(num)&np.isfinite(den)&(num>=0)&(den>(meta[-1] if rts=='vrs' else 0));expected_ce[valid]=num[valid]/den[valid]
 assert np.array_equal(np.isfinite(ce),np.isfinite(expected_ce)),('CE mask',ce,expected_ce)
 assert np.allclose(ce,expected_ce,rtol=1e-8,atol=3e-6,equal_nan=True)
 means=[];counts=[]
 for j in range(n):
  values=ce[np.arange(n)!=j,j];values=values[np.isfinite(values)];counts.append(len(values));means.append(values.mean() if len(values) else np.nan)
 means=np.array(means);assert np.array_equal(peers,np.array(counts));assert np.allclose(scores[:,0],means,rtol=1e-8,atol=1e-10,equal_nan=True)
 mav=np.full(n,np.nan);ok=np.isfinite(means)&(means>0);mav[ok]=scores[ok,1]/means[ok]
 assert np.allclose(scores[:,2],mav,rtol=1e-8,atol=1e-10,equal_nan=True)
 avg=lambda a:np.mean(a[np.isfinite(a)]) if np.isfinite(a).any() else np.nan
 expected_meta=[n,avg(scores[:,0]),avg(scores[:,1]),np.sum(st==0),np.sum(st==3),np.sum(~np.isfinite(means)),n*(n-1)-sum(counts),meta[-1]]
 assert np.allclose(meta,expected_meta,equal_nan=True,rtol=1e-9,atol=1e-10),(meta,expected_meta)
 if rts=='crs':assert np.allclose(np.diag(ce),scores[:,1],rtol=1e-8,atol=3e-6,equal_nan=True)
 (d/'independent.json').write_text(json.dumps(dict(self_scores=oracle,secondary=records),indent=2)+'\n')
 return dict(primary_LPs=n,secondary_LPs=len(records),self_max_gap=max_self,secondary_max_relative_gap=max_sec,max_weight_residual=max_residual,secondary_unbounded=unbounded,secondary_feasible_fallback=fallback_feasible,missing_peer_cells=int(np.sum(~np.isfinite(ce))-np.sum(~np.isfinite(np.diag(ce)))),matrix_cells=n*n)

R_SCRIPT=r'''
suppressPackageStartupMessages(library(deaR))
a=commandArgs(TRUE);jobs=read.csv(a[1]);ans=list();ii=0
for(j in seq_len(nrow(jobs))) {
 q=jobs[j,];z=read.csv(q$file);z=z[z$selected==1,];x=paste0('x',1:q$m);y=paste0('y',1:q$s)
 for(v in c(x,y)) {mx=max(z[[v]]);if(mx>0) z[[v]]=z[[v]]/mx}
 dat=make_deadata(z,dmus=1,inputs=match(x,names(z)),outputs=match(y,names(z)))
 f=cross_efficiency(dat,epsilon=0,orientation=if(q$ori=='input')'io' else 'oo',rts=q$rts,selfapp=FALSE,M2=q$rts=='crs'&&q$ori=='input',M3=FALSE)
 for(k in seq_along(z$id)) {
  primary=f$Arbitrary$efficiency[k];if(q$ori=='output')primary=1/primary
  for(method in c('arbitrary','benevolent','aggressive')) {
   obj=if(method=='arbitrary') f$Arbitrary else if(q$rts=='crs'&&q$ori=='input') f[[if(method=='benevolent')'M2_ben' else 'M2_agg']] else NULL
   if(is.null(obj))next
   v=as.numeric(obj$multiplier_input[k,]);u=as.numeric(obj$multiplier_output[k,]);xx=as.matrix(z[x]);yy=as.matrix(z[y])
   sec=if(method=='arbitrary')NA else sum(v*colSums(xx[seq_len(nrow(xx))!=k,,drop=FALSE]))-sum(u*colSums(yy[seq_len(nrow(yy))!=k,,drop=FALSE]))
   if(method=='aggressive')sec=-sec
   ce=as.numeric(obj$cross_eff[k,]);if(q$ori=='output')ce=1/ce
   for(l in seq_along(z$id)) {ii=ii+1;ans[[ii]]=data.frame(test=q$test,method=method,row=k,col=l,self=primary,secondary=sec,ce=ce[l])}
  }
 }
}
write.csv(do.call(rbind,ans),a[2],row.names=FALSE,na='NA')
writeLines(c(paste('deaR',packageVersion('deaR')),capture.output(sessionInfo())),a[3])
'''

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--source-ado',type=Path);ap.add_argument('--only');ap.add_argument('--rscript');a=ap.parse_args()
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir();hashes={}
 paths=list(a.source_ado.iterdir()) if a.source_ado else [REPO/n for n in package_files(REPO)]
 for p in paths:
  if p.suffix in ['.ado','.do','.sthlp']:raw=p.read_bytes();(pkg/p.name).write_bytes(raw);hashes[p.name]=hashlib.sha256(raw).hexdigest()
 (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,stata=a.stata),indent=2)+'\n');rows=[];jobs=[]
 def run(name,body,check=None):
  if a.only and not re.search(a.only,name):return
  d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset processors 1\nset type double\nset linesize 160\n'
  for f in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/f).mkdir();head+=f'sysdir set {f} "{d/f}"\n'
  head+=f'adopath ++ "{pkg}"\n'+WRITER;(d/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n');r=dict(test=name,passed=False);start=time.monotonic()
  try:
   proc=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=180);log=(d/'test.log').read_text(errors='replace')
   assert proc.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2000:]
   if check:r.update(check(d))
   r['passed']=True
  except Exception as e:r['error']=repr(e)
  r['seconds']=round(time.monotonic()-start,3);rows.append(r);(out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(r),flush=True)
 variants=['base','units','tiny','zero','zerocol','duplicate','sample','permuted','near','singleton']
 for m,s,variant in itertools.product([1,2],[1,2],variants):
  z,cols=fixture(m,s,variant);file=out/f'data_{m}x{s}_{variant}.csv'
  with file.open('w') as f:w=csv.writer(f);w.writerow(cols);w.writerows(['' if np.isnan(v) else format(v,'.17g') for v in row] for row in z)
  xs=' '.join(f'x{i+1}' for i in range(m));ys=' '.join(f'y{i+1}' for i in range(s))
  for ci,(rts,ori,method) in enumerate(itertools.product(['crs','vrs'],['input','output'],['arbitrary','benevolent','aggressive'])):
   name=f'{m}x{s}_{variant}_{rts}_{ori}_{method}'
   if a.only and not re.search(a.only,name):continue
   cmd=f'dea_crosseff {xs}={ys} if selected,' if ci%2 else f'dea_crosseff if selected,inputs({xs}) outputs({ys})'
   tol='1e-9' if ci%2 else '1e-7'
   body=f'import delimited using "{file}",clear varnames(1) asdouble\n{cmd} dmu(id) rts({rts}) ort({ori}) method({method}) tol({tol}) gen(res) matrix detail\n'
   body+='matrix meta=(r(n),r(mean_ce),r(mean_self),r(n_failed_primary),r(n_failed_secondary),r(n_failed_ce),r(n_missing_peer_pairs),r(tolerance))\n'
   for f,mat in [('ce','ce_matrix'),('v','weights_input'),('u','weights_output'),('w','intercept'),('status','stage_status'),('peers','peer_count')]:body+=f'mata: savebin("{f}.bin",st_matrix("r({mat})"))\n'
   body+='mata: savebin("meta.bin",st_matrix("meta"))\nassert missing(res,res_self,res_maverick) if !selected\n'
   body+=f'keep if selected\nmata: savebin("data.bin",st_data(.,tokens("id {xs} {ys}")))\nmata: savebin("scores.bin",st_data(.,tokens("res res_self res_maverick")))\n'
   run(name,body,lambda d,m=m,s=s,rts=rts,ori=ori,method=method:verify(d,m,s,rts,ori,method))
   if variant=='base' and method=='arbitrary':jobs.append(dict(test=name,file=str(file),m=m,s=s,rts=rts,ori=ori))
 data='input double(id x y)\n1 1 1\n2 2 1\n3 3 2\nend\n'
 for bad in ['0','-1','.','.a','.z']:
  run('invalid_tol_'+bad.replace('.','dot'),data+'gen str12 res="keep"\n'+unchanged(f'dea_crosseff x=y,gen(res) replacevars tol({bad})',198))
 for key,setup in [('negative_x','replace x=-1 in 1'),('negative_y','replace y=-1 in 1'),('zero_x','replace x=0 in 1'),('zero_y','replace y=0 in 1')]:run('invalid_data_'+key,data+setup+'\ngen str12 res="keep"\n'+unchanged('dea_crosseff x=y,gen(res) replacevars',459))
 run('empty_sample',data+'gen str12 res="keep"\n'+unchanged('dea_crosseff x=y if id==99,gen(res) replacevars',2000))
 run('reload',f'run "{pkg}/dea_crosseff.ado"\nrun "{pkg}/dea_crosseff.ado"\n'+data+'dea_crosseff x=y,gen(res)\nassert !missing(res)\n')
 run('matrix_state',data+'dea_crosseff x=y,gen(res) matrix\nassert colsof(r(weights_input))==1\nlocal cn: colnames r(weights_input)\nassert "`cn\'"=="x"\ndea_crosseff x=y,gen(res) replacevars\ncapture matrix list r(ce_matrix)\nassert _rc==111\ncapture matrix list r(weights_input)\nassert _rc==111\n')
 # Positive denominator smaller than tol must remain usable; the secondary LP
 # fixes v*x0=1 at a much larger DMU. A feasibility tolerance is not a ratio cutoff.
 run('small_peer_denominator','input double(id x y)\n1 1 1\n2 1e8 1e8\nend\ndea_crosseff x=y,gen(res) tol(1e-7) matrix\nassert r(n_missing_peer_pairs)==0\nassert abs(res-1)<1e-6\n')
 if a.rscript and jobs:
  with (out/'r_jobs.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(jobs[0]));w.writeheader();w.writerows(jobs)
  (out/'compare.R').write_text(R_SCRIPT);p=subprocess.run([a.rscript,str(out/'compare.R'),str(out/'r_jobs.csv'),str(out/'r_scores.csv'),str(out/'R_session.txt')],capture_output=True,env=dict(os.environ,RGL_USE_NULL='TRUE'),timeout=300);(out/'R.log').write_bytes(p.stdout+p.stderr);r=dict(test='R_deaR',passed=False)
  try:
   assert p.returncode==0,p.stderr.decode()[-2000:]
   checks=0;maxgap=0.;direct=0;maxce=0.;sec_checks=0;maxsec=0
   with (out/'r_scores.csv').open() as f:
    for rr in csv.DictReader(f):
     name=rr['test'].replace('_arbitrary','_'+rr['method']);d=out/name
     if not d.exists():continue
     k=int(rr['row'])-1;l=int(rr['col'])-1
     if l==0:
      score=readbin(d/'scores.bin');gap=abs(score[k,1]-float(rr['self']));assert gap<3e-6;checks+=1;maxgap=max(maxgap,gap)
     if rr['secondary']!='NA' and l==0:
      parts=name.split('_');m,s=map(int,parts[0].split('x'));z=readbin(d/'data.bin');x=z[:,1:1+m];y=z[:,1+m:1+m+s];sx=x.max(0);sy=y.max(0);q=secondary(x/sx,y/sy,k,'crs','input',rr['method']);a0=np.r_[readbin(d/'v.bin')[k]*sx,readbin(d/'u.bin')[k]*sy];delta=abs(q@a0-float(rr['secondary']))/max(1,abs(float(rr['secondary'])));assert delta<3e-6,(name,delta);sec_checks+=1;maxsec=max(maxsec,delta)
     if name.startswith('1x1_') and '_crs_' in name:
      delta=abs(readbin(d/'ce.bin')[k,l]-float(rr['ce']));assert delta<3e-6,(name,delta);direct+=1;maxce=max(maxce,delta)
   r.update(passed=True,self_checks=checks,self_max_gap=maxgap,secondary_checks=sec_checks,secondary_max_relative_gap=maxsec,unique_CRS_CE_cells=direct,unique_CRS_CE_max_gap=maxce)
  except Exception as e:r['error']=repr(e)
  rows.append(r);print(json.dumps(r),flush=True)
 (out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');summary=dict(total=len(rows),passed=sum(r['passed'] for r in rows),failed=sum(not r['passed'] for r in rows));(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary);return bool(summary['failed'])
if __name__=='__main__':raise SystemExit(main())
