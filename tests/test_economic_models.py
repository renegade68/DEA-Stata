#!/usr/bin/env python3
"""Economic DEA: independent lambda-only LPs, Benchmarking, and price validation.

Requires NumPy/SciPy and Stata 16+; --rscript compares cost/revenue and radial
decompositions with unmodified Benchmarking. Profit uses the documented
restricted-improvement model, not Benchmarking::profit.opt.
"""
from _support import package_files
from pathlib import Path
import argparse, csv, hashlib, itertools, json, re, subprocess, time
import numpy as np
import scipy
from scipy.optimize import linprog
from verify_bootstrap_lp import solve as radial

REPO=Path(__file__).resolve().parents[1]
RTS=['crs','vrs','drs','irs']; ENGINES=['simplex','ipm','auto']
MODELS=['cost','revenue','profit','alloc_in','alloc_out']
ALIASES=dict(crs='ccr',vrs='bcc',drs='nirs',irs='ndrs')

def fixture(m,s,variant):
 x=np.array([[2,12],[2,8],[5,5],[10,4],[10,6],[3,13]],float)[:,:m]
 y=np.array([[1,4],[2,1],[2,3],[4,2],[3,3],[1.5,2]],float)[:,:s]
 rng=np.random.default_rng(93412+10*m+s)
 w=rng.uniform(.2,.6,(6,m));p=rng.uniform(10,15,(6,s))
 if variant=='random':x=rng.uniform(.5,9,x.shape);y=rng.uniform(.5,9,y.shape);p*=4
 if variant=='zero':
  if m>1:x[0,0]=0
  if s>1:y[1,0]=0
 if variant=='units':
  sx=np.array([.01,100])[:m];sy=np.array([10,.1])[:s]
  x*=sx;w/=sx;y*=sy;p/=sy
 if variant=='prices':
  f=np.array([.1,2,10,.5,5,1])[:,None];w*=f;p*=f
 a=np.c_[np.arange(1,7),x,y,w,p]
 cols=['id']+[f'{v}{j+1}' for v,n in [('x',m),('y',s),('w',m),('p',s)] for j in range(n)]
 if variant=='permuted':a=a[[5,2,0,4,1,3]]
 if variant=='duplicate':a=np.vstack([a,a[1]]);a[-1,0]=7
 return a,cols

def oracle(x,y,w,p,rts):
 # Positive prices eliminate x*,y*: cost/revenue LPs have lambda variables only.
 n=len(x);xx=x/np.maximum(x.max(axis=0),1e-300);yy=y/np.maximum(y.max(axis=0),1e-300)
 ae=np.ones((1,n)) if rts=='vrs' else None;be=np.ones(1) if rts=='vrs' else None;ans=[]
 for k in range(n):
  def opt(c,a,b):
   if rts in ['drs','irs']:
    sign=1 if rts=='drs' else -1;a=np.vstack([a,sign*np.ones(n)]);b=np.r_[b,sign]
   scale=max(abs(c).max(),1e-300)
   f=linprog(c/scale,A_ub=a,b_ub=b,A_eq=ae,b_eq=be,bounds=(0,None),method='highs',options=dict(primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9))
   assert f.success,f.message
   assert np.max(a@f.x-b)<1e-7
   if ae is not None:assert abs(f.x.sum()-1)<1e-7
   return f.fun*scale
  ce=opt(x@w[k],-yy.T,-yy[k])/(x[k]@w[k])
  rev=(y[k]@p[k])/-opt(-y@p[k],xx.T,xx[k])
  pr=(y[k]@p[k]-x[k]@w[k])/-opt(x@w[k]-y@p[k],np.vstack([xx.T,-yy.T]),np.r_[xx[k],-yy[k]])
  ti=radial(xx[k],yy[k],xx,yy,'input',rts);to=radial(xx[k],yy[k],xx,yy,'output',rts)
  ans.append([ce,rev,pr,ce/ti,rev/to])
 return np.array(ans)

def read(path):
 with path.open() as f:rows=list(csv.DictReader(f))
 return {k:np.array([float(r[k]) if r[k] not in ['', '.'] and not re.fullmatch(r'\.[a-z]',r[k]) else np.nan for r in rows]) for k in rows[0]}

R_SCRIPT=r'''
suppressPackageStartupMessages(library(Benchmarking))
args=commandArgs(TRUE);jobs=read.csv(args[1],stringsAsFactors=FALSE);ans=list();i=0
for(j in seq_len(nrow(jobs))) {
 z=jobs[j,];d=read.csv(z$file);d=d[d$selected==1,]
 mat=function(cols) { a=as.matrix(d[cols]);storage.mode(a)='double';a }
 x=mat(paste0('x',1:z$m));y=mat(paste0('y',1:z$s))
 w=mat(paste0('w',1:z$m));p=mat(paste0('p',1:z$s))
 for(rts in c('crs','vrs','drs','irs')) {
  co=cost.opt(x,y,w,RTS=rts);ro=revenue.opt(x,y,p,RTS=rts)
  ti=dea(x,y,RTS=rts,ORIENTATION='in',SLACK=FALSE)$eff
  phi=dea(x,y,RTS=rts,ORIENTATION='out',SLACK=FALSE)$eff
  ce=co$cost/rowSums(x*w);re=rowSums(y*p)/ro$rev
  i=i+1;ans[[i]]=data.frame(test=z$test,rts=rts,id=d$id,cost=ce,revenue=re,alloc_in=ce/ti,alloc_out=re*phi)
 }
}
write.csv(do.call(rbind,ans),args[2],row.names=FALSE)
writeLines(c(paste('Benchmarking',packageVersion('Benchmarking')),capture.output(sessionInfo())),args[3])
for(fun in c('cost.opt','revenue.opt','profit.opt','dea'))
 writeLines(deparse(get(fun,asNamespace('Benchmarking'))),file.path(dirname(args[2]),paste0('R_original_',fun,'.txt')))
'''

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True)
 ap.add_argument('--baseline-ado',type=Path);ap.add_argument('--rscript');ap.add_argument('--only');a=ap.parse_args()
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir();hashes={}
 for name in package_files(REPO):
  raw=(REPO/name).read_bytes();(pkg/Path(name).name).write_bytes(raw);hashes[name]=hashlib.sha256(raw).hexdigest()
 (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,stata=a.stata,numpy=np.__version__,scipy=scipy.__version__,baseline=str(a.baseline_ado)),indent=2)+'\n')
 rows=[];jobs=[];expectations={};variants=['base','permuted','duplicate','units','prices','sample','zero','random']
 def chosen(name):return not a.only or re.search(a.only,name)
 def run(name,body,check=None,source=pkg):
  if not chosen(name):return
  d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset processors 1\nset type double\nset linesize 160\n'
  for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/key).mkdir();head+=f'sysdir set {key} "{d/key}"\n'
  head+=f'adopath ++ "{source}"\n';(d/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n')
  row=dict(test=name,passed=False);t=time.monotonic()
  try:
   p=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=240)
   log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else p.stderr.decode(errors='replace')
   assert p.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2500:]
   if check:row.update(check(d))
   row['passed']=True
  except Exception as ex:row['error']=repr(ex)
  row['seconds']=round(time.monotonic()-t,3);rows.append(row);(out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(row),flush=True)
 for m,s,variant in itertools.product([1,2],[1,2],variants):
  name=f'valid_{m}x{s}_{variant}'
  if not chosen(name):continue
  data,cols=fixture(m,s,variant);f=out/(name+'_input.csv');np.savetxt(f,data,delimiter=',',header=','.join(cols),comments='',fmt='%.17g')
  xs=' '.join(f'x{j+1}' for j in range(m));ys=' '.join(f'y{j+1}' for j in range(s));ws=' '.join(f'w{j+1}' for j in range(m));ps=' '.join(f'p{j+1}' for j in range(s))
  q='if id!=2 in 1/5' if variant=='sample' else ''
  body=f'import delimited using "{f}",clear varnames(1) asdouble\ngen byte selected=0\nreplace selected=1 {q}\n'
  if variant=='sample':body+='replace w1=.z if id==4\nreplace selected=0 if id==4\n'
  body+='count if selected\nlocal n=r(N)\n'
  for ei,eng in enumerate(ENGINES):
   for ri,rts in enumerate(RTS):
    for mi,model in enumerate(MODELS):
     ori='output' if model in ['revenue','alloc_out'] else 'input';mod='allocative' if model.startswith('alloc_') else model
     spec=ALIASES[rts].upper() if (ei+mi)%2 else rts
     h=f'dea {xs}={ys} {q},' if (ei+mi)%2 else f'dea {ys} {q}, inputs({xs})'
     v=f'e{ei}_{ri}_{mi}';tol='1e-9' if ri%2 else '1e-7'
     body+=f'{h} model({mod}) w({ws}) p({ps}) ort({ori}) rts({spec}) engine({eng}) tol({tol}) stage({1+(ei+mi)%2}) dmu(id) gen({v})\n'
     body+=f'assert r(N)==`n\' & r(n_failed)==0\nassert !missing({v}) & {v}>0 & {v}<1+2e-6 if selected\nassert missing({v}) if !selected\n'
  body+='format _all %24.17g\nexport delimited using "scores.csv",replace datafmt\n'
  def check(d,m=m,s=s,name=name):
   got=read(d/'scores.csv');sel=got['selected']==1
   x,y,w,p=[np.column_stack([got[f'{v}{j+1}'][sel] for j in range(n)]) for v,n in [('x',m),('y',s),('w',m),('p',s)]]
   assert np.all(np.sum(y*p,axis=1)-np.sum(x*w,axis=1)>0)
   maxgap=0.;expected={}
   for ri,rts in enumerate(RTS):
    exp=oracle(x,y,w,p,rts);expected[rts]=exp
    for ei in range(3):
     actual=np.column_stack([got[f'e{ei}_{ri}_{mi}'][sel] for mi in range(5)])
     gap=float(np.max(abs(actual-exp)));assert gap<2e-6,(name,ei,rts,gap,actual,exp);maxgap=max(maxgap,gap)
   expectations[name]=dict(ids=got['id'][sel],values=expected);jobs.append(dict(test=name,file=str(d/'scores.csv'),m=m,s=s))
   return dict(commands=60,independent_LPs=20*len(x),score_checks=60*len(x),max_score_gap=maxgap)
  run(name,body,check)
  if a.baseline_ado and variant in ['base','units']:
   def compare(d,name=name):
    new=read(out/name/'scores.csv');old=read(d/'scores.csv');count=0
    for k in new:
     if k.startswith('e'):assert np.array_equal(new[k],old[k],equal_nan=True),(name,k);count+=int(np.isfinite(new[k]).sum())
    return dict(unchanged_scores=count,commands=60)
   run(name+'_baseline',body,compare,a.baseline_ado.resolve())
 for eng,rts in itertools.product(ENGINES,RTS):
  body='''input double(id x1 x2 y w1 w2)
1 1 4 1 1 2
2 4 1 1 1 2
3 3 3 1 1 2
end
'''
  body+=f'''dea x1 x2=y, model(cost) w(w1 w2) rts({rts}) engine({eng}) gen(ce) dmu(id)
assert abs(ce-cond(id==2,1,2/3))<2e-7
dea x1 x2=y, model(allocative) w(w1 w2) rts({rts}) engine({eng}) gen(ae) dmu(id)
assert abs(ae-cond(id==1,2/3,cond(id==2,1,.8)))<2e-7
clear
input double(id x y1 y2 p1 p2)
1 1 1 4 1 2
2 1 4 1 1 2
3 1 2 2 1 2
end
dea x=y1 y2, model(revenue) p(p1 p2) rts({rts}) engine({eng}) gen(rev) dmu(id)
assert abs(rev-cond(id==1,1,2/3))<2e-7
dea x=y1 y2, model(allocative) ort(output) p(p1 p2) rts({rts}) engine({eng}) gen(ae) dmu(id)
assert abs(ae-cond(id==1,1,cond(id==2,2/3,5/6)))<2e-7
clear
input double(id x y w p)
1 1 2 1 2
2 3 10 1 2
end
dea x=y, model(profit) w(w) p(p) rts({rts}) engine({eng}) gen(pe) dmu(id)
assert abs(pe-cond(id==2,1,{9/17 if rts in ['crs','drs'] else 1}))<2e-7
'''
  run(f'analytic_{eng}_{rts}',body,lambda d:dict(commands=5,analytic_scores=14))
 data,cols=fixture(2,2,'base');f=out/'invalid_input.csv';np.savetxt(f,data,delimiter=',',header=','.join(cols),comments='',fmt='%.17g')
 for eng,ori in itertools.product(ENGINES,['input','output']):
  active='w' if ori=='input' else 'p';other='p' if ori=='input' else 'w'
  cases={'none':'','wrongside':f'{other}({other}1 {other}2)','short':f'{active}({active}1)','long':f'{active}({active}1 {active}2 id)','other_short':f'{active}({active}1 {active}2) {other}({other}1)'}
  for issue,opts in cases.items():
   for autoid in [False,True]:
    tag=f'invalid_{eng}_{ori}_{issue}'+('_autoid' if autoid else '_id');idopt='' if autoid else 'dmu(id)'
    head='dea y1 y2, inputs(x1 x2)' if autoid else 'dea x1 x2=y1 y2,'
    body=f'''import delimited using "{f}",clear varnames(1) asdouble
gen double result=123
label variable result "caller result"
save "protected.dta",replace
datasignature
local sig `"`r(datasignature)'"'
capture noisily {head} model(allocative) ort({ori}) {opts} engine({eng}) {idopt} gen(result) replacevars saving("protected.dta") replace
assert _rc==198
datasignature
assert `"`r(datasignature)'"'==`"`sig'"'
assert result==123
assert `"`:variable label result'"'=="caller result"
capture confirm variable _dmu
assert _rc==111
cf _all using "protected.dta",all
'''
    run(tag,body,lambda d:dict(commands=1))
 # Match identities across observationally equivalent specifications.
 invariants=0;max_invariant_gap=0.
 for m,s in itertools.product([1,2],[1,2]):
  base=out/f'valid_{m}x{s}_base'/'scores.csv'
  if not base.exists():continue
  v=read(base)
  for variant in ['permuted','duplicate','units','prices']:
   path=out/f'valid_{m}x{s}_{variant}'/'scores.csv'
   if not path.exists():continue
   z=read(path)
   for j,identity in enumerate(z['id']):
    i=int(np.where(v['id']==(2 if identity==7 else identity))[0][0])
    for k in v:
     if not k.startswith('e'):continue
     gap=abs(z[k][j]-v[k][i]);assert gap<2e-6,(path,k,gap)
     invariants+=1;max_invariant_gap=max(max_invariant_gap,gap)
 if invariants:
  rows.append(dict(test='unit_price_order_duplicate_invariance',passed=True,score_checks=invariants,max_score_gap=max_invariant_gap))
  (out/'results.json').write_text(json.dumps(rows,indent=2)+'\n')
 if a.rscript and jobs:
  with (out/'r_jobs.csv').open('w') as f:writer=csv.DictWriter(f,list(jobs[0]));writer.writeheader();writer.writerows(jobs)
  (out/'reference.R').write_text(R_SCRIPT)
  with (out/'reference.log').open('w') as f:proc=subprocess.run([a.rscript,str(out/'reference.R'),str(out/'r_jobs.csv'),str(out/'r_scores.csv'),str(out/'R_session.txt')],stdout=f,stderr=subprocess.STDOUT,timeout=600)
  row=dict(test='R_Benchmarking',passed=False)
  try:
   assert proc.returncode==0,(out/'reference.log').read_text()[-2000:]
   count=0;maxgap=0.
   with (out/'r_scores.csv').open() as f:
    for r in csv.DictReader(f):
     exp=expectations[r['test']];i=int(np.where(exp['ids']==float(r['id']))[0][0]);v=exp['values'][r['rts']][i]
     for k,mi in [('cost',0),('revenue',1),('alloc_in',3),('alloc_out',4)]:
      gap=abs(float(r[k])-v[mi]);assert gap<1e-7,(r,gap);count+=1;maxgap=max(maxgap,gap)
   row.update(passed=True,reference_scores=count,max_score_gap=maxgap)
  except Exception as ex:row['error']=repr(ex)
  rows.append(row);(out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(row),flush=True)
 raise SystemExit(any(not r['passed'] for r in rows))

if __name__=='__main__':main()
