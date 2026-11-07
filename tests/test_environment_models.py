#!/usr/bin/env python3
"""Environmental model conventions: eliminated-slack SBM oracle, DDF multipliers,
legacy equality-ratio counterexamples, domain/RTS/backend contracts and original R.
Requires licensed Stata 16+, Python with NumPy/SciPy, optionally R/deaR.
"""
from _support import package_files
from pathlib import Path
import argparse,csv,hashlib,itertools,json,re,shutil,subprocess,time
import numpy as np
from scipy.optimize import linprog
from test_core_retries import WRITER,readbin
from test_output_protection import unchanged
from test_fuzzy_models import dual as radial_dual
REPO=Path(__file__).resolve().parents[1]

def solve(c,U=None,b=None,E=None,e=None,bounds=(0,None)):
 f=linprog(c,A_ub=U,b_ub=b,A_eq=E,b_eq=e,bounds=bounds,method='highs-ds',options=dict(presolve=False,primal_feasibility_tolerance=1e-10,dual_feasibility_tolerance=1e-10,time_limit=3.))
 if not f.success:f=linprog(c,A_ub=U,b_ub=b,A_eq=E,b_eq=e,bounds=bounds,method='highs-ipm',options=dict(presolve=False,primal_feasibility_tolerance=1e-10,dual_feasibility_tolerance=1e-10,time_limit=3.))
 assert f.success,f.message
 if U is not None:assert np.max(U@f.x-b,initial=0)<1e-7
 if E is not None:assert np.max(abs(E@f.x-e),initial=0)<1e-7
 return f.fun

def oracle(data,m,s,q,model,rts,direction):
 x=data[:,1:1+m].copy();y=data[:,1+m:1+m+s].copy();b=data[:,1+m+s:1+m+s+q].copy();n=len(x)
 for z in [x,y,b]:sc=z.max(0);sc[sc==0]=1;z/=sc
 ans=[]
 for k in range(n):
  if model=='sbmundesirable':
   # Eliminate all slacks first; use only (z=t*lambda,t), not production's full LP.
   c=np.r_[np.mean(x/x[k],axis=1),0];normal=np.r_[(np.sum(y/y[k],axis=1)-np.sum(b/b[k],axis=1))/(s+q),2*q/(s+q)]
   U=np.vstack([np.c_[x.T,-x[k]],np.c_[-y.T,y[k]],np.c_[b.T,-b[k]]]);rhs=np.zeros(m+s+q);E=normal[None,:];e=np.ones(1);rt=np.r_[np.ones(n),-1]
   if rts=='vrs':E=np.vstack([E,rt]);e=np.r_[e,0]
   if rts in ['drs','irs']:U=np.vstack([U,rt*(1 if rts=='drs' else -1)]);rhs=np.r_[rhs,0]
   ans.append([solve(c,U,rhs,E,e)])
  elif model=='undesirable':
   gx=x[k] if direction in ['x','xyb'] else np.zeros(m);gy=y[k] if direction!='x' else np.zeros(s);gb=b[k] if direction!='x' else np.zeros(q)
   # Multiplier dual of max beta, with free bad-output multipliers.
   c=np.r_[x[k],-y[k],b[k]];U=-np.c_[x,-y,b];directionrow=-np.r_[gx,gy,gb];bounds=[(0,None)]*(m+s)+[(None,None)]*q
   if rts!='crs':c=np.r_[c,1];U=np.c_[U,-np.ones(n)];directionrow=np.r_[directionrow,0];bounds+=[(None,None) if rts=='vrs' else (0,None) if rts=='drs' else (None,0)]
   ans.append([solve(c,np.vstack([U,directionrow]),np.r_[np.zeros(n),-1],bounds=bounds)])
  else:
   # Equality LP dual, including the upper bound theta<=1.
   A=np.vstack([np.c_[x.T,-x[k]],np.c_[y.T,np.zeros(s)],np.r_[np.ones(n),0]]);rhs=np.r_[np.zeros(m),y[k],1]
   U=np.c_[A.T,np.r_[np.zeros(n),-1]];c=np.r_[-rhs,1];ce=-solve(c,U,np.r_[np.zeros(n),1],bounds=[(None,None)]*len(rhs)+[(0,None)])
   cr=radial_dual(x,y,k,'crs','input');vr=radial_dual(x,y,k,'vrs','input');ans.append([vr/ce,cr,vr,ce,vr/ce,cr/vr])
 return np.array(ans)

def fixture(m,s,q,variant):
 rng=np.random.default_rng(7721+m);z=rng.uniform(.7,8,(5,m+s+q));z[-1,:m]=z[0,:m]*1.6;z[-1,m:m+s]=z[0,m:m+s]*.65;z[-1,m+s:]=z[0,m+s:]*1.7
 if variant=='units':z*=10.**np.linspace(-10,10,m+s+q)
 if variant in ['duplicate','near']:z[1]=z[0]*(1+1e-8 if variant=='near' else 1)
 if variant=='zero':
  if m>1:z[1,0]=0
  if s>1:z[2,m]=0
  if q>1:z[3,m+s]=0
 a=np.c_[np.arange(1,6),z,np.ones(5)]
 if variant=='sample':a[1,-1]=0
 if variant=='permuted':a=a[[4,2,0,3,1]]
 cols=['id']+[f'x{i+1}' for i in range(m)]+[f'y{i+1}' for i in range(s)]+[f'b{i+1}' for i in range(q)]+['selected']
 return a,cols

R_SCRIPT=r'''
Sys.setenv(RGL_USE_NULL='TRUE');options(rgl.useNULL=TRUE);suppressPackageStartupMessages(library(deaR));a<-commandArgs(TRUE)
d<-read.csv(a[1]);m<-as.integer(a[2]);s<-as.integer(a[3]);q<-as.integer(a[4]);d<-d[d$selected==1,];z<-as.matrix(d[,2:(1+m+s+q)]);z<-sweep(z,2,apply(z,2,max),'/');df<-data.frame(id=d$id,z)
f<-make_deadata(df,inputs=2:(m+1),outputs=(m+2):(m+s+q+1),ud_outputs=(s+1):(s+q));ans<-list()
for(rt in c('crs','vrs','nirs','ndrs')) {
 fit<-model_sbmeff(f,orientation='no',rts=rt);e<-vapply(fit$DMU,function(z)z$efficiency,numeric(1));ans[[length(ans)+1]]<-data.frame(model='sbmundesirable',rts=rt,direction='none',id=d$id,score=e)
 for(di in c('yb','xyb','x')) {
  dx<-if(di=='yb') matrix(0,m,nrow(z)) else t(z[,1:m,drop=FALSE]);dy<-if(di=='x') matrix(0,s+q,nrow(z)) else t(z[,(m+1):(m+s+q),drop=FALSE])
  fit<-model_dir(f,dir_input=dx,dir_output=dy,rts=rt,maxslack=FALSE);e<-vapply(fit$DMU,function(z)z$beta,numeric(1));ans[[length(ans)+1]]<-data.frame(model='undesirable',rts=rt,direction=di,id=d$id,score=e)
 }
}
write.csv(do.call(rbind,ans),a[5],row.names=FALSE);cat('deaR',as.character(packageVersion('deaR')),'\n')
'''

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',required=True,type=Path);ap.add_argument('--only');ap.add_argument('--baseline-ado',type=Path);ap.add_argument('--rscript');a=ap.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir();h={};rows=[]
 for n in package_files(REPO):raw=(REPO/n).read_bytes();(pkg/Path(n).name).write_bytes(raw);h[n]=hashlib.sha256(raw).hexdigest()
 (out/'manifest.json').write_text(json.dumps(dict(source_sha256=h),indent=2)+'\n')
 def run(name,body,source=pkg,check=None):
  if a.only and not re.search(a.only,name):return
  d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset type double\nset processors 1\nset linesize 180\n'
  for f in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/f).mkdir();head+=f'sysdir set {f} "{d/f}"\n'
  (d/'test.do').write_text(head+f'adopath ++ "{source}"\n'+WRITER+body+'\ndi "TEST_PASS"\n');r=dict(test=name,passed=False);t=time.monotonic()
  try:
   proc=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=180);log=(d/'test.log').read_text(errors='replace');assert proc.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2500:]
   if check:r.update(check(d))
   r['passed']=True
  except Exception as e:r['error']=repr(e)
  r['seconds']=round(time.monotonic()-t,3);rows.append(r);(out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(r),flush=True)
 specs=[]
 for m,s,q in [(1,1,1),(2,2,2)]:
  xs=' '.join(f'x{i+1}' for i in range(m));ys=' '.join(f'y{i+1}' for i in range(s));bs=' '.join(f'b{i+1}' for i in range(q))
  for variant in ['base','units','duplicate','near','sample','permuted','zero']:
   data,cols=fixture(m,s,q,variant);key=f'{m}{s}{q}_{variant}';f=out/(key+'.csv');np.savetxt(f,data,delimiter=',',header=','.join(cols),comments='',fmt='%.17g');specs.append((key,m,s,q,f))
   for model in ['congestion','sbmundesirable','undesirable']:
    if model=='sbmundesirable' and variant=='zero':continue
    for rts,direction,engine in itertools.product(['vrs'] if model=='congestion' else ['crs','vrs','drs','irs'],['yb','xyb','x'] if model=='undesirable' else ['none'],['simplex','ipm','auto']):
     opts=f'model({model}) rts({rts}) engine({engine}) dmu(id) gen(res)'+(f' bad({bs})' if model!='congestion' else '')+(f' direction({direction})' if model=='undesirable' else '')
     body=f'import delimited using "{f}", clear asdouble\ndea {xs}={ys} if selected, {opts}\nassert "`r(lp_backend)\'"=="stata_linearprogram"\n'
     want={'congestion':'input','sbmundesirable':'non-oriented','undesirable':'directional'}[model];body+=f'assert "`r(orientation)\'"=="{want}"\n'
     vs='res res_crs res_vrs res_cong_eff res_congestion res_se' if model=='congestion' else 'res';body+=f'matrix meta=(r(N),r(n_failed))\nmata: savebin("meta.bin",st_matrix("meta"))\nmata: savebin("scores.bin",st_data(.,tokens("{vs}"),"selected"))\nmata: savebin("data.bin",st_data(.,tokens("'+ ' '.join(cols[:-1])+'"),"selected"))\n'
     body+='foreach v of varlist '+vs+' {\nassert missing(`v\') if !selected\n}\n'
     def check(d,m=m,s=s,q=q,model=model,rts=rts,direction=direction):
      actual=readbin(d/'scores.bin');ref=oracle(readbin(d/'data.bin'),m,s,q,model,rts,direction);ok=np.isfinite(actual);gap=float(np.max(abs(actual[ok]-ref[ok]),initial=0));assert gap<2e-6,(gap,d.name);meta=readbin(d/'meta.bin').ravel();assert meta[1]==sum((~ok).any(1));return dict(independent_LPs=len(ref)*(3 if model=='congestion' else 1),missing=int((~ok).sum()),max_gap=gap)
     run('_'.join([key,model,rts,direction,engine]),body,check=check)
 # Constant-output false positive and a nonradial false negative for the old label.
 for case,data,expected in [('constant_output','1 1 1 1\n2 2 3 1',.5),('nonradial','1 1 1 2\n2 2 1 1',1.)]:
  body='input double(id x1 x2 y)\n'+data+'\nend\ndea_congestion x1 x2=y,dmu(id) gen(res)\n'+f'assert abs(res-{expected})<1e-6 in 2\n'
  run('interpretation_'+case,body)
  if a.baseline_ado:run('old_interpretation_'+case,body,source=a.baseline_ado.resolve())
 D='input double(id x y b)\n1 1 1 1\n2 3 2 4\n3 2 3 2\nend\n'
 for model in ['congestion','sbmundesirable']:
  for field in ['x','y']+(['b'] if model=='sbmundesirable' else []):
   for value in [0,-1]:
    cmd=f'dea x=y,model({model}) gen(res) replacevars'+(' bad(b)' if model=='sbmundesirable' else '')
    run(f'domain_{model}_{field}_{value}',D+f'replace {field}={value} in 2\ngen str12 res="old result"\n'+unchanged(cmd,459))
 for model in ['congestion','sbmundesirable','undesirable']:
  opt=f'model({model})'+(' bad(b)' if model!='congestion' else '')
  run('orientation_stage_'+model,D+f'dea x=y,{opt} ort(input) stage(1) dmu(id) gen(a)\ndea x=y,{opt} ort(output) stage(2) dmu(id) gen(c)\nassert abs(a-c)<1e-9 | missing(a,c)\n')
 # SBM's requested RTS must be honored, including aliases and wrappers.
 for rt in ['crs','vrs','drs','irs','nirs','ndrs','ccr','bcc']:
  run('wrapper_rts_'+rt,D+f'dea x=y,model(sbmundesirable) bad(b) rts({rt}) dmu(id) gen(a)\ndea_sbmundesirable x=y,bad(b) rts({rt}) dmu(id) gen(c)\nassert a==c\n')
 if a.baseline_ado:
  run('old_ignored_sbm_rts',D+'dea x=y,model(sbmundesirable) bad(b) rts(crs) dmu(id) gen(a)\ndea x=y,model(sbmundesirable) bad(b) rts(drs) dmu(id) gen(c)\nassert a==c\n',source=a.baseline_ado.resolve())
 if a.rscript:
  (out/'compare.R').write_text(R_SCRIPT)
  def rcheck(d):
   count=missing=0;gap=0
   for key,m,s,q,f in specs:
    if not key.endswith('_base'):continue
    dest=d/(key+'.csv');p=subprocess.run([a.rscript,str(out/'compare.R'),str(f),str(m),str(s),str(q),str(dest)],capture_output=True,text=True,timeout=180);(d/(key+'.console')).write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr
    with dest.open() as fh:rr=list(csv.DictReader(fh))
    for model,rts,di in itertools.product(['sbmundesirable','undesirable'],['crs','vrs','drs','irs'],['none','yb','xyb','x']):
     if (model=='sbmundesirable')!=(di=='none'):continue
     path=out/'_'.join([key,model,rts,di,'ipm'])/'scores.bin'
     if not path.exists():continue
     tag={'drs':'nirs','irs':'ndrs'}.get(rts,rts);ref=np.array([float(z['score']) for z in rr if z['model']==model and z['rts']==tag and z['direction']==di]);actual=readbin(path).ravel();assert ref.shape==actual.shape and np.isfinite(ref).all();ok=np.isfinite(actual);g=float(np.max(abs(actual[ok]-ref[ok]),initial=0));assert g<2e-6,(model,rts,di,g);gap=max(gap,g);count+=int(ok.sum());missing+=int((~ok).sum())
   assert count>0;return dict(R_scores=count,R_missing=missing,R_max_gap=gap)
  run('R_deaR','',check=rcheck)
 summary=dict(total=len(rows),passed=sum(x['passed'] for x in rows),failed=sum(not x['passed'] for x in rows));(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary);return bool(summary['failed'])
if __name__=='__main__':raise SystemExit(main())
