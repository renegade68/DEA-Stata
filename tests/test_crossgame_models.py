#!/usr/bin/env python3
"""Crossgame: independent primary/secondary LPs, every complete game iteration,
public matrices, failed-iteration isolation and paper-example diagnostics."""
from _support import package_files
from pathlib import Path
import argparse,csv,hashlib,itertools,json,os,re,shutil,subprocess,time
import numpy as np
from scipy.optimize import linprog
from test_core_retries import WRITER,readbin
from test_crosseff_models import fixture,problem,fit
from test_output_protection import unchanged
REPO=Path(__file__).resolve().parents[1]
PAPER=np.array([[7,7,7,4,4],[5,9,7,7,7],[4,6,5,5,7],[5,9,8,6,2],[6,8,5,3,6]],float)
PAPER_GAME=np.array([.6384,.9766,1,.7988,.6670])
MATS=['self_eff','weights_input','weights_output','ce_initial']
for ab,word in [('agg','aggressive'),('ben','benevolent')]:MATS += ['ce_'+ab,'weights_'+ab+'_input','weights_'+ab+'_output','status_'+ab,'count_'+ab]
MATS+=['ce_nash','nash_history','nash_trace']

def game(x,y,d,j,theta):
 c,e,b,a,h,bd=problem(x,y,j,'crs','input');a=np.vstack([a,np.r_[theta*x[d],-y[d]]]);h=np.r_[h,0.]
 return linprog(c,A_eq=e,b_eq=b,A_ub=a,b_ub=h,bounds=bd,method='highs',options=dict(primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9))

def verify(d,m,s,strategy,tol,limit):
 z=readbin(d/'data.bin');x=z[:,1:1+m];y=z[:,1+m:];n=len(z);sx=x.max(0);sy=y.max(0);sx[sx==0]=1;sy[sy==0]=1;xx=x/sx;yy=y/sy
 A={p.stem:readbin(p) for p in d.glob('*.bin')};meta=A['meta'].ravel();self=A['self_eff'].ravel();scores=A['scores'];self_gap=sec_gap=game_gap=0.;game_lps=secondary_lps=feasible_fallback=unbounded=0
 for k in range(n):
  p=problem(xx,yy,k,'crs','input');f=fit(p);assert f.success
  assert np.isfinite(self[k]);self_gap=max(self_gap,abs(self[k]+f.fun));assert abs(self[k]+f.fun)<3e-6
  for ab in ['','agg_','ben_']:
   vk='weights_'+ab+'input'
   if vk not in A:continue
   w=np.r_[A[vk][k]*sx,A['weights_'+ab+'output'][k]*sy];c,e,b,u,h,bd=p
   assert np.max(u@w)<3e-6 and np.max(abs(e@w-b))<3e-6 and min(w)>-3e-6
   assert abs(c@w+self[k])<3e-6
   if ab:
    st=A['status_'+ab[:-1]][k,0];q=np.r_[np.zeros(m),yy.sum(0)-yy[k]]
    if ab=='ben_':q=-q
    g=fit(p,q,lock=-self[k]);secondary_lps+=1
    if st==2:
     assert g.success;gap=abs(q@w-g.fun)/max(1,abs(g.fun));sec_gap=max(sec_gap,gap);assert gap<3e-6
    else:
     assert st==3;feasible_fallback+=int(g.success);unbounded+=int(g.status==3)
     assert np.array_equal(A[vk][k],A['weights_input'][k]) and np.array_equal(A['weights_'+ab+'output'][k],A['weights_output'][k])
 for ab in ['initial','agg','ben']:
  ck='ce_'+ab
  if ck not in A:continue
  tag='' if ab=='initial' else ab+'_';num=A['weights_'+tag+'output']@y.T;den=A['weights_'+tag+'input']@x.T
  with np.errstate(divide='ignore',invalid='ignore'):q=num/den
  mask=np.isfinite(q)&np.isfinite(den)&np.isfinite(num)&(den>0)&(num>=0)&(q<=1+tol);expected=np.where(mask,q,np.nan)
  assert np.array_equal(np.isfinite(A[ck]),np.isfinite(expected)),(ab,A[ck],expected)
  assert np.allclose(A[ck],expected,equal_nan=True,atol=3e-6,rtol=1e-8)
  if ab=='initial':assert np.allclose(A[ck].diagonal(),self,atol=3e-6)
  else:
   counts=[];means=[]
   for j in range(n):
    v=A[ck][np.arange(n)!=j,j];v=v[np.isfinite(v)];counts.append(len(v));means.append(v.mean() if len(v) else np.nan)
   assert np.array_equal(A['count_'+ab].ravel(),counts)
   si=0 if ab=='agg' else (1 if strategy=='all' else 0)
   assert np.allclose(scores[:,si],means,equal_nan=True,atol=1e-10)
   mi=1 if ab=='agg' else 2;mpi=9 if ab=='agg' else 10
   assert meta[mi]==np.sum(A['status_'+ab]==3) and meta[mpi]==n*(n-1)-sum(counts)
 if strategy in ['all','nash']:
  ns,ni,na,nf,change,missing=meta[3:9];hist=A['nash_history'];trace=A['nash_trace'];ce=A['ce_nash'];assert ns in [0,1,2,3]
  assert missing==np.sum(~np.isfinite(A['ce_initial']))
  if ns==2:
   assert missing>0 and ni==na==0 and np.isnan(scores[:,-1]).all() and np.isnan(ce).all()
  else:
   assert missing==0 and len(hist)==ni+1 and len(trace)==na+1
   assert np.allclose(hist[0],A['ce_initial'].mean(0),atol=1e-12)
   last=None
   for t in range(1,len(hist)):
    new=np.empty((n,n))
    for k,j in itertools.product(range(n),repeat=2):
     f=game(xx,yy,k,j,hist[t-1,k]);assert f.success,(t,k,j,f.message);new[k,j]=-f.fun;game_lps+=1
    gap=np.max(abs(new.mean(0)-hist[t]));game_gap=max(game_gap,gap);assert gap<3e-6,(t,gap)
    delta=np.max(abs(hist[t]-hist[t-1]));assert abs(trace[t,1]-delta)<1e-12 and trace[t,2]==0 and trace[t,3]==1
    if t<len(hist)-1:assert delta>=tol
    last=new
   if ni>0:assert np.allclose(ce,last,atol=3e-6);assert abs(change-np.max(abs(hist[-1]-hist[-2])))<1e-12
   else:assert np.array_equal(ce,A['ce_initial'])
   assert np.allclose(scores[:,-1],hist[-1],atol=1e-12) and np.allclose(ce.mean(0),hist[-1],atol=1e-12)
   if ns==0:assert change<tol and ni==na and nf==0
   if ns==1:assert ni==na==limit and change>=tol and nf==0
   if ns==3:assert na==ni+1 and nf>0 and trace[-1,3]==0
 return dict(primary_LPs=n,secondary_LPs=secondary_lps,game_LPs=game_lps,self_max_gap=self_gap,secondary_max_relative_gap=sec_gap,game_mean_max_gap=game_gap,secondary_feasible_fallback=feasible_fallback,secondary_unbounded=unbounded,nash_status=int(meta[3]) if strategy in ['all','nash'] else None,nash_iterations=int(meta[4]) if strategy in ['all','nash'] else 0,matrix_cells=sum(A[k].size for k in ['ce_initial','ce_agg','ce_ben','ce_nash'] if k in A))

R_SCRIPT=r'''
library(lpSolve)
a=commandArgs(TRUE);z=as.matrix(read.csv(a[1]));x=z[,2:4];y=z[,5:6];n=nrow(x);m=ncol(x);s=ncol(y)
solve1=function(j,d=NULL,theta=NULL) {
 obj=c(rep(0,m),y[j,]);A=rbind(c(x[j,],rep(0,s)),cbind(-x,y));dir=c('=',rep('<=',n));rhs=c(1,rep(0,n))
 if(!is.null(d)) {A=rbind(A,c(theta*x[d,],-y[d,]));dir=c(dir,'<=');rhs=c(rhs,0)}
 f=lp('max',obj,A,dir,rhs);stopifnot(f$status==0);f
}
f=lapply(1:n,solve1);v=t(sapply(f,function(q)q$solution[1:m]));u=t(sapply(f,function(q)q$solution[(m+1):(m+s)]));theta=colMeans((u%*%t(y))/(v%*%t(x)))
for(iter in 1:200) {ce=matrix(0,n,n);for(d in 1:n)for(j in 1:n)ce[d,j]=solve1(j,d,theta[d])$objval;nextscore=colMeans(ce);delta=max(abs(theta-nextscore));theta=nextscore;if(delta<1e-10)break}
write.csv(data.frame(self=sapply(f,function(q)q$objval),game=theta),a[2],row.names=FALSE)
writeLines(c(paste('iterations',iter,'max_change',delta),capture.output(sessionInfo())),a[3])
'''

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--source-ado',type=Path);ap.add_argument('--only');ap.add_argument('--rscript');a=ap.parse_args()
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir();h={}
 paths=list(a.source_ado.iterdir()) if a.source_ado else [REPO/n for n in package_files(REPO)]
 for p in paths:
  if p.suffix in ['.ado','.do','.sthlp']:raw=p.read_bytes();(pkg/p.name).write_bytes(raw);h[p.name]=hashlib.sha256(raw).hexdigest()
 (out/'manifest.json').write_text(json.dumps(dict(source_sha256=h,stata=a.stata),indent=2)+'\n');rows=[]
 def run(name,body,check=None,source=pkg):
  if a.only and not re.search(a.only,name):return
  d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset processors 1\nset type double\nset linesize 160\n'
  for f in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/f).mkdir();head+=f'sysdir set {f} "{d/f}"\n'
  head+=f'adopath ++ "{source}"\n'+WRITER;(d/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n');r=dict(test=name,passed=False);t=time.monotonic()
  try:
   p=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=180);log=(d/'test.log').read_text(errors='replace');assert p.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2300:]
   if check:r.update(check(d))
   r['passed']=True
  except Exception as e:r['error']=repr(e)
  r['seconds']=round(time.monotonic()-t,3);rows.append(r);(out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(r),flush=True)
 specs=[]
 for m,s,var in itertools.product([1,2],[1,2],['base','units','tiny','zero','zerocol','duplicate','sample','permuted','near']):
  z,cols=fixture(m,s,var);specs.append((f'{m}x{s}_{var}',m,s,z,cols))
 specs.append(('paper',3,2,np.c_[np.arange(1,6),PAPER,np.ones(5)],['id','x1','x2','x3','y1','y2','selected']))
 for key,m,s,z,cols in specs:
  file=out/('data_'+key+'.csv')
  with file.open('w') as f:w=csv.writer(f);w.writerow(cols);w.writerows(['' if np.isnan(v) else format(v,'.17g') for v in row] for row in z)
  xs=' '.join(f'x{i+1}' for i in range(m));ys=' '.join(f'y{i+1}' for i in range(s))
  params=list(itertools.product(['all','aggressive','benevolent','nash'],['1e-6','1e-9'],[25]))
  if key=='paper':params += [('nash','1e-3',100),('nash','1e-10',200),('nash','1e-10',1)]
  for ci,(strategy,tol,limit) in enumerate(params):
   cmd=f'dea_crossgame {xs}={ys} if selected,' if ci%2 else f'dea_crossgame if selected,inputs({xs}) outputs({ys})'
   cmd+=f' strategy({strategy}) dmu(id) gen(res) tol({tol}) maxiter({limit}) matrix detail'
   body=f'import delimited using "{file}",clear varnames(1) asdouble\n'+cmd+'\n'
   # Unrequested scalars are explicitly missing when constructing this test record.
   body+='matrix meta=J(1,11,.)\n'
   for k,key0 in enumerate(['n_failed_primary','n_failed_agg','n_failed_ben','nash_status','nash_iterations','nash_attempts','nash_failed_pairs','nash_last_change','nash_initial_missing','n_missing_agg','n_missing_ben']):
    if (key0.endswith('agg') and strategy not in ['all','aggressive']) or (key0.endswith('ben') and strategy not in ['all','benevolent']) or (key0.startswith('nash_') and strategy not in ['all','nash']):continue
    body+=f'matrix meta[1,{k+1}]=r({key0})\n'
   body+='mata: savebin("meta.bin",st_matrix("meta"))\n'
   for name in MATS:
    if ('agg' in name and strategy not in ['all','aggressive']) or ('ben' in name and strategy not in ['all','benevolent']) or ('nash' in name and strategy not in ['all','nash']):continue
    body+=f'mata: savebin("{name}.bin",st_matrix("r({name})"))\n'
   outs='res_agg res_ben res_nash' if strategy=='all' else 'res_'+dict(aggressive='agg',benevolent='ben',nash='nash')[strategy]
   body+=f'foreach v of varlist {outs} {{\nassert missing(`v\') if !selected\n}}\nkeep if selected\nmata: savebin("data.bin",st_data(.,tokens("id {xs} {ys}")))\nmata: savebin("scores.bin",st_data(.,tokens("{outs}")))\n'
   run(f'{key}_{strategy}_{tol}_{limit}',body,lambda d,m=m,s=s,strategy=strategy,tol=tol,limit=limit:verify(d,m,s,strategy,float(tol),limit))
 data='input double(id x y)\n1 1 1\n2 2 1\n3 3 2\nend\n'
 for opt in ['tol(0)','tol(-1)','tol(.)','tol(.a)','maxiter(0)','maxiter(-1)','maxiter(.)','maxiter(1.5)']:
  run('invalid_'+re.sub(r'\W','_',opt),data+'gen str12 res_nash="keep"\n'+unchanged('dea_crossgame x=y,gen(res) strategy(nash) replacevars '+opt,198))
 for key,setup in [('negative_x','replace x=-1 in 1'),('negative_y','replace y=-1 in 1'),('zero_x','replace x=0 in 1'),('zero_y','replace y=0 in 1')]:run('domain_'+key,data+setup+'\ngen str12 res_nash="keep"\n'+unchanged('dea_crossgame x=y,gen(res) strategy(nash) replacevars',459))
 run('string_quantity',data+'tostring x,replace\ngen str12 res_nash="keep"\n'+unchanged('dea_crossgame x=y,gen(res) strategy(nash) replacevars',109))
 run('missing_ids',data+'replace id=. in 1\ndea_crossgame x=y,dmu(id) gen(res) matrix detail\nassert r(n)==2\nassert missing(res_nash) in 1\n')
 run('singleton',data+'keep in 1\ndea_crossgame x=y,strategy(nash) gen(res) matrix\nassert r(nash_status)==0 & r(nash_iterations)==1\nassert abs(res_nash-1)<1e-10\n')
 run('small_denominator','input double(x y)\n1 1\n1e8 1e8\nend\ndea_crossgame x=y,gen(res) matrix\nassert r(nash_status)==0 & r(nash_initial_missing)==0\nassert abs(res_nash-1)<1e-6\n')
 run('workspace_matrices',data+'matrix sentinel=(991,992)\nlocal before : all matrices\ndea_crossgame x=y,gen(res) matrix\nlocal after : all matrices\nassert "`before\'"=="`after\'"\nassert rowsof(r(ce_nash))==3\ndea_crossgame x=y,gen(res) replacevars\nlocal after : all matrices\nassert "`before\'"=="`after\'"\n')
 run('reload_state',data+f'run "{pkg}/dea_crossgame.ado"\nrun "{pkg}/dea_crossgame.ado"\ndea_crossgame x=y,gen(res) matrix\nassert rowsof(r(ce_agg))==3\ndea_crossgame x=y,gen(res) strategy(nash) replacevars\ncapture matrix list r(ce_agg)\nassert _rc==111\ncapture matrix list r(nash_history)\nassert _rc==111\n')
 if a.rscript and (out/'paper_nash_1e-10_200/scores.bin').exists():
  (out/'compare.R').write_text(R_SCRIPT);p=subprocess.run([a.rscript,str(out/'compare.R'),str(out/'data_paper.csv'),str(out/'r_scores.csv'),str(out/'R_session.txt')],capture_output=True,env=dict(os.environ,RGL_USE_NULL='TRUE'),timeout=120);(out/'R.log').write_bytes(p.stdout+p.stderr);r=dict(test='paper_and_R',passed=False)
  try:
   assert p.returncode==0,p.stderr.decode();rr=np.genfromtxt(out/'r_scores.csv',delimiter=',',names=True);d=out/'paper_nash_1e-10_200';native=readbin(d/'scores.bin').ravel();self=readbin(d/'self_eff.bin').ravel();gap=float(max(abs(native-rr['game'])));assert gap<3e-6;assert max(abs(self-rr['self']))<3e-6
   r.update(passed=True,R_game_max_gap=gap,paper_rounded_column_max_gap=float(max(abs(native-PAPER_GAME))),native_game=native.tolist(),R_game=rr['game'].tolist(),paper_stopping_epsilon=.001)
  except Exception as e:r['error']=repr(e)
  rows.append(r);print(json.dumps(r),flush=True)
 (out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');summary=dict(total=len(rows),passed=sum(r['passed'] for r in rows),failed=sum(not r['passed'] for r in rows));(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary);return bool(summary['failed'])
if __name__=='__main__':raise SystemExit(main())
