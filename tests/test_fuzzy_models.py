#!/usr/bin/env python3
"""Triangular alpha-cut model contracts, independent multiplier duals, interval
extrema, complete-cohort summaries, and optional original deaR comparison.
Requires Stata 16+, Python 3.10+, NumPy/SciPy. No installed source is modified.
"""
from _support import package_files
from pathlib import Path
import argparse,csv,hashlib,itertools,json,re,shutil,subprocess,time
import numpy as np
from scipy.optimize import linprog
from test_core_retries import WRITER,readbin
from test_output_protection import unchanged
REPO=Path(__file__).resolve().parents[1]
ALPHAS=np.array([0,.001,.004,.37,.8,1.])
DATA='input double(id x y)\n1 1 1\n2 2 1\n3 4 1\nend\n'

def dual(x,y,k,rts,orientation):
 n,m=x.shape;ss=y.shape[1];U=np.c_[-x,y];bounds=[(0,None)]*(m+ss)
 c=np.r_[np.zeros(m),-y[k]] if orientation=='input' else np.r_[x[k],np.zeros(ss)]
 E=np.r_[x[k],np.zeros(ss)][None,:] if orientation=='input' else np.r_[np.zeros(m),y[k]][None,:]
 if rts!='crs':
  U=np.c_[U,-np.ones(n)];E=np.c_[E,0];c=np.r_[c,1]
  bounds.append((None,None) if rts=='vrs' else (0,None) if rts=='drs' else (None,0))
 f=linprog(c,A_ub=U,b_ub=np.zeros(n),A_eq=E,b_eq=[1],bounds=bounds,method='highs-ds',options=dict(presolve=False,primal_feasibility_tolerance=1e-10,dual_feasibility_tolerance=1e-10,time_limit=2.))
 if not f.success:
  # Use a bounded-time independent alternate algorithm on the same dual.
  # Neither presolve nor relaxed constraints is used.
  f=linprog(c,A_ub=U,b_ub=np.zeros(n),A_eq=E,b_eq=[1],bounds=bounds,method='highs-ipm',options=dict(presolve=False,primal_feasibility_tolerance=1e-10,dual_feasibility_tolerance=1e-10,time_limit=2.))
 assert f.success,f.message
 assert max(U@f.x,default=0)<1e-7 and abs((E@f.x)[0]-1)<1e-7
 return -f.fun if orientation=='input' else 1/f.fun

def fixture(m,s,variant):
 if variant=='leon':
  # Leon2003 as distributed in deaR; R comparison also checks these source values.
  d=np.array([[3,2,3,1],[4,.5,2.5,1],[4.5,1.5,6,1],[6.5,.5,4,1.25],[7,2,5,.5],[8,.5,3.5,.5],[10,1,6,.5],[6,.5,2,1.5]])
  c=d[:,[0,2]];spread=d[:,[1,3]]
  return np.c_[np.arange(1,9),c,c-spread,c+spread,np.ones(8)],['id','x1','y1','lo_x1','lo_y1','hi_x1','hi_y1','selected']
 rng=np.random.default_rng(6637+10*m+s);c=rng.uniform(.7,9,(5,m+s));c[4,:m]=c[0,:m]*2;c[4,m:]=c[0,m:]*.5
 lo=c*rng.uniform(.55,.9,c.shape);hi=c*rng.uniform(1.1,1.65,c.shape)
 if variant=='crisp':lo=c.copy();hi=c.copy()
 if variant=='onesided':lo=c.copy()
 if variant=='inputs_only':lo[:,m:]=c[:,m:];hi[:,m:]=c[:,m:]
 if variant=='outputs_only':lo[:,:m]=c[:,:m];hi[:,:m]=c[:,:m]
 if variant=='zero':
  if m>1:c[1,0]=lo[1,0]=hi[1,0]=0
  if s>1:c[2,m]=lo[2,m]=hi[2,m]=0
 if variant=='zerocol':
  if m>1:c[:,m-1]=lo[:,m-1]=hi[:,m-1]=0
  if s>1:c[:,-1]=lo[:,-1]=hi[:,-1]=0
 if variant in ['duplicate','near']:
  for v in [c,lo,hi]:v[1]=v[0]*(1+1e-8 if variant=='near' else 1)
 if variant=='units':
  for v in [c,lo,hi]:v*=10.**np.linspace(-12,12,m+s)
 a=np.c_[np.arange(1,6),c,lo,hi,np.ones(5)]
 if variant=='permuted':a=a[[4,2,0,3,1]]
 if variant=='sample':a[1,-1]=0;a=np.vstack([a,a[0]]);a[-1,0]=6;a[-1,-1]=0;a[-1,1:1+m]*=.001
 vs=[f'x{i+1}' for i in range(m)]+[f'y{i+1}' for i in range(s)]
 return a,['id']+vs+['lo_'+v for v in vs]+['hi_'+v for v in vs]+['selected']

def oracle(a,m,s,rts,orientation,alphas=ALPHAS):
 p=m+s;c=a[:,1:1+p].copy();lo=a[:,1+p:1+2*p].copy();hi=a[:,1+2*p:1+3*p].copy()
 scale=hi.max(0);scale[scale==0]=1;c/=scale;lo/=scale;hi/=scale
 ref=np.zeros((len(a),len(alphas),2))
 for j,alpha in enumerate(alphas):
  L=(1-alpha)*lo+alpha*c;H=(1-alpha)*hi+alpha*c
  for k in range(len(a)):
   for sense in [0,1]:
    x=(L if sense==0 else H)[:,:m].copy();y=(H if sense==0 else L)[:,m:].copy()
    x[k]=(H if sense==0 else L)[k,:m];y[k]=(L if sense==0 else H)[k,m:]
    ref[k,j,sense]=dual(x,y,k,rts,orientation)
 return ref

def verify(d,m,s,rts,orientation):
 a=readbin(d/'data.bin');ref=oracle(a,m,s,rts,orientation);n=len(a)
 if not (d/'pess.bin').exists():return dict(independent_LPs=int(ref.size),numerical_aborts=1,unavailable=ref.size,max_gap=0.,scored=0)
 P=readbin(d/'pess.bin');O=readbin(d/'opt.bin');S=readbin(d/'summary.bin');complete=readbin(d/'complete.bin').ravel();M=readbin(d/'meta.bin').ravel();R=np.stack([P,O],axis=2);good=np.isfinite(R);gap=np.max(abs(R[good]-ref[good]),initial=0);assert gap<2e-6,(gap,d.name)
 co=good.all((1,2));assert np.array_equal(co,complete==1);assert M[0]==n and M[1]==co.sum() and M[2]==n-co.sum() and M[3]==(~good).sum()
 assert np.array_equal(S[:,0],ALPHAS)
 for j in range(len(ALPHAS)):
  expected=[sum(np.isfinite(P[:,j])),sum(np.isfinite(O[:,j])),sum(np.isfinite(P[:,j])&np.isfinite(O[:,j])),sum(co)];assert np.array_equal(S[j,1:5],expected)
  if co.any():assert np.allclose(S[j,5:],[np.mean(P[co,j]),np.mean(O[co,j]),min(P[co,j]),max(P[co,j]),min(O[co,j]),max(O[co,j])],rtol=0,atol=1e-14)
  else:assert np.isnan(S[j,5:]).all()
 if co.any():
  assert np.max(P[co]-O[co])<=1e-8;assert np.max(-np.diff(P[co],axis=1))<=1e-8;assert np.max(np.diff(O[co],axis=1))<=1e-8
 assert np.array_equal(P[:,-1],O[:,-1],equal_nan=True)
 return dict(independent_LPs=int(ref.size),numerical_aborts=0,unavailable=int((~good).sum()),max_gap=float(gap),scored=int(good.sum()))

R_SCRIPT=r'''
Sys.setenv(RGL_USE_NULL='TRUE');options(rgl.useNULL=TRUE)
suppressPackageStartupMessages(library(deaR));a<-commandArgs(TRUE)
d<-read.csv(a[1]);m<-as.integer(a[2]);s<-as.integer(a[3]);d<-d[d$selected==1,];p<-m+s
cent<-as.matrix(d[,2:(p+1)]);low<-as.matrix(d[,(p+2):(2*p+1)]);high<-as.matrix(d[,(2*p+2):(3*p+1)])
if(grepl('11_leon',a[1])) {data('Leon2003');stopifnot(max(abs(cent-cbind(Leon2003$x,Leon2003$y)))==0,max(abs(low-cbind(Leon2003$x-Leon2003$alpha,Leon2003$y-Leon2003$beta)))==0)}
scale<-apply(high,2,max);scale[scale==0]<-1;cent<-sweep(cent,2,scale,'/');low<-sweep(low,2,scale,'/');high<-sweep(high,2,scale,'/')
z<-data.frame(id=d$id,cent,cent-low,high-cent)
f<-make_deadata_fuzzy(z,inputs.mL=2:(m+1),outputs.mL=(m+2):(p+1),inputs.dL=(p+2):(p+m+1),outputs.dL=(p+m+2):(2*p+1),inputs.dR=(2*p+2):(2*p+m+1),outputs.dR=(2*p+m+2):(3*p+1))
al<-c(0,.001,.004,.37,.8,1);ans<-list()
for(rt in c('crs','vrs','nirs','ndrs')) for(ori in c('io','oo')) {
 res<-modelfuzzy_kaoliu(f,kaoliu_modelname='basic',alpha=al,orientation=ori,rts=rt,maxslack=FALSE)
 # efficiencies.dea_fuzzy rounds to five decimals; compare raw model values.
 e<-list(Worst=sapply(res$alphacut,function(cut) vapply(cut$DMU$Worst,function(d) d$efficiency,numeric(1))),Best=sapply(res$alphacut,function(cut) vapply(cut$DMU$Best,function(d) d$efficiency,numeric(1))))
 printed<-efficiencies(res);stopifnot(all(printed$Worst==round(e$Worst,5)),all(printed$Best==round(e$Best,5)))
 for(j in seq_along(al)) ans[[length(ans)+1]]<-data.frame(rts=rt,orientation=ori,alpha=al[j],id=d$id,pess=if(ori=='io') e$Worst[,j] else 1/e$Worst[,j],opt=if(ori=='io') e$Best[,j] else 1/e$Best[,j])
}
write.csv(do.call(rbind,ans),a[4],row.names=FALSE);cat('deaR',as.character(packageVersion('deaR')),'\n')
'''

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--only');ap.add_argument('--rscript');ap.add_argument('--baseline-ado',type=Path);a=ap.parse_args()
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir();h={}
 for n in package_files(REPO):
  raw=(REPO/n).read_bytes();(pkg/Path(n).name).write_bytes(raw);h[n]=hashlib.sha256(raw).hexdigest()
 (out/'manifest.json').write_text(json.dumps(dict(source_sha256=h,stata=a.stata),indent=2)+'\n');rows=[]
 def run(name,body,source=pkg,check=None):
  if a.only and not re.search(a.only,name):return
  d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset type double\nset processors 1\nset linesize 180\n'
  for f in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/f).mkdir();head+=f'sysdir set {f} "{d/f}"\n'
  head+=f'adopath ++ "{source}"\n'+WRITER;(d/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n');r=dict(test=name,passed=False);start=time.monotonic()
  try:
   p=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=180);log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else ''
   if p.returncode or not log:(d/'process.json').write_text(json.dumps(dict(returncode=p.returncode,stdout=p.stdout.decode(errors='replace'),stderr=p.stderr.decode(errors='replace'))))
   assert p.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2400:]
   if check:r.update(check(d))
   r['passed']=True
  except Exception as e:r['error']=repr(e)
  r['seconds']=round(time.monotonic()-start,3);rows.append(r);(out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(r),flush=True);return d
 specs=[]
 for m,s in [(1,1),(2,1),(1,2),(2,2)]:
  for var in ['base','units','crisp','onesided','inputs_only','outputs_only','zero','zerocol','duplicate','near','sample','permuted']+(['leon'] if (m,s)==(1,1) else []):
   z,cols=fixture(m,s,var);key=f'{m}{s}_{var}';f=out/('data_'+key+'.csv');np.savetxt(f,z,delimiter=',',header=','.join(cols),comments='',fmt='%.17g');specs.append((key,m,s,f))
   xs=' '.join(f'x{i+1}' for i in range(m));ys=' '.join(f'y{i+1}' for i in range(s));bs=' '.join(f'{side}{end}('+ ' '.join(f'{pre}_{v}' for v in (xs if side=='x' else ys).split())+')' for side,end,pre in [('x','lo','lo'),('x','hi','hi'),('y','lo','lo'),('y','hi','hi')])
   if var=='crisp':bs=''
   if var=='onesided':bs=re.sub(r'[xy]lo\([^)]*\)','',bs)
   if var=='inputs_only':bs=re.sub(r'y(lo|hi)\([^)]*\)','',bs)
   if var=='outputs_only':bs=re.sub(r'x(lo|hi)\([^)]*\)','',bs)
   for rts,ori in itertools.product(['crs','vrs','drs','irs'],['input','output']):
    body=f'import delimited using "{f}",clear asdouble\nquietly datasignature\nlocal before "`r(datasignature)\'"\n'
    body+=f'capture noisily dea_fuzzy if selected, inputs({xs}) outputs({ys}) {bs} dmu(id) rts({rts}) ort({ori}) alpha(0 .001 .004 .37 .8 1) gen(res) matrix detail\nlocal rc=_rc\n'
    body+='if `rc\'==0 {\nmatrix meta=(r(N),r(N_complete),r(N_partial),r(lp_missing),r(lp_count),r(lp_retries),r(grid_retries))\nmata: savebin("meta.bin",st_matrix("meta"))\n'
    for k,rn in [('pess','eff_pess'),('opt','eff_opt'),('complete','complete'),('summary','summary')]:body+=f'mata: savebin("{k}.bin",st_matrix("r({rn})"))\n'
    body+='foreach v of varlist res_* {\nassert missing(`v\') if !selected\n}\n}\nelse {\nassert `rc\'==498\nquietly datasignature\nassert "`before\'"=="`r(datasignature)\'"\n}\nkeep if selected\n'
    body+='mata: savebin("data.bin",st_data(.,tokens("'+' '.join(cols[:-1])+'")))\n'
    run(key+'_'+rts+'_'+ori,body,check=lambda p,m=m,s=s,rts=rts,ori=ori:verify(p,m,s,rts,ori))
 cmd='dea_fuzzy, inputs(x) outputs(y) dmu(id) alpha(0 .5 1) gen(res) matrix detail'
 for opt in ['ort(invalid)','ort(orange)','ort(outward)','rts(bad)','tol(0)','tol(.)']:
  run('invalid_'+re.sub(r'\W','_',opt),DATA+unchanged(cmd+' '+opt,198))
 for side,alpha in itertools.product(['x','y'],['0','0 .5 1','1']):
  run('zero_row_'+side+'_'+alpha.replace(' ','_'),DATA+f'replace {side}=0 in 2\n'+unchanged(cmd.replace('alpha(0 .5 1)',f'alpha({alpha})'),459))
 for side in ['x','y']:
  b=DATA+f'gen lo={side}\nreplace lo=0 in 2\n'
  run('zero_support_'+side,b+unchanged(cmd+f' {side}lo(lo)',459))
  run('positive_selected_cut_'+side,b+cmd.replace('alpha(0 .5 1)','alpha(.1 1)')+f' {side}lo(lo)\nassert r(N_complete)==3\n')
 run('workspace_reload',DATA+'matrix user=(9,8)\nlocal before : all matrices\n'+cmd+'\nlocal after : all matrices\nassert "`before\'"=="`after\'"\ndiscard\n'+cmd+' replacevars\nlocal after : all matrices\nassert "`before\'"=="`after\'"\n')
 labels=DATA+'gen str80 sid="Unit "+string(id)+char(34)+" quote "+char(96)+"probe"+char(39)\nlocal probe "EXPANDED_BAD"\n'+cmd.replace('dmu(id)','dmu(sid)')+'\n'
 run('literal_string_ids',labels,check=lambda d:dict(literal_ids=3) if 'quote `probe\'' in (d/'test.log').read_text() else (_ for _ in ()).throw(AssertionError('literal ID altered')))
 run('large_numeric_ids',DATA+'replace id=1234567890123+_n\n'+cmd+'\n',check=lambda d:dict(numeric_ids=3) if all(str(1234567890123+i) in (d/'test.log').read_text() for i in range(1,4)) else (_ for _ in ()).throw(AssertionError('numeric ID precision')))
 run('label_missing_kept',DATA+'replace id=. in 2\n'+cmd+'\nassert r(N)==3 & r(N_complete)==3\n')
 # Force failure after a temporary matrix exists: data and user matrices survive.
 target=out/'fault_summary_error';shutil.copytree(pkg,target);p=target/'dea_fuzzy.ado'
 p.write_text(p.read_text().replace('    st_matrix(destinations[1],T)','    st_matrix(destinations[1],T)\n    _error(3498)'))
 run('injected_summary_error',DATA+'matrix user=(9,8)\nlocal beforemat : all matrices\n'+unchanged(cmd,3498)+'local aftermat : all matrices\nassert "`beforemat\'"=="`aftermat\'"\n',source=target)
 for fault in ['partial','all','gap','crisp']:
  target=out/('fault_'+fault);shutil.copytree(pkg,target);p=target/'dea_fuzzy.ado';src=p.read_text()
  code={'partial':'eff_pess=(.9\\.\\.5); eff_opt=(.\\.1\\.5)','all':'eff_pess=J(n,1,.); eff_opt=J(n,1,.)','gap':'eff_pess[1]=(alpha==0 ? .9 : (alpha==.5 ? . : .1)); eff_opt[1]=(alpha==1 ? .1 : 1)','crisp':'eff_pess[1]=.1; eff_opt[1]=.9'}[fault]
  p.write_text(src.replace('    st_matrix(stats_name,lpstats)',code+'\n    st_matrix(stats_name,lpstats)'))
  # Corrupt the final temporary center outputs as well: the new exact center
  # route otherwise repairs the earlier solver-level injection. These cases
  # test summary/nesting checks, independently of recovery correctness.
  center=target/'_dea_fuzzy_exact.ado'
  if center.exists():
   cs=center.read_text();hook=' st_store(obs,st_varindex(pess),P);st_store(obs,st_varindex(opt),O)'
   assert hook in cs
   late=code.replace('eff_pess','P').replace('eff_opt','O').replace('alpha','1')
   center.write_text(cs.replace(hook,late+'\n'+hook))

  if fault in ['gap','crisp']:
   body=DATA+'gen str12 res_pess_a0="old value"\n'+unchanged(cmd+' replacevars',498)
  else:
   body=DATA+cmd+'\nmatrix S=r(summary)\n'+('assert r(N_complete)==1 & r(N_partial)==2\nmata: assert(all(st_matrix("S")[,2..3]:==2) & all(st_matrix("S")[,4..5]:==1) & all(st_matrix("S")[,6..11]:==.5))\n' if fault=='partial' else 'assert r(N_complete)==0 & r(N_partial)==3\nmata: assert(all(st_matrix("S")[,2..5]:==0) & all(st_matrix("S")[,6..11]:>=.))\n')
  run('injected_'+fault,body,source=target)
 if a.baseline_ado:
  baseline=out/'before_ado';shutil.copytree(a.baseline_ado,baseline)
  for fault in ['partial','gap','crisp']:
   target=out/('old_'+fault);shutil.copytree(baseline,target);p=target/'dea_fuzzy.ado';src=p.read_text()
   code={'partial':'eff_pess=(.9\\.\\.5); eff_opt=(.\\.1\\.5)','gap':'eff_pess[1]=(alpha==0 ? .9 : (alpha==.5 ? . : .1)); eff_opt[1]=(alpha==1 ? .1 : 1)','crisp':'eff_pess[1]=.1; eff_opt[1]=.9'}[fault];p.write_text(src.replace('    st_matrix(stats_name,lpstats)',code+'\n    st_matrix(stats_name,lpstats)'))
   body=DATA+cmd.replace(' matrix','')+'\n'
   if fault=='partial':body+='quietly summarize res_pess_a0\nscalar meanp=r(mean)\nquietly summarize res_opt_a0\nassert meanp>r(mean)\n'
   if fault=='gap':body+='assert res_pess_a0>res_pess_a100 in 1\n'
   if fault=='crisp':body+='assert res_opt_a100-res_pess_a100>.5 in 1\n'
   run('old_accepted_'+fault,body,source=target)
 # Tiny full-box enumeration is independent of the best/worst construction.
 def corners(d):
  c=np.array([[1.,1.],[1.,1.]]);lo=c*.7;hi=c*1.4;calls=0;strict=0
  for rts,ori,alpha in itertools.product(['crs','vrs','drs','irs'],['input','output'],[0,.4,1]):
   aa=np.c_[np.arange(1,3),c,lo,hi];ref=oracle(aa,1,1,rts,ori,np.array([alpha]));L=(1-alpha)*lo+alpha*c;H=(1-alpha)*hi+alpha*c;vals=[]
   for bits in itertools.product([0,1],repeat=4):
    z=np.where(np.array(bits).reshape(2,2),H,L);vals.append([dual(z[:,:1],z[:,1:],k,rts,ori) for k in range(2)]);calls+=2
   vals=np.array(vals);assert np.allclose(vals.min(0),ref[:,0,0],atol=1e-8,rtol=0) and np.allclose(vals.max(0),ref[:,0,1],atol=1e-8,rtol=0)
   # Marginal endpoint means need not be attainable together in one dataset.
   lower=ref[:,0,0].mean();upper=ref[:,0,1].mean();joint=vals.mean(1)
   assert joint.min()>=lower-1e-8 and joint.max()<=upper+1e-8
   strict+=int(joint.min()>lower+1e-6 or joint.max()<upper-1e-6)
  assert strict>0
  return dict(corner_LPs=calls,marginal_mean_not_joint_interval_cases=strict)
 run('independent_box_extrema','',check=corners)
 if a.rscript:
  (out/'compare.R').write_text(R_SCRIPT)
  def rcheck(d):
   count=0;missing=0;gap=0
   for key,m,s,f in specs:
    if not (key.endswith('_base') or key.endswith('_leon')):continue
    dest=d/(key+'.csv');p=subprocess.run([a.rscript,str(out/'compare.R'),str(f),str(m),str(s),str(dest)],capture_output=True,text=True,timeout=180);(d/(key+'.console')).write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr
    with dest.open() as fh:ref=list(csv.DictReader(fh))
    for rts,ori in itertools.product(['crs','vrs','drs','irs'],['input','output']):
     q=out/(key+'_'+rts+'_'+ori)
     if not q.exists():continue
     tag={'drs':'nirs','irs':'ndrs'}.get(rts,rts);rr=[row for row in ref if row['rts']==tag and row['orientation']==('io' if ori=='input' else 'oo')]
     for v in ['pess','opt']:
      got=readbin(q/(v+'.bin'));want=np.array([float(row[v]) for row in rr]).reshape(len(ALPHAS),-1).T;ok=np.isfinite(got);missing+=int((~ok).sum());g=np.max(abs(got[ok]-want[ok]),initial=0);assert g<2e-6,(key,rts,ori,v,g);gap=max(gap,g);count+=int(ok.sum())
   assert count>0
   return dict(R_scores=count,R_missing=missing,R_max_gap=float(gap))
  run('R_deaR','',check=rcheck)
 summary=dict(total=len(rows),passed=sum(r['passed'] for r in rows),failed=sum(not r['passed'] for r in rows));(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary);return bool(summary['failed'])
if __name__=='__main__':raise SystemExit(main())
