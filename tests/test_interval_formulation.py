#!/usr/bin/env python3
"""Consistent-realization interval DEA: independent crisp LPs and deaR.

Requires Stata 16+, NumPy/SciPy; --rscript requires installed deaR. Uses isolated
ado copies. The oracle constructs complete realized datasets and solves ordinary
radial DEA; it never reads production LP coefficients. Optional baseline is the
pre-R09 (2.0.55) package. R output-oriented phi is inverted before comparison.
"""
from _support import package_files
from pathlib import Path
import argparse, csv, hashlib, itertools, json, os, re, shutil, subprocess, time
import numpy as np
import scipy
from verify_bootstrap_lp import solve
REPO=Path(__file__).resolve().parents[1]
SPECS={'crs':'crs','ccr':'crs','vrs':'vrs','bcc':'vrs','drs':'drs','nirs':'drs','irs':'irs','ndrs':'irs'}
BASE=np.array([[1,1,1,1],[2,3,1,2],[4,2,3,1],[3,4,2,4],[5,5,2,2],[1,2,2,3]],float)

def data(dim,variant):
 a=BASE.copy()
 if variant=='singleton':a=np.ones((1,4))
 if variant=='two':a=np.array([[1,1,1,1],[2,3,1,2]],float)
 if variant.startswith('random'):a=np.random.default_rng(7300+int(variant[6:])).uniform(.2,8,(7,4))
 if variant=='zero':a[0,0]=0;a[1,3]=0
 if variant=='units':a*=np.array([.001,1e3,100,.01])
 if variant=='tiny':a*=1e-9
 if variant=='large':a*=1e9
 if variant=='duplicate':a=np.vstack([a,a[1]])
 if variant=='permuted':a=a[[5,2,0,4,1,3]]
 cols=['x'+str(i+1) for i in range(dim)]+['y'+str(i+1) for i in range(dim)]
 a=a[:,list(range(dim))+list(range(2,2+dim))]
 script='input double(id '+' '.join(cols)+')\n'+'\n'.join(str(j+1)+' '+' '.join(format(v,'.17g') for v in row) for j,row in enumerate(a))+'\nend\n'
 q='if id!=2 in 1/5' if variant=='sample' else ''
 script+='gen byte selected=0\nreplace selected=1 '+q+'\n'
 # Unequal widths expose accidental swaps, and a visibly dominated DMU is included.
 for v in cols:
  script+=f'gen double l_{v}=.8*{v}\ngen double h_{v}=1.2*{v}\n'
 return script,cols,q

def options(cols,mode):
 opts=[]
 for side in ['x','y']:
  vs=[v for v in cols if v.startswith(side)]
  if mode=='point':opts += [side+'lo('+' '.join(vs)+')',side+'hi('+' '.join(vs)+')'];continue
  if mode in ['full','lower',side]:opts.append(side+'lo('+' '.join('l_'+v for v in vs)+')')
  if mode in ['full','upper',side]:opts.append(side+'hi('+' '.join('h_'+v for v in vs)+')')
 return ' '.join(opts)

def read(p):
 with p.open() as f:rr=list(csv.DictReader(f))
 return {k:np.array([float(r[k]) if r[k] not in ['', '.'] and not re.fullmatch(r'\.[a-z]',r[k]) else np.nan for r in rr]) for k in rr[0] if k.lower()!='rts'}

def corners(xl,xu,yl,yu,ori,rts):
 # Column normalization is an exact unit change in the independent ordinary LP.
 sx=np.maximum(xu.max(axis=0),1e-300);sy=np.maximum(yu.max(axis=0),1e-300)
 xl,xu,yl,yu=xl/sx,xu/sx,yl/sy,yu/sy
 ans=[]
 for worst in [True,False]:
  vals=[]
  for k in range(len(xl)):
   x,y=(xl.copy(),yu.copy()) if worst else (xu.copy(),yl.copy())
   x[k],y[k]=(xu[k],yl[k]) if worst else (xl[k],yu[k])
   vals.append(solve(x[k],y[k],x,y,ori,rts))
  ans.append(vals)
 return np.array(ans).T

R_SCRIPT=r'''
suppressPackageStartupMessages(library(deaR))
args=commandArgs(TRUE); jobs=read.csv(args[1],stringsAsFactors=FALSE)
ans=list(); i=0
for (j in seq_len(nrow(jobs))) {
 job=jobs[j,]; d=read.csv(job$file); p=job$dim
 f=make_deadata_fuzzy(d,dmus=1,inputs.mL=2:(1+p),inputs.mR=(2+p):(1+2*p),
     outputs.mL=(2+2*p):(1+3*p),outputs.mR=(2+3*p):(1+4*p))
 for (rts in c('crs','vrs','nirs','ndrs')) {
  z=modelfuzzy_kaoliu(f,alpha=0,orientation=if(job$ori=='input') 'io' else 'oo',rts=rts,maxslack=FALSE)
  # efficiencies() rounds to five digits by default; compare raw solver results.
  e=z$alphacut[[1]]$DMU
  lo=vapply(e$Worst,function(v) v$efficiency,numeric(1))
  hi=vapply(e$Best,function(v) v$efficiency,numeric(1))
  if(job$ori=='output') {lo=1/lo;hi=1/hi}
  i=i+1;ans[[i]]=data.frame(test=job$test,rts=rts,id=d$id,pess=lo,opt=hi)
 }
}
write.csv(do.call(rbind,ans),args[2],row.names=FALSE)
writeLines(c(paste('deaR',packageVersion('deaR')),capture.output(sessionInfo())),args[3])
'''

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True)
 ap.add_argument('--rscript');ap.add_argument('--baseline-ado',type=Path);ap.add_argument('--only');args=ap.parse_args()
 out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir();hashes={}
 for name in package_files(REPO):
  raw=(REPO/name).read_bytes();(pkg/Path(name).name).write_bytes(raw);hashes[name]=hashlib.sha256(raw).hexdigest()
 manifest=dict(source_sha256=hashes,numpy=np.__version__,scipy=scipy.__version__,stata=args.stata,baseline=str(args.baseline_ado))
 if args.baseline_ado:manifest['baseline_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in args.baseline_ado.iterdir() if p.is_file()}
 (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');rows=[];jobs=[];expected={}
 def chosen(n):return not args.only or re.search(args.only,n)
 def run(name,body,check=None,source=pkg):
  if not chosen(name):return
  d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 120\n'
  for k in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/k).mkdir();head+=f'sysdir set {k} "{d/k}"\n'
  head+=f'adopath ++ "{source}"\n';(d/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n')
  r=dict(test=name,passed=False);start=time.monotonic()
  try:
   p=subprocess.run([args.stata,'-b','do','test.do'],cwd=d,timeout=240,capture_output=True)
   if p.returncode or not (d/'test.log').exists():
    (d/'process.json').write_text(json.dumps(dict(returncode=p.returncode,stdout=p.stdout.decode(errors='replace'),stderr=p.stderr.decode(errors='replace')),indent=2))
   log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else p.stderr.decode(errors='replace')
   assert p.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2200:]
   if check:r.update(check(d))
   r['passed']=True
  except Exception as e:r['error']=repr(e)
  r['seconds']=round(time.monotonic()-start,3);rows.append(r);(out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(r),flush=True)

 variants=['regular','sample','permuted','duplicate']
 grid=list(itertools.product([1,2],['input','output'],['point','full','x','y','lower','upper'],variants))
 grid+=list(itertools.product([2],['input','output'],['full'],['zero','units','tiny','large','random1','random2','random3','unused_missing']))
 grid+=list(itertools.product([1],['input','output'],['full'],['singleton','two']))
 for dim,ori,mode,variant in grid:
  name=f'{dim}d_{ori}_{mode}_{variant}'
  if not chosen(name):continue
  body,cols,q=data(dim,variant);xs=' '.join(v for v in cols if v[0]=='x');ys=' '.join(v for v in cols if v[0]=='y')
  opts=options(cols,mode)
  if variant=='unused_missing':body+='foreach v of varlist '+ ' '.join(cols) +" {\n replace `v'=.z\n}\n"
  body+='count if selected\nlocal n=r(N)\n'
  for i,(spec,norm) in enumerate(SPECS.items()):
   engine=['simplex','ipm','auto'][i%3];tol=['1e-7','1e-9'][i%2]
   head=f'dea {xs}={ys} {q},' if i%2==0 else f'dea {ys} {q}, inputs({xs})'
   body+=f'{head} {opts} rts({spec.upper() if i%2 else spec}) ort({ori}) engine({engine}) stage({1+i%2}) tol({tol}) dmu(id) gen(e_{spec})\n'
   body+=f'''assert "`r(idea_method)'"=="interval" & "`r(rts)'"=="{norm}"
assert r(N)==`n' & r(n_failed)==0 & r(n_failed_primary)==0 & r(n_failed_opt)==0
assert !missing(e_{spec},e_{spec}_opt) if selected
assert e_{spec}>=-1e-8 & e_{spec}<=e_{spec}_opt+1e-8 & e_{spec}_opt<=1+1e-8 if selected
assert missing(e_{spec},e_{spec}_opt) if !selected
'''
  if variant=='regular' and mode in ['point','full']:
   for spec in ['crs','vrs','drs','irs']:
    body+=f'dea_fuzzy {q}, inputs({xs}) outputs({ys}) {opts} rts({spec}) ort({ori}) alpha(0) dmu(id) gen(f_{spec})\n'
    body+=f'assert abs(f_{spec}_pess_a0-e_{spec})<1e-7 & abs(f_{spec}_opt_a0-e_{spec}_opt)<1e-7 if selected\n'
    if mode=='point':
     for engine in ['simplex','ipm','auto']:
      body+=f'dea {xs}={ys}, rts({spec}) ort({ori}) engine({engine}) tolerance(1e-9) stage(1) dmu(id) gen(c_{spec}_{engine}) replacevars\n'
      body+=f'assert abs(c_{spec}_{engine}-e_{spec})<1e-7 & abs(c_{spec}_{engine}-e_{spec}_opt)<1e-7\n'
   if mode=='full':
    body+='preserve\n'
    for v in cols:body+=f'replace l_{v}=(l_{v}+{v})/2\nreplace h_{v}=(h_{v}+{v})/2\n'
    for spec in ['crs','vrs','drs','irs']:
     body+=f'dea {xs}={ys}, {opts} rts({spec}) ort({ori}) dmu(id) gen(n_{spec})\n'
     body+=f'assert !missing(n_{spec},n_{spec}_opt) & n_{spec}>=e_{spec}-1e-7 & n_{spec}_opt<=e_{spec}_opt+1e-7\n'
    body+='restore\n'
    # Save the exact numerical intervals, then compare the legacy optimistic path.
    body+='preserve\n'
    for v in cols:body+=f'gen strL b_{v}="["+strtrim(string(l_{v},"%24.17g"))+","+strtrim(string(h_{v},"%24.17g"))+"]"\ndrop {v}\nrename b_{v} {v}\n'
    for spec,engine in itertools.product(['crs','vrs','drs','irs'],['simplex','ipm','auto']):
     body+=f'dea {xs}={ys}, bounded({xs} {ys}) rts({spec}) ort({ori}) engine({engine}) tolerance(1e-9) dmu(id) gen(b_{spec}_{engine}) replacevars\n'
     body+=f'assert abs(b_{spec}_{engine}-e_{spec}_opt)<1e-7\n'
    body+='restore\n'
  body+="ds, has(type numeric)\nformat `r(varlist)' %24.17g\nexport delimited using \"scores.csv\",replace datafmt\n"
  def check(d,dim=dim,ori=ori,mode=mode,variant=variant,cols=cols,name=name):
   got=read(d/'scores.csv');sel=got['selected']==1
   def bound(side,end):
    use=mode in ['full',side, 'lower' if end=='l' else 'upper']
    return np.column_stack([got[(end+'_' if use else '')+v][sel] for v in cols if v[0]==side])
   xl,xu,yl,yu=bound('x','l'),bound('x','h'),bound('y','l'),bound('y','h')
   gap=0.;count=0
   for rts in ['crs','vrs','drs','irs']:
    ref=corners(xl,xu,yl,yu,ori,rts);count+=len(xl)*2
    for spec,norm in SPECS.items():
     if norm!=rts:continue
     vals=np.column_stack([got['e_'+spec][sel],got['e_'+spec+'_opt'][sel]])
     err=float(abs(vals-ref).max());assert err<1e-7,(spec,err,vals,ref);gap=max(gap,err)
    expected[name,rts]=np.column_stack([got['e_'+rts][sel],got['e_'+rts+'_opt'][sel]])
   if variant=='singleton':assert np.max(abs(expected[name,'crs']-1))<1e-8
   prop=0
   # For two DMUs, enumerate every corner and interior draws; this independently
   # establishes min/max over all corner realizations without the endpoint recipe.
   if variant=='two':
    lower=np.c_[xl,yl];upper=np.c_[xu,yu];random=np.random.default_rng(42)
    realizations=[lower+(upper-lower)*np.array(z).reshape(lower.shape) for z in itertools.product([0,1],repeat=lower.size)]
    realizations += [lower+(upper-lower)*random.uniform(size=lower.shape) for j in range(25)]
    for rts in ['crs','vrs','drs','irs']:
     vals=[]
     for z in realizations:
      vals.append([solve(z[k,:dim],z[k,dim:],z[:,:dim],z[:,dim:],ori,rts) for k in range(len(z))]);prop+=len(z)
     vals=np.array(vals);ref=expected[name,rts]
     assert np.max(abs(vals.min(axis=0)-ref[:,0]))<1e-7 and np.max(abs(vals.max(axis=0)-ref[:,1]))<1e-7
   # Nested half-width bounds, compared against the observed full-width scores.
   if mode=='full':
    cx=(xl+xu)/2;cy=(yl+yu)/2
    for rts in ['crs','vrs','drs','irs']:
     narrow=corners((xl+cx)/2,(xu+cx)/2,(yl+cy)/2,(yu+cy)/2,ori,rts);prop+=2*len(xl)
     ref=expected[name,rts];assert np.all(ref[:,0]<=narrow[:,0]+1e-7) and np.all(ref[:,1]>=narrow[:,1]-1e-7)
   # R uses a common unit change to avoid comparing unrelated solver scaling limits.
   sx=np.maximum(xu.max(axis=0),1e-300);sy=np.maximum(yu.max(axis=0),1e-300)
   mat=np.c_[got['id'][sel],xl/sx,xu/sx,yl/sy,yu/sy];header=['id']+[f'{v}{i}' for v in ['xl','xu','yl','yu'] for i in range(dim)]
   with (d/'r_input.csv').open('w') as f:w=csv.writer(f);w.writerow(header);w.writerows(mat)
   jobs.append(dict(test=name,file=str(d/'r_input.csv'),dim=dim,ori=ori))
   return dict(lp_checks=count,property_lp_checks=prop,max_gap=gap,command_grid=8)
  run(name,body,check)

 # Exhaustive public syntax combinations on a shared analytical LP oracle.
 for ori in ['input','output']:
  body,cols,q=data(2,'regular');opts=options(cols,'full');idx=0
  ref={rts:corners(.8*BASE[:,:2],1.2*BASE[:,:2],.8*BASE[:,2:],1.2*BASE[:,2:],ori,rts) for rts in ['crs','vrs','drs','irs']}
  for spec,engine,form,stage in itertools.product(SPECS,['simplex','ipm','auto'],['A','B'],[1,2]):
   head='dea x1 x2=y1 y2,' if form=='A' else 'dea y1 y2, inputs(x1 x2)'
   model=['','model(te)','model(imprecise)'][idx%3]
   body+=f'{head} {opts} {model} rts({spec}) ort({ori}) engine({engine}) stage({stage}) dmu(id) gen(g{idx})\n'
   body+='assert r(n_failed)==0\n'
   for k,(lo,hi) in enumerate(ref[SPECS[spec]],1):body+=f'assert abs(g{idx}-{lo:.17g})<1e-7 & abs(g{idx}_opt-{hi:.17g})<1e-7 if id=={k}\n'
   idx+=1
  run('syntax_grid_'+ori,body,check=lambda d:dict(command_grid=96,score_checks=1152))

 # Observe tolerance at the actual nested radial engine, for bounded-only,
 # ordinal-only and mixed numeric/ordinal transformations. Test-only copies.
 trace=out/'trace_ado';shutil.copytree(pkg,trace)
 for engine in ['ipm','simplex']:
  p=trace/f'_dea_{engine}.ado';s=p.read_text();start=s.index(f'void _dea_{engine}_envelopment(');at=s.index('    X = st_data',start)
  s=s[:at]+'    assert(tol==st_numscalar("expected_tol"))\n'+s[at:];p.write_text(s)
 for kind,engine,tol in itertools.product(['bounded','ordinal','mixed'],['simplex','ipm','auto'],['1e-7','1e-9']):
  name=f'legacy_tolerance_{kind}_{engine}_{tol}'
  body,cols,q=data(2,'regular');body+=f'scalar expected_tol={tol}\n'
  opts={'bounded':'bounded(x1 y1)','ordinal':'ordinal(x2 y2)','mixed':'xlo(l_x1 x2) xhi(h_x1 x2) ordinal(x2 y2)'}[kind]
  body+=f'dea x1 x2=y1 y2, {opts} rts(vrs) ort(output) engine({engine}) tolerance({tol}) dmu(id) gen(e)\nassert !missing(e)\n'
  run(name,body,source=trace)

 if args.baseline_ado:
  # Freeze all eight old singleton cases; this is an error reproduction, not an oracle.
  for ori in ['input','output']:
   body,cols,q=data(1,'singleton');opts=options(cols,'full')
   for rts in ['crs','vrs','drs','irs']:
    body+=f'dea x1=y1, {opts} rts({rts}) ort({ori}) dmu(id) gen(e_{rts})\n'
   body+="ds, has(type numeric)\nformat `r(varlist)' %24.17g\nexport delimited using \"scores.csv\",replace datafmt\n"
   run('old_singleton_'+ori,body,source=args.baseline_ado.resolve())
  body,cols,q=data(1,'singleton')
  body=f'quietly run "{args.baseline_ado.resolve()}/_dea_ipm.ado"\nmata: assert(_dea_ipm_audit_ready()==2)\n'+body
  body+='dea x1=y1, '+options(cols,'full')+' rts(vrs) ort(output) gen(e)\nmata: assert(_dea_ipm_audit_ready()==4)\nassert abs(e-1)<1e-8 & abs(e_opt-1)<1e-8\n'
  run('reload_from_r08',body)
 if args.rscript and jobs:
  with (out/'r_jobs.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=['test','file','dim','ori']);w.writeheader();w.writerows(jobs)
  (out/'compare.R').write_text(R_SCRIPT);start=time.monotonic();row=dict(test='external_deaR',passed=False)
  try:
   p=subprocess.run([args.rscript,str(out/'compare.R'),str(out/'r_jobs.csv'),str(out/'r_scores.csv'),str(out/'R_session.txt')],env=dict(os.environ,RGL_USE_NULL='TRUE'),capture_output=True,text=True,timeout=600)
   (out/'R.log').write_text(p.stdout+'\n'+p.stderr);assert p.returncode==0,p.stderr[-1800:]
   with (out/'r_scores.csv').open() as f:rr=list(csv.DictReader(f))
   offsets={};gap=0
   for r in rr:
    key=r['test'],{'nirs':'drs','ndrs':'irs'}.get(r['rts'],r['rts']);i=offsets.get(key,0);offsets[key]=i+1
    err=float(np.max(abs(expected[key][i]-np.array([float(r['pess']),float(r['opt'])]))));assert err<1e-7,(key,i,err);gap=max(gap,err)
   assert sum(len(z) for z in expected.values())==len(rr)
   row.update(passed=True,score_checks=2*len(rr),max_gap=gap)
  except Exception as e:row['error']=repr(e)
  row['seconds']=round(time.monotonic()-start,3);rows.append(row);(out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(row),flush=True)
 summary=dict(total=len(rows),passed=sum(r['passed'] for r in rows),failed=sum(not r['passed'] for r in rows),lp_checks=sum(r.get('lp_checks',0) for r in rows),property_lp_checks=sum(r.get('property_lp_checks',0) for r in rows))
 (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
 raise SystemExit(summary['failed']>0)
if __name__=='__main__':main()
