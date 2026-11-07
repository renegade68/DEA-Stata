#!/usr/bin/env python3
"""Input VRS Lim-Zhu peer formula: independent LPs, paper weights, R and compatibility.

The correction changes peer appraisal, not primary or secondary LP objectives.
R comparisons separate solver-selected weights from same-weight formula checks.
"""
from _support import package_files
from pathlib import Path
import argparse,csv,hashlib,itertools,json,os,re,shutil,subprocess,time
import numpy as np
from test_core_retries import WRITER,readbin
from test_crosseff_models import fixture,verify,problem,fit
from test_output_protection import unchanged
from test_nonte_tolerance import instrument,CAPTURE,verify as verify_lp
REPO=Path(__file__).resolve().parents[1]
PAPER=np.c_[np.arange(1,8),[5/3,2,4,3,6,5,7],np.arange(1,8),np.ones(7)]
DATA='input double(id x y)\n1 1 1\n2 2 1\n3 3 2\nend\n'
R_SCRIPT=r'''
suppressPackageStartupMessages(library(deaR))
a=commandArgs(TRUE);jobs=read.csv(a[1]);ans=list();ii=0
for(j in seq_len(nrow(jobs))) {
 q=jobs[j,];z=read.csv(q$file);z=z[z$selected==1,];xs=paste0('x',seq_len(q$m));ys=paste0('y',seq_len(q$s))
 dat=make_deadata(z,dmus=1,inputs=match(xs,names(z)),outputs=match(ys,names(z)))
 sol=model_multiplier(dat,orientation='io',rts='vrs',epsilon=0,compute_lambda=FALSE)
 independent=cross_efficiency(sol,rts='vrs',correction=TRUE,selfapp=FALSE,M2=FALSE,M3=FALSE)$Arbitrary
 # Public dea input object with Stata-selected weights; original R function unchanged.
 supplied=read.csv(q$weights)
 for(k in seq_len(nrow(z))) {
  sol$DMU[[k]]$multiplier_input=as.numeric(supplied[k,paste0('v',seq_len(q$m)),drop=TRUE])
  sol$DMU[[k]]$multiplier_output=as.numeric(supplied[k,paste0('u',seq_len(q$s)),drop=TRUE])
  sol$DMU[[k]]$multiplier_rts=supplied$w[k]
  sol$DMU[[k]]$efficiency=supplied$self[k]
 }
 fixed=cross_efficiency(sol,rts='vrs',correction=TRUE,selfapp=FALSE,M2=FALSE,M3=FALSE)$Arbitrary
 for(k in seq_len(nrow(z)))for(l in seq_len(nrow(z))) {
  ii=ii+1;ans[[ii]]=data.frame(test=q$test,row=k,col=l,self=independent$efficiency[k],fixed=fixed$cross_eff[k,l],independent=independent$cross_eff[k,l],mean=fixed$e[l])
 }
}
write.csv(do.call(rbind,ans),a[2],row.names=FALSE,na='NA')
writeLines(c(paste('deaR',packageVersion('deaR')),capture.output(sessionInfo())),a[3])
'''

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--source-ado',type=Path);ap.add_argument('--baseline-ado',type=Path);ap.add_argument('--rscript');ap.add_argument('--only');a=ap.parse_args()
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir();h={}
 paths=list(a.source_ado.iterdir()) if a.source_ado else [REPO/n for n in package_files(REPO)]
 for p in paths:
  if p.suffix in ['.ado','.do','.sthlp']:raw=p.read_bytes();(pkg/p.name).write_bytes(raw);h[p.name]=hashlib.sha256(raw).hexdigest()
 (out/'manifest.json').write_text(json.dumps(dict(source_sha256=h,stata=a.stata),indent=2)+'\n');rows=[];jobs=[]
 audit=out/'audit_ado';instrument(pkg,audit)
 def run(name,body,check=None,source=pkg):
  if a.only and not re.search(a.only,name):return
  d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset processors 1\nset type double\nset linesize 160\n'
  for f in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/f).mkdir();head+=f'sysdir set {f} "{d/f}"\n'
  head+=f'adopath ++ "{source}"\n'+WRITER
  (d/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n');r=dict(test=name,passed=False);t=time.monotonic()
  try:
   proc=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=180);log=(d/'test.log').read_text(errors='replace');assert proc.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2500:]
   if check:r.update(check(d))
   r['passed']=True
  except Exception as e:r['error']=repr(e)
  r['seconds']=round(time.monotonic()-t,3);rows.append(r);(out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(r),flush=True)
 def save_fixture(name,z,cols):
  file=out/(name+'.csv')
  with file.open('w') as f:w=csv.writer(f);w.writerow(cols);w.writerows(['' if np.isnan(v) else format(v,'.17g') for v in row] for row in z)
  return file
 specs=[]
 for m,s,var in itertools.product([1,2],[1,2],['base','units','tiny','zero','zerocol','duplicate','sample','permuted','near','singleton']):
  z,cols=fixture(m,s,var);specs.append((f'{m}x{s}_{var}',m,s,save_fixture(f'data_{m}x{s}_{var}',z,cols)))
 specs.append(('paper7',1,1,save_fixture('data_paper7',PAPER,['id','x1','y1','selected'])))
 if a.rscript:
  script=out/'export_data.R';script.write_text("suppressPackageStartupMessages(library(deaR))\ndata('Lim_Zhu_2015')\nz=as.data.frame(Lim_Zhu_2015)\nnames(z)=c('id','x1',paste0('y',1:5))\nz$id=seq_len(nrow(z))\nz$selected=1\nwrite.csv(z,commandArgs(TRUE)[1],row.names=FALSE)\n")
  file=out/'data_paper37.csv';p=subprocess.run([a.rscript,str(script),str(file)],capture_output=True,env=dict(os.environ,RGL_USE_NULL='TRUE'),timeout=60);assert p.returncode==0,p.stderr.decode();specs.append(('paper37',1,5,file))
 for key,m,s,file in specs:
  xs=' '.join(f'x{i+1}' for i in range(m));ys=' '.join(f'y{i+1}' for i in range(s))
  for ci,(method,tol) in enumerate(itertools.product(['arbitrary','benevolent','aggressive'],['1e-4','1e-7','1e-10'])):
   name=f'{key}_{method}_{tol}'
   cmd=f'dea_crosseff {xs}={ys} if selected,' if ci%2 else f'dea_crosseff if selected, in({xs}) out({ys})'
   cmd+=f' dmu(id) r(v) ort(i) method({method}) tol({tol}) gen(res) vrsc(limzhu) matrix detail'
   body=f'import delimited using "{file}",clear varnames(1) asdouble\n'+cmd+'\n'
   body+='assert "`r(peer_rule)\'"=="limzhu_input"\nassert "`r(vrscross)\'"=="limzhu"\nmata: assert(st_numscalar("r(n_missing_diagonal)")==sum(diagonal(st_matrix("r(ce_matrix)")):>=.))\n'
   body+='matrix meta=(r(n),r(mean_ce),r(mean_self),r(n_failed_primary),r(n_failed_secondary),r(n_failed_ce),r(n_missing_peer_pairs),r(tolerance))\n'
   for f,mat in [('ce','ce_matrix'),('v','weights_input'),('u','weights_output'),('w','intercept'),('status','stage_status'),('peers','peer_count')]:body+=f'mata: savebin("{f}.bin",st_matrix("r({mat})"))\n'
   body+='mata: savebin("meta.bin",st_matrix("meta"))\nassert missing(res,res_self,res_maverick) if !selected\n'
   body+=f'keep if selected\nmata: savebin("data.bin",st_data(.,tokens("id {xs} {ys}")))\nmata: savebin("scores.bin",st_data(.,tokens("res res_self res_maverick")))\n'
   run(name,body,lambda d,m=m,s=s,method=method:verify(d,m,s,'vrs','input',method,'limzhu'))
   if (key.endswith('_base') or key.startswith('paper')) and method=='arbitrary' and tol=='1e-7' and (out/name/'scores.bin').exists():
    d=out/name;weights=d/'weights.csv';mat=np.c_[readbin(d/'v.bin'),readbin(d/'u.bin'),readbin(d/'w.bin'),readbin(d/'scores.bin')[:,1]]
    with weights.open('w') as f:w=csv.writer(f);w.writerow([f'v{i+1}' for i in range(m)]+[f'u{i+1}' for i in range(s)]+['w','self']);w.writerows(mat)
    jobs.append(dict(test=name,file=str(file),weights=str(weights),m=m,s=s))
 # Paper Table 1 supplies explicit optimal weights; do not assume a solver chooses them.
 body=f'run "{pkg}/dea_crosseff.ado"\nmata:\nx=(5/3\\2\\4\\3\\6\\5\\7); y=(1::7)\nv=(3/5\\1/2\\1/4\\1/3\\1/6\\1/5\\1/7); u=(0\\1/6\\1/8\\1/6\\1/6\\1/5\\2/7); w=(1\\2/3\\1/4\\1/3\\-1/6\\-1/5\\-1)\nce=J(7,7,.)\nfor (i=1;i<=7;i++) for(j=1;j<=7;j++) ce[i,j]=_dea_crosseff_lzratio(v[i]*x[j],u[i]*y[j],w[i],1e-10)\nsavebin("paper.bin",ce)\nend\n'
 def paper_check(d):
  from fractions import Fraction as F
  x=list(map(F,['5/3','2','4','3','6','5','7']));y=list(map(F,range(1,8)));v=list(map(F,['3/5','1/2','1/4','1/3','1/6','1/5','1/7']));u=list(map(F,['0','1/6','1/8','1/6','1/6','1/5','2/7']));w=list(map(F,['1','2/3','1/4','1/3','-1/6','-1/5','-1']))
  expected=np.full((7,7),np.nan)
  for i,j in itertools.product(range(7),repeat=2):
   num=u[i]*y[j];den=v[i]*x[j]-w[i]
   if den>0:expected[i,j]=float(num/den)
  got=readbin(d/'paper.bin');assert np.allclose(got,expected,equal_nan=True,atol=2e-15,rtol=0);assert abs(got[6,0]-3/13)<1e-15;assert got[2,2]==.5
  # C's ordinary input score is 5/8, whereas translated diagonal is 1/2.
  for k in range(7):
   p=problem(np.array(x,dtype=float)[:,None],np.array(y,dtype=float)[:,None],k,'vrs','input');f=fit(p);assert f.success and abs(-f.fun-float(u[k]*y[k]+w[k]))<1e-12
  return dict(exact_rational_cells=49,primary_LPs=7,diagonal_C=.5,ordinary_self_C=.625,G_to_A=float(got[6,0]))
 run('paper_fixed_weights',body,paper_check)
 # Relative guard must be homogeneous; no arbitrary 0/0 replacement or clipping.
 body=f'run "{pkg}/dea_crosseff.ado"\nmata:\n'
 for scale in ['1e-150','1','1e150']:
  body+=f'a={scale}\nassert(_dea_crosseff_lzratio(2*a,a,a,1e-7)==1)\nassert(_dea_crosseff_lzratio(a,0,a,1e-7)==.)\nassert(_dea_crosseff_lzratio(2*a,0,a,1e-7)==0)\nassert(_dea_crosseff_lzratio(a,1e-12*a,(1-1e-12)*a,1e-7)==.)\nassert(_dea_crosseff_lzratio(a,2*a,0,1e-7)==.)\n'
 body+='assert(_dea_crosseff_lzratio(.,1,0,1e-7)==.)\nassert(_dea_crosseff_lzratio(1,.,0,1e-7)==.)\nassert(_dea_crosseff_lzratio(1,1,.,1e-7)==.)\nassert(_dea_crosseff_lzratio(-1,0,0,1e-7)==.)\nassert(_dea_crosseff_lzratio(1,-1,0,1e-7)==.)\nassert(_dea_crosseff_lzratio(1,1+1e-8,0,1e-7)>1)\nend\n';run('ratio_boundaries',body)
 for key,options in [('bad_value','rts(vrs) vrscross(unknown)'),('crs','rts(crs) vrscross(limzhu)'),('default_crs','vrscross(limzhu)'),('output_vrs','rts(vrs) ort(output) vrscross(limzhu)'),('output_alias','r(v) or(o) vrsc(limzhu)')]:
  for form in ['A','B']:
   cmd='dea_crosseff x=y,' if form=='A' else 'dea_crosseff,inputs(x) outputs(y)'
   run(f'invalid_{key}_{form}',DATA+'gen str12 res="keep"\n'+unchanged(cmd+' '+options+' gen(res) replacevars',198))
 run('reload_state',DATA+f'run "{pkg}/dea_crosseff.ado"\nrun "{pkg}/dea_crosseff.ado"\ndea_crosseff x=y,rts(vrs) vrsc(limzhu) gen(res) matrix\ndea_crosseff x=y,rts(vrs) gen(res) replacevars\nassert "`r(peer_rule)\'"=="ratio_no_intercept"\nassert "`r(vrscross)\'"=="legacy"\ncapture matrix list r(ce_matrix)\nassert _rc==111\n')
 # Compatibility: all old LP selections and outputs remain bit-identical.
 if a.baseline_ado:
  for rts,ori,method in itertools.product(['crs','vrs'],['input','output'],['arbitrary','benevolent','aggressive']):
   cmd=f'dea_crosseff x=y,rts({rts}) ort({ori}) method({method}) matrix'
   body=DATA+f'adopath ++ "{a.baseline_ado.resolve()}"\nrun "{a.baseline_ado.resolve()}/dea_crosseff.ado"\n'+cmd+' gen(old)\n'
   for k in ['ce_matrix','weights_input','weights_output','intercept','stage_status','peer_count']:body+=f'matrix old_{k}=r({k})\n'
   body+=f'adopath ++ "{pkg}"\nrun "{pkg}/dea_crosseff.ado"\n'+cmd+' gen(new) vrscross(legacy)\nassert old==new & old_self==new_self & old_maverick==new_maverick\n'
   for k in ['ce_matrix','weights_input','weights_output','intercept','stage_status','peer_count']:body+=f'mata: assert(all(st_matrix("old_{k}"):==st_matrix("r({k})")))\n'
   body+=cmd+' gen(default)\nassert new==default & new_self==default_self & new_maverick==default_maverick\n'
   run(f'legacy_{rts}_{ori}_{method}',body)
 # Audit real LP submissions for the new option; it must not alter the LP models.
 for method,tol in itertools.product(['arbitrary','benevolent','aggressive'],['1e-4','1e-7','1e-10']):
  body=CAPTURE+'scalar r08calls=0\n'+DATA+f'dea_crosseff x=y,rts(vrs) vrsc(limzhu) method({method}) tol({tol}) matrix\nassert r08calls>0\n'
  run(f'lp_{method}_{tol}',body,lambda d,tol=tol:verify_lp(d,float(tol)),source=audit)
 # Explicit primary failure and secondary fallback with the correction active.
 for fault in ['primary','secondary']:
  target=out/('fault_ado_'+fault);shutil.copytree(pkg,target);p=target/'_dea_nonte_lp.ado';s=p.read_text();cond='rows(S.getEquality())=='+('1' if fault=='primary' else '2');s=s.replace('    c=S.getCoefficients(); n=cols(c)','    if ('+cond+') return(1)\n    c=S.getCoefficients(); n=cols(c)');p.write_text(s)
  body=DATA+'dea_crosseff x=y,rts(vrs) vrsc(limzhu) gen(res) matrix method(benevolent)\n'
  if fault=='primary':body+='assert r(n_failed_primary)==3 & r(n_missing_peer_pairs)==6 & r(n_missing_diagonal)==3\nassert missing(res,res_self,res_maverick)\n'
  else:body+='assert r(n_failed_primary)==0 & r(n_failed_secondary)==3\nassert !missing(res_self)\nmatrix old=r(ce_matrix)\ndea_crosseff x=y,rts(vrs) vrsc(limzhu) gen(res) replacevars matrix\nmata: assert(all(st_matrix("old"):==st_matrix("r(ce_matrix)")))\n'
  run('fault_'+fault,body,source=target)
 if a.rscript and jobs:
  with (out/'r_jobs.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(jobs[0]));w.writeheader();w.writerows(jobs)
  (out/'compare.R').write_text(R_SCRIPT);p=subprocess.run([a.rscript,str(out/'compare.R'),str(out/'r_jobs.csv'),str(out/'r_scores.csv'),str(out/'R_session.txt')],capture_output=True,env=dict(os.environ,RGL_USE_NULL='TRUE'),timeout=300);(out/'R.log').write_bytes(p.stdout+p.stderr);r=dict(test='R_deaR_input_correction',passed=False)
  try:
   assert p.returncode==0,p.stderr.decode()[-2500:];nself=ncells=ndiag=0;maxself=maxcell=0.;selected_gap=0.;different_diag=0;guarded=0
   with (out/'r_scores.csv').open() as f:
    for rr in csv.DictReader(f):
     d=out/rr['test'];k=int(rr['row'])-1;l=int(rr['col'])-1;ce=readbin(d/'ce.bin');scores=readbin(d/'scores.bin');fixed=float(rr['fixed']) if rr['fixed']!='NA' else np.nan
     if l==0:
      delta=abs(scores[k,1]-float(rr['self']));assert delta<3e-6,(rr['test'],k,delta);maxself=max(maxself,delta);nself+=1
     if k==l:
      assert abs(fixed-scores[k,1])<1e-10;ndiag+=1
      if not np.isfinite(ce[k,l]) or abs(ce[k,l]-fixed)>1e-6:different_diag+=1
     elif np.isfinite(ce[k,l]):
      delta=abs(ce[k,l]-fixed);assert delta<3e-6,(rr,ce[k,l]);ncells+=1;maxcell=max(maxcell,delta)
      if rr['independent']!='NA':selected_gap=max(selected_gap,abs(ce[k,l]-float(rr['independent'])))
     else:guarded+=1
   r.update(passed=True,independent_primary_scores=nself,self_max_gap=maxself,same_weight_off_diagonal_cells=ncells,same_weight_max_gap=maxcell,R_diagonals_replaced=ndiag,diagonals_differing_from_translated=different_diag,guarded_missing_cells=guarded,alternative_weight_CE_max_difference=selected_gap)
  except Exception as e:r['error']=repr(e)
  rows.append(r);print(json.dumps(r),flush=True)
 (out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');summary=dict(total=len(rows),passed=sum(r['passed'] for r in rows),failed=sum(not r['passed'] for r in rows));(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary);return bool(summary['failed'])
if __name__=='__main__':raise SystemExit(main())
