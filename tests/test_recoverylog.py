#!/usr/bin/env python3
"""Public recovery context/lifecycle and byte parity; Stata 16+, NumPy/SciPy.
Fault injections are restricted to private copies, never the tested package.
"""
from pathlib import Path
import argparse,hashlib,itertools,json,re,shutil,subprocess,sys,time
import numpy as np
from test_core_retries import WRITER,readbin
from test_simplex_stages import DATA, instrument
REPO=Path(__file__).resolve().parents[1]
BASE=DATA+'gen double b1=1+id/3\ngen double b2=2+mod(id,3)\ngen str30 label=cond(mod(id,2),"same label",char(34)+"quoted"+char(34))\n'
SNAP='''local mats : r(matrices)
local nums : r(scalars)
foreach m of local mats {
    mata: savebin("PREFIX_m_`m'.bin",st_matrix("r(`m')"))
}
foreach n of local nums {
    mata: savebin("PREFIX_s_`n'.bin",st_numscalar("r(`n')"))
}
mata: savebin("PREFIX_scores.bin",st_data(.,"eff"))
'''
def stage_linkage_additions(before, after, prefix):
 """Allow only the three documented new returns; verify their row/count contract."""
 old={p.name for p in before.glob(prefix+'*.bin')};new={p.name for p in after.glob(prefix+'*.bin')}
 suffixes=['m_stage_linkage.bin','s_n_stage_link_adjusted.bin','s_n_stage_link_retries.bin']
 allowed={prefix+'_'+s for s in suffixes};extra=new-old
 assert old<=new and (not extra or extra==allowed),(prefix,old-new,extra)
 if allowed & new:
  assert allowed<=new
  L=readbin(after/(prefix+'_'+suffixes[0]));E=readbin(after/(prefix+'_m_efficiency.bin'))
  assert L.shape==(len(E),30) and np.isfinite(L[:,0]).all()
  assert (L[:,0]==np.floor(L[:,0])).all() and (L[:,0]>0).all() and (np.diff(L[:,0])>0).all()
  assert readbin(after/(prefix+'_'+suffixes[1]))[0,0]==sum(L[:,14]==1)
  assert readbin(after/(prefix+'_'+suffixes[2]))[0,0]==sum(np.isfinite(L[:,9]))
 return sorted(extra)
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',required=True,type=Path);ap.add_argument('--baseline-ado',required=True,type=Path);ap.add_argument('--fixture',type=Path);args=ap.parse_args()
 out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';shutil.copytree(REPO/'ado',pkg);results=[];cases=[]
 def run(name,body,source=pkg):
  d=out/name;d.mkdir();script=f'clear all\nset more off\nset matsize 400\nadopath ++ "{source}"\n'+WRITER+body+'\ndi "TEST_PASS"\nexit, clear\n';(d/'test.do').write_text(script)
  t=time.monotonic();p=subprocess.run([args.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=240);log=(d/'test.log').read_text() if (d/'test.log').exists() else '';ok=p.returncode==0 and 'TEST_PASS' in log.splitlines();results.append(dict(test=name,passed=ok,seconds=time.monotonic()-t));(out/'runs.json').write_text(json.dumps(results,indent=2)+'\n')
  assert ok,(name,p.returncode,log[-4500:]);return d
 def check(d,prefix,model,engine,stage,rts,sample,expect_source=None):
  c=readbin(d/f'{prefix}_m_recovery_stats.bin').ravel();assert c.shape==(8,);assert c[0]==c[6]+c[7] and c[2]<=c[1]<=c[0] and c[3]+c[4]==c[5]
  if c[7]==0:
   assert engine=='ipm' and not (d/f'{prefix}_m_recovery_log.bin').exists();assert (c==0).all();return dict(lp_calls=0,attempted=0,recovered=0)
  a=readbin(d/f'{prefix}_m_recovery_log.bin');assert a.shape==(int(c[7]),25)
  assert np.array_equal(a[:,6],np.arange(1,len(a)+1));np.testing.assert_array_equal(a[:,11]+a[:,12],a[:,13]);np.testing.assert_array_equal(a[:,[7,8,11,12,13]].sum(0),c[1:6])
  assert np.isin(a[:,0],sample).all();assert np.isfinite(a[:,19:21]).all();assert (a[:,17:19]>0).all()
  score=readbin(d/f'{prefix}_scores.bin').ravel();np.testing.assert_array_equal(a[:,23],score[a[:,0].astype(int)-1]);np.testing.assert_array_equal(a[:,22],np.isfinite(a[:,23]));assert np.isin(a[:,24],[0,1,2,3]).all()
  if expect_source is not None: assert (a[a[:,22]==1,24]==expect_source).all()
  if model=='te':
   assert (a[:,2]==1).all() and np.isin(a[:,4],[1,2] if stage==2 else [1]).all();assert np.isfinite(a[:,21]).all()
   if rts=='vrs':assert set(a[:,1])=={1,2} and set(a[a[:,1]==2,3])=={1,2,3} and (a[a[:,1]==2,4]==1).all()
   else:assert (a[:,1]==1).all()
  elif model=='sbm': assert (a[:,2]==2).all() and (a[:,4:6]==1).all()
  elif model=='supersbm':assert np.isin(a[:,2],[2,3]).all() and (a[:,4:6]==1).all()
  else:assert (a[:,2]==4).all() and np.isin(a[:,5],[1,2,3,4]).all()
  return dict(lp_calls=int(c[0]),attempted=int(c[1]),recovered=int(c[2]))
 # Full combinations for TE; SBM orientations and RTS; weak directions and RTS.
 for model in ['te','sbm','supersbm','weak']:
  for rts,ori in itertools.product(['crs','vrs','drs','irs'],(['input','output'] if model=='te' else ['none','input','output'] if model in ['sbm','supersbm'] else ['x','yb','xyb'])):
   for stage in ([1,2] if model=='te' else [2]):
    for engine in ['simplex','auto','ipm']:
     i=len(cases);qual='if mod(id,3)!=0 in 2/8' if (i//3)%3==0 else '';sample=[x for x in range(1,9) if (2<=x<=8 and x%3!=0)] if qual else list(range(1,9))
     io=f'x1 x2=y1 y2 {qual},' if i%2==0 else f'y1 y2 {qual}, inputs(x1 x2)'
     opts=f'model({model}) '+('' if ori=='none' else f'ort({ori}) ') if model!='weak' else f'model(undesirable) bad(b1 b2) disp(weak) dir({ori}) weakmatrix '
     cmd=f'dea {io} dmu(label) engine({engine}) rts({rts}) {opts} stage({stage}) gen(eff)'
     cases.append(dict(case=f'c{i:03}',model=model,engine=engine,stage=stage,rts=rts,sample=sample,cmd=cmd))
 # All commands are grouped into a single session per mode, exercising repeated lifecycles.
 for mode,source in [('baseline',args.baseline_ado.resolve()),('off',pkg),('on',pkg)]:
  body=''
  for x in cases:
   body+='clear\n'+BASE+x['cmd']+(' recoveryl' if mode=='on' else '')+'\n'+SNAP.replace('PREFIX',x['case'])
   if mode!='baseline':body+='mata: assert(_dea_rec_on==0)\n'
  run('grid_'+mode,body,source)
 changes=[];counts=[];parity=0
 for x in cases:
  prefix=x['case'];b=out/'grid_baseline';off=out/'grid_off';on=out/'grid_on'
  files=sorted(p.name for p in b.glob(prefix+'*.bin'))
  added=stage_linkage_additions(b,off,prefix)
  for name in added:assert (off/name).read_bytes()==(on/name).read_bytes(),(prefix,name)
  for name in files:
   assert (b/name).read_bytes()==(off/name).read_bytes()==(on/name).read_bytes(),(prefix,name);parity+=1
  counts.append(dict(case=prefix,**check(on,prefix,**{k:x[k] for k in ['model','engine','stage','rts','sample']})))
 # Recovery case from the sealed weak-disposal investigation.
 if args.fixture:
  body=f'import delimited using "{args.fixture.resolve()}", clear asdouble\ndea x1 x2=y1 y2, model(undesirable) bad(b1 b2) disp(weak) dir(x) rts(vrs) dmu(id) engine(simplex) weakmatrix gen(eff) recoverylog\n'+SNAP.replace('PREFIX','target')
  d=run('actual_recovery',body);c=check(d,'target','weak','simplex',2,'vrs',list(range(1,readbin(d/'target_scores.bin').shape[0]+1)));assert c['attempted']>0 and c['recovered']>0;assert np.isfinite(readbin(d/'target_scores.bin')).all();counts.append(dict(case='actual_recovery',**c))
 # Errors, default-off, nonnumeric DMU labels, internal raw-ledger coexistence and cap.
 body=BASE+'''run "PKG/_dea_simplex.ado"
mata: _rsm_audit_reset(1)
dea x1 x2=y1 y2, dmu(label) engine(simplex) recoverylog gen(eff)
matrix C=r(recovery_stats)
mata: assert(_rsm_audit_counts()[1]==st_matrix("C")[1] & _dea_rec_on==0)
drop eff
mata: _rsm_audit_reset(0)
'''.replace('PKG',str(pkg))
 for i,extra in enumerate(['model(supereff)','fdh','bootstrap reps1(60)','bounded','contractcheck','unknownoption','tolerance(-1)','saving("missing/subdir/result.dta")']):
  body+=f'capture noisily dea x1 x2=y1 y2, dmu(label) engine(simplex) recoverylog gen(eff) {extra}\nassert _rc!=0\nmata: assert(_dea_rec_on==0 & rows(_dea_rec_rows)==0)\ncapture matrix list r(recovery_log)\nassert _rc==111\ncapture drop eff\nquietly dea x1 x2=y1 y2, dmu(label) engine(simplex) gen(eff)\nassert r(N)==8\ncapture matrix list r(recovery_stats)\nassert _rc==111\ndrop eff\n'
 body+='''di "QUIET_START"
quietly dea x1 x2=y1 y2, dmu(label) engine(simplex) recoverylog gen(eff)
di "QUIET_END"
assert r(recovery_stats)[1,1]>0
drop eff
replace x1=-1 in 1
capture noisily dea x1 x2=y1 y2, model(sbm) engine(simplex) recoverylog gen(eff)
assert _rc==459
mata: assert(_dea_rec_on==0 & rows(_dea_rec_rows)==0)
replace x1=1 in 1
set matsize 40
expand 4
replace id=_n
dea x1 x2=y1 y2, dmu(id) engine(simplex) recoverylog gen(eff)
assert r(recovery_stats)[1,1]==64 & r(recovery_stats)[1,7]==24 & r(recovery_stats)[1,8]==40
assert rowsof(r(recovery_log))==40
mata: assert(_dea_rec_on==0)
clear mata
capture drop eff
dea x1 x2=y1 y2, dmu(id) engine(simplex) recoverylog gen(eff)
assert r(recovery_stats)[1,1]==64
mata: assert(_dea_rec_on==0)
'''
 d=run('lifecycle',body);log=(d/'test.log').read_text();quiet=log.split('\nQUIET_START\n')[1].split('\nQUIET_END\n')[0];assert 'Data Envelopment' not in quiet and 'DMU' not in quiet
 # Existing fault hooks applied to a private source copy.
 faulty=out/'fault_ado';instrument(pkg,faulty)
 init='scalar audit_s1=0\nscalar audit_s2=0\nscalar audit_dual=0\nscalar fault_primal=0\nscalar fault_dual=0\nscalar fault_stage2=0\n'
 for name,settings,eng,checks in [
  ('stage2','scalar fault_stage2=1','simplex','assert all(A[,21]:==1) & all(A[,22]:==0) & all(A[,23]:==1)\nassert any(A[,20]:==3)'),
  ('dual','scalar fault_primal=1','simplex','assert any(A[,6]:==2) & all(A[,23]:==1)'),
  ('missing','scalar fault_primal=1\nscalar fault_dual=1','simplex','assert all(A[,21]:==0) & all(A[,23]:==0) & all(A[,25]:==0)'),
  ('auto_native','scalar fault_primal=1\nscalar fault_dual=1','auto','assert all(A[,21]:==0) & all(A[,23]:==1) & all(A[,25]:==2)')]:
  checks='\n'.join('assert('+line[7:]+')' if line.startswith('assert ') else line for line in checks.splitlines())
  body=BASE+init+settings+f'\ndea x1 x2=y1 y2, dmu(id) engine({eng}) recoverylog gen(eff)\n'+SNAP.replace('PREFIX','fault')+'mata:\nA=st_matrix("r(recovery_log)")\n'+checks+'\nend\n'
  d=run('fault_'+name,body,faulty);check(d,'fault','te',eng,2,'crs',list(range(1,9)))
 # Corrupt a candidate before the non-TE original-constraint check: core success is not model success.
 reject=out/'reject_ado';shutil.copytree(pkg,reject);p=reject/'_dea_simplex.ado';s=p.read_text();a=s.index('struct RSM_Result scalar _dea_nonte_simplex(');b=s.index('// #',a);part=s[a:b];anchor='\n    if (res.status == 0) {';assert anchor in part;part=part.replace(anchor,anchor+'\n        res.x=J(rows(res.x),1,0);',1);p.write_text(s[:a]+part+s[b:])
 body=BASE+'dea x1 x2=y1 y2, dmu(id) model(sbm) engine(simplex) recoverylog gen(eff)\n'+SNAP.replace('PREFIX','reject')+'mata:\nA=st_matrix("r(recovery_log)")\nassert(all(A[,15]:==0) & all(A[,20]:==4) & all(A[,21]:==0) & all(A[,23]:==0))\nend\n'
 run('fault_original_residual',body,reject)
 summary=dict(grid_commands=len(cases),grid_executions=3*len(cases),parity_numeric_files=parity,runs=len(results),all_passed=True,lp_calls=sum(c['lp_calls'] for c in counts),attempted=sum(c['attempted'] for c in counts),recovered=sum(c['recovered'] for c in counts),fault_scenarios=5)
 (out/'grid.json').write_text(json.dumps(cases,indent=2)+'\n');(out/'counts.json').write_text(json.dumps(counts,indent=2)+'\n');(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
