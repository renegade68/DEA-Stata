#!/usr/bin/env python3
"""Explicit IPM projection recovery: unchanged scores, declared RHS, saved metadata.
Requires immediately preceding dea 2.0.48 snapshot and NumPy/SciPy/Stata 16+.
"""
from _support import package_files
from pathlib import Path
import argparse,itertools,json,hashlib,subprocess,shutil,re,time
import numpy as np
from scipy.optimize import linprog
from test_core_retries import WRITER,readbin,same_bits,check_trace
from test_ipm_stages import FIXTURES
from test_lp_residuals import residual
REPO=Path(__file__).resolve().parents[1]
BASE=FIXTURES['p03_vrs']


def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--baseline-ado',type=Path,required=True);ap.add_argument('--only');args=ap.parse_args()
 out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir();hashes={}
 for n in package_files(REPO):
  raw=(REPO/n).read_bytes();(pkg/Path(n).name).write_bytes(raw);hashes[n]=hashlib.sha256(raw).hexdigest()
 (out/'manifest.json').write_text(json.dumps(hashes,indent=2)+'\n');results=[]
 baseline=args.baseline_ado.resolve();(out/'baseline_manifest.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in baseline.iterdir() if p.is_file()},indent=2)+'\n')
 def run(name,body,source=pkg,check=None):
  if args.only and not re.search(args.only,name):return None
  d=out/name;d.mkdir();s='version 16\nclear all\nset more off\nset linesize 160\n'
  for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/key).mkdir();s+=f'sysdir set {key} "{d/key}"\n'
  s+=f'adopath ++ "{source}"\n'+WRITER+body+'\ndi "TEST_PASS"\n';(d/'test.do').write_text(s);r=dict(test=name,passed=False);start=time.monotonic()
  try:
   p=subprocess.run([args.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=240)
   log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else ''
   assert p.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-3000:] or p.stderr.decode(errors='replace')
   if check:r.update(check(d))
   r['passed']=True
  except Exception as e:r['error']=repr(e)
  r['seconds']=round(time.monotonic()-start,3);results.append(r);(out/'results.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(r),flush=True);return d
 def save(mats,suffix):return ''.join(f'mata: savebin("{m}{suffix}.bin",st_matrix("r({m})"))\n' for m in mats)
 mats=['efficiency','lambda','slack_in','slack_out','stage_status','ipm_trace']
 for fixture,rts,ori,tol in itertools.product(FIXTURES,['crs','vrs','drs','irs'],['input','output'],[1e-7,1e-10]):
  name=f'grid_{fixture}_{rts}_{ori}_{tol}'
  if args.only and not re.search(args.only,name):continue
  q='if id>=2 in 2/8' if fixture=='duplicate' else ''
  cmd=f'dea x1 x2=y1 y2 {q}, dmu(id) engine(ipm) rts({rts}) ort({ori}) tolerance({tol}) stage(2) reference'
  old=run(name+'_before',FIXTURES[fixture]+cmd+' gen(e)\n'+save(mats,''),baseline)
  s=FIXTURES[fixture]+cmd+' gen(e0)\n'+save(mats,'0')
  s+=cmd+' gen(ez) radialtolerance(0) replacevars\nassert ez==e0\n'+save(mats,'z')
  s+='capture matrix list r(projection)\nassert _rc==111\n'
  s+=f'dea y1 y2 {q}, inputs(x1 x2) dmu(id) engine(ipm) rts({rts}) ort({ori}) tolerance({tol}) stage(2) reference gen(e) radialtolerance(1e-7) replacevars\n'
  s+='assert e==e0\n'+save(mats+['projection','projection_trace'],'1')
  s+='matrix C=(r(n_projection_attempted),r(n_projection_adjusted),r(n_projection_solves),r(n_failed),r(n_stage2_failed))\nmata: savebin("counts.bin",st_matrix("C"))\n'
  s+=f'gen byte selected=1\n'+('replace selected=inrange(_n,2,8) & id>=2\n' if q else '')
  s+='assert missing(_dea_radial,_dea_radial2,_dea_radial_adj) if !selected\n'
  s+='mata: savebin("data.bin",st_data(.,tokens("x1 x2 y1 y2"),"selected"))\n'
  s+='mata: savebin("stored.bin",st_data(.,tokens("_dea_radial _dea_radial2 _dea_radial_adj"),"selected"))\n'
  def check(d,old=old,ori=ori,rts=rts,tol=tol):
   for m in mats:
    assert (d/f'{m}0.bin').read_bytes()==(old/f'{m}.bin').read_bytes(),m
    assert (d/f'{m}z.bin').read_bytes()==(old/f'{m}.bin').read_bytes(),m
   assert (d/'efficiency1.bin').read_bytes()==(old/'efficiency.bin').read_bytes()
   assert (d/'ipm_trace1.bin').read_bytes()==(old/'ipm_trace.bin').read_bytes()
   p=readbin(d/'projection1.bin');h=readbin(d/'projection_trace1.bin');status=readbin(d/'stage_status1.bin');oldstatus=readbin(old/'stage_status.bin');counts=readbin(d/'counts.bin')[0]
   np.testing.assert_array_equal(p[:,:3],readbin(d/'stored.bin'))
   np.testing.assert_array_equal(p[:,4],oldstatus[:,1]);np.testing.assert_array_equal(p[:,6],status[:,1])
   use=p[:,7]==1;tried=np.isfinite(p[:,5]);basegood=oldstatus[:,1]==0
   assert np.all(~use|tried) and not np.any(basegood&tried)
   assert np.all(np.isnan(p[~tried,3:4])) and np.isnan(h[~tried]).all()
   eff=readbin(d/'efficiency1.bin')[:,0];rad=eff if ori=='input' else 1/eff
   np.testing.assert_array_equal(rad,p[:,0]);np.testing.assert_array_equal(p[basegood,1],p[basegood,0]);assert np.all(p[basegood,2]==0)
   assert np.isnan(p[status[:,1]!=0,1:3]).all()
   if tried.any():
    expected=rad+1e-7*(1+abs(rad)) if ori=='input' else np.maximum(1,rad-1e-7*(1+abs(rad)))
    np.testing.assert_array_equal(p[tried,3],expected[tried]);np.testing.assert_array_equal(p[use,1],p[use,3])
    np.testing.assert_array_equal(p[use,2],p[use,1]-p[use,0])
   hc=check_trace(np.c_[h,np.full((len(h),50),np.nan)],tol)
   np.testing.assert_array_equal(counts[:3],[tried.sum(),use.sum(),tried.sum()+hc[0]])
   assert counts[3]==counts[4]==np.count_nonzero(status[:,1]!=0)
   for m in ['lambda','slack_in','slack_out']:same_bits(readbin(d/f'{m}1.bin')[basegood],readbin(old/f'{m}.bin')[basegood])
   data=readbin(d/'data.bin');x=data[:,:2];y=data[:,2:];xs=abs(x).max(axis=0);ys=abs(y).max(axis=0);x=x/xs;y=y/ys;n=len(x);L=readbin(d/'lambda1.bin');si=readbin(d/'slack_in1.bin')/xs;so=readbin(d/'slack_out1.bin')/ys
   maxres=maxgap=0.;nonoptimal=0
   for k in np.where(use)[0]:
    ae=np.block([[x.T,np.eye(2),np.zeros((2,2))],[y.T,np.zeros((2,2)),-np.eye(2)]]);b=np.r_[p[k,1]*x[k] if ori=='input' else x[k],y[k] if ori=='input' else p[k,1]*y[k]];au=np.empty((0,n+4));bu=np.empty(0);r=np.r_[np.ones(n),np.zeros(4)]
    if rts=='vrs':ae=np.vstack([ae,r]);b=np.r_[b,1]
    if rts in ['irs','drs']:au=(r*(1 if rts=='drs' else -1))[None,:];bu=np.array([1 if rts=='drs' else -1])
    z=np.r_[L[k],si[k],so[k]];diag=residual(z,ae,b,au,bu,np.zeros(n+4),np.full(n+4,np.inf));maxres=max(maxres,max(diag[3:]));assert max(diag[3:])<=tol
    c=np.r_[np.zeros(n),-np.ones(4)];f=linprog(c,A_eq=ae,b_eq=b,A_ub=au if len(au) else None,b_ub=bu if len(bu) else None,bounds=(0,None),method='highs-ds',options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    if f.success:maxgap=max(maxgap,abs(c@z-f.fun));assert abs(c@z-f.fun)<1e-6
    else:nonoptimal+=1
   return dict(dmus=n,unchanged_initial_projection_rows=int(basegood.sum()),attempted=int(tried.sum()),adjusted=int(use.sum()),unavailable_before=int((oldstatus[:,1]!=0).sum()),unavailable_after=int((status[:,1]!=0).sum()),max_scaled_residual=float(maxres),max_slack_gap=float(maxgap),independent_nonoptimal=nonoptimal)
  run(name+'_after',s,check=check)
 # Rejection occurs before replacing public outputs or generating automatic IDs.
 for form,value in itertools.product(['equal','standard'],['-1','.','.a','.z','1','1e308','abc']):
  io='x1 x2=y1 y2,' if form=='equal' else 'y1 y2, inputs(x1 x2)'
  s=BASE+'gen double keepme=42\ngen double _dea_radial=71\nquietly datasignature\nlocal sig "`r(datasignature)\'"\n'
  s+=f'capture noisily dea {io} engine(ipm) radialtolerance({value}) gen(keepme) replacevars\nassert _rc==198\n'
  s+='quietly datasignature\nassert "`r(datasignature)\'"=="`sig\'"\ncapture confirm variable _dmu\nassert _rc==111\n'
  run(f'invalid_{form}_{value}',s)
 for scope in ['engine(simplex)','engine(auto)','engine(ipm) stage(1)','engine(ipm) model(supereff) rts(vrs)','engine(ipm) model(sbm)','engine(ipm) fdh','engine(ipm) xlo(x1 x2) xhi(x1 x2)']:
  s=BASE+'gen double keepme=42\ncapture noisily dea x1 x2=y1 y2, '+scope+' radialtolerance(1e-7) gen(keepme) replacevars\nassert _rc==198\nassert keepme==42\n'
  run('scope_'+str(len(results)),s)
 for fault in ['negative_input','negative_output','zero_inputs','zero_outputs','missing_excluded','invalid_outside']:
  changes={'negative_input':'replace x1=-1 in 1','negative_output':'replace y1=-1 in 1','zero_inputs':'replace x1=0 in 1\nreplace x2=0 in 1','zero_outputs':'replace y1=0 in 1\nreplace y2=0 in 1','missing_excluded':'replace x1=. in 1','invalid_outside':'replace x1=-1 in 1'}[fault]
  q='if id>1' if fault=='invalid_outside' else ''
  s=BASE+changes+'\ngen double keepme=42\n'+f'capture noisily dea x1 x2=y1 y2 {q}, engine(ipm) rts(irs) radialtolerance(1e-7) gen(keepme) replacevars\n'
  s+=('assert _rc==0\nassert missing(keepme[1])\n' if fault in ['missing_excluded','invalid_outside'] else 'assert _rc==459\nassert keepme==42\n')
  run(f'quantity_{fault}',s)
 for ori,rts,mode in itertools.product(['input','output'],['irs','vrs'],['saving','bootstrap','contractcheck','no_reference']):
  opts={'saving':'reference saveref saving(result)','bootstrap':'reference bootstrap reps1(50)','contractcheck':'reference contractcheck','no_reference':''}[mode]
  cmd=f'dea x1 x2=y1 y2, dmu(id) engine(ipm) rts({rts}) ort({ori}) tolerance(1e-10) radialtolerance(1e-7)'
  s=BASE+cmd+' gen(e0) reference\nmatrix P0=r(projection)\nmatrix H0=r(projection_trace)\n'
  s+=cmd+f' gen(e) {opts} replacevars\nassert e==e0\nmatrix P=r(projection)\nmatrix H=r(projection_trace)\n'
  s+='mata: assert(all(st_matrix("P"):==st_matrix("P0")) & all(st_matrix("H"):==st_matrix("H0")))\n'
  if mode=='saving':s+='local refs "`r(refvars)\'"\nmatrix L=r(lambda)\nuse result,clear\nmata: assert(all(st_data(.,tokens("_dea_radial _dea_radial2 _dea_radial_adj")):==st_matrix("P")[.,1..3]))\nmata: assert(all(st_data(.,tokens(st_local("refs"))):==st_matrix("L")))\n'
  run(f'preserve_{ori}_{rts}_{mode}',s)
 for collision in ['_dea_radial','_dea_radial2','_dea_radial_adj','output_alias','input_alias']:
  if collision in ['output_alias','input_alias']:
   extra='gen(_dea_radial)' if collision=='output_alias' else 'inputs(_dea_radial x2)'
   body=('dea x1 x2=y1 y2,' if collision=='output_alias' else 'dea y1 y2,')
   setup='' if collision=='output_alias' else 'gen double _dea_radial=x1\n'
  else:setup=f'gen double {collision}=73\n';body='dea x1 x2=y1 y2,';extra=''
  s=BASE+setup+f'capture noisily {body} engine(ipm) radialtolerance(1e-7) {extra}\nassert _rc==110\n'
  if collision not in ['output_alias','input_alias']:s+=f'assert {collision}==73\n'
  run('collision_'+collision,s)
 for band in ['1e-300','0.1']:
  s=BASE+f'dea x1 x2=y1 y2, engine(ipm) rts(irs) tolerance(1e-10) radialtolerance({band}) reference\nmatrix P=r(projection)\n'
  if band=='1e-300':s+='assert r(n_projection_attempted)==0\nassert r(n_projection_adjusted)==0\nassert r(n_stage2_failed)==1\n'
  run('band_'+band,s)
 # A shifted RHS that is deliberately made infeasible must remain unavailable.
 target=out/'failed_shift_ado';shutil.copytree(pkg,target);p=target/'_dea_ipm.ado';s=p.read_text();a='    rc = S2.optimize()';assert s.count(a)==1;s=s.replace(a,'    if (args()>=12) { b_eq[1]=-1;S2.setEquality(A_eq,b_eq); }\n'+a);p.write_text(s)
 run('shift_still_fails',BASE+'''dea x1 x2=y1 y2, dmu(id) engine(ipm) rts(irs) tolerance(1e-10) radialtolerance(1e-7) reference
assert r(n_projection_attempted)==1 & r(n_projection_adjusted)==0 & r(n_stage2_failed)==1
matrix P=r(projection)
matrix L=r(lambda)
assert !missing(P[3,4]) & P[3,6]!=0 & missing(P[3,2],P[3,3])
mata: assert(all(st_matrix("L")[3,.]:>=.))
''',target)
 # Contradictory/unbounded fixed-face evidence must prevent radial movement.
 target=out/'unbounded_fixed_ado';shutil.copytree(pkg,target);p=target/'_dea_ipm.ado';s=p.read_text()
 a='status=_dea_ipm_finish(status,S2.errorcode(),tol,S2.iterations(),';assert s.count(a)==1
 s=s.replace(a,'status=_dea_ipm_finish(2,2,tol,S2.iterations(),');p.write_text(s)
 for ori in ['input','output']:
  run('unbounded_fixed_'+ori,BASE+f"dea x1 x2=y1 y2, dmu(id) engine(ipm) rts(irs) ort({ori}) tolerance(1e-10) radialtolerance(1e-7) reference\n"+"""assert r(n_projection_attempted)==0 & r(n_projection_adjusted)==0
assert r(n_stage2_failed)==3
matrix P=r(projection)
mata: assert(all(st_matrix("P")[.,2..4]:>=.))
""",target)
 from test_ipm_stages import instrument, INIT
 target=out/'dual_only_ado';instrument(pkg,target)
 for ori in ['input','output']:
  run('dual_only_'+ori,INIT+'scalar fault_primal=1\n'+BASE+f"dea x1 x2=y1 y2, dmu(id) engine(ipm) rts(irs) ort({ori}) tolerance(1e-10) radialtolerance(1e-7) reference\n"+"""assert r(n_projection_attempted)==0 & r(n_projection_adjusted)==0
assert r(n_radial_failed)==0 & r(n_stage2_failed)==3
assert audit_s2==0
""",target)
 for form in ['equal','standard']:
  io='x1 x2=y1 y2,' if form=='equal' else 'y1 y2, inputs(x1 x2)'
  run('abbreviation_'+form,BASE+f"dea {io} method(ipm) rts(irs) tol(1e-10) radialt(1e-7)\nassert r(radial_tolerance)==1e-7 & r(n_projection_adjusted)==1\n")
 run('individual_zeros',BASE+"""replace x1=0 in 1
replace y1=0 in 2
dea x1 x2=y1 y2, dmu(id) engine(ipm) radialtolerance(1e-7)
assert rowsof(r(projection))==3
""")
 run('longer_existing_name',BASE+"""gen double _dea_radial_long=713
dea x1 x2=y1 y2, dmu(id) engine(ipm) radialtolerance(1e-7)
assert _dea_radial_long==713
confirm variable _dea_radial, exact
""")
 # A previously loaded engine must be refreshed by the existing loader path.
 run('cached_before_engine',f'run "{baseline/"_dea_ipm.ado"}"\n'+BASE+'''dea x1 x2=y1 y2, dmu(id) engine(ipm) rts(irs) tolerance(1e-10) radialtolerance(1e-7)
assert colsof(r(projection))==8 & r(n_projection_adjusted)==1
''')
 summary=dict(cases=len(results),failed=[r['test'] for r in results if not r['passed']])
 for key in ['dmus','unchanged_initial_projection_rows','attempted','adjusted','unavailable_before','unavailable_after','independent_nonoptimal']:summary[key]=sum(r.get(key,0) for r in results)
 (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');assert not summary['failed'],summary

if __name__=='__main__':main()
