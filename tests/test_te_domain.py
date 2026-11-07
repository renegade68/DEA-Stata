#!/usr/bin/env python3
"""Ordinary TE quantity domain, atomic failure, sample and numeric parity.
Independent radial LPs check allowed zero components/columns. Stata 16+, NumPy/SciPy.
"""
from pathlib import Path
import argparse,itertools,json,subprocess,shutil,hashlib,time
import numpy as np
from test_recoverylog import BASE,SNAP,stage_linkage_additions
from test_core_retries import WRITER,readbin
from verify_bootstrap_lp import solve
REPO=Path(__file__).resolve().parents[1]
DATA=BASE+'recast double x1 x2 y1 y2\n'

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',required=True,type=Path);ap.add_argument('--baseline-ado',required=True,type=Path);a=ap.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';shutil.copytree(REPO/'ado',pkg);runs=[]
 def run(name,body,source=pkg):
  d=out/name;d.mkdir();(d/'test.do').write_text(f'clear all\nset more off\nset linesize 160\nadopath ++ "{source}"\n'+WRITER+body+'\ndi "TEST_PASS"\nexit,clear\n');t=time.monotonic();p=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=300);log=(d/'test.log').read_text() if (d/'test.log').exists() else '';ok=p.returncode==0 and 'TEST_PASS' in log.splitlines();runs.append(dict(test=name,passed=ok,seconds=time.monotonic()-t));(out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n');assert ok,(name,log[-4500:],p.stderr.decode(errors='replace'));return d
 changes={'negative_input':'replace x1=-1e-250 in 1','negative_output':'replace y2=-1 in 1','zero_inputs':'replace x1=0 in 1\nreplace x2=0 in 1','zero_outputs':'replace y1=0 in 1\nreplace y2=0 in 1'}
 body='file open guard using "protected.dta", write replace\nfile write guard "DO NOT OVERWRITE"\nfile close guard\n';invalid=[]
 for engine,rts,ori,stage,form,rec,fault in itertools.product(['simplex','ipm','auto'],['crs','vrs','nirs','ndrs'],['input','output'],[1,2],['A','B'],['','recoverylog'],changes):
  rr={'crs':'ccr','vrs':'bcc','nirs':'drs','ndrs':'irs'}[rts] if form=='B' else rts;io='x1 x2=y1 y2,' if form=='A' else 'y1 y2, inputs(x1 x2)';tol='1e-3' if len(invalid)%2 else '1e-12'
  body+='clear\n'+DATA+changes[fault]+'\ngen double eff=42\nlabel var eff "keep result"\nchar eff[review] "preserve metadata"\nquietly datasignature\nlocal sig "`r(datasignature)\'"\n'
  body+=f'capture noisily dea {io} dmu(label) engine({engine}) rts({rr}) ort({ori}) stage({stage}) tol({tol}) {rec} gen(eff) replacevars reference saveref saving("protected.dta") replace\n'
  body+='assert _rc==459\nquietly datasignature\nassert "`r(datasignature)\'"=="`sig\'"\nassert "`:variable label eff\'"=="keep result"\nassert "`:char eff[review]\'"=="preserve metadata"\ncapture confirm variable _dmu\nassert _rc==111\ncapture matrix list r(recovery_log)\nassert _rc==111\n'
  # The helper may not have been loaded in an initial native OFF call.
  if rec:body+='mata: assert(_dea_rec_on==0 & rows(_dea_rec_rows)==0)\n'
  invalid.append(dict(engine=engine,rts=rr,ori=ori,stage=stage,form=form,recording=bool(rec),fault=fault))
 d=run('invalid_grid',body);assert (d/'protected.dta').read_bytes()==b'DO NOT OVERWRITE'
 fixtures={
  'positive':('', ''),
  'zero_cells':('replace x1=0 in 1\nreplace y2=0 in 3\n',''),
  'zero_columns':('replace x2=0\nreplace y2=0\n',''),
  'missing':('replace x1=.a in 1\nreplace y2=.z in 3\n',''),
  'invalid_outside':('replace x1=-1 in 1\nreplace x1=0 in 8\nreplace x2=0 in 8\n','if id>1 in 2/7'),
  'missing_hides_invalid':('replace x1=. in 1\nreplace y1=-1 in 1\nreplace x1=0 in 3\nreplace x2=0 in 3\nreplace y2=.a in 3\n',''),
  'shuffled':('gsort -id\nreplace x1=. in 4\n','if mod(id,2) in 2/8')}
 cases=[]
 for fixture,engine,rts,ori,stage in itertools.product(fixtures,['simplex','ipm','auto'],['crs','vrs','drs','irs'],['input','output'],[1,2]):
  i=len(cases);change,qual=fixtures[fixture];io=f'x1 x2=y1 y2 {qual},' if i%2==0 else f'y1 y2 {qual}, inputs(x1 x2)'
  cases.append(dict(prefix=f'c{i:03}',fixture=fixture,engine=engine,rts=rts,ori=ori,stage=stage,change=change,command=f'dea {io} dmu(label) engine({engine}) rts({rts}) ort({ori}) stage({stage}) gen(eff)'))
 for mode,source in [('baseline',a.baseline_ado.resolve()),('off',pkg),('on',pkg)]:
  body=''
  for x in cases:
   prefix=x['prefix'];body+='clear\n'+DATA+x['change']+x['command']+(' recoverylog' if mode=='on' else '')+'\n'+SNAP.replace('PREFIX',prefix)
   body+='assert r(N)>0\n'
   body+=f'mata: savebin("{prefix}_data.bin",st_data(.,tokens("x1 x2 y1 y2")))\n'
   if mode!='baseline':body+='mata: assert(_dea_rec_on==0)\n'
  run('valid_'+mode,body,source)
 numeric=oracle=0;maxgap=0.
 for x in cases:
  prefix=x['prefix'];base=out/'valid_baseline';off=out/'valid_off';on=out/'valid_on';files=sorted(p.name for p in base.glob(prefix+'*.bin'))
  added=stage_linkage_additions(base,off,prefix)
  for name in added:assert (off/name).read_bytes()==(on/name).read_bytes(),(prefix,name)
  for name in files:assert (base/name).read_bytes()==(off/name).read_bytes()==(on/name).read_bytes(),(prefix,name);numeric+=1
  if x['fixture'] in ['zero_cells','zero_columns']:
   z=readbin(on/(prefix+'_data.bin'));e=readbin(on/(prefix+'_scores.bin')).ravel();assert np.isfinite(e).all()
   for k in range(len(z)):
    v=solve(z[k,:2],z[k,2:],z[:,:2],z[:,2:],x['ori'],x['rts']);gap=abs(v-e[k]);assert gap<2e-6,(prefix,k,v,e[k]);maxgap=max(maxgap,gap);oracle+=1
 # A native positive allowance retains the same domain rule and original sample.
 body=DATA+'gen double eff=42\nreplace x1=-1 in 1\ncapture noisily dea x1 x2=y1 y2, engine(ipm) radialtolerance(1e-7) gen(eff) replacevars recoverylog\nassert _rc==459\nassert eff==42\nmata: assert(_dea_rec_on==0 & rows(_dea_rec_rows)==0)\nreplace x1=1 in 1\n'
 body+='dea x1 x2=y1 y2, dmu(id) engine(ipm) radialtolerance(1e-7) gen(eff) replacevars recoverylog\nassert r(N)==8\n'
 body+='capture noisily dea x1 x2=y1 y2 if id>99, dmu(id) gen(eff) replacevars recoverylog\nassert _rc==2000\nmata: assert(_dea_rec_on==0 & rows(_dea_rec_rows)==0)\n'
 body+='replace x1=. in 1/8\ncapture noisily dea x1 x2=y1 y2, gen(eff) replacevars\nassert _rc==2000\n'
 run('allowance_and_empty',body)
 # No automatic ID or destination is created by rejected quantities.
 body=DATA+'replace y1=0 in 1\nreplace y2=0 in 1\ncapture noisily dea x1 x2=y1 y2, gen(eff) recoverylog saving("new.dta")\nassert _rc==459\ncapture confirm variable _dmu\nassert _rc==111\ncapture confirm variable eff\nassert _rc==111\ncapture confirm file "new.dta"\nassert _rc==601\n'
 run('atomic_new_outputs',body)
 # Strict sign check is independent of scale; max-based row checks do not sum quantities.
 body=''
 for val in ['1e-300','8e307']:
  body+='clear\n'+DATA+f'replace x1={val}\nreplace x2={val}\nreplace y1={val}\nreplace y2={val}\n'
  body+='dea x1 x2=y1 y2, dmu(id) engine(simplex) gen(eff) recoverylog\nassert r(N)==8 & r(n_radial_failed)==0\nassert abs(eff-1)<1e-10\n'
 for ori in ['input','output']:
  body+='clear\ninput double id x y\n'+('1 1 1\n2 1e11 1\n' if ori=='input' else '1 1 1e11\n2 1 1\n')+'end\n'
  body+=f'dea x=y, dmu(id) engine(simplex) ort({ori}) tol(1e-12) gen(eff) recoverylog\nassert eff[2]>0 & abs(eff[2]-1e-11)<1e-18\n'
 run('extremes_and_small_scores',body)
 # Repeated ON/OFF, command errors, clear mata, and default/non-default stages.
 body=DATA+'dea x1 x2=y1 y2, dmu(id) gen(eff) recoverylog\nassert "`r(recovery_scope)\'"=="package_simplex_calls"\nmata: assert(_dea_rec_on==0)\nclear mata\ndea x1 x2=y1 y2, dmu(id) gen(eff) replacevars\nassert "`r(recovery_scope)\'"==""\ncapture matrix list r(recovery_log)\nassert _rc==111\nreplace x1=-1 in 1\ncapture noisily dea x1 x2=y1 y2, dmu(id) gen(eff) replacevars recoverylog\nassert _rc==459\nreplace x1=1 in 1\ndea x1 x2=y1 y2, dmu(id) stage(1) gen(eff) replacevars recoverylog\nassert r(recovery_stats)[1,1]==8\nmata: assert(_dea_rec_on==0)\n'
 run('lifecycle',body)
 summary=dict(invalid_commands=len(invalid),valid_combinations=len(cases),valid_executions=3*len(cases),numeric_files_identical=numeric,independent_zero_domain_LPs=oracle,max_radial_gap=maxgap,batches=len(runs),all_passed=True)
 (out/'invalid_cases.json').write_text(json.dumps(invalid,indent=2)+'\n');(out/'valid_cases.json').write_text(json.dumps(cases,indent=2)+'\n');(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)
if __name__=='__main__':main()
