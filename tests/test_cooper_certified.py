#!/usr/bin/env python3
"""Portable Cooper engine, weight-domain and explicit selection-cap regression.
Requires licensed Stata16+ and NumPy. Full endpoint Fraction checks are also
run by test_cooper_congestion.py via cooper_certified/exact.py.
"""
from pathlib import Path
import argparse,json,subprocess,time,shutil
import numpy as np
from test_core_retries import WRITER,readbin

REPO=Path(__file__).resolve().parents[1]
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',required=True,type=Path);a=ap.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';shutil.copytree(REPO/'ado',pkg);runs=[]
 def run(name,body):
  d=out/name;d.mkdir();s='version 16\nclear all\nset more off\nset type double\nset processors 1\nset linesize 180\n'
  for k in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/k).mkdir();s+=f'sysdir set {k} "{d/k}"\n'
  s+=f'adopath ++ "{pkg}"\n'+WRITER+body+'\ndi "CERTIFIED_PASS"\nexit,clear\n';(d/'test.do').write_text(s);t=time.monotonic();p=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=1200);log=(d/'test.log').read_text();assert p.returncode==0 and 'CERTIFIED_PASS' in log.splitlines(),log[-3500:];runs.append(dict(name=name,passed=True,seconds=time.monotonic()-t));(out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n');return d
 for engine in ['simplex','ipm']:
  for weight in ['raw','normalized']:
   s='input double id x y\n1 1e300 1e-300\n2 1e300 2e-300\nend\ngen double old=765\n'
   if weight=='raw':s+=f'capture dea x=y,model(congestion) congm(cooper) engine({engine}) congw(raw) dmu(id) gen(old) replacevars\nassert _rc==459 & old[1]==765\ncapture confirm variable old_phi\nassert _rc==111\n'
   else:s+=f'dea x=y,model(congestion) congm(cooper) engine({engine}) congw(normalized) dmu(id) gen(new)\nassert r(n_failed)==0\nassert "`r(actual_engine)\'"=="{engine}"\n'
   if engine=='simplex' and weight=='normalized':s+='assert "`r(integration_status)\'"=="integrated"\n'
   run('units_'+engine+'_'+weight,s)
  run('reciprocal_'+engine,'input double id x y\n1 1e-320 1\n2 1 2\nend\ngen double old=765\n'+f'capture dea x=y,model(congestion) congm(cooper) engine({engine}) congw(normalized) dmu(id) gen(old) replacevars\nassert _rc==459 & old[1]==765\n')
 # IPM does not inherit the bounded certificate's64variable restriction.
 run('ipm_large','set obs 31\ngen double id=_n\ngen double x=1\ngen double y=1\ndea x=y,model(congestion) congm(cooper) engine(ipm) dmu(id) gen(e)\nassert r(N)==31 & r(n_failed)==0\nassert e_status==0\nassert "`r(actual_engine)\'"=="ipm"\nmata: assert(findexternal("_rsm_solve_lp()")==NULL)\n')
 fixtures=json.loads((REPO/'tests/fixtures/cooper_certified.json').read_text())['fixtures']
 for f in fixtures:
  X=np.array(f['X']);Y=np.array(f['Y']);n,m=X.shape;s=Y.shape[1];xs=' '.join('x'+str(j+1) for j in range(m));ys=' '.join('y'+str(j+1) for j in range(s));data='input double(id '+xs+' '+ys+')\n'+'\n'.join(' '.join(repr(float(v)) for v in [i+1,*X[i],*Y[i]]) for i in range(n))+'\nend\n'
  def literal(v):return '.' if v is None else repr(float(v))
  mats='matrix A=('+'\\'.join(','.join(map(literal,row)) for row in f['anchor'])+')\nmatrix M=('+'\\'.join(','.join(map(literal,row)) for row in f['reference'])+')\n'
  for cap in [16,128]:
   k=f['index'];status=3 if cap==16 else 1
   d=run(f'cap_{f["case"]}_{cap}',data+mats+f'dea {xs}={ys},model(congestion) congm(cooper) engine(simplex) congw(normalized) conga(A) congp(stage2) congref(M) congs({cap}) dmu(id) gen(e) congmat\nassert r(n_failed)==0\nassert "`r(integration_status)\'"=="integrated"\nmatrix C=r(certificate_status)\nassert C[{k},3]=={status}\nassert C[{k},1]==1 & C[{k},8]=={cap}\nassert C[{k},4]=={int(cap==128)}\n')
 summary=dict(passed=True,batches=len(runs),weight_parity_engines=['simplex','ipm'],native_large_sample_supported=True,default16_fallback_cases=[93,109],explicit128_certified_cases=[93,109],runs=runs)
 (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k!='runs'}),flush=True)
if __name__=='__main__':main()
