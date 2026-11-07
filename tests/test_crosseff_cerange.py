#!/usr/bin/env python3
"""Explicit CE ranges against immutable, independent binary-rational references."""
from pathlib import Path
from fractions import Fraction as F
import argparse,csv,hashlib,json,subprocess,time
import numpy as np
from test_core_retries import WRITER,readbin
from test_crosseff_optima import saved
REPO=Path(__file__).resolve().parents[1]
COUNTERS=['requested','verified','fixed','variable','undefined','unknown','unsupported','face_unavailable','limit_only','unbounded','den_zero_possible','cutoff_variable','display_loss','display_unknown','offface','no_weights','link_unknown','raw_link_unknown','lp_calls','integer_ops']
def recount(a):
 return [len(a),sum(np.isin(a[:,4],[1,2])),sum(a[:,4]==1),sum(a[:,4]==2),sum(a[:,4]==0),sum(a[:,4]==-1),sum(a[:,4]==-2),sum(a[:,4]==-3),sum(np.any(a[:,9:11]==2,axis=1)),sum(np.any(a[:,9:11]==3,axis=1)),sum(a[:,5]==3),sum(a[:,6]==3),sum(np.any(np.isin(a[:,11:13],[3,4]),axis=1)),sum(np.any(a[:,11:13]==5,axis=1)),sum(a[:,15]==2),sum(a[:,15]==3),sum(a[:,15]==-1),sum(np.any(a[:,13:15]==-1,axis=1))]
def review(a,refs,unavailable,faces,counts):
 assert a.shape==(len(refs)+len(unavailable),20)
 assert len({tuple(x[:3]) for x in a})==len(a)
 np.testing.assert_array_equal(counts[:18],recount(a))
 for q in refs:
  face=int(q['face'][1]);r=a[(a[:,0]==q['evaluator'])&(a[:,1]==face)&(a[:,2]==q['peer'])][0]
  assert r[3]==1
  assert r[4]=={'fixed_on_positive_domain':1,'variable_on_positive_domain':2,'undefined_everywhere':0}[q['state']],(q,r)
  assert r[5]=={'strictly_positive':1,'identically_zero':2,'zero_possible':3}[q['denominator_state']]
  assert r[6]=={'always_included':1,'always_excluded':2,'mask_variable':3}[q['cutoff_state']]
  fr=faces[(faces[:,0]==q['evaluator'])&(faces[:,1]==face)][0]
  assert r[15]==fr[6] and r[16]==fr[7] and r[17]==0
  if not q['endpoints']:
   assert all(np.isnan(r[7:9])) and all(r[9:11]==-1) and all(r[11:15]==0)
  else:
   for j,e in enumerate(q['endpoints']):
    st={'attained':1,'limit_only':2,'unbounded':3}[e['state']]
    assert r[9+j]==st
    if st==3:assert np.isnan(r[7+j]) and r[11+j]==0 and r[13+j]==0
    else:
     value=F(e['value']);assert r[7+j]==float(value),(q,r,'bound')
     assert r[11+j]==(1 if F(float(value))==value else 2)
     assert r[13+j]==(1 if st==1 else 0)
 for q in unavailable:
  r=a[(a[:,0]==q['evaluator'])&(a[:,1]==int(q['face'][1]))&(a[:,2]==q['peer'])][0]
  assert r[3]==2 and r[4]==-3 and r[17]==8 and all(np.isnan(r[7:9]))
 assert counts[5]==0 and counts[18]<=4096 and counts[19]<=1e9
 return dict(cells=len(a),finite_face_cells=len(refs),face_unavailable=len(unavailable),fixed=int(counts[2]),variable=int(counts[3]),undefined=int(counts[4]),unknown=int(counts[5]),lp_calls=int(counts[18]),integer_ops=int(counts[19]))
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',required=True,type=Path);ap.add_argument('--baseline-ado',type=Path);args=ap.parse_args()
 out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=REPO/'ado';baseline=args.baseline_ado.resolve() if args.baseline_ado else pkg
 specs=json.loads((REPO/'tests/fixtures/crosseff_cerange.json').read_text());runs=[]
 manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in pkg.iterdir() if p.is_file()}
 for z in specs['cases']:
  name=z['name'];xs=' '.join('x'+str(i+1) for i in range(z['m']));ys=' '.join('y'+str(i+1) for i in range(z['s']))
  for engine in ['ipm','simplex']:
   d=out/(name+'_'+engine);d.mkdir();f=d/'data.csv'
   with f.open('w') as h:
    w=csv.writer(h);w.writerow(['id']+xs.split()+ys.split());w.writerows([[format(v,'.17g') for v in r] for r in z['data']])
   head='version 16\nclear all\nset more off\nset type double\nset processors 1\nset linesize 200\n'
   for k in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/k).mkdir();head+=f'sysdir set {k} "{d/k}"\n'
   head+=WRITER+f'adopath ++ "{baseline}"\n'
   data=f'import delimited using "{f}",clear varnames(1) asdouble\n'
   cmd=f'dea_crosseff {xs}={ys},dmu(id) rts({z["rts"]}) ort({z["orientation"]}) method({z["method"]}) tol(1e-7) gen(res) engine({engine}) matrix'
   body=data+cmd+' optimadiag\n'+saved('before',True)
   if baseline!=pkg:body+='clear all\n'+WRITER+f'adopath ++ "{pkg}"\n'
   body+=data+cmd+' cerange\n'+saved('after',True)
   body+='assert r(cerange_version)==7701 & r(optima_diagnostic_version)==7301\nassert r(optima_ce_ranges_available)==1 & r(cerange_rank_ranges_available)==0 & r(cerange_mean_ranges_available)==0\n'
   body+='matrix CR=('+','.join('r(n_cerange_'+k+')' for k in COUNTERS)+')\nmata: savebin("ranges.bin",st_matrix("r(ce_ranges)"));savebin("counts.bin",st_matrix("CR"))\ndi "CERANGE77_PASS"\nexit,clear\n'
   (d/'test.do').write_text(head+body);t=time.monotonic();p=subprocess.run([args.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=2400)
   log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else '';assert p.returncode==0 and 'CERANGE77_PASS' in log.splitlines(),(name,engine,'\n'.join(log.splitlines()[-40:]))
   legacy=list(d.glob('before_*.bin'));assert len(legacy)==18
   for f in legacy:assert f.read_bytes()==(d/('after_'+f.name[7:])).read_bytes(),(name,engine,f.name)
   r=review(readbin(d/'ranges.bin'),z['expected'],z['unavailable'],readbin(d/'after_faces.bin'),readbin(d/'counts.bin').ravel())
   runs.append(dict(test=name,engine=engine,passed=True,seconds=time.monotonic()-t,legacy_numeric_pairs=len(legacy),**r))
   (out/'results.json').write_text(json.dumps(runs,indent=2)+'\n');print(name,engine,'passed',flush=True)
 assert manifest=={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in pkg.iterdir() if p.is_file()}
 summary=dict(passed=True,batches=len(runs),cells=sum(r['cells'] for r in runs),finite_face_cells=sum(r['finite_face_cells'] for r in runs),face_unavailable=sum(r['face_unavailable'] for r in runs),unknown=sum(r['unknown'] for r in runs),legacy_numeric_pairs=sum(r['legacy_numeric_pairs'] for r in runs),source_manifest=manifest,scope='12 overlapping datasets/configurations x 2 native engines; diagnostics use the shared independent exact path. No universal accuracy or native IPM repair claim.')
 (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
if __name__=='__main__':main()
