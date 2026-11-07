#!/usr/bin/env python3
"""Opt-in optimal-face diagnostics: immutable independent rational references.
No reference solver or research-directory dependency at run time.
"""
from pathlib import Path
import argparse,csv,hashlib,json,subprocess,time
import numpy as np
from fractions import Fraction as F
from test_core_retries import WRITER,readbin
REPO=Path(__file__).resolve().parents[1]
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def point(p,w):
 return all(dot(list(map(F,a)),w)==F(b) for a,b in zip(p['E'],p['b'])) and all(dot(list(map(F,a)),w)<=F(b) for a,b in zip(p['U'],p['ub'])) and all(v is None or w[i]>=F(v) for i,v in enumerate(p['lo']))
def saved(prefix,diag=False):
 code=''
 for name,key in [('ce','ce_matrix'),('v','weights_input'),('u','weights_output'),('w','intercept'),('status','stage_status'),('peers','peer_count'),('pc','primary_cert_diagnostics'),('sc','secondary_cert_diagnostics'),('pr','primary_recovery_diagnostics'),('sr','secondary_recovery_diagnostics'),('objective','secondary_objective_diagnostics'),('native','native_numerical_diagnostics')]+([('faces','optima_faces'),('coords','optima_coordinates'),('denoms','optima_denominators')] if diag else []):
  code+=f'mata: savebin("{prefix}_{name}.bin",st_matrix("r({key})"))\n'
 code+='matrix nums=(r(n),r(mean_ce),r(mean_self),r(n_failed_primary),r(n_failed_secondary),r(n_failed_ce),r(n_missing_peer_pairs),r(n_missing_diagonal),r(tolerance),r(has_numerical_warning),r(n_objective_warning))\n'
 code+=f'mata: savebin("{prefix}_nums.bin",st_matrix("nums"));savebin("{prefix}_scores.bin",st_data(.,tokens("res res_self res_maverick")))\n'
 if diag:
  code+='matrix od=(r(n_optima_requested),r(n_optima_certified),r(n_optima_unique),r(n_optima_multiple),r(n_optima_unknown),r(n_optima_unsupported),r(n_optima_unbounded),r(n_optima_linked),r(n_optima_offface),r(n_optima_no_weights),r(n_optima_den_positive),r(n_optima_den_zero),r(n_optima_den_possible_zero),r(n_optima_den_unknown),r(n_optima_cutoff_variable),r(n_optima_lp_calls),r(n_optima_integer_ops),r(n_optima_link_unknown))\n'
  code+=f'mata: savebin("{prefix}_counts.bin",st_matrix("od"))\n'
 return code

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',required=True,type=Path);ap.add_argument('--baseline-ado',type=Path);args=ap.parse_args()
 out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=REPO/'ado';baseline=args.baseline_ado.resolve() if args.baseline_ado else pkg
 specs=json.loads((REPO/'tests/fixtures/crosseff_optima.json').read_text())['cases'];results=[]
 manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in pkg.iterdir() if p.is_file()}
 for z in specs:
  name=z['name'];m=z['m'];s=z['s'];xs=' '.join(f'x{j+1}' for j in range(m));ys=' '.join(f'y{j+1}' for j in range(s))
  for engine in ['ipm','simplex']:
   d=out/(name+'_'+engine);d.mkdir();f=d/'data.csv'
   with f.open('w') as h:
    w=csv.writer(h);w.writerow(['id']+xs.split()+ys.split());w.writerows([[format(v,'.17g') for v in rr] for rr in z['data']])
   head='version 16\nclear all\nset more off\nset type double\nset processors 1\nset linesize 200\n'
   for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/key).mkdir();head+=f'sysdir set {key} "{d/key}"\n'
   head+=WRITER+f'adopath ++ "{baseline}"\n'
   data=f'import delimited using "{f}",clear varnames(1) asdouble\n'
   cmd=f'dea_crosseff {xs}={ys},dmu(id) rts({z["rts"]}) ort({z["orientation"]}) method({z["method"]}) tol(1e-7) gen(res) engine({engine}) matrix'
   body=data+cmd+'\n'+saved('before')
   if baseline!=pkg:body+='clear all\n'+WRITER+f'adopath ++ "{pkg}"\n'
   body+=data+cmd+' optimadiag\n'+saved('after',True)+'di "OPTIMA73_PASS"\nexit,clear\n'
   (d/'test.do').write_text(head+body);start=time.monotonic();p=subprocess.run([args.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=1800)
   log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else '';assert p.returncode==0 and 'OPTIMA73_PASS' in log.splitlines(),'\n'.join(log.splitlines()[-40:])
   for f in d.glob('before_*.bin'):assert f.read_bytes()==(d/('after_'+f.name[7:])).read_bytes(),(name,engine,f.name)
   fd=readbin(d/'after_faces.bin');cd=readbin(d/'after_coords.bin');dd=readbin(d/'after_denoms.bin');counts=readbin(d/'after_counts.bin').ravel()
   ww=np.c_[readbin(d/'after_v.bin'),readbin(d/'after_u.bin')]
   if z['rts']=='vrs':ww=np.c_[ww,readbin(d/'after_w.bin')]
   assert len(fd)==2*len(z['data']) and counts[4]==0 and counts[13]==0 and counts[17]==0
   for ref in z['expected']:
    i=ref['evaluator'];w=list(map(lambda v:F(float(v)),ww[i-1]))
    for face in [1,2]:
     key='F'+str(face);q=ref[key];fr=fd[(fd[:,0]==i)&(fd[:,1]==face)][0]
     if q is None:assert fr[2]==2 and fr[3]==0;continue
     assert fr[2]==1 and fr[3]==(1 if q['weights']['state']=='unique' else 2)
     assert fr[6]==(1 if point(ref['faces'][key],w) else 2),(name,engine,i,face,'stored-weight link')
     for k,c in enumerate(q['weights']['coordinates']):
      row=cd[(cd[:,0]==i)&(cd[:,1]==face)&(cd[:,2]==k+1)][0]
      assert row[5]==(2 if c['lower_unbounded'] or c['upper_unbounded'] or c['lower']!=c['upper'] else 1)
      for e,key2 in enumerate(['lower','upper']):
       assert row[3+e]==(2 if c[key2+'_unbounded'] else 1)
       if not c[key2+'_unbounded']:assert row[6+e]==float(F(c[key2]))
     for j,c in enumerate(q['peers']):
      row=dd[(dd[:,0]==i)&(dd[:,1]==face)&(dd[:,2]==j+1)][0];den=c['denominator']
      assert row[5]=={'strictly_positive':1,'identically_zero':2,'zero_possible':3}[c['denominator_state']]
      cutoff=F(float(1e-7)) if z['rts']=='vrs' else F(0)
      expected=1 if F(den['lower'])>cutoff else 2 if not den['upper_unbounded'] and F(den['upper'])<=cutoff else 3
      assert row[6]==expected
      for e,key2 in enumerate(['lower','upper']):
       assert row[3+e]==(2 if den[key2+'_unbounded'] else 1)
       if not den[key2+'_unbounded']:assert row[7+e]==float(F(den[key2]))
   results.append(dict(test=name,engine=engine,passed=True,evaluators=len(z['data']),finite_faces=int(counts[1]),unique=int(counts[2]),multiple=int(counts[3]),unbounded=int(counts[6]),legacy_numeric_files=14,seconds=time.monotonic()-start))
   (out/'results.json').write_text(json.dumps(results,indent=2)+'\n');print(name,engine,'passed',flush=True)
 assert manifest=={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in pkg.iterdir() if p.is_file()}
 (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 (out/'summary.json').write_text(json.dumps(dict(passed=True,cases=len(results),evaluators=sum(z['evaluators'] for z in results),independent_reference=True,scope='12commands×2engines,same-data no-option/option pairs,exact states/bounds/cutoff and stored-weight links;not independent datasets'),indent=2)+'\n')
if __name__=='__main__':main()
