#!/usr/bin/env python3
"""Portable native prepared-array/warning checks using private observation copies."""
from pathlib import Path
import argparse,csv,hashlib,itertools,json,math,shutil,subprocess,time
import numpy as np
from test_core_retries import WRITER,readbin
REPO=Path(__file__).resolve().parents[1]
OBSERVER='''
real scalar _dea_nonte_optimize(class LinearProgram scalar S,real scalar tol,|real scalar cap) {
 if(args()<3) cap=1e-8
 assert(tol==st_numscalar("PREFLIGHT_REQUESTED_TOL"))
 st_numscalar("PREFLIGHT_CALLS",st_numscalar("PREFLIGHT_CALLS")+1)
 return(_c66_observed_optimize(S,tol,cap))
}
'''
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';shutil.copytree(REPO/'ado',pkg)
 source={str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (REPO/'ado').iterdir() if p.is_file()}
 f=pkg/'_dea_nonte_lp.ado';s=f.read_text();assert s.count('real scalar _dea_nonte_optimize(')==1
 s=s.replace('real scalar _dea_nonte_optimize(','real scalar _c66_observed_optimize(',1);f.write_text(s.rsplit('\nend',1)[0]+OBSERVER+'\nend\n')
 code='version 16\nclear all\nset more off\nset processors 1\nset type double\nset linesize 200\n'
 for k in ['PLUS','PERSONAL','SITE','OLDPLACE']:(out/k).mkdir();code+=f'sysdir set {k} "{out/k}"\n'
 code+=f'adopath ++ "{pkg}"\n'+WRITER;cases=[]
 for fixture in json.loads((REPO/'tests/fixtures/crosseff_preflight.json').read_text())['fixtures']:
  file=out/(fixture['name']+'.csv')
  with file.open('w') as fh:
   w=csv.writer(fh);w.writerow(['id','x','y'])
   for i,e in enumerate(fixture['exponents']):
    z=math.ldexp(1.,e);assert math.isfinite(z) and z>0;w.writerow([i+1,format(z,'.17g'),format(z,'.17g')])
  for rts,ori,method,requested_tol in itertools.product(['crs','vrs'],['input','output'],['arbitrary','benevolent','aggressive'],['1e-4','1e-7','1e-12']):
   name='_'.join([fixture['name'],rts,ori,method,requested_tol]);d=out/name;d.mkdir()
   p=fixture['primary_blocked'];sec=fixture['secondary_blocked_when_requested']*(method!='arbitrary')
   code+=f'import delimited using "{file}",clear varnames(1) asdouble\nassert x>0 & y>0\nscalar PREFLIGHT_REQUESTED_TOL={requested_tol}\nscalar PREFLIGHT_CALLS=0\n'
   cmd=f'dea_crosseff x=y,dmu(id) rts({rts}) ort({ori}) method({method}) tol({requested_tol}) gen(res)'
   code+=cmd+' matrix detail\n'+f'assert r(tolerance)=={requested_tol}\nassert r(diagnostic_version)==6201 & r(numerical_diagnostic_version)==6601\nassert r(n_native_nonfinite_primary)=={p}\n'
   code+='assert r(has_numerical_warning)==(r(n_failed_primary)>0 | r(n_failed_secondary)>0 | r(n_missing_peer_pairs)>0 | r(n_failed_ce)>0 | r(n_normalization_lost_positive)>0 | r(n_objective_warning)>0)\n'
   if p:code+='assert r(has_numerical_warning)==1\n'
   code+='mata: assert(findexternal("_rsm_solve_lp()")==NULL)\nmatrix ND=r(native_numerical_diagnostics)\nmatrix ST=r(stage_status)\nmatrix WI=r(weights_input)\nmatrix WO=r(weights_output)\n'
   code+='mata: assert(st_numscalar("PREFLIGHT_CALLS")==sum(st_matrix("ND")[.,3]:<.)+sum(st_matrix("ND")[.,5]:<.))\n'
   if p:code+='assert missing(res_self) in 1\nmata: assert(st_matrix("ND")[1,2]==2 & st_matrix("ND")[1,3]>=.);assert(st_matrix("ST")[1,1]==0);assert(st_matrix("WI")[1,1]>=. & st_matrix("WO")[1,1]>=.);assert(all(st_matrix("r(ce_matrix)")[1,.]:>=.))\n'
   code+='assert r(n_native_nonfinite_secondary)==0\nassert r(objective_diagnostic_version)==7001\nmatrix OD=r(secondary_objective_diagnostics)\n'
   if sec:code+='local expect_secondary=(ND[1,2]==1)\nif `expect_secondary\' {\nassert abs(res_self-1)<1e-7 in 1\nmata: assert(st_matrix("OD")[1,2]==1 & st_matrix("OD")[1,3]==1)\n}\n'
   code+=f'mata: savebin("{d}/native.bin",st_matrix("ND"));savebin("{d}/stage.bin",st_matrix("ST"));savebin("{d}/data.bin",st_data(.,tokens("id x y")));savebin("{d}/objective.bin",st_matrix("OD"))\n'
   code+='matrix CC=(r(n_failed_primary),r(n_failed_secondary),r(n_missing_peer_pairs),r(n_failed_ce),r(n_normalization_lost_positive),r(n_native_nonfinite_primary),r(n_native_nonfinite_secondary),r(has_numerical_warning))\nmata: M66=st_data(.,tokens("res res_self res_maverick"))\nscalar CALLS66=PREFLIGHT_CALLS\nscalar PREFLIGHT_CALLS=0\n'
   code+=cmd+' replacevars\n'
   code+='matrix NC=(r(n_failed_primary),r(n_failed_secondary),r(n_missing_peer_pairs),r(n_failed_ce),r(n_normalization_lost_positive),r(n_native_nonfinite_primary),r(n_native_nonfinite_secondary),r(has_numerical_warning))\nmata: asserteq(st_matrix("CC"),st_matrix("NC"));asserteq(M66,st_data(.,tokens("res res_self res_maverick")))\nassert CALLS66==PREFLIGHT_CALLS\ncapture matrix list r(native_numerical_diagnostics)\nassert _rc==111\n'
   if sec:
    code+=f'if `expect_secondary\' & ST[1,1]==3 {{\ndea_crosseff x=y,gen(base) rts({rts}) ort({ori}) method(arbitrary) tol({requested_tol}) matrix\nassert base_self==res_self\nmata: assert(st_matrix("WI")[1,1]==st_matrix("r(weights_input)")[1,1]);assert(st_matrix("WO")[1,1]==st_matrix("r(weights_output)")[1,1])\n}}\n'
   cases.append(dict(test=name,primary_blocked=p,secondary_blocked=sec,requested_tolerance=float(requested_tol)))
 code+=f'import delimited using "{out/"ray_p512.csv"}",clear varnames(1) asdouble\nscalar PREFLIGHT_REQUESTED_TOL=1e-7\n'
 code+='di "PREFLIGHT_QUIET_BEGIN"\nquietly dea_crosseff x=y,gen(quietres) tol(1e-7)\nlocal quiet_flag=r(has_numerical_warning)\ndi "PREFLIGHT_QUIET_END"\nassert `quiet_flag\'==1\n'
 code+='di "CROSSEFF_PREFLIGHT_PASS"\nexit,clear\n';(out/'test.do').write_text(code);t=time.monotonic()
 q=subprocess.run([a.stata,'-b','do','test.do'],cwd=out,capture_output=True,timeout=1800)
 (out/'process.txt').write_bytes(q.stdout+q.stderr);log=(out/'test.log').read_text(errors='replace')
 assert q.returncode==0 and 'CROSSEFF_PREFLIGHT_PASS' in log.splitlines(),'\n'.join(log.splitlines()[-35:])
 assert 'Warning: ' in log and 'optimizer not called' in log
 quiet=log.split('\nPREFLIGHT_QUIET_BEGIN\n',1)[1].split('\nPREFLIGHT_QUIET_END\n',1)[0]
 assert 'Warning: 1 primary LP(s) skipped because prepared coefficients exceeded the finite numeric range.' in quiet
 assert 'primary weight solve(s) unavailable.' in quiet and 'averages use available peers.' in quiet
 for case in cases:
  d=out/case['test'];D=readbin(d/'native.bin');S=readbin(d/'stage.bin');data=readbin(d/'data.bin')
  assert np.all(np.isfinite(data)) and np.all(data[:,1:]>0)
  np.testing.assert_array_equal(D[:,0],data[:,0])
  reached=case['secondary_blocked']*int(D[0,1]==1)
  assert np.sum(D[:,1]==2)==case['primary_blocked'] and np.sum(D[:,3]==2)==0
  O=readbin(d/'objective.bin');assert not reached or (O[0,1]==1 and O[0,2]==1)
  case['old_secondary_overflow_reached']=reached
  case['secondary_blocked']=0
  assert np.all(S[D[:,1]==2,0]==0) and np.all(S[D[:,3]==2,0]==3)
 assert sum(z['old_secondary_overflow_reached'] for z in cases)>0,'No old secondary-overflow case reached;coverage is incomplete'
 for n,h in source.items():assert hashlib.sha256((REPO/n).read_bytes()).hexdigest()==h
 (out/'manifest.json').write_text(json.dumps(dict(source_sha256=source),indent=2)+'\n');(out/'results.json').write_text(json.dumps(cases,indent=2)+'\n')
 summary=dict(passed=True,cases=len(cases),matrix_and_nomatrix_commands=2*len(cases),conditional_arbitrary_fallback_comparisons=True,old_secondary_overflow_cases_processed=sum(z['old_secondary_overflow_reached'] for z in cases),quiet_warning_commands=1,full_quiet_warning_text_asserted=True,actual_tolerance_delivery_asserted=True,preflight_before_optimizer_asserted=True,cold_native_no_simplex=True,source_files_preserved=True,seconds=time.monotonic()-t)
 (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
