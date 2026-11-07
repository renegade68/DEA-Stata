#!/usr/bin/env python3
"""Same-pseudo-sample MLE recovery, unchanged accepted results, bounded failure."""
from _support import package_files
import argparse,difflib,hashlib,json,subprocess
from pathlib import Path
from test_sw2007 import DATA,REPO
from test_output_protection import unchanged

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',required=True,type=Path);a=ap.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
 names=package_files(REPO);hashes={n:hashlib.sha256((REPO/n).read_bytes()).hexdigest() for n in names};(out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes),indent=2)+'\n');results=[]
 for case in ['control','short_failure','large_limit','persistent_failure']:
  d=out/case;d.mkdir();pkg=d/'ado';pkg.mkdir()
  for n in names:(pkg/Path(n).name).write_bytes((REPO/n).read_bytes())
  f=pkg/'_dea_sw2007_seq.ado';original=f.read_text();modified=original
  if case!='control':
   needle="        local fit_ok = (_rc == 0)";assert original.count(needle)==2
   injection="        if `reg_try'==1 ereturn scalar converged=0\n" if case in ['short_failure','large_limit'] else '        ereturn scalar converged=0\n'
   modified=original.replace(needle,injection+needle,1 if case in ['short_failure','large_limit'] else 2)
   if case=='large_limit':modified=modified.replace('local maxiter = max(25,3*`iterations\')','local maxiter = 300')
   f.write_text(modified)
   (d/'instrumentation.patch').write_text(''.join(difflib.unified_diff(original.splitlines(True),modified.splitlines(True))))
  cmd='dea_sw2007 x=y, envars(z) rts(crs) reps(20) bcreps(30) seed(127) gen(a) replacevars'
  body=DATA+'gen a_eff=991\ngen a_eff_bc=992\n'
  if case=='persistent_failure':body+=unchanged(cmd,498)+'\nassert "`e(cmd)\'"==""\n'
  else:
   body+=cmd+'''
assert e(boot_success)==20
assert e(boot_fit_attempts)==20+e(boot_failed)+e(boot_retried)
'''
   if case=='control':
    body+='''assert e(boot_failed)==0 & e(boot_retried)==0 & e(boot_recovered)==0
estimates save expected,replace
keep a_eff a_eff_bc
save expected,replace
file open fh using rng.txt,write replace
file write fh "`c(rngstate)'"
file close fh
'''
   elif case=='large_limit':
    body+=f'''assert e(boot_failed)==1 & e(boot_retried)==0 & e(boot_recovered)==0
cf a_eff a_eff_bc using "{out}/control/expected.dta",all
matrix got=e(b_reps)
estimates use "{out}/control/expected.ster"
mata: assert(max(abs(st_matrix("got")[1..19,.]-st_matrix("e(b_reps)")[2..20,.]))==0)
'''
   else:
    body+=f'''assert e(boot_failed)==0 & e(boot_retried)==1 & e(boot_recovered)==1
cf a_eff a_eff_bc using "{out}/control/expected.dta",all
local current_rng "`c(rngstate)'"
file open fh using "{out}/control/rng.txt",read
file read fh expected_rng
file close fh
assert "`current_rng'"=="`expected_rng'"
local names b V b_initial V_initial V_modelbased b_reps b_boot sigma_boot b_mean ci_percentile ci_pivot b_bc beta_bc ci_lo ci_hi
foreach name of local names {{
 matrix test_`name'=e(`name')
}}
estimates use "{out}/control/expected.ster"
foreach name of local names {{
 mata: assert(max(abs(st_matrix("test_`name'")-st_matrix("e(`name')")))==0)
}}
'''
  head='version 16\nclear all\nset more off\nset processors 1\n'
  for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/key).mkdir();head+=f'sysdir set {key} "{d/key}"\n'
  (d/'test.do').write_text(head+f'adopath ++ "{pkg}"\n'+body+'\ndi "TEST_PASS"\n');p=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=180)
  (d/'process.json').write_text(json.dumps(dict(returncode=p.returncode,stderr=p.stderr.decode(errors='replace')),indent=2)+'\n');log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else '';row=dict(test=case,passed=p.returncode==0 and 'TEST_PASS' in log.splitlines())
  if not row['passed']:row['tail']=log[-2000:]
  results.append(row);(out/'results.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(row),flush=True)
 raise SystemExit(any(not r['passed'] for r in results))
if __name__=='__main__':main()
