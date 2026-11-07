#!/usr/bin/env python3
"""Independent SW2007 quantile formulas and initial boundary-sample checks.

--reference-ado extracts the original cipsimarwilson program without editing
its installation or requiring its plugin. Boundary tests intercept after the
real initial regression, before bootstrap; they do not test coverage.
"""
from _support import package_files
import argparse,csv,hashlib,json,math,subprocess
from pathlib import Path
REPO=Path(__file__).resolve().parents[1]

def endpoints(values,level,public=False):
 v=sorted(values);n=len(v);a=(100-level)/200
 if not public:return v[max(1,math.ceil(n*a))-1],v[min(n,math.ceil(n*(1-a)))-1]
 lo=1+math.floor(n*a);hi=math.ceil(n-n*a)
 w=min(max((1-level/100)*n-((n-hi)+(lo-1)),0),2)
 if lo>1 and hi<n:
  return tuple((.5+.25*w)*x+(.5-.25*w)*y for x,y in [(v[lo-1],v[lo-2]),(v[hi-1],v[hi])])
 return v[lo-1],v[hi-1]

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',required=True,type=Path);ap.add_argument('--reference-ado',type=Path);a=ap.parse_args()
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir();hashes={}
 for name in package_files(REPO):
  raw=(REPO/name).read_bytes();(pkg/Path(name).name).write_bytes(raw);hashes[name]=hashlib.sha256(raw).hexdigest()
 manifest=dict(source_sha256=hashes);reference=''
 if a.reference_ado:
  raw=a.reference_ado.read_bytes();text=raw.decode();start=text.index('program cipsimarwilson, nclass');stop=text.index('\nend',start)+len('\nend');reference=text[start:stop]+'\n'
  (out/'public_original.ado').write_bytes(raw);(out/'public_quantile_program.ado').write_text(reference)
  manifest['reference']=dict(path=str(a.reference_ado.resolve()),sha256=hashlib.sha256(raw).hexdigest())
 (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 head='version 16\nclear all\nset more off\nset processors 1\nset type double\nset linesize 255\n'
 for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:
  (out/key).mkdir();head+=f'sysdir set {key} "{out/key}"\n'
 head+=f'adopath ++ "{pkg}"\nquietly run "{pkg}/_dea_simplex.ado"\nquietly run "{pkg}/_dea_sw2007_seq.ado"\n'+reference
 body='''tempname handle
postfile `handle' int n double level byte method col double lower upper pivot_lo pivot_hi using "quantiles.dta", double replace
'''
 expected={};cases=0
 for n in [2,3,9,10,39,40,41,99,100,199,399]:
  values=[[math.sin(i)*3+i/n,round((i%7)/3,2),math.exp(i/n)] for i in range(1,n+1)]
  # Read through a CSV to avoid Stata's command-line length limit at large B.
  with (out/f'draws_{n}.csv').open('w') as f:
   w=csv.writer(f);w.writerow(['a','b','c']);w.writerows(values)
  body+=f'import delimited using "draws_{n}.csv",clear varnames(1) asdouble\nmkmat a b c,matrix(D)\nmata: B=st_matrix("D")\nmatrix point=(.2,1.5,.3)\n'
  for level in [.1,50,80,90,95,99,99.9]:
   body+=f'mata: _sw07_seq_summary("D","point",{level},"V","M","C","P")\n'
   for method in range(1+bool(reference)):
    if method:
     body+=f"cipsimarwilson {level} B publicC\nmatrix C=publicC'\nmatrix P=(2*point-C[2,1...])\\(2*point-C[1,1...])\n"
    for col,point in enumerate([.2,1.5,.3],1):
     lo,hi=endpoints([v[col-1] for v in values],level,bool(method));expected[(n,level,method,col)]=(lo,hi,2*point-hi,2*point-lo)
     body+=f"post `handle' ({n}) ({level}) ({method}) ({col}) (C[1,{col}]) (C[2,{col}]) (P[1,{col}]) (P[2,{col}])\n"
   cases+=1
 body+='''postclose `handle'
use quantiles.dta,clear
format _all %24.17g
export delimited using "quantiles.csv",replace datafmt
// Intercept after the actual boundary handling and initial MLE.
capture program drop _dea_sw2007_seq
program define _dea_sw2007_seq, eclass
 syntax, SCORES(varname) TOUSE(varname) INITIAL(name) TAU(real) GENerate(name) *
 quietly gen double `generate'_eff=`scores' if `touse'
 quietly count if `touse' & `tau'/`scores'>`tau'
 local n_initial=r(N)
 ereturn clear
 ereturn matrix b_initial=`initial'
 ereturn scalar N_initial=`n_initial'
end
clear
set seed 72931
set obs 80
gen x=1+runiform()*4
gen z=2*runiform()-1
gen mu=1.5+.2*z
gen delta=mu-.3*invnormal((1-runiform())*normal((mu-1)/.3))
gen theta=1/delta
gen y=x/delta
'''
 bc=0
 for tol in [1e-7,1e-5]:
  for tau in [1,2]:
   for ori in ['input','output']:
    body+=f'''replace theta=1/delta
replace theta=1 in 1
replace theta=1-{tol}/2 in 2
replace theta=1-2*{tol} in 3
gen expected_theta=theta
replace expected_theta=1 if abs(expected_theta-1)<={tol}
gen expected_delta={tau}/expected_theta
quietly truncreg expected_delta z if expected_delta>{tau},ll({tau}) iterate(200)
matrix expected_b=e(b)
local expected_n=e(N)
quietly dea_sw2007 x=y, envars(z) scores(theta) ort({ori}) tau({tau}) tolerance({tol}) reps(2) bcreps(2) gen(test)
assert e(N_initial)==`expected_n'
assert test_eff==expected_theta
assert test_eff==1 in 1/2
assert test_eff<1 in 3
mata: assert(max(abs(st_matrix("e(b_initial)")-st_matrix("expected_b")))<1e-12)
drop test_eff expected_theta expected_delta
'''
    bc+=1
 body+='di "TEST_PASS"\n';(out/'test.do').write_text(head+body)
 proc=subprocess.run([a.stata,'-b','do','test.do'],cwd=out,capture_output=True,timeout=240)
 (out/'process.json').write_text(json.dumps(dict(returncode=proc.returncode,stdout=proc.stdout.decode(errors='replace'),stderr=proc.stderr.decode(errors='replace')),indent=2)+'\n')
 log=(out/'test.log').read_text(errors='replace') if (out/'test.log').exists() else ''
 assert proc.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-4000:]
 with (out/'quantiles.csv').open() as f:rows=list(csv.DictReader(f))
 assert len(rows)==len(expected)
 maxgap=0;diff=0
 for row in rows:
  key=tuple(float(row[x]) for x in ['n','level','method','col']);got=[float(row[x]) for x in ['lower','upper','pivot_lo','pivot_hi']];want=expected[key]
  gap=max(abs(x-y) for x,y in zip(got,want));maxgap=max(maxgap,gap);assert gap<1e-12,(key,got,want)
  if key[2] and want!=expected[(key[0],key[1],0,key[3])]:diff+=1
 result=dict(passed=True,quantile_grid_cases=cases,boundary_cases=bc,endpoint_checks=4*len(rows),maximum_gap=maxgap,reference_rule_differences=diff,reference_checked=bool(reference))
 (out/'results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
