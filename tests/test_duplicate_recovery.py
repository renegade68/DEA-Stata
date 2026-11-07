#!/usr/bin/env python3
"""Exact duplicate recovery: stored-double witnesses, priority, accounting, guards.

Requires NumPy/SciPy and licensed Stata16+. Only packaged fixtures are used.
Private instrumentation/fault copies are confined to the requested new output.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,shutil,subprocess,time
import numpy as np
from test_core_retries import WRITER,readbin
from test_simplex_stability import mata
from test_precision_recovery import decode
REPO=Path(__file__).resolve().parents[1]
dot=lambda a,b:sum((x*y for x,y in zip(a,b)),F(0))
fm=lambda A:tuple(tuple(F(float(v)) for v in row) for row in A)
fv=lambda A:tuple(F(float(v)) for v in np.ravel(A))


def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args()
 r=args.output.resolve();r.mkdir(parents=True,exist_ok=False);pkg=r/'ado';shutil.copytree(REPO/'ado',pkg)
 (r/'checks').mkdir();h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();fp=REPO/'tests/fixtures/duplicate_recovery.json';fixture=json.loads(fp.read_text())['cases'];assert len(fixture)==1215
 (r/'manifest.json').write_text(json.dumps(dict(source_sha256={p.name:h(p) for p in pkg.iterdir() if p.is_file()},fixtures_sha256={fp.name:h(fp)}),indent=2)+'\n')
 source=(pkg/'_dea_simplex.ado').read_text();runs=[]
 def inputs(index):
  f=fixture[index]['problem'];n=len(f['xstar']);return f,[decode(f[k],nc) for k,nc in [('E',n),('b',1),('U',n),('u',1),('c',n)]]
 def assign(values):return ''.join(k+'='+mata(a)+'\n' for k,a in zip(['E','b','U','u','c'],values))
 def run(name,source,body):
  out=r/'checks'/name;out.mkdir();head='version 16\nclear all\nset more off\nset processors 1\nset linesize 160\n'
  for k in ['PLUS','PERSONAL','SITE','OLDPLACE']:(out/k).mkdir();head+=f'sysdir set {k} "{out/k}"\n'
  head+=WRITER+f'adopath ++ "{source}"\n';(out/'test.do').write_text(head+body+'\ndi "DUPLICATE_PASS"\nexit,clear\n');start=time.monotonic()
  p=subprocess.run([args.stata,'-b','do','test.do'],cwd=out,capture_output=True,timeout=900)
  (out/'process.json').write_text(json.dumps(dict(returncode=p.returncode,seconds=time.monotonic()-start,stdout=p.stdout.decode(errors='replace'),stderr=p.stderr.decode(errors='replace')),indent=2)+'\n')
  log=(out/'test.log').read_text(errors='replace') if (out/'test.log').exists() else '';assert p.returncode==0 and 'DUPLICATE_PASS' in log.splitlines(),(name,log[-4500:])
  runs.append(dict(name=name,passed=True,seconds=time.monotonic()-start));return out
 # Observers sum actual inner exits; no numerical operations are replaced.
 s=source
 for name in ['_rsm_two_phase_once','_rsm_fp_once','_rsm_precision_once']:
  a=s.index('struct RSM_Result scalar '+name+'(');b=s.index('\n}',a)+2;part=s[a:b]
  part=part.replace('    real scalar i, max_iter, p1_iter','    external real scalar _test_dup_work\n    real scalar i, max_iter, p1_iter',1).replace('return(result)','_test_dup_work=_test_dup_work+result.iters;return(result)');s=s[:a]+part+s[b:]
 a=s.index('struct RSM_Result scalar _rsm_duplicate_retry(');b=s.index('\n}',a)+2;part=s[a:b]
 part=part.replace('    real scalar work','    external real scalar _test_dup_depth,_test_dup_max,_test_dup_calls\n    real scalar work',1).replace('    result.status=3;result.objval=.;result.x=J(original.n_orig,1,.)','    _test_dup_calls++\n    result.status=3;result.objval=.;result.x=J(original.n_orig,1,.)',1)
 part=part.replace('    result=_rsm_two_phase_core(reduced,tol,trace,0,check_success,0)','    _test_dup_depth++\n    _test_dup_max=max((_test_dup_max,_test_dup_depth))\n    result=_rsm_two_phase_core(reduced,tol,trace,0,check_success,0)\n    _test_dup_depth--')
 s=s[:a]+part+s[b:];observed=r/'observed_ado';shutil.copytree(pkg,observed);(observed/'_dea_simplex.ado').write_text(s)
 body='mata:\n_test_dup_work=0\n_test_dup_depth=0\n_test_dup_max=0\n_test_dup_calls=0\n'+r'''
void probe(string scalar name,real matrix E,real colvector b,real matrix U,real colvector u,real rowvector c,
 real scalar tol,real scalar expected,real scalar newly,real scalar observed) {
 external real scalar _test_dup_work,_test_dup_depth,_test_dup_max,_test_dup_calls
 struct RSM_Prob scalar P,T
 struct RSM_Result scalar R,S
 real scalar enabled,f
 real matrix log
 _rsm_std_form(P,U,u,E,b,c)
 for(enabled=0;enabled<=1;enabled++) {
  _test_dup_work=0;_test_dup_depth=0;_test_dup_max=0;_test_dup_calls=0
  _rsm_audit_reset(enabled);T=P;R=_rsm_two_phase(T,tol,"")
  assert(R.status==0 & R.iters==expected)
  assert(R.recovery[5]+R.recovery[6]==R.iters & R.recovery[7]==R.iters & R.recovery[8]==0)
  if(observed) assert(_test_dup_work==R.iters & _test_dup_depth==0 & _test_dup_max==newly & _test_dup_calls==newly)
  log=_rsm_audit_read()
  if(enabled) {
   assert(rows(log)==1 & cols(log)==12 & all(log[1,1..10]:==R.recovery))
   assert(R.status==S.status & R.objval==S.objval & all(R.x:==S.x) & R.msg==S.msg & all(R.recovery:==S.recovery))
   savebin(name+"_result.bin",(R.status,R.objval,R.iters));savebin(name+"_x.bin",R.x);savebin(name+"_record.bin",log)
   f=fopen(name+"_message.txt","w");fput(f,R.msg);fclose(f)
  }
  else {assert(rows(log)==0);S=R;}
 }
 _rsm_audit_reset(0)
 savebin(name+"_eqmat.bin",E);savebin(name+"_eqrhs.bin",b);savebin(name+"_ineqmat.bin",U);savebin(name+"_ineqrhs.bin",u);savebin(name+"_cost.bin",c)
}
''';data=[]
 for i,z in enumerate(fixture):
  f,values=inputs(i);E,b,U,u,c=values;x=tuple(map(F,f['xstar']));y=tuple(map(F,f['ystar']))
  assert all(np.all(A>0) for A in values) and min(x)>=0
  assert all(dot(ar,x)==rhs for ar,rhs in zip(fm(E),fv(b))) and all(dot(ar,x)<=rhs for ar,rhs in zip(fm(U),fv(u)))
  assert min(ci-dot(y,col) for ci,col in zip(fv(c),zip(*fm(E))))>=0 and dot(fv(c),x)==dot(y,fv(b))==6
  data.append(values);body+=assign(values)+f'probe("p{i:04}",E,b,U,u,c,{f["tol"]:.17g},{z["expected_iterations"]},{int(z["previous_status"]==3)},OBSERVED)\n'
 body+='end\n';plain=run('plain',pkg,f'run "{pkg}/_dea_simplex.ado"\n'+body.replace('OBSERVED','0'));obs=run('observed',observed,f'run "{observed}/_dea_simplex.ado"\n'+body.replace('OBSERVED','1'))
 files=[p for p in plain.iterdir() if p.suffix=='.bin' or p.name.endswith('_message.txt')]
 for p in files:assert p.read_bytes()==(obs/p.name).read_bytes(),p.name
 exact=preserved=corrected=0
 for i,(z,values) in enumerate(zip(fixture,data)):
  pre=f'p{i:04}';E,b,U,u,c=values;n=E.shape[1]
  for key,A in zip(['eqmat','eqrhs','ineqmat','ineqrhs','cost'],values):assert np.array_equal(readbin(plain/(pre+'_'+key+'.bin')),A)
  result=readbin(plain/(pre+'_result.bin')).ravel();xf=readbin(plain/(pre+'_x.bin')).ravel()
  assert result[0]==0 and result[2]==z['expected_iterations'] and [float(t).hex() for t in xf]==z['expected_x']
  assert (plain/(pre+'_message.txt')).read_text()==z['expected_message']
  assert [float(t).hex() for t in result]==z['expected_result']
  assert [float(t).hex() for t in readbin(plain/(pre+'_record.bin')).ravel()]==z['expected_record']
  if z['previous_status']==3:
   x=fv(xf[:n]);assert min(xf)>=0 and all(dot(ar,x)==rhs for ar,rhs in zip(fm(E),fv(b))) and all(dot(ar,x)<=rhs for ar,rhs in zip(fm(U),fv(u))) and dot(fv(c),x)==result[1]==6;exact+=1
  else:
   assert z['previous_status']==0 and z['previous_iterations']==result[2]
   if 'pre_quality_x' in z:
    x=fv(xf[:n]);assert min(xf)>=0 and all(dot(ar,x)==rhs for ar,rhs in zip(fm(E),fv(b))) and all(dot(ar,x)<=rhs for ar,rhs in zip(fm(U),fv(u))) and dot(fv(c),x)==result[1]==6;corrected+=1
   else:preserved+=1
 assert (exact,preserved,corrected)==(388,784,43)
 body=f'run "{pkg}/_dea_simplex.ado"\nmata:\n'+r'''
 void maps(real matrix E,real colvector b,real matrix U,real colvector u,real rowvector c,real scalar nr,real scalar nc,real scalar eligible) {
  struct RSM_Prob scalar P
  real colvector kr,kc,rep
  real scalar work,ok
  _rsm_std_form(P,U,u,E,b,c)
  ok=_rsm_duplicate_map(P,kr,kc,rep,work)
  assert(ok==eligible & rows(kr)==nr & rows(kc)==nc & work<=1000000)
 }
 void basicguards() {
  struct RSM_Prob scalar P,S
  struct RSM_Result scalar R
  real colvector kr,kc,rep
  real scalar work,ok
  assert(_rsm_duplicate_ratio((2,-4,0),(1,-2,0))==2)
  assert(_rsm_duplicate_ratio((1,-2,0),(2,-4,0))==.5)
  assert(_rsm_duplicate_ratio((1,2),(1,2))==1)
  assert(_rsm_duplicate_ratio(1,3.8725919148493183e-121)==2.5822498780869086e+120)
  assert(_rsm_duplicate_ratio(1.499696813895631e-241,3.8725919148493183e-121)==3.8725919148493183e-121)
  assert(missing(_rsm_duplicate_ratio((1+2^-52,2),(1,2))))
  assert(missing(_rsm_duplicate_ratio((0,0),(0,0))))
  assert(missing(_rsm_duplicate_ratio((-1,-2),(1,2))))
  assert(missing(_rsm_duplicate_ratio((3,6),(1,2))))
  assert(missing(_rsm_duplicate_ratio(2^401,1)))
  assert(missing(_rsm_duplicate_ratio(2^400,1))) // Approximate Mata power is not an exact power of two.
  assert(missing(_rsm_duplicate_ratio(2^-1022,2^-1022)))
  assert(missing(_rsm_duplicate_ratio((.,1),(1,1))))
  _rsm_std_form(P,J(0,2,0),J(0,1,0),(1,1\2,2),(1\2),(1,1))
  ok=_rsm_duplicate_map(P,kr,kc,rep,work,0);assert(!ok & work==0)
  P.m_con=.;assert(!_rsm_duplicate_map(P,kr,kc,rep,work));P.m_con=2
  P.A=J(2,4093,1);P.n_orig=4093;P.c=J(1,4093,1)
  ok=_rsm_duplicate_map(P,kr,kc,rep,work);assert(!ok & work==0)
  // Full original bound rows block unsafe column merging.
  maps((1,1\2,2),(1\2),I(2),(2\2),(1,1),3,4,1)
  // Exact zero rows/columns remain: this transform is intentionally narrower.
  maps((0,0\1,0\2,0),(0\1\2),J(0,2,0),J(0,1,0),(1,0),2,2,1)
  // Signed coefficients and negative RHS preserve equivalent standard rows.
  maps((1,-1\2,-2),(-1\-2),J(0,2,0),J(0,1,0),(1,2),1,2,1)
  // Direct helper never exports a new infeasible/unbounded classification.
  _rsm_std_form(P,J(0,2,0),J(0,1,0),(1,1\2,2),(1\3),(1,1))
  R=_rsm_duplicate_retry(P,1e-12,"");assert(R.status==3)
  _rsm_std_form(P,J(0,2,0),J(0,1,0),(1,-1\2,-2),(0\0),(-1,-1))
  R=_rsm_duplicate_retry(P,1e-12,"");assert(R.status==3)
  // A degenerate support without a certified basis must be rejected.
  _rsm_std_form(P,J(0,2,0),J(0,1,0),(1,0\0,1),(1\0),(1,1));S=P
  R.status=0;R.x=(1\0);R.objval=1;kr=kc=(1\2)
  assert(!_rsm_duplicate_accept(P,S,R,kr,kc,1e-12))
 }
 basicguards()
 '''
 v=inputs(81)[1];E,b,U,u,c=v;E=E.copy();E[3,2]=np.nextafter(E[3,2],np.inf);body+=assign([E,b,U,u,c])+'maps(E,b,U,u,c,4,7,0)\n'
 E,b,U,u,c=inputs(81)[1];b=b.copy();b[3,0]=np.nextafter(b[3,0],np.inf);body+=assign([E,b,U,u,c])+'maps(E,b,U,u,c,4,7,0)\n'
 for i,kind in [(162,'cost'),(891,'ineq')]:
  E,b,U,u,c=inputs(i)[1];n=E.shape[1]
  if kind=='cost':c[0,-1]=np.nextafter(c[0,-1],np.inf)
  else:U[0,-1]=np.nextafter(U[0,-1],np.inf)
  body+=assign([E,b,U,u,c])+f'maps(E,b,U,u,c,{len(E)+len(U)},{n+len(U)},0)\n'
 body+='end\n';run('map_guards',pkg,body)
 f,values=inputs(83);old=fixture[83]['previous_iterations'];guardbody='mata:\nvoid test() {\nstruct RSM_Prob scalar P\nstruct RSM_Result scalar R\nreal matrix E,U\nreal colvector b,u\nreal rowvector c\n'+assign(values)+r'''
  _rsm_std_form(P,U,u,E,b,c);_rsm_audit_reset(1);R=_rsm_two_phase(P,1e-12,"")
  assert(R.status==3 & R.recovery[2]==0 & R.recovery[8]==3)
  assert(R.recovery[5]+R.recovery[6]==R.iters & rows(_rsm_audit_read())==1)
  savebin("result.bin",(R.status,R.objval,R.iters));savebin("recovery.bin",R.recovery)
 }
 test()
 end
 '''
 faults=[]
 for fault in ['budget','lifted_primal','lifted_dual','quality','support','iteration_budget']:
  d=r/('fault_'+fault);shutil.copytree(pkg,d);s=source.replace('if(args()<6) allow_rows=1','if(args()<6) allow_rows=0') # Isolate legacy duplicate-path fault controls; T35 separately faults its new path.
  if fault=='budget':s=s.replace('_rsm_duplicate_map(original,kr,kc,representative,work))','_rsm_duplicate_map(original,kr,kc,representative,work,0))')
  elif fault=='lifted_primal':s=s.replace('    Y[kr,.]=YP','    x=J(original.n_orig,1,0)\n    Y[kr,.]=YP')
  elif fault=='lifted_dual':s=s.replace('    Y[kr,.]=YP','    Y[kr,.]=YP\n    Y[,1]=Y[,1]:+1')
  elif fault=='quality':s=s.replace('    if(!_dea_lpcheck_ok(check,tol) | hasmissing(quality)', '    quality=J(1,16,.)\n    if(!_dea_lpcheck_ok(check,tol) | hasmissing(quality)')
  elif fault=='support':s=s.replace('    basic=selectindex(candidate.x:>0)','    candidate.x=J(n,1,0)\n    basic=selectindex(candidate.x:>0)')
  else:
   for name in ['_rsm_two_phase_once','_rsm_fp_once','_rsm_precision_once']:
    a=s.index('struct RSM_Result scalar '+name+'(');z=s.index('\n}',a)+2;part=s[a:z];part=part.replace('    real scalar i, max_iter, p1_iter','    external real scalar _fault_reduced\n    real scalar i, max_iter, p1_iter',1);anchor='    max_iter = max((10000, 20 * prob.m_con))';assert part.count(anchor)==1;part=part.replace(anchor,anchor+'\n    if(_fault_reduced) max_iter=1');s=s[:a]+part+s[z:]
   a=s.index('struct RSM_Result scalar _rsm_duplicate_retry(');z=s.index('\n}',a)+2;part=s[a:z].replace('    real scalar work','    external real scalar _fault_reduced\n    real scalar work').replace('    result=_rsm_two_phase_core(reduced,tol,trace,0,check_success,0)','    _fault_reduced=1\n    result=_rsm_two_phase_core(reduced,tol,trace,0,check_success,0)\n    _fault_reduced=0');s=s[:a]+part+s[z:]
  (d/'_dea_simplex.ado').write_text(s);out=run('guard_'+fault,d,f'run "{d}/_dea_simplex.ado"\nmata: _fault_reduced=0\n'+guardbody);result=readbin(out/'result.bin').ravel();assert result[2]>=old
  if fault=='budget':assert result[2]==old
  faults.append(dict(name=fault,passed=True,iterations=int(result[2]),old_iterations=old))
 summary=dict(all_passed=True,conditions=1215,new_exact_recoveries=exact,previous_outputs_preserved=preserved,successful_bases_corrected=corrected,already_exact_preserved=preserved+exact,LP_grid_executions=4860,observed_numeric_message_files_identical=len(files),actual_iteration_sums_match=True,single_record_per_call=True,max_reduction_depth=1,actual_Mata_one_ULP_controls=4,exact_ratio_assertions=13,other_map_and_certificate_controls=9,false_recovery_controls=faults,batches=runs)
 (r/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
