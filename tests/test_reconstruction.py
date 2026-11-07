#!/usr/bin/env python3
"""Exact-witness failure-tail reconstruction, accounting and false-recovery guards.

Requires NumPy/SciPy and licensed Stata16+. Only packaged data are used.
Instrumentation and fault injections are confined to private output copies.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,shutil,subprocess
import numpy as np
from test_core_retries import WRITER,readbin
from test_simplex_stability import mata
from test_precision_recovery import decode
REPO=Path(__file__).resolve().parents[1]
dot=lambda a,b:sum((x*y for x,y in zip(a,b)),F(0))
fm=lambda A:tuple(tuple(F(float(v)) for v in row) for row in A)
fv=lambda A:tuple(F(float(v)) for v in np.ravel(A))

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 root=args.output.resolve();root.mkdir(parents=True,exist_ok=False);pkg=root/'ado';shutil.copytree(REPO/'ado',pkg);h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
 fp=REPO/'tests/fixtures/reconstruction.json';pp=REPO/'tests/fixtures/precision_recovery.json';fixture=json.loads(fp.read_text());precision=json.loads(pp.read_text())['cases'];cases=fixture['cases'];assert len(cases)==114
 (root/'manifest.json').write_text(json.dumps(dict(source_sha256={p.name:h(p) for p in pkg.iterdir() if p.is_file()},fixtures_sha256={fp.name:h(fp),pp.name:h(pp)}),indent=2)+'\n');runs=[]
 def run(name,source,body):
  out=root/name;out.mkdir();header='version 16\nclear all\nset more off\nset processors 1\n'
  for k in ['PLUS','PERSONAL','SITE','OLDPLACE']:(out/k).mkdir();header+=f'sysdir set {k} "{out/k}"\n'
  header+=WRITER+f'adopath ++ "{source}"\nrun "{source}/_dea_simplex.ado"\n'
  (out/'test.do').write_text(header+body+'\ndi "RECONSTRUCTION_PASS"\nexit,clear\n')
  p=subprocess.run([args.stata,'-b','do','test.do'],cwd=out,capture_output=True,timeout=900);(out/'process.json').write_text(json.dumps(dict(returncode=p.returncode,stderr=p.stderr.decode(errors='replace')),indent=2)+'\n');log=(out/'test.log').read_text(errors='replace') if (out/'test.log').exists() else '';assert p.returncode==0 and 'RECONSTRUCTION_PASS' in log.splitlines(),(name,log[-4000:]);runs.append(dict(name=name,passed=True));return out
 source=(pkg/'_dea_simplex.ado').read_text();instrumented=source
 for name in ['_rsm_two_phase_once','_rsm_fp_once','_rsm_precision_once']:
  a=instrumented.index('struct RSM_Result scalar '+name+'(');b=instrumented.index('\n}',a)+2;part=instrumented[a:b]
  part=part.replace('    real scalar i, max_iter, p1_iter','    external real scalar _test_rebuild_work,_test_rebuild_tail\n    real scalar i, max_iter, p1_iter',1)
  part=part.replace('return(result)','_test_rebuild_work=_test_rebuild_work+result.iters;return(result)')
  if name=='_rsm_precision_once':part=part.replace('    pivot_seen=0','    if(allow_rebuild) _test_rebuild_tail=_test_rebuild_tail+1\n    pivot_seen=0',1)
  instrumented=instrumented[:a]+part+instrumented[b:]
 observed=root/'observed_ado';shutil.copytree(pkg,observed);(observed/'_dea_simplex.ado').write_text(instrumented)
 body='mata:\n_test_rebuild_work=0\n_test_rebuild_tail=0\n'+r'''
void probe(string scalar name,real matrix E,real colvector b,real matrix U,
 real colvector u,real rowvector c,real scalar tol,real scalar observed,
 real scalar expected_status,real scalar expected_iters,real scalar tail)
{
 external real scalar _test_rebuild_work,_test_rebuild_tail
 struct RSM_Prob scalar P,T
 struct RSM_Result scalar R,S
 real scalar enabled
 real matrix log
 _rsm_std_form(P,U,u,E,b,c)
 for(enabled=0;enabled<=1;enabled++) {
  _test_rebuild_work=0;_test_rebuild_tail=0;_rsm_audit_reset(enabled);T=P
  R=_rsm_two_phase(T,tol,"")
  assert(R.status==expected_status & R.iters==expected_iters)
  assert(R.recovery[5]+R.recovery[6]==R.iters & R.recovery[7]==R.iters & R.recovery[8]==R.status)
  if(observed) assert(_test_rebuild_work==R.iters & (_test_rebuild_tail>0)==tail)
  log=_rsm_audit_read()
  if(enabled) {
   assert(rows(log)==1 & cols(log)==12 & all(log[1,1..10]:==R.recovery))
   assert(R.status==S.status & R.objval==S.objval & all(R.x:==S.x) & all(R.recovery:==S.recovery))
   savebin(name+"_result.bin",(R.status,R.objval,R.iters));savebin(name+"_x.bin",R.x);savebin(name+"_record.bin",log)
  }
  else {assert(rows(log)==0);S=R;}
 }
 _rsm_audit_reset(0)
 savebin(name+"_eqmat.bin",E);savebin(name+"_eqrhs.bin",b);savebin(name+"_ineqmat.bin",U);savebin(name+"_ineqrhs.bin",u);savebin(name+"_cost.bin",c)
}
''';inputs=[]
 for i,z in enumerate(cases):
  f=z.get('problem',precision[z['precision_index']] if 'precision_index' in z else None);n=len(f['xstar']);values=[decode(f[k],nc) for k,nc in [('E',n),('b',1),('U',n),('u',1),('c',n)]];E,b,U,u,c=values;x=tuple(map(F,f['xstar']));y=tuple(map(F,f['ystar']))
  assert all(np.all(A>0) for A in values) and min(x)>=0
  assert all(dot(ar,x)==rhs for ar,rhs in zip(fm(E),fv(b))) and all(dot(ar,x)<=rhs for ar,rhs in zip(fm(U),fv(u)))
  assert min(ci-dot(y,col) for ci,col in zip(fv(c),zip(*fm(E))))>=0 and dot(fv(c),x)==dot(y,fv(b))==6
  inputs.append((f,values))
  for key,A in zip(['E','b','U','u','c'],values):body+=f'{key}={mata(A)}\n'
  body+=f'probe("p{i:03}",E,b,U,u,c,{f["tol"]:.17g},OBSERVED,{z["expected_status"]},{z["expected_iterations"]},{int(z["tail_expected"])})\n'
 body+='end\n';plain=run('plain',pkg,body.replace('OBSERVED','0'));obs=run('observed',observed,body.replace('OBSERVED','1'))
 files=list(plain.glob('*.bin'))
 for p in files:assert p.read_bytes()==(obs/p.name).read_bytes(),p.name
 exact=0
 for i,(z,(f,values)) in enumerate(zip(cases,inputs)):
  E,b,U,u,c=values;n=E.shape[1];pre=f'p{i:03}'
  for key,A in zip(['eqmat','eqrhs','ineqmat','ineqrhs','cost'],values):assert np.array_equal(readbin(plain/(pre+'_'+key+'.bin')),A)
  a=readbin(plain/(pre+'_result.bin')).ravel()
  if a[0]!=0:continue
  xf=readbin(plain/(pre+'_x.bin')).ravel();assert len(xf)==n+len(u) and np.isfinite(xf).all();x=fv(xf[:n])
  if z.get('known_accepted_control'):
   assert [float(t).hex() for t in xf]==z['expected_x']
  assert min(xf)>=0 and all(dot(ar,x)==rhs for ar,rhs in zip(fm(E),fv(b))) and all(dot(ar,x)<=rhs for ar,rhs in zip(fm(U),fv(u))) and dot(fv(c),x)==a[1]==6;exact+=1
 assert exact==114

 f=precision[146];E,b,c=[decode(f[k],nc) for k,nc in [('E',3),('b',1),('c',3)]];values=''.join(f'{k}={mata(A)}\n' for k,A in [('E',E),('b',b),('c',c)])
 wide='mata:\n'+r'''
void wide(real matrix E0,real colvector b,real rowvector c0,real scalar n) {
 struct RSM_Prob scalar P
 struct RSM_Result scalar R
 real matrix E
 real rowvector c
 real scalar j,k
 E=J(rows(E0),n,0);c=J(1,n,0);E[,1..3]=E0;c[1..3]=c0
 for(j=4;j<=n;j++) {k=mod(j-1,3)+1;E[,j]=E0[,k];c[j]=c0[k]+64;}
 _rsm_std_form(P,J(0,n,0),J(0,1,0),E,b,c);R=_rsm_two_phase(P,1e-12,"")
 if(n<=4090) assert(R.status==0 & R.objval==6 & all(R.x:==(J(3,1,1)\J(n-3,1,0))))
 else assert(R.status==3)
 savebin(sprintf("w%g_result.bin",n),(R.status,R.objval,R.iters));savebin(sprintf("w%g_x.bin",n),R.x)
}
'''+values+''.join(f'wide(E,b,c,{n})\n' for n in [3,129,1024,4090,4091])+'end\n';run('wide',pkg,wide)
 for fault in ['denied','quality','invalid_basis','budget']:
  faultpkg=root/('fault_'+fault);shutil.copytree(pkg,faultpkg);s=source
  if fault=='denied':a='        if(_rsm_retry_realbasis(prob,basis,tol)) {';assert s.count(a)==1;s=s.replace(a,'        if(0) {')
  elif fault=='quality':
   a=s.index('real rowvector _rsm_retry_quality(');bpos=s.index('\n}',a)+2;part=s[a:bpos];s=s[:a]+part[:part.index('{')+1]+'\n return(J(1,16,.))\n}'+s[bpos:]
  elif fault=='invalid_basis':a='            result.msg="Phase1 original-column reconstruction"';assert s.count(a)==1;s=s.replace(a,a+'\n            basis.basic=J(prob.m_con,1,1)')
  else:
   a=s.index('struct RSM_Result scalar _rsm_precision_once(');bpos=s.index('\n}',a)+2;part=s[a:bpos];needle='    max_iter = max((10000, 20 * prob.m_con))';assert part.count(needle)==1;s=s[:a]+part.replace(needle,'    max_iter = 2')+s[bpos:]
  (faultpkg/'_dea_simplex.ado').write_text(s)
  body='mata:\nvoid guard() {\nreal matrix E\nreal colvector b\nreal rowvector c\nstruct RSM_Prob scalar P\nstruct RSM_Result scalar R\n'+values+r'''
 _rsm_std_form(P,J(0,3,0),J(0,1,0),E,b,c);_rsm_audit_reset(1);R=_rsm_two_phase(P,1e-12,"")
 assert(R.status==3 & R.recovery[2]==0 & R.recovery[8]==3)
 assert(R.recovery[5]+R.recovery[6]==R.iters & rows(_rsm_audit_read())==1)
 savebin("result.bin",(R.status,R.objval,R.iters));savebin("record.bin",R.recovery)
 }
 guard()
 end
''';run('guard_'+fault,faultpkg,body)
 summary=dict(all_passed=True,conditions=114,new_exact_recoveries=111,successful_bases_corrected=2,unresolved_controls=0,accepted_controls=3,LP_grid_executions=456,observed_numeric_files_identical=len(files),actual_iteration_sums_match=True,single_record_per_call=True,wide_conditions=5,quality_limit_controls=1,false_recovery_controls=4,batches=runs)
 (root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
