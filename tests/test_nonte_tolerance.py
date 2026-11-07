#!/usr/bin/env python3
"""Core non-TE tolerance delivery, original-LP validation and independent HiGHS checks.

Requires Stata 16+, NumPy/SciPy. Isolated copies record submitted LPs without
changing solves; separate copies inject bad candidates. Mathematical formula
review is distinct from checking the submitted LP against another solver.
"""
from _support import package_files
from pathlib import Path
import argparse, hashlib, itertools, json, re, shutil, subprocess, time
import numpy as np
from scipy.optimize import linprog
from test_core_retries import WRITER, readbin
from test_lp_residuals import residual
from test_specialized_failures import DATA, lp as specialized_lp
REPO=Path(__file__).resolve().parents[1]
PROFILES={
 'cost':('cost','w(w1 w2)','input'), 'revenue':('revenue','p(p1 p2)','output'),
 'profit':('profit','w(w1 w2) p(p1 p2)','input'),
 'alloc_in':('allocative','w(w1 w2)','input'), 'alloc_out':('allocative','p(p1 p2)','output'),
 'ddf':('ddf','gx(gx1 gx2) gy(gy1 gy2)','input'), 'additive':('additive','','input'),
 'sbm':('sbm','','input'), 'sbmin':('sbm','ort(sbm_input)','input'), 'sbmout':('sbm','ort(sbm_output)','output'),
 'superin':('supereff','','input'), 'superout':('supereff','','output'),
 'supersbm':('supersbm','','input'), 'multin':('multiplier','saveweights','input'), 'multout':('multiplier','saveweights','output'),
 'supersbmin':('supersbm','','input'), 'supersbmout':('supersbm','','output'),
 'bad_yb':('undesirable','bad(b1 b2) direction(yb)','input'),
 'bad_x':('undesirable','bad(b1 b2) direction(x)','input'),
 'bad_xyb':('undesirable','bad(b1 b2) direction(xyb)','input'),
 'intervalin':('imprecise','xlo(x1l x2l) xhi(x1h x2h) ylo(y1l y2l) yhi(y1h y2h)','input'),
 'intervalout':('imprecise','xlo(x1l x2l) xhi(x1h x2h) ylo(y1l y2l) yhi(y1h y2h)','output'),
 'congestion':('congestion','','input'), 'sbmu':('sbmundesirable','bad(b1 b2)','input'),
}
PROFILES['sbm']=('sbm','','')
PROFILES['supersbm']=('supersbm','','')
PROFILES['sbmin']=('sbm','','input');PROFILES['sbmout']=('sbm','','output')
CAPTURE='''mata:
void r08save(real rowvector c,real matrix E,real matrix U,real matrix B,
 real rowvector x,real scalar status,real scalar obj,real scalar tol,real scalar actual,real scalar sign) {
 string scalar prefix
 real scalar i
 i=st_numscalar("r08calls")+1;st_numscalar("r08calls",i)
 prefix="lp"+strofreal(i,"%06.0f")
 savebin(prefix+"_c.bin",c);savebin(prefix+"_E.bin",E);savebin(prefix+"_U.bin",U)
 savebin(prefix+"_B.bin",B);savebin(prefix+"_x.bin",x)
 savebin(prefix+"_meta.bin",(status,obj,tol,actual,sign))
}
end
'''

def instrument(src,dst):
 shutil.copytree(src,dst)
 p=dst/'_dea_nonte_lp.ado';s=p.read_text().replace('real scalar _dea_nonte_optimize(', 'real scalar _r08_native(',1)
 s=s.replace('        (void) S.optimize()',
     '        (void) S.optimize()\n        savebin("attempt"+strofreal(st_numscalar("r08calls")+1,"%06.0f")+"_"+strofreal(attempt)+".bin",(tol,use_tol,S.converged(),S.errorcode()))')
 at=s.rindex('\nend');s=s[:at]+'''
real scalar _dea_nonte_optimize(class LinearProgram scalar S,real scalar tol, | real scalar cap) {
 real scalar status,n,sign
 real rowvector c,x
 real matrix E,U,B
 if (args()<3) cap=1e-8
 c=S.getCoefficients();n=cols(c);E=S.getEquality();U=S.getInequality();B=S.getBounds()
 if (rows(E)==0) E=J(0,n+1,0)
 if (rows(U)==0) U=J(0,n+1,0)
 sign=(S.getMaxOrMin()=="max" ? -1 : 1)
 status=_r08_native(S,tol,cap)
 x=(status==0 ? S.parameters() : J(1,n,.))
 r08save(c,E,U,B,x,status,(status==0 ? S.value() : .),tol,S.getTol(),sign)
 savebin("lp"+strofreal(st_numscalar("r08calls"),"%06.0f")+"_native.bin",(S.converged(),S.errorcode(),S.value()))
 if (S.converged() & S.errorcode()==0) savebin("lp"+strofreal(st_numscalar("r08calls"),"%06.0f")+"_candidate.bin",S.parameters())
 return(status)
}
'''+s[at:];p.write_text(s)
 p=dst/'_dea_simplex.ado';s=p.read_text().replace('struct RSM_Result scalar _dea_nonte_simplex(', 'struct RSM_Result scalar _r08_simplex(',1)
 at=s.rindex('\nend');s=s[:at]+'''
struct RSM_Result scalar _dea_nonte_simplex(real matrix U,real colvector bu,
 real matrix E,real colvector be,real rowvector c,real scalar tol) {
 struct RSM_Result scalar res
 res=_r08_simplex(U,bu,E,be,c,tol)
 r08save(c,(E,be),(U,bu),(J(1,cols(c),0)\\J(1,cols(c),.)),res.x',res.status,res.objval,tol,tol,1)
 return(res)
}
'''+s[at:];p.write_text(s)
 # Verify tolerance in composite/radial parts as well (does not alter calculation).
 for engine in ['ipm','simplex']:
  p=dst/f'_dea_{engine}.ado';s=p.read_text()
  targets=['envelopment'] if engine=='ipm' else ['multiplier']
  for name in targets:
   start=s.index(f'void _dea_{engine}_{name}(');pos=s.index('    X = st_data',start)
   s=s[:pos]+'    assert(tol==st_numscalar("expected_tol"))\n'+s[pos:]
  p.write_text(s)

def fixture(variant):
 s=DATA
 if variant=='duplicate':s+='expand 2 if id==2\nreplace id=_n\n'
 if variant=='units':s+='replace x1=x1*.01\nreplace y2=y2*100\nreplace b1=b1*.1\n'
 if variant=='tiny':
  # An explicitly dominated DMU exposes omitted denominator coefficients;
  # an all-efficient fixture alone would not detect the old score error.
  s+="foreach v of varlist x1 x2 b1 b2 {\n replace `v'=2*`v'[1] if id==3\n}\n"
  s+="foreach v of varlist y1 y2 {\n replace `v'=.5*`v'[1] if id==3\n}\n"
  s+="foreach v of varlist x1 x2 y1 y2 b1 b2 {\n replace `v'=1e-9*`v'\n}\n"
 s+='gen double w1=.5\ngen double w2=.3\ngen double p1=2\ngen double p2=1.5\n'
 s+='gen double gx1=x1\ngen double gx2=x2\ngen double gy1=y1\ngen double gy2=y2\n'
 s+='foreach v of varlist x1 x2 y1 y2 {\n gen double `v\'l=.8*`v\'\n gen double `v\'h=1.2*`v\'\n}\n'
 return s

def command(profile,engine,rts,tol,form='A',variant='regular'):
 model,opts,ori=PROFILES[profile];q='if id!=2 in 1/5' if variant=='sample' else ''
 head=f'dea x1 x2=y1 y2 {q},' if form=='A' else f'dea y1 y2 {q}, inputs(x1 x2)'
 orient=f'orientation({ori})' if ori else ''
 return f'{head} model({model}) {opts} {orient} rts({rts}) engine({engine}) tolerance({tol}) dmu(id) gen(res)'

def verify(d,tol):
 counts=dict(restarts=0,lp_checks=0,accepted=0,rejected=0,rejected_feasible=0,infeasible=0,unbounded=0,max_scaled=0.,max_relative_objective_gap=0.)
 for p in sorted(d.glob('lp*_meta.bin')):
  base=p.name[:-9];meta=readbin(p).ravel();status,obj,requested,actual,sign=meta
  assert requested==tol and actual<=tol,(meta,tol)
  c=readbin(d/(base+'_c.bin')).ravel();E=readbin(d/(base+'_E.bin'));U=readbin(d/(base+'_U.bin'));B=readbin(d/(base+'_B.bin'));x=readbin(d/(base+'_x.bin')).ravel()
  n=len(c);assert E.shape[1]==U.shape[1]==n+1
  attempts=sorted(d.glob('attempt'+base[2:]+'_*.bin'))
  if attempts:
   assert len(attempts)<=2
   first=readbin(attempts[0]).ravel();assert first[0]==tol and first[1] in (min(tol,1e-8),min(tol,1e-12))
   if len(attempts)==2:
    counts['restarts']+=1;last=readbin(attempts[1]).ravel()
    assert first[2]==1 and first[3]==0 and last[0]==tol and last[1]==min(tol,1e-12)
    assert first[1]>last[1]
  # HiGHS receives positive row rescaling, preserving each original constraint.
  def normalized(M):
   if not len(M):return None,None
   z=np.maximum(abs(M).max(axis=1),1e-300);return M[:,:-1]/z[:,None],M[:,-1]/z
  ae,be=normalized(E);au,bu=normalized(U)
  lb=np.nan_to_num(B[0],nan=-np.inf);ub=np.nan_to_num(B[1],nan=np.inf)
  cs=max(1,np.max(abs(c)))
  fit=linprog(sign*c/cs,A_eq=ae,b_eq=be,A_ub=au,b_ub=bu,bounds=list(zip(lb,ub)),method='highs',options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
  assert fit.status in (0,2,3),fit.message
  counts['lp_checks']+=1
  if fit.status==2:counts['infeasible']+=1
  if fit.status==3:counts['unbounded']+=1
  if status==0:
   counts['accepted']+=1;assert fit.success,(p.name,meta,fit.message)
   assert np.isfinite(x).all() and np.isfinite(obj)
   r=residual(x,E[:,:-1],E[:,-1],U[:,:-1],U[:,-1],lb,ub)
   scaled=max(r[3:]);assert scaled<=tol+5e-16,(p.name,scaled,tol)
   gap=abs(obj-sign*fit.fun*cs)/max(1,abs(obj),abs(fit.fun*cs))
   assert gap<=1e-6,(p.name,obj,fit.fun*cs,gap)
   assert abs(obj-c@x)<=1e-8*max(1,abs(obj)),(p.name,obj,c@x)
   counts['max_scaled']=max(counts['max_scaled'],float(scaled));counts['max_relative_objective_gap']=max(counts['max_relative_objective_gap'],float(gap))
  else:
   counts['rejected']+=1;counts['rejected_feasible']+=int(fit.success)
   if fit.success and (d/(base+'_native.bin')).exists():
    info=readbin(d/(base+'_native.bin')).ravel().tolist()
    if (d/(base+'_candidate.bin')).exists():
     raw=readbin(d/(base+'_candidate.bin')).ravel();info+=residual(raw,E[:,:-1],E[:,-1],U[:,:-1],U[:,-1],lb,ub).tolist()
    counts.setdefault('rejected_feasible_evidence',[]).append(dict(lp=base,native=info))
 return counts

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',required=True,type=Path);ap.add_argument('--only');a=ap.parse_args()
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir();hashes={}
 for f in package_files(REPO):
  raw=(REPO/f).read_bytes();(pkg/Path(f).name).write_bytes(raw);hashes[f]=hashlib.sha256(raw).hexdigest()
 audit=out/'audit_ado';instrument(pkg,audit)
 (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,stata=a.stata),indent=2)+'\n');rows=[]
 def run(name,body,source=audit,check=None):
  if a.only and not re.search(a.only,name):return
  d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 120\n'
  for f in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/f).mkdir();head+=f'sysdir set {f} "{d/f}"\n'
  head+=f'adopath ++ "{source}"\n'+WRITER+CAPTURE+'scalar r08calls=0\n'
  (d/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n');r=dict(test=name,passed=False);t=time.monotonic()
  try:
   proc=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=180)
   if proc.returncode or not (d/'test.log').exists():
    (d/'process.json').write_text(json.dumps(dict(returncode=proc.returncode,stdout=proc.stdout.decode(errors='replace'),stderr=proc.stderr.decode(errors='replace')),indent=2)+'\n')
   log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else proc.stderr.decode(errors='replace')
   assert proc.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-1800:]
   if check:r.update(check(d))
   r['passed']=True
  except Exception as e:r['error']=str(e)
  r['seconds']=round(time.monotonic()-t,3);rows.append(r);(out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(r),flush=True)
 def modelcase(profile,engine,rts,tol,form='A',variant='regular'):
  name=f'lp_{profile}_{engine}_{rts}_{tol}_{form}_{variant}'
  body=fixture(variant)+f'scalar expected_tol={tol}\n'+command(profile,engine,rts,tol,form,variant)+'\n'
  body+=f'assert r(tolerance)=={tol}\nscalar reported_failed=r(n_failed)\n'
  if profile!='congestion':
   if profile.startswith('interval'):body+='quietly count if '+('id!=2 & inrange(_n,1,5) & ' if variant=='sample' else '')+'missing(res,res_opt)\n'
   elif variant=='sample':body+='quietly count if id!=2 & inrange(_n,1,5) & missing(res)\n'
   else:body+='quietly count if missing(res)\n'
   body+='assert reported_failed==r(N)\n'
  body+='mata: savebin("scores.bin",st_data(.,"res"))\n'
  # All new LP routes must be observed; simplex multiplier is already validated by TE.
  if not (engine!='ipm' and profile.startswith('mult')):body+='assert scalar(r08calls)>0\n'
  run(name,body,check=lambda d:verify(d,float(tol)))
 for profile,engine,rts,tol in itertools.product(PROFILES,['simplex','ipm','auto'],['crs','vrs','drs','irs'],['1e-4','1e-7','1e-10']):
  if profile=='congestion' and rts!='vrs':continue
  if profile in ('profit','sbmu') and rts not in ('crs','vrs'):continue
  if profile=='profit' and rts=='crs':continue
  modelcase(profile,engine,rts,tol)
 for profile,engine,variant in itertools.product(PROFILES,['simplex','ipm'],['units','duplicate','sample']):
  modelcase(profile,engine,'vrs','1e-9','B',variant)
 # A formerly tolerance-dependent formula: use exact same positive ratios in
 # much smaller units. The independent model builder uses every positive value.
 for rts,tol in itertools.product(['crs','vrs'],['1e-7','1e-10']):
  body=fixture('tiny')+f'scalar expected_tol={tol}\n'+command('sbmu','ipm',rts,tol)+'\n'
  body+='mata: savebin("scores.bin",st_data(.,"res"))\nmata: savebin("data.bin",st_data(.,tokens("x1 x2 y1 y2 b1 b2")))\n'
  def smallcheck(d,rts=rts,tol=tol):
   # Independent LP in original ratios but well-conditioned units.
   data=readbin(d/'data.bin')*1e9;actual=readbin(d/'scores.bin').ravel();expected=np.array([specialized_lp(data,'sbmundesirable',rts,k).fun for k in range(len(data))])
   assert np.isfinite(actual).all() and np.max(abs(actual-expected))<1e-6,(actual,expected)
   return dict(independent_formula_checks=len(data),max_score_gap=float(np.max(abs(actual-expected))))
  run(f'tiny_sbmu_{rts}_{tol}',body,check=smallcheck)
 # Primal rejection is exercised independently of native convergence flags.
 for fault in ['bound','equality','nonfinite']:
  target=out/('fault_'+fault);shutil.copytree(pkg,target);p=target/'_dea_nonte_lp.ado';s=p.read_text();anchor='        x=S.parameters()'
  change={'bound':'x[1]=-1000','equality':'x[1]=x[1]+1000','nonfinite':'x[1]=.'}[fault]
  assert s.count(anchor)==1;s=s.replace(anchor,anchor+'\n    '+change);p.write_text(s)
  for profile in ['cost','revenue','profit','ddf','additive','sbm','sbmin','sbmout','superin','superout','supersbm','bad_yb','intervalin','intervalout','sbmu']:
   run(f'reject_{fault}_{profile}',fixture('regular')+command(profile,'ipm','vrs','1e-7')+'\nassert r(n_failed)==_N\nassert missing(res)\n',source=target)
 # Prove that passive LP recording did not change ordinary output values.
 for profile,engine in itertools.product(PROFILES,['simplex','ipm']):
  regular=f'lp_{profile}_{engine}_vrs_1e-7_A_regular'
  def parity(d,regular=regular):
   left=readbin(d/'scores.bin');right=readbin(out/regular/'scores.bin')
   np.testing.assert_array_equal(left,right)
   return dict(recording_parity=True)
  if not a.only:
   run(f'parity_{profile}_{engine}',fixture('regular')+command(profile,engine,'vrs','1e-7')+'\nmata: savebin("scores.bin",st_data(.,"res"))\n',source=pkg,check=parity)
 # A corrupted simplex result must fail the original-row/bound gate as well.
 for fault,assignment in [('bound','res.x[1]=-1000'),('equality','res.x[1]=res.x[1]+1000'),('nonfinite','res.x[1]=.')]:
  target=out/('simplex_fault_'+fault);shutil.copytree(pkg,target)
  p=target/'_dea_simplex.ado';s=p.read_text();anchor='    res = _rsm_solve_lp(A_iq,b_iq,A_eq,b_eq,c,tol)'
  assert s.count(anchor)==1;s=s.replace(anchor,anchor+'\n    if (res.status==0) '+assignment);p.write_text(s)
  for profile in ['cost','revenue','profit','ddf','additive','sbm','sbmin','sbmout','superin','superout','supersbm','alloc_in','alloc_out']:
   run(f'simplex_reject_{fault}_{profile}',fixture('regular')+command(profile,'simplex','vrs','1e-7')+'\nassert r(n_failed)==_N\nassert missing(res)\n',source=target)
 # A 1e-9 multiplier is essential when data are 1e9. Cleanup must not zero it.
 for engine,rts,ori in itertools.product(['simplex','ipm','auto'],['crs','vrs','drs','irs'],['input','output']):
  body='input double id x1 y1\n1 1e9 1e9\n2 2e9 1e9\nend\n'
  body+=f'dea x1=y1,model(multiplier) orientation({ori}) rts({rts}) engine({engine}) tolerance(1e-9) saveweights gen(res) dmu(id)\n'
  if ori=='input':body+='assert abs(_dea_v1*x1-1)<1e-7\n'
  else:body+='assert abs(_dea_u1*y1-1)<1e-7\n'
  body+='assert !missing(_dea_v1,_dea_u1,_dea_u0)\n'
  body+="forvalues k=1/2 {\n assert _dea_u1[`k']*y1-_dea_v1[`k']*x1-_dea_u0[`k']<1e-7\n}\n"
  run(f'multiplier_small_{engine}_{rts}_{ori}',body,source=pkg)
 target=out/'weights_fail';shutil.copytree(pkg,target);p=target/'_dea_nonte_lp.ado';s=p.read_text();s=s.replace('    c=S.getCoefficients(); n=cols(c)','    return(1) // test-only failure of multiplier LP\n    c=S.getCoefficients(); n=cols(c)');p.write_text(s)
 for rts,ori in itertools.product(['crs','vrs','drs','irs'],['input','output']):
  body=fixture('regular')+command('multin' if ori=='input' else 'multout','ipm',rts,'1e-7')+'\n'
  body+='assert !missing(res)\nassert missing(_dea_v1,_dea_v2,_dea_u1,_dea_u2,_dea_u0)\n'
  run(f'multiplier_failure_{rts}_{ori}',body,source=target)
 # Existing cached engine signatures must reload when the core sees version 1.
 for engine in ['simplex','ipm']:
  body=fixture('regular')+f'quietly run "{pkg}/_dea_{engine}.ado"\n'
  body+=f'mata: mata drop _dea_{engine}_audit_ready()\nmata: real scalar _dea_{engine}_audit_ready() return(1)\n'
  marker=re.search(r'real scalar _dea_'+engine+r'_audit_ready\(\) return\((\d+)\)', (pkg/('_dea_'+engine+'.ado')).read_text())
  assert marker, 'Engine readiness marker is missing'
  expected_ready=int(marker.group(1))
  body+=command('cost',engine,'vrs','1e-9')+f'\nmata: assert(_dea_{engine}_audit_ready()=={expected_ready})\nassert !missing(res)\n'
  run('reload_'+engine,body,source=pkg)
 summary=dict(total=len(rows),passed=sum(r['passed'] for r in rows),restarts=sum(r.get('restarts',0) for r in rows),lp_checks=sum(r.get('lp_checks',0) for r in rows),accepted=sum(r.get('accepted',0) for r in rows),rejected_feasible=sum(r.get('rejected_feasible',0) for r in rows),infeasible=sum(r.get('infeasible',0) for r in rows),max_scaled=max((r.get('max_scaled',0) for r in rows),default=0),max_relative_objective_gap=max((r.get('max_relative_objective_gap',0) for r in rows),default=0))
 (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary));raise SystemExit(summary['passed']!=summary['total'])
if __name__=='__main__':main()
