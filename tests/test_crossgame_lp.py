#!/usr/bin/env python3
"""Capture crossgame LPs and verify failures cannot become partial Nash updates."""
from _support import package_files
from pathlib import Path
import argparse,csv,hashlib,itertools,json,re,shutil,subprocess,time
from test_nonte_tolerance import instrument,verify,CAPTURE
from test_crosseff_models import fixture
from test_core_retries import WRITER
REPO=Path(__file__).resolve().parents[1]
DATA='input double(id x y)\n1 1 1\n2 2 1\n3 3 2\nend\n'
PAPER='input double(id x1 x2 x3 y1 y2)\n1 7 7 7 4 4\n2 5 9 7 7 7\n3 4 6 5 5 7\n4 5 9 8 6 2\n5 6 8 5 3 6\nend\n'
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--only');a=ap.parse_args()
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir();h={}
 for n in package_files(REPO):
  raw=(REPO/n).read_bytes();(pkg/Path(n).name).write_bytes(raw);h[n]=hashlib.sha256(raw).hexdigest()
 audit=out/'audit_ado';instrument(pkg,audit);(out/'manifest.json').write_text(json.dumps(dict(source_sha256=h,stata=a.stata),indent=2)+'\n');rows=[]
 def run(name,body,source=audit,check=None):
  if a.only and not re.search(a.only,name):return
  d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset processors 1\nset type double\nset linesize 160\n'
  for f in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/f).mkdir();head+=f'sysdir set {f} "{d/f}"\n'
  head+=f'adopath ++ "{source}"\n'+WRITER+CAPTURE+'scalar r08calls=0\nscalar gamecalls=0\n';(d/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n');r=dict(test=name,passed=False);t=time.monotonic()
  try:
   p=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=180);log=(d/'test.log').read_text(errors='replace');assert p.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2400:]
   if check:r.update(check(d))
   r['passed']=True
  except Exception as e:r['error']=repr(e)
  r['seconds']=round(time.monotonic()-t,3);rows.append(r);(out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(r),flush=True)
 for strategy,tol,variant in itertools.product(['all','aggressive','benevolent','nash'],['1e-4','1e-7','1e-10'],['base','units','zero']):
  body=DATA+'gen double z=0\n'
  if variant=='units':body+='replace x=x*1e9\nreplace y=y*1e-10\n'
  xs='x z' if variant=='zero' else 'x'
  body+=f'dea_crossgame {xs}=y,gen(res) strategy({strategy}) tol({tol}) matrix detail\nassert r08calls>0\n'
  run(f'lp_{strategy}_{tol}_{variant}',body,check=lambda d,tol=tol:verify(d,float(tol)))
 # Natural multivariate boundaries complement the simple analytic submission grid.
 for m,s,var,tol in itertools.product([2],[1,2],['base','tiny'],['1e-6','1e-9']):
  z,cols=fixture(m,s,var);f=out/f'boundary_{m}x{s}_{var}.csv'
  with f.open('w') as stream:w=csv.writer(stream);w.writerow(cols);w.writerows(z)
  xs=' '.join(f'x{k+1}' for k in range(m));ys=' '.join(f'y{k+1}' for k in range(s))
  body=f'import delimited using "{f}",clear asdouble\ndea_crossgame {xs}={ys},gen(res) tol({tol}) matrix detail maxiter(25)\nassert r08calls>0\n'
  run(f'boundary_{m}x{s}_{var}_{tol}',body,check=lambda d,tol=tol:verify(d,float(tol)))
 for fault in ['primary','first_primary','secondary','first_game','second_game','bound','equality','nonfinite']:
  target=out/('fault_ado_'+fault);shutil.copytree(pkg,target)
  if fault in ['bound','equality','nonfinite']:
   p=target/'_dea_nonte_lp.ado';s=p.read_text();change={'bound':'x[1]=-1000','equality':'x=x:+1000','nonfinite':'x[1]=.'}[fault];s=s.replace('        x=S.parameters()','        x=S.parameters()\n        '+change);p.write_text(s)
  else:
   p=target/'_dea_crossgame_core.ado';s=p.read_text();needle='    v=J(1,m,.);u=J(1,s,.);value=.'
   cond={'primary':'mode==0','first_primary':'mode==0 & j==1','secondary':'mode==1 | mode==2','first_game':'mode==3 & d==1 & j==1','second_game':'mode==3 & st_numscalar("gamecalls")==rows(X0)^2+1'}[fault]
   assert needle in s;s=s.replace(needle,needle+'\n    if(mode==3) st_numscalar("gamecalls",st_numscalar("gamecalls")+1)\n    if('+cond+') return(1)',1);p.write_text(s)
  # These cases test final unavailability, so make the recovery unavailable too.
  # Native-only failures and successful certified recovery have a separate suite.
  if fault in ['secondary','first_game','second_game']:
   p=target/'_dea_crossgame_recovery.ado';s=p.read_text()
   needle=' n=rows(X);m=cols(X);s=cols(Y);primary_value=.;rc=1';assert s.count(needle)==1
   s=s.replace(needle,needle+'\n diag[15]=0;diag[16]=3;return(1);',1);p.write_text(s)
  body=PAPER+'dea_crossgame x1 x2 x3=y1 y2,gen(res) matrix detail tol(1e-8)\n'
  if fault in ['primary','bound','equality','nonfinite']:
   body+='assert r(n_failed_primary)==5 & r(nash_status)==2 & r(nash_converged)==0 & r(nash_attempts)==0\nassert r(nash_initial_missing)==25\nassert missing(res_agg,res_ben,res_nash)\n'
  elif fault=='first_primary':body+='assert r(n_failed_primary)==1 & r(nash_status)==2 & r(nash_attempts)==0\nassert missing(res_nash)\nmata: assert(st_matrix("r(status_agg)")[1]==0)\n'
  elif fault=='secondary':body+='assert r(n_failed_agg)==5 & r(n_failed_ben)==5\nmata: assert(all(st_matrix("r(weights_input)"):==st_matrix("r(weights_agg_input)")))\nmata: assert(all(st_matrix("r(weights_output)"):==st_matrix("r(weights_ben_output)")))\nmata: assert(all(st_matrix("r(status_agg)"):==3) & all(st_matrix("r(status_ben)"):==3))\n'
  else:
   it=0 if fault=='first_game' else 1
   body+=f'assert r(nash_status)==3 & r(nash_converged)==0 & r(nash_iterations)=={it} & r(nash_attempts)=={it+1} & r(nash_failed_pairs)==1\n'
   body+='mata: H=st_matrix("r(nash_history)"); T=st_matrix("r(nash_trace)"); C=st_matrix("r(ce_nash)" )\nmata: assert(max(abs(st_data(.,"res_nash")-H[rows(H),.]\'))<1e-12)\nmata: assert(max(abs(colsum(C)/rows(C)-H[rows(H),.]))<1e-12)\nmata: assert(T[rows(T),3]==1 & T[rows(T),4]==0)\n'
   if fault=='first_game':body+='mata: assert(all(C:==st_matrix("r(ce_initial)")))\n'
   else:body+='mata: assert(max(abs(H[1,.]-H[2,.]))>1e-6)\n'
  run('fault_'+fault,body,source=target)
 summary=dict(total=len(rows),passed=sum(r['passed'] for r in rows),failed=sum(not r['passed'] for r in rows),LPs=sum(r.get('lp_checks',0) for r in rows),rejected_feasible=sum(r.get('rejected_feasible',0) for r in rows));(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary);return bool(summary['failed'])
if __name__=='__main__':raise SystemExit(main())
