#!/usr/bin/env python3
"""Audit cross-efficiency submitted LPs and exercise explicit primary fallback."""
from _support import package_files
from pathlib import Path
import argparse,hashlib,itertools,json,re,shutil,subprocess,time
from test_nonte_tolerance import instrument,verify,CAPTURE
from test_core_retries import WRITER
REPO=Path(__file__).resolve().parents[1]
DATA='''input double(id x1 x2 y1 y2)
1 2 4 2 1
2 4 2 1 2
3 3 3 1 1
4 4 5 3 2
5 6 4 2 3
end
'''
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--only');a=ap.parse_args()
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir();h={}
 for n in package_files(REPO):
  raw=(REPO/n).read_bytes();(pkg/Path(n).name).write_bytes(raw);h[n]=hashlib.sha256(raw).hexdigest()
 audit=out/'audit_ado';instrument(pkg,audit);(out/'manifest.json').write_text(json.dumps(dict(source_sha256=h,stata=a.stata),indent=2)+'\n');rows=[]
 def run(name,body,source=audit,check=None):
  if a.only and not re.search(a.only,name):return
  d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset processors 1\nset type double\n'
  for f in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/f).mkdir();head+=f'sysdir set {f} "{d/f}"\n'
  head+=f'adopath ++ "{source}"\n'+WRITER+CAPTURE+'scalar r08calls=0\nscalar fault_calls=0\n';(d/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n');r=dict(test=name,passed=False);t=time.monotonic()
  try:
   p=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=180);log=(d/'test.log').read_text(errors='replace');assert p.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2000:]
   if check:r.update(check(d))
   r['passed']=True
  except Exception as e:r['error']=repr(e)
  r['seconds']=round(time.monotonic()-t,3);rows.append(r);(out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(r),flush=True)
 for rts,ori,method,tol,variant in itertools.product(['crs','vrs'],['input','output'],['arbitrary','benevolent','aggressive'],['1e-4','1e-7','1e-10'],['base','units','zero']):
  body=DATA+f'scalar expected_tol={tol}\n'
  if variant=='units':body+='replace x1=x1*1e9\nreplace x2=x2*1e-9\nreplace y1=y1*1e-8\nreplace y2=y2*1e10\n'
  if variant=='zero':body+='replace x1=0 if id==1\nreplace y1=0 if id==2\n'
  body+=f'dea_crosseff x1 x2=y1 y2,dmu(id) gen(res) rts({rts}) ort({ori}) method({method}) tol({tol}) matrix\nassert r(tolerance)=={tol}\nassert r08calls>0\n'
  run(f'lp_{rts}_{ori}_{method}_{tol}_{variant}',body,check=lambda d,tol=tol:verify(d,float(tol)))
 for fault in ['primary','first_primary','secondary','first_secondary','bound','equality','nonfinite']:
  target=out/('fault_'+fault);shutil.copytree(pkg,target);p=target/'_dea_nonte_lp.ado';s=p.read_text()
  if fault in ['bound','equality','nonfinite']:
   change={'bound':'x[1]=-1000','equality':'x=x:+1000','nonfinite':'x[1]=.'}[fault]
   s=s.replace('        x=S.parameters()','        x=S.parameters()\n        '+change)
  else:
   cond={'primary':'rows(S.getEquality())==1','first_primary':'st_numscalar("fault_calls")==1','secondary':'rows(S.getEquality())==2','first_secondary':'st_numscalar("fault_calls")==2'}[fault]
   s=s.replace('    c=S.getCoefficients(); n=cols(c)','    st_numscalar("fault_calls",st_numscalar("fault_calls")+1)\n    if ('+cond+') return(1)\n    c=S.getCoefficients(); n=cols(c)')
  p.write_text(s)
  for ori in ['input','output']:
   cmd=f'dea_crosseff x1 x2=y1 y2,dmu(id) gen(res) rts(crs) ort({ori}) matrix'
   body=DATA
   if fault in ['secondary','first_secondary']:
    body+='scalar fault_calls=-1000\n'
    body+=cmd+'\nmatrix expected_v=r(weights_input)\nmatrix expected_u=r(weights_output)\nmatrix expected_self=J(5,1,.)\nmata: st_matrix("expected_self",st_data(.,"res_self"))\nscalar fault_calls=0\n'
    body+=cmd+' method(benevolent) replacevars\n'
    if fault=='secondary':
     body+='assert r(n_failed_primary)==0 & r(n_failed_secondary)==5\nmata: assert(all(st_matrix("r(stage_status)"):==3))\nmata: assert(all(st_matrix("expected_v"):==st_matrix("r(weights_input)")))\nmata: assert(all(st_matrix("expected_u"):==st_matrix("r(weights_output)")))\n'
    else:body+='assert r(n_failed_secondary)>=1\nmata: assert(st_matrix("r(stage_status)")[1]==3)\n'
    body+='mata: assert(all(st_matrix("expected_self"):==st_data(.,"res_self")))\nassert !missing(res)\n'
   else:
    body+=cmd+' method(arbitrary)\n'
    if fault=='first_primary':body+='assert r(n_failed_primary)==1 & r(n_failed_ce)==0 & r(n_missing_peer_pairs)==4\nassert missing(res_self,res_maverick) if id==1\nassert !missing(res)\nmata: assert(st_matrix("r(peer_count)")[1]==4 & all(st_matrix("r(peer_count)")[2..5]:==3))\n'
    else:body+='assert r(n_failed_primary)==5 & r(n_failed_secondary)==0 & r(n_failed_ce)==5 & r(n_missing_peer_pairs)==20\nassert missing(res,res_self,res_maverick)\nassert missing(r(mean_ce),r(mean_self))\n'
   run(f'fault_{fault}_{ori}',body,source=target)
 summary=dict(total=len(rows),passed=sum(r['passed'] for r in rows),failed=sum(not r['passed'] for r in rows),LPs=sum(r.get('lp_checks',0) for r in rows),rejected_feasible=sum(r.get('rejected_feasible',0) for r in rows));(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary);return bool(summary['failed'])
if __name__=='__main__':raise SystemExit(main())
