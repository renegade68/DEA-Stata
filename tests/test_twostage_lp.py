#!/usr/bin/env python3
"""Actual two-stage LP capture and forced-stage-failure data/return contracts."""
from _support import package_files
from pathlib import Path
import argparse,hashlib,itertools,json,re,shutil,subprocess,time
import numpy as np
from test_nonte_tolerance import instrument,verify,CAPTURE
from test_core_retries import WRITER
from test_twostage_models import fixture,PUBLISHED_DATA
REPO=Path(__file__).resolve().parents[1]
DATA='input double(id x z1 z2 y)\n1 1 1 2 1\n2 2 3 1 2\n3 3 2 4 3\n4 2 2 2 1\nend\n'
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
  head+=f'adopath ++ "{source}"\n'+WRITER+CAPTURE+'scalar r08calls=0\nscalar stagecalls=0\n';(d/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n');r=dict(test=name,passed=False);t=time.monotonic()
  try:
   p=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=180);log=(d/'test.log').read_text(errors='replace');assert p.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2400:]
   if check:r.update(check(d))
   r['passed']=True
  except Exception as e:r['error']=repr(e)
  r['seconds']=round(time.monotonic()-t,3);rows.append(r);(out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(r),flush=True)
 profiles=[('network',method,rts,ort,'stage1') for method,rts,ort in itertools.product(['multiplicative','independent','additive'],['crs','vrs'],['input','output'])]
 profiles += [('stackelberg',method,'crs','input',l) for method,l in itertools.product(['independent','noncooperative'],['stage1','stage2'])]+[('stackelberg','independent','vrs','input',l) for l in ['stage1','stage2']]
 for profile,tol,var in itertools.product(profiles,['1e-4','1e-7','1e-10'],['base','units','zerocol']):
  kind,method,rts,ort,l=profile;z,cols=fixture((2,2,2),var);f=out/('data_'+var+'.csv')
  if not f.exists():np.savetxt(f,z,delimiter=',',header=','.join(cols),comments='',fmt='%.17g')
  body=f'import delimited using "{f}",clear asdouble\ndea_{kind},inputs(x1 x2) intermediates(z1 z2) outputs(y1 y2) method({method}) rts({rts}) '+(f'ort({ort})' if kind=='network' else f'leader({l})')+f' gen(res) tol({tol}) matrix detail\nassert r08calls>0\n'
  run('lp_'+'_'.join(profile)+'_'+tol+'_'+var,body,check=lambda d,tol=tol:verify(d,float(tol)))
 (out/'published.csv').write_text(PUBLISHED_DATA)
 for kind,method,l in [('network','multiplicative','stage1'),('stackelberg','noncooperative','stage1'),('stackelberg','noncooperative','stage2')]:
  for tol in ['1e-4','1e-7','1e-10']:
   body=f'import delimited using "{out/"published.csv"}",clear asdouble\ndea_{kind},inputs(x1 x2) intermediates(z1 z2) outputs(y1 y2) method({method}) '+(f'leader({l})' if kind=='stackelberg' else '')+f' gen(res) tol({tol}) matrix detail\nassert r08calls>0\n'
   run('published_'+kind+'_'+l+'_'+tol,body,check=lambda d,tol=tol:verify(d,float(tol)))
 for fault in ['primary','decomp','follower','stage1','stage2','bound','equality','nonfinite']:
  target=out/('fault_ado_'+fault);shutil.copytree(pkg,target)
  if fault in ['bound','equality','nonfinite']:
   p=target/'_dea_nonte_lp.ado';s=p.read_text();change={'bound':'x[1]=-1000','equality':'x=x:+1000','nonfinite':'x[1]=.'}[fault];s=s.replace('        x=S.parameters()','        x=S.parameters()\n        '+change);p.write_text(s)
  else:
   p=target/'_dea_twostage.ado';s=p.read_text()
   if fault in ['primary','decomp','follower']:
    needle='    weights=J(1,p,.)';cond={'primary':'q==0','decomp':'q==1','follower':'q>=2'}[fault];s=s.replace(needle,needle+'\n    if('+cond+') return(1)')
   else:
    needle='n=rows(X);m=cols(X);s=cols(Y);lambda=J(1,n,.)';s=s.replace(needle,needle+'\n    st_numscalar("stagecalls",st_numscalar("stagecalls")+1)\n    if(mod(st_numscalar("stagecalls"),2)=='+('1' if fault=='stage1' else '0')+') return(.)')
   p.write_text(s)
  for kind,method,leader in [('network','multiplicative','stage1'),('network','independent','stage1'),('stackelberg','independent','stage1'),('stackelberg','noncooperative','stage1'),('stackelberg','noncooperative','stage2')]:
   rel=kind=='network' and method=='multiplicative';nc=method=='noncooperative'
   body=DATA+f'dea_{kind},inputs(x) intermediates(z1 z2) outputs(y) method({method}) gen(res) matrix detail '+(f'leader({leader})' if kind=='stackelberg' else '')+'\n'
   if fault in ['bound','equality','nonfinite'] or fault=='primary' and rel:
    body+='assert missing(res_overall)\nassert r(n_failed_overall)==4\n'
   elif fault=='decomp' and rel:body+='assert r(n_decomp_fallback)==4\nassert !missing(res_overall)\nmata: assert(all(st_matrix("r(decomposition_status)"):==1))\n'
   elif fault=='follower' and nc:body+='assert r(n_failed_follower_lp)==4\nassert missing(res_overall,res_follower)\nassert !missing(res_leader,res_indep1,res_indep2)\nmata: assert(all(st_matrix("r(follower_status)"):==1))\n'
   elif fault in ['stage1','stage2'] and not rel:
    j=1 if fault=='stage1' else 2
    if kind=='network':body+=f'assert missing(res_stage{j},res_overall)\nassert !missing(res_stage{3-j})\n'
    else:
     body+=f'assert missing(res_indep{j})\nassert !missing(res_indep{3-j})\n'
     if nc and leader!=fault:body+='assert !missing(res_leader,res_follower,res_overall)\n'
     elif nc:body+='assert missing(res_leader,res_follower,res_overall)\nmata: assert(all(st_matrix("r(follower_status)"):==0))\n'
     else:body+='assert missing(res_overall)\n'
   else:body+='assert !missing(res_overall)\n'
   run('fault_'+fault+'_'+kind+'_'+method+'_'+leader,body,source=target)
 summary=dict(total=len(rows),passed=sum(r['passed'] for r in rows),failed=sum(not r['passed'] for r in rows),LPs=sum(r.get('lp_checks',0) for r in rows),rejected_feasible=sum(r.get('rejected_feasible',0) for r in rows));(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary);return bool(summary['failed'])
if __name__=='__main__':raise SystemExit(main())
