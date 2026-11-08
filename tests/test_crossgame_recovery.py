#!/usr/bin/env python3
"""Crossgame recovery contract, output protection and shared-state regression.

Fixed synthetic boundary references; raw logs stay in the requested output folder.
Numerical and native-failure references were verified on Stata/MP16.1 macOS.
"""
from pathlib import Path
import argparse,csv,hashlib,json,shutil,subprocess,time
import numpy as np
from _support import package_files
from test_core_retries import WRITER,readbin
from test_output_protection import unchanged
from test_crossgame_models import MATS
REPO=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--stata',required=True)
ap.add_argument('--output',type=Path,required=True)
a=ap.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
pkg=out/'ado';pkg.mkdir()
for name in package_files(REPO):shutil.copy2(REPO/name,pkg/Path(name).name)
STATA=a.stata;rows=[]
spec=json.loads((REPO/'tests/fixtures/crossgame_recovery.json').read_text())['cases']
def expected(case,name):return np.array(spec[case][name],float)
def dump(name,value):(out/name).write_text(json.dumps(value,indent=2)+'\n')
def csvdata(case):
 import io
 z=expected(case,'data');f=io.StringIO();w=csv.writer(f)
 w.writerow(['id']+(['x1','x2','y1'] if case=='game' else ['x1','y1']))
 w.writerows([[format(v,'.17g') for v in row] for row in z]);return f.getvalue()
data=csvdata('game')
check=r'''
mata:
void check_state(transmorphic scalar saved) {
 string rowvector names
 pointer scalar p
 real scalar i
 names=_c57_context_names()
 for(i=1;i<=cols(names);i++) if(asarray(saved,names[i]+":existed")) {
  p=findexternal(names[i]);assert(p!=NULL)
  assert(rows(*p)==rows(asarray(saved,names[i])) & cols(*p)==cols(asarray(saved,names[i])))
  if(rows(*p)*cols(*p)>0) assert(all(*p:==asarray(saved,names[i])))
 }
}
end
mata: _dea_rec_abort();_rsm_audit_reset(1);_rsm_success_reset(1,1);_dea_link_reset(2,1)
mata: snapshot=_c57_context_save()
'''
cmd='dea_crossgame x1 x2=y1,dmu(id) gen(res) tol(1e-9) maxiter(25) matrix detail'
get='''matrix C=r(recovery_diagnostics)
matrix N=(r(n_recovery_attempted),r(n_recovered),r(n_recovery_unresolved),r(n_recovery_limit),r(n_recovery_rows),r(n_recovery_rows_omitted))
mata: savebin("C.bin",st_matrix("C"));savebin("N.bin",st_matrix("N"))
'''
def run(name,body,checker=None,package=pkg,init='import delimited data.csv,clear asdouble\n'):
 d=out/name;d.mkdir();(d/'data.csv').write_text(data);head='version 16\nclear all\nset more off\nset processors 1\nset type double\nset linesize 240\n'
 for k in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/k).mkdir();head+=f'sysdir set {k} "{d/k}"\n'
 head+=f'adopath ++ "{package}"\n'+WRITER
 (d/'test.do').write_text(head+init+body+'\ndi "CONTRACT_PASS"\n');t=time.monotonic();p=subprocess.run([STATA,'-b','do','test.do'],cwd=d,capture_output=True,timeout=180);log=(d/'test.log').read_text(errors='replace');r=dict(test=name,passed=False)
 try:
  assert p.returncode==0 and 'CONTRACT_PASS' in log.splitlines(),log[-2400:]
  if checker:checker(d)
  assert not list(d.glob('recover*.bin')) and not list(d.glob('lp*.bin'))
  r.update(passed=True,seconds=round(time.monotonic()-t,3))
 except Exception as e:r['error']=str(e)
 rows.append(r);dump('results.json',rows);print(json.dumps(r),flush=True)
 return r['passed']
# Actual native error provenance and complete Nash acceptance.
def parity(d):
 np.testing.assert_array_equal(readbin(d/'C.bin'),expected('game','recovery_diagnostics'))
 np.testing.assert_array_equal(readbin(d/'N.bin'),[[1,1,0,0,1,0]])
run('cold',cmd+'\nassert r(nash_status)==0 & r(n_recovered_game)==1\n'+get,parity)
# Return contract exists independently of display switches, and never leaves work matrices.
run('without_matrix',cmd.replace(' matrix detail','')+'\nassert r(recovery_version)==11401 & r(n_recovered_game)==1\nassert colsof(r(recovery_diagnostics))==32\ncapture matrix list r(ce_nash)\nassert _rc==111\n'+get,parity)
run('warm_context',cmd+'\n'+check+'scalar g113calls=991\nmatrix sentinel=(111,222)\nlocal before : all matrices\n'+cmd+' replacevars\n'+get+'mata: check_state(snapshot)\nassert g113calls==991\nmatrix drop C N\nlocal after : all matrices\nassert "`before\'"=="`after\'"\n',parity)
# A second strategy cannot retain prior game records or output matrices.
run('strategy_reset',cmd+'\n'+cmd.replace(' matrix detail',' strategy(aggressive)')+' replacevars\nassert r(n_recovered)==0 & r(n_recovery_attempted)==0 & r(n_recovery_rows)==0\ncapture matrix list r(recovery_diagnostics)\nassert _rc==111\ncapture matrix list r(ce_nash)\nassert _rc==111\n')
run('string_ids_subset','tostring id,replace\nreplace id="" in 1\n'+cmd.replace(' x1 x2=y1,',' in 2/6,inputs(x1 x2) outputs(y1)')+'\nassert r(n)==5\nassert missing(res_nash) in 1\n')
# Real sample-index and original-observation mapping after exclusions, including a recovered LP.
run('observation_mapping','expand 2\nreplace x1=-1 if _n>6\n'+cmd.replace(' x1 x2=y1,',' x1 x2=y1 if _n<=6,')+'\nassert r(n)==6 & r(n_recovered_game)==1\nassert missing(res_nash) if _n>6\n'+get,parity)
# Repeated loading and recovery after partial helper loss must preserve numerical results.
for name,action in [('main_reload',f'run "{pkg}/dea_crossgame.ado"\nrun "{pkg}/dea_crossgame.ado"\n'),('mata_clear','mata: mata clear\n'+WRITER)]+[(x,'mata: mata drop '+x+'()\n') for x in ['_dea_cg_lp','_g114_guard','_c58_fit','_c70_build','_c84_dot','_rsm_solve_lp']]:
 run('reload_'+name,cmd+'\n'+action+cmd+' replacevars\n'+get,parity)
# Core and crosseff coexistence; the saved shared diagnostic state survives crossgame.
run('crosseff_coexist',cmd+'\ndea_crosseff x1 x2=y1,dmu(id) eng(simplex) method(benevolent) gen(ce) matrix\n'+check+cmd+' replacevars\n'+get+'mata: check_state(snapshot)\ndea_crosseff x1 x2=y1,dmu(id) eng(simplex) method(aggressive) gen(ce) replacevars matrix\n',parity)
# Fault injection into a real recovery: same error/Break handling, old string outputs intact.
sig='real scalar _g114_recover(real matrix X,real matrix Y,real rowvector sx,real rowvector sy,real scalar j,real scalar mode,real scalar own,real scalar d,real scalar theta,real scalar tol,real rowvector v,real rowvector u,real scalar value,real rowvector diag)'
for rc in [3498,1]:
 inj=f'mata: mata drop _g114_recover()\nmata:\n{sig} {{ _error({rc});return(1); }}\nend\n'
 body='_dea_crossgame_load\n'+check+inj+'gen str12 res_agg="keepA"\ngen str12 res_ben="keepB"\ngen str12 res_nash="keepN"\nmatrix sentinel=(991,992)\nlocal before : all matrices\n'+unchanged(cmd+' replacevars',rc)+'mata: check_state(snapshot)\nlocal after : all matrices\nassert "`before\'"=="`after\'"\n'
 body+='mata: mata drop _g114_recover()\n'+cmd+' replacevars\n'+get
 run('compute_error_'+str(rc),body,parity)
# A missing or mismatched dependency fails before outputs are replaced.
for tag,kind in [('missing_recovery','missing'),('wrong_readiness','version')]:
 p=out/(tag+'_ado');shutil.copytree(pkg,p);f=p/'_dea_crossgame_recovery.ado'
 if kind=='missing':f.unlink();rc=601
 else:f.write_text(f.read_text().replace('return(11401)','return(11400)',1));rc=3498
 run(tag,'gen str12 res_agg="keepA"\ngen str12 res_ben="keepB"\ngen str12 res_nash="keepN"\n'+unchanged(cmd+' replacevars',rc),package=p)
# Internal call: actual attempt cap, diagnostic truncation and unresolved reasons.
helper=r'''
mata:
void boundaries() {
 real matrix XR,YR,X,Y,D
 real rowvector sx,sy,C,v,u,diag
 real colvector obs
 real scalar value,source,rc,k
 XR=st_data(.,tokens("x1 x2"));YR=st_data(.,tokens("y1"));sx=colmax(XR);sy=colmax(YR);X=XR:/sx;Y=YR:/sy;obs=1::rows(X)
 D=J(0,32,.);C=J(1,13,0);C[1]=256
 // Stored canonical failed game theta, passed with exact binary precision.
 rc=_g114_dispatch(X,Y,XR,YR,sx,sy,2,3,.,4,THETA,1e-9,v,u,value,1,obs,D,C,source)
 assert(rc==1 & source==0 & C[4]==1 & rows(D)==0)
 _g114_report(D,C,"cap",0)
 savebin("cap_counts.bin",st_matrix("cap_recovery_counts"))
 diag=J(1,32,.)
 rc=_g114_recover(J(33,1,1),J(33,1,1),1,1,1,1,1,.,.,1e-9,v,u,value,diag);assert(rc==1 & diag[16]==1)
 rc=_g114_recover(J(3,9,1),J(3,8,1),J(1,9,1),J(1,8,1),1,1,1,.,.,1e-9,v,u,value,diag);assert(rc==1 & diag[16]==1)
 rc=_g114_recover(J(3,1,1e300),J(3,1,1),1e300,1,1,1,1,.,.,1e-9,v,u,value,diag);assert(rc==1 & diag[16]==2)
 rc=_g114_recover(J(3,1,1),J(3,1,1),1,1,1,3,.,2,1.01,1e-9,v,u,value,diag);assert(rc==1 & diag[16]==3 & diag[15]==0)
 savebin("infeasible_diag.bin",diag)
 D=J(50,32,.);D[.,15]=J(50,1,0);C=J(1,13,0);C[1]=C[3]=C[7]=50
 _g114_report(D,C,"trunc",0)
 savebin("trunc_counts.bin",st_matrix("trunc_recovery_counts"))
}
boundaries()
end
'''
theta=expected('game','recovery_diagnostics')[0,13];helper=helper.replace('THETA',format(theta,'.17g'))
for size in [10,40]:
 run('limits_matsize'+str(size),'_dea_crossgame_load\nset matsize '+str(size)+'\n'+helper,lambda d,size=size:np.testing.assert_array_equal(readbin(d/'trunc_counts.bin')[0,13:],[0,50] if size==10 else [40,10]))

run('core_coexist',cmd+'\ndea x1 x2=y1, dmu(id) engine(simplex) rts(crs) ort(in) stage(2) gen(te)\n'+check+cmd+' replacevars\n'+get+'mata: check_state(snapshot)\ndea x1 x2=y1, dmu(id) engine(simplex) rts(crs) ort(in) stage(1) gen(te1)\n',parity)
run('shifted_observations','expand 2\ntostring id,replace\nreplace id="" in 1/6\n'+cmd.replace(' x1 x2=y1,',' x1 x2=y1 in 7/12,')+'\nassert r(n)==6 & r(n_recovered_game)==1\nassert missing(res_nash) in 1/6\n'+get+'assert C[1,2]==2 & C[1,3]==8 & C[1,5]==4 & C[1,6]==10\n')
run('matsize10_real','set matsize 10\n'+cmd.replace(' matrix detail','')+'\nassert r(n_recovered_game)==1 & r(n_recovery_rows)==0 & r(n_recovery_rows_omitted)==1\ncapture matrix list r(recovery_diagnostics)\nassert _rc==111\n')
broken=out/'missing_after_warm_ado';shutil.copytree(pkg,broken);(broken/'_dea_crossgame_recovery.ado').unlink()
run('load_error_state',cmd+'\n'+check+f'adopath - "{pkg}"\nadopath ++ "{broken}"\nmata: mata drop _g114_recover()\n'+unchanged(cmd+' replacevars',601)+'mata: check_state(snapshot)\n')
fault=out/'discard_ado';shutil.copytree(pkg,fault);f=fault/'_dea_crossgame_recovery.ado';t=f.read_text()
anchor=' rc=_dea_cg_lp(X,Y,j,mode,own,d,theta,tol,v,u,value,nd)';assert t.count(anchor)==1
# Inject an extra failure only in the second iteration, after the real recovery.
t=t.replace(anchor,anchor+'\n if(mode==3 & iteration==2 & j==2 & d==5) {rc=1;nd=(0,99,.,1);}')
anchor=' n=rows(X);m=cols(X);s=cols(Y);primary_value=.;rc=1';assert t.count(anchor)==1
t=t.replace(anchor,anchor+'\n if(mode==3 & j==2 & d==5) {diag[15]=0;diag[16]=3;return(1);}');f.write_text(t)
body=cmd+'\nassert r(nash_status)==3 & r(nash_iterations)==1 & r(nash_failed_pairs)==1\nassert r(n_recovered_game)==1 & r(n_recovery_unresolved)==1\n'+get
for name in MATS:body+=f'mata: savebin("{name}.bin",st_matrix("r({name})"))\n'
def discarded(d):
 C=readbin(d/'C.bin');assert C.shape==(2,32) and C[0,14]==1 and C[1,14]==0 and np.all(C[:,31]==0) and np.all(C[:,0]==2)
 for name in MATS:np.testing.assert_array_equal(readbin(d/(name+'.bin')),expected('discarded_game',name))
run('discard_complete_iteration',body,discarded,package=fault)
near=out/'near.csv';near.write_text(csvdata('near'))
for reason in [4,5,6,7]:
 bad=out/('reason'+str(reason)+'_ado');shutil.copytree(pkg,bad)
 f=bad/'_dea_crossgame_recovery.ado';t=f.read_text();anchor='  diag[12]=primary_value;diag[13]=primary_value-own';assert t.count(anchor)==1
 t=t.replace(anchor,anchor+f'\n  diag[16]={reason};return(1);');f.write_text(t)
 c='dea_crossgame x1=y1,dmu(id) gen(res) strategy(all) tol(1e-9) maxiter(25) matrix detail'
 body=f'import delimited "{near}",clear asdouble\n'+c+'\nassert r(n_recovered)==0 & r(n_recovery_unresolved)==2 & r(n_failed_agg)==1 & r(n_failed_ben)==1\n'+get+f'assert C[1,16]=={reason} & C[2,16]=={reason}\nassert C[1,17]==1 & C[2,17]==1\n'
 body+='matrix A=r(source_agg)\nmatrix B=r(source_ben)\nassert A[1,1]==3 & B[1,1]==3\n'
 run('unresolved_reason'+str(reason),body,package=bad)
summary=dict(passed=all(r['passed'] for r in rows),cases=len(rows),passed_cases=sum(r['passed'] for r in rows),failed_cases=sum(not r['passed'] for r in rows),simulated_break=True,actual_user_break_not_tested=True)
dump('summary.json',summary);print(json.dumps(summary),flush=True)
raise SystemExit(0 if summary['passed'] else 1)
