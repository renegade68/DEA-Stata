#!/usr/bin/env python3
"""Research objective-unit flow, inherited resource limits and reload checks."""
from pathlib import Path
import argparse,json,subprocess
import numpy as np
from test_core_retries import WRITER,readbin
from test_crosseff_recovery import SAVE,check,integer
from test_crosseff_costunits import LOAD
REPO=Path(__file__).resolve().parents[1]
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);o=ap.parse_args();out=o.output.resolve();out.mkdir(parents=True,exist_ok=False)
 # LP families: max C*x, x in [0,10]; x>=0 unbounded; inconsistent bounds.
 cases=[]
 for family in ['positive','negative','zero','unbounded','infeasible']:
  for seed in ['lower','empty']:
   for used in [0,2,15]:
    for spent in [0,12345]:cases.append(dict(family=family,seed=seed,used=used,spent=spent))
 for spent in [79999999,80000000,-1,1.5]:cases.append(dict(family='positive',seed='lower',used=0,spent=spent))
 for used in [16,17,-1,.5]:cases.append(dict(family='positive',seed='lower',used=used,spent=0))
 code='version 16\nclear all\nset more off\nset processors 1\n'
 for n in ['PLUS','PERSONAL','SITE','OLDPLACE']:(out/n).mkdir();code+=f'sysdir set {n} "{out/n}"\n'
 helper=REPO/'ado/_dea_crosseff_costunits.ado'
 code+=f'adopath ++ "{REPO/"ado"}"\ninput id x y\n1 1 1\n2 2 1\nend\nquietly dea_crosseff x=y,eng(simplex) replacevars\nrun "{helper}"\n'+WRITER+'mata:\n'+SAVE
 code+='''
 void flow(string scalar path,real scalar family,real scalar empty,real scalar spent,real scalar used) {
  struct _C107_Result scalar result
  struct _C58_Work scalar W
  transmorphic scalar C
  real matrix G
  real colvector h,ids
  real rowvector a
  real scalar signum
  W.ops=0;W.limit=10000000;W.limbs=128;W.maxlimbs=0;W.error=0
  C=asarray_create("real",2);signum=(family==2 ? -1 : (family==3 ? 0 : 1))
  a=_c58f_restore(_c58b_int(signum*3),1531,W);asarray(C,(1,1),a)
  G=(1\\-1);h=(10\\0);ids=(empty ? J(0,1,.) : 2)
  if(family==4) {G=-1;h=0;ids=(empty ? J(0,1,.) : 1);}
  if(family==5) h=(-1\\0)
  // Opposite proposal sign deliberately tests exact-cost provenance.
  result=_c107_attempt(-signum,G,h,0,ids,C,-1500,spent,used,0)
  savebin(path+"_context.bin",(result.attempted,result.complete,result.ops,result.maxlimbs,result.original_value,result.cost.valid,result.cost.error,result.prior_ops,result.prior_pivots,result.new_pivots,result.total_pivots))
  savebin(path+"_factor.bin",result.cost.factor)
  _c63_record(path,-signum,G,h,0,result.fit)
  if(result.cost.valid) _c70_cost_record(path,result.cost.working,-1500,1)
  assert(all(asarray(C,(1,1)):==a))
 }
 '''
 code+='\n void invalid_models() {\n  struct _C107_Result scalar result\n  transmorphic scalar C\n  real matrix G,h,c\n  real scalar i,ne,wide\n  C=asarray_create("real",2);asarray(C,(1,1),(1,3))\n  for(i=1;i<=10;i++) {\n   G=(1\\-1);h=(10\\0);c=1;ne=0;wide=0\n   if(i==1) G[1,1]=.\n   if(i==2) h[1]=.\n   if(i==3) c=.\n   if(i==4) ne=.\n   if(i==5) ne=.5\n   if(i==6) ne=3\n   if(i==7) wide=2\n   if(i==8) G=J(2,2,1)\n   if(i==9) h=J(3,1,0)\n   if(i==10) G=J(257,1,1)\n   result=_c107_attempt(c,G,h,ne,2,C,0,0,0,wide)\n   assert(result.attempted==0 & result.complete==0 & result.ops==0 & result.cost.valid==0)\n  }\n }\n invalid_models()\n'
 for i,z in enumerate(cases):code+=f'flow("{out/str(i)}",{["positive","negative","zero","unbounded","infeasible"].index(z["family"])+1},{int(z["seed"]=="empty")},{z["spent"]},{z["used"]})\n'
 code+=f'end\nrun "{helper}"\nmata: assert(_c107_ready()==10702)\nmata: flow("{out/"reload"}",1,0,12345,2)\nquietly dea_crosseff x=y,eng(simplex) replacevars\ndi "T107_FLOW_PASS"\nexit,clear\n'
 (out/'test.do').write_text(code);p=subprocess.run([o.stata,'-b','do','test.do'],cwd=out,capture_output=True,timeout=1800)
 log=(out/'test.log').read_text(errors='replace') if (out/'test.log').exists() else '';assert p.returncode==0 and 'T107_FLOW_PASS' in log.splitlines(),log[-2500:]
 records=[]
 for i,z in enumerate(cases):
  a=readbin(out/(str(i)+'_context.bin'))[0];meta=readbin(out/(str(i)+'_meta.bin'))[0]
  assert a[7]==z['spent'] and a[8]==z['used']
  invalid=z['spent']<0 or z['spent']>=80000000 or z['spent']!=int(z['spent']) or z['used']<0 or z['used']>=16 or z['used']!=int(z['used'])
  if invalid:assert a[0]==a[1]==0 and a[2]==z['spent']
  elif z['spent']==79999999:assert a[0]==1 and a[1]==a[5]==0 and a[6]==1
  else:
   assert a[2]<=80000000 and a[10]<=16 and a[10]==z['used']+a[9]
   proof,ray=check(out/str(i))
   if z['family']=='infeasible':assert not proof and not ray and a[1]==0
   elif z['used']==15 and z['seed']=='empty' and z['family']=='positive':assert a[1] in [0,1] # auxiliary pivots may exhaust the remaining single pivot
   elif z['family']=='unbounded':assert ray==1 and proof==0 and a[1]==1
   else:
    assert proof==1 and a[1]==1
    x=integer(readbin(out/(str(i)+'_num.bin'))[0])/integer(readbin(out/(str(i)+'_den.bin'))[0]);assert x==(10 if z['family']=='positive' else 0)
    # _c70_value reports the minimizing objective, -C*x, in original units.
    assert a[4]==(-3*2**31*10 if z['family']=='positive' else 0)
  records.append(dict(**z,attempted=int(a[0]),complete=int(a[1]),ops=a[2],new_pivots=a[9],total_pivots=a[10],status=meta[0]))
 # Reload the new helper with a compiled caller still resident.
 same=next(str(i) for i,z in enumerate(cases) if z==dict(family='positive',seed='lower',used=2,spent=12345))
 for f in out.glob('reload_*.bin'):assert f.read_bytes()==(out/(same+f.name[len('reload'):])).read_bytes()
 summary=dict(passed=True,cases=len(cases),invalid_models=10,reload_repeated_case=True,existing_command_after_reload=True,certificate_unit_proofs=True,records=records)
 (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('FLOW_PASS',len(cases),'cases plus reload')
if __name__=='__main__':main()
