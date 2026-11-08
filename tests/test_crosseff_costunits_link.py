#!/usr/bin/env python3
"""Bounded secondary-cost retry eligibility and cumulative journal contracts."""
from pathlib import Path
import argparse,json,subprocess
import numpy as np
from test_core_retries import WRITER,readbin
REPO=Path(__file__).resolve().parents[1]
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);o=ap.parse_args();out=o.output.resolve();out.mkdir(parents=True,exist_ok=False)
 code='version 16\nclear all\nset more off\nset processors 1\n'
 for n in ['PLUS','PERSONAL','SITE','OLDPLACE']:(out/n).mkdir();code+=f'sysdir set {n} "{out/n}"\n'
 code+=f'adopath ++ "{REPO/"ado"}"\ninput id x y\n1 1 1\n2 2 1\nend\nquietly dea_crosseff x=y,eng(simplex)\n'+WRITER
 code+='''mata:
 void links() {
  struct _C58_Cert scalar prior,changed,joined
  struct _C107_Result scalar retry
  struct _C58_Work scalar W
  transmorphic scalar C
  real rowvector a,t
  real matrix summary
  real scalar i,kind,used,spent,want,n,kept
  W.ops=0;W.limit=10000000;W.limbs=128;W.maxlimbs=0;W.error=0
  a=_c58f_restore((1,3),1531,W);C=asarray_create("real",2);asarray(C,(1,1),a)
  prior=_c70b_certify(1,(1\\-1),(10\\0),0,J(0,1,.),0,C,0)
  prior.ops=12345;prior.maxlimbs=138;prior.solver_status=4
  t=_c62_initial(prior);t[1..4]=(1,3,0,2);t[10]=1;t[11]=t[13]=2;t[17]=10000;t[20]=2345;t[22]=12345;t[23]=138;t[32..34]=(-1,12,.)
  prior.trace62=t;assert(_c108_needed(prior))
  for(i=1;i<=12;i++) {
   changed=prior;want=0
   if(i==1) changed.status=1
   if(i==2) changed.status=0
   if(i==3) changed.trace62[16]=1
   if(i==4) changed.ops=80000000
   if(i==5) changed.trace62[13]=16
   if(i==6) changed.trace62[33]=11
   if(i==7) changed.trace62[33]=4
   if(i==8) changed.trace62[33]=.
   if(i==9) {changed.reason=12;want=1;}
   if(i==10) {changed.trace62[33]=.;changed.trace62[27]=12;want=1;}
   if(i==11) {changed.trace62[33]=.;changed.trace62[27]=11;}
   if(i==12) want=1
   assert(_c108_needed(changed)==want)
  }
  summary=J(16,9,.)
  for(kind=1;kind<=4;kind++) for(used=0;used<=3;used++) {
   changed=prior;changed.trace62[11]=changed.trace62[13]=used
   if(kind==3) {changed.ops=79999999;changed.trace62[20]=79989999;changed.trace62[22]=79999999;}
   if(kind==2) retry=_c107_attempt(1,-1,0,0,J(0,1,.),C,0,changed.ops,used,0)
   else if(kind==4) retry=_c107_attempt(1,(1\\-1),(10\\0),0,2,C,0,changed.ops,used,0)
   else retry=_c107_attempt(1,(1\\-1),(10\\0),0,J(0,1,.),C,0,changed.ops,used,0)
   joined=_c108_join(changed,retry);t=joined.trace62
   assert(joined.solver_status==4 & joined.ops==retry.ops & joined.maxlimbs>=138)
   assert(t[13]==sum(t[(7,9,11,12,25)]) & t[13]<=16 & t[22]==sum(t[17..21]))
   assert(t[5..6]==changed.trace62[5..6])
   assert(t[7]==changed.trace62[7] & t[17]==changed.trace62[17])
   assert(t[3]==5 & t[22]==joined.ops & t[23]==joined.maxlimbs)
   if(kind==1 | kind==4) assert(joined.status==1 & t[2]==1 & retry.cost.valid)
   if(kind==2) assert(joined.status!=1 & t[2]==2 & t[16]==1 & retry.cost.valid)
   if(kind==3) assert(joined.status==changed.status & joined.reason==changed.reason & t[2]==3 & !retry.cost.valid)
   i=(kind-1)*4+used+1;summary[i,.]=(kind,used,joined.status,t[2],changed.ops,joined.ops,t[13],retry.new_pivots,retry.complete)
   // Mismatched inherited accounting must not be adopted.
   retry.prior_ops=retry.prior_ops+1;joined=_c108_join(changed,retry)
   assert(joined.ops==changed.ops & joined.status==changed.status & all(joined.trace62:==changed.trace62))
  }
  savebin("links.bin",summary)
 }
 links()
 end
 di "T108_LINK_PASS"
 exit,clear
 '''
 (out/'test.do').write_text(code);p=subprocess.run([o.stata,'-b','do','test.do'],cwd=out,capture_output=True,timeout=1800);log=(out/'test.log').read_text(errors='replace');assert p.returncode==0 and 'T108_LINK_PASS' in log.splitlines(),log[-2000:]
 rows=readbin(out/'links.bin');assert rows.shape==(16,9)
 np.testing.assert_array_equal(rows[:,0],np.repeat([1,2,3,4],4))
 summary=dict(passed=True,eligibility_cases=12,joined_flows=16,accounting_mismatch_refusals=16,finite_ray_budget_failure=True,prior_solver_status_preserved=True)
 (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))
if __name__=='__main__':main()
