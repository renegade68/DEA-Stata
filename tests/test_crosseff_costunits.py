#!/usr/bin/env python3
"""Research checks for bounded exact positive objective factors (fresh Stata)."""
from pathlib import Path
import argparse,json,math,random,struct,subprocess
import numpy as np
from test_core_retries import WRITER,readbin
from test_crosseff_recovery import integer
REPO=Path(__file__).resolve().parents[1]
LOAD='''
real matrix loadbin(string scalar file) {
 real matrix B,A
 real scalar f,n,m
 B=bufio();bufbyteorder(B,2);f=fopen(file,"r");n=fbufget(B,f,"%4b",1,1);m=fbufget(B,f,"%4b",1,1)
 A=fbufget(B,f,"%8z",n,m);fclose(f);return(A)
}
'''
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);o=ap.parse_args();out=o.output.resolve();out.mkdir(parents=True,exist_ok=False)
 def pack(x):
  s=(x>0)-(x<0);a=[s];x=abs(x)
  while x:a.append(x%(2**26));x//=2**26
  assert len(a)<=130
  return a+[0]*(130-len(a))
 def write(n,A):
  a=np.asarray(A,dtype='<f8');(out/n).write_bytes(struct.pack('<ii',*a.shape)+a.tobytes())
 rng=random.Random(10720261108);pairs=[]
 def pair(a,b,tag='random',prior=0,limbs=128,limit=10000000,spent=0):pairs.append(dict(a=a,b=b,tag=tag,prior=prior,limbs=limbs,limit=limit,spent=spent))
 for i in range(256):
  g=rng.getrandbits(rng.choice([1,26,100,510,1531,3000]))|1
  a=g*(rng.getrandbits(rng.randrange(1,200))|1);b=g*(rng.getrandbits(rng.randrange(1,200))|1)
  if a.bit_length()>3328 or b.bit_length()>3328:continue
  if i%2:a=-a
  if i%3:b=-b
  pair(a,b,'common_factor');pair(a*b//g if (a*b//g).bit_length()<=3328 else a,b,'division')
 for a,b in [(0,0),(0,7),(9,0),(1,1),(-1,1),(1,-1),(-9,-3),(3,7),(7,3),(2**3327,2**3326),(2**3328-1,1),(2**3328-1,2**3328-1),(2**3328-2,(2**3328-2)//2),(2**3328,1),(1,2**3328)]:pair(a,b,'boundary')
 for p in range(1,6):pair(3,7,'prior',prior=p,spent=31)
 for a,b in [(3,7),(1<<3000,1<<2999),(2**3328-1,1)]:
  pair(a,b,'zero_budget',limit=0);pair(a,b,'spent_budget',limit=10,spent=10);pair(a,b,'small_cap',limbs=1)
 write('pairs.bin',[pack(v) for z in pairs for v in [z['a'],z['b']]]);write('pair_meta.bin',[[z[k] for k in ['prior','limbs','limit','spent']] for z in pairs])
 costs=[]
 def cost(v,tag='random',prior=0,limbs=128,limit=10000000,spent=0,mutation=0,n=None,ce=-1020):costs.append(dict(v=v,tag=tag,prior=prior,limbs=limbs,limit=limit,spent=spent,mutation=mutation,n=len(v) if n is None else n,ce=ce))
 for i in range(256):
  g=(rng.getrandbits(rng.choice([1,52,511,1531,3000]))|1)<<rng.randrange(0,31)
  v=[g*rng.randrange(-10,11) for _ in range(rng.randrange(1,9))];cost(v)
 for v in [[0]*6,[0],[1,-1,0],[6,10,15],[2**3328-1],[-(2**3328-1),0],[2**3328-2,(2**3328-2)//2],[2**3328],[9,10],[2**3000,3*2**3000]]:cost(v,'boundary')
 for p in range(1,6):cost([3,7], 'prior',prior=p,spent=19)
 for change in range(1,7):cost([15,-30,0], 'malformed',mutation=change)
 for n in [0,65,1.5]:cost([1], 'invalid_n',n=n)
 cost([1],'invalid_ce',ce=.5)
 cost([3,6],'zero_budget',limit=0)
 cost([3,6],'spent_budget',limit=12,spent=12)
 cost([1<<40,1<<39],'small_cap',limbs=1)
 starts=[];flat=[]
 for z in costs:starts.append(len(flat)+1);flat.extend(pack(v) for v in z['v'])
 write('costs.bin',flat);write('cost_meta.bin',[[s,len(z['v'])]+[z[k] for k in ['prior','limbs','limit','spent','mutation','n','ce']] for s,z in zip(starts,costs)])
 code='version 16\nclear all\nset more off\nset processors 1\n'
 for n in ['PLUS','PERSONAL','SITE','OLDPLACE']:(out/n).mkdir();code+=f'sysdir set {n} "{out/n}"\n'
 code+=f'adopath ++ "{REPO/"ado"}"\ninput id x y\n1 1 1\n2 2 1\nend\nquietly dea_crosseff x=y,eng(simplex)\nrun "{REPO/"ado/_dea_crosseff_costunits.ado"}"\n'+WRITER+'mata:\n'+LOAD
 code+='''
 void arithmetic_tests() {
  struct _C58_Work scalar W,W2
  struct _C107_Cost scalar unit
  real matrix M,K,Q,V,D
  real rowvector a,b,v,copya,copyb
  real scalar i,j,k,s,n,mode,ops,exact,shortage,unchanged
  transmorphic scalar C
  M=loadbin("pairs.bin");K=loadbin("pair_meta.bin");Q=J(2*rows(K),137,0)
  for(i=1;i<=rows(K);i++) for(mode=1;mode<=2;mode++) {
   a=_c58b_norm(M[2*i-1,.]);b=_c58b_norm(M[2*i,.]);copya=a;copyb=b
   W.ops=K[i,4];W.limit=K[i,3];W.limbs=K[i,2];W.maxlimbs=0;W.error=K[i,1]
   v=(mode==1 ? _c107_gcd(a,b,W) : _c107_div(a,b,W));k=2*i-2+mode
   Q[k,1..7]=(W.error,W.ops,W.maxlimbs,all(a:==copya)&all(b:==copyb),W.limit,W.limbs,mode);Q[k,8..(7+cols(v))]=v
  }
  savebin("pair_results.bin",Q)
  M=loadbin("costs.bin");K=loadbin("cost_meta.bin");Q=J(rows(K),139,0);D=J(rows(M),129,0)
  for(i=1;i<=rows(K);i++) {
   s=K[i,1];n=K[i,2];C=asarray_create("real",2)
   for(j=1;j<=n;j++) asarray(C,(1,j),_c58b_norm(M[s+j-1,.]))
   a=asarray(C,(1,1))
   if(K[i,7]==1) a[2]=a[2]+.5
   if(K[i,7]==2) a[1]=2
   if(K[i,7]==3) a[2]=67108864
   if(K[i,7]==4) a=(a,0)
   if(K[i,7]==5) a[2]=.
   asarray(C,(1,1),a)
   if(K[i,7]==6) asarray_remove(C,(1,1))
   W.ops=K[i,6];W.limit=K[i,5];W.limbs=K[i,4];W.maxlimbs=0;W.error=K[i,3]
   unit=_c107_cost(C,K[i,8],K[i,9],W);ops=W.ops;exact=shortage=0
   if(unit.valid) {
    for(j=1;j<=n;j++) {
     v=asarray(unit.working,(1,j));D[s+j-1,1..cols(v)]=v
     assert(all(asarray(unit.original,(1,j)):==_c58b_norm(M[s+j-1,.])))
    }
    // Repeat at the measured exact remaining budget, and one operation less.
    W2.ops=K[i,6];W2.limit=ops;W2.limbs=K[i,4];W2.maxlimbs=0;W2.error=0
    unit=_c107_cost(C,n,K[i,9],W2);exact=unit.valid & W2.ops==ops
    W2.ops=K[i,6];W2.limit=ops-1;W2.maxlimbs=0;W2.error=0
    unit=_c107_cost(C,n,K[i,9],W2);shortage=!unit.valid & W2.error==1
    // Return the successful context, then prove its copies do not alias C.
    W2.ops=K[i,6];W2.limit=K[i,5];W2.maxlimbs=0;W2.error=0
    unit=_c107_cost(C,n,K[i,9],W2)
    asarray(unit.working,(1,1),(1,999));asarray(unit.original,(1,1),(1,998))
   }
   unchanged=1
   if(K[i,7]==0) for(j=1;j<=n;j++) unchanged=unchanged & all(asarray(C,(1,j)):==_c58b_norm(M[s+j-1,.]))
   Q[i,1..10]=(unit.valid,W.error,W.ops,W.maxlimbs,exact,shortage,unchanged,W.limit,W.limbs,asarray_elements(unit.working))
   Q[i,11..(10+cols(unit.factor))]=unit.factor
  }
  savebin("cost_results.bin",Q);savebin("working.bin",D)
 }
 arithmetic_tests()
 end
 di "T107_ARITHMETIC_PASS"
 exit,clear
 '''
 (out/'test.do').write_text(code);p=subprocess.run([o.stata,'-b','do','test.do'],cwd=out,capture_output=True,timeout=1800)
 log=(out/'test.log').read_text(errors='replace') if (out/'test.log').exists() else ''
 assert p.returncode==0 and 'T107_ARITHMETIC_PASS' in log.splitlines(),log[-2500:]
 Q=readbin(out/'pair_results.bin');errors=0
 for i,z in enumerate(pairs):
  for mode in [1,2]:
   row=Q[2*i+mode-1];err=int(row[0]);a=z['a'];b=z['b'];value=integer(row[7:]);assert row[3]==1 and row[4]==z['limit'] and row[5]==z['limbs']
   want=math.gcd(a,b) if mode==1 else (a//b if b and a%b==0 else None)
   if z['prior']:assert err==z['prior'] and row[1]==z['spent']
   elif max(abs(a).bit_length(),abs(b).bit_length())>26*z['limbs']:assert err==2
   elif z['limit']<=z['spent']:assert err==1
   elif want is None:assert err==3
   else:assert err==0 and value==want,(i,mode,z,err,value,want)
   if err:errors+=1
 Q=readbin(out/'cost_results.bin');D=readbin(out/'working.bin');success=0
 for i,z in enumerate(costs):
  row=Q[i];err=int(row[1]);g=integer(row[10:]);assert row[6]==1 and row[7]==z['limit'] and row[8]==z['limbs']
  expected=z['prior'] or (4 if z['mutation'] or z['n']<1 or z['n']>64 or z['n']!=int(z['n']) or z['ce']!=int(z['ce']) else 2 if max(map(lambda x:abs(x).bit_length(),z['v']))>26*z['limbs'] else 1 if z['limit']<=z['spent'] else 0)
  assert err==expected,(i,z,row[:10])
  if err:assert row[0]==0 and row[9]==0
  else:
   want=math.gcd(*z['v']) or 1;assert g==want and row[0]==1 and row[4]==row[5]==1
   v=[integer(a) for a in D[starts[i]-1:starts[i]-1+len(z['v'])]];assert [g*x for x in v]==z['v'];success+=1
 summary=dict(passed=True,pairs=len(pairs),pair_evaluations=2*len(pairs),pair_expected_refusals=errors,cost_cases=len(costs),cost_successes=success,cost_budget_threshold_repeats=2*success,independent_integer_identity=True,inputs_immutable=True)
 (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))
if __name__=='__main__':main()
