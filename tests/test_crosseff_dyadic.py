#!/usr/bin/env python3
"""Portable exact dyadic arithmetic: signed quotients, residuals and fixed resource limits."""
from pathlib import Path
import argparse,json,subprocess,time
import random,struct
import numpy as np
from fractions import Fraction as F
from test_core_retries import WRITER,readbin
from test_crosseff_recovery import integer
REPO=Path(__file__).resolve().parents[1]

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata',required=True)
    ap.add_argument('--output',required=True,type=Path)
    options=ap.parse_args()
    out=options.output.resolve();out.mkdir(parents=True,exist_ok=False);rng=random.Random(20261111)
    def pack(x):
     s=(x>0)-(x<0);x=abs(x);a=[s]
     while x:a.append(x%(2**26));x//=2**26
     assert len(a)<=129
     return a+[0]*(129-len(a))
    def write(n,A):
     A=np.array(A,dtype='<f8');(out/n).write_bytes(struct.pack('<ii',*A.shape)+A.tobytes())
    cases=[]
    def add(vals,tag='random',error=None,limit=10000000,limbs=128,prior=0):
     a,b,c,d,den=vals;want=F(a*b-c*d,den) if den else None
     if error is None:error=3 if want is None or want.denominator!=1 else (2 if abs(want.numerator).bit_length()>3328 else 0)
     cases.append(dict(vals=vals,tag=tag,error=error,want=want,limit=limit,limbs=limbs,prior=prior))
    for i in range(2048):
     factor=rng.getrandbits(31)|1;ea,eb,ec,ed=[rng.randrange(1700,1901) for _ in range(4)]
     a=(rng.getrandbits(64)|1)*factor<<ea;b=(rng.getrandbits(64)|1)<<eb;c=(rng.getrandbits(64)|1)*factor<<ec;d=(rng.getrandbits(64)|1)<<ed
     if i%2:a=-a
     if i%3:c=-c
     if i%5:b=-b
     if i%7:d=-d
     if i%17==0:a=0
     if i%19==0:c=0
     if i%23==0:c=a;d=b
     den=factor<<rng.randrange(0,min(ea+eb,ec+ed,3200)+1)
     if i%2:den=-den
     add([a,b,c,d,den])
    edge=[([1<<1663,1<<1664,0,0,1],'last128'),([1<<1664,1<<1664,0,0,1],'result129'),([3<<1500,1<<1500,1<<1500,1<<1500,1<<3001],'cancellation_exponent'),([3<<1500,1<<1500,0,0,1<<3001],'fractional_power2'),([5<<1800,7<<1800,0,0,3<<1000],'fractional_odd'),([0,0,0,0,1<<3327],'zero'),([0,0,0,0,0],'zero_den'),([3,7,2,6,-3],'negative_den')]
    for v,t in edge:add(v,t)
    add([(1<<2500)-1]*4+[1],'odd_intermediate_over_limit',2)
    for setting,value,err in [('limit',0,1),('limbs',1,2),('prior',3,3),('prior',1,1)]:
     args={setting:value};add([3<<1700,5<<1700,1<<1700,1<<1700,1<<1000],setting,err,**args)
    for i in range(40):add([i+3,21,i+2,6,3],'ordinary')
    flat=[pack(v) for z in cases for v in z['vals']];write('quotient_inputs.bin',flat);write('quotient_meta.bin',[[z['limit'],z['limbs'],z['prior']] for z in cases])
    residuals=[]
    for i in range(512):
     p=rng.randrange(1300,1600);q=rng.randrange(1300,1600);a=(rng.getrandbits(120)|1)<<p;b=(rng.getrandbits(120)|1)<<q
     c=(rng.getrandbits(32)|1)<<p;d=(rng.getrandbits(32)|1)<<q
     rhs=c if i%3==0 else (c-2**p if i%3==1 else c+2**p)
     residuals.append(dict(a=[a,-a,c],n=[b,b,d],rhs=rhs,den=d,error=0,limit=10000000,limbs=128,prior=0,tag='signed_cancellation'))
    for setting,value,err in [('limit',0,1),('limbs',1,2),('prior',3,3)]:
     z=dict(a=[1<<1700],n=[1<<1700],rhs=0,den=1,error=err,limit=10000000,limbs=128,prior=0,tag=setting);z[setting]=value;residuals.append(z)
    flat=[];meta=[]
    for z in residuals:
     start=len(flat)+1;flat.extend(pack(v) for v in z['a']+z['n']+[z['rhs'],z['den']]);meta.append([start,len(z['a']),z['limit'],z['limbs'],z['prior']])
    write('residual_inputs.bin',flat);write('residual_meta.bin',meta)
    for name in ['PLUS','PERSONAL','SITE','OLDPLACE']:(out/name).mkdir()
    code='version 16\nclear all\nset more off\nset processors 1\n'
    for name in ['PLUS','PERSONAL','SITE','OLDPLACE']:code+=f'sysdir set {name} "{out/name}"\n'
    code+=f'adopath ++ "{REPO/"ado"}"\ninput id x y\n1 1 1\n2 2 1\nend\nquietly dea_crosseff x=y,eng(simplex)\n'+WRITER
    code+='''mata:
    real matrix loadbin(string scalar file) {
     real matrix B,A
     real scalar f,n,m
     B=bufio();bufbyteorder(B,2);f=fopen(file,"r");n=fbufget(B,f,"%4b",1,1);m=fbufget(B,f,"%4b",1,1)
     A=fbufget(B,f,"%8z",n,m);fclose(f);return(A)
    }
    void run_tests() {
     struct _C58_Work scalar W,W0
     real matrix M,K,Q,D,old
     real rowvector a,b,c,d,den,v,rhs,v0
     real scalar i,j,k,mode,wide,s,n,signum,untouched
     transmorphic scalar A,N
     M=loadbin("quotient_inputs.bin");K=loadbin("quotient_meta.bin");Q=J(rows(K)*2,136,0)
     for(i=1;i<=rows(K);i++) for(mode=1;mode<=2;mode++) {
      s=(i-1)*5;a=_c58b_norm(M[s+1,.]);b=_c58b_norm(M[s+2,.]);c=_c58b_norm(M[s+3,.]);d=_c58b_norm(M[s+4,.]);den=_c58b_norm(M[s+5,.]);old=M[(s+1)..(s+5),.]
      W.ops=0;W.limit=K[i,1];W.limbs=K[i,2];W.maxlimbs=0;W.error=K[i,3];W0=W
      wide=_c58f_wide_product(a,b,W)|_c58f_wide_product(c,d,W)
      if(mode==1) v=_c58f_quotient(a,b,c,d,den,W)
      else v=_c58f_update_quotient(a,b,c,d,den,W)
      untouched=all(a:==_c58b_norm(old[1,.])) & all(b:==_c58b_norm(old[2,.])) & all(c:==_c58b_norm(old[3,.])) & all(d:==_c58b_norm(old[4,.])) & all(den:==_c58b_norm(old[5,.]))
      k=(i-1)*2+mode;Q[k,1..7]=(W.error,W.maxlimbs,W.ops,wide,untouched,W.limit,W.limbs);Q[k,8..(7+cols(v))]=v
      if(mode==2 & !wide) {
       a=_c58b_mul(a,b,W0);b=_c58b_mul(c,d,W0);v0=_c58b_div(_c58b_add(a,_c58b_neg(b),W0),den,W0)
       assert(all(v0:==v) & W.ops==W0.ops & W.error==W0.error & W.maxlimbs==W0.maxlimbs)
      }
     }
     savebin("quotient_results.bin",Q)
     M=loadbin("residual_inputs.bin");K=loadbin("residual_meta.bin");Q=J(rows(K),7,0)
     for(i=1;i<=rows(K);i++) {
      s=K[i,1];n=K[i,2];A=asarray_create("real",2);N=asarray_create("real",2)
      for(j=1;j<=n;j++) {asarray(A,(1,j),_c58b_norm(M[s+j-1,.]));asarray(N,(j,1),_c58b_norm(M[s+n+j-1,.]));}
      rhs=_c58b_norm(M[s+2*n,.]);den=_c58b_norm(M[s+2*n+1,.])
      W.ops=0;W.limit=K[i,3];W.limbs=K[i,4];W.maxlimbs=0;W.error=K[i,5]
      signum=_c58f_residual_sign(A,1,N,rhs,den,n,W);untouched=1
      for(j=1;j<=n;j++) untouched=untouched & all(asarray(A,(1,j)):==_c58b_norm(M[s+j-1,.])) & all(asarray(N,(j,1)):==_c58b_norm(M[s+n+j-1,.]))
      untouched=untouched & all(rhs:==_c58b_norm(M[s+2*n,.])) & all(den:==_c58b_norm(M[s+2*n+1,.]))
      Q[i,.]=(signum,W.error,W.maxlimbs,W.ops,untouched,W.limit,W.limbs)
     }
     savebin("residual_results.bin",Q)
    }
    run_tests()
    end
    di "CROSSEFF_DYADIC_PASS"
    exit,clear
    '''
    (out/'test.do').write_text(code)
    p=subprocess.run([options.stata,'-b','do','test.do'],cwd=out,capture_output=True,timeout=180)
    log=(out/'test.log').read_text(errors='replace');assert p.returncode==0 and 'CROSSEFF_DYADIC_PASS' in log.splitlines(),log[-3000:]
    Q=readbin(out/'quotient_results.bin')
    for i,z in enumerate(cases):
     for row in Q[2*i:2*i+2]:
      assert row[0]==z['error'],(i,z['tag'],row[:7],z['error']);assert row[4]==1 and row[5]==z['limit'] and row[6]==z['limbs']
      if not z['error']:assert integer(row[7:])==z['want'].numerator and row[1]<=128
    R=readbin(out/'residual_results.bin')
    for z,row in zip(residuals,R):
     want=sum(a*n for a,n in zip(z['a'],z['n']))-z['rhs']*z['den'];assert row[1]==z['error'] and row[4]==1 and row[5]==z['limit'] and row[6]==z['limbs']
     if not z['error']:assert row[0]==(want>0)-(want<0) and row[2]<=128
    summary=dict(passed=True,quotient_cases=len(cases),quotient_interfaces=2,nonzero_expected_quotients=sum(z["want"] is not None and z["want"]!=0 for z in cases),zero_expected_quotients=sum(z["want"]==0 for z in cases),quotient_evaluations=len(Q),wide_quotient_cases=int(sum(Q[::2,3])),quotient_errors_refused=sum(bool(z['error']) for z in cases),residual_cases=len(residuals),residual_error_controls=3,exact_fraction_identities=True,ordinary_results_and_costs_preserved=True,inputs_and_limits_unchanged=True,Stata_batches=1)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary),flush=True)

if __name__=='__main__':main()
