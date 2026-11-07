#!/usr/bin/env python3
"""Exact binary64 decoder boundaries and unchanged default coefficient domain."""
from pathlib import Path
from fractions import Fraction
import argparse,json,math,random,struct,subprocess,time
import numpy as np
from test_core_retries import WRITER
REPO=Path(__file__).resolve().parents[1]

def readbin(path):
    data=path.read_bytes();n,m=struct.unpack('<ii',data[:8])
    return np.frombuffer(data[8:],dtype='<f8').reshape(n,m)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata',required=True)
    ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    values={0.}
    for exponent in range(-1074,1023):
        v=math.ldexp(1.,exponent)
        for q in [v,np.nextafter(v,0),np.nextafter(v,math.inf)]:
            if 0<q<2.**1023:values.update([float(q),-float(q)])
    rng=random.Random(20261108)
    for _ in range(1024):
        v=struct.unpack('<d',struct.pack('<Q',rng.randrange(1,0x7fe0000000000000)))[0]
        values.update([v,-v])
    values=sorted(values);(out/'inputs.bin').write_bytes(np.array(values,dtype='<f8').tobytes())
    (out/'controls.bin').write_bytes(np.array([-2.**1023,-np.finfo(float).max,np.nextafter(0.,1.),np.nextafter(2.**1023,0.)],dtype='<f8').tobytes())
    code='version 16\nclear all\nset more off\nset processors 1\n'
    for name in ['PLUS','PERSONAL','SITE','OLDPLACE']:
        (out/name).mkdir();code+=f'sysdir set {name} "{out/name}"\n'
    code+=f'adopath ++ "{REPO/"ado"}"\ninput id x y\n1 1 1\n2 2 1\nend\nquietly dea_crosseff x=y,dmu(id) eng(simplex)\n'
    code+='mata: assert(_c58_cert_ready()==6202 & _c70_exact_ready()==7002 & _c84_ready()==8402 & _c57_simplex_ready()==8402)\n'+WRITER
    code+='''mata:
void decoder_check(real colvector values,real scalar mode,string scalar file) {
 struct _C58_Work scalar W
 real matrix Q
 real rowvector d
 real scalar i,e
 Q=J(rows(values),9,0)
 for(i=1;i<=rows(values);i++) {
  W.ops=0;W.limit=10000000;W.limbs=128;W.maxlimbs=0;W.error=0
  if(mode==-1) d=_c58b_decode(values[i],e,W)
  else d=_c58b_decode(values[i],e,W,mode)
  Q[i,1..6]=(values[i],e,W.error,W.ops,W.maxlimbs,d[1])
  Q[i,7..(5+cols(d))]=d[2..cols(d)]
 }
 savebin(file,Q)
}
void decoder_control(real scalar x,real scalar limit,real scalar limbs,real scalar prior) {
 struct _C58_Work scalar W
 real rowvector a
 real scalar e
 W.ops=0;W.limit=limit;W.limbs=limbs;W.maxlimbs=0;W.error=prior
 a=_c58b_decode(x,e,W,1)
 printf("DECODER_CONTROL %g %g %g %g\\n",W.error,W.ops,W.maxlimbs,e)
}
'''
    code+=f'B=bufio();bufbyteorder(B,2);f=fopen("inputs.bin","r");V=fbufget(B,f,"%8z",{len(values)},1);fclose(f)\n'
    code+='decoder_check(V,-1,"default.bin");decoder_check(V,0,"narrow.bin");decoder_check(V,1,"wide.bin")\n'
    code+='f=fopen("controls.bin","r");C=fbufget(B,f,"%8z",4,1);fclose(f)\n'
    for x in ['.','.a','.z','C[1]','C[2]']:code+=f'decoder_control({x},10000000,128,0)\n'
    for x in ['C[3]','C[4]']:
        code+=f'decoder_control({x},0,128,0)\ndecoder_control({x},10000000,2,0)\ndecoder_control({x},10000000,128,2)\n'
    code+='end\ndi "CROSSEFF_DECODER_PASS"\nexit,clear\n';(out/'test.do').write_text(code)
    start=time.monotonic();p=subprocess.run([args.stata,'-b','do','test.do'],cwd=out,capture_output=True,timeout=300)
    log=(out/'test.log').read_text(errors='replace');assert p.returncode==0 and 'CROSSEFF_DECODER_PASS' in log.splitlines(),log[-3000:]
    a=readbin(out/'default.bin');b=readbin(out/'narrow.bin');c=readbin(out/'wide.bin');assert a.tobytes()==b.tobytes()
    eligible=np.array([v==0 or 2.**-900<=abs(v)<=2.**900 for v in values])
    assert np.all(a[~eligible,2]==4) and np.all(c[:,2]==0)
    assert a[eligible].tobytes()==c[eligible].tobytes()
    assert c[:,0].tobytes()==np.array(values,dtype='<f8').tobytes()
    for value,row in zip(values,c):
        assert row[5] in [-1,0,1] and all(d==int(d) and 0<=d<2**26 for d in row[6:])
        n=int(row[5])*sum(int(d)*2**(26*j) for j,d in enumerate(row[6:]))
        assert Fraction(n)*Fraction(2)**int(row[1])==Fraction(value)
        if value:assert abs(n)%2==1 and abs(n)<2**53
    controls=[[float(v) for v in line.split()[1:]] for line in log.splitlines() if line.startswith('DECODER_CONTROL ')]
    assert [int(v[0]) for v in controls]==[4]*5+[1,2,2]*2
    summary=dict(passed=True,operands=len(values),expanded_operands=int(sum(~eligible)),subnormal_operands=sum(0<abs(v)<2.**-1022 for v in values),exact_fraction_identities=True,default_domain_unchanged=True,default_and_explicit_narrow_identical=True,old_domain_costs_identical=True,invalid_budget_controls=len(controls),seconds=time.monotonic()-start)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)

if __name__=='__main__':main()
