#!/usr/bin/env python3
"""Portable failure-only raw recovery/peer regression. Fixed independent references."""
from pathlib import Path
from fractions import Fraction as F
import argparse,csv,json,struct,subprocess,time,hashlib
import numpy as np
from test_core_retries import WRITER
REPO=Path(__file__).resolve().parents[1]
MISS=2.**1023;MAX=np.nextafter(MISS,0)
def readbin(p):
 b=Path(p).read_bytes();r,c=struct.unpack('<ii',b[:8]);return np.frombuffer(b[8:],dtype='<f8').reshape(r,c)
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
 spec=json.loads((REPO/'tests/fixtures/crosseff_raw.json').read_text());tol=spec['tolerance'];tick=chr(96)
 code='version 16\nclear all\nset more off\nset processors 1\nset type double\nset linesize 240\n'
 for k in ['PLUS','PERSONAL','SITE','OLDPLACE']:(out/k).mkdir();code+=f'sysdir set {k} "{out/k}"\n'
 code+=f'adopath ++ "{REPO/"ado"}"\n'+WRITER
 for z in spec['cases']:
  d=out/z['tag'];d.mkdir()
  with (d/'data.csv').open('w') as f:
   w=csv.writer(f);w.writerow(['id','x1','x2','y1','y2']);w.writerows(z['data'])
  code+=f'import delimited using "{d/"data.csv"}",clear varnames(1) asdouble\n'
  code+=f'dea_crosseff x1 x2=y1 y2,dmu(id) eng(simplex) rts({z["rts"]}) ori({z["orientation"]}) method({z["method"]}) vrscross({z["rule"]}) tol({tol}) gen(res) matrix\nassert r(raw_recovery_version)==8401\n'
  for k in ['ce_matrix','weights_input','weights_output','intercept','stage_status','peer_count','raw_recovery_diagnostics','raw_peer_status']:
   code+=f'mata: savebin("{d/k}.bin",st_matrix("r({k})"))\n'
  code+=f'mata: savebin("{d}/data.bin",st_data(.,tokens("id x1 x2 y1 y2 res res_self res_maverick")))\n'
 for i,z in enumerate(spec['boundaries']):
  vs=[ '('+','.join(format(v,'.17g') for v in z[k])+')' for k in ['v','x','u','y']]
  args=','.join(vs)+','+','.join(format(z[k],'.17g') for k in ['w','free','lz','tol'])
  code+=f'mata: D=J(1,3,.);v=_c84_peer({args},D);savebin("{out}/boundary{i}.bin",(v,D))\n'
 code+='di "CROSSEFF_RAW_PASS"\nexit,clear\n';(out/'test.do').write_text(code)
 start=time.monotonic();p=subprocess.run([a.stata,'-b','do','test.do'],cwd=out,capture_output=True,timeout=1800)
 log=(out/'test.log').read_text(errors='replace');assert p.returncode==0 and 'CROSSEFF_RAW_PASS' in log.splitlines(),log[-4000:]
 records=[]
 for z in spec['cases']:
  d=out/z['tag'];data=readbin(d/'data.bin');W=np.c_[readbin(d/'weights_input.bin'),readbin(d/'weights_output.bin'),readbin(d/'intercept.bin')]
  np.testing.assert_array_equal(W,z['weights']);np.testing.assert_array_equal(data[:,6],z['self'])
  np.testing.assert_array_equal(readbin(d/'stage_status.bin').ravel(),z['stages'])
  np.testing.assert_array_equal(data[:,:5],z['data'])
  C=readbin(d/'ce_matrix.bin');PS=readbin(d/'raw_peer_status.bin');RD=readbin(d/'raw_recovery_diagnostics.bin')
  for i,recovered in enumerate(z['recovered']):
   assert RD[i,6]==recovered and RD[i,1]==(z['initial_stage'][i]==0)
   if not recovered:
    np.testing.assert_array_equal(C[i],z['initial_ce'][i]);assert np.all(PS[i]==0);continue
   x=[F(float(v)) for v in W[i]]
   for j,row in enumerate(z['data']):
    vx=dot(x[:2],[F(float(v)) for v in row[1:3]]);uy=dot(x[2:4],[F(float(v)) for v in row[3:]])
    den=vx;cut=F(tol) if z['rts']=='vrs' else F(0);state=2;val=MISS
    if z['rule']=='limzhu':den=vx-x[-1];cut=F(tol)*max(abs(vx),abs(x[-1]),uy)
    if den>0 and den>cut:
     q=uy/den
     if z['rule']=='limzhu' and q>F(1+tol):state=6
     elif q>F(float(MAX)) or (q and float(q)==0):state=3
     else:state=1;val=float(q)
    assert PS[i,j]==state and C[i,j]==val
  counts=[]
  for j in range(len(C)):
   total=0.;n=0
   for i in range(len(C)):
    if i!=j and C[i,j]<MISS:total+=float(C[i,j]);n+=1
   counts.append(n);value=total/n if n and total<MISS else MISS
   assert data[j,5]==value
  np.testing.assert_array_equal(readbin(d/'peer_count.bin').ravel(),counts)
  records.append(dict(test=z['tag'],passed=True,evaluators=len(C),recovered=sum(z['recovered'])))
 for i,z in enumerate(spec['boundaries']):
  got=readbin(out/f'boundary{i}.bin').ravel();want=z['expected']
  assert got[1]==want['state'] and got[0]==(MISS if want['value'] is None else want['value'])
 files={str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in out.rglob('*.bin')}
 (out/'numeric_manifest.json').write_text(json.dumps(files,indent=2)+'\n')
 (out/'results.json').write_text(json.dumps(records,indent=2)+'\n')
 summary=dict(passed=True,cases=len(records),peer_boundaries=len(spec['boundaries']),evaluators=sum(z['evaluators'] for z in records),
  recovered_records=sum(z['recovered'] for z in records),accepted_weights_preserved=True,exact_peer_verification=True,
  scope='8fixed datasets x3methods x5RTS-orientation-peer-rule combinations at1e-10;not independent datasets',
  seconds=time.monotonic()-start)
 (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
