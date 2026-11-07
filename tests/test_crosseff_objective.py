#!/usr/bin/env python3
"""Verify intended peer sums and bounded exact witnesses on private copies."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,math,shutil,subprocess,time
import numpy as np
from test_core_retries import WRITER,readbin
from test_crosseff_recovery import integer
REPO=Path(__file__).resolve().parents[1]
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def fl(x):
 n=x.numerator.bit_length()-x.denominator.bit_length();return n-(x<F(2)**n)
def literal(A):return '('+'\\'.join(','.join(format(float(x),'.17g') for x in r) for r in A)+')'
SAVE='''
void witness70(string scalar tag,transmorphic scalar C,real scalar ce,real rowvector od,
 real rowvector c,real matrix G,real colvector h,real scalar ne,struct _C58_Cert scalar fit) {
 real matrix packed
 real rowvector a
 real scalar j
 packed=J(cols(c),129,0)
 for(j=1;j<=cols(c);j++) {a=asarray(C,(1,j));packed[j,1..cols(a)]=a;}
 savebin(tag+"_exact_cost.bin",packed);savebin(tag+"_ce.bin",ce);savebin(tag+"_od.bin",od)
 savebin(tag+"_c.bin",c);savebin(tag+"_G.bin",G);savebin(tag+"_h.bin",h)
 savebin(tag+"_meta.bin",(fit.status,fit.reason,ne,fit.ops,fit.maxlimbs))
 savebin(tag+"_num.bin",fit.numerator);savebin(tag+"_den.bin",fit.denominator)
 savebin(tag+"_basis.bin",fit.basis);savebin(tag+"_dual.bin",fit.dual_numerator)
 savebin(tag+"_basis_den.bin",fit.basis_denominator);savebin(tag+"_exponent.bin",fit.row_exponent)
 savebin(tag+"_face.bin",fit.face);savebin(tag+"_trace.bin",fit.trace62)
 savebin(tag+"_ray.bin",fit.ray);savebin(tag+"_ray_num.bin",fit.ray_num);savebin(tag+"_ray_den.bin",fit.ray_den)
}
'''
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',required=True,type=Path);a=ap.parse_args()
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';shutil.copytree(REPO/'ado',pkg)
 manifest={str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (REPO/'ado').iterdir() if p.is_file()}
 specs=json.loads((REPO/'tests/fixtures/crosseff_objective.json').read_text())['cases'];cases=[]
 code='version 16\nclear all\nset more off\nset processors 1\nset type double\nset linesize 200\n'
 for k in ['PLUS','PERSONAL','SITE','OLDPLACE']:(out/k).mkdir();code+=f'sysdir set {k} "{out/k}"\n'
 code+=f'adopath ++ "{pkg}"\n'+WRITER+'input double x y\n1 1\n2 2\nend\ndea_crosseff x=y,eng(simplex) gen(warm)\nmata:\n'+SAVE+'\nend\n'
 for z in specs:
  data=[[F(float(v)) for v in row[1:]] for row in z['data']];X=[r[:2] for r in data];Y=[r[2:] for r in data];i=z['index']-1;free=z['rts']=='vrs';output=z['orientation']=='output';n=4+free
  E=[([F(0)]*2+Y[i] if output else X[i]+[F(0)]*2)+([F(0)] if free else [])];b=[F(1)]
  c=(X[i]+[F(0)]*2 if output else [F(0)]*2+[-v for v in Y[i]])+([F(1 if output else -1)] if free else [])
  U=[[-v for v in x]+y+([F(-1 if output else 1)] if free else []) for x,y in zip(X,Y)]
  xp=list(map(F,z['primary_x']));dual=z['primary_dual'];eq=list(map(F,dual['equality']));iq=list(map(F,dual['inequality']));lo={int(k):F(v) for k,v in dual['lower'].items()}
  assert dot(E[0],xp)==1 and all(dot(r,xp)<=0 for r in U) and all(v>=0 for v in xp[:4])
  assert all(v<=0 for v in iq) and all(v>=0 for v in lo.values())
  assert [eq[0]*E[0][k]+sum((r[k]*v for r,v in zip(U,iq)),F(0))+lo.get(k,F(0)) for k in range(n)]==c
  assert dot(c,xp)==eq[0]==F(z['primary_objective'])
  Fm=E+ [r for r,v in zip(U,iq) if v<0]+[[-F(k==j) for k in range(n)] for j,v in lo.items() if v>0];fb=b+[F(0)]*(len(Fm)-1)
  ex=[-fl(v) if v else 0 for v in data[i]]+([0] if free else []);D=[F(2)**e for e in ex]
  T=[[D[j] if k==j else F(0) for j in range(n)] for k in range(n)]
  rows=Fm+U;rhs=fb+[F(0)]*len(U);new=[];nr=[]
  for row,v in zip(rows,rhs):
   r=[x*d for x,d in zip(row,D)];es=[fl(abs(x)) for x in r if x]+([fl(abs(v))] if v else []);e=-((max(es)+min(es))//2) if es else 0
   new.append([x*F(2)**e for x in r]);nr.append(v*F(2)**e)
  for A in [T,new]:assert all(F(float(v))==v and math.isfinite(float(v)) and abs(float(v))<float.fromhex('0x1.0000000000000p+1023') for row in A for v in row)
  path=out/z['name'];path.mkdir();prefix=str(path/'fit')
  body='mata:\nvoid probe70() {\nreal matrix X,Y,T,E,U,G\nreal colvector b,ub,h\nreal rowvector c,od\ntransmorphic scalar C\nreal scalar ce\nstruct _C58_Cert scalar fit\n'
  for k,A in [('X',X),('Y',Y),('T',T),('E',new[:len(Fm)]),('U',new[len(Fm):]),('b',[[v] for v in nr[:len(Fm)]]),('ub',[[v] for v in nr[len(Fm):]])]:body+=k+'='+literal(A)+'\n'
  body+=f'_c70_build(X,Y,{i+1},{int(output)},{int(z["method"]=="aggressive")},T,c,C,ce,od)\nassert(od[2]==1)\n_c70_fit(c,E,b,U,ub,1e-10,C,ce,od,fit,G,h)\nwitness70("{prefix}",C,ce,od,c,G,h,rows(E),fit)\n}}\nend\nmata: probe70()\nmata: mata drop probe70()\n';code+=body
  peer=[sum((r[k] for j,r in enumerate(data) if j!=i),F(0)) for k in range(4)];minc=peer[:2]+([-v for v in peer[2:]] if not output else [F(0)]*2)+([F(0)] if free else [])
  if z['method']=='aggressive':minc=[-v for v in minc]
  cases.append(dict(name=z['name'],D=D,target=[-v*d for v,d in zip(minc,D)],E=new[:len(Fm)],b=nr[:len(Fm)],U=new[len(Fm):],ub=nr[len(Fm):],primary_E=E,primary_U=U,primary_c=c,primary_opt=F(z['primary_objective'])))
 code+='di "OBJECTIVE70_PASS"\nexit,clear\n';(out/'test.do').write_text(code);start=time.monotonic();q=subprocess.run([a.stata,'-b','do','test.do'],cwd=out,capture_output=True,timeout=1800)
 log=(out/'test.log').read_text(errors='replace');assert q.returncode==0 and 'OBJECTIVE70_PASS' in log.splitlines(),'\n'.join(log.splitlines()[-40:])
 rows=[];verified=rounded=0
 for z in cases:
  p=out/z['name']/'fit';read=lambda k:readbin(Path(str(p)+'_'+k+'.bin'));meta=read('meta').ravel();ce=int(read('ce')[0,0]);C=[F(integer(r))*F(2)**ce for r in read('exact_cost')];assert C==z['target']
  od=read('od').ravel();assert od[1]==1 and od[4]<=10000000 and od[5]<=128;rounded+=int(od[2]>1)
  G=[[F(float(v)) for v in row] for row in read('G')];h=[F(float(v)) for v in read('h').ravel()];ne=int(meta[2]);n=len(C)
  assert G==z['E']+z['U']+[[-F(k==j) for k in range(n)] for j in range(n)] and h==z['b']+z['ub']+[F(0)]*n
  trace=read('trace').ravel();assert trace[21]==sum(trace[16:21]) and trace[12]==sum(trace[[6,8,10,11,24]]) and trace[12]<=16 and meta[3]<=80000000 and meta[4]<=128
  assert trace[15]==0,'All fixture intended secondary models are finite'
  if meta[0]==1:
   dx=integer(read('den').ravel());x=[F(integer(r),dx) for r in read('num')];assert all(dot(r,x)==v for r,v in zip(G[:ne],h[:ne])) and all(dot(r,x)<=v for r,v in zip(G[ne:],h[ne:]))
   ids=[int(v)-1 for v in read('basis').ravel()];ex=list(map(int,read('exponent').ravel()));db=integer(read('basis_den').ravel());dual=[F(0)]*len(G)
   for k,row in enumerate(read('dual')):dual[ids[k]]=F(integer(row),db)*F(2)**(ex[-1]-ex[k])
   assert all(v>=0 for v in dual[ne:]) and [sum((r[j]*v for r,v in zip(G,dual)),F(0)) for j in range(n)]==C and dot(C,x)==dot(h,dual)
   original=[v*d for v,d in zip(x,z['D'])];assert dot(z['primary_E'][0],original)==1 and all(dot(r,original)<=0 for r in z['primary_U']) and dot(z['primary_c'],original)==z['primary_opt'];verified+=1
  rows.append(dict(test=z['name'],passed=True,exact_status=int(meta[0]),reason=int(meta[1]),numeric_state=int(od[2]),build_ops=float(od[4]),total_ops=float(meta[3])))
 assert len(rows)==len(specs) and all(z['exact_status']==1 for z in rows)
 for n,h in manifest.items():assert hashlib.sha256((REPO/n).read_bytes()).hexdigest()==h
 (out/'manifest.json').write_text(json.dumps(dict(source_sha256=manifest),indent=2)+'\n');(out/'results.json').write_text(json.dumps(rows,indent=2)+'\n')
 summary=dict(passed=True,cases=len(specs),exact_peer_sum_vectors_verified=len(specs),finite_certificates=verified,numerical_rounding_vectors=rounded,all_certificates_use_intended_objective=True,unverified=len(rows)-verified,budget_ledger_verified=True,source_files_preserved=True,seconds=time.monotonic()-start)
 (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
