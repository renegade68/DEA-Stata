#!/usr/bin/env python3
"""Wide failure-only recovery and exact reversible initial scaling.
Requires NumPy/SciPy, Stata16+ and only the packaged fixtures. Tests independent
rational primal/dual witnesses, coefficient preservation and safe rejection.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,math,shutil,struct,subprocess
import numpy as np
from test_core_retries import WRITER,readbin
from test_simplex_stability import mata
from test_precision_recovery import decode,dot,fm
REPO=Path(__file__).resolve().parents[1]
def rank(A):
 T=[list(row) for row in A];m=len(T);n=len(T[0]) if m else 0;k=0
 for col in range(n):
  p=next((i for i in range(k,m) if T[i][col]),None)
  if p is None:continue
  T[k],T[p]=T[p],T[k];v=T[k][col];T[k]=[x/v for x in T[k]]
  for i in range(k+1,m):
   v=T[i][col]
   if v:T[i]=[a-v*b for a,b in zip(T[i],T[k])]
  k+=1
  if k==m:break
 return k

def boundaries(out,pkg,stata,fixture):
 out.mkdir();rows=[]
 # Use binary input, so decimal parser underflow cannot disguise a helper error.
 reader='''mata:
 real matrix loadbin(string scalar file) {
  real matrix B,A
  real scalar f,m,n
  B=bufio();bufbyteorder(B,2);f=fopen(file,"r")
  m=fbufget(B,f,"%4b",1,1);n=fbufget(B,f,"%4b",1,1)
  A=fbufget(B,f,"%8z",m,n);fclose(f);return(A)
 }
 end
 '''
 def write(name,A):
  A=np.atleast_2d(A).astype('<f8');(out/name).write_bytes(struct.pack('<ii',*A.shape)+A.tobytes())
 body='version 16\nclear all\nset more off\nset processors 1\nset linesize 180\n'
 for k in ['PLUS','PERSONAL','SITE','OLDPLACE']:(out/k).mkdir();body+=f'sysdir set {k} "{out/k}"\n'
 body+=WRITER+reader+f'adopath ++ "{pkg}"\nrun "{pkg}/_dea_lpcontract.ado"\nmata:\n'
 vectors=[(x['name'],np.array([float.fromhex(y) for y in x['values']]).reshape(1,-1)) for x in fixture['scale_vectors']]
 for i,(name,a) in enumerate(vectors):
  write(f'p{i:03}_input.bin',a)
  body+=f'v=loadbin("p{i:03}_input.bin");v0=v;mult=_dea_lpcontract_pow2(v)\nif(cols(v)>0) assert(all(v:==v0))\nsavebin("p{i:03}_mult.bin",mult);savebin("p{i:03}_scaled.bin",v*mult);savebin("p{i:03}_back.bin",(v*mult)/mult)\n'
 body+='assert(_dea_lpcontract_pow2((1,.,2))==1)\n'
 # Deterministic, well conditioned thin matrices, duplicates/zeros/ties, fixed columns.
 rng=np.random.default_rng(2026110314);selects=[]
 for m in [1,2,3,8,16,32]:
  for n in [129,256,1024,10000]:
   for mode in ['identity','dyadic']:
    A=np.zeros((m,n));A[:,:m]=np.eye(m)*(2 if mode=='identity' else 16)
    if mode=='identity':
     A[:,m:]=np.eye(m)[:,np.arange(n-m)%m]
    else:
     A[:,:m]+=rng.integers(-2,3,(m,m));A[:,m:]=rng.integers(-2,3,(m,n-m))/8
    fixed=np.zeros(n);indices=np.arange(m)[::2];fixed[indices]=1
    selects.append((f'{mode}_{m}_{n}',A,fixed,True))
 # Many fixed columns must move through swap positions; input dimension/missing controls.
 A=np.c_[np.eye(8)[:,np.arange(248)%8]*.5,np.eye(8)*2];f=np.zeros(256);f[248::2]=1;selects.append(('late_fixed',A,f,True))
 selects += [('zero_rank',np.zeros((3,129)),np.zeros(129),False),('too_many_fixed',np.eye(3,129),np.ones(129),False),('too_few_columns',np.eye(3,2),np.zeros(2),False),('empty_rows',np.zeros((0,129)),np.zeros(129),False),('wrong_fixed_size',np.eye(3,129),np.zeros(128),False)]
 for i,(name,A,fixed,expected) in enumerate(selects):
  write(f's{i:03}_A.bin',A);write(f's{i:03}_fixed.bin',fixed)
  body+=f'A=loadbin("s{i:03}_A.bin");A0=A;fixed=loadbin("s{i:03}_fixed.bin");f0=fixed\np=_rsm_retry_columns(A,fixed);if(rows(A)>0 & cols(A)>0) assert(all(A:==A0));assert(all(fixed:==f0));savebin("s{i:03}_p.bin",p)\n'
 body+='assert(rows(_rsm_retry_columns((1,.,3),(0,0,0)))==0)\nassert(rows(_rsm_retry_columns((1,2,3),(0,.,0)))==0)\nend\ndi "BOUNDARY_PASS"\nexit,clear\n'
 (out/'test.do').write_text(body);p=subprocess.run([stata,'-b','do','test.do'],cwd=out,capture_output=True,timeout=240);(out/'process.json').write_text(json.dumps(dict(returncode=p.returncode,stderr=p.stderr.decode(errors='replace')),indent=2)+'\n');log=(out/'test.log').read_text();assert p.returncode==0 and 'BOUNDARY_PASS' in log.splitlines(),log[-4500:]
 scale_results=[]
 for i,(name,a) in enumerate(vectors):
  v=np.ravel(a);mult=float(readbin(out/f'p{i:03}_mult.bin')[0,0]);z=readbin(out/f'p{i:03}_scaled.bin').ravel();back=readbin(out/f'p{i:03}_back.bin').ravel()
  assert mult>0 and math.isfinite(1/mult) and math.frexp(mult)[0]==.5
  assert np.array_equal(back,v),(name,mult,v,back)
  assert all(F(float(y))==F(float(x))*F(mult) for x,y in zip(v,z)),name
  scale_results.append(dict(name=name,multiplier=mult,exact_forward=True,exact_reverse=True,unscaled=mult==1))
 selection_results=[]
 for i,(name,A,fixed,expected) in enumerate(selects):
  p=readbin(out/f's{i:03}_p.bin').ravel();assert bool(len(p))==expected,name
  if len(p):
   ids=p.astype(int)-1;assert np.array_equal(ids+1,p) and len(set(ids))==len(A) and min(ids)>=0 and max(ids)<A.shape[1]
   assert set(ids[:int(fixed.sum())])==set(np.flatnonzero(fixed));assert rank(fm(A[:,ids]))==len(A),name
  selection_results.append(dict(name=name,accepted=bool(len(p)),exact_selected_rank=bool(len(p)),unchanged_input=True))
 summary=dict(power_two_cases=len(vectors)+1,exact_forward_and_reverse_cases=len(vectors),missing_control=1,unscaled=sum(x['unscaled'] for x in scale_results),selector_cases=len(selects)+2,selector_valid=sum(x['accepted'] for x in selection_results),invalid_or_zero_controls=7,max_rows=32,max_columns=10000,all_passed=True)
 (out/'results.json').write_text(json.dumps(dict(scale=scale_results,selection=selection_results),indent=2)+'\n');(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');return summary

def header(out,pkg):
 out.mkdir();body='version 16\nclear all\nset more off\nset processors 1\nset linesize 180\n'
 for k in ['PLUS','PERSONAL','SITE','OLDPLACE']:(out/k).mkdir();body+=f'sysdir set {k} "{out/k}"\n'
 return body+WRITER+f'adopath ++ "{pkg}"\nrun "{pkg}/_dea_simplex.ado"\n'

def run(out,body,stata,marker):
 (out/'test.do').write_text(body+f'\ndi "{marker}"\nexit,clear\n')
 p=subprocess.run([stata,'-b','do','test.do'],cwd=out,capture_output=True,timeout=900)
 (out/'process.json').write_text(json.dumps(dict(returncode=p.returncode,stderr=p.stderr.decode(errors='replace')),indent=2)+'\n')
 log=(out/'test.log').read_text(errors='replace');assert p.returncode==0 and marker in log.splitlines(),log[-4000:]

def witness(f):
 n=len(f['xstar']);E=decode(f['E'],n);b=decode(f['b'],1).ravel();c=decode(f['c'],n).ravel();x=list(map(F,f['xstar']));y=list(map(F,f['ystar']))
 assert np.all(E>0) and np.all(b>0) and np.all(c>0)
 assert all(dot(ar,x)==rhs for ar,rhs in zip(fm(E),map(F,b)))
 assert min(F(ci)-dot(y,list(map(F,col))) for ci,col in zip(c,E.T))>=0
 assert dot(list(map(F,c)),x)==dot(y,list(map(F,b)))==6
 return E,b,c,x,y

def wide(out,pkg,stata,fixture,cases):
 body=header(out,pkg)+'mata:\n'+r'''
void wide_lp(string scalar name,real matrix E0,real colvector b,real rowvector c0,real scalar n,real scalar tol) {
 struct RSM_Prob scalar P
 struct RSM_Result scalar R
 real matrix E
 real rowvector c
 real scalar j,k
 E=J(rows(E0),n,0);c=J(1,n,0);E[,1..cols(E0)]=E0;c[1..cols(c0)]=c0
 for(j=cols(E0)+1;j<=n;j++) {k=mod(j-1,cols(E0))+1;E[,j]=E0[,k];c[j]=c0[k]+64;}
 _rsm_std_form(P,J(0,n,0),J(0,1,0),E,b,c);R=_rsm_two_phase(P,tol,"")
 savebin(name+"_result.bin",(R.status,R.objval,R.iters));savebin(name+"_x.bin",R.x)
 savebin(name+"_E.bin",E);savebin(name+"_b.bin",b);savebin(name+"_c.bin",c)
}
'''
 for i,s in enumerate(fixture['wide_seeds']):
  E,b,c,xx,y=witness(cases[s['fixture_index']])
  for n in fixture['widths']:body+=f'wide_lp("w{i}_{n}",{mata(E)},{mata(b.reshape(-1,1))},{mata(c.reshape(1,-1))},{n},{s["tol"]})\n'
 run(out,body+'end\n',stata,'WIDE_PASS');rows=[]
 for i,s in enumerate(fixture['wide_seeds']):
  E0,b0,c0,target,y=witness(cases[s['fixture_index']])
  for n in fixture['widths']:
   label=f'w{i}_{n}';E=readbin(out/(label+'_E.bin'));b=readbin(out/(label+'_b.bin')).ravel();c=readbin(out/(label+'_c.bin')).ravel();a=readbin(out/(label+'_result.bin')).ravel();xf=readbin(out/(label+'_x.bin')).ravel();xx=target+[F(0)]*(n-len(target))
   assert np.array_equal(E[:,:len(target)],E0) and np.array_equal(b,b0) and np.array_equal(c[:len(target)],c0)
   assert np.all(E>0) and np.all(b>0) and np.all(c>0)
   A=fm(E);bf=list(map(F,b));cf=list(map(F,c));assert all(dot(ar,xx)==rhs for ar,rhs in zip(A,bf));assert min(ci-dot(y,col) for ci,col in zip(cf,zip(*A)))>=0
   assert dot(cf,xx)==dot(y,bf)==6
   expected=3 if n+2*len(b)>4096 else 0;assert int(a[0])==expected,(i,n,a)
   if expected==0:assert list(map(F,xf))==xx and F(a[1])==6
   rows.append(dict(seed=i,width=n,status=int(a[0]),exact_primal=expected==0,passed=True))
 s=dict(cases=len(rows),available=sum(x['status']==0 for x in rows),exact_primals=sum(x['exact_primal'] for x in rows),size_limit_controls=3,all_passed=True)
 (out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');(out/'summary.json').write_text(json.dumps(s,indent=2)+'\n');return s

def public(out,pkg,stata,fixture,cases):
 body=header(out,pkg)
 for z in fixture['public']:
  f=cases[z['fixture_index']];E,b,c,x,y=witness(f);assert len(x)==3
  label=f'c{z["case"]:03}_{z["scale"]}'
  body+='clear\ninput double x1 double x2 double x3 str2 rel double rhs\n'+' '.join(format(v,'.17g') for v in c)+' "" .\n'
  for ar,rhs in zip(E,b):body+=' '.join(format(v,'.17g') for v in ar)+f' "=" {rhs:.17g}\n'
  body+=f'end\nmata: savebin("{label}_input.bin",st_data(.,tokens("x1 x2 x3 rhs")))\nlp x1 x2 x3,min engine(simplex) scaling({z["scale"]}) tolerance({f["tol"]})\n'
  body+=f'mata: savebin("{label}_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)")))\nmata: savebin("{label}_x.bin",st_matrix("r(x)"))\n'
 run(out,body,stata,'PUBLIC_PASS');rows=[]
 for z in fixture['public']:
  f=cases[z['fixture_index']];E,b,c,target,y=witness(f);label=f'c{z["case"]:03}_{z["scale"]}';inp=readbin(out/(label+'_input.bin'));assert np.array_equal(inp[0,:3],c) and np.array_equal(inp[1:,:3],E) and np.array_equal(inp[1:,3],b)
  a=readbin(out/(label+'_result.bin')).ravel();row=dict(case=z['case'],scale=z['scale'],status=int(a[0]),passed=True);assert a[0] in [0,3]
  if z['required_available']:assert a[0]==0,label
  if a[0]==0:
   x=list(map(F,readbin(out/(label+'_x.bin')).ravel()));assert len(x)==3 and min(x)>=0;obj=dot(list(map(F,c)),x)
   for ar,rhs in zip(fm(E),map(F,b)):
    den=max(abs(rhs),dot(list(map(abs,ar)),[max(F(1),abs(v)) for v in x]));assert abs(dot(ar,x)-rhs)<=F(f['tol'])*den
   assert abs(obj-6)/7<=F(1e-8) and abs(obj-F(a[1]))<=F(1e-12)*(1+abs(obj))
   if z['required_exact']:assert x==target and obj==6,label
   row.update(exact_primal=x==target,exact_objective=str(obj),relative_objective_error=float(abs(obj-6)/7))
  rows.append(row)
 s=dict(cases=len(rows),available=sum(x['status']==0 for x in rows),exact_primals=sum(x.get('exact_primal',False) for x in rows),required_available=sum(x['required_available'] for x in fixture['public']),required_exact=sum(x['required_exact'] for x in fixture['public']),all_passed=True)
 (out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');(out/'summary.json').write_text(json.dumps(s,indent=2)+'\n');return s

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',required=True,type=Path);a=ap.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';shutil.copytree(REPO/'ado',pkg)
 paths={n:REPO/'tests/fixtures'/n for n in ['precision_recovery.json','wide_scaling.json']};fixture=json.loads(paths['wide_scaling.json'].read_text());cases=json.loads(paths['precision_recovery.json'].read_text())['cases']
 manifest=dict(source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in pkg.iterdir() if p.is_file()},fixtures_sha256={n:hashlib.sha256(p.read_bytes()).hexdigest() for n,p in paths.items()},numpy=np.__version__,stata=a.stata)
 (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 s=dict(boundaries=boundaries(out/'boundaries',pkg,a.stata,fixture),wide=wide(out/'wide',pkg,a.stata,fixture,cases),public=public(out/'public',pkg,a.stata,fixture,cases),all_passed=True)
 s['checks']=s['boundaries']['power_two_cases']+s['boundaries']['selector_cases']+s['wide']['cases']+s['public']['cases'];assert s['checks']==269
 (out/'summary.json').write_text(json.dumps(s,indent=2)+'\n');print(json.dumps(s,indent=2))
if __name__=='__main__':main()
