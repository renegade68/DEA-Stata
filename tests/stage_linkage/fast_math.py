from runlib import *
from exact_utils import *
cases=[]
for z in j(r/"FAST_FIXTURE.json"):
 cases.append((np.array([[float.fromhex(v) for v in row] for row in z["A"]]),np.array([float.fromhex(v) for v in z["b"]]),np.array([float.fromhex(v) for v in z["x"]]),z["kind"]))
body=f'run "{r}/candidate/ado/_dea_simplex.ado"\nmata:\n'+r'''
void fastprobe(string scalar id,real matrix A,real colvector b,real colvector x) {
 struct RSM_Prob scalar P
 P.A=A;P.b=b;P.m_con=rows(A);P.n_orig=cols(A)
 savebin(id+"_result.bin",(_dea_link_small(P,x),_dea_link_residual(P,x)))
}
end
'''
meta=[]
for i,(A,b,x,kind) in enumerate(cases):
 id=f'f{i:04}';body+='mata:\n'
 for name,V in [('AA',A),('bb',b.reshape(-1,1)),('xx',x.reshape(-1,1))]:
  body+=f'{name}=J({V.shape[0]},{V.shape[1]},0)\n'
  for ri in range(V.shape[0]):
   for cj in range(0,V.shape[1],32):
    end=min(cj+32,V.shape[1]);body+=f'{name}[{ri+1},{cj+1}..{end}]={mata(V[ri:ri+1,cj:end])}\n'
 body+=f'fastprobe("{id}",AA,bb,xx)\nend\n';meta.append(dict(id=id,kind=kind,A=[[v.hex() for v in row] for row in A],b=[v.hex() for v in b],x=[v.hex() for v in x]))
dump('FAST_FIXTURE.json',meta);d=run('fast_math',body,r/'candidate/ado');accepted=available=0;rows=[]
for (A,b,x,kind),z in zip(cases,meta):
 result=readbin(d/(z['id']+'_result.bin')).ravel();af=fm(A);bf=fv(b);xf=fv(x)
 exact=max(abs(bb-dot(a,xf))/(1+abs(bb)+sum(abs(u*v) for u,v in zip(a,xf))) for a,bb in zip(af,bf))
 gate=F(64*(A.shape[0]+A.shape[1]+4))*F(float(np.finfo(float).eps))
 if result[0]==1:assert exact<gate,z['id'];accepted+=1
 if np.isfinite(result[1:]).all():
  assert F(float(result[1]))<=exact<=F(float(result[2])),(z['id'],result,float(exact));available+=1
 rows.append(dict(id=z['id'],kind=kind,fast=bool(result[0]),bounded=bool(np.isfinite(result[1:]).all()),exact_ratio=str(exact),gate=str(gate)))
out=dict(passed=True,cases=len(cases),fast_accepted=accepted,accepted_outside_trigger_bound=0,residual_enclosures_verified=available,rows=rows)
dump('FAST_REVIEW.json',out);print({k:v for k,v in out.items() if k!='rows'},flush=True)
