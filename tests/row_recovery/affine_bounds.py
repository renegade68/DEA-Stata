from runlib import *
from exact_utils import fm,fv,dot
import itertools
body=f'run "{r}/candidate/ado/_dea_simplex.ado"\n';cases=[];milps=[]
for k,dep,shift,scale,sense,policy in itertools.product([5,25],['sum','odd'],[(0,0,0),(-6,0,0),(2,-4,1),(-2,-1,-3)],['none','row'],['min','max'],[0,1,2]):
 d=2.**-k;E=np.array([[1,2,4],[1,2+d,4-d],[1+d,2,4-d]],dtype=float);E=np.vstack([E,E[0]+E[1] if dep=='sum' else 3*E[0]])
 L=np.array(shift,dtype=float).reshape(1,-1);xstar=np.ones(3)+L.ravel();b=E@xstar;c=np.array([[1.,2.,3.]]);ub=L+4
 if sense=='max':c=-c
 assert all(dot(ar,fv(xstar))==bb for ar,bb in zip(fm(E),fv(b)))
 obj=dot(fv(c),fv(xstar));id=f'a{len(cases):03}'
 cases.append(dict(id=id,k=k,dependency=dep,lb=L.tolist(),ub=ub.tolist(),E=E.tolist(),b=b.tolist(),c=c.tolist(),xstar=xstar.tolist(),objective=str(obj),scale=scale,sense=sense,policy=policy))
 body+=f'mata: _rsm_success_reset({policy},1);st_matrix("E",{mata(E)});st_matrix("B",{mata(b.reshape(-1,1))});st_matrix("C",{mata(c)});st_matrix("LB",{mata(L)});st_matrix("UB",{mata(ub)})\n'
 body+=f'_dea_lpcontract, c(C) aeq(E) beq(B) lb(LB) ub(UB) sense({sense}) engine(simplex) scaling({scale}) tolerance(1e-10)\n'
 body+=f'mata: savebin("{id}_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")));savebin("{id}_x.bin",st_matrix("r(x)"));savebin("{id}_success.bin",_rsm_success_read())\n'
for policy,kind in itertools.product([0,1,2],['mixed','integer','limited']):
 id=f'm_{kind}_p{policy}';integer='1 2' if kind=='integer' else '1';limit='maxnodes(1)' if kind=='limited' else ''
 body+=f'mata: _rsm_success_reset({policy},1)\nmatrix C=(-1,0)\nmatrix E=(1,1\\3,3)\nmatrix B=(1.5\\4.5)\nmatrix LB=(-1,-1)\nmatrix UB=(3,3)\n'
 body+=f'milp_bb, c(C) aeq(E) beq(B) lb(LB) ub(UB) sense(min) integer({integer}) lpengine(simplex) lptol(1e-10) {limit}\n'
 body+=f'mata: savebin("{id}_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(nodes)")));savebin("{id}_x.bin",st_matrix("r(x)"));savebin("{id}_success.bin",_rsm_success_read())\n'
 milps.append(dict(id=id,policy=policy,kind=kind))
dump('AFFINE_FIXTURE.json',cases);out=run('affine_bounds',body);rows=[]
for z in cases:
 res=readbin(out/(z['id']+'_result.bin')).ravel();x=fv(readbin(out/(z['id']+'_x.bin')));E=fm(z['E']);b=fv(z['b']);c=fv(z['c']);obj=dot(c,x);known=F(z['objective']);err=abs(obj-known)/(1+abs(obj));p=[bb-dot(ar,x) for ar,bb in zip(E,b)];meta=readbin(out/(z['id']+'_success.bin'))
 assert res[0]==0 and all(lo<=v<=hi for lo,v,hi in zip(fv(z['lb']),x,fv(z['ub']))),(z['id'],res)
 if z['policy']:
  assert err<=F(1e-10),(z['id'],err)
  assert all(abs(pp)<=F(1e-10)*max(F(1),abs(bb),dot(tuple(map(abs,ar)),tuple(max(F(1),abs(v)) for v in x))) for ar,bb,pp in zip(E,b,p))
  assert all(v in [1,2] for v in meta[:,0]),(z['id'],meta)
 rows.append(dict(id=z['id'],policy=z['policy'],status=0,exact_primal=all(v==0 for v in p),relative_objective_error=float(err),iterations=int(res[2]),outcomes=meta[:,0].tolist()))
for z in milps:
 res=readbin(out/(z['id']+'_result.bin')).ravel();x=readbin(out/(z['id']+'_x.bin')).ravel();assert res[0]==dict(mixed=0,integer=1,limited=4)[z['kind']],(z,res)
 if z['kind']=='mixed':assert res[1]==-2 and list(x)==[2,-.5] and res[2]>1,(z,res,x)
 rows.append(dict(**z,status=int(res[0]),nodes=int(res[2])))
dump('AFFINE_REVIEW.json',dict(passed=True,dependent_row_LP_commands=len(cases),bounded_MILP_commands=len(milps),actual_branching_verified=True,original_units_verified=True,rows=rows))
print('Dependent rows with affine bounds and MILP passed',len(cases),len(milps),flush=True)
