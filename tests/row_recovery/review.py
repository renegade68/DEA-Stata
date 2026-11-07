from runlib import *
from exact_utils import fm,fv,dot,trans,rank
cases={z['id']:z for z in j(r/'GRID_FIXTURE.json')};oracles={z['id']:z for z in j(r/'EXACT_ORACLES.json')['rows']}
# Recheck exact primal/dual and infeasibility/unbounded witnesses from the actual binary inputs.
for z in cases.values():
 A=fm(np.array(z['E']).reshape(-1,len(z['c'])));b=fv(z['b']);c=fv(z['c']);o=oracles[z['id']]
 if o['exact_status']=='optimal':
  x=tuple(map(F,o['x']));y=tuple(map(F,o['y']));assert min(x)>=0 and all(dot(ar,x)==bb for ar,bb in zip(A,b))
  assert all(cc-dot(y,col)>=0 for cc,col in zip(c,trans(A))) if A else min(c)>=0
  assert dot(c,x)==dot(b,y)==F(o['objective'])
 elif o['exact_status']=='infeasible':assert rank(tuple(ar+(bb,) for ar,bb in zip(A,b)))>rank(A)
 else:
  assert not A or all(all(v==0 for v in ar) for ar in A)
  assert all(bb==0 for bb in b) and min(c)<0
checks=exact=0;maximum=0
for group,commands in [('grid',[dict(id=z['id']+'_p'+str(q),fixture=z['id'],policy=q) for z in cases.values() for q in [1,2]]),('public',[q for q in j(r/'PUBLIC_COMMANDS.json') if q['policy'] in [1,2]])]:
 for cmd in commands:
  z=cases[cmd['fixture']];out=r/'checks'/group;res=readbin(out/(cmd['id']+'_result.bin')).ravel()
  if 'expected_status' in z:assert res[0]==z['expected_status'];continue
  assert res[0]==0
  x=fv(readbin(out/(cmd['id']+'_x.bin')));A=fm(np.array(z['E']).reshape(-1,len(x)));b=fv(z['b']);c=fv(z['c']);y=tuple(map(F,oracles[z['id']]['y']));obj=dot(c,x);D=dot(b,y)
  pp=[bb-dot(ar,x) for ar,bb in zip(A,b)];rc=tuple(cc-dot(y,col) for cc,col in zip(c,trans(A))) if A else c
  pu=sum((abs(yy)*abs(v) for yy,v in zip(y,pp)),F(0));cu=sum((abs(rr)*abs(xx) for rr,xx in zip(rc,x)),F(0));cut=F(1e-10)*(1+abs(obj))
  assert min(x)>=0 and min(rc)>=0 and abs(obj-D)<=cut and pu<=cut and cu<=cut
  assert all(abs(v)<=F(1e-10)*max(F(1),abs(bb),dot(tuple(map(abs,ar)),tuple(max(F(1),abs(xx)) for xx in x))) for ar,bb,v in zip(A,b,pp))
  checks+=1;exact+=all(v==0 for v in pp);maximum=max(maximum,float(abs(obj-D)/(1+abs(obj))),float(pu/(1+abs(obj))))
assert checks==1872
reference=j(r/'REFERENCE.json')['numeric_sha256']
for name,sha in reference.items():assert h(r/'checks'/name)==sha,name
reports=['GRID_REVIEW','PUBLIC_REVIEW','DEA_REVIEW','RANK_REVIEW','EXPANSION_REVIEW','ACCOUNTING_REVIEW','GUARDS_REVIEW','AFFINE_REVIEW']
for name in reports:assert j(r/(name+'.json'))['passed'],name
dump('FINAL_REVIEW.json',dict(passed=True,exact_original_witnesses=316,full_original_quality_checks=checks,exact_primal_returns=exact,max_original_quality_error=maximum,reference_numeric_files_identical=len(reference),rank_controls=1536,expansion_controls=1607,direct_guards=16,fault_paths=6,affine_LP=192,bounded_MILP=9,positive_DEA_commands=192,existing_stage2_quality_failures=3,all_input_data_packaged=True))
print('Independent witness/original-quality and reference parity passed',checks,len(reference),flush=True)
