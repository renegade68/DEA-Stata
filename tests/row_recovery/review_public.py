from runlib import *
from exact_utils import fm,fv,dot
from collections import Counter
cases={z['id']:z for z in j(r/'GRID_FIXTURE.json')};out=r/'checks/public';obs=r/'checks/public_observed';rows=[];parity=0
for cmd in j(r/'PUBLIC_COMMANDS.json'):
 id=cmd['id'];z=cases[cmd['fixture']];rc=int(readbin(out/(id+'_rc.bin'))[0,0]);row=dict(**cmd,rc=rc)
 if cmd['policy']==1:
  files=sorted(p.name for p in out.glob(id+'_*.bin'))
  for n in files:assert h(out/n)==h(obs/n),n;parity+=1
 if rc==0:
  res=readbin(out/(id+'_result.bin')).ravel();meta=readbin(out/(id+'_success.bin'));status=int(res[0])
  row.update(status=status,iterations=int(res[2]) if np.isfinite(res[2]) else None,outcomes=meta[:,0].astype(int).tolist(),reasons=meta[:,1].astype(int).tolist())
  if status==0:
   x=fv(readbin(out/(id+'_x.bin')));E=fm(np.array(z['E']).reshape(-1,len(z['c'])));b=fv(z['b']);c=fv(z['c']);known=F(z['objective']);obj=dot(c,x)
   row.update(exact_primal=all(dot(ar,x)==bb for ar,bb in zip(E,b)),objective=str(obj),objective_relative_error=float(abs(obj-known)/(1+abs(known))))
   if cmd['policy']==1:
    observed=list(obs.glob(id+'_obs*_meta.bin'))
    assert len(observed)==len(meta)
    for p in observed:
     prefix=p.name.removesuffix('_meta.bin');M=readbin(p);A=fm(readbin(obs/(prefix+'_A.bin')));bb=fv(readbin(obs/(prefix+'_b.bin')));cc=fv(readbin(obs/(prefix+'_c.bin')))
     assert A==E and bb==b and cc==c,(id,'observer original units differ')
 else:
  assert z['kind'].startswith('no_rows') or z['kind']=='unbounded',row
 rows.append(row)
for z in cases.values():
 for scale in ['none','row']:
  triple=[q for q in rows if q['fixture']==z['id'] and q['scale']==scale]
  assert len(triple)==3
  if triple[1].get('outcomes')==[-1]:
   assert triple[0]['status']==triple[1]['status']==0 and triple[2]['status']==3
   for suffix in ['result.bin','x.bin']:assert h(out/(triple[0]['id']+'_'+suffix))==h(out/(triple[1]['id']+'_'+suffix))
summary=dict(passed=True,commands=len(rows),observed_commands=632,observed_numeric_files_identical=parity,rc_counts=dict(Counter(str(z['rc']) for z in rows)),policy_statuses={str(p):dict(Counter(str(z.get('status','rc'+str(z['rc']))) for z in rows if z['policy']==p)) for p in [0,1,2]},retained_unavailable=sum(z['policy']==1 and z.get('outcomes')==[-1] for z in rows),retained_unavailable_objective_exceeds=sum(z['policy']==1 and z.get('outcomes')==[-1] and z.get('objective_relative_error',0)>1e-10 for z in rows),rows=rows)
dump('PUBLIC_REVIEW.json',summary);print(json.dumps({k:v for k,v in summary.items() if k!='rows'},indent=2),flush=True)
