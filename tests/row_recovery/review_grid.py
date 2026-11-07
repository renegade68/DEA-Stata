from runlib import *
from exact_utils import fm,fv,dot,rank,solve,trans,basis
from collections import Counter

def certificate(A,b,c,basic):
 m=len(A);n=len(c);k=len(basic)
 if k==0:
  if all(bb==0 for bb in b) and all(cc>=0 for cc in c):return dict(available=True,rows=[],x=[F(0)]*n,y=[F(0)]*m)
  return dict(available=False,reason='empty_basis_not_trivial')
 if len(set(basic))!=k or min(basic)<0 or max(basic)>=n:return dict(available=False,reason='invalid_basis')
 Ab=basis(A,basic);keep=[];rr=0
 for i,ar in enumerate(Ab):
  if rank(tuple(Ab[j] for j in keep+[i]))>rr:keep.append(i);rr+=1
  if rr==k:break
 if rr!=k:return dict(available=False,reason='basis_column_rank')
 B=tuple(Ab[i] for i in keep);xx=solve(B,tuple(b[i] for i in keep));yy=solve(trans(B),tuple(c[i] for i in basic))
 x=[F(0)]*n;y=[F(0)]*m
 for i,val in zip(basic,xx):x[i]=val
 for i,val in zip(keep,yy):y[i]=val
 pf=min(x)>=0 and all(dot(ar,x)==bb for ar,bb in zip(A,b));df=all(cc-dot(y,col)>=0 for cc,col in zip(c,trans(A)))
 return dict(available=pf and df,primal_feasible=pf,dual_feasible=df,rows=keep,x=x,y=y)

cases=j(r/'GRID_FIXTURE.json');out=r/'checks/grid';rows=[];oracles=[];corrections=[]
for z in cases:
 E=np.array(z['E']).reshape(-1,len(z['c']));b=np.array(z['b']);c=np.array(z['c']);A=fm(E);bb=fv(b);cc=fv(c)
 # Independent floating-point optimizer is comparison evidence,not the exact oracle.
 # Exact witnesses are packaged; no external floating optimizer is called.
 for policy in [0,1,2]:
  id=z['id']+'_p'+str(policy);rc=int(readbin(out/(id+'_rc.bin'))[0,0]);assert rc==0,(id,rc)
  res=readbin(out/(id+'_result.bin')).ravel();meta=readbin(out/(id+'_success.bin'));basic=(readbin(out/(id+'_basis.bin')).ravel()-1).astype(int).tolist()
  row=dict(id=z['id'],kind=z['kind'],policy=policy,status=int(res[0]),iterations=int(res[2]),rows=len(E),columns=len(c),exact_rank=rank(A),basis_size=len(basic),basis=basic,outcomes=meta[:,0].astype(int).tolist(),reasons=meta[:,1].astype(int).tolist(),message=(out/(id+'_message.txt')).read_text().strip())
  if 'expected_status' in z:
   assert row['status']==z['expected_status'],row
  if row['status']==0:
   x=fv(readbin(out/(id+'_x.bin')));obj=dot(cc,x);known=F(z['objective']);p=[bb0-dot(ar,x) for ar,bb0 in zip(A,bb)]
   row.update(exact_primal=all(v==0 for v in p),objective=str(obj),known_objective=str(known),objective_relative_error=float(abs(obj-known)/(1+abs(known))),max_exact_primal_residual=str(max(map(abs,p),default=F(0))))
   assert min(x)>=0
   cert=certificate(A,bb,cc,basic)
   row['offline_certificate_available']=cert['available']
   if cert['available']:
    assert dot(cc,cert['x'])==dot(bb,cert['y'])==known
    pu=sum(abs(yy)*abs(pp) for yy,pp in zip(cert['y'],p));row['weighted_primal_relative_error']=float(pu/(1+abs(obj)))
    row['offline_quality_exceeds']=abs(obj-known)>F(1e-10)*(1+abs(obj)) or pu>F(1e-10)*(1+abs(obj))
    if policy==1:corrections.append(dict(id=z['id'],basis=basic,selected_rows=cert['rows'],x=[str(v) for v in cert['x']],y=[str(v) for v in cert['y']],unchanged_candidate=row['exact_primal'] and obj==known))
  rows.append(row)
 # Recording and certification policies must not change solver iteration budget.
 a,b1,c1=rows[-3:];assert len({z0['iterations'] for z0 in [a,b1,c1]})==1,z['id']
 if b1['outcomes']==[-1]:
  assert a['status']==b1['status']==0 and c1['status']==3,z['id']
  for suffix in ['result.bin','x.bin','message.txt','basis.bin']:assert h(out/(z['id']+'_p0_'+suffix))==h(out/(z['id']+'_p1_'+suffix))
 if a['status']!=0:assert a['status']==b1['status']==c1['status']
summary=dict(passed=True,conditions=len(cases),executions=len(rows),policy_statuses={str(p):dict(Counter(str(z['status']) for z in rows if z['policy']==p)) for p in [0,1,2]},retained_unavailable=sum(z['policy']==1 and z['outcomes']==[-1] for z in rows),retained_unavailable_quality_exceeds=sum(z['policy']==1 and z['outcomes']==[-1] and z.get('offline_quality_exceeds',False) for z in rows),raw_numerical_failures=sum(z['policy']==0 and z['status']==3 for z in rows),offline_certificate_records=len(corrections),rows=rows)
dump('GRID_REVIEW.json',summary);dump('OFFLINE_CERTIFICATES.json',corrections)
print(json.dumps({k:v for k,v in summary.items() if k!='rows'},indent=2),flush=True)
print('Largest unavailable gaps',sorted([z for z in rows if z['policy']==1 and z['outcomes']==[-1]],key=lambda z:z.get('objective_relative_error',0),reverse=True)[:3],flush=True)
