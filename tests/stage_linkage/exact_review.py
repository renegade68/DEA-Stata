from runlib import *
from oracle import *
from collections import Counter
d=r/'checks/observed';data={z['kind']:z for z in j(r/'GRID_DATA.json')};grid=j(r/'GRID_COMMANDS.json');bounds_checked=0;fallback=0;proof_count=0;records=[]
def certificate(A,b,c,basic):
 global fallback,proof_count
 ids=tuple(int(v)-1 for v in np.ravel(basic));p=None
 if len(ids)==len(A) and len(set(ids))==len(ids) and all(0<=v<len(c) for v in ids):
  B=basis(A,ids);xb=solve(B,b);y=solve(trans(B),tuple(c[i] for i in ids))
  if xb is not None and y is not None and min(xb)>=0 and min(cc-dot(y,col) for cc,col in zip(c,trans(A)))>=0:
   x=[F(0)]*len(c)
   for i,v in zip(ids,xb):x[i]=v
   p=dict(feasible=True,x=list(map(str,x)),y=list(map(str,y)),objective=str(dot(c,x)),basis=list(ids))
 if p is None:p=lp(A,b,c);fallback+=1
 if p['feasible']:
  x=tuple(map(F,p['x']));y=tuple(map(F,p['y']))
  assert min(x)>=0 and all(dot(a,x)==bb for a,bb in zip(A,b))
  assert all(cc-dot(y,col)>=0 for cc,col in zip(c,trans(A))) and dot(c,x)==dot(b,y)==F(p['objective'])
 else:
  y=tuple(map(F,p['farkas']));assert dot(b,y)<0 and all(dot(y,col)>=0 for col in trans(A))
 proof_count+=1;return p
def saved(id,phase):
 return tuple(fn(readbin(d/f'{id}_{phase}_{name}.bin')) for name,fn in [('A',fm),('b',fv),('c',fv),('x',fv),('basis',lambda x:x)])
for z in grid:
 if z['stage']!=2:continue
 id=z['id'];meta=readbin(d/f'{id}_link.bin')[:z['DMUs']];v=data[z['dataset']];raw=np.array(v['data']);xi=[v['names'].index(n) for n in v['inputs']];yi=[v['names'].index(n) for n in v['outputs']];XR=raw[:,xi];YR=raw[:,yi];xs=XR.max(axis=0);ys=YR.max(axis=0);X=XR/xs;Y=YR/ys;n=len(raw);mi=len(xi);so=len(yi)
 L=readbin(d/f'{id}_m_lambda.bin');SI=readbin(d/f'{id}_m_slack_in.bin');SO=readbin(d/f'{id}_m_slack_out.bin');D=readbin(d/f'{id}_m_numerical_diagnostics.bin');scores=readbin(d/f'{id}_scores.bin').ravel()
 for k,m in enumerate(meta):
  key=f'{id}_d{k+1:02}'
  for phase in ['old','retry']:
   f=d/f'{key}_{phase}_bounds.bin'
   if not f.exists():continue
   A,b,c,x,basic=saved(key,phase);ids=tuple(int(a)-1 for a in basic.ravel());xx=solve(basis(A,ids),b);assert xx is not None
   bo=readbin(f);xp=readbin(d/f'{key}_{phase}_xp.bin');detail=readbin(d/f'{key}_{phase}_detail.bin').ravel()
   assert all(F(float(lo))<=true<=F(float(hi)) for true,(lo,hi) in zip(xx,bo)),key
   assert max(abs(true-F(float(p[0]))-F(float(p[1]))) for true,p in zip(xx,xp))<=F(float(detail[1])),key
   bounds_checked+=1
   if phase=='old':assert sum(a<0 for a in xx)>=m[6]
  if m[14]!=1:continue
  A,b,c,x,basic=saved(key,'new');p1=certificate(A,b,c,basic);assert p1['feasible']
  optimum=F(p1['objective']);radial=F(float(m[16]));score=F(float(scores[k]));targetscore=optimum if z['ort']=='input' else -1/optimum
  assert abs(dot(c,x)-optimum)/(1+abs(optimum))<=F('1e-10'),key
  if z['ort']=='input':assert radial>=optimum
  else:assert radial<=-optimum
  A2,b2,c2,x2,basic2=saved(key,'s2');p2=certificate(A2,b2,c2,basic2)
  assert p2['feasible'],(key,'candidate target still infeasible')
  assert np.isfinite(L[k]).all(),(key,'stage2 unavailable')
  # Actual builder must reuse the same normalized reference data and new radial.
  assert A2[:mi]==fm(np.column_stack([X.T,np.eye(mi),np.zeros((mi,so+(z['rts'] in ['drs','irs'])))]))
  assert all(row[:n]==col for row,col in zip(A2[mi:mi+so],fm(Y.T)))
  rhs=np.r_[float(radial)*X[k] if z['ort']=='input' else X[k],Y[k] if z['ort']=='input' else float(radial)*Y[k]]
  assert b2[:mi+so]==fv(rhs) and D[k,0]==m[16] and D[k,1]==m[16]
  l=fv(L[k]);si=fv(SI[k]);ss=fv(SO[k]);px=tuple(dot(col,l) for col in trans(fm(XR)));py=tuple(dot(col,l) for col in trans(fm(YR)))
  tx=tuple(radial*t for t in fv(XR[k])) if z['ort']=='input' else fv(XR[k]);ty=fv(YR[k]) if z['ort']=='input' else tuple(radial*t for t in fv(YR[k]))
  res=tuple(a+s-t for a,s,t in zip(px,si,tx))+tuple(a-s-t for a,s,t in zip(py,ss,ty))
  den=tuple(1+abs(a)+abs(s)+abs(t) for a,s,t in zip(px,si,tx))+tuple(1+abs(a)+abs(s)+abs(t) for a,s,t in zip(py,ss,ty))
  max_scaled=max(abs(a)/b for a,b in zip(res,den));assert max_scaled<=F('1e-10'),key
  tech=max([F(0)]+[a-t for a,t in zip(px,tx)]+[t-a for t,a in zip(ty,py)])
  if z['rts']=='vrs':assert abs(sum(l)-1)<=F('1e-10')
  if z['rts']=='drs':assert sum(l)<=1+F('1e-10')
  if z['rts']=='irs':assert sum(l)>=1-F('1e-10')
  assert min(l+si+ss)>=0
  obj=-sum(s/F(float(scale)) for s,scale in zip(si,xs))-sum(s/F(float(scale)) for s,scale in zip(ss,ys));objerror=abs(obj-F(p2['objective']))/(1+abs(F(p2['objective'])))
  assert objerror<=F('1e-10'),(key,'slack objective',float(objerror))
  # Transform the normalized dual to physical units, then check against raw
  # binary64 observations and the unrounded radial*observation expression.
  scale=fv(np.r_[xs,ys,np.ones(len(b2)-mi-so)])
  yr=tuple(a/s for a,s in zip(map(F,p2['y']),scale));ar=[list(row) for row in A2]
  for i in range(mi):ar[i][:n]=fv(XR[:,i])
  for i in range(so):ar[mi+i][:n]=fv(YR[:,i])
  br=tx+ty+b2[mi+so:];cr=tuple([F(0)]*n)+tuple(-1/F(float(t)) for t in np.r_[xs,ys])+tuple([F(0)]*(len(c2)-n-mi-so))
  xr=l+si+ss
  if z['rts']=='drs':xr+= (max(F(0),1-sum(l)),)
  if z['rts']=='irs':xr+= (max(F(0),sum(l)-1),)
  gap=abs(dot(cr,xr)-dot(yr,br))/(1+abs(dot(cr,xr)))
  weighted=sum(abs(w)*abs(dot(row,xr)-b) for row,b,w in zip(ar,br,yr))/(1+abs(dot(cr,xr)))
  prices=max([F(0)]+[max(F(0),dot(yr,col)-cc)/(1+abs(cc)+sum(abs(y*a) for y,a in zip(yr,col))) for cc,col in zip(cr,trans(ar))])
  assert max(gap,weighted,prices)<=F('1e-10'),(key,'raw primal/dual quality')
  rec=dict(case=id,DMU=k+1,stage1_proof=p1,stage2_proof=p2,score_absolute_error=float(abs(score-targetscore)),score=float(score),radial=float(radial),max_original_residual=float(max(map(abs,res))),max_original_scaled_residual=float(max_scaled),max_original_technology_violation=float(tech),relative_slack_objective_error=float(objerror),raw_relative_dual_gap=float(gap),raw_relative_weighted_primal=float(weighted),raw_scaled_dual_violation=float(prices))
  records.append(rec)
  if len(records)%30==0:print('Exact candidate checks',len(records),'fallbacks',fallback,flush=True)
out=dict(passed=True,repaired_main_stage2_records=len(records),basis_error_enclosures_verified=bounds_checked,exact_witness_checks=proof_count,enumeration_fallbacks=fallback,max_score_error=max(x['score_absolute_error'] for x in records),max_original_scaled_residual=max(x['max_original_scaled_residual'] for x in records),max_relative_slack_objective_error=max(x['relative_slack_objective_error'] for x in records),max_raw_dual_gap=max(x['raw_relative_dual_gap'] for x in records),max_raw_weighted_primal=max(x['raw_relative_weighted_primal'] for x in records),max_raw_scaled_dual_violation=max(x['raw_scaled_dual_violation'] for x in records),records=records)
dump('EXACT_REVIEW.json',out);print(json.dumps({k:v for k,v in out.items() if k!='records'},indent=2),flush=True)
