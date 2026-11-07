from runlib import *
from exact_utils import fm,fv,rank
from collections import Counter
rng=np.random.default_rng(2026110535);cases=[]
for i in range(1200):
 m=int(rng.integers(1,9));n=int(rng.integers(1,9));A=rng.integers(-8,9,(m,n)).astype(float)*np.exp2(rng.integers(-12,13,(m,n)));b=rng.integers(-8,9,m).astype(float)
 # Exact dependent and exact inconsistent augmented rows, signed odd factors.
 if i%3!=0:
  weights=rng.integers(-3,4,m);ar=weights@A;br=weights@b
  A=np.vstack([A,ar]);b=np.append(b,br)
  if i%3==2:b[-1]=np.nextafter(b[-1],np.inf)
 if i%5==0:A=np.vstack([A,np.zeros(n)]);b=np.append(b,0)
 if i%7==0:A=A[::-1].copy();b=b[::-1].copy()
 cases.append(dict(id=f'r{i:04}',kind='random',A=A.tolist(),b=b.tolist(),budget=1000000))
for z in j(r/'GRID_FIXTURE.json'):
 cases.append(dict(id=z['id'],kind='grid',A=z['E'],b=z['b'],n=len(z['c']),budget=1000000))
for k in [0,1,20,40,50,52]:
 d=2.**-k;A=[[1,1],[1,1+d],[2,2+d]]
 for bad in [False,True]:
  b=[2,2+d,4+d]
  if bad:b[-1]=np.nextafter(b[-1],np.inf)
  cases.append(dict(id=f'ulp{k}_{int(bad)}',kind='one_ulp',A=A,b=b,budget=1000000))
# Explicit bounds, including work failure before a numerical decision.
for name,A,b,budget in [('rows257',np.zeros((257,1)),np.zeros(257),1000000),('cols257',np.zeros((1,257)),[0],1000000),('rank33',np.eye(33),np.ones(33),1000000),('budget0',np.eye(2),[1,1],0),('budget100',np.eye(2),[1,1],100),('range',[[2.**400,2.**-400]],[0],1000000),('tiny_rhs',[[0,0]],[2.**-1000],1000000),('inconsistent',[[1,1],[3,3]],[1,4],1000000)]:
 cases.append(dict(id=name,kind='guard',A=np.asarray(A).tolist(),b=list(b),budget=budget))
dump('RANK_FIXTURE.json',cases)
body=f'run "{r}/candidate/ado/_dea_simplex.ado"\nmata:\n'+r'''
void mapcase(string scalar id,real matrix A,real colvector b,real scalar budget) {
 struct RSM_Prob scalar P
 real colvector kr,kc
 real rowvector detail
 real scalar ok
 P=RSM_Prob();P.A=A;P.b=b;P.m_con=rows(A);P.n_orig=cols(A)
 ok=_rsm_rowmap(P,kr,kc,detail,budget)
 savebin(id+"_map.bin",(ok,detail));savebin(id+"_rows.bin",kr);savebin(id+"_cols.bin",kc)
}
end
'''
for z in cases:
 A=np.array(z['A']).reshape(-1,z.get('n',len(z['A'][0]) if z['A'] else 0));b=np.array(z['b']).reshape(-1,1)
 ae,be=mata(A),mata(b)
 if z['id']=='rows257':ae,be='J(257,1,0)','J(257,1,0)'
 if z['id']=='cols257':ae,be='J(1,257,0)','0'
 if z['id']=='rank33':ae,be='I(33)','J(33,1,1)'
 body+=f'mata: mapcase("{z["id"]}",{ae},{be},{z["budget"]})\n'
out=r/'checks/rank' if '--reuse-rank' in sys.argv else run('rank',body)
assert j(out/'summary.json')['source_sha256']==j(r/'BUILD.json')['source_sha256']
records=[]
for z in cases:
 A=np.array(z['A']).reshape(-1,z.get('n',len(z['A'][0]) if z['A'] else 0));b=np.array(z['b']).reshape(-1,1);Q=fm(A);aug=fm(np.column_stack([A,b]));rr=rank(Q);ra=rank(aug)
 meta=readbin(out/(z['id']+'_map.bin')).ravel();kr=(readbin(out/(z['id']+'_rows.bin')).ravel()-1).astype(int).tolist();kc=(readbin(out/(z['id']+'_cols.bin')).ravel()-1).astype(int).tolist();ok=bool(meta[0])
 assert meta[4]<=z['budget'],(z['id'],'work')
 if ok:
  assert rr==ra==len(kr)==len(kc)==meta[3],(z,meta,rr,ra)
  assert rank(tuple(Q[i] for i in kr))==rr and rank(tuple(aug[i] for i in kr))==ra
  assert rank(tuple(tuple(Q[i][k] for k in kc) for i in kr))==rr
 if meta[2]==2:assert ra>rr,z['id']
 records.append(dict(id=z['id'],kind=z['kind'],available=ok,reason=int(meta[2]),rank=rr,augmented_rank=ra,work=int(meta[4]),expansion=int(meta[5])))
assert all(not q['available'] for q in records if q['kind']=='guard')
assert all(q['available']==(q['rank']==q['augmented_rank']) for q in records if q['kind'] in ['grid','one_ulp'])
assert all(q['reason'] in [0,2,12,13,14,15] for q in records if q['kind']=='random')
dump('RANK_REVIEW.json',dict(passed=True,conditions=len(records),available=sum(q['available'] for q in records),reasons=dict(Counter(str(q['reason']) for q in records)),no_false_deletion=True,rows=records));print('Exact rank validation passed',len(records),flush=True)
# Expansions formed by growing arbitrary signed finite terms; compare all components to Fraction.
ac=[]
for i in range(1600):
 count=int(rng.integers(1,13));a=np.ldexp(rng.uniform(-2,2,count),rng.integers(-180,181,count));b=np.ldexp(rng.uniform(-2,2,count),rng.integers(-180,181,count))
 if i%4==0:a=np.r_[a,-a]
 if i%7==0:b=np.r_[b,-b,np.nextafter(1.,2.)-1.]
 ac.append(dict(id=f'a{i:04}',a=a.tolist(),b=b.tolist()))
for power in [-400,-399,-100,0,100,399,400]:
 ac.append(dict(id=f'edge{power}',a=[2.**power],b=[-2.**power]))
dump('ARITHMETIC_FIXTURE.json',ac)
body=f'run "{r}/candidate/ado/_dea_simplex.ado"\nmata:\n'+r'''
void arith(string scalar id,real rowvector a,real rowvector b) {
 real scalar w,e,wa,ea,wm,em
 real rowvector x,y,z,p
 w=e=0;x=_rsm_xadd(0,a,w,e,1000000);y=_rsm_xadd(0,b,w,e,1000000)
 wa=w;ea=e;z=_rsm_xadd(x,y,wa,ea,1000000)
 wm=w;em=e;p=_rsm_xmul(x,y,wm,em,1000000)
 savebin(id+"_x.bin",x);savebin(id+"_y.bin",y);savebin(id+"_add.bin",z);savebin(id+"_mul.bin",p)
 savebin(id+"_meta.bin",(w,e,wa,ea,wm,em))
}
end
'''
for z in ac:body+=f'mata: arith("{z["id"]}",{mata(np.array(z["a"]).reshape(1,-1))},{mata(np.array(z["b"]).reshape(1,-1))})\n'
body+='''mata:
w=e=0;out=_rsm_xadd(0,(1,2,3),w,e,0);savebin("guard_work.bin",(w,e))
w=e=0;out=_rsm_xmul(2.0^401,1,w,e,1000000);savebin("guard_range.bin",(w,e))
w=e=0;out=0
for(i=0;i<=33;i++) out=_rsm_xgrow(out,2.0^(-900+54*i),w,e,1000000)
savebin("guard_capacity.bin",(w,e))
end
'''
out=r/'checks/expansion' if '--reuse-expansion' in sys.argv else run('expansion',body)
assert j(out/'summary.json')['source_sha256']==j(r/'BUILD.json')['source_sha256']
maxlen=verified=unavailable=0
for z in ac:
 a=sum(map(F,z['a']));b=sum(map(F,z['b']));meta=readbin(out/(z['id']+'_meta.bin')).ravel();assert np.all(meta[[1,3]]==0) and meta[5] in [0,3],(z['id'],meta)
 for suffix,expected in [('x',a),('y',b),('add',a+b),('mul',a*b)]:
  if suffix=='mul' and meta[5]!=0:
   assert readbin(out/(z['id']+'_mul.bin')).size==0;unavailable+=1;continue
  verified+=1
  terms=fv(readbin(out/(z['id']+'_'+suffix+'.bin')));assert sum(terms)==expected,(z['id'],suffix)
  assert (all(v==0 for v in terms))==(expected==0),(z['id'],'zero detection')
  maxlen=max(maxlen,len(terms))
for guard,code in [('work',4),('range',2),('capacity',3)]:assert readbin(out/('guard_'+guard+'.bin'))[0,1]==code,guard
dump('EXPANSION_REVIEW.json',dict(passed=True,conditions=len(ac),exact_component_checks=verified,multiplications_unavailable_capacity=unavailable,guards=3,max_components=maxlen,no_false_zero=True))
print('Exact expansion validation passed',len(ac),flush=True)
