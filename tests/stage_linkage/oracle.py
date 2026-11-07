from exact_utils import *
from functools import lru_cache
from itertools import combinations
@lru_cache(None)
def bases(A,c):
 m=len(A);n=len(c);out=[]
 for ids in combinations(range(n),m):
  B=basis(A,ids);cols=[solve(B,tuple(F(int(i==j)) for i in range(m))) for j in range(m)]
  if cols[0] is None:continue
  inv=trans(cols);y=tuple(dot(tuple(c[k] for k in ids),col) for col in trans(inv))
  rc=tuple(cc-dot(y,col) for cc,col in zip(c,trans(A)))
  out.append((ids,inv,y,rc))
 return tuple(out)
def lp(A,b,c,proof=True):
 feasible=[]
 for ids,inv,y,rc in bases(A,c):
  xB=tuple(dot(row,b) for row in inv)
  if min(xB)<0:continue
  x=[F(0)]*len(c)
  for k,z in zip(ids,xB):x[k]=z
  assert all(dot(row,x)==rhs for row,rhs in zip(A,b))
  if not proof:return True
  feasible.append((dot(c,x),ids,tuple(x),y,rc))
 if not proof:return False
 if feasible:
  val=min(z[0] for z in feasible);opts=[z for z in feasible if z[0]==val];z=next(z for z in opts if min(z[4])>=0)
  assert dot(b,z[3])==val
  return dict(feasible=True,objective=str(val),x=list(map(str,z[2])),y=list(map(str,z[3])),basis=list(z[1]),optimal_bases=len(opts),distinct_optimal_vertices=len({z[2] for z in opts}),feasible_bases=len(feasible))
 for active in combinations(range(len(c)),len(A)-1):
  M=tuple(tuple(A[i][col] for i in range(len(A))) for col in active)+(b,)
  y=solve(M,tuple([F(0)]*(len(A)-1)+[F(-1)]))
  if y is not None and all(dot(y,col)>=0 for col in trans(A)):
   assert dot(b,y)==-1
   return dict(feasible=False,farkas=list(map(str,y)),feasible_bases=0)
 raise AssertionError('No primal or Farkas proof found')
