from fractions import Fraction as F
from functools import lru_cache
import numpy as np
def fm(A):return tuple(tuple(F(float(v)) for v in row) for row in np.atleast_2d(A))
def fv(a):return tuple(F(float(v)) for v in np.ravel(a))
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def trans(A):return tuple(zip(*A))
@lru_cache(None)
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
@lru_cache(None)
def solve(A,b):
 n=len(A)
 if any(len(row)!=n for row in A) or len(b)!=n:return None
 T=[list(row)+[rhs] for row,rhs in zip(A,b)]
 for col in range(n):
  p=next((i for i in range(col,n) if T[i][col]),None)
  if p is None:return None
  T[col],T[p]=T[p],T[col];v=T[col][col];T[col]=[x/v for x in T[col]]
  for i in range(n):
   if i==col:continue
   v=T[i][col]
   if v:T[i]=[a-v*b for a,b in zip(T[i],T[col])]
 return tuple(row[-1] for row in T)
def basis(A,basic):return tuple(tuple(row[k] for k in basic) for row in A)
