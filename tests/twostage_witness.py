"""Exact original-data optimality witnesses, independent of the Mata solver.

The compact fixture stores basis row indices, not expected efficiency values.
Fraction Gaussian elimination reconstructs primal and dual optima and verifies
strong duality before deriving the noncooperative stage ratios.
"""
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path

def dot(a,b):
    return sum((x*y for x,y in zip(a,b)),F(0))

def solve(a,b):
    n=len(b);v=[list(row)+[rhs] for row,rhs in zip(a,b)]
    assert n and all(len(row)==n+1 for row in v)
    for j in range(n):
        i=next(i for i in range(j,n) if v[i][j])
        v[i],v[j]=v[j],v[i];q=v[j][j];v[j]=[x/q for x in v[j]]
        for i in range(n):
            if i!=j and v[i][j]:
                q=v[i][j];v[i]=[x-q*y for x,y in zip(v[i],v[j])]
    return [row[-1] for row in v]

def optimum(G,h,c,ne,basis):
    p=len(c);assert len(basis)==len(set(basis))==p
    assert all(0<=i<len(G) for i in basis)
    A=[G[i] for i in basis];x=solve(A,[h[i] for i in basis])
    dual=solve(list(map(list,zip(*A))),c)
    residual=[dot(row,x)-rhs for row,rhs in zip(G,h)]
    assert all(v==0 for v in residual[:ne]) and all(v<=0 for v in residual[ne:])
    assert all(v>=0 for i,v in zip(basis,dual) if i>=ne)
    assert all(dot([G[i][j] for i in basis],dual)==c[j] for j in range(p))
    assert all(residual[i]*v==0 for i,v in zip(basis,dual))
    assert dot(c,x)==dot([h[i] for i in basis],dual)
    face=[i for i,v in zip(basis,dual) if i>=ne and v>0]
    return x,face

@lru_cache(maxsize=1)
def fixtures():
    return json.loads((Path(__file__).parent/'fixtures/twostage_bases.json').read_text())

@lru_cache(maxsize=None)
def noncooperative(name):
    pack=fixtures();case=pack['cases'][name];data=pack['data'][case['data']]
    a=[[F(float.fromhex(v)) for v in row] for row in data]
    m,d,s=case['shape'];n=len(a);p=m+d+s
    X=[row[1:1+m] for row in a];Z=[row[1+m:1+m+d] for row in a];Y=[row[1+m+d:] for row in a]
    zero=F(0);one=F(1)
    U=[[-v for v in x]+z+[zero]*s for x,z in zip(X,Z)]+[[zero]*m+[-v for v in z]+y for z,y in zip(Z,Y)]
    bounds=[[-F(i==j) for j in range(p)] for i in range(p)]
    answers=[]
    for k,witness in enumerate(case['rows']):
        ex=X[k]+[zero]*(d+s);ez=[zero]*m+Z[k]+[zero]*s;ey=[zero]*(m+d)+Y[k]
        c=[-v for v in ex] if case['leader']==1 else ey
        G=[ez]+U+bounds;h=[one]+[zero]*(2*n+p)
        primary,face=optimum(G,h,c,1,witness['primary'])
        assert all(h[i]==0 for i in face)
        E=[ez if case['leader']==1 else ex]+[G[i] for i in face]
        H=[one]+[zero]*len(face);secondary,_=optimum(E+U+bounds,H+[zero]*(2*n+p),ey if case['leader']==1 else ez,len(E),witness['secondary'])
        den,mid,num=[dot(q,secondary) for q in [ex,ez,ey]]
        assert den>0 and mid>0 and num>=0
        if case['leader']==1:assert mid==1 and den==dot(ex,primary)
        else:assert den==1 and num==dot(ey,primary)*mid
        values=[num/den,mid/den,num/mid];assert all(0<=v<=1 for v in values)
        answers.append(values)
    return a,answers

def check(name,data):
    stored,values=noncooperative(name)
    assert stored==[[F(float(v)) for v in row] for row in data], 'Witness input differs from executed data'
    return values
