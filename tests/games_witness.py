"""Original binary64 input certificates and independent exact game formulas.

Basis row indices are hints only: Fraction elimination proves primal/dual
feasibility, complementarity and strong duality before using an efficiency.
Rounded exported scores are never inputs to the exact derived oracle.
"""
from fractions import Fraction as F
from functools import lru_cache
from decimal import Decimal, localcontext
from pathlib import Path
import json
import math
import re
import numpy as np
from twostage_witness import dot, solve, optimum
from test_core_retries import readbin

@lru_cache(maxsize=1)
def fixtures():
    return json.loads((Path(__file__).parent/'fixtures/games_bases.json').read_text())

@lru_cache(maxsize=None)
def certified(name):
    pack=fixtures();case=pack['cases'][name]
    a=[[F(float.fromhex(v)) for v in row] for row in pack['data'][case['data']]]
    m,s=case['shape'];n=len(a);vrs=case['rts']=='vrs';p=m+s+2*vrs
    X=[row[1:1+m] for row in a];Y=[row[1+m:] for row in a]
    U=[[-v for v in x]+y+([F(1),F(-1)] if vrs else []) for x,y in zip(X,Y)]
    bounds=[[-F(i==j) for j in range(p)] for i in range(p)];h=[F(1)]+[F(0)]*(n+p);answers=[]
    for k,basis in enumerate(case['basis']):
        E=X[k]+[F(0)]*(s+2*vrs);c=[F(0)]*m+Y[k]+([F(1),F(-1)] if vrs else [])
        G=[E]+U+bounds;x,_=optimum(G,h,c,1,basis);score=dot(c,x)
        dual=solve(list(map(list,zip(*[G[i] for i in basis]))),c)
        lam=[F(0)]*n
        for i,v in zip(basis,dual):
            if 1<=i<=n:lam[i-1]=v
        assert 0<score<=1 and all(v>=0 for v in lam)
        if vrs:assert sum(lam)==1
        assert all(dot([r[j] for r in X],lam)<=X[k][j]*score for j in range(m))
        assert all(dot([r[j] for r in Y],lam)>=Y[k][j] for j in range(s))
        answers.append(score)
    return a,answers

def check_input(name,data):
    stored,answers=certified(name)
    assert stored==[[F(float(v)) for v in row] for row in data], 'Certificate data differ from executed inputs'
    return answers

def export(q,tol):
    try:v=float(q)
    except OverflowError:return 3,1102,None
    if not math.isfinite(v):return 3,1102,None
    if q and not v:return 3,1101,None
    err=abs(F(v)-q)/abs(q) if q else F(0)
    return (3,1103,None) if err>F(tol) else (2,0,v)

def allocation(d,tol):
    c=[F(1) if v>=1-F(tol) else v/(F(4,5)+v/5) for v in d]
    total=sum((a-b for a,b in zip(c,d)),F(0));w=[max(F(1,100),1-v) for v in d]
    shares=[total*v/sum(w) for v in w];b=[min(F(1),v+q) for v,q in zip(d,shares)]
    g=[a-v for a,v in zip(b,d)];allocated=sum(g,F(0));product=math.prod(g)
    zeros=sum(v==0 for v in g);caps=sum(v+q>1 for v,q in zip(d,shares));log=None
    if not zeros:
        with localcontext() as ctx:
            ctx.prec=90
            log=float(sum((ctx.ln(Decimal(v.numerator))-ctx.ln(Decimal(v.denominator)) for v in g),Decimal(0)))
    return dict(rows=list(zip(d,c,b,g)),scalars=[total,allocated,total-allocated,allocated/len(d),product],zeros=zeros,caps=caps,log=log)

def check_derived(d,truth,mode,tol,detail=True):
    R=readbin(d/'results.bin');Q=readbin(d/'quality.bin');S=readbin(d/'reason.bin');A=readbin(d/'aggregate.bin');M=readbin(d/'meta.bin').ravel();checked=0;withheld=0;coalitions=0
    def one(value,quality,reason,expected):
        nonlocal checked,withheld
        checked+=1;want=export(expected,tol)
        assert (quality,reason)==want[:2],(d.name,quality,reason,want)
        if want[0]==3:assert np.isnan(value);withheld+=1
        else:assert value==want[2],(d.name,value,want)
    assert Q.shape==S.shape==R.shape
    if mode=='cooperative':
        for i,v in enumerate(truth):
            for j,q in enumerate([v,v,F(0),F(1+sum(w>v for w in truth))]):one(R[i,j],Q[i,j],S[i,j],q)
        one(A[0,2],A[0,0],A[0,1],sum(truth,F(0)))
        assert M[6]==A[0,2] or np.isnan(M[6]) and np.isnan(A[0,2])
        if detail and len(truth)<=15:
            masks={int(m.group(1)):(float(m.group(2)),int(m.group(3))) for line in (d/'test.log').read_text().splitlines() if (m:=re.match(r'^\s*(\d+)\s*\|\s*([\d.eE+-]+)\s*\|\s*(\d)\s*$',line))}
            assert len(masks)==2**len(truth)
            for mask,(v,q) in masks.items():
                want=sum((v for i,v in enumerate(truth) if mask&(1<<i)),F(0));assert q==2
                assert abs(v-float(want))<=max(1,abs(float(want)))*6e-15
            coalitions=len(masks)
    else:
        want=allocation(truth,tol)
        for i,row in enumerate(want['rows']):
            # Synthetic strict-export rule tests intentionally inject certified
            # scores; they test only columns 2..4. Public score tests check col 1.
            for j in range(1,4):one(R[i,j],Q[i,j],S[i,j],row[j])
        for j,expected in enumerate(want['scalars']):one(A[j,2],A[j,0],A[j,1],expected)
        if want['zeros']:assert A[5,0]==0 and A[5,1]==1004 and np.isnan(A[5,2])
        else:assert A[5,0]==1 and A[5,1]==1005 and abs(A[5,2]-want['log'])<=2e-12*max(1,abs(want['log']))
        assert M[12]==want['zeros'] and M[13]==want['caps'] and M[14]==(export(want['scalars'][4],tol)[1]==1101)
        assert M[15]==(np.isfinite(R).all() and np.isfinite(A[:4,2]).all())
        assert np.array_equal(M[6:12],A[:,2],equal_nan=True)
    return dict(exact_derived_fields=checked,withheld_export_fields=withheld,coalition_values=coalitions)

def residual(E,b,U,h,c,x,opt):
    eq=max((abs(dot(a,x)-v)/(max(abs(v),sum(abs(q*t) for q,t in zip(a,x))) or F(1)) for a,v in zip(E,b)),default=F(0))
    iq=max((max(F(0),dot(a,x)-v)/(max(abs(v),sum(abs(q*t) for q,t in zip(a,x))) or F(1)) for a,v in zip(U,h)),default=F(0))
    gap=abs(dot(c,x)-opt)/max(F(1),abs(opt))
    return max(eq,iq,gap)

def check_scores(d,truth,m,s,rts,tol):
    a=readbin(d/'data.bin');X=[[F(float(v)) for v in row[1:1+m]] for row in a];Y=[[F(float(v)) for v in row[1+m:]] for row in a];n=len(a);vrs=rts=='vrs'
    Q=readbin(d/'score_quality.bin');D=readbin(d/'lp_diagnostics.bin');N=readbin(d/'native_lp.bin');R=readbin(d/'results.bin');L=readbin(d/'lambda.bin');W=readbin(d/'weights.bin');C=readbin(d/'counts.bin').ravel();M=readbin(d/'meta.bin').ravel()
    assert Q.shape==(n,36) and np.all(Q[:,6]==2) and np.all(Q[:,7]==0)
    assert C[0]==n and C[1]==C[2]==0 and C[3]==sum(N[:,2]==0) and C[4]==sum(N[:,3])
    assert np.array_equal(D[:,0],D[:,1]) and np.array_equal(D[:,1],R[:,0]) and np.all(D[:,2]==1) and np.all(D[:,3]==0)
    assert np.array_equal(Q[:,5],N[:,1],equal_nan=True) and np.array_equal(Q[:,8],D[:,1])
    assert M[0]==M[1]==n and M[2]==M[3]==0 and M[4]==1 and M[5]==tol
    maxerr=F(0)
    for k,score in enumerate(truth):
        assert R[k,0]==float(score) and abs(F(float(score))-score)/score<=F(tol)
        lam=list(map(lambda x:F(float(x)),L[k]));w=list(map(lambda x:F(float(x)),W[k]))
        assert all(v>=0 for v in lam+w)
        E=[[F(1)]*n+[F(0)]] if vrs else [];h=[F(1)] if vrs else []
        env=residual(E,h,[[r[j] for r in X]+[-X[k][j]] for j in range(m)]+[[-r[j] for r in Y]+[F(0)] for j in range(s)],[F(0)]*m+[-v for v in Y[k]],[F(0)]*n+[F(1)],lam+[F(float(score))],score)
        mult=residual([X[k]+[F(0)]*(s+2*vrs)],[F(1)],[[-v for v in x]+y+([F(1),F(-1)] if vrs else []) for x,y in zip(X,Y)],[F(0)]*n,[F(0)]*m+Y[k]+([F(1),F(-1)] if vrs else []),w,score)
        assert max(env,mult)<=F(tol);maxerr=max(maxerr,env,mult)
    return dict(original_data_certificates=n,max_export_scaled_residual=float(maxerr),native_unavailable=int(C[3]))

COMMON=['n','n_lp','n_failed','n_boundary','complete','tolerance']
SPEC={'cooperative':['v_grand','detail_shown'],'nash':['total_gain','allocated_gain','unallocated_gain','mean_gain','nash_product','log_product','n_zero_gain','n_capped','product_underflow','allocation_complete']}
EXTRA=['n_certified','n_unverified','n_unavailable','n_native_failed','n_native_boundary','rank_quality','rank_reason','allocation_quality']

def capture_results(mode):
    body='matrix meta=('+','.join('r('+k+')' for k in COMMON+SPEC[mode])+')\nmata: savebin("meta.bin",st_matrix("meta"))\n'
    body+='matrix counts=('+','.join('r('+k+')' for k in EXTRA)+')\nmata: savebin("counts.bin",st_matrix("counts"))\n'
    for name,ret in dict(results='results',lp_diagnostics='lp_diagnostics',score_quality='score_diagnostics',native_lp='native_lp_diagnostics',**{'lambda':'lambda'},weights='weights',quality='result_quality',reason='result_reason',aggregate='aggregate_diagnostics',numerical='numerical_results',numaggregate='numerical_aggregates').items():
        body+=f'mata: savebin("{name}.bin",st_matrix("r({ret})"))\n'
    return body

def input_data(a,cols):
    return 'input double('+' '.join(cols)+')\n'+'\n'.join(' '.join(format(float(v),'.17g') for v in row) for row in a)+'\nend\n'
