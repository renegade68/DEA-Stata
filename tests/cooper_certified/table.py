"""Exact checks of packaged Cooper source-example witnesses, without an LP solver.

The three data representations are distinct models. Reconstruct their equations
and exact stage faces before checking sparse primal/dual witnesses. Numeric
proposals and returned doubles do not authorize these proofs.
"""
from fractions import Fraction as F
import copy

REPRESENTATIONS=('decimal','binary64','normalized_binary64')
STEPS=('radial','slack','projection_0_min','projection_0_max',
       'projection_1_min','projection_1_max','projection_2_min','projection_2_max',
       'decomposition','target_mu_min','target_mu_max')

def dot(a,b):
    assert len(a)==len(b)
    return sum((x*y for x,y in zip(a,b)),F(0))

def vector(sparse,size):
    assert isinstance(sparse,dict)
    result=[F(0)]*size
    for key,value in sparse.items():
        index=int(key);assert str(index)==key and 0<=index<size
        result[index]=F(value)
    return result

def certify(A,b,c,proof):
    assert len(A)==len(b) and all(len(row)==len(c) for row in A)
    x=vector(proof['x'],len(c));y=vector(proof['y'],len(b));value=F(proof['objective'])
    assert all(v>=0 for v in x)
    assert [dot(row,x) for row in A]==b
    assert all(dot(col,y)<=v for col,v in zip(zip(*A),c))
    assert dot(c,x)==dot(b,y)==value
    return x,value

def case_data(fixture,case,representation):
    ds=fixture['datasets'][case['industry']]
    assert ds['years']==list(range(1981,1998))
    assert ds['columns']==['labor','capital','output'] and len(ds['data'])==17
    data=[[F(x) for x in row] for row in ds['data']]
    assert all(len(row)==3 and min(row)>0 for row in data)
    units=[F(1)]*3;weights=[F(1)]*3
    if representation=='binary64':data=[[F(float(x)) for x in row] for row in data]
    elif representation=='normalized_binary64':
        maxima=[max(float(row[h]) for row in data) for h in range(3)]
        data=[[F(float(row[h])/maxima[h]) for h in range(3)] for row in data]
        units=[F(x) for x in maxima]
        weights=[F(float.fromhex(x)) for x in case['normalized_slack_weights']]
        assert weights==[F(x/max(maxima)) for x in maxima]
    else:assert representation=='decimal'
    return data,units,weights

def evaluate(fixture,mutations=True):
    assert fixture['schema']==1 and len(fixture['cases'])==4
    assert set(fixture['datasets'])=={'textile','automobile'}
    assert {(c['industry'],c['year'],c['input']) for c in fixture['cases']}=={
        ('textile',1989,'capital'),('textile',1993,'labor'),
        ('automobile',1990,'labor'),('automobile',1996,'labor')}
    models=[];records=[];rejected=0
    def check(tag,A,b,c,proof):
        nonlocal rejected
        x,value=certify(A,b,c,proof)
        if mutations:
            for field in ['x','y','objective','rhs']:
                bad=copy.deepcopy(proof);bb=b[:];eps=F(1,10**50)
                if field=='rhs':bb[0]+=eps
                elif field=='objective':bad[field]=str(F(bad[field])+eps)
                else:bad[field]['0']=str(F(bad[field].get('0','0'))+eps)
                try:certify(A,bb,c,bad)
                except AssertionError:rejected+=1
                else:raise AssertionError('Changed witness was accepted: '+field)
        models.append(dict(id=tag,A=A,b=b,c=c,objective=value))
        return x,value
    for case in fixture['cases']:
        assert case['id']==case['industry']+'_'+str(case['year'])
        k=case['year']-1981;j=['labor','capital'].index(case['input'])
        for representation in REPRESENTATIONS:
            data,units,weights=case_data(fixture,case,representation)
            proofs=case['representations'][representation]['proofs'];assert set(proofs)==set(STEPS)
            prefix=case['id']+'_'+representation+'_'
            def run(name,A,b,c):return check(prefix+name,A,b,c,proofs[name])
            # t=phi-1 >=0; omitted lambda/slack upper bounds follow from
            # convexity, positive quantities and the original equalities.
            E=[[F(0)]*21 for _ in range(4)]
            for h in range(3):
                for i in range(17):E[h][i]=data[i][h]
            E[0][17]=E[1][18]=1;E[2][19]=-1;E[2][20]=-data[k][2]
            E[3][:17]=[F(1)]*17;b=data[k]+[F(1)]
            radial,rf=run('radial',E,b,[F(0)]*20+[F(-1)])
            trow=[F(0)]*20+[F(1)];As=E+[trow];bs=b+[radial[-1]]
            cs=[F(0)]*17+[-w for w in weights]+[F(0)]
            slack,sf=run('slack',As,bs,cs);assert slack[-1]==radial[-1]
            Af=As+[cs];bf=bs+[sf];projection=[]
            for h,col in enumerate(zip(*data)):
                q=list(col)+[F(0)]*4
                _,lo=run(f'projection_{h}_min',Af,bf,q)
                _,negative_hi=run(f'projection_{h}_max',Af,bf,[-v for v in q])
                assert lo==-negative_hi==dot(list(col),slack[:17])
                projection.append(lo)
            D=[[F(0)]*21 for _ in range(6)]
            for h in range(2):
                for i in range(17):D[h][i]=D[h+2][i]=data[i][h]
                D[h][17+h]=1;D[h+2][19+h]=-1
            for i in range(17):D[4][i]=data[i][2];D[5][i]=1
            bd=data[k][:2]+projection[:2]+[projection[2],F(1)]
            ct=([F(float.fromhex(x)) for x in case['normalized_technical_objective']]
                if representation=='normalized_binary64' else [row[0]+row[1] for row in data])
            assert len(ct)==17 and min(ct)>0
            cd=[-x for x in ct]+[F(0)]*4
            mu,df=run('decomposition',D,bd,cd);Dm=D+[cd];bm=bd+[df]
            q=[row[j] for row in data]+[F(0)]*4
            _,lo=run('target_mu_min',Dm,bm,q)
            _,negative_hi=run('target_mu_max',Dm,bm,[-x for x in q])
            assert lo==-negative_hi
            component=(data[k][j]-lo)*units[j]
            assert component==F(case['representations'][representation]['component'])
            table=F(fixture['datasets'][case['industry']]['manuscript_table'][k][j])
            assert component<max(F(0),table-F('0.005')) or component>table+F('0.005')
            records.append(dict(id=case['id'],representation=representation,component=str(component),
                projection=[str(v*u) for v,u in zip(projection,units)],phi=str(1+radial[-1]),
                projection_unique=True,component_unique=True,printed_interval_disjoint=True))
    assert len(models)==132 and len(records)==12
    assert rejected==(528 if mutations else 0)
    return records,models,rejected
