#!/usr/bin/env python3
"""Verify precision-recovery arithmetic with exact stored-double oracles.

Requires NumPy/SciPy and Stata16+. All inputs are deterministic synthetic values;
no external files or private research directories are required.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,json,subprocess,itertools,argparse,hashlib,shutil
import numpy as np
from test_core_retries import WRITER,readbin
from test_simplex_stability import mata
REPO=Path(__file__).resolve().parents[1]

def check_arithmetic(root,pkg,stata):
    out=root/'arithmetic';out.mkdir();rng=np.random.default_rng(202611038)
    body='version 16\nclear all\nset more off\nset processors 1\nset linesize 160\n'
    for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:(out/key).mkdir();body+=f'sysdir set {key} "{out/key}"\n'
    body+=WRITER+f'adopath ++ "{pkg}"\nrun "{pkg}/_dea_simplex.ado"\nmata:\n'
    def number(x):return format(float(x),'.17g')
    def mat(label,A):
     nonlocal body
     A=np.array(A,dtype=float);body+=f'{label}=J({A.shape[0]},{A.shape[1]},0)\n'
     for i,row in enumerate(A):
      for k in range(0,len(row),64):
       last=min(k+64,len(row));body+=f'{label}[{i+1},{k+1}..{last}]={mata(row[None,k:last])}\n'
    products=[];adds=[];dots=[]
    for i in range(600):
     a=float(rng.uniform(-2,2)*2.**rng.integers(-390,390));b=float(rng.uniform(-2,2)*2.**rng.integers(-390,390))
     if i%10==0:b=-a
     products.append((a,b));body+=f'savebin("p{i:04}.bin",_rsm_q_product({number(a)},{number(b)}))\n'
    for i in range(800):
     a=float(rng.uniform(-2,2)*2.**rng.integers(-800,800));b=float(rng.uniform(-2,2)*2.**rng.integers(-800,800))
     if i%4==0:b=-float(np.nextafter(a,np.inf))
     if i%8==0:b=-a
     adds.append((a,b));body+=f'savebin("a{i:04}.bin",(_rsm_q_add({number(a)},{number(b)},0),_rsm_q_add({number(a)},{number(b)},1)))\n'
    for n,kind,rep in itertools.product([0,1,2,3,8,32,128,512,4096],['random','cancel','flat'],range(4)):
     if kind=='random':
      a=rng.uniform(-2,2,n)*2.**rng.integers(-390,390,n);b=rng.uniform(-2,2,n)*2.**rng.integers(-390,390,n)
     elif kind=='flat':
      a=rng.uniform(-2,2,n)*2.**(rep*80-120);b=rng.uniform(-2,2,n)*2.**(120-rep*80)
     else:
      a=np.zeros(n);b=np.zeros(n)
      for i in range(0,n-1,2):
       a[i]=float(rng.uniform(-2,2)*2.**rng.integers(-390,390));a[i+1]=-a[i];b[i]=b[i+1]=float(rng.uniform(-2,2)*2.**rng.integers(-390,390))
      if n%2:a[-1]=1.;b[-1]=2.**(-100-rep)
     D=0. if kind=='cancel' else float(rng.uniform(-2,2)*2.**rng.integers(-200,200))
     dots.append((a,b,D,kind));i=len(dots)-1;mat('a',a[None,:]);mat('b',b[:,None]);body+=f'savebin("d{i:03}.bin",_rsm_q_dot(a,b,{number(D)}))\n'
    # Product envelope violations and missing/dimension/length rejection controls.
    body+='assert(hasmissing(_rsm_q_product(1e-130,1)))\nassert(hasmissing(_rsm_q_product(1e130,1)))\nassert(hasmissing(_rsm_q_product(.,0)))\nassert(hasmissing(_rsm_q_dot((1,2),1,0)))\nassert(hasmissing(_rsm_q_dot(J(1,4097,1),J(4097,1,1),0)))\nassert(hasmissing(_rsm_q_dot(1,1,.)))\nassert(hasmissing(_rsm_q_add(1e300,0,1)))\nassert(hasmissing(_rsm_q_add(1e-300,0,1)))\nassert(hasmissing(_rsm_q_mulup(-1,1)))\n'
    body+='end\ndi "ARITHMETIC_PASS"\nexit,clear\n';(out/'test.do').write_text(body)
    p=subprocess.run([stata,'-b','do','test.do'],cwd=out,capture_output=True,timeout=600)
    (out/'process.json').write_text(json.dumps(dict(returncode=p.returncode,stderr=p.stderr.decode(errors='replace')),indent=2)+'\n')
    log=(out/'test.log').read_text();assert p.returncode==0 and 'ARITHMETIC_PASS' in log.splitlines(),log[-3000:]
    for i,(a,b) in enumerate(products):
     q=readbin(out/f'p{i:04}.bin').ravel();assert np.isfinite(q).all() and F(a)*F(b)==F(q[0])+F(q[1]),('product',i)
    for i,(a,b) in enumerate(adds):
     q=readbin(out/f'a{i:04}.bin').ravel();exact=F(a)+F(b);assert np.isfinite(q).all() and F(q[0])<=exact<=F(q[1]),('add',i)
    rows=[]
    for i,(a,b,D,kind) in enumerate(dots):
     q=readbin(out/f'd{i:03}.bin').ravel();exact=F(D)-sum(F(x)*F(y) for x,y in zip(a,b));assert np.isfinite(q).all() and F(q[1])>=0 and abs(F(q[0])-exact)<=F(q[1]),('dot',i)
     rows.append(dict(case=i,n=len(a),kind=kind,center=float(q[0]),radius=float(q[1]),exact_value=str(exact),enclosed=True))
    summary=dict(exact_products=len(products),outward_additions=len(adds),bounded_dot_products=len(dots),largest_dot_dimension=4096,rejection_controls=9,all_passed=True)
    (root/'ARITHMETIC_CHECKS.json').write_text(json.dumps(dict(summary=summary,dots=rows),indent=2)+'\n');print(json.dumps(summary,indent=2))
    return summary

def check_edges(root,pkg,stata):
    out=root/'edges';out.mkdir();lo=2.**-400;hi=2.**400
    values=[lo,np.nextafter(lo,0),np.nextafter(lo,np.inf),hi,np.nextafter(hi,0),np.nextafter(hi,np.inf)]
    products=[(sign*a,b) for a,b,sign in itertools.product(values,[lo,hi,1.,0.],[-1,1])]
    adds=[(sign*a,sign*b) for a,b,sign in itertools.product([2.**-970,np.nextafter(2.**-970,0),2.**900,np.nextafter(2.**900,np.inf)],[0.,2.**-970],[-1,1])]
    body='version 16\nclear all\nset more off\nset processors 1\n'
    for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:(out/key).mkdir();body+=f'sysdir set {key} "{out/key}"\n'
    body+=WRITER+f'adopath ++ "{pkg}"\nrun "{pkg}/_dea_simplex.ado"\nmata:\n'
    for i,(a,b) in enumerate(products):body+=f'savebin("p{i:03}.bin",_rsm_q_product({a:.17g},{b:.17g}))\n'
    for i,(a,b) in enumerate(adds):body+=f'savebin("a{i:03}.bin",(_rsm_q_add({a:.17g},{b:.17g},0),_rsm_q_add({a:.17g},{b:.17g},1)))\n'
    for i,exp in enumerate([-400,400]):
     v=2.**exp
     body+=f'a=J(1,4096,{v:.17g});b=J(4096,1,{v:.17g})\nsavebin("d{i}_positive.bin",_rsm_q_dot(a,b,0))\na[1,2049..4096]=-a[1,2049..4096]\nsavebin("d{i}_cancel.bin",_rsm_q_dot(a,b,0))\n'
    body+='end\ndi "EDGE_PASS"\nexit,clear\n';(out/'test.do').write_text(body)
    p=subprocess.run([stata,'-b','do','test.do'],cwd=out,capture_output=True,timeout=180)
    (out/'process.json').write_text(json.dumps(dict(returncode=p.returncode,stderr=p.stderr.decode(errors='replace')),indent=2)+'\n');assert p.returncode==0 and 'EDGE_PASS' in (out/'test.log').read_text().splitlines()
    rejected=0
    for i,(a,b) in enumerate(products):
     q=readbin(out/f'p{i:03}.bin').ravel();expected=a==0 or b==0 or (lo<=abs(a)<=hi and lo<=abs(b)<=hi)
     assert np.isfinite(q).all()==expected,(i,q,expected)
     if expected:assert F(a)*F(b)==F(q[0])+F(q[1])
     else:rejected+=1
    for i,(a,b) in enumerate(adds):
     q=readbin(out/f'a{i:03}.bin').ravel();expected=all(z==0 or 2.**-970<=abs(z)<=2.**900 for z in [a,b]);assert np.isfinite(q).all()==expected
     if expected:assert F(q[0])<=F(a)+F(b)<=F(q[1])
     else:rejected+=1
    for i,exp in enumerate([-400,400]):
     for kind in ['positive','cancel']:
      exact=-4096*F(2.**exp)**2 if kind=='positive' else F(0);q=readbin(out/f'd{i}_{kind}.bin').ravel();assert np.isfinite(q).all() and abs(F(q[0])-exact)<=F(q[1])
    summary=dict(product_edges=len(products),addition_edges=len(adds),dot_edges=4,largest_dot_dimension=4096,expected_range_rejections=rejected,all_passed=True)
    (root/'EDGE_CHECKS.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
    return summary

def check_quality(root,pkg,stata):
    out=root/'quality';out.mkdir();rng=np.random.default_rng(202611039);cases=[]
    def add(name,A,b,c,x,h,l,tol):
     cases.append(dict(name=name,A=np.array(A).reshape(len(b),len(x)).tolist(),b=list(b),c=list(c),x=list(x),yh=list(h),yl=list(l),tol=tol))
    for kind,exp,tol,shift in itertools.product(['gap','dual','complementarity','weighted_primal'],[-80,0,80],[1e-8,1e-10,1e-12],range(-3,4)):
     s=2.**exp
     if kind=='gap':
      t=tol/(1-tol)
     elif kind=='dual':t=tol*min(1.,s)
     elif kind=='complementarity':t=2*tol
     else:t=tol/2
     for k in range(abs(shift)):t=float(np.nextafter(t,0. if shift<0 else np.inf))
     name=f'{kind}_{exp}_{tol}_{shift}'
     if kind=='gap':add(name,[[s]],[0.],[1.],[t],[1/s],[0.],tol)
     elif kind=='dual':add(name,[[1.]],[0.],[s],[0.],[s],[t],tol)
     elif kind=='complementarity':add(name,[[s]],[s],[1.],[1.],[1/s],[-t/s],tol)
     else:add(name,[[s],[s]],[0.,0.],[0.],[t],[1/s,-1/s],[0.,0.],tol)
    for m,rep,exp in itertools.product([1,2,3,8,16,32,64],range(4),[-40,0,40]):
     n=m+4;A=rng.integers(-4,5,size=(m,n)).astype(float);x=np.r_[rng.integers(1,5,m)/4,np.zeros(4)]
     rs=2.**rng.integers(-8,9,m);cs=2.**rng.integers(-8,9,n);A=A*rs[:,None]*cs;x=x/cs
     y=rng.integers(-5,6,m)/rs*2.**exp;b=A@x;c=y@A;c[-4:]+=np.arange(1,5)*cs[-4:]*2.**exp
     if rep==1 and m>=2:
      # Exact dyadic dual with zero objective; inequalities encoded by nonbasic cost.
      y=np.zeros(m);y[0]=b[1];y[1]=-b[0];c=y@A;c[-4:]+=np.arange(1,5)*cs[-4:]
     if rep==2:c[0]+=2.**-40
     if rep==3:x[0]=np.nextafter(x[0],np.inf)
     for tol in [1e-8,1e-12]:add(f'dense_{m}_{rep}_{exp}_{tol}',A,b,c,x,y,np.zeros(m),tol)
    body='version 16\nclear all\nset more off\nset processors 1\nset linesize 160\n'
    for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:(out/key).mkdir();body+=f'sysdir set {key} "{out/key}"\n'
    body+=WRITER+f'adopath ++ "{pkg}"\nrun "{pkg}/_dea_simplex.ado"\nmata:\n'
    body+=r'''
    void quality_probe(string scalar name,real matrix A,real colvector b,real rowvector c,real colvector x,real rowvector h,real rowvector l,real scalar tol) {
     struct RSM_Prob scalar P
     P.A=A;P.b=b;P.c=c;P.n_orig=cols(A);P.m_con=rows(A)
     savebin(name+"_q.bin",_rsm_retry_quality(P,x,h,l,tol))
     savebin(name+"_A.bin",A);savebin(name+"_b.bin",b);savebin(name+"_c.bin",c)
     savebin(name+"_x.bin",x);savebin(name+"_h.bin",h);savebin(name+"_l.bin",l)
    }
    '''
    for i,case in enumerate(cases):
     for label in ['A','b','c','x','yh','yl']:
      A=np.array(case[label]);A=A.reshape(-1,1) if label in ['b','x'] else np.atleast_2d(A)
      body+=f'{label}=J({A.shape[0]},{A.shape[1]},0)\n'
      for k,row in enumerate(A):
       for col in range(0,len(row),32):
        last=min(col+32,len(row));body+=f'{label}[{k+1},{col+1}..{last}]={mata(row[None,col:last])}\n'
     body+=f'quality_probe("q{i:03}",A,b,c,x,yh,yl,{case["tol"]:.17g})\n'
    body+='end\ndi "QUALITY_PASS"\nexit,clear\n';(out/'test.do').write_text(body)
    (out/'cases.json').write_text(json.dumps(cases,indent=2)+'\n')
    p=subprocess.run([stata,'-b','do','test.do'],cwd=out,capture_output=True,timeout=600)
    (out/'process.json').write_text(json.dumps(dict(returncode=p.returncode,stderr=p.stderr.decode(errors='replace')),indent=2)+'\n');log=(out/'test.log').read_text();assert p.returncode==0 and 'QUALITY_PASS' in log.splitlines(),log[-2000:]
    rows=[];dot=lambda a,b:sum(x*y for x,y in zip(a,b))
    for i,case in enumerate(cases):
     for key,label in [('A','A'),('b','b'),('c','c'),('x','x'),('yh','h'),('yl','l')]:assert np.array_equal(np.array(case[key]).ravel(),readbin(out/f'q{i:03}_{label}.bin').ravel())
     A=[[F(z) for z in row] for row in case['A']];b,c,x,h,l=[list(map(F,case[k])) for k in ['b','c','x','yh','yl']];y=[v+w for v,w in zip(h,l)];tol=F(case['tol']);cutoff=min(F(1e-8),tol)
     obj=dot(c,x);dobj=dot(y,b);rp=[rhs-dot(row,x) for row,rhs in zip(A,b)];rd=[ci-dot(y,col) for ci,col in zip(c,zip(*A))]
     nums=[abs(obj-dobj),dot([abs(v)+abs(w) for v,w in zip(h,l)],list(map(abs,rp))),dot(list(map(abs,rd)),list(map(abs,x)))];den=1+abs(obj);du=max(F(0),-min(rd));threshold=cutoff*min(F(1),max(map(abs,c)))
     expected=all(v<=cutoff*den for v in nums) and du<=threshold
     q=readbin(out/f'q{i:03}_q.bin').ravel();available=np.isfinite(q).all();accept=bool(q[0]==1)
     assert not accept or expected,(i,case['name'])
     if available:
      assert F(q[8])<=den and all(F(z)>=t for z,t in zip(q[9:12],nums)) and F(q[6])>=du,(i,'bounds')
      expr=[v-cutoff*den for v in nums]+[du-threshold];assert all(F(v)>=e for v,e in zip(q[12:16],expr)),(i,'comparison bound')
     rows.append(dict(case=i,name=case['name'],expected=expected,accepted=accept,available=bool(available),false_rejection=bool(expected and not accept),exact_objective=str(obj),enclosures_valid=bool(available)))
    summary=dict(cases=len(rows),boundary_cases=252,dense_cases=len(rows)-252,max_rows=64,max_cols=68,all_stored_data_verified=True,false_acceptances=0,false_rejections=sum(x['false_rejection'] for x in rows),unavailable=sum(not x['available'] for x in rows),exact_zero_objectives=sum(x['exact_objective']=='0' for x in rows),all_finite_bounds_enclose_exact=True)
    (root/'QUALITY_CHECKS.json').write_text(json.dumps(dict(summary=summary,results=rows),indent=2)+'\n');print(json.dumps(summary,indent=2))
    return summary

def check_helpers(root,pkg,stata):
    out=root/'helpers';out.mkdir()
    probes=[]
    for a in [3.,7.,11.,1.5]:probes.append((f'scalar_{a}',np.array([[a]]),np.ones((1,1))))
    for k,scale,transpose in itertools.product([25,28,30,35,40],['base','row','col','both'],[False,True]):
        d=2.**-k;B=np.array([[1.,2,4],[1,2+d,4-d],[1+d,2,4-d]])
        if scale in ['row','both']:B*=np.array([2.**-20,1,2.**20])[:,None]
        if scale in ['col','both']:B*=np.array([2.**12,2.**-12,1])[None,:]
        if transpose:B=B.T
        D=np.array([[1.],[2.],[3.]]) if transpose else np.full((3,1),7.)
        probes.append((f'near_{k}_{scale}_t{int(transpose)}',B,D))
    rng=np.random.default_rng(202611037)
    for n,rep,scale in itertools.product([2,3,5,10],range(5),['base','global','mixed']):
        B=rng.integers(-3,4,size=(n,n)).astype(float)+np.eye(n)*(5*n);D=rng.integers(-5,6,size=(n,1)).astype(float)
        if scale=='global':B*=2.**(-80+40*rep);D*=2.**(-80+40*rep)
        if scale=='mixed':B*=2.**rng.integers(-20,21,size=(n,1));B*=2.**rng.integers(-10,11,size=(1,n))
        probes.append((f'dense_{n}_{rep}_{scale}',B,D))
    body='version 16\nclear all\nset more off\nset processors 1\n'
    for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:(out/key).mkdir();body+=f'sysdir set {key} "{out/key}"\n'
    body+=WRITER+f'adopath ++ "{pkg}"\nrun "{pkg}/_dea_simplex.ado"\nmata:\n'
    for i,(name,B,D) in enumerate(probes):
        for label,A in [('B',B),('D',D)]:body+=f'{label}=J({A.shape[0]},{A.shape[1]},0)\n'+''.join(f'{label}[{j+1},.]={mata(row[None,:])}\n' for j,row in enumerate(A))
        body+=f'savebin("h{i:03}_mat.bin",B)\nsavebin("h{i:03}_rhs.bin",D)\nsavebin("h{i:03}_single.bin",_rsm_factor_solve(B,D))\nsavebin("h{i:03}_pair.bin",_rsm_retry_refine(B,D))\n'
    body+=r'''
    assert(hasmissing(_rsm_retry_refine((1,2\2,4),(1\2))))
    assert(hasmissing(_rsm_retry_refine(J(2,3,1),(1\2))))
    assert(hasmissing(_rsm_retry_refine(I(2),(1\2\3))))
    assert(hasmissing(_rsm_retry_refine((.,0\0,1),(1\2))))
    assert(hasmissing(_rsm_retry_refine(I(2),(1\.))))
    void quality_controls() {
        struct RSM_Prob scalar P
        real rowvector q
        _rsm_std_form(P,J(0,3,0),J(0,1,0),(1,0,1\0,1,1),(1\2),(3,4,8))
        q=_rsm_retry_quality(P,(1\2\0),(3,4),(0,0),1e-10);assert(q[1]==1)
        P.c=(3,4,6)
        q=_rsm_retry_quality(P,(1\2\0),(3,4),(0,0),1e-10);assert(q[1]==0)
        P.c=(3,4,8)
        q=_rsm_retry_quality(P,(1.01\2\0),(3,4),(0,0),1e-10);assert(q[1]==0)
        q=_rsm_retry_quality(P,(1\2\0),(3,4),(.,0),1e-10);assert(hasmissing(q))
        q=_rsm_retry_quality(P,(1\2\0),(3,4),(0,0,0),1e-10);assert(hasmissing(q))
    }
    quality_controls()
    void cutoff_probe(string scalar name,real scalar x,real scalar tol) {
        struct RSM_Prob scalar P
        _rsm_std_form(P,J(0,1,0),J(0,1,0),1,0,1)
        savebin(name+"_quality.bin",_rsm_retry_quality(P,x,1,0,tol))
        savebin(name+"_input.bin",(x,tol))
    }
    '''
    cutoffs=[]
    for tol,shift in itertools.product([1e-8,1e-10,1e-12],range(-3,4)):
        x=tol/(1-tol)
        for _ in range(abs(shift)):x=np.nextafter(x,0. if shift<0 else np.inf)
        cutoffs.append(dict(x=float(x),tol=tol,shift=shift));i=len(cutoffs)-1;body+=f'cutoff_probe("b{i:03}",{x:.17g},{tol:.17g})\n'
    for i,exp in enumerate([-900,-500,500,900]):
        body+=f'B={3*2.**exp:.17g}\nD={2.**exp:.17g}\nsavebin("extreme_{i}.bin",_rsm_retry_refine(B,D))\n'
    body+='end\ndi "HELPER_PASS"\nexit,clear\n';(out/'test.do').write_text(body)
    p=subprocess.run([stata,'-b','do','test.do'],cwd=out,capture_output=True,timeout=240)
    (out/'process.json').write_text(json.dumps(dict(returncode=p.returncode,stderr=p.stderr.decode(errors='replace')),indent=2)+'\n');log=(out/'test.log').read_text();assert p.returncode==0 and 'HELPER_PASS' in log.splitlines(),log[-2500:]
    def solve(B,D):
        T=[[F(float(v)) for v in row]+[F(float(rhs))] for row,rhs in zip(B,D.ravel())];n=len(T)
        for j in range(n):
            k=next(k for k in range(j,n) if T[k][j]);T[j],T[k]=T[k],T[j];d=T[j][j];T[j]=[v/d for v in T[j]]
            for i in range(n):
                if i!=j:d=T[i][j];T[i]=[a-d*b for a,b in zip(T[i],T[j])]
        return [row[-1] for row in T]
    rows=[]
    for i,(name,B,D) in enumerate(probes):
        assert np.array_equal(B,readbin(out/f'h{i:03}_mat.bin')) and np.array_equal(D,readbin(out/f'h{i:03}_rhs.bin'))
        X=readbin(out/f'h{i:03}_pair.bin');single=readbin(out/f'h{i:03}_single.bin').ravel();assert np.isfinite(X).all()
        exact=solve(B,D);pair=[F(float(h))+F(float(l)) for h,l in X]
        old=max(abs(F(float(s))-t)/(1+abs(t)) for s,t in zip(single,exact));new=max(abs(s-t)/(1+abs(t)) for s,t in zip(pair,exact))
        rows.append(dict(name=name,single_forward_error=float(old),twofold_forward_error=float(new),exact_twofold=pair==exact,improved=new<old))
        assert new<=max(old,F(1e-30)),(name,float(old),float(new))
    decisions=[]
    for i,c in enumerate(cutoffs):
        x,tol=readbin(out/f'b{i:03}_input.bin').ravel();q=readbin(out/f'b{i:03}_quality.bin').ravel();assert x==c['x'] and tol==c['tol']
        metric=F(float(x))/(1+F(float(x)));expected=metric<=F(float(tol));accepted=q[0]==1
        decisions.append(dict(**c,exact_metric=str(metric),expected=expected,accepted=bool(accepted),false_accept=bool(accepted and not expected),false_reject=bool(expected and not accepted)))
    summary=dict(linear_systems=len(rows),improved=sum(x['improved'] for x in rows),max_twofold_forward_error=max(x['twofold_forward_error'] for x in rows),malformed_singular_controls=5,quality_controls=5,boundary_controls=len(decisions),boundary_false_acceptances=sum(x['false_accept'] for x in decisions),boundary_false_rejections=sum(x['false_reject'] for x in decisions),extreme_scalars_unavailable=sum(not np.isfinite(readbin(out/f'extreme_{i}.bin')).all() for i in range(4)))
    (root/'HELPER_CHECKS.json').write_text(json.dumps(dict(summary=summary,systems=rows,boundaries=decisions),indent=2)+'\n');print(json.dumps(summary,indent=2))
    return summary

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();root=args.output.resolve();root.mkdir(parents=True,exist_ok=False)
    pkg=root/'ado';shutil.copytree(REPO/'ado',pkg)
    manifest=dict(source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in pkg.iterdir() if p.is_file()},numpy=np.__version__,stata=args.stata)
    (root/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    results={}
    for name,fn in [('arithmetic',check_arithmetic),('edges',check_edges),('quality',check_quality),('helpers',check_helpers)]:
        results[name]=fn(root,pkg,args.stata)
    assert results['quality']['false_acceptances']==results['quality']['false_rejections']==results['quality']['unavailable']==0
    assert results['helpers']['boundary_false_acceptances']==results['helpers']['boundary_false_rejections']==0
    summary=dict(all_passed=True,stages=results,boundary_sets_overlap=True)
    (root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
