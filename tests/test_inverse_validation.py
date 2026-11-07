#!/usr/bin/env python3
"""Check inverse identities and cutoff decisions against exact input-double arithmetic.

Requires NumPy/SciPy and Stata 16+. Synthetic/captured basis fixtures are distributed
with this test. No research-directory paths, downloads, or external R packages.
Failed logs and source hashes remain in the new output directory.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,shutil,subprocess
import numpy as np
from test_core_retries import WRITER,readbin
from test_simplex_stability import mata

REPO=Path(__file__).resolve().parents[1]
def decode(A):return np.array([[float.fromhex(v) if v is not None else np.nan for v in row] for row in A])
def assign(name,A):
    return f'{name}=J({A.shape[0]},{A.shape[1]},.)\n'+''.join(f'{name}[{i+1},.]={mata(row[None,:]).replace("nan",".")}\n' for i,row in enumerate(A))
def header(out,pkg):
    body='version 16\nclear all\nset more off\nset processors 1\n'
    for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:
        (out/key).mkdir();body+=f'sysdir set {key} "{out/key}"\n'
    return body+f'adopath ++ "{pkg}"\nrun "{pkg/"_dea_simplex.ado"}"\n'
def run(out,body,marker,stata):
    (out/'test.do').write_text(body)
    proc=subprocess.run([stata,'-b','do','test.do'],cwd=out,capture_output=True,timeout=240)
    (out/'process.json').write_text(json.dumps(dict(returncode=proc.returncode,stderr=proc.stderr.decode(errors='replace')),indent=2)+'\n')
    log=(out/'test.log').read_text(errors='replace')
    assert proc.returncode==0 and marker in log.splitlines(),log[-4000:]
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();root=args.output.resolve();root.mkdir(parents=True,exist_ok=False)
    pkg=root/'ado';shutil.copytree(REPO/'ado',pkg)
    path=REPO/'tests/fixtures/inverse_validation.json';fixture=json.loads(path.read_text())
    (root/'manifest.json').write_text(json.dumps(dict(source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in pkg.iterdir() if p.is_file()},fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),stata=args.stata,numpy=np.__version__),indent=2)+'\n')
    out=root/'guard';out.mkdir();tests=[]
    def add(name,B,X,expect=None):
        tests.append((name,np.asarray(B,dtype=float),np.asarray(X,dtype=float),expect))
    for row in fixture['pairs']:add(row['name'],decode(row['B']),decode(row['X']),row['expect_fast'])
    for boundary in [1e-8/2,1e-8]:
        for delta in [np.nextafter(boundary,0.),boundary,np.nextafter(boundary,np.inf)]:
            X=np.eye(3);X[0,1]=delta
            add('cutoff_'+delta.hex(),np.eye(3),X,int(delta<=1e-8/2))
    for exp in [500,900,1000]:
        B=np.array([[2.**exp,2.**exp],[2.**-exp,-2.**-exp]])
        X=np.array([[2.**(-exp-1),2.**(exp-1)],[2.**(-exp-1),-2.**(exp-1)]])
        add(f'wide_products_{exp}',B,X,0)
    for n in [1,2,3,10,25,64]:add(f'identity_{n}',np.eye(n),np.eye(n),1)
    for v in [np.nextafter(2.**-200,0.),2.**-200,np.nextafter(2.**-200,np.inf),np.nextafter(2.**200,0.),2.**200,np.nextafter(2.**200,np.inf)]:
        add('range_scalar_'+v.hex(),[[v]],[[1/v]],int(2.**-200<=v<=2.**200 and 2.**-200<=1/v<=2.**200))
    for v in [np.nextafter(2.**-200,0.),2.**-200,np.nextafter(2.**-200,np.inf),np.nextafter(0.,1.)]:
        add('range_offdiag_'+v.hex(),[[1,v],[0,1]],[[1,-v],[0,1]],int(v>=2.**-200))
    edge=np.sqrt((1e-8/4)/(64*np.finfo(float).eps*4)-1)
    for v in [128.,256.,np.nextafter(edge,0.),edge,np.nextafter(edge,np.inf)]:
        add('radius_'+v.hex(),[[1,v],[0,1]],[[1,-v],[0,1]],1 if v==128. else 0 if v==256. else None)
    rng=np.random.default_rng(202611035)
    for n in [2,3,5,10,20]:
        for rep in range(5):
            base=rng.integers(-3,4,size=(n,n)).astype(float)+np.eye(n)*(4*n)
            for scaling in ['none','rows','global']:
                B=base.copy()
                if scaling=='rows':B*=2.**rng.integers(-4,5,size=(n,1))
                if scaling=='global':B*=2.**([-180,-100,0,100,180][rep])
                X=np.linalg.inv(B)
                for delta in [0.,1e-10,2e-8]:
                    Y=X.copy();Y[0,0]+=delta/np.max(abs(B))
                    add(f'random_n{n}_r{rep}_{scaling}_d{delta}',B,Y)
    # Exact power-of-two reciprocal pairs with severe cancellation must fall back.
    for exp in [5,10,15,20,30,50,100,200]:
        v=2.**exp;add(f'triangular_{exp}',[[1,v],[0,1]],[[1,-v],[0,1]],None)

    head=header(out,pkg)+WRITER+'mata:\n'
    for i,(name,B,X,expect) in enumerate(tests):
        for label,A in [('B',B),('X',X)]:
            head+=f'{label}=J({len(A)},{A.shape[1]},.)\n'
            head+=''.join(f'{label}[{j+1},.]={mata(row[None,:])}\n' for j,row in enumerate(A))
        head+=f'savebin("g{i}_B.bin",B)\nsavebin("g{i}_X.bin",X)\n'
        head+=f'savebin("g{i}_fast.bin",_rsm_inverse_fast(B,X))\n'
        head+=f'savebin("g{i}_left.bin",B*X-I(rows(B)))\nsavebin("g{i}_right.bin",X*B-I(rows(B)))\n'
        head+=f'savebin("g{i}_radius.bin",64*epsilon(1)*rows(B)*rows(B)*(max(abs(B))*max(abs(X))+1))\n'
        if name.startswith(('wide_products_','cutoff_')):
            head+=f'savebin("g{i}_strictleft.bin",_rsm_inverse_residual(B,X,I(rows(B))))\nsavebin("g{i}_strictright.bin",_rsm_inverse_residual(X,B,I(rows(B))))\n'
    head+=r'''assert(_rsm_inverse_fast(J(0,0,0),J(0,0,0))==0)
    assert(_rsm_inverse_fast(J(2,3,1),I(2))==0)
    assert(_rsm_inverse_fast(I(2),I(3))==0)
    assert(_rsm_inverse_fast((.,0\0,1),I(2))==0)
    assert(_rsm_inverse_fast(I(2),(.,0\0,1))==0)
    assert(_rsm_inverse_fast(J(2,2,0),I(2))==0)
    assert(_rsm_inverse_fast(I(2),J(2,2,0))==0)
    assert(_rsm_inverse_fast(J(4097,0,0),J(0,0,0))==0)
    end
    di "GUARD_PASS"
    exit,clear
    '''
    (out/'test.do').write_text(head)
    p=subprocess.run([args.stata,'-b','do','test.do'],cwd=out,capture_output=True,timeout=240)
    log=(out/'test.log').read_text();assert p.returncode==0 and 'GUARD_PASS' in log.splitlines(),log[-4000:]
    def exact(A,B):
        af=[[F(float(v)) for v in row] for row in A];bf=[[F(float(v)) for v in row] for row in B]
        return [[sum(x*y for x,y in zip(row,col))-int(i==j) for j,col in enumerate(zip(*bf))] for i,row in enumerate(af)]
    rows=[]
    for i,(name,B,X,expect) in enumerate(tests):
        assert np.array_equal(B,readbin(out/f'g{i}_B.bin')) and np.array_equal(X,readbin(out/f'g{i}_X.bin'))
        accepted=bool(readbin(out/f'g{i}_fast.bin')[0,0]);radius=float(readbin(out/f'g{i}_radius.bin')[0,0]);emax=[];errors=[]
        if expect is not None:assert accepted==bool(expect),(name,accepted,expect)
        for side,A,C in [('left',B,X),('right',X,B)]:
            E=exact(A,C);emax.append(max(abs(v) for row in E for v in row))
            observed=readbin(out/f'g{i}_{side}.bin')
            errors.append(max(abs(F(float(v))-z) for row,erow in zip(observed,E) for v,z in zip(row,erow)) if np.isfinite(observed).all() else None)
        eligible=(len(B)<=4096 and np.max(abs(B))>0 and np.max(abs(X))>0 and
            all(np.all((A==0)|((abs(A)>=2.**-200)&(abs(A)<=2.**200))) for A in [B,X]))
        bound_verified=None
        if eligible:
            assert np.isfinite(radius) and all(e is not None and e<=F(radius) for e in errors),name
            a=F(float(np.max(abs(B))));b=F(float(np.max(abs(X))))
            assert F(radius)>=64*F(2.**-53)*len(B)**2*(a*b+1),name
            bound_verified=True
        if accepted:
            assert eligible and radius<=1e-8/4 and max(emax)<F(1e-8),name
        if name.startswith('wide_products_'):
            assert not accepted and any(not np.isfinite(readbin(out/f'g{i}_strict{s}.bin')).all() for s in ['left','right'])
        if name.startswith('cutoff_'):
            strict=max(np.max(abs(readbin(out/f'g{i}_strict{s}.bin'))) for s in ['left','right'])
            assert (strict<=1e-8)==(max(emax)<=F(1e-8)),name
        rows.append(dict(name=name,fast_accepted=accepted,exact_residual=[float(v) for v in emax],evaluation_error=[float(v) if v is not None else None for v in errors],radius=radius if np.isfinite(radius) else None,bound_verified=bound_verified))
    summary=dict(pairs=len(rows),fast_accepted=sum(x['fast_accepted'] for x in rows),fallback=sum(not x['fast_accepted'] for x in rows),invalid_controls=8,exact_false_acceptances=0,bound_verified=sum(x['bound_verified'] is True for x in rows),T2_false_acceptances_blocked=8,extreme_unavailable_preserved=3,cutoff_neighbors=6,random_pairs=225)
    (out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)

    out=root/'matrices';out.mkdir()
    head=header(out,pkg)+WRITER+'mata:\n'
    for i,row in enumerate(fixture['matrices']):
        B=decode(row['B']);head+=assign('B',B)+f'X=_rsm_inv(B)\nsavebin("m{i}_B.bin",B)\nsavebin("m{i}_X.bin",X)\n'
    head+='end\ndi "INVERSE_MATRIX_PASS"\nexit,clear\n'
    run(out,head,'INVERSE_MATRIX_PASS',args.stata)
    rows=[]
    for i,row in enumerate(fixture['matrices']):
        B=readbin(out/f'm{i}_B.bin');X=readbin(out/f'm{i}_X.bin')
        assert np.array_equal(B,decode(row['B']),equal_nan=True)
        finite=bool(np.isfinite(X).all());assert finite==row['expected_finite'],row['name']
        record=dict(name=row['name'],finite=finite)
        if finite:
            left=exact(B,X);right=exact(X,B);maximum=max(abs(v) for E in [left,right] for rr in E for v in rr)
            assert maximum<=F(1e-8),row['name'];record['exact_max_residual']=float(maximum)
        rows.append(record)
    (out/'results.json').write_text(json.dumps(rows,indent=2)+'\n')
    (out/'summary.json').write_text(json.dumps(dict(matrices=len(rows),finite=sum(x['finite'] for x in rows),false_acceptances=0),indent=2)+'\n')

    out=root/'scaling';out.mkdir();cases=[]
    for exp in [-1000,-900,-500,500,900,1000]:
        for multiplier in [1.,1.1,1.3,1.5,1.7,1.9]:
            b=multiplier*2.**exp
            for sign in [-1,1]:
                center=(1+sign*1e-8)/b
                for shift in range(-3,4):
                    x=center
                    for _ in range(abs(shift)):x=np.nextafter(x,0. if shift<0 else np.inf)
                    cases.append(dict(exp=exp,multiplier=multiplier,sign=sign,shift=shift,b=b,x=float(x)))
    head=header(out,pkg)+WRITER+'mata:\n'
    head+='savebin("powers.bin",(6.223015277861142e-61,1.6069380442589903e60,1.157920892373162e77))\n'
    for i,c in enumerate(cases):
        head+=f'B={c["b"]:.17g}\nX={c["x"]:.17g}\nsavebin("s{i}.bin",(B,X,_rsm_inverse_residual(B,X,1)))\n'
    head+='end\ndi "INVERSE_SCALING_PASS"\nexit,clear\n'
    run(out,head,'INVERSE_SCALING_PASS',args.stata)
    rows=[]
    for i,c in enumerate(cases):
        b,x,v=readbin(out/f's{i}.bin').ravel();assert b==c['b'] and x==c['x'] and np.isfinite(v)
        value=1-F(c['b'])*F(c['x'])
        assert (abs(v)<=1e-8)==(abs(value)<=F(1e-8)) and v==float(value),(i,v,str(value))
        rows.append(dict(case=i,**c,exact_residual=str(value),returned_residual=float(v),passed=True))
    assert np.array_equal(readbin(out/'powers.bin').ravel(),[2.**-200,2.**200,2.**256])
    (out/'results.json').write_text(json.dumps(rows,indent=2)+'\n')
    (out/'summary.json').write_text(json.dumps(dict(cases=len(rows),wrong_classifications=0,correctly_rounded=len(rows)),indent=2)+'\n')
    summaries={n:json.loads((root/n/'summary.json').read_text()) for n in ['guard','matrices','scaling']}
    (root/'summary.json').write_text(json.dumps(dict(all_passed=True,**summaries),indent=2)+'\n');print(json.dumps(summaries,indent=2))

if __name__=='__main__':main()
