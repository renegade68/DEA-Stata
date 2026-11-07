#!/usr/bin/env python3
"""Validate failure-only LP recovery against exact feasible/dual witnesses.

Requires NumPy/SciPy and Stata16+. Hexadecimal fixtures preserve stored binary64 data.
The single known preexisting objective-boundary observation is explicitly reported;
it is not counted as an accurate new recovery. Improvements to that case are allowed.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,hashlib,shutil,subprocess
import numpy as np
from test_core_retries import WRITER,readbin
from test_simplex_stability import mata
REPO=Path(__file__).resolve().parents[1]
def decode(A,n):return np.array([[float.fromhex(v) for v in row] for row in A]).reshape(-1,n)
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def fm(A):return [[F(float(v)) for v in row] for row in A]
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    root=args.output.resolve();root.mkdir(parents=True,exist_ok=False);pkg=root/'ado';shutil.copytree(REPO/'ado',pkg)
    path=REPO/'tests/fixtures/precision_recovery.json';cases=json.loads(path.read_text())['cases'];assert len(cases)==396
    manifest=dict(source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in pkg.iterdir() if p.is_file()},fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),numpy=np.__version__,stata=args.stata)
    (root/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');body='version 16\nclear all\nset more off\nset processors 1\nset linesize 180\n'
    for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:(root/key).mkdir();body+=f'sysdir set {key} "{root/key}"\n'
    body+=WRITER+f'adopath ++ "{pkg}"\nrun "{pkg}/_dea_simplex.ado"\nmata:\n'
    body+=r'''
void precision_probe(string scalar name,real matrix E,real colvector b,real matrix U,real colvector u,real rowvector c,real scalar tol) {
 struct RSM_Prob scalar P
 struct RSM_Result scalar R
 _rsm_std_form(P,U,u,E,b,c);R=_rsm_two_phase(P,tol,"")
 savebin(name+"_result.bin",(R.status,R.objval,R.iters));savebin(name+"_x.bin",R.x)
 savebin(name+"_eqmat.bin",E);savebin(name+"_eqrhs.bin",b)
 savebin(name+"_ineqmat.bin",U);savebin(name+"_ineqrhs.bin",u);savebin(name+"_cost.bin",c)
}
'''
    inputs=[]
    for i,case in enumerate(cases):
        n=len(case['xstar']);E=decode(case['E'],n);b=decode(case['b'],1);U=decode(case['U'],n);u=decode(case['u'],1);c=decode(case['c'],n)
        x=list(map(F,case['xstar']));y=list(map(F,case['ystar']));cost=list(map(F,c.ravel()));ef=fm(E);uf=fm(U);rhs=list(map(F,b.ravel()));irhs=list(map(F,u.ravel()));obj=F(case['objective'])
        assert min(x)>=0 and all(z>0 for row in E for z in row) and (b>0).all() and (c>0).all()
        assert all(dot(row,x)==v for row,v in zip(ef,rhs)) and all(dot(row,x)<=v for row,v in zip(uf,irhs))
        rc=[ci-dot(y,col) for ci,col in zip(cost,zip(*ef))];assert min(rc)>=0 and dot(cost,x)==dot(y,rhs)==obj==6
        active=[k for k,v in enumerate(x) if v>0];assert len(active)==3 and all(rc[k]>0 for k in range(n) if k not in active)
        B=[[row[k] for k in active] for row in ef];det=B[0][0]*(B[1][1]*B[2][2]-B[1][2]*B[2][1])-B[0][1]*(B[1][0]*B[2][2]-B[1][2]*B[2][0])+B[0][2]*(B[1][0]*B[2][1]-B[1][1]*B[2][0]);assert det!=0
        inputs.append((E,b,U,u,c))
        for label,A in zip(['E','b','U','u','c'],inputs[-1]):
            body+=f'{label}=J({A.shape[0]},{A.shape[1]},0)\n'+''.join(f'{label}[{k+1},.]={mata(row[None,:])}\n' for k,row in enumerate(A))
        body+=f'precision_probe("p{i:03}",E,b,U,u,c,{case["tol"]:.17g})\n'
    body+='end\ndi "PRECISION_LP_PASS"\nexit,clear\n';(root/'test.do').write_text(body)
    proc=subprocess.run([args.stata,'-b','do','test.do'],cwd=root,capture_output=True,timeout=900)
    (root/'process.json').write_text(json.dumps(dict(returncode=proc.returncode,stderr=proc.stderr.decode(errors='replace')),indent=2)+'\n');log=(root/'test.log').read_text(errors='replace');assert proc.returncode==0 and 'PRECISION_LP_PASS' in log.splitlines(),log[-4000:]
    rows=[]
    for i,case in enumerate(cases):
        E,b,U,u,c=inputs[i];n=len(case['xstar'])
        for suffix,A in zip(['eqmat','eqrhs','ineqmat','ineqrhs','cost'],inputs[i]):assert np.array_equal(A,readbin(root/f'p{i:03}_{suffix}.bin'))
        out=readbin(root/f'p{i:03}_result.bin').ravel();xf=readbin(root/f'p{i:03}_x.bin').ravel();status=int(out[0]);assert status in [0,3]
        if case['required_available']:assert status==0,case['name']
        row=dict(name=case['name'],passed=True,status=status,required_exact=case['required_exact'])
        if status==0:
            assert len(xf)==n+len(u) and np.isfinite(xf).all();x=list(map(F,xf[:n]));cost=list(map(F,c.ravel()));obj=dot(cost,x);target=list(map(F,case['xstar']));tol=F(case['tol']);err=abs(obj-6)/7
            assert min(x)>=0 and abs(F(float(out[1]))-obj)<=F(1e-12)*(1+abs(obj))
            for A,bb,eq in [(E,b,True),(U,u,False)]:
                for a,v in zip(fm(A),map(F,bb.ravel())):
                    residual=abs(dot(a,x)-v) if eq else max(F(0),dot(a,x)-v)
                    scale=max(abs(v),dot(list(map(abs,a)),[max(F(1),abs(z)) for z in x]));assert scale>0 and residual<=tol*scale
            known=err>F(1e-8)
            if known:assert case.get('known_objective_boundary')==str(obj),(case['name'],float(err))
            if case['required_exact']:assert x==target and obj==6
            row.update(exact_primal=x==target,exact_objective=str(obj),relative_objective_error=float(err),known_objective_boundary=known)
        rows.append(row)
    summary=dict(cases=len(rows),passed=len(rows),all_passed=True,available=sum(x['status']==0 for x in rows),exact_recovery_requirements=sum(x['required_exact'] for x in rows),exact_primals=sum(x.get('exact_primal',False) for x in rows),known_preexisting_objective_boundary=sum(x.get('known_objective_boundary',False) for x in rows),exact_primal_dual_witnesses=len(cases),stored_inputs_verified=True)
    (root/'results.json').write_text(json.dumps(rows,indent=2)+'\n');(root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
