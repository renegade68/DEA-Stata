#!/usr/bin/env python3
"""Bounded native LP retries, status evidence, exact structural guards and MILP.
Uses isolated ado copies. Fault injection tests policy, not mathematical status.
"""
from _support import package_files
import argparse
import csv
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import random
import re
import shutil
import subprocess
import time

import numpy as np
from scipy.optimize import linprog
from test_lp_residuals import residual

REPO = Path(__file__).resolve().parents[1]
TIGHT = '''matrix C=(0,0,0,0,-1,-1)
matrix A=(1,2,4,5,1,0\\1,3,4,3,0,-1)
matrix B=(2+718293/2251799813685248\\3)
'''
BOUND = 'matrix C=(1,2)\nmatrix A=(1,1)\nmatrix B=(2)\n'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata', required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--only')
    args = ap.parse_args(); out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    pkg = out/'ado'; pkg.mkdir(); hashes = {}
    for name in package_files(REPO):
        raw=(REPO/name).read_bytes(); (pkg/Path(name).name).write_bytes(raw)
        hashes[name]=hashlib.sha256(raw).hexdigest()
    (out/'manifest.json').write_text(json.dumps(hashes,indent=2)+'\n')
    results=[]; oracle=[]

    def run(name, script, source=pkg, check=None):
        if args.only and not re.search(args.only,name): return
        d=out/name;d.mkdir();start=time.monotonic()
        head='version 16\nclear all\nset more off\nset linesize 120\n'
        for s in ['PLUS','PERSONAL','SITE','OLDPLACE']:
            (d/s).mkdir();head+=f'sysdir set {s} "{d/s}"\n'
        head+=f'adopath ++ "{source}"\n'
        (d/'test.do').write_text(head+script+'\ndi "TEST_PASS"\n')
        row=dict(test=name,passed=False)
        try:
            p=subprocess.run([args.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=240)
            log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else ''
            assert p.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2200:] or p.stderr.decode(errors='replace')
            if check: row['evidence']=check(d)
            row['passed']=True
        except Exception as e: row['error']=repr(e)
        row['seconds']=round(time.monotonic()-start,3);results.append(row)
        print(json.dumps(row),flush=True);(out/'results.json').write_text(json.dumps(results,indent=2)+'\n')

    def clone(name):
        target=out/('source_'+name);shutil.copytree(pkg,target);return target

    # Exact binary-rational primal/dual proof for the tight feasible LP.
    b=Fraction(float(2+718293/2251799813685248))
    x=[Fraction(0),b/2,Fraction(0),Fraction(0),Fraction(0),3*b/2-3]
    c=[0,0,0,0,-1,-1]; ae=[[1,2,4,5,1,0],[1,3,4,3,0,-1]]; y=[Fraction(-3,2),Fraction(1)]
    assert all(v>=0 for v in x)
    assert [sum(a*v for a,v in zip(row,x)) for row in ae]==[b,Fraction(3)]
    assert all(sum(y[i]*ae[i][j] for i in range(2))<=c[j] for j in range(6))
    exact=sum(a*v for a,v in zip(c,x)); assert exact==y[0]*b+y[1]*3
    (out/'exact_certificate.json').write_text(json.dumps(dict(x=list(map(str,x)),dual=list(map(str,y)),objective=str(exact)),indent=2)+'\n')
    run('tight_exact_recovery',TIGHT+'''_dea_lpcontract, c(C) aeq(A) beq(B) engine(ipm) tolerance(1e-10)
assert r(status)==0 & r(initial_status)==1 & r(native_retries)==1
assert r(fallback_used)==0 & r(status_conflict)==0
assert "`r(engine)'"=="ipm" & "`r(certificate)'"=="none"
matrix H=r(attempts)
assert rowsof(H)==2 & colsof(H)==8
assert H[1,1]==2 & H[1,2]==0 & H[1,3]==1 & H[1,4]==1
assert H[2,2]==1 & H[2,3]==0 & H[2,4]==0 & H[2,5]==1e-10
assert H[2,7]<=1e-10 & abs(H[2,8]-r(objective))<1e-20
assert abs(r(objective)-(-4.784792562162465e-10))<1e-10
matrix R=(r(x),r(max_scaled_eq_violation))
clear
svmat double R
export delimited using result.csv,replace
''',check=lambda d: check_tight(d,b))

    # Sign/zero proofs do not depend on tolerance, engine or row normalization.
    cases=[
        ('eq_zero','matrix C=(1)\nmatrix A=(0)\nmatrix B=(1)\n','aeq(A) beq(B)',1,'zero_row'),
        ('ub_zero','matrix C=(1)\nmatrix A=(0)\nmatrix B=(-1)\n','aub(A) bub(B)',1,'zero_row'),
        ('eq_positive','matrix C=(1)\nmatrix A=(2)\nmatrix B=(-1)\n','aeq(A) beq(B)',1,'nonnegative_row'),
        ('eq_negative','matrix C=(1)\nmatrix A=(-2)\nmatrix B=(1)\n','aeq(A) beq(B)',1,'nonnegative_row'),
        ('ub_positive','matrix C=(1)\nmatrix A=(2)\nmatrix B=(-1)\n','aub(A) bub(B)',1,'nonnegative_row'),
        ('ray','matrix C=(-1)\nmatrix A=(-1)\nmatrix B=(0)\n','aub(A) bub(B)',2,'zero_origin_ray'),
        ('ray_empty','matrix C=(-1e-15)\n','',2,'zero_origin_ray'),
        ('ray_eq','matrix C=(-1,0)\nmatrix A=(0,1)\nmatrix B=(0)\n','aeq(A) beq(B)',2,'zero_origin_ray'),
    ]
    for (name,setup,opts,status,cert),engine,scale,tol in itertools.product(cases,['ipm','simplex','auto'],['none','row'],['1e-7','1e-12']):
        run(f'proof_{name}_{engine}_{scale}_{tol}',setup+f'''_dea_lpcontract, c(C) {opts} engine({engine}) scaling({scale}) tolerance({tol})
assert r(status)=={status} & r(fallback_used)==0 & r(native_retries)==0
assert "`r(certificate)'"=="{cert}"
assert missing(r(objective))
matrix X=r(x)
mata: assert(all(st_matrix("X"):>=.))
''')
    guards=[
        ('tiny_eq','matrix A=(1e-15)\nmatrix B=(0)\n','aeq(A) beq(B)'),
        ('tiny_ub','matrix A=(1e-15)\nmatrix B=(0)\n','aub(A) bub(B)'),
        ('finite_upper','matrix U=(2)\n','ub(U)'),
        ('positive_lower','matrix L=(1)\n','lb(L)'),
        ('negative_rhs','matrix A=(-1)\nmatrix B=(-1)\n','aub(A) bub(B)'),
        ('nonzero_eq','matrix A=(0,1)\nmatrix B=(1)\nmatrix C=(-1,0)\n','aeq(A) beq(B)'),
        ('negative_lower','matrix A=(1)\nmatrix B=(-1)\nmatrix L=(-2)\n','aub(A) bub(B) lb(L)'),
        ('zero_objective','matrix C=(0)\n',''),
    ]
    for (name,setup,opts),engine in itertools.product(guards,['ipm','simplex','auto']):
        run(f'proof_guard_{name}_{engine}','matrix C=(-1)\n'+setup+f'''_dea_lpcontract, c(C) {opts} engine({engine}) tolerance(1e-12)
assert "`r(certificate)'"=="none"
''')

    # Compare successful initial results bit-for-bit, and verify their objective
    # independently with HiGHS. Includes signs, nonzero lower bounds and row units.
    rng=random.Random(73102)
    for i in range(18):
        n=3; a=np.array([[rng.randint(1,8) for _ in range(n)] for _ in range(4)],float)
        lower=np.array([rng.randint(-2,1) for _ in range(n)],float)
        upper=lower+np.array([rng.randint(2,5) for _ in range(n)],float)
        feasible=lower+.4*(upper-lower); rhs=a@feasible+np.array([rng.randint(1,4) for _ in range(4)])
        cost=np.array([rng.randint(1,7) for _ in range(n)],float)*(-1 if i%2 else 1)
        sense='max' if i%3 else 'min'; sign=-1 if sense=='max' else 1
        opt=linprog(sign*cost,A_ub=a,b_ub=rhs,bounds=list(zip(lower,upper)),method='highs-ds');assert opt.success
        answer=sign*opt.fun; oracle.append(dict(case=i,objective=answer))
        def mat(v): return '('+'\\'.join(','.join(f'{z:.17g}' for z in row) for row in np.atleast_2d(v))+')'
        setup=f'matrix C={mat(cost)}\nmatrix A={mat(a)}\nmatrix B={mat(rhs[:,None])}\nmatrix L={mat(lower)}\nmatrix U={mat(upper)}\n'
        for engine,scale in itertools.product(['ipm','simplex','auto'],['none','row']):
            s=setup+f'''_dea_lpcontract, c(C) aub(A) bub(B) lb(L) ub(U) sense({sense}) engine({engine}) scaling({scale}) tolerance(1e-9)
assert r(status)==0 & abs(r(objective)-({answer:.17g}))<1e-6
assert r(native_retries)==0 & r(fallback_used)==0
matrix X=r(x)
matrix H=r(attempts)
assert rowsof(H)==1 & H[1,2]==0 & H[1,4]==0 & H[1,7]<=1e-9
'''
            if scale=='none':
                backend='ipm' if engine=='ipm' else 'simplex'; tail=',""' if backend=='ipm' else ''
                s+=f'''mata:
first=_dea_lpcontract_{backend}(st_matrix("C"),"{sense}",J(0,3,0),J(0,1,0),st_matrix("A"),st_matrix("B"),st_matrix("L"),st_matrix("U"),1e-9{tail})
assert(first.status==0 & all(first.x:==st_matrix("X")))
end
'''
            run(f'initial_success_{i}_{engine}_{scale}',s)
    (out/'highs_oracle.json').write_text(json.dumps(oracle,indent=2)+'\n')

    # Test only: alter initial reports/candidates after the native solver call.
    def initial_fault(name, backend, code=None, corrupt=False, target=None):
        target=target or clone(name);p=target/'_dea_lpcontract_impl.do';s=p.read_text()
        lo=s.index('struct _DEA_LPResult scalar _dea_lpcontract_'+backend+'(')
        hi=s.index('\nstruct _DEA_LPResult scalar _dea_lpcontract_simplex(',lo) if backend=='ipm' else s.index('\nvoid _dea_lpcontract_verify(',lo)
        part=s[lo:hi];pos=part.rindex('    return(out)')
        if corrupt: inj='    if (out.status==0) out.x[1]=out.x[1]+10\n'
        else: inj=f'''    out.status={code};out.reported_status={code};out.status_text="injected";out.x=J(1,cols(c),.);out.objective=.
    out.attempts[1,3]={code};out.attempts[1,4]={code}
'''
        p.write_text(s[:lo]+part[:pos]+inj+part[pos:]+s[hi:]);return target

    for code in [1,3]:
        target=initial_fault('native_'+str(code),'ipm',code)
        # Explicit ipm must never call the package simplex solver.
        p=target/'_dea_lpcontract_impl.do';s=p.read_text().replace('    raw=_rsm_solve_lp(', '    assert(0)\n    raw=_rsm_solve_lp(');p.write_text(s)
        run('native_initial_'+str(code),BOUND+f'''_dea_lpcontract, c(C) aeq(A) beq(B) sense(max) engine(ipm) tolerance(1e-10)
assert r(status)==0 & r(initial_status)=={code} & r(native_retries)==1
assert r(fallback_used)==0 & r(status_conflict)==0 & abs(r(objective)-4)<1e-7
''',target)
    target=initial_fault('native_corrupt','ipm',corrupt=True)
    run('native_residual_recovery',BOUND+'''_dea_lpcontract, c(C) aeq(A) beq(B) sense(max) engine(ipm) tolerance(1e-10)
assert r(status)==0 & r(initial_status)==0 & r(native_retries)==1
matrix H=r(attempts)
assert H[1,4]==90 & H[1,7]>1e-10 & H[2,4]==0 & H[2,7]<=1e-10
''',target)
    for code in [1,2,3]:
        target=initial_fault('simplex_'+str(code),'simplex',code)
        status=3 if code==2 else 0
        run('auto_simplex_'+str(code),BOUND+f'''_dea_lpcontract, c(C) aeq(A) beq(B) sense(max) engine(auto) tolerance(1e-10)
assert r(status)=={status} & r(fallback_used)==1 & r(initial_status)=={code}
assert r(status_conflict)=={int(code==2)}
matrix H=r(attempts)
assert rowsof(H)==2 & H[1,1]==1 & H[2,1]==2 & H[2,4]==0
''',target)
        if code==2:
            run('milp_status_conflict',BOUND+'''milp_bb, c(C) aeq(A) beq(B) integer(1 2) lpengine(auto) lptol(1e-10)
assert r(status)==5 & r(status_conflicts)==1 & r(hasinc)==0
assert r(lp_fallbacks)==1 & r(simplex_solves)==1 & r(ipm_solves)==1
''',target)
    target=initial_fault('native_unbounded','ipm',2)
    run('native_status_conflict',BOUND+'''_dea_lpcontract, c(C) aeq(A) beq(B) sense(max) engine(ipm) tolerance(1e-10)
assert r(status)==3 & r(status_conflict)==1 & r(initial_status)==2
assert missing(r(objective))
''',target)

    # All retries fail: cap remains three and every attempt stays inspectable.
    target=initial_fault('all_failed','ipm',3)
    p=target/'_dea_native_retry.ado';s=p.read_text();s=s.replace('if (S.converged() & S.errorcode()==0)', 'if (0 & S.converged() & S.errorcode()==0)')
    s=s.replace('status=(S.errorcode()==0 ? 3 : S.errorcode())','status=3')
    p.write_text(s)
    run('bounded_exhaustion',BOUND+'''_dea_lpcontract, c(C) aeq(A) beq(B) sense(max) engine(ipm) tolerance(1e-8)
assert r(status)==3 & r(native_retries)==3
matrix H=r(attempts)
assert rowsof(H)==4 & H[2,2]==1 & H[3,2]==2 & H[4,2]==3
assert H[2,5]==1e-8 & H[3,5]==1e-8 & H[4,5]==1e-12
assert H[2,4]==3 & H[3,4]==3 & H[4,4]==3
''',target)
    run('tight_skips_duplicate',BOUND+'''_dea_lpcontract, c(C) aeq(A) beq(B) sense(max) engine(ipm) tolerance(1e-12)
assert r(status)==3 & r(native_retries)==2
matrix H=r(attempts)
assert rowsof(H)==3
''',target)

    # Contradictory infeasible/unbounded reports and incomplete checks must
    # remain failures. Agreement alone is explicitly not a structural proof.
    for first,second,expected,conflict in [(1,2,3,1),(2,1,3,1),(1,1,1,0),(2,2,2,0),(1,3,3,0),(2,3,3,0)]:
        name=f'auto_reports_{first}_{second}'
        target=initial_fault(name,'simplex',first)
        initial_fault(name,'ipm',second,target=target)
        p=target/'_dea_native_retry.ado';s=p.read_text()
        anchor='out.status=3;out.x=J(1,cols(c),.);out.history=J(3,6,.)'
        s=s.replace(anchor,anchor+f'\n    out.history=J(3,1,1)*({second},{second},tol,0,.,.);return(out)')
        p.write_text(s)
        run(name,BOUND+f'''_dea_lpcontract, c(C) aeq(A) beq(B) sense(max) engine(auto) tolerance(1e-8)
assert r(status)=={expected} & r(status_conflict)=={conflict}
assert r(fallback_used)==1 & r(initial_status)=={first} & r(native_retries)==3
assert "`r(certificate)'"=="none" & missing(r(objective))
matrix H=r(attempts)
assert rowsof(H)==5 & H[1,3]=={first} & H[2,3]=={second}
''',target)
    target=initial_fault('native_inconsistent','ipm',1)
    p=target/'_dea_native_retry.ado';s=p.read_text()
    anchor='out.status=3;out.x=J(1,cols(c),.);out.history=J(3,6,.)'
    p.write_text(s.replace(anchor,anchor+'\n    out.history=J(3,1,1)*(2,2,tol,0,.,.);return(out)'))
    run('native_inconsistent_reports',BOUND+'''_dea_lpcontract, c(C) aeq(A) beq(B) engine(ipm) tolerance(1e-8)
assert r(status)==3 & r(status_conflict)==1 & r(initial_status)==1
''',target)

    # Every retry candidate is independently rejected in original coordinates.
    target=initial_fault('all_bad','ipm',corrupt=True)
    p=target/'_dea_native_retry.ado';s=p.read_text().replace('x=S.parameters();d=', 'x=S.parameters();x[1]=x[1]+10;d=');p.write_text(s)
    run('all_candidates_bad',BOUND+'''_dea_lpcontract, c(C) aeq(A) beq(B) sense(max) engine(ipm) tolerance(1e-8)
assert r(status)==3 & r(native_retries)==3
matrix H=r(attempts)
mata: assert(all(st_matrix("H")[.,4]:==90) & all(st_matrix("H")[.,7]:>1e-8))
''',target)

    # Force native retry at every otherwise-successful MILP node. The optimum
    # is independently enumerated, and counters must include those extra solves.
    target=initial_fault('milp_retry','ipm',1)
    run('milp_retry_counters','''matrix C=(3,2)
matrix A=(2,1\\1,2)
matrix B=(4\\4)
milp_bb, c(C) aub(A) bub(B) integer(1 2) lpengine(ipm) lptol(1e-10)
assert r(status)==0 & r(objective)==6 & r(hasinc)==1
assert r(native_retries)>0 & r(status_conflicts)==0
assert r(ipm_solves)==r(lp_solves) & r(simplex_solves)==0 & r(lp_fallbacks)==0
''',target)
    run('milp_real_boundary_recovery',TIGHT+'''milp_bb, c(C) aeq(A) beq(B) integer(1) sense(min) lpengine(ipm) lptol(1e-10)
assert r(status)==0 & r(hasinc)==1 & r(native_retries)==1
assert r(status_conflicts)==0 & r(lp_solves)==1 & r(ipm_solves)==1
assert abs(r(objective)-(-4.784792562162465e-10))<1e-10
''')
    # Preserve the translator's explicit radial band after removing its
    # duplicate auto dispatch. The injected failure affects one fixed face only.
    target=clone('explicit_radial_band')
    p=target/'_dea_lpcontract_impl.do';text=p.read_text()
    lo=text.index('struct _DEA_LPResult scalar _dea_lpcontract_dispatch(')
    hi=text.index('real rowvector _dea_lpcontract_diagnostics(',lo)
    part=text[lo:hi];pos=part.rindex('    return(out)')
    part=part[:pos]+'''    if (rows(Aeq)==2 & cols(c)==6 & st_numscalar("fail_face_once")==1) {
        st_numscalar("fail_face_once",0)
        out.status=3;out.x=J(1,cols(c),.);out.objective=.
    }
'''+part[pos:];p.write_text(text[:lo]+part+text[hi:])
    for engine,band in itertools.product(['ipm','auto'],['0','1e-7']):
        s=f'''scalar fail_face_once=1
input double x y
1 1
2 3
4 4
5 3
end
_dea_tecontract, inputs(x) outputs(y) gen(eff) status(status) rts(crs) orientation(input) stage(2) engine({engine}) tolerance(1e-10) radialtolerance({band})
assert scalar(fail_face_once)==0
matrix R=r(radial)
matrix A=r(radial_adjustment)
'''
        if band=='0':
            s+='assert r(n_failed)==1 & r(stage2_failures)==1\nassert missing(eff)==(_n==1)\nassert status==3*(_n==1)\n'
        elif engine=='ipm':
            # This feasible shifted face remains a native numerical failure in
            # Stata 16.1. A witness is lambda2=1/3, sin=b-2/3, sout=0.
            s+='''assert r(n_failed)==1 & r(stage2_failures)==1
assert missing(eff)==(_n==1) & status==3*(_n==1)
assert missing(A[1,1]) & R[1,1]+1e-7*(1+abs(R[1,1]))>2/3
assert r(native_retries)==4
'''
        else:
            s+='''assert r(n_failed)==0 & r(stage2_failures)==0
assert status==0 & !missing(eff)
assert abs(A[1,1]-1e-7*(1+abs(R[1,1])))<1e-14
assert abs(eff[1]-R[1,1]-A[1,1])<1e-14
mata: assert(all(st_matrix("A")[2..4,.]:==0))
'''
        run(f'explicit_radial_band_{engine}_{band}',s,target)

    run('public_history_survives_saving','''input double a b c d e f str2 rel double rhs
0 0 0 0 -1 -1 "" 0
1 2 4 5 1 0 "=" 2
1 3 4 3 0 -1 "=" 3
end
replace rhs=2+718293/2251799813685248 in 2
lp a b c d e f, min engine(ipm) tolerance(1e-10) saving(solution)
assert r(status)==0 & r(native_retries)==1 & r(initial_status)==1
matrix H=r(attempts)
assert H[1,3]==1 & H[2,4]==0
assert _N==3 & rhs[2]==2+718293/2251799813685248
''')
    summary=dict(cases=len(results),failed=[r['test'] for r in results if not r['passed']])
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');assert not summary['failed'],summary


def check_tight(d,b):
    with (d/'result.csv').open() as f: rows=list(csv.reader(f))
    a=np.array(list(map(float,rows[1])));x=a[:6]
    diag=residual(x,np.array([[1,2,4,5,1,0],[1,3,4,3,0,-1]],float),np.array([float(b),3]),np.empty((0,6)),np.empty(0),np.zeros(6),np.full(6,np.inf))
    assert max(diag[3:])<=1e-10
    assert abs(a[6]-diag[3])<1e-15
    return dict(max_scaled_residual=float(max(diag[3:])),objective=float(-x[-2:].sum()))


if __name__=='__main__': main()
