#!/usr/bin/env python3
"""Residual acceptance: exact boundary witnesses, rescaling, fault injection,
dual efficiency and integer rounding. Requires NumPy/SciPy and licensed Stata.
All instrumentation is confined to copied ado files in --output.
"""
from _support import package_files
import argparse
import csv
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np
from scipy.optimize import linprog

REPO = Path(__file__).resolve().parents[1]
DATA = 'input double id x y\n1 1 1\n2 2 4\n3 4 2\nend\n'


def residual(x, ae, be, au, bu, lb, ub):
    """Independent NumPy computation, directly in original coordinates."""
    result = []
    scaled = []
    for a, b, equality in [(ae, be, True), (au, bu, False)]:
        v = abs(a @ x - b) if equality else np.maximum(a @ x - b, 0)
        den = np.maximum(abs(b), abs(a) @ np.maximum(1, abs(x)))
        result.append(float(v.max(initial=0)))
        scaled.append(float(np.divide(v, den, out=np.zeros_like(v), where=den != 0).max(initial=0)))
    v = np.r_[np.maximum(lb - x, 0), np.maximum(x - ub, 0)]
    den = np.r_[np.maximum.reduce([np.ones(len(x)), abs(lb), abs(x)]),
                np.maximum.reduce([np.ones(len(x)), abs(ub), abs(x)])]
    # Exclude only absent bounds, never arithmetic overflow at finite bounds.
    use = np.r_[np.isfinite(lb), np.isfinite(ub)]
    return np.r_[result, v[use].max(initial=0), scaled, (v[use] / den[use]).max(initial=0)]


def check_recovered_projection(d, ori):
    with (d / 'projection.csv').open() as f:
        rows=list(csv.reader(f))
    a=np.array([[float(v) for v in row] for row in rows[1:]])
    x=np.array([1,2,4,5.]);y=np.array([1,3,4,3.])
    ae=np.array([np.r_[x,1,0],np.r_[y,0,-1]])
    maximum=0.; optimal=0;nonoptimal=0
    np.testing.assert_array_equal(a[:,6],a[:,7]); assert np.all(a[:,8]==0)
    for k,row in enumerate(a):
        b=np.array([row[6]*x[k],y[k]]) if ori=='input' else np.array([x[k],row[6]*y[k]])
        diag=residual(row[:6],ae,b,np.empty((0,6)),np.empty(0),np.zeros(6),np.full(6,np.inf))
        maximum=max(maximum,max(diag[3:])); assert max(diag[3:])<=1e-4+1e-15
        opt=linprog([0,0,0,0,-1,-1],A_eq=ae,b_eq=b,bounds=(0,None),method='highs-ds')
        if opt.success:
            optimal+=1; assert abs(opt.fun+row[4:6].sum())<2e-3
        else: nonoptimal+=1
    return dict(max_scaled_residual=maximum,highs_optimal=optimal,highs_nonoptimal=nonoptimal)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata', required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    pkg = out / 'ado'; pkg.mkdir(); hashes = {}
    for name in package_files(REPO):
        raw = (REPO / name).read_bytes(); (pkg / Path(name).name).write_bytes(raw)
        hashes[name] = hashlib.sha256(raw).hexdigest()
    (out / 'manifest.json').write_text(json.dumps(hashes, indent=2) + '\n')
    results = []

    def run(name, script, source=pkg, check=None):
        d = out / name; d.mkdir(); start = time.monotonic()
        head = 'version 16\nclear all\nset more off\nset linesize 120\n'
        for s in ['PLUS', 'PERSONAL', 'SITE', 'OLDPLACE']:
            (d / s).mkdir(); head += f'sysdir set {s} "{d/s}"\n'
        head += f'adopath ++ "{source}"\n'
        (d / 'test.do').write_text(head + script + '\ndi "TEST_PASS"\n')
        row = dict(test=name, passed=False)
        try:
            proc = subprocess.run([args.stata, '-b', 'do', 'test.do'], cwd=d, capture_output=True, timeout=240)
            log = (d / 'test.log').read_text(errors='replace') if (d / 'test.log').exists() else ''
            assert proc.returncode == 0 and 'TEST_PASS' in log.splitlines(), log[-1800:] or proc.stderr.decode(errors='replace')
            if check: row['evidence'] = check(d)
            row['passed'] = True
        except Exception as e:
            row['error'] = repr(e)
        row['seconds'] = round(time.monotonic() - start, 3)
        results.append(row); print(json.dumps(row), flush=True)
        (out / 'results.json').write_text(json.dumps(results, indent=2) + '\n')

    # Exact separator for the infeasible fixture: (3/2)*row1-row2 has
    # nonnegative coefficients but negative RHS. A feasible control has
    # lambda2=1, slack_in=delta, slack_out=0 (all other variables zero).
    ae0 = np.array([[1, 2, 4, 5, 1, 0], [1, 3, 4, 3, 0, -1]], dtype=float)
    assert np.all(1.5 * ae0[0] - ae0[1] >= 0)
    certificates = []
    for delta, scales in itertools.product([-1e-8, 718293 / 2251799813685248, .125],
                                           [(1, 1), (1e-8, 1e8), (1e8, 1e-8)]):
        ae = [[Fraction(float(v)) for v in row] for row in ae0 * np.array(scales)[:, None]]
        be = [Fraction(float(v)) for v in np.array([2 + delta, 3]) * scales]
        if delta < 0:
            alpha = max(ae[1][j] / ae[0][j] for j in range(4))
            assert all(alpha * ae[0][j] - ae[1][j] >= 0 for j in range(6))
            margin = alpha * be[0] - be[1]; assert margin < 0
            certificate = dict(kind='infeasible_separator', row1_multiplier=str(alpha), negative_rhs=str(margin))
        else:
            witness = [Fraction(0), Fraction(1), Fraction(0), Fraction(0), (be[0] - ae[0][1]) / ae[0][4], Fraction(0)]
            assert all(v >= 0 for v in witness)
            assert all(sum(a * v for a, v in zip(row, witness)) == b for row, b in zip(ae, be))
            certificate = dict(kind='feasible_witness', x=[str(v) for v in witness])
        certificates.append(dict(delta=delta, row_multipliers=scales, **certificate))
    (out / 'exact_certificates.json').write_text(json.dumps(certificates, indent=2) + '\n')
    for delta, tol, engine, scaling, scales in itertools.product(
            [-1e-8, 718293 / 2251799813685248, .125], [1e-7, 1e-10, 1e-12],
            ['simplex', 'ipm', 'auto'], ['none', 'row'], [(1, 1), (1e-8, 1e8), (1e8, 1e-8)]):
        index = len(results)
        s = f'''matrix C=(0,0,0,0,-1,-1)
matrix A=(1,2,4,5,1,0\\1,3,4,3,0,-1)
matrix B=(2+({delta:.17g})\\3)
matrix W=({scales[0]:.17g},0\\0,{scales[1]:.17g})
matrix A=W*A
matrix B=W*B
_dea_lpcontract, c(C) aeq(A) beq(B) engine({engine}) scaling({scaling}) tolerance({tol})
matrix R=(r(status),r(objective),r(max_eq_violation),r(max_ineq_violation),r(max_bound_violation),r(max_scaled_eq_violation),r(max_scaled_ineq_violation),r(max_scaled_bound_violation),r(reported_status),r(fallback_used),r(x),r(candidate_x))
assert r(feasibility_tolerance)=={tol}
if r(status)==0 assert max(r(max_scaled_eq_violation),r(max_scaled_ineq_violation),r(max_scaled_bound_violation))<={tol}
clear
svmat double R
export delimited using result.csv, replace
'''
        def check(d, delta=delta, tol=tol, engine=engine, scales=scales):
            with (d / 'result.csv').open() as f:
                values = next(csv.reader(f)); values = next(csv.reader(f))
            a = np.array([float(v) if v not in ('.', '') else np.nan for v in values])
            status, obj = a[:2]; x = a[10:16]; candidate = a[16:22]
            ae = ae0 * np.array(scales)[:, None]; be = np.array([2 + delta, 3]) * scales
            if delta < 0 and tol <= 1e-10: assert status != 0, (status, a)
            if delta == .125 and scales == (1, 1): assert status == 0
            if np.isfinite(candidate).all():
                diag = residual(candidate, ae, be, np.empty((0, 6)), np.empty(0), np.zeros(6), np.full(6, np.inf))
                np.testing.assert_allclose(a[2:8], diag, rtol=2e-6, atol=2e-14 * max(1, max(scales)))
                np.testing.assert_allclose(a[5:8], diag[3:], rtol=2e-6, atol=3e-16)
                if status == 0:
                    np.testing.assert_array_equal(x, candidate)
                    assert max(diag[3:]) <= tol + 3e-16
                    assert abs(obj + x[-2:].sum()) < 1e-12
            if status != 0: assert np.isnan(x).all() and np.isnan(obj)
            return dict(status=int(status), candidate_available=bool(np.isfinite(candidate).all()), tolerance=tol)
        run(f'boundary_{index:03}', s, check=check)

    # Public table API preserves rejected candidates and residual diagnostics.
    run('public_lp_rejection', '''input double a b c d e f str2 rel double rhs
0 0 0 0 -1 -1 "" 0
1 2 4 5 1 0 "=" 1.99999999
1 3 4 3 0 -1 "=" 3
end
lp a b c d e f, min engine(simplex) tolerance(1e-12)
assert r(status)==3
assert "`r(status_text)'"=="residual check failed"
assert r(max_scaled_eq_violation)>1e-12
matrix X=r(x)
matrix Q=r(candidate_x)
mata: assert(all(st_matrix("X"):>=.) & !hasmissing(st_matrix("Q")))
''')

    # Force an invalid candidate after a backend's success return. The
    # dispatcher must reject it, and auto must try its already-defined fallback.
    for engine in ['simplex', 'ipm']:
        target = out / ('fault_' + engine); shutil.copytree(pkg, target)
        p = target / '_dea_lpcontract_impl.do'; text = p.read_text()
        lo = text.index('struct _DEA_LPResult scalar _dea_lpcontract_' + engine + '(')
        hi = text.index('\nvoid _dea_lpcontract_verify(', lo) if engine == 'simplex' else text.index('\nstruct _DEA_LPResult scalar _dea_lpcontract_simplex(', lo)
        part = text[lo:hi]; pos = part.rindex('    return(out)')
        part = part[:pos] + '    if (out.status==0) out.x[1]=out.x[1]+10\n' + part[pos:]
        p.write_text(text[:lo] + part + text[hi:])
        if engine == 'ipm':
            # This fixture corrupts every candidate, including native retries.
            # Initial-only corruption/recovery is tested in test_native_retries.
            p = target / '_dea_native_retry.ado'
            p.write_text(p.read_text().replace('x=S.parameters();d=', 'x=S.parameters();x[1]=x[1]+10;d='))
        s = 'matrix C=(1)\nmatrix A=(1)\nmatrix B=(2)\n'
        s += f'_dea_lpcontract, c(C) aub(A) bub(B) sense(max) engine({engine}) tolerance(1e-10)\n'
        s += 'assert r(status)==3 & r(reported_status)==0\nassert r(max_scaled_ineq_violation)>1e-10\n'
        s += 'assert "`r(status_text)\'"=="residual check failed"\n'
        if engine == 'simplex':
            s += '_dea_lpcontract, c(C) aub(A) bub(B) sense(max) engine(auto) tolerance(1e-10)\nassert r(status)==0 & r(fallback_used)==1\nassert "`r(engine)\'"=="ipm"\nassert abs(r(objective)-2)<1e-8\n'
        run('reject_fault_' + engine, s, target)

    # The same exact output-oriented dual objective is phi, not efficiency.
    # Analytic frontier: DMU 3 has phi=4 under CRS/IRS and phi=2 under VRS/DRS.
    for rts, ori, tol in itertools.product(['crs', 'vrs', 'drs', 'irs'], ['input', 'output'], [1e-7, 1e-10, 1e-12]):
        x = np.array([[1.], [2.], [4.]]); y = np.array([[1.], [4.], [2.]])
        if ori == 'output':
            au = np.r_[np.c_[x.T, [0]], np.c_[-y.T, [2]]]; bu = np.array([4, 0]); c = np.array([0, 0, 0, -1])
        else:
            au = np.r_[np.c_[x.T, [-4]], np.c_[-y.T, [0]]]; bu = np.array([0, -2]); c = np.array([0, 0, 0, 1])
        ae = be = None
        if rts == 'vrs': ae = np.array([[1, 1, 1, 0]]); be = [1]
        if rts in ['drs', 'irs']:
            sign = 1 if rts == 'drs' else -1
            au = np.vstack([au, sign * np.array([1, 1, 1, 0])]); bu = np.r_[bu, sign]
        opt = linprog(c, A_ub=au, b_ub=bu, A_eq=ae, b_eq=be, bounds=(0, None), method='highs-ds')
        assert opt.success; answer = -1 / opt.fun if ori == 'output' else opt.fun
        run(f'dual_{rts}_{ori}_{tol}', f'run "{pkg/"_dea_ipm.ado"}"\nmata: st_numscalar("answer",_dea_ipm_dual_{ori}((1\\2\\4),(1\\4\\2),3,"{rts}",{tol}))\nassert abs(answer-({answer:.17g}))<1e-6\n')

    fault = out / 'fault_stage2'; shutil.copytree(pkg, fault)
    p = fault / '_dea_ipm.ado'; text = p.read_text()
    assert text.count('sol = S2.parameters()') == 1
    p.write_text(text.replace('sol = S2.parameters()', 'sol = S2.parameters()\n        sol[1]=sol[1]+10'))
    p = fault / '_dea_native_retry.ado'
    p.write_text(p.read_text().replace('x=S.parameters();d=', 'x=S.parameters();x[1]=x[1]+10;d='))
    run('reject_stage2_candidate', DATA + '''dea x=y, dmu(id) engine(ipm) rts(vrs) tolerance(1e-10) reference
matrix S=r(stage_status)
matrix L=r(lambda)
assert r(n_stage2_failed)==3 & r(n_radial_failed)==0
mata: assert(all(st_matrix("S")[.,2]:==90) & all(st_matrix("L"):>=.))
''', fault)

    # Every native TE path reads the actual requested setting, including dual
    # recovery and the auxiliary VRS frontier. These assertions are test-only.
    audit = out / 'audit_tol'; shutil.copytree(pkg, audit)
    p = audit / '_dea_ipm.ado'; text = p.read_text()
    end = text.index('// 5. UTILITY:')
    part = text[:end]
    for obj in ['S', 'S2', 'lp']:
        anchor = obj + '.setTol(tol)'
        part = part.replace(anchor, anchor + f';assert({obj}.getTol()==tol)')
    p.write_text(part + text[end:])
    for ori, tol in itertools.product(['input', 'output'], [1e-7, 1e-10, 1e-12]):
        run(f'tolerance_delivery_{ori}_{tol}', DATA + f'dea x=y, dmu(id) engine(ipm) rts(vrs) ort({ori}) tolerance({tol}) reference\nrun "{audit/"_dea_ipm.ado"}"\nmata: assert(_dea_ipm_dual_{ori}((1\\2\\4),(1\\4\\2),3,"crs",{tol})<.)\n', audit)

    fault = out / 'fault_radial'; shutil.copytree(pkg, fault)
    p = fault / '_dea_ipm.ado'; text = p.read_text()
    lo = text.index('void _dea_ipm_envelopment(')
    hi = text.index('real scalar _dea_ipm_slacks(', lo)
    part = text[lo:hi]; assert part.count('sol = S.parameters()') == 2
    part = part.replace('sol = S.parameters()', 'sol = S.parameters()\n                sol[1]=sol[1]+100')
    p.write_text(text[:lo] + part + text[hi:])
    p = fault / '_dea_native_retry.ado'
    p.write_text(p.read_text().replace('x=S.parameters();d=', 'x=S.parameters();x[1]=x[1]+100;d='))
    for rts, ori in itertools.product(['crs', 'vrs', 'drs', 'irs'], ['input', 'output']):
        answer = (.25 if rts in ['crs', 'irs'] else .5) if ori == 'output' else (.25 if rts in ['crs', 'drs'] else 1/3)
        run(f'reject_radial_recover_dual_{rts}_{ori}', DATA + f'''dea x=y, dmu(id) engine(ipm) rts({rts}) ort({ori}) tolerance(1e-10) reference gen(eff)
matrix S=r(stage_status)
matrix L=r(lambda)
assert r(n_radial_failed)==0 & r(n_stage2_failed)==3
assert abs(eff[3]-({answer}))<1e-6
mata: assert(all(st_matrix("S")[.,2]:>=.) & all(st_matrix("L"):>=.))
''', fault)

    run('zero_rows_bounds_nonfinite', f'''run "{pkg/"_dea_lpcheck.ado"}"
mata:
empty=J(0,2,0); rhs=J(0,1,0)
assert(_dea_lpcheck_ok(_dea_lpcheck((-2,3),J(1,2,0),0,empty,rhs,(-2,.),(.,3)),1e-12))
assert(!_dea_lpcheck_ok(_dea_lpcheck((-2,3),J(1,2,0),1,empty,rhs,(-2,.),(.,3)),1e-12))
d=_dea_lpcheck((-2,3.25),empty,rhs,empty,rhs,(-2,.),(.,3))
assert(abs(d[3]-.25)<1e-15 & abs(d[6]-1/13)<1e-15)
assert(!_dea_lpcheck_ok(_dea_lpcheck((.,3),empty,rhs,empty,rhs,(.,.),(.,.)),1e-12))
assert(!_dea_lpcheck_ok(J(1,6,0),.a))
// Mata max() ignores missing: reject arithmetic overflow before taking maxima.
assert(!_dea_lpcheck_ok(_dea_lpcheck((-8e307,0),empty,rhs,empty,rhs,(8e307,.),(.,.)),1e-12))
assert(!_dea_lpcheck_ok(_dea_lpcheck((8e307,0),empty,rhs,empty,rhs,(.,.),(-8e307,.)),1e-12))
assert(!_dea_lpcheck_ok(_dea_lpcheck((2,0),(8e307,0),0,empty,rhs,(.,.),(.,.)),1e-12))
assert(!_dea_lpcheck_ok(_dea_lpcheck((8e307,-8e307),(1,1),0,empty,rhs,(.,.),(.,.)),1e-12))
end
''')

    audit = out / 'audit_contract_rejections'; shutil.copytree(pkg, audit)
    p = audit / '_dea_lpcontract_impl.do'; text = p.read_text()
    anchor = 'out.status=3; out.status_text="residual check failed"'
    assert text.count(anchor) == 1
    text = text.replace(anchor, '''assert(max(out.residual[4..6])>tol)
            st_numscalar("audit_rejects",st_numscalar("audit_rejects")+1)
            ''' + anchor)
    p.write_text(text)
    for ori in ['input', 'output']:
        run('contract_loose_rejections_' + ori, f'''scalar audit_rejects=0
input double x y
1 1
2 3
4 4
5 3
end
_dea_tecontract, inputs(x) outputs(y) gen(eff) status(status) rts(crs) orientation({ori}) stage(2) engine(ipm) tolerance(1e-4)
assert audit_rejects==3 & r(stage2_failures)==0 & r(n_failed)==0
assert status==0 & !missing(eff)
assert r(native_retries)==3 & r(status_conflicts)==0
matrix R=(r(lambda),r(slack_in),r(slack_out),r(radial),r(radial_stage2),r(radial_adjustment))
clear
svmat double R
export delimited using projection.csv,replace
''', audit, check=lambda d, ori=ori: check_recovered_projection(d,ori))

    # Near-integer LP optimum must not turn into an infeasible rounded incumbent.
    # Enumeration is exact: the only feasible nonnegative integer is zero.
    for engine, scaling in itertools.product(['simplex', 'ipm', 'auto'], ['none', 'row']):
        run(f'milp_rounding_{engine}_{scaling}', f'''matrix C=(1)
matrix U=(.99999999)
milp_bb, c(C) ub(U) integer(1) lpengine({engine}) lpscaling({scaling}) inttol(1e-7) lptol(1e-12)
assert r(status)==0 & r(hasinc)==1 & r(objective)==0
matrix X=r(x)
assert X[1,1]==0
''')

    summary = dict(cases=len(results), failed=[r['test'] for r in results if not r['passed']])
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    assert not summary['failed'], summary


if __name__ == '__main__':
    main()
