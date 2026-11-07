#!/usr/bin/env python3
"""Audit every pseudo-frontier LP using SciPy HiGHS on identical random draws.

Requires NumPy, SciPy, Git and Stata. A test-only source copy records the actual
normalized data, initial scores and resampled scores; production source is never
instrumented. This audits LPs and summaries, not the sampling model's coverage.
"""
from _support import package_files
import argparse
import csv
import difflib
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
import scipy
from scipy.optimize import linprog

from test_core_bootstrap import DATA, REPO


def solve(x0, y0, xref, yref, orientation, rts):
    n = len(xref)
    if orientation == 'input':
        cost = np.r_[1., np.zeros(n)]
        a = np.vstack([np.column_stack([-x0, xref.T]), np.column_stack([np.zeros(len(y0)), -yref.T])])
        b = np.r_[np.zeros(len(x0)), -y0]
    else:
        cost = np.r_[-1., np.zeros(n)]
        a = np.vstack([np.column_stack([np.zeros(len(x0)), xref.T]), np.column_stack([y0, -yref.T])])
        b = np.r_[x0, np.zeros(len(y0))]
    convexity = np.r_[0., np.ones(n)]
    ae = be = None
    if rts == 'vrs': ae, be = convexity[None, :], np.array([1.])
    if rts in ['drs', 'irs']:
        sign = 1 if rts == 'drs' else -1
        a = np.vstack([a, sign*convexity]); b = np.r_[b, sign]
    fit = linprog(cost, A_ub=a, b_ub=b, A_eq=ae, b_eq=be, bounds=(0, None), method='highs',
                  options={'primal_feasibility_tolerance':1e-9, 'dual_feasibility_tolerance':1e-9})
    if not fit.success: raise RuntimeError(fit.message)
    return fit.x[0] if orientation == 'input' else 1/fit.x[0]


def read_matrix(path):
    with path.open() as f:
        rows = csv.reader(f); next(rows)
        return np.array([[float(v) if v != '.' else np.nan for v in row] for row in rows])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata', required=True)
    ap.add_argument('--output', required=True, type=Path)
    args = ap.parse_args(); out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    pkg = out/'instrumented_ado'; pkg.mkdir()
    hashes = {}
    names = package_files(REPO)
    for name in names:
        data = (REPO/name).read_bytes(); (pkg/Path(name).name).write_bytes(data)
        hashes[name] = hashlib.sha256(data).hexdigest()
    f = pkg/'_dea_bootstrap_core.ado'; original = f.read_text()
    modified = original.replace('real matrix    X, Y, eff_boot, X_pseudo, Y_pseudo',
                                'real matrix    X, Y, eff_boot, X_pseudo, Y_pseudo, audit_beta', 1)
    modified = modified.replace('eff_boot = J(n, B, .)', 'eff_boot = J(n, B, .)\n    audit_beta = J(n, B, .)', 1)
    modified = modified.replace('beta_b = _sw_sample(eff_orig, h)', 'beta_b = _sw_sample(eff_orig, h)\n        audit_beta[.,b] = beta_b', 1)
    needle = '    // ── Display'
    audit = '''    st_matrix("AUDIT_X", X)
    st_matrix("AUDIT_Y", Y)
    st_matrix("AUDIT_ORIG", eff_orig)
    st_matrix("AUDIT_DRAWS", audit_beta)
    st_matrix("AUDIT_BOOT", eff_boot)

'''
    assert needle in modified
    modified = modified.replace(needle, audit+needle, 1); f.write_text(modified)
    (out/'instrumentation.patch').write_text(''.join(difflib.unified_diff(original.splitlines(True), modified.splitlines(True),
                                                                    fromfile='production', tofile='test-only')))
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes, numpy=np.__version__, scipy=scipy.__version__,
        instrumented_sha256=hashlib.sha256(f.read_bytes()).hexdigest(), stata=args.stata,
        method='SciPy HiGHS independently solves fixed-data LPs; empirical inverse-CDF quantiles via NumPy'), indent=2))
    results = []
    for dimensions in [1, 2]:
        for engine in ['simplex', 'ipm']:
            for rts in ['crs', 'vrs', 'drs', 'irs']:
                for orientation in ['input', 'output']:
                    name = f'{dimensions}d_{engine}_{rts}_{orientation}'; work = out/name; work.mkdir()
                    header = 'version 16\nclear all\nset more off\nset linesize 120\n'
                    for folder in ['plus', 'personal', 'site', 'oldplace']:
                        (work/folder).mkdir(); header += f'sysdir set {folder.upper()} "{work/folder}"\n'
                    header += f'adopath ++ "{pkg}"\n'
                    extra = 'gen double x2=1+mod(id*7,5)\ngen double y2=1+mod(id*3,4)\n' if dimensions==2 else ''
                    io = 'x x2=y y2' if dimensions==2 else 'x=y'
                    body = DATA+extra+f'''set seed 714
    dea {io} if selected, dmu(id) rts({rts}) ort({orientation}) engine({engine}) tolerance(1e-8) bootstrap reps1(51) level(90) gen(eff)
    matrix AUDIT_RESULT=r(bootstrap)
    '''
                    for matrix in ['X', 'Y', 'ORIG', 'DRAWS', 'BOOT', 'RESULT']:
                        body += f'''preserve
    clear
    svmat double AUDIT_{matrix}, names(v)
    export delimited using "{matrix}.csv", replace
    restore
    '''
                    (work/'test.do').write_text(header+body+'\ndi "TEST_PASS"\n')
                    start = time.monotonic()
                    row = dict(test=name, passed=False)
                    try:
                        subprocess.run([args.stata, '-b', 'do', 'test.do'], cwd=work, capture_output=True, timeout=180)
                        log = (work/'test.log').read_text(errors='replace')
                        if 'TEST_PASS' not in log.splitlines(): raise RuntimeError(log[-1400:])
                        x, y, orig, draws, boot, result = [read_matrix(work/f'{m}.csv') for m in ['X','Y','ORIG','DRAWS','BOOT','RESULT']]
                        orig = orig[:,0]; n, bcount = draws.shape
                        if not all(np.isfinite(m).all() for m in [x,y,orig,draws,boot,result]): raise RuntimeError('Missing audit values')
                        initial = np.array([solve(x[i], y[i], x, y, orientation, rts) for i in range(n)])
                        expected = np.empty_like(boot)
                        for b in range(bcount):
                            xr = (orig/draws[:,b])[:,None]*x if orientation == 'input' else x
                            yr = (draws[:,b]/orig)[:,None]*y if orientation == 'output' else y
                            for i in range(n): expected[i,b] = solve(x[i], y[i], xr, yr, orientation, rts)
                        bias = expected.mean(axis=1)-orig
                        low, high = np.quantile(expected-orig[:,None], [.05,.95], axis=1, method='inverted_cdf')
                        summaries = np.column_stack([orig-bias, bias, expected.std(axis=1, ddof=1), orig-high, orig-low])
                        row.update(lp_count=n*(bcount+1), max_initial_gap=float(np.max(abs(initial-orig))),
                                   max_replication_gap=float(np.max(abs(expected-boot))), max_summary_gap=float(np.max(abs(summaries-result))))
                        row['passed'] = row['max_initial_gap']<2e-5 and row['max_replication_gap']<1e-7 and row['max_summary_gap']<1e-7
                        np.savetxt(work/'independent_replications.csv', expected, delimiter=',', fmt='%.17g')
                        np.savetxt(work/'independent_summaries.csv', summaries, delimiter=',', fmt='%.17g')
                    except Exception as exc:
                        row['error'] = str(exc)
                    row['seconds'] = round(time.monotonic()-start, 3)
                    results.append(row); (out/'results.json').write_text(json.dumps(results, indent=2)); print(json.dumps(row), flush=True)
    summary = dict(total=len(results), passed=sum(r['passed'] for r in results), lp_count=sum(r.get('lp_count',0) for r in results))
    (out/'summary.json').write_text(json.dumps(summary, indent=2)); print(json.dumps(summary), flush=True)
    raise SystemExit(not all(r['passed'] for r in results))


if __name__ == '__main__':
    main()
