#!/usr/bin/env python3
"""Audit undesirable DDF against independent HiGHS LPs and primal residuals.

Requires NumPy, SciPy, git and Stata 16+. Production runs are compared with
private instrumented engine copies that record lambda/beta without altering
the LP. Optional pre-F25 sources reproduce the defect and check the x path.
"""
from _support import package_files
import argparse
import csv
import difflib
import hashlib
import itertools
import json
from pathlib import Path
import re
import shutil
import subprocess
import time

import numpy as np
import scipy
from scipy.optimize import linprog

from test_output_protection import unchanged
from test_score_classification import check_display, capture
from test_saveref import instrument as force_retry

REPO = Path(__file__).resolve().parents[1]
CANONICAL = dict(crs='crs', ccr='crs', vrs='vrs', bcc='vrs', drs='drs', nirs='drs', irs='irs', ndrs='irs')


def read(path):
    with path.open() as f:
        rows = list(csv.DictReader(f))
    return {key: np.array([float(r[key]) if r[key] and not re.fullmatch(r'\.[a-z]?', r[key]) else np.nan for r in rows])
            for key in rows[0] if key != 'label'}


def oracle(x, y, b, k, direction, rts):
    """Build beta-first constraints from projected observations (not ado rows)."""
    n = len(x); x0, y0, b0 = x[k], y[k], b[k]
    gx = x0 if direction in ['x', 'xyb'] else np.zeros_like(x0)
    gy = y0 if direction in ['yb', 'xyb'] else np.zeros_like(y0)
    gb = b0 if direction in ['yb', 'xyb'] else np.zeros_like(b0)
    # Normalize rows for HiGHS using original units. Positive row scaling does
    # not change the feasible technology, lambda or the optimum beta.
    sx = np.maximum(np.max(abs(x), axis=0), 1e-300)
    sy = np.maximum(np.max(abs(y), axis=0), 1e-300)
    sb = np.maximum(np.max(abs(b), axis=0), 1e-300)
    aub = np.vstack([np.column_stack([gx/sx, (x/sx).T]),
                     np.column_stack([gy/sy, -(y/sy).T])])
    bub = np.r_[x0/sx, -y0/sy]
    aeq = np.column_stack([gb/sb, (b/sb).T]); beq = b0/sb
    convexity = np.r_[0., np.ones(n)]
    if rts == 'vrs': aeq = np.vstack([aeq, convexity]); beq = np.r_[beq, 1.]
    if rts in ['drs', 'irs']:
        sign = 1 if rts == 'drs' else -1
        aub = np.vstack([aub, sign*convexity]); bub = np.r_[bub, sign]
    fit = linprog(np.r_[-1., np.zeros(n)], A_ub=aub, b_ub=bub,
                  A_eq=aeq, b_eq=beq, bounds=(0, None), method='highs',
                  options={'primal_feasibility_tolerance': 1e-9, 'dual_feasibility_tolerance': 1e-9})
    if not fit.success: raise AssertionError(fit.message)
    return fit.x[0]


def residual(x, y, b, sol, direction, rts):
    lam, beta = sol[:, :-1], sol[:, -1:]
    gx = x if direction in ['x', 'xyb'] else np.zeros_like(x)
    gy = y if direction in ['yb', 'xyb'] else np.zeros_like(y)
    gb = b if direction in ['yb', 'xyb'] else np.zeros_like(b)
    sx = np.maximum(np.max(abs(x), axis=0), 1e-300)
    sy = np.maximum(np.max(abs(y), axis=0), 1e-300)
    sb = np.maximum(np.max(abs(b), axis=0), 1e-300)
    vals = [0., float(-np.min(sol)), float(np.max((lam@x+beta*gx-x)/sx)),
            float(np.max((y+beta*gy-lam@y)/sy)), float(np.max(abs(lam@b+beta*gb-b)/sb))]
    if rts == 'vrs': vals.append(float(np.max(abs(lam.sum(axis=1)-1))))
    if rts == 'drs': vals.append(float(np.max(lam.sum(axis=1)-1)))
    if rts == 'irs': vals.append(float(np.max(1-lam.sum(axis=1))))
    return max(vals)


def fixture(dims=(2, 2, 2), seed=551, factors=None, zeros=False):
    rng = np.random.default_rng(seed); n = 18
    x = rng.uniform(1, 5, (n, dims[0])); y = rng.uniform(1, 5, (n, dims[1])); b = rng.uniform(1, 5, (n, dims[2]))
    # Two duplicate frontier candidates and a dominated observation.
    x[1] = x[0]; y[1] = y[0]; b[1] = b[0]
    x[2] = x[0]*2; y[2] = y[0]*.5; b[2] = b[0]*2
    data = np.column_stack([x, y, b])
    if zeros:
        data[::3, 0] = 0; data[1::4, dims[0]] = 0; data[2::4, dims[0]+dims[1]] = 0
        # Include an identically zero bad-output column as a redundant equality.
        if dims[2] > 1: data[:, -1] = 0
    if factors is not None: data *= np.array(factors)
    names = [f'{prefix}{j+1}' for prefix, count in zip(['x', 'y', 'b'], dims) for j in range(count)]
    s = 'input double(id '+' '.join(names)+')\n'
    s += '\n'.join(str(i+1)+' '+' '.join(f'{v:.17g}' for v in row) for i, row in enumerate(data))+'\nend\n'
    return s+'gen str12 label="unit"+string(id)\ngen unrelated=991\n', names


def instrument(source, dest):
    dest.mkdir(); patches = []
    for engine in ['ipm', 'simplex']:
        p = source/f'_dea_{engine}.ado'; original = p.read_text()
        start = original.index(f'void _dea_{engine}_ddf_undesirable(')
        end = original.index('\n}\n', start)+3
        function = original[start:end]
        function = function.replace('real matrix X, Y, B', 'real matrix X, Y, B, audit_solution', 1)
        anchor = '    for (k = 1; k <= n; k++) {'
        assert function.count(anchor) == 1
        function = function.replace(anchor, '    audit_solution = J(n, n+1, .)\n'+anchor)
        anchor = '            beta = sol[n+1]'
        assert function.count(anchor) == 1
        function = function.replace(anchor, '            audit_solution[k,.] = sol\n'+anchor)
        function = function[:-2]+f'    st_matrix("F25_{engine}", audit_solution)\n}}\n'
        modified = original[:start]+function+original[end:]
        (dest/p.name).write_text(modified)
        patches.extend(difflib.unified_diff(original.splitlines(True), modified.splitlines(True), fromfile=p.name, tofile='test-only/'+p.name))
    (dest/'instrumentation.patch').write_text(''.join(patches))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata', required=True); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--baseline-ado', type=Path); ap.add_argument('--only')
    args = ap.parse_args()
    if args.baseline_ado and not (args.baseline_ado/'dea.ado').is_file(): ap.error('Invalid pre-F25 ado snapshot')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    pkg = out/'current_ado'; pkg.mkdir(); hashes = {}
    names = package_files(REPO)
    for name in names:
        data = (REPO/name).read_bytes(); (pkg/Path(name).name).write_bytes(data); hashes[name] = hashlib.sha256(data).hexdigest()
    audit = out/'instrumented_engines'; instrument(pkg, audit)
    retry = out/'forced_retry_ado'; force_retry(pkg, retry)
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes, stata=args.stata,
        baseline=str(args.baseline_ado), numpy=np.__version__, scipy=scipy.__version__,
        score_tolerance=1e-6, normalized_residual_tolerance=1e-7), indent=2))
    results = []; errors = []; residuals = []; lp_count = 0

    def wanted(name): return not args.only or re.search(args.only, name)

    def run(name, script, source=pkg, verify=None):
        d = out/name; d.mkdir(); row = dict(test=name, passed=False); start = time.monotonic()
        head = 'version 16\nclear all\nset more off\nset linesize 120\n'
        for folder in ['PLUS', 'PERSONAL', 'SITE', 'OLDPLACE']:
            (d/folder).mkdir(); head += f'sysdir set {folder} "{d/folder}"\n'
        head += f'adopath ++ "{source.resolve()}"\n'; (d/'test.do').write_text(head+script+'\ndi "TEST_PASS"\n')
        try:
            subprocess.run([args.stata, '-b', 'do', 'test.do'], cwd=d, capture_output=True, timeout=180)
            log = (d/'test.log').read_text(errors='replace')
            if 'TEST_PASS' not in log.splitlines(): raise RuntimeError(log[-1800:])
            if verify: row.update(verify(d) or {})
            row['passed'] = True
        except Exception as e: row['error'] = repr(e)
        row['seconds'] = round(time.monotonic()-start, 3); results.append(row)
        (out/'results.json').write_text(json.dumps(results, indent=2)); print(json.dumps(row), flush=True)
        return d

    def numeric_case(name, dims, engine, rts, direction, variant=0, setup=None, names=None, historical=False):
        nonlocal lp_count
        if setup is None: setup, names = fixture(dims)
        xs = ' '.join(n for n in names if n.startswith('x')); ys = ' '.join(n for n in names if n.startswith('y')); bs = ' '.join(n for n in names if n.startswith('b'))
        q = ['', 'if id<=16', 'in 2/17', 'if id<=16 in 2/17'][variant%4]
        selected = ['1', 'id<=16', 'inrange(_n,2,17)', 'id<=16 & inrange(_n,2,17)'][variant%4]
        setup += f'gen byte selected={selected}\n'
        for v in names: setup += f'replace selected=0 if missing({v})\n'
        dmu = 'label' if variant%2 else 'id'
        prefix = f'dea {xs}={ys} {q},' if variant%2 else f'dea {ys} {q}, inputs({xs})'
        opts = f'model(undesirable) bad({bs}) direction({direction}) engine({engine}) rts({rts}) dmu({dmu})'
        cmd = prefix+' '+opts
        s = setup+'quietly ds\nlocal original_vars "`r(varlist)\'"\ntempfile original\nsave `original\'\n'
        s += cmd+f' gen(res) stage({1+variant%2}) orientation({["input","output"][variant%2]}) saving("saved.dta")\n'+capture()
        s += 'assert r(n_failed)==0\nmatrix scores=r(efficiency)\ncf `original_vars\' using `original\', all\nassert missing(res) if !selected\n'
        s += 'mata: assert(all(st_matrix("scores"):==select(st_data(.,"res"),st_data(.,"selected"))))\n'
        s += 'export delimited using "scores.csv", replace\n'
        s += 'preserve\nkeep if selected\nkeep '+dmu+' '+xs+' '+ys+' res\ntempfile expected\nsave `expected\'\nuse "saved.dta", clear\ncf _all using `expected\', all\nrestore\n'
        if not historical:
            # Engine loaders may reload through findfile during a command;
            # keep that lookup in the private copy as well as rebinding callers.
            s += f'program drop _all\nadopath ++ "{audit}"\n'
            for e in ['ipm', 'simplex']: s += f'quietly run "{audit}/_dea_{e}.ado"\n'
            s += cmd+f' gen(audited) stage({2-variant%2}) orientation({["output","input"][variant%2]})\n'
            s += 'local final_engine "`r(engine)\'"\nassert res==audited\nmatrix S=F25_`final_engine\'\n'
            s += 'preserve\nclear\nsvmat double S, names(s)\nexport delimited using "solution.csv", replace\nrestore\n'
        def verify(d):
            nonlocal lp_count
            data = read(d/'scores.csv'); chosen = data['selected']==1
            if historical: return {}
            check_display(d, 'undesirable', 1e-7, dmu)
            x, y, b = [np.column_stack([data[n][chosen] for n in group.split()]) for group in [xs, ys, bs]]
            predicted = np.array([oracle(x, y, b, k, direction, CANONICAL[rts.lower()]) for k in range(len(x))]); lp_count += len(x)
            actual = data['res'][chosen]; error = float(np.max(abs(actual-predicted)))
            sol = np.column_stack(list(read(d/'solution.csv').values()))
            resid = residual(x, y, b, sol, direction, CANONICAL[rts.lower()])
            assert np.isfinite(actual).all() and np.isfinite(sol).all()
            assert error < 1e-6, ('independent score error', error, actual.tolist(), predicted.tolist())
            assert resid < 1e-7, ('normalized primal residual', resid)
            assert np.max(abs(sol[:, -1]-actual)) < 1e-12
            extra = {}
            if args.baseline_ado and (out/('old_'+name)/'scores.csv').exists():
                old = read(out/('old_'+name)/'scores.csv')
                # Scaling and tighter solver tolerance deliberately improve
                # numerical precision; the x LP itself is unchanged.
                delta = float(np.max(abs(actual-old['res'][chosen])))
                olderror = float(np.max(abs(predicted-old['res'][chosen])))
                assert delta < 1e-6 or (olderror >= 1e-6 and error < olderror), ('Unexplained historical x change', delta)
                extra.update(historical_x_max_delta=delta, historical_x_oracle_error=olderror)
            if name.startswith('units_') and not name.startswith('units_base_'):
                base = out/('units_base_'+engine+'_'+rts+'_'+direction)/'scores.csv'
                if base.exists():
                    delta = float(np.max(abs(actual-read(base)['res'])))
                    assert delta < 1e-6, ('Unit invariance', delta)
                    extra['unit_change_max_delta'] = delta
            errors.append(error); residuals.append(resid)
            return dict(lp_count=len(x), max_score_error=error, max_primal_residual=resid, **extra)
        return run(('old_' if historical else '')+name, s, args.baseline_ado if historical else pkg, verify)

    grid = itertools.product([(1,1,1),(2,2,2),(3,2,2)], ['simplex','ipm','auto'], ['crs','vrs','nirs','ndrs'], ['yb','xyb','x'])
    for i, (dims, engine, rts, direction) in enumerate(grid):
        name = f'grid_{"".join(map(str,dims))}_{engine}_{rts}_{direction}'
        if not wanted(name): continue
        if args.baseline_ado and direction=='x': numeric_case(name,dims,engine,rts,direction,i,historical=True)
        numeric_case(name,dims,engine,rts,direction,i)

    # One-coordinate zeros, a zero equality, unit changes and independent DGPs.
    for seed, engine, rts, direction in itertools.product([7,921], ['simplex','ipm'], ['crs','vrs','nirs','ndrs'], ['yb','xyb','x']):
        name = f'zeros_{seed}_{engine}_{rts}_{direction}'
        if wanted(name):
            setup, names = fixture(seed=seed, zeros=True)
            numeric_case(name,(2,2,2),engine,rts,direction,setup=setup,names=names)
    for factors, tag in [([1]*6,'base'),([1e-6,1e6,1e3,1e-3,1e-4,1e4],'mixed6'),([1e-10,1e10,1e8,1e-8,1e-6,1e6],'mixed10'),([1e-8]*6,'small'),([1e8]*6,'large')]:
        for engine,rts,direction in itertools.product(['simplex','ipm'],['crs','vrs','nirs','ndrs'],['yb','xyb','x']):
            name=f'units_{tag}_{engine}_{rts}_{direction}'
            if wanted(name):
                setup,names=fixture(factors=factors)
                numeric_case(name,(2,2,2),engine,rts,direction,setup=setup,names=names)

    for engine,direction in itertools.product(['simplex','ipm','auto'],['yb','xyb','x']):
        name=f'identical_{engine}_{direction}'
        if not wanted(name): continue
        setup='input id x y b\n1 1 1 1\n2 1 1 1\nend\n'
        cmd=f'dea x=y, dmu(id) model(undesirable) bad(b) direction({direction}) engine({engine}) gen(res)\n'
        run(name,setup+cmd+'assert r(N_eff)==2\nassert r(n_failed)==0\nassert abs(res)<1e-8\n')
        if args.baseline_ado: run('old_'+name,setup+cmd+f'assert abs(res-{0 if direction=="x" else 1})<1e-8\n',args.baseline_ado)

    # Active-variable missingness must agree with a physically restricted sample.
    for v,missing,engine in itertools.product(['x1','y1','b1','b2'],['.','.a','.z'],['simplex','ipm']):
        name=f'missing_{v}_{missing.replace(".","dot")}_{engine}'
        if not wanted(name): continue
        setup,names=fixture();setup+=f'replace {v}={missing} in 3\n'
        numeric_case(name,(2,2,2),engine,'vrs','yb',3,setup,names)
    for role,direction,engine in itertools.product(['x1','y1','b1'],['yb','xyb','x'],['simplex','ipm']):
        name=f'negative_{role}_{direction}_{engine}'
        if not wanted(name):continue
        setup,_=fixture();setup+=f'replace {role}=-1 in 1\ngen str20 res="preserve"\n'
        cmd=f'dea x1 x2=y1 y2, model(undesirable) bad(b1 b2) direction({direction}) engine({engine}) gen(res) replacevars'
        run(name,setup+unchanged(cmd,459))
        # Invalid values outside if/in do not invalidate a valid selected sample.
        run('excluded_'+name,setup+cmd.replace('=y1 y2,','=y1 y2 if id>1,')+'\nassert missing(res) in 1\n')
    for direction,engine in itertools.product(['yb','xyb','x'],['simplex','ipm']):
        name=f'zero_direction_{direction}_{engine}'
        if not wanted(name):continue
        setup,names=fixture();active=[n for n in names if (direction=='xyb' or n.startswith('x')==(direction=='x'))]
        setup+=''.join(f'replace {v}=0 in 1\n' for v in active)+'gen str20 res="preserve"\n'
        cmd=f'dea x1 x2=y1 y2, model(undesirable) bad(b1 b2) direction({direction}) engine({engine}) gen(res) replacevars'
        run(name,setup+unchanged(cmd,459))
    for form,engine in itertools.product(['A','B'],['simplex','ipm']):
        name=f'no_observations_{form}_{engine}'
        if not wanted(name):continue
        setup,_=fixture();setup+='replace b1=.\ngen str20 res="preserve"\n'
        cmd=('dea x1 x2=y1 y2,' if form=='A' else 'dea y1 y2, inputs(x1 x2)')+f' model(undesirable) bad(b1 b2) engine({engine}) gen(res) replacevars'
        run(name,setup+unchanged(cmd,2000))
    for engine, rts, direction in itertools.product(['simplex','ipm'], ['CCR','BCC','DRS','IRS'], ['yb','xyb','x']):
        name = f'alias_{engine}_{rts}_{direction}'
        if not wanted(name): continue
        setup,_=fixture()
        opts=f'inputs(x1 x2) model(undesirable) bad(b1 b2) engine({engine}) dmu(id)'
        canonical=CANONICAL[rts.lower()]
        s=setup+f'dea y1 y2, {opts} rts({canonical}) direction({direction}) gen(control)\n'
        s+='local count=r(N_eff)\nlocal avg=r(mean_eff)\n'
        # Full names, uppercase aliases, and documented minimum option spellings.
        short=f'i(x1 x2) mod(undesirable) bad(b1 b2) engine({engine}) dmu(id)'
        s+=f'dea y1 y2, {short} r({rts}) dir({direction.upper()}) gen(res)\n'
        s+=f'assert "`r(rts)\'"=="{canonical}"\nassert r(N_eff)==`count\'\nassert r(mean_eff)==`avg\'\nassert res==control\n'
        if direction=='yb':
            s+=f'dea y1 y2, {opts} rts({canonical}) gen(default)\nassert default==control\n'
        run(name,s)
    for direction,rts in itertools.product(['yb','xyb','x'],['crs','vrs']):
        name=f'retry_{direction}_{rts}'
        if not wanted(name):continue
        setup,_=fixture()
        opts=f'inputs(x1 x2) model(undesirable) bad(b1 b2) direction({direction}) rts({rts}) dmu(id)'
        s=setup+f'dea y1 y2, {opts} engine(ipm) gen(control)\n'
        s+='local count=r(N_eff)\nlocal avg=r(mean_eff)\n'
        s+=f'dea y1 y2, {opts} engine(auto) gen(res) saving("saved.dta")\n'
        s+='assert "`r(engine)\'"=="ipm"\nassert r(n_failed)==0\nassert r(N_eff)==`count\'\nassert r(mean_eff)==`avg\'\nassert res==control\nuse "saved.dta",clear\nassert !missing(res)\n'
        run(name,s,retry)
    for value,tag in [('0','zero'),('-1','negative'),('.','missing'),('.a','extended')]:
        name=f'invalid_tolerance_{tag}'
        if not wanted(name):continue
        setup,_=fixture();setup+='gen str20 res="preserve"\n'
        cmd=f'dea x1 x2=y1 y2, model(undesirable) bad(b1 b2) tolerance({value}) gen(res) replacevars'
        run(name,setup+unchanged(cmd,198))
    for direction,engine,rts in itertools.product(['yb','xyb','x'],['simplex','ipm'],['crs','vrs','nirs','ndrs']):
        name=f'singleton_{direction}_{engine}_{rts}'
        if not wanted(name):continue
        s='input id x y b\n1 1 1 1\nend\n'
        s+=f'dea x=y, model(undesirable) bad(b) direction({direction}) engine({engine}) rts({rts}) dmu(id) gen(res)\n'
        run(name,s+'assert r(N)==1\nassert r(N_eff)==1\nassert r(n_failed)==0\nassert abs(res)<1e-8\n')
    for direction,engine in itertools.product(['yb','xyb'],['simplex','ipm','auto']):
        name=f'unbounded_{direction}_{engine}'
        if not wanted(name):continue
        s='input id x y b\n1 0 1 0\n2 0 1 0\nend\n'
        s+=f'dea x=y, model(undesirable) bad(b) direction({direction}) engine({engine}) dmu(id) gen(res)\n'
        run(name,s+'assert r(N_eff)==0\nassert r(n_failed)==2\nassert missing(res)\n')
    # These are controlled invalid-solution injections, not empirical failure rates.
    mutations=dict(nonfinite='x[n]=.;',width='x=J(1,n-1,0);',negative='x[1]=-1;',
                   inequality='x[n]=2;',equality='x=J(1,n,0);')
    for fault,mutation in mutations.items():
        for engine in ['simplex','ipm','auto']:
            name=f'reject_solution_{fault}_{engine}'
            if not wanted(name):continue
            package=out/(name+'_ado');shutil.copytree(pkg,package)
            # R08 moved acceptance into the shared native helper. Corrupt the
            # candidate before that gate, not after a successful validation.
            p=package/'_dea_nonte_lp.ado';text=p.read_text()
            anchor='        x=S.parameters()';assert text.count(anchor)==1
            p.write_text(text.replace(anchor,anchor+'\n        '+mutation))
            setup,_=fixture()
            s=setup+f'dea x1 x2=y1 y2, model(undesirable) bad(b1 b2) engine({engine}) dmu(id) gen(res)\n'
            run(name,s+'assert r(N_eff)==0\nassert r(n_failed)==_N\nassert missing(res)\n',package)
    summary=dict(total=len(results),passed=sum(r['passed'] for r in results),failed=sum(not r['passed'] for r in results),
                 independent_lp_count=lp_count,max_score_error=max(errors,default=None),max_primal_residual=max(residuals,default=None))
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
    raise SystemExit(summary['failed']>0)


if __name__=='__main__':main()
