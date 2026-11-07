#!/usr/bin/env python3
"""Check core score counts/ranks with real Stata and explicit score fixtures.

Standard library only. An optional pre-F07 ado snapshot proves historical score
preservation and the original defects. Solver optimality is a separate test.
"""
from _support import package_files
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
import subprocess
import time

from test_saveref import BASIC, instrument

REPO = Path(__file__).resolve().parents[1]
ZERO = {'ddf', 'undesirable', 'additive'}
SUPER = {'supereff', 'supersbm'}
COMMON = ['te', 'multiplier', 'cost', 'revenue', 'profit', 'ddf', 'sbm',
          'additive', 'supereff', 'supersbm', 'allocative', 'undesirable']
FIXTURE = '''input id x y bad
1 1 1 1
2 2 3 2
3 4 4 3
4 5 3 4
5 2 3 2
6 8 4 5
7 3 2 3
8 6 5 4
end
gen double w=1
gen double p=10
gen double gx=1
gen double gy=1
gen str12 label="unit"+string(id)
'''


def records(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def number(s):
    try:
        return float(s)
    except ValueError:
        return None


def rank_rows(log):
    # Only compact score tables, not the reference tables or echoed commands.
    rows = re.findall(r'^\s*(\S+)\s*\|\s*(FAILED|[-+\d.eE]+)\s+(\d+|\.)\s*$', log, re.M)
    return {key: None if rank == '.' else int(rank) for key, _, rank in rows}


def check_display(d, model, tolerance, dmu='id'):
    rows = records(d/'scores.csv')
    selected = [r for r in rows if r['selected'] == '1']
    finite = [r for r in selected if number(r['res']) is not None]
    frontier = 0 if model in ZERO else 1
    values = [number(r['res']) for r in finite]
    expected_eff = sum(v >= 1-tolerance for v in values) if model in SUPER else sum(abs(v-frontier) <= tolerance for v in values)
    expected_super = sum(v > 1+tolerance for v in values) if model in SUPER else None
    meta = records(d/'metadata.csv')[0]
    assert int(meta['N']) == len(selected)
    assert int(meta['N_eff']) == expected_eff, (meta, expected_eff)
    assert int(meta['n_failed']) == len(selected)-len(finite)
    assert number(meta['N_super']) == expected_super, (meta, expected_super)
    if values:
        for name, expected in [('mean_eff', sum(values)/len(values)), ('min_eff', min(values)), ('max_eff', max(values))]:
            assert math.isclose(float(meta[name]), expected, abs_tol=1e-12, rel_tol=1e-12), (name, meta[name], expected)
    # Independently compare each printed rank with the number of strictly better
    # tolerance-grid keys. Half-grid points are avoided by the fixed fixtures.
    keys = {}
    for row in finite:
        value = number(row['res'])
        key = frontier if abs(value-frontier) <= tolerance else math.copysign(math.floor(abs(value)/tolerance+.5)*tolerance, value)
        keys[row[dmu]] = key if model in ZERO else -key
    actual = rank_rows((d/'test.log').read_text())
    expected = {r[dmu]: None if number(r['res']) is None else 1+sum(k < keys[r[dmu]] for k in keys.values()) for r in selected}
    assert actual == expected, (actual, expected)
    log = (d/'test.log').read_text()
    assert ('DMU Inefficiency Scores' in log) == (model in ZERO)


def capture():
    # Save r() before any command can clear it, including the new optional scalar.
    names = ['N', 'N_eff', 'N_super', 'n_failed', 'mean_eff', 'min_eff', 'max_eff']
    return ''.join(f'local m_{n}=r({n})\n' for n in names)+'''
file open metadata using "metadata.csv", write replace
file write metadata "N,N_eff,N_super,n_failed,mean_eff,min_eff,max_eff" _n
'''+ 'file write metadata '+ ' "," '.join(f'%24.17g (`m_{n}\')' for n in names)+' _n\nfile close metadata\n'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata', required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--baseline-ado', type=Path)
    ap.add_argument('--only', help='Case-name regex, also includes needed baseline controls')
    args = ap.parse_args()
    if args.baseline_ado and not (args.baseline_ado/'dea.ado').is_file():
        ap.error('Invalid pre-F07 ado snapshot')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    pkg = out/'current_ado'; pkg.mkdir(); hashes = {}
    names = package_files(REPO)
    for name in names:
        data = (REPO/name).read_bytes(); (pkg/Path(name).name).write_bytes(data)
        hashes[name] = hashlib.sha256(data).hexdigest()
    forced = out/'forced_retry_ado'; instrument(pkg, forced)
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes, baseline=str(args.baseline_ado), stata=args.stata), indent=2))
    results = []

    def wanted(name):
        return not args.only or re.search(args.only, name)

    def run(name, script, source=pkg, verify=None):
        d = out/name; d.mkdir()
        head = 'version 16\nclear all\nset more off\nset linesize 120\n'
        for folder in ['PLUS', 'PERSONAL', 'SITE', 'OLDPLACE']:
            (d/folder).mkdir(); head += f'sysdir set {folder} "{d/folder}"\n'
        head += f'adopath ++ "{source.resolve()}"\n'
        (d/'test.do').write_text(head+script+'\ndi "TEST_PASS"\n')
        start = time.monotonic(); row = dict(test=name, passed=False)
        try:
            subprocess.run([args.stata, '-b', 'do', 'test.do'], cwd=d, capture_output=True, timeout=180)
            log = (d/'test.log').read_text(errors='replace')
            if 'TEST_PASS' not in log.splitlines():
                raise RuntimeError(log[-2400:])
            if verify: verify(d)
            row['passed'] = True
        except Exception as e:
            row['error'] = repr(e)
        row['seconds'] = round(time.monotonic()-start, 3); results.append(row)
        (out/'results.json').write_text(json.dumps(results, indent=2))
        print(json.dumps(row), flush=True)
        return d

    # Explicit binary-exact boundaries, ties, near-boundary roundoff, extended
    # missings, exclusions, and all-missing solver results for every common model.
    # Both endpoints also fit Stata syntax's 12-significant-digit real macro;
    # 2^-20 would be rounded by syntax itself before the score comparison.
    for model, orientation, exponent in itertools.product(COMMON, ['input', 'output'], [3, 16]):
        name = f'boundary_{model}_{orientation}_{exponent}'
        if not wanted(name): continue
        tol = 2.0**-exponent; frontier = 0 if model in ZERO else 1
        offsets = [0, -1, 1, -.5, .5, 2, 2, 4]
        ranks = [1]*5+[6, 6, 8] if model in ZERO else [4]*5+[2, 2, 1]
        if model not in ZERO:
            offsets += [-2]; ranks += [9]
        n_finite = len(offsets); n = n_finite+4
        lines = [f'{i} {frontier+v*tol:.17g} {r}' for i, (v, r) in enumerate(zip(offsets, ranks), 1)]
        lines += [f'{n_finite+i} {miss} .' for i, miss in enumerate(['.', '.a', '.z'], 1)]
        lines += [f'{n} 99 .']
        s = f'run "{pkg}/dea.ado"\ninput id double(res expected)\n'+'\n'.join(lines)+'\nend\n'
        s += f'gen byte selected=id<{n}\nclonevar original=res\n'
        s += f'_dea_score_summary, model({model}) score(res) touse(selected) tolerance({tol:.17g}) rankvar(rank)\n'
        s += f'assert r(n_efficient)=={8 if model in SUPER else 5}\nassert r(n_failed)==3\nassert r(n_super)=={3 if model in SUPER else "."}\nassert rank==expected\nassert res==original\n'
        display = f'_dea_display, model({model}) rts(crs) orientation({orientation}) generate(res) touse(selected) tolerance({tol:.17g}) stage(2) dmu(id)\n'
        s += display+'assert r(n_efficient)=='+str(8 if model in SUPER else 5)+'\nassert res==original\n'
        s += 'export delimited id res selected expected using "scores.csv", replace\n'
        def synthetic(d):
            rows = records(d/'scores.csv'); actual = rank_rows((d/'test.log').read_text())
            expected = {r['id']: number(r['expected']) for r in rows if r['selected']=='1'}
            assert actual == expected, (actual, expected)
        run(name, s, verify=synthetic)
        # The empty-finite-set result must remain missing, with no false efficient DMUs.
        if exponent == 3:
            s = f'run "{pkg}/dea.ado"\ninput id res\n1 .\n2 .a\n3 .z\nend\ngen selected=1\n'+display
            s += 'assert r(n_efficient)==0\nassert r(n_failed)==3\nassert missing(r(mean_efficiency),r(min_eff),r(max_eff))\n'
            if model in SUPER: s += 'assert r(n_super)==0\n'
            run(f'all_missing_{model}_{orientation}', s)

    def real(name, model, options, variant=0, source=pkg, old=False, fixture=FIXTURE):
        quals = ['', 'if id<=6', 'in 2/7', 'if id<=6 in 2/7']
        conditions = ['1', 'id<=6', 'inrange(_n,2,7)', 'id<=6 & inrange(_n,2,7)']
        v = variant % 4; dmu = 'label' if variant % 2 else 'id'
        setup = fixture+f'gen byte selected={conditions[v]}\n'
        setup += 'tempfile original\nquietly ds\nlocal original_vars "`r(varlist)\'"\nsave `original\'\n'
        extra = dict(cost='w(w)', revenue='p(p)', profit='w(w) p(p)', allocative='w(w) p(p)', undesirable='bad(bad)').get(model, '')
        if model == 'ddf' and variant % 2: extra += ' gx(gx) gy(gy)'
        prefix = f'dea x=y {quals[v]},' if variant % 2 == 0 else f'dea y {quals[v]}, inputs(x)'
        generate = 'ddfbeta(res)' if model == 'ddf' and variant % 2 else 'generate(res)'
        cmd = f'{prefix} dmu({dmu}) model({model}) {extra} {generate} {options}'
        if not old and variant % 5 == 0:
            setup += 'gen double res=991\n'; cmd += ' replacevars'
        saving = variant % 3 == 0
        if saving: cmd += ' saving("saved.dta") replace'
        s = setup+cmd+'\n'+capture()
        if saving and model=='te' and not old:
            s+='matrix SavedStage=r(stage_status)\nmatrix SavedSource=r(solution_source)\n'
        s += 'assert r(N)>0\nmatrix scores=r(efficiency)\n'
        s += 'assert missing(res) if !selected\ncf `original_vars\' using `original\', all\n'
        s += 'mata: assert(all(st_matrix("scores"):==select(st_data(.,"res"),st_data(.,"selected"))))\n'
        s += 'export delimited id label res selected using "scores.csv", replace\n'
        if saving:
            s += 'preserve\nkeep if selected\nkeep '+dmu+' x y res\ntempfile expected\nsave `expected\'\nuse "saved.dta", clear\n'
            if model=='te' and not old:
                s+='mata: assert(all(st_data(.,tokens("_dea_radial_status _dea_stage2_status _dea_score_source _dea_ref_stage")):==(st_matrix("SavedStage"),st_matrix("SavedSource"))))\n'
                s+='drop _dea_radial_status _dea_stage2_status _dea_score_source _dea_ref_stage\n'
            s+='cf _all using `expected\', all\ncapture ds __*\nassert _rc==111\nrestore\n'
        def verify(d):
            if not old:
                check_display(d, model, 1e-7, dmu)
                if args.baseline_ado:
                    a = records(d/'scores.csv'); b = records(out/('old_'+name)/'scores.csv')
                    if model == 'undesirable':
                        # F25 deliberately changes yb/xyb LPs and improves
                        # precision. Their pre-F07 scores are not an oracle.
                        assert [{k:v for k,v in r.items() if k!='res'} for r in a] == [{k:v for k,v in r.items() if k!='res'} for r in b]
                        if 'direction(x)' in options:
                            for current, prior in zip(a,b):
                                av,bv=number(current['res']),number(prior['res'])
                                assert av==bv if av is None or bv is None else abs(av-bv)<1e-6
                    else:
                        assert a == b, 'Raw scores or selected sample changed from pre-F07'
        return run(('old_' if old else '')+name, s, source, verify)

    # Complete 5 x 3 x 4 x 2 grid, with rotating stages, syntax forms, samples,
    # IDs, saving/replacement and custom DDF directions. Missings are legitimate
    # for infeasible leave-one-out problems; the table and returns must agree.
    grid = itertools.product(['ddf', 'additive', 'undesirable', 'supereff', 'supersbm'], ['simplex', 'ipm', 'auto'], ['crs', 'vrs', 'nirs', 'ndrs'], ['input', 'output'])
    for i, (model, engine, rts, ori) in enumerate(grid):
        name = f'grid_{model}_{engine}_{rts}_{ori}'
        if not wanted(name): continue
        opts = f'engine({engine}) rts({rts}) orientation({ori}) stage({1+i%2})'
        if model == 'undesirable': opts += f' direction({["yb", "xyb", "x"][i%3]})'
        if args.baseline_ado: real(name, model, opts, i, args.baseline_ado, old=True)
        real(name, model, opts, i)
    # Other callers of the same display, including non-oriented/oriented SBM.
    for i, (model, engine, ori) in enumerate(itertools.product([m for m in COMMON if m not in ZERO|SUPER], ['simplex', 'ipm'], ['input', 'output'])):
        name = f'control_{model}_{engine}_{ori}'
        if not wanted(name): continue
        opts = f'engine({engine}) orientation({ori}) rts(vrs)'
        if args.baseline_ado: real(name, model, opts, i, args.baseline_ado, old=True)
        real(name, model, opts, i)

    # Original four-DMU counterexamples, known analytic CRS DDF and super scores.
    for model, ori in itertools.product(['ddf', 'supereff', 'supersbm'], ['input', 'output']):
        name = f'original_{model}_{ori}'
        if not wanted(name): continue
        common = BASIC+f'dea x=y, model({model}) orientation({ori}) dmu(id) gen(res)\n'
        s = common+'assert r(N_eff)==1\n'
        if model in SUPER: s += 'assert r(N_super)==1\n'
        if model == 'ddf': s += 'assert abs(res-cond(id==2,0,cond(id==4,3/7,1/5)))<1e-7\n'
        if model == 'supereff': s += 'assert abs(res-cond(id==2,1.5,cond(id==4,.4,2/3)))<1e-7\n'
        s += 'export delimited id res using "scores.csv", replace\n'
        def original_rank(d):
            actual = rank_rows((d/'test.log').read_text())
            assert actual == {'1': 2, '2': 1, '3': 2, '4': 4}, actual
        run(name, s, verify=original_rank)
        if args.baseline_ado:
            oldcount = 0 if model == 'ddf' else 3 if ori == 'output' else 1
            run('old_'+name, common+f'assert r(N_eff)=={oldcount}\n', args.baseline_ado)

    # Exact duplicate frontier observations remain efficient, but not strictly
    # super-efficient; they are not dropped from N_eff merely because score=1.
    for model, engine, ori in itertools.product(SUPER, ['simplex', 'ipm'], ['input', 'output']):
        name = f'duplicate_frontier_{model}_{engine}_{ori}'
        if not wanted(name): continue
        s = 'input id x y\n1 1 1\n2 1 1\n3 2 1\nend\n'
        cmd = f'dea x=y, model({model}) engine({engine}) orientation({ori}) dmu(id) gen(res)\n'
        run(name, s+cmd+'assert r(N_eff)==2\nassert r(N_super)==0\nassert abs(res-1)<1e-7 in 1/2\n')
        if args.baseline_ado:
            run('old_'+name, s+cmd+f'assert r(N_eff)=={0 if ori=="input" else 1}\n', args.baseline_ado)

    # Invalidating simplex in a private package forces a real IPM retry.
    # Check that final counts and the added scalar survive display/save commands.
    for model, ori, refs in itertools.product(['ddf', 'supereff', 'supersbm'], ['input', 'output'], [False, True]):
        if refs and model != 'supereff': continue
        name = f'retry_{model}_{ori}_ref{int(refs)}'
        if not wanted(name): continue
        opts = f'model({model}) orientation({ori}) dmu(id)'
        s = BASIC+f'dea x=y, {opts} engine(ipm) gen(control)\n'
        s += 'local n_eff=r(N_eff)\nlocal n_super=r(N_super)\n'
        s += f'dea x=y, {opts} engine(auto) gen(res) saving("saved.dta") '+('reference saveref ' if refs else '')+'replacevars\n'
        s += 'assert "`r(engine)\'"=="ipm"\nassert r(N_eff)==`n_eff\'\nassert r(N_super)==`n_super\'\nassert r(n_failed)==0\nassert res==control\n'
        s += 'use "saved.dta", clear\nassert !missing(res)\n'
        run(name, s, forced)
    failed = [r for r in results if not r['passed']]
    print(json.dumps(dict(total=len(results), passed=len(results)-len(failed), failed=len(failed))), flush=True)
    raise SystemExit(bool(failed))


if __name__ == '__main__':
    main()
