#!/usr/bin/env python3
"""F11 alpha/name mapping, independent LPs, compatibility and atomic rejection.

Requires NumPy/SciPy and Stata 16+. Each case uses a fresh isolated process.
Optional --baseline-ado is the pre-F11 package, not the original Git commit.
"""
from _support import package_files
import argparse
from collections import Counter
from decimal import Decimal, ROUND_HALF_UP
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time

import numpy as np
import scipy
from test_fuzzy_rts import DATA, SPECS, SAMPLES, read
from test_output_protection import unchanged
from verify_bootstrap_lp import solve

REPO = Path(__file__).resolve().parents[1]
GRIDS = [
    ('original', '.004 .001', ['.001', '.004']),
    ('middle', '.501 .499 .5', ['.499', '.5', '.501']),
    ('mixed', '1 .999 .004 .001 .25 .251 .5', ['.001', '.004', '.25', '.251', '.5', '.999', '1']),
    ('tiny', '1 0 .999999999999 .999999999998 2e-11 1e-11', ['0', '1e-11', '2e-11', '.999999999998', '.999999999999', '1']),
    ('step', '0(.001).012', [str(Decimal(i)/1000) for i in range(13)]),
    ('dense', '0(.0001).004', [str(Decimal(i)/10000) for i in range(41)]),
    ('boundary', '.005 .004999 .994999 .995 .995001', ['.004999', '.005', '.994999', '.995', '.995001']),
    ('precise', '.0010000002 .001 .0010000001', ['.001', '.0010000001', '.0010000002']),
]
LEGACY = [('', ['0', '.25', '.5', '.75', '1']),
          ('.751 .001 .999 .249 .501', ['.001', '.249', '.501', '.751', '.999']),
          ('0(.01).1', [str(Decimal(i)/100) for i in range(11)])]


def names(alphas, prefix='f'):
    buckets = [str((Decimal(a)*100).quantize(Decimal(1), rounding=ROUND_HALF_UP)) for a in alphas]
    totals = Counter(buckets)
    indices = Counter()
    suffixes = []
    for b in buckets:
        indices[b] += 1
        suffixes.append(b if totals[b] == 1 else f'{b}_{indices[b]}')
    return {s: [f'{prefix}_{s}_a{b}' for b in suffixes] for s in ['pess', 'opt']}


def command(dim=1, ori='input', mode='full', sample='', rts='crs', alpha='.004 .001', prefix='f', extra='', dmu='id'):
    x = ' '.join(f'x{i}' for i in range(1, dim+1))
    y = ' '.join(f'y{i}' for i in range(1, dim+1))
    bounds = []
    for opt, side, variables, modes in [
        ('xlo', 'lo', x, ['full', 'inputs', 'lower']),
        ('xhi', 'hi', x, ['full', 'inputs', 'upper']),
        ('ylo', 'lo', y, ['full', 'outputs', 'lower']),
        ('yhi', 'hi', y, ['full', 'outputs', 'upper']),
    ]:
        if mode in modes:
            bounds.append(f'{opt}('+ ' '.join(f'{side}_{v}' for v in variables.split())+')')
    return (f'dea_fuzzy {sample}, inputs({x}) outputs({y}) '+ ' '.join(bounds)
            + f' ort({ori}) rts({rts}) gen({prefix}) '
            + (f'dmu({dmu}) ' if dmu else '') + (f'alpha({alpha}) ' if alpha else '') + extra)


def mapping_checks(alphas, prefix='f', selected='selected'):
    expected = names(alphas, prefix)
    s = f'assert r(n_alpha)=={len(alphas)}\n'
    s += 'local levels "`r(alpha)\'"\n'
    for sense in ['pess', 'opt']:
        s += f'''assert "`r({sense}vars)'"=="{' '.join(expected[sense])}"
'''
    for i, a in enumerate(alphas, 1):
        s += f'''local a : word {i} of `levels'
assert `a'=={a}
'''
        for sense in ['pess', 'opt']:
            v = expected[sense][i-1]
            s += f'''confirm numeric variable {v}
assert "`:type {v}'"=="double"
local has_alpha=regexm("`:variable label {v}'", "α=([^,]+),")
local label_alpha=real(regexs(1))
assert `has_alpha' & `label_alpha'==`a'
assert !missing({v}) if {selected}
assert missing({v}) if !{selected}
'''
    return s


def check_lp(d, alphas, dim, ori, mode, rts):
    data = read(d/'results.csv'); chosen = data['selected'] == 1
    def columns(kind, side=''):
        return np.column_stack([data[f'{side}{kind}{i}'][chosen] for i in range(1, dim+1)])
    x, y = columns('x'), columns('y')
    xl = columns('x', 'lo_') if mode in ['full', 'inputs', 'lower'] else x
    xu = columns('x', 'hi_') if mode in ['full', 'inputs', 'upper'] else x
    yl = columns('y', 'lo_') if mode in ['full', 'outputs', 'lower'] else y
    yu = columns('y', 'hi_') if mode in ['full', 'outputs', 'upper'] else y
    expected = names(alphas); gap = 0.; count = 0
    for i, a in enumerate(alphas):
        a = float(a)
        lowx, highx = xl+a*(x-xl), xu-a*(xu-x)
        lowy, highy = yl+a*(y-yl), yu-a*(yu-y)
        for sense in ['pess', 'opt']:
            values = []
            for k in range(len(x)):
                if sense == 'pess':
                    x0, y0, xx, yy = highx[k], lowy[k], lowx.copy(), highy.copy()
                else:
                    x0, y0, xx, yy = lowx[k], highy[k], highx.copy(), lowy.copy()
                xx[k], yy[k] = x0, y0
                values.append(solve(x0, y0, xx, yy, ori, SPECS[rts])); count += 1
            actual = data[expected[sense][i]][chosen]
            assert np.isfinite(actual).all()
            delta = float(np.max(np.abs(actual-values))); gap = max(gap, delta)
            assert delta < 1e-6, (a, sense, delta)
    return dict(lp_count=count, max_oracle_gap=gap)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stata', required=True); p.add_argument('--output', required=True, type=Path)
    p.add_argument('--baseline-ado', type=Path); p.add_argument('--only')
    args = p.parse_args(); out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    pkg = out/'current_ado'; pkg.mkdir(); hashes = {}
    paths = package_files(REPO)
    for name in paths:
        data = (REPO/name).read_bytes(); (pkg/Path(name).name).write_bytes(data)
        hashes[name] = hashlib.sha256(data).hexdigest()
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes, numpy=np.__version__, scipy=scipy.__version__, baseline=str(args.baseline_ado)), indent=2))
    rows = []
    def run(name, script, source=pkg, check=None):
        if args.only and not re.search(args.only, name): return
        d = out/name; d.mkdir(); head = 'version 16\nclear all\nset more off\nset linesize 120\n'
        for folder in ['PLUS', 'PERSONAL', 'SITE', 'OLDPLACE']:
            (d/folder).mkdir(); head += f'sysdir set {folder} "{d/folder}"\n'
        head += f'adopath ++ "{source.resolve()}"\n'
        (d/'test.do').write_text(head+script+'\ndi "TEST_PASS"\n')
        row = dict(test=name, passed=False); start = time.monotonic()
        try:
            proc = subprocess.run([args.stata, '-b', 'do', 'test.do'], cwd=d, capture_output=True, timeout=240)
            log = (d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else proc.stderr.decode(errors='replace')
            assert 'TEST_PASS' in log.splitlines(), log[-2200:]
            if check: row.update(check(d))
            row['passed'] = True
        except Exception as e: row['error'] = str(e)
        row['seconds'] = round(time.monotonic()-start, 3); rows.append(row)
        (out/'results.json').write_text(json.dumps(rows, indent=2)); print(json.dumps(row), flush=True)

    # Every RTS spelling, orientation and bound mode; stagger dimensions,
    # selected samples and grids so close levels reach every solver branch.
    case = 0
    for rts in SPECS:
        for ori in ['input', 'output']:
            for mode in ['crisp', 'full', 'inputs', 'outputs', 'lower', 'upper']:
                grid, spec, alphas = GRIDS[case % len(GRIDS)]
                dim = 1+case % 2; sample = list(SAMPLES.values())[(case//2) % 3]
                case += 1
                opts = dict(dim=dim, ori=ori, mode=mode, sample=sample, rts=rts)
                selected = f'gen byte selected=0\nreplace selected=1 {sample}\n'
                s = DATA+selected+command(**opts, alpha=spec, extra='detail debug')+'\n'+mapping_checks(alphas)
                # Check association survives reversed alpha input and output replacement.
                s += command(**opts, alpha=' '.join(reversed(alphas)), prefix='g')+'\n'+mapping_checks(alphas, 'g')
                for sense in ['pess', 'opt']:
                    for v, w in zip(names(alphas)[sense], names(alphas, 'g')[sense]):
                        s += f'assert {v}=={w}\n'
                s += 'assert unrelated=="preserve me"\nexport delimited using "results.csv", replace\n'
                run(f'grid_{rts}_{ori}_{mode}_{grid}', s,
                    check=lambda d, a=alphas, dim=dim, ori=ori, mode=mode, rts=rts: check_lp(d, a, dim, ori, mode, rts))

    if args.baseline_ado:
        # A separate old process calculates each close alpha individually. The
        # new grid must match those unrounded calculations exactly, not merely
        # agree with an oracle at the rounded suffix value.
        for rts in SPECS:
            for ori in ['input', 'output']:
                for grid, spec, alphas in [('close', '.004 .001 .501 .499 1 .999', ['.001', '.004', '.499', '.501', '.999', '1'])]+[(f'legacy{i}', s, a) for i, (s, a) in enumerate(LEGACY)]:
                    key = f'compat_{rts}_{ori}_{grid}'
                    opts = dict(dim=2, ori=ori, mode='full', rts=rts)
                    if grid == 'close':
                        old = DATA
                        for i, a in enumerate(alphas):
                            old += command(**opts, alpha=a, prefix=f'c{i}')+'\n'
                            for sense in ['pess', 'opt']:
                                old += f'rename {names([a],f"c{i}")[sense][0]} {names(alphas)[sense][i]}\n'
                    else: old = DATA+command(**opts, alpha=spec)+'\n'
                    run(key+'_before', old+'keep f_*\nsave "scores.dta", replace\n', source=args.baseline_ado)
                    s = DATA+'gen byte selected=1\n'+command(**opts, alpha=spec)+'\n'+mapping_checks(alphas)
                    s += f'cf f_* using "{out/(key+"_before")/"scores.dta"}", all\n'
                    run(key, s)

    # Full dataset/metadata must survive invalid requests and late name conflicts.
    for i, duplicate in enumerate(['.001 .001', '.004 .001 .004', '0 0', '1 1', '.001 .00100000000001', '.1 .10', '0(.001).004 .002']):
        for replace in ['', 'replacevars']:
            run(f'duplicate_{i}_{bool(replace)}', DATA+'gen str8 f_pess_a0_1="keep"\n'+unchanged(command(alpha=duplicate, extra=replace), 198))
    for i, (bad, rc) in enumerate([('-.001 .004',125), ('0 1.0001',125), ('0 .',127), ('word',121)]):
        run(f'invalid_{i}', DATA+unchanged(command(alpha=bad), rc))
    for protected in ['x1', 'y1', 'lo_x1', 'hi_x1', 'lo_y1', 'hi_y1', 'id']:
        s = DATA+f'rename {protected} f_opt_a0_2\n'
        cmd = command()
        cmd = re.sub(r'\b'+protected+r'\b', 'f_opt_a0_2', cmd)
        run('protected_'+protected, s+unchanged(cmd+' replacevars'))
    for target in ['f_pess_a0_1', 'f_opt_a0_2']:
        run('existing_'+target, DATA+f'gen str8 {target}="keep"\n'+unchanged(command()))
    for i, (prefix, alpha, length) in enumerate([('p'*22, '.004 .001', 32), ('p'*20, '.999 1', 32), ('p'*21, '.001 .002 .003 .004 .0001 .0002 .0003 .0004 .0005 .0006', 32)]):
        run(f'prefix_ok_{i}', DATA+'gen byte selected=1\n'+command(prefix=prefix, alpha=alpha)+'\nassert r(n_alpha)>=2\n')
        run(f'prefix_long_{i}', DATA+unchanged(command(prefix=prefix+'p', alpha=alpha, extra='replacevars'), 7))
    alphas = ['.001', '.004']
    s = DATA+'gen byte selected=id>=2 & id<=7\n'
    for sense in ['pess', 'opt']:
        for v in names(alphas)[sense]: s += f'gen str8 {v}="old"\n'
    s += 'gen f_pess_a0=991\ngen f_opt_a0=992\ngen f_other=993\n'
    s += command(sample='if selected', extra='replacevars', dmu='')+'\n'+mapping_checks(alphas)
    s += 'assert _dmu=="dmu"+string(_n) if selected\nassert _dmu=="" if !selected\nassert f_pess_a0==991 & f_opt_a0==992 & f_other==993\n'
    run('replace_subset_autoid', s)
    run('autoid_collision', DATA+'gen _dmu=991\n'+unchanged(command(dmu='', extra='replacevars')))
    # Automatic-ID preflight and invalid bounds occur before any output is published.
    run('bad_bounds_atomic', DATA+'replace lo_x1=x1+1 in 1\ngen f_pess_a0_1=991\n'+unchanged(command(extra='replacevars'), 459))
    run('empty_sample_atomic', DATA+unchanged(command(sample='if id<0', extra='replacevars'), 2000))
    run('minimum_abbreviations', DATA+'gen byte selected=1\n'
        + 'dea_fuzzy if id>0 in 1/8, in(x1) out(y1) xlo(lo_x1) xhi(hi_x1) ylo(lo_y1) yhi(hi_y1) alp(.004 .001) r(crs) or(input) dmu(id) gen(f) tol(1e-10) det debug\n'
        + mapping_checks(['.001', '.004']))
    run('dense_1001', DATA+'keep in 1\n'+command(mode='crisp', alpha='0(.000001).001')+'\n'+"""
assert r(n_alpha)==1001
local pv "`r(pessvars)'"
local ov "`r(optvars)'"
assert `:word count `pv''==1001 & `:word count `ov''==1001
assert "`:word 1 of `pv''"=="f_pess_a0_1"
assert "`:word 1001 of `pv''"=="f_pess_a0_1001"
foreach v of varlist f_* {
 assert abs(`v'-1)<1e-8
}
""")
    summary = dict(total=len(rows), passed=sum(r['passed'] for r in rows), lp_count=sum(r.get('lp_count',0) for r in rows), max_oracle_gap=max((r.get('max_oracle_gap',0) for r in rows), default=0))
    (out/'summary.json').write_text(json.dumps(summary, indent=2)); print(json.dumps(summary), flush=True)
    raise SystemExit(summary['total'] != summary['passed'])


if __name__ == '__main__': main()
