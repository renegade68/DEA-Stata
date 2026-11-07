#!/usr/bin/env python3
"""Compare qualified DEA calls with physically restricted estimation datasets.

Requires Python 3, Git and Stata 16+. No third-party Python packages are needed.
Every case runs in a fresh Stata process with isolated ado paths. The fixture
includes an excluded, highly productive DMU so reference-set leakage matters.
"""
from _support import package_files, require_git_history
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time

from test_output_protection import unchanged, fuzzy_historical_check, games_historical_check

REPO = Path(__file__).resolve().parents[1]
DATA = '''input id x1 x2 y1 y2
1 .05 .05 80 100
2 2 3 2 1
3 3 2 3 2
4 4 5 3 4
5 2 4 1 2
6 5 3 4 3
7 6 5 5 4
8 3 6 2 5
end
gen byte group=mod(id,2)==0
gen str20 tag=cond(group,"keep = this","omit")
gen str20 sid="Unit "+string(id)
gen double unused=.
gen double xl1=.9*x1
gen double xl2=.9*x2
gen double xh1=1.1*x1
gen double xh2=1.1*x2
gen double yl1=.9*y1
gen double yl2=.9*y2
gen double yh1=1.1*y1
gen double yh2=1.1*y2
'''
PROFILES = {
    'bootstrap': ('reps(50) seed(731)', ['res', 'res_bc', 'res_se', 'res_lo', 'res_hi'],
                  ['n', 'reps', 'level', 'mean_eff', 'mean_bc', 'mean_se'], ['orientation', 'rts']),
    'cooperative': ('', ['res_indiv', 'res_shapley', 'res_contrib', 'res_rank'],
                    ['n', 'v_grand'], ['rts']),
    'crossgame': ('', ['res_agg', 'res_ben', 'res_nash'], ['n'], ['strategy']),
    'fuzzy': ('alpha(0 .5 1)', [f'res_{s}_a{a}' for a in [0, 50, 100] for s in ['pess', 'opt']],
              ['N', 'm_in', 'm_out'], ['cmd', 'orientation', 'rts', 'alpha', 'generate']),
}
BOUNDS = 'xlo(xl1 xl2) xhi(xh1 xh2) ylo(yl1 yl2) yhi(yh1 yh2)'


def command(key, form='B', qualifier='', options='', dmu='id', inputs='x1 x2', outputs='y1 y2'):
    default = PROFILES[key][0]
    io = f'{inputs} = {outputs} {qualifier},' if form == 'A' else f'{qualifier}, inputs({inputs}) outputs({outputs})'
    return f'dea_{key} {io} {default} {options} {"dmu("+dmu+")" if dmu else ""} gen(res)'


def pair(key, form='B', qualifier='if group', selection='group', options='', setup='',
         required='x1 x2 y1 y2', dmu='id', outputs=None, inputs='x1 x2', rhs='y1 y2'):
    outs = outputs or PROFILES[key][1]
    scalars, macros = PROFILES[key][2:]
    checks = ''.join(f'local expected_{s} = r({s})\n' for s in scalars)
    checks += ''.join(f'local expected_{s} `"`r({s})\'"\'\n' for s in macros)
    verify = ''.join(f'assert reldif(r({s}),`expected_{s}\')<1e-12\n' for s in scalars)
    verify += ''.join(f'assert `"`r({s})\'"\'==`"`expected_{s}\'"\'\n' for s in macros)
    nscalar = scalars[0]
    c0 = command(key, form, '', options, dmu, inputs, rhs)
    c1 = command(key, form, qualifier, options, dmu, inputs, rhs)
    body = DATA + setup + f'gen byte selected=(({selection})!=0)\n'
    body += ''.join(f'replace selected=0 if missing({v})\n' for v in required.split())
    if key in ['bootstrap','cooperative'] and dmu:
        body += f'replace selected=0 if missing({dmu})\n'
    body += '''quietly ds
local original_vars "`r(varlist)'"
tempfile original expected
save `original'
count if selected
local selected_n=r(N)
'''
    peer_only = key=='crossgame' and 'strategy(nash)' not in options
    if peer_only:
        # R04: a singleton has no peer-excluding average. Reject it before
        # replacing any result instead of accepting an all-missing score.
        body += "if `selected_n'==1 {\n" + unchanged(c1,2001) + '\n}\nelse {\n'
    body += 'preserve\nkeep if selected\n'
    body += c0 + '\n' + checks + f'assert r({nscalar})==`selected_n\'\n'
    body += 'keep id ' + ' '.join(outs) + '\nsave `expected\'\nrestore\n'
    body += c1 + '\n' + verify + f'assert r({nscalar})==`selected_n\'\n'
    body += ''.join(f'assert missing({v}) if !selected\n' for v in outs)
    for v in outs:
        body += f'count if selected & !missing({v})\nassert r(N)>0\n'
    if key == 'fuzzy' and not dmu:
        body += 'assert _dmu=="dmu"+string(_n) if selected\nassert _dmu=="" if !selected\n'
    body += 'preserve\nkeep if selected\ncf id ' + ' '.join(outs) + ' using `expected\', all\nrestore\n'
    # Compare all original data in original observation order, except explicitly replaced outputs.
    body += 'local replaced "' + ' '.join(outs) + '"\nlocal original_vars : list original_vars - replaced\n'
    body += 'cf `original_vars\' using `original\', all\n'
    if peer_only: body += '}\n'
    return body


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata', required=True)
    ap.add_argument('--output', required=True, type=Path)
    ap.add_argument('--baseline-ado', type=Path, help='Optional source snapshot; otherwise use Git HEAD')
    ap.add_argument('--only')
    args = ap.parse_args()
    if not args.baseline_ado:
        try:
            require_git_history(REPO, 'HEAD')
        except ValueError as error:
            ap.error(str(error))
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    pkg = out/'current_ado'; pkg.mkdir()
    names = package_files(REPO)
    hashes = {}
    for name in names:
        data = (REPO/name).read_bytes(); (pkg/Path(name).name).write_bytes(data)
        hashes[name] = hashlib.sha256(data).hexdigest()
    baseline = out/'baseline_ado'; baseline.mkdir()
    if args.baseline_ado:
        for src in args.baseline_ado.glob('*'):
            if src.is_file(): (baseline/src.name).write_bytes(src.read_bytes())
        baseline_label = str(args.baseline_ado.resolve())
    else:
        baseline_label = subprocess.check_output(['git', '-C', str(REPO), 'rev-parse', 'HEAD'], text=True).strip()
        names = subprocess.check_output(['git', '-C', str(REPO), 'ls-tree', '-r', '--name-only', 'HEAD', 'ado/'], text=True).splitlines()
        for name in names:
            (baseline/Path(name).name).write_bytes(subprocess.check_output(['git', '-C', str(REPO), 'show', f'HEAD:{name}']))
    baseline_hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in baseline.iterdir()}
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes, baseline=baseline_label,
                                                   baseline_sha256=baseline_hashes, stata=args.stata), indent=2))
    results = []

    def run(name, body, source=pkg, force=False):
        if not force and args.only and not re.search(args.only, name): return None
        work = out/name; work.mkdir()
        head = 'version 16\nclear all\nset more off\nset linesize 120\n'
        for folder in ['plus', 'personal', 'site', 'oldplace']:
            (work/folder).mkdir(); head += f'sysdir set {folder.upper()} "{work/folder}"\n'
        head += f'adopath ++ "{source}"\n'
        (work/'test.do').write_text(head + body + '\ndi "TEST_PASS"\n')
        start = time.monotonic()
        try:
            proc = subprocess.run([args.stata, '-b', 'do', 'test.do'], cwd=work, capture_output=True, timeout=180)
            log = (work/'test.log').read_text(errors='replace') if (work/'test.log').exists() else proc.stderr.decode(errors='replace')
            row = dict(test=name, passed='TEST_PASS' in log.splitlines(), seconds=round(time.monotonic()-start, 3))
            if not row['passed']: row['tail'] = log[-2200:]
        except subprocess.TimeoutExpired:
            row = dict(test=name, passed=False, error='timeout')
        results.append(row); (out/'results.json').write_text(json.dumps(results, indent=2))
        print(json.dumps(row), flush=True)
        return work

    qualifiers = {
        'boolean': ('if group', 'group'),
        'equality': ('if group==1', 'group==1'),
        'inequality': ('if group!=0', 'group!=0'),
        'string': ('if tag=="keep = this"', 'tag=="keep = this"'),
        'function': ('if inlist(id,2,4,6,8)', 'inlist(id,2,4,6,8)'),
        'nested': ('if inrange(id,2,8) & (group==1 | id==3)', 'inrange(id,2,8) & (group==1 | id==3)'),
        'in': ('in 2/6', '_n>=2 & _n<=6'),
        'if_in': ('if group==1 in 2/7', 'group==1 & _n>=2 & _n<=7'),
        'in_last': ('in -4/l', '_n>=5'),
        'paren': ('if(group==1)', 'group==1'),
        'all': ('if id>=1 in 1/l', '1'),
        'less_equal': ('if id<=6 & group~=0', 'id<=6 & group~=0'),
        'one': ('in 4/4', '_n==4'),
        'two': ('if inlist(id,4,8)', 'inlist(id,4,8)'),
    }
    for key in PROFILES:
        forms = ['B'] if key == 'fuzzy' else ['A', 'B']
        for form in forms:
            for label, (qualifier, selection) in qualifiers.items():
                run(f'qual_{key}_{form}_{label}', pair(key, form, qualifier, selection))
            for label, setup, required in [
                ('missing_input', 'replace x1=. in 4\n', 'x1 x2 y1 y2'),
                ('missing_output', 'replace y2=. in 6\n', 'x1 x2 y1 y2'),
            ]:
                run(f'{label}_{key}_{form}', pair(key, form, setup=setup, required=required))
            run(f'string_id_{key}_{form}', pair(key, form, dmu='sid', options='detail'))
            run(f'auto_id_{key}_{form}', pair(key, form, dmu=''))
            run(f'missing_condition_{key}_{form}', pair(key, form, setup='replace group=. in 4\n'))
            run(f'missing_equality_{key}_{form}', pair(key, form, 'if group==.', 'group==.', setup='replace group=. in 4\n'))
            for varlist, label in [('x1-x2', 'range'), ('x?', 'wildcard')]:
                run(f'varlist_{key}_{form}_{label}', pair(key, form, inputs=varlist, rhs='y1-y2'))
            first = PROFILES[key][1][0]
            run(f'replace_qualifier_{key}_{form}', pair(key, form, f'if {first}==1 in 2/8', 'group',
                                                       setup=f'gen {first}=group\n', options='replacevars'))
        for label, qualifier in [('if', 'if 0'), ('if_in', 'if group==0 in 2/2')]:
            for replace in [False, True]:
                setup = ''.join(f'gen str12 {v}="old value"\n' for v in PROFILES[key][1]) if replace else ''
                cmd = command(key, qualifier=qualifier, options='replacevars' if replace else '', dmu='')
                run(f'empty_{key}_{label}_{replace}', DATA + setup + unchanged(cmd, 2000))
        run(f'all_missing_{key}', DATA + 'replace x1=.\n' + unchanged(command(key), 2000))
        # Baseline controls preserve qualifier-free numerical behavior exactly.
        compare_name = f'unchanged_full_{key}'
        baseline_form = 'A' if key == 'bootstrap' else 'B'
        control = run(f'baseline_full_{key}', DATA + command(key, baseline_form) + '\nkeep id res*\nsave "expected.dta"\n',
                      baseline, force=bool(args.only and re.search(args.only, compare_name)))
        if control:
            # F03 repairs the smoother and CI quantiles. Preserve point scores as
            # the historical bootstrap control; inference has its own audit.
            compare_outputs = 'res' if key == 'bootstrap' else ' '.join(PROFILES[key][1])
            comparison=f'cf id {compare_outputs} using "{control}/expected.dta", all\n'
            if key=='bootstrap':
                # Refactorization changes final bits in the old point scores.
                # DMU 1 dominates this CRS fixture: verify that frontier as
                # well as a tight historical comparison. Qualified-vs-restricted
                # current fits above retain their exact dataset comparisons.
                comparison=f'''cf id using "{control}/expected.dta", all
tempname observed historical
mkmat res, matrix(`observed')
preserve
use "{control}/expected.dta", clear
mkmat res, matrix(`historical')
restore
mata: assert(!hasmissing(st_matrix("`observed'")) & !hasmissing(st_matrix("`historical'")))
mata: assert(max(abs(st_matrix("`observed'")-st_matrix("`historical'")))<1e-12)
tempvar intensity expected
gen double `intensity'=max(y1/y1[1],y2/y2[1])
assert x1[1]*`intensity'<=x1 & x2[1]*`intensity'<=x2
gen double `expected'=max(x1[1]*`intensity'/x1,x2[1]*`intensity'/x2)
assert abs(res-`expected')<1e-12
drop `intensity' `expected'
'''
            if key == 'cooperative':
                comparison=games_historical_check(key,PROFILES[key][1],control/'expected.dta').replace('assert note=="keep me"\n','')
            if key == 'fuzzy':
                comparison=fuzzy_historical_check(PROFILES[key][1],control/'expected.dta')
            run(compare_name, DATA + command(key) + '\n'+comparison)
        nscalar = PROFILES[key][2][0]
        first = PROFILES[key][1][0]
        run(f'old_ignored_if_{key}', DATA + command(key, baseline_form, qualifier='if group') +
            f'\nassert r({nscalar})==8\nassert !missing({first}) if !group\n', baseline)
        if key != 'fuzzy':
            run(f'old_equality_parse_{key}', DATA + unchanged(command(key, qualifier='if group==1'), 198), baseline)

    for rts in ['crs', 'vrs', 'irs', 'drs']:
        for orient in ['input', 'output']:
            for form in ['A', 'B']:
                run(f'options_bootstrap_{rts}_{orient}_{form}', pair('bootstrap', form, options=f'rts({rts}) ort({orient}) level(90) detail'))
            for fuzzy in [False, True]:
                run(f'options_fuzzy_{rts}_{orient}_{fuzzy}', pair('fuzzy', options=f'rts({rts}) orientation({orient}) detail '+(BOUNDS if fuzzy else '')))
    for rts in ['crs', 'vrs']:
        for form in ['A', 'B']:
            run(f'options_cooperative_{rts}_{form}', pair('cooperative', form, options=f'rts({rts}) detail'))
    for strategy, suffix in [('aggressive', 'agg'), ('benevolent', 'ben'), ('nash', 'nash'), ('all', '')]:
        for form in ['A', 'B']:
            run(f'options_crossgame_{strategy}_{form}', pair('crossgame', form, options=f'strategy({strategy}) matrix detail maxiter(30)',
                                                           outputs=[f'res_{suffix}'] if suffix else None))
    for bound in ['xl1', 'xh2', 'yl2', 'yh1']:
        run(f'missing_fuzzy_{bound}', pair('fuzzy', options=BOUNDS, setup=f'replace {bound}=. in 4\n', required=f'x1 x2 y1 y2 {bound}'))
    run('missing_bootstrap_id', pair('bootstrap', setup='replace id=. in 4\n'))
    # Additive scores have no n=20 cap; selected-sample control and old limit remain.
    large = DATA + 'set obs 25\nforeach v in id x1 x2 y1 y2 {\nreplace `v\'=1 if _n>8\n}\n'
    run('cooperative_limit_selected', large + command('cooperative', qualifier='in 2/5') + '\nassert r(n)==4\nassert missing(res_indiv) if _n<2 | _n>5\n')
    run('cooperative_additive_large',large+command('cooperative')+'\nassert r(n)==25 & r(n_lp)==25\nassert res_indiv==res_shapley & res_contrib==0\n')
    run('missing_cooperative_id',pair('cooperative',setup='replace id=. in 4\n'))
    run('missing_cooperative_string_id',pair('cooperative',dmu='sid',setup='replace sid="" in 4\n'))
    run('old_cooperative_limit', large + unchanged(command('cooperative', qualifier='in 2/5'), 198), baseline)
    run('known_baseline_crossgame_one', DATA + 'keep in 4\n' + command('crossgame') +
        '\nassert r(n)==1\nassert missing(res_agg) & missing(res_ben)\nassert abs(res_nash-1)<1e-6\n', baseline)
    # Compact Form A separator and bootstrap's output-varlist Form B.
    for key in ['bootstrap', 'cooperative', 'crossgame']:
        run(f'compact_{key}', DATA + command(key, 'A', 'if group==1').replace('x1 x2 = y1 y2', 'x1 x2=y1 y2') + '\nassert r(n)==4\n')
        for form in ['A', 'B']:
            run(f'tabs_{key}_{form}', DATA + command(key, form, 'if\tgroup!=0', inputs='x1\tx2') + '\nassert r(n)==4\n')
        run(f'no_comma_{key}', DATA + f'dea_{key} x1 x2 = y1 y2 if group==1\nassert r(n)==4\n')
    run('bootstrap_output_varlist', pair('bootstrap').replace('dea_bootstrap if group, inputs(x1 x2) outputs(y1 y2)',
                                                            'dea_bootstrap y1 y2 if group, inputs(x1 x2)'))
    run('bootstrap_requires_outputs', DATA + unchanged('dea_bootstrap, inputs(x1 x2) reps(50)', 198))
    for key in ['bootstrap', 'cooperative', 'crossgame']:
        run(f'invalid_if_{key}', DATA + unchanged(command(key, 'A', 'if unknown_flag==1'), 111))
    # Independent one-input CRS oracle confirms the excluded frontier DMU is absent.
    for key, first in [('bootstrap', 'res'), ('cooperative', 'res_indiv'), ('fuzzy', 'res_pess_a100')]:
        run(f'analytic_frontier_{key}', DATA + command(key, qualifier='if group', inputs='x1', outputs='y1') +
            f'\nassert abs({first} - (y1/x1))<1e-7 if group\nassert missing({first}) if !group\n')
    if not results: raise SystemExit('No cases matched --only')
    failures = sum(not row['passed'] for row in results)
    print(json.dumps(dict(total=len(results), passed=len(results)-failures, failed=failures)), flush=True)
    raise SystemExit(bool(failures))


if __name__ == '__main__':
    main()
