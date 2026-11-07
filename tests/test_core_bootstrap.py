#!/usr/bin/env python3
"""Validate core bootstrap integration, sample mapping, outputs and rejected options."""
from _support import package_files
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time

from test_output_protection import unchanged

REPO = Path(__file__).resolve().parents[1]
DATA = '''input id x y
1 .1 80
2 1 1
3 2 3
4 4 4
5 5 3
6 3 2
7 2 2
8 6 4
end
gen byte selected=id>=2 & id<=7
gen str20 sid="Unit "+string(id)
gen double unused=.
gen double z=id/10
gen double w=1
gen double lo=.9*x
gen double hi=1.1*x
gen double bad=1
'''
OUTS = ['eff', 'eff_bc', 'eff_bias', 'eff_se', 'eff_lo', 'eff_hi']


def command(form='A', qualifier='if selected', options=''):
    io = f'x = y {qualifier},' if form == 'A' else f'y {qualifier}, inputs(x)'
    return f'dea {io} dmu(id) gen(eff) bootstrap reps1(60) {options}'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata', required=True)
    ap.add_argument('--output', required=True, type=Path)
    ap.add_argument('--baseline-ado', required=True, type=Path)
    ap.add_argument('--only')
    args = ap.parse_args(); out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    pkg = out/'current_ado'; pkg.mkdir()
    names = package_files(REPO)
    hashes = {}
    for name in names:
        data = (REPO/name).read_bytes(); (pkg/Path(name).name).write_bytes(data)
        hashes[name] = hashlib.sha256(data).hexdigest()
    baseline = out/'baseline_ado'; baseline.mkdir()
    for path in args.baseline_ado.glob('*'):
        if path.is_file(): (baseline/path.name).write_bytes(path.read_bytes())
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes, baseline=str(args.baseline_ado.resolve()),
        baseline_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in baseline.iterdir()}, stata=args.stata), indent=2))
    results = []

    def run(name, body, source=pkg):
        if args.only and not re.search(args.only, name): return
        work = out/name; work.mkdir()
        head = 'version 16\nclear all\nset more off\nset linesize 120\n'
        for folder in ['plus', 'personal', 'site', 'oldplace']:
            (work/folder).mkdir(); head += f'sysdir set {folder.upper()} "{work/folder}"\n'
        head += f'adopath ++ "{source}"\n'
        (work/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n')
        start = time.monotonic()
        try:
            p = subprocess.run([args.stata, '-b', 'do', 'test.do'], cwd=work, capture_output=True, timeout=180)
            log = (work/'test.log').read_text(errors='replace') if (work/'test.log').exists() else p.stderr.decode(errors='replace')
            row = dict(test=name, passed='TEST_PASS' in log.splitlines(), seconds=round(time.monotonic()-start, 3))
            if not row['passed']: row['tail'] = log[-2200:]
        except subprocess.TimeoutExpired:
            row = dict(test=name, passed=False, error='timeout')
        results.append(row); (out/'results.json').write_text(json.dumps(results, indent=2)); print(json.dumps(row), flush=True)

    for form in ['A', 'B']:
        for engine in ['simplex', 'ipm', 'auto']:
            for rts in ['crs', 'vrs', 'nirs', 'ndrs']:
                for ort in ['input', 'output']:
                    options = f'engine({engine}) rts({rts}) ort({ort}) tolerance(1e-8) reference saveref'
                    cmd = command(form, options=options)
                    body = DATA + 'set seed 871\n' + cmd + '''
assert r(N)==6 & r(bootstrap_reps)==60 & r(bootstrap_level)==95
assert "`r(bootstrap_engine)'"=="simplex"
local refs "`r(refvars)'"
matrix first_L=r(lambda)
matrix B=r(bootstrap)
mata: assert(max(abs(st_matrix("B")-st_data(.,tokens("eff_bc eff_bias eff_se eff_lo eff_hi"),"selected")))<1e-12)
assert abs(eff-eff_bc-eff_bias)<1e-12 if selected
assert eff_lo<=eff_hi if selected
foreach v of varlist eff* {
    assert missing(`v') if !selected
    assert !missing(`v') if selected
}
tempfile results
save `results'
'''
                    body += cmd.replace(' bootstrap reps1(60)', '') + ' replacevars\n'
                    body += '''mata: assert(max(abs(st_matrix("first_L")-st_matrix("r(lambda)")))<1e-12)
cf eff `refs' using `results', all
set seed 871
'''
                    normalized = {'nirs':'drs', 'ndrs':'irs'}.get(rts, rts)
                    body += f'dea_bootstrap if selected, inputs(x) outputs(y) dmu(id) rts({normalized}) ort({ort}) reps(60) tolerance(1e-8) gen(bs)\n'
                    tol = '2e-5' if engine == 'ipm' else '1e-7'
                    for suffix in ['', '_bc', '_se', '_lo', '_hi']:
                        body += f'assert abs(eff{suffix}-bs{suffix})<{tol} if selected\n'
                    run(f'shared_{form}_{engine}_{rts}_{ort}', body)

    qualifiers = [('if', 'if selected'), ('eq', 'if selected==1'), ('string', 'if sid!="Unit 1" & id<8'),
                  ('in', 'in 2/7'), ('combined', 'if inrange(id,2,7) in 2/8')]
    for form in ['A', 'B']:
        for label, qualifier in qualifiers:
            for ort in ['input', 'output']:
                cmd = command(form, qualifier, f'ort({ort})')
                body = DATA + '''tempfile original expected
save `original'
preserve
keep if selected
set seed 331
''' + command(form, '', f'ort({ort})') + '''
keep id eff*
save `expected'
restore
set seed 331
''' + cmd + '''
preserve
keep if selected
cf id eff* using `expected', all
restore
cf id x y sid z w lo hi bad selected unused using `original', all
'''
                run(f'sample_{form}_{label}_{ort}', body)

    for form in ['A', 'B']:
        for option in ['envars(z) reps2(1)', 'reps2(2000)', 'reps2(10)', 'reps2(-1)', 'reps1(0)', 'reps1(49)',
                       'level(0)', 'level(100)', 'level(.)', 'tolerance(0)', 'tolerance(-1)', 'tolerance(.)']:
            cmd = command(form).replace('reps1(60)', '') if option.startswith('reps1') else command(form)
            name = re.sub(r'\W+', '_', option).strip('_')
            run(f'reject_option_{form}_{name}', DATA + 'gen str12 eff="old data"\n' + unchanged(cmd+' '+option+' replacevars', 198))
        for model in ['fdh', 'ddf', 'undesirable', 'cost', 'revenue', 'profit', 'allocative', 'sbm',
                      'supereff', 'supersbm', 'additive', 'multiplier', 'imprecise', 'congestion', 'sbmundesirable']:
            run(f'reject_model_{form}_{model}', DATA + unchanged(command(form, options=f'model({model}) w(w) p(w) gx(x) gy(y) bad(bad)'), 198))
        for opt in ['fdh', 'bounded(x)', 'ordinal(x)', 'xlo(lo) xhi(hi)', 'xord(x)', 'yord(y)']:
            name = re.sub(r'\W+', '_', opt).strip('_')
            run(f'reject_imprecise_{form}_{name}', DATA + unchanged(command(form, options=opt), 198))
        for target in OUTS:
            run(f'reject_output_{form}_{target}', DATA + f'gen str12 {target}="keep"\n' + unchanged(command(form), 110))
            run(f'replace_output_{form}_{target}', DATA + f'gen str12 {target}="replace"\n' + command(form, options='replacevars') + '\nassert !missing(eff_hi) if selected\n')
        run(f'protect_ci_input_{form}', DATA + 'rename x eff_hi\n' + unchanged(command(form).replace('x', 'eff_hi')+' replacevars', 110))
        run(f'replace_if_{form}', DATA + 'gen eff_lo=selected\n' + command(form, 'if eff_lo==1', 'replacevars') + '\nassert missing(eff_lo) if !selected\n')
        run(f'empty_{form}', DATA + 'gen eff=991\n' + unchanged(command(form, 'if 0', 'replacevars'), 2000))

    for level in [1, 90, 95, 99, 99.5]:
        body = DATA + 'set seed 838\n' + command(options=f'level({level})').replace('reps1(60)', 'reps1(51)') + '\nassert eff_lo<=eff_hi if selected\n'
        run(f'level_order_{level}', body)
    run('level_changes_only_ci', DATA + 'set seed 839\n' + command(options='level(90)') + '''
local state1 `"`c(rngstate)'"'
foreach v of varlist eff* {
    rename `v' old_`v'
}
set seed 839
''' + command(options='level(99)') + '''
assert `"`state1'"'==`"`c(rngstate)'"'
assert eff==old_eff & eff_bc==old_eff_bc & eff_bias==old_eff_bias & eff_se==old_eff_se
assert eff_lo<=old_eff_lo & eff_hi>=old_eff_hi if selected
count if selected & (eff_lo<old_eff_lo | eff_hi>old_eff_hi)
assert r(N)>0
''')
    run('reps_changes_draws', DATA + 'set seed 839\n' + command() + '''
local state1 `"`c(rngstate)'"'
rename eff_bc old_bc
set seed 839
''' + command(options='replacevars').replace('reps1(60)', 'reps1(80)') + '''
assert `"`state1'"'!=`"`c(rngstate)'"'
count if selected & abs(eff_bc-old_bc)>1e-10
assert r(N)>0
''')
    run('saving_all_outputs', DATA + command(options='reference saveref saving("results.dta")') + '''
local refs "`r(refvars)'"
matrix B=r(bootstrap)
matrix L=r(lambda)
use "results.dta", clear
assert _N==6
confirm variable eff eff_bc eff_bias eff_se eff_lo eff_hi
mata: assert(max(abs(st_matrix("B")-st_data(.,tokens("eff_bc eff_bias eff_se eff_lo eff_hi"))))<1e-12)
mata: assert(max(abs(st_matrix("L")-st_data(.,tokens(st_local("refs")))))<1e-12)
''')
    for dmu in ['sid', '']:
        cmd = command().replace('dmu(id)', f'dmu({dmu})' if dmu else '')
        run(f'identifier_{dmu or "auto"}', DATA + cmd + '\nassert !missing(eff_hi) if selected\n')
    run('missing_input', DATA + 'replace x=. in 4\n' + command() + '\nassert r(N)==5\nassert missing(eff_bc) in 4\n')
    run('output_not_all_one', DATA + 'set seed 839\n' + command(options='ort(output)') + '\ncount if selected & eff_bc<.9\nassert r(N)>0\n')
    run('dedicated_then_core', DATA + '''set seed 871
dea_bootstrap if selected, inputs(x) outputs(y) dmu(id) reps(60) gen(bs) tolerance(1e-8)
set seed 871
''' + command(options='tolerance(1e-8)') + '''
foreach suffix in bc se lo hi {
    assert abs(eff_`suffix'-bs_`suffix')<1e-7 if selected
}
''')
    for opt in ['envars(z) reps2(1)', 'reps2(2000)']:
        run('reject_without_bootstrap_'+opt.split('(')[0], DATA + unchanged('dea x=y, dmu(id) '+opt, 198))

    # Reflection must be performed AFTER variance correction. The old ordering
    # has a deterministic upper support bound below 1 for a nonconstant sample.
    sampler = '''set seed 771
mata:
theta=(.4\\.5\\.6\\.7\\.8\\.9\\1)
h=.2
cap=mean(theta)+(1-mean(theta))/sqrt(1+h^2/variance(theta))
draws=J(7,4000,.)
for (b=1;b<=cols(draws);b++) draws[.,b]=_sw_sample(theta,h)
assert(min(draws)>0 & max(draws)<=1)
'''
    run('sampler_reaches_boundary', '''quietly findfile _dea_bootstrap_core.ado
quietly run `"`r(fn)'"'
''' + sampler + 'assert(max(draws)>cap+.01 & max(draws)>.999)\nend\n')
    run('old_sampler_excludes_boundary', '''quietly findfile dea_bootstrap.ado
quietly run `"`r(fn)'"'
''' + sampler + 'assert(max(draws)<=cap)\nend\n', baseline)

    # Defects must still be observable in the frozen pre-F03 source.
    run('old_output_all_one', DATA + 'set seed 839\n' + command(options='ort(output)') + '\nassert eff_bc==1 if selected\n', baseline)
    run('old_sample_misalignment', DATA + 'set seed 839\n' + command() + '''
keep if selected
rename eff_bc old_bc
set seed 839
''' + command(qualifier='', options='replacevars') + '''
count if abs(eff_bc-old_bc)>1e-4
assert r(N)>0
''', baseline)
    run('old_inert_level_reps2', '''set obs 80
set seed 27811
gen id=_n
gen double z=rnormal()
gen double x=1+5*runiform()
gen double y=x/(2+exp(.3*z)+abs(rnormal()))
set seed 839
dea x=y, dmu(id) bootstrap envars(z) reps1(20) reps2(10) level(90) gen(a)
local state `"`c(rngstate)'"'
set seed 839
dea x=y, dmu(id) bootstrap envars(z) reps1(20) reps2(100) level(99) gen(b)
assert a_bc==b_bc & a_bias==b_bias
assert `"`state'"'==`"`c(rngstate)'"'
capture confirm variable a_lo
assert _rc==111
''', baseline)
    if not results: raise SystemExit('No cases matched --only')
    failures = sum(not r['passed'] for r in results)
    print(json.dumps(dict(total=len(results), passed=len(results)-failures, failed=failures)), flush=True)
    raise SystemExit(bool(failures))


if __name__ == '__main__':
    main()
