#!/usr/bin/env python3
"""Reproduce reference persistence/retry defects and test the corrected command.

Requires Python 3, git, and a licensed Stata 16+ batch executable. All generated
packages, datasets and logs go into a new --output directory, never ado/.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time

REPO = Path(__file__).resolve().parents[1]
BASIC = """input id x y
1 1 1
2 2 3
3 4 4
4 5 3
end
"""
RANDOM = """set seed 7351
set obs 12
gen long id=_n
gen double x1=1+10*runiform()
gen double x2=1+10*runiform()
gen double y1=1+5*runiform()
gen double y2=1+5*runiform()
"""
HOOK = """    // =================================================================
    // Display results
    // =================================================================
    // lambda_vars and slack_in/out_vars already live in this scope"""
INJECT = """    // TEST COPY ONLY: invalidate the first attempt; IPM still solves normally.
    if "`engine'" == "simplex" {
        quietly replace `generate' = . if `touse'
        foreach v of local lambda_vars {
            quietly replace `v' = . if `touse'
        }
        foreach v of local slack_in_vars {
            quietly replace `v' = 991 if `touse'
        }
        foreach v of local slack_out_vars {
            quietly replace `v' = 991 if `touse'
        }
        if "`weight_u0_var'" != "" quietly replace `weight_u0_var' = 991 if `touse'
    }

"""


def instrument(source, dest):
    shutil.copytree(source, dest)
    path = dest / 'dea.ado'
    code = path.read_text()
    assert code.count(HOOK) == 1, 'Retry instrumentation anchor changed'
    path.write_text(code.replace(HOOK, INJECT + HOOK))


def checks(inputs='x1 x2', outputs='y1 y2', orientation='input', rts='crs'):
    """Check primal feasibility independently of the storage implementation."""
    inequality = ('L*X :- X:*E' if orientation == 'input' else 'L*X :- X')
    out = ('Y :- L*Y' if orientation == 'input' else 'Y:/E :- L*Y')
    convexity = {
        'crs': '', 'vrs': 'assert(max(abs(rowsum(L):-1))<1e-6)',
        'nirs': 'assert(max(rowsum(L):-1)<1e-6)',
        'ndrs': 'assert(min(rowsum(L):-1)>-1e-6)',
    }[rts]
    return f'''
local refs "`r(refvars)'"
local failures=r(n_failed)
local average=r(mean_eff)
matrix L=r(lambda)
assert `failures'==0
assert !missing(eff)
quietly summarize eff, meanonly
assert abs(r(mean)-`average')<1e-9
mata:
L=st_matrix("L")
X=st_data(.,tokens("{inputs}"))
Y=st_data(.,tokens("{outputs}"))
E=st_data(.,"eff")
assert(sum(L:>=.)==0)
assert(min(L)>-1e-7)
assert(max({inequality})<1e-5)
assert(max({out})<1e-5)
{convexity}
assert(max(abs(st_data(.,tokens(st_local("refs")))-L))<1e-12)
end
foreach v of local refs {{
    assert "`:type `v''"=="double"
}}
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stata', required=True)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--baseline-ref', help='Optional pre-fix git commit for old-defect reproduction and numerical comparisons')
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    package = out / 'current_ado'
    package.mkdir()
    tracked = subprocess.check_output(['git', '-C', str(REPO), 'ls-files', 'ado/'], text=True).splitlines()
    hashes = {}
    for name in tracked:
        data = (REPO/name).read_bytes()
        (package/Path(name).name).write_bytes(data)
        hashes[name] = hashlib.sha256(data).hexdigest()
    forced = out / 'forced_retry_ado'
    instrument(package, forced)
    baseline = None
    if args.baseline_ref:
        baseline = out / 'baseline_ado'
        baseline.mkdir()
        for name in tracked:
            data = subprocess.check_output(['git', '-C', str(REPO), 'show', args.baseline_ref+':'+name])
            (baseline/Path(name).name).write_bytes(data)
        instrument(baseline, out/'baseline_forced_retry_ado')
    (out/'manifest.json').write_text(json.dumps({
        'baseline_ref': args.baseline_ref, 'stata': args.stata,
        'source_sha256': hashes, 'retry_method': 'test-only first-attempt invalidation; unchanged real IPM',
    }, indent=2))
    results = []

    def run(name, code, ado=package):
        work = out/name
        work.mkdir()
        header = 'version 16\nclear all\nset more off\nset linesize 120\n'
        for folder in ['plus', 'personal', 'site', 'oldplace']:
            (work/folder).mkdir()
            header += f'sysdir set {folder.upper()} "{work/folder}"\n'
        header += f'adopath ++ "{ado}"\n'
        (work/'test.do').write_text(header+code+'\ndisplay "TEST_PASS"\n')
        start = time.monotonic()
        try:
            process = subprocess.run([args.stata, '-b', 'do', 'test.do'], cwd=work,
                                     timeout=180, capture_output=True, text=True)
            logpath = work/'test.log'
            log = logpath.read_text(errors='replace') if logpath.exists() else process.stderr
            passed = 'TEST_PASS' in log.splitlines()
            result = dict(test=name, passed=passed, seconds=round(time.monotonic()-start, 3))
            if not passed:
                result['tail'] = log[-2200:]
        except subprocess.TimeoutExpired:
            result = dict(test=name, passed=False, error='timeout')
        results.append(result)
        (out/'results.json').write_text(json.dumps(results, indent=2))
        print(json.dumps(result), flush=True)
        return work

    # A 32-cell grid, plus optional controls from the committed pre-fix source.
    for engine in ['simplex', 'ipm']:
        for rts in ['crs', 'vrs', 'nirs', 'ndrs']:
            for orientation in ['input', 'output']:
                for stage in [1, 2]:
                    name = f'grid_{engine}_{rts}_{orientation}_{stage}'
                    opts = f'dmu(id) engine({engine}) rts({rts}) orientation({orientation}) stage({stage}) generate(eff) reference'
                    setup = RANDOM
                    compare = ''
                    if baseline:
                        control = run('control_'+name, RANDOM+f'dea x1 x2 = y1 y2, {opts}\n'+'''
matrix C=r(lambda)
svmat double C, names(base_lambda)
rename eff base_eff
save "control.dta", replace
''', baseline)
                        setup = f'use "{control}/control.dta", clear\n'
                        compare = '''
assert abs(eff-base_eff)<1e-10
mata: assert(max(abs(st_matrix("L")-st_data(.,"base_lambda":+strofreal(1..st_nobs()))))<1e-10)
'''
                    run(name, setup+f'dea x1 x2 = y1 y2, {opts} saveref saving("final.dta")\n'+
                        checks(orientation=orientation, rts=rts)+compare+'''
preserve
use "final.dta", clear
confirm variable id x1 x2 y1 y2 eff
capture ds __*
assert _rc==111
capture ds base_*
assert _rc==111
mata: assert(max(abs(st_matrix("L")-st_data(.,tokens(st_local("refs")))))<1e-12)
restore
''')

    # Forced retries test lambda, statistics, stale-slack removal and both file modes.
    for rts in ['crs', 'vrs']:
        for orientation in ['input', 'output']:
            for ref in [False, True]:
                for replace in [False, True]:
                    opts = f'dmu(id) rts({rts}) orientation({orientation}) reference'
                    prep = '' if not replace else '''preserve
clear
set obs 1
gen sentinel=991
save "final.dta", replace
restore
'''
                    run(f'retry_{rts}_{orientation}_ref{int(ref)}_replace{int(replace)}', BASIC+f'''
dea x = y, {opts} engine(ipm) generate(control)
matrix C=r(lambda)
local mean0=r(mean_eff)
local min0=r(min_eff)
local max0=r(max_eff)
local count0=r(N_eff)
'''+prep+f'''
dea x = y, {opts} engine(auto) generate(eff) {'saveref' if ref else ''} saving("final.dta") {'replace' if replace else ''}
local refs "`r(refvars)'"
matrix L=r(lambda)
assert "`r(engine)'"=="ipm"
assert r(n_failed)==0
assert abs(r(mean_eff)-`mean0')<1e-10
assert abs(r(min_eff)-`min0')<1e-10
assert abs(r(max_eff)-`max0')<1e-10
assert r(N_eff)==`count0'
assert abs(eff-control)<1e-10
mata: assert(max(abs(st_matrix("L")-st_matrix("C")))<1e-10)
capture matrix list r(slack_in)
assert _rc==111
capture matrix list r(slack_out)
assert _rc==111
use "final.dta", clear
assert !missing(eff)
capture ds __*
assert _rc==111
capture confirm variable sentinel
assert _rc==111
'''+('''mata: assert(max(abs(st_matrix("L")-st_data(.,tokens(st_local("refs")))))<1e-12)
''' if ref else '''capture ds _ref_*
assert _rc==111
'''), forced)

    run('implicit_id_success', BASIC+'''
dea x = y in 2/4, rts(crs) generate(eff) reference saveref saving("final.dta")
assert r(generated_dmu)==1
assert "`r(refvars)'"=="_ref_dmu2 _ref_dmu3 _ref_dmu4"
assert _dmu=="dmu"+string(_n) in 2/4
assert "`:char _ref_dmu2[dea_peer_dmu]'"=="_dmu"
assert missing(eff,_ref_dmu2,_ref_dmu3,_ref_dmu4) in 1
use "final.dta", clear
assert _N==3
confirm variable _dmu
''')
    run('retry_without_reference', BASIC+'''
dea x = y, dmu(id) engine(ipm) generate(control)
local mean0=r(mean_eff)
dea x = y, dmu(id) engine(auto) generate(eff) saving("final.dta")
assert "`r(engine)'"=="ipm"
assert r(n_failed)==0
assert abs(r(mean_eff)-`mean0')<1e-10
assert "`r(refvars)'"==""
capture matrix list r(lambda)
assert _rc==111
assert abs(eff-control)<1e-10
use "final.dta", clear
assert !missing(eff)
capture ds __*
assert _rc==111
''', forced)
    for saveref in ['', 'reference saveref']:
        run('existing_auto_id_'+str(bool(saveref)), BASIC+'''
gen str20 _dmu="original"+string(_n)
datasignature
local before "`r(datasignature)'"
'''+f'capture noisily dea x = y, generate(eff) {saveref}\n'+'''
assert _rc==110
datasignature
assert "`before'"=="`r(datasignature)'"
''')
    for pre, opts in [
        ('gen _ref_dmu1=991', 'generate(eff) reference saveref'),
        ('', 'generate(_dmu) reference saveref'),
        ('', 'generate(_ref_dmu1) reference saveref'),
    ]:
        run('auto_preflight_'+str(len(results)), BASIC+pre+'''
datasignature
local before "`r(datasignature)'"
'''+f'capture noisily dea x = y, {opts}\n'+'''
assert _rc==110
datasignature
assert "`before'"=="`r(datasignature)'"
capture confirm variable _dmu
assert _rc==111
''')
    run('explicit_existing_dmu', BASIC+'''
gen str20 _dmu="original"+string(_n)
dea x = y, dmu(_dmu) generate(eff) reference saveref
assert r(generated_dmu)==0
assert _dmu=="original"+string(_n)
confirm variable _ref_original1 _ref_original2 _ref_original3 _ref_original4
''')

    for opts in ['model(cost) w(x)', 'model(sbm)', 'fdh', 'model(supereff) fdh',
                 'xlo(x) xhi(x)', 'model(imprecise) xlo(x) xhi(x)', 'ordinal(x)']:
        run('unsupported_'+str(len(results)), BASIC+'''
datasignature
local before "`r(datasignature)'"
'''+f'capture noisily dea x = y, generate(eff) reference saveref {opts}\n'+'''
assert _rc==198
datasignature
assert "`before'"=="`r(datasignature)'"
''')

    run('sample_alignment', '''input id x y selected
40 1 1 1
99 2 3 0
30 3 . 1
10 2 3 1
77 5 3 0
20 4 4 1
end
gen long original_order=_n
dea x = y if selected in 1/6, dmu(id) rts(vrs) generate(eff) reference saveref
local refs "`r(refvars)'"
matrix L=r(lambda)
assert "`refs'"=="_ref_40 _ref_10 _ref_20"
assert original_order==_n
foreach v of local refs {
    assert missing(`v') if !selected | missing(y)
}
mata: obs=selectindex(st_data(.,"selected"):==1 :& st_data(.,"y"):<.)
mata: assert(max(abs(st_data(obs,tokens(st_local("refs")))-st_matrix("L")))<1e-12)
''')
    run('name_metadata', '''input str80 id x y
"A B" 1 1
"A-B" 2 3
"A_B_2" 4 4
"abcdefghijklmnopqrstuvwxyz0123456789 A" 5 3
"abcdefghijklmnopqrstuvwxyz0123456789 B" 6 4
"서울" 7 3
"서울" 8 4
end
dea x = y, dmu(id) rts(crs) generate(eff) reference saveref
local refs "`r(refvars)'"
local count: word count `refs'
assert `count'==7
local j=0
foreach v of local refs {
    local ++j
    assert strlen("`v'")<=32
    mata: assert(st_global("`v'[dea_peer_id]")==st_sdata(`j',"id"))
    assert "`:char `v'[dea_peer_index]'"=="`j'"
    assert "`:char `v'[dea_peer_dmu]'"=="id"
}
''')
    run('repeated_reference_collision', BASIC+'''
dea x = y, dmu(id) generate(eff) reference saveref
datasignature
local before "`r(datasignature)'"
capture noisily dea x = y, dmu(id) generate(second) reference saveref replace
assert _rc==110
datasignature
assert "`before'"=="`r(datasignature)'"
''')
    for engine in ['simplex', 'ipm']:
        for orientation in ['input', 'output']:
            run(f'super_{engine}_{orientation}', BASIC+f'''
dea x = y, dmu(id) model(supereff) engine({engine}) orientation({orientation}) generate(eff) reference saveref
'''+checks('x', 'y', orientation)+'''
mata: assert(max(abs(diagonal(st_matrix("L"))))<1e-7)
''')
            run(f'multiplier_{engine}_{orientation}', BASIC+f'''
dea x = y, dmu(id) model(multiplier) rts(vrs) engine({engine}) orientation({orientation}) generate(eff) saveweights saving("final.dta")
matrix U=r(multiplier_u0)
preserve
use "final.dta", clear
confirm variable id x y eff _dea_v1 _dea_u1 _dea_u0
mata: assert(max(abs(st_data(.,"_dea_u0")-st_matrix("U")))<1e-10)
restore
''')
    run('multiplier_retry', BASIC+'''
dea x = y, dmu(id) model(multiplier) rts(vrs) engine(ipm) generate(control) saveweights
matrix C=r(multiplier_u0)
dea x = y, dmu(id) model(multiplier) rts(vrs) engine(auto) generate(eff) saveweights saving("final.dta")
matrix U=r(multiplier_u0)
assert "`r(engine)'"=="ipm"
assert r(n_failed)==0
assert abs(eff-control)<1e-10
mata: assert(max(abs(st_matrix("U")-st_matrix("C")))<1e-10)
use "final.dta", clear
mata: assert(max(abs(st_data(.,"_dea_u0")-st_matrix("U")))<1e-10)
''', forced)
    # VRS super-efficiency can be infeasible after excluding the evaluated DMU.
    # Exercise a natural retry with partial results, as well as explicit engines.
    for engine in ['simplex', 'ipm', 'auto']:
        for orientation in ['input', 'output']:
            run(f'super_vrs_partial_{engine}_{orientation}', BASIC+f'''
dea x = y, dmu(id) model(supereff) rts(vrs) engine(ipm) orientation({orientation}) generate(control) reference
matrix C=r(lambda)
dea x = y, dmu(id) model(supereff) rts(vrs) engine({engine}) orientation({orientation}) generate(eff) reference saveref saving("final.dta")
local refs "`r(refvars)'"
local failures=r(n_failed)
matrix L=r(lambda)
assert missing(eff)==missing(control)
assert abs(eff-control)<1e-6 if !missing(eff)
quietly count if missing(eff)
assert r(N)>0 & r(N)<_N
assert r(N)==`failures'
mata:
S=st_data(.,tokens(st_local("refs")))
L=st_matrix("L")
assert(all((S:==L):|(abs(S-L):<1e-12)))
end
use "final.dta", clear
assert _N==4
quietly count if missing(eff)
assert r(N)==`failures'
''')
    for model in ['fdh', 'imprecise', 'cost', 'sbm', 'additive']:
        extra = {'imprecise':'xlo(x) xhi(x)', 'cost':'w(x)'}.get(model, '')
        run('saving_'+model, BASIC+f'''
dea x = y, dmu(id) model({model}) {extra} generate(control)
dea x = y, dmu(id) model({model}) {extra} generate(eff) saving("final.dta")
assert abs(eff-control)<1e-7 | (missing(eff) & missing(control))
preserve
use "final.dta", clear
confirm variable id x y eff
capture ds __*
assert _rc==111
restore
''')
    run('existing_file_preserved', BASIC+'''
preserve
keep id
gen sentinel=991
save "existing.dta"
restore
capture noisily dea x = y, dmu(id) generate(eff) reference saveref saving("existing.dta")
assert _rc==602
assert _N==4
confirm variable x y
use "existing.dta", clear
assert sentinel==991
capture confirm variable eff
assert _rc==111
''')
    run('save_failure_restores_sample', BASIC+'''
gen unrelated=991
capture noisily dea x = y in 2/4, dmu(id) generate(eff) reference saveref saving("nonexistent/subdir/out.dta")
assert _rc==603
assert _N==4
assert unrelated==991
confirm variable id x y
''')
    for seed in [1, 22, 301]:
        run('stress_'+str(seed), f'''set seed {seed}
set obs 120
gen id=_n
gen double x1=1+10*runiform()
gen double x2=1+10*runiform()
gen double y1=1+5*runiform()
gen double y2=1+5*runiform()
dea x1 x2 = y1 y2, dmu(id) stage(1) engine(auto) generate(eff) reference saveref
'''+checks())
    run('auto_preemptive_ipm', BASIC+'''
replace x=x*1e12
dea x = y, dmu(id) engine(auto) generate(eff) reference saveref
local refs "`r(refvars)'"
matrix L=r(lambda)
assert r(auto_preemptive_ipm)==1
assert "`r(engine)'"=="ipm"
assert !missing(eff)
mata: assert(max(abs(st_matrix("L")-st_data(.,tokens(st_local("refs")))))<1e-12)
''')
    if baseline:
        run('old_defect_missing_reference', BASIC+'''
dea x = y, dmu(id) generate(eff) reference saveref
matrix L=r(lambda)
assert !missing(eff)
capture ds _ref_*
assert _rc==111
''', baseline)
        run('old_defect_stale_returns', BASIC+'''
dea x = y, dmu(id) engine(auto) generate(eff) reference
matrix L=r(lambda)
assert "`r(engine)'"=="ipm"
assert r(n_failed)==4
assert missing(r(mean_eff))
assert !missing(eff)
mata: assert(all(st_matrix("L"):>=.))
''', out/'baseline_forced_retry_ado')
        run('old_defect_double_save', BASIC+'''
capture noisily dea x = y, dmu(id) engine(auto) generate(eff) reference saving("final.dta")
assert _rc==602
assert !missing(eff)
use "final.dta", clear
assert missing(eff)
''', out/'baseline_forced_retry_ado')
        run('old_defect_auto_id_overwrite', BASIC+'''
gen str20 _dmu="original"+string(_n)
dea x = y, generate(eff)
assert _dmu=="dmu"+string(_n)
''', baseline)
    failed = [r['test'] for r in results if not r['passed']]
    print(json.dumps(dict(total=len(results), passed=len(results)-len(failed), failed=failed)))
    raise SystemExit(bool(failed))


if __name__ == '__main__':
    main()
