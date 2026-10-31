#!/usr/bin/env python3
"""Execute literal help examples and public DEA syntax combinations in Stata."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time

REPO = Path(__file__).resolve().parents[1]
DATA = '''input id x1 x2 y1 y2 group
1 1 2 3 2 1
2 2 1 4 3 0
3 3 4 5 4 1
4 4 2 6 2 1
5 2 3 4 5 0
6 3 2 5 3 1
end
gen str4 tag=cond(group,"A B","C")
gen unrelated=991
'''


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stata', required=True)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--only', help='Run case names matching this regex')
    args = p.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    package = out/'ado'
    package.mkdir()
    tracked = subprocess.check_output(['git', '-C', str(REPO), 'ls-files', 'ado/'], text=True).splitlines()
    hashes = {}
    for name in tracked:
        data = (REPO/name).read_bytes()
        (package/Path(name).name).write_bytes(data)
        hashes[name] = hashlib.sha256(data).hexdigest()
    cases = {}
    example_map = {}
    for name in tracked:
        if not name.endswith('.sthlp'):
            continue
        code = []
        mapping = []
        active = False
        for lineno, line in enumerate((REPO/name).read_text().splitlines(), 1):
            if 'DEAS_HELP_EXAMPLES_BEGIN' in line:
                active = True
            elif 'DEAS_HELP_EXAMPLES_END' in line:
                active = False
            elif active:
                match = re.fullmatch(r'\{phang2\}\{cmd:(.*)\}\{p_end\}', line)
                if match:
                    command = re.sub(r'^\.\s*', '', match[1])
                    command = command.replace('{c -(}', '{').replace('{c )-}', '}')
                    command = re.sub(r'\{space (\d+)\}', lambda m: ' '*int(m[1]), command)
                    code.append(command)
                    mapping.append({'help_line': lineno, 'command': command})
        assert code, f'No marked executable examples in {name}'
        key = 'help_'+Path(name).stem
        example_map[key] = mapping
        cases[key] = '''set obs 2
gen original_sentinel=991
datasignature
local signature_before "`r(datasignature)'"
'''+ '\n'.join(code)+'''
datasignature
assert "`signature_before'"=="`r(datasignature)'"
'''

    # The fuzzy help also has an executable alpha=1 check in Remarks.
    fuzzy_text=(REPO/'ado/dea_fuzzy.sthlp').read_text()
    remarks=fuzzy_text.split('{it:Verifying results at α = 1.}',1)[1].split('{marker examples}',1)[0]
    remark_commands=[re.sub(r'^\.\s*','',x) for x in re.findall(r'^\{phang2\}\{cmd:(.*)\}\{p_end\}$',remarks,re.M)]
    prelude=[]
    for row in example_map['help_dea_fuzzy']:
        if row['command'].startswith('dea_fuzzy'):
            break
        prelude.append(row['command'])
    cases['remarks_fuzzy_alpha1']='\n'.join(prelude)+'''\ngen str20 _dmu="existing user ID"\n'''+ '\n'.join(remark_commands)+'''\nrestore\n'''

    # Equivalent commands in both public syntax forms, with sample restrictions
    # containing equality signs that must not be mistaken for the input separator.
    qualifiers = {'all':'', 'if':'if group==1', 'in':'in 2/6',
                  'if_in':'if tag=="A B" in 2/6'}
    for form in ['equal', 'inputs']:
        for engine in ['simplex', 'ipm', 'auto']:
            for rts in ['crs', 'vrs']:
                for qname, qualifier in qualifiers.items():
                    for mode in ['score', 'reference', 'saveref']:
                        key=f'syntax_{form}_{engine}_{rts}_{qname}_{mode}'
                        extra={'score':'', 'reference':'reference', 'saveref':'reference saveref'}[mode]
                        command=(f'dea x1 x2 = y1 y2 {qualifier},' if form=='equal'
                                 else f'dea y1 y2 {qualifier}, inputs(x1 x2)')
                        opts=f'dmu(id) rts({rts}) engine({engine}) orientation(output) stage(2)'
                        cases[key]=DATA+f'''
dea x1 x2 = y1 y2 {qualifier}, {opts} generate(control) reference
matrix C=r(lambda)
local N0=r(N)
{command} {opts} generate(eff) {extra} saving("saved results.dta")
local refs "`r(refvars)'"
assert r(N)==`N0'
assert missing(eff)==missing(control)
assert abs(eff-control)<1e-10 if !missing(eff)
'''+('''matrix L=r(lambda)
mata: assert(max(abs(st_matrix("C")-st_matrix("L")))<1e-10)
''' if mode!='score' else '''capture matrix list r(lambda)
assert _rc==111
''')+'''
assert unrelated==991
preserve
use "saved results.dta", clear
assert _N==`N0'
assert !missing(eff)
capture ds __*
assert _rc==111
capture confirm variable unrelated
assert _rc==111
'''+('''mata: assert(max(abs(st_matrix("L")-st_data(.,tokens(st_local("refs")))))<1e-12)
''' if mode=='saveref' else '''capture ds _ref_*
assert _rc==111
''')+'''restore
'''

    aliases = [
        ('documented_minimum', 'mod(te) r(crs) or(input) met(simplex) gen(eff) tol(1e-7)', 'model(te) rts(crs) orientation(input) engine(simplex)'),
        ('model_minimum', 'mod(te) gen(eff)', 'model(te)'),
        ('rts_minimum', 'r(crs) gen(eff)', 'rts(crs)'),
        ('orientation_minimum', 'or(input) gen(eff)', 'orientation(input)'),
        ('method_minimum', 'met(simplex) gen(eff)', 'engine(simplex)'),
        ('ccr', 'rts(ccr) ort(i) gen(eff)', 'rts(crs) orientation(input)'),
        ('bcc', 'rts(bcc) ort(o) gen(eff)', 'rts(vrs) orientation(output)'),
        ('drs', 'rts(drs) gen(eff)', 'rts(nirs)'),
        ('irs', 'rts(irs) gen(eff)', 'rts(ndrs)'),
        ('orientation_inp', 'ori(inp) gen(eff)', 'orientation(input)'),
        ('orientation_outp', 'ori(outp) gen(eff)', 'orientation(output)'),
        ('engine_minimum', 'eng(ipm) gen(eff)', 'engine(ipm)'),
        ('engine_precedence', 'engine(simplex) method(ipm) gen(eff)', 'engine(simplex)'),
        ('reference_save_minimum', 'ref saver gen(eff) sav("short path.dta")', 'reference'),
        ('uppercase_values', 'model(TE) rts(BCC) orientation(OUT) engine(SIMPLEX) gen(eff)', 'model(te) rts(vrs) orientation(output) engine(simplex)'),
        ('debug_trace', 'debug trace gen(eff)', ''),
    ]
    for form in ['equal', 'inputs']:
        prefix='dea x1 x2 = y1 y2,' if form=='equal' else 'dea y1 y2, i(x1 x2)'
        for name, opts, canonical in aliases:
            cases[f'alias_{form}_{name}']=DATA+f'''
dea x1 x2 = y1 y2, dmu(id) {canonical} generate(control)
{prefix} dmu(id) {opts}
assert abs(eff-control)<1e-10
'''

    setup_prices='''gen w1=.2
gen w2=.3
gen p1=10
gen p2=12
gen bad=1+id/10
gen xlo1=x1*.9
gen xlo2=x2*.9
gen xhi1=x1*1.1
gen xhi2=x2*1.1
gen gx1=1
gen gx2=1
gen gy1=1
gen gy2=1
'''
    models={'te':'reference', 'multiplier':'saveweights', 'cost':'w(w1 w2)',
            'revenue':'p(p1 p2)', 'profit':'w(w1 w2) p(p1 p2)',
            'ddf':'gx(gx1 gx2) gy(gy1 gy2)', 'sbm':'', 'additive':'',
            'supereff':'reference saveref', 'supersbm':'', 'allocative':'w(w1 w2)',
            'fdh':'', 'imprecise':'xlo(xlo1 xlo2) xhi(xhi1 xhi2)',
            'undesirable':'bad(bad) direction(xyb)', 'congestion':'', 'sbmundesirable':'bad(bad)'}
    for model, opts in models.items():
        # Both forms are called in the same session to catch persistent-variable interactions.
        canonical_opts=opts.replace(' saveref','')
        cases['model_forms_'+model]=DATA+setup_prices+f'''
dea x1 x2 = y1 y2, dmu(id) model({model}) {canonical_opts} generate(control)
dea y1 y2, inputs(x1 x2) dmu(id) model({model}) {opts} generate(eff) saving("model.dta")
assert missing(eff)==missing(control)
assert abs(eff-control)<1e-8 if !missing(eff)
use "model.dta", clear
confirm variable id x1 x2 y1 y2 eff
capture ds __*
assert _rc==111
'''
    cases['ddfbeta_alias']=DATA+'''
dea y1 y2, inputs(x1 x2) model(ddf) dmu(id) ddfbeta(beta)
confirm variable beta
'''
    cases['wildcard_compact_equal']=DATA+'''
dea x1 x2 = y1 y2, dmu(id) generate(control)
quietly dea x*=y*, gen(eff) dmu(id)
assert abs(eff-control)<1e-10
'''
    cases['range_varlists']=DATA+'''
dea x1 x2 = y1 y2 if inlist(group,0,1) in 2/5, dmu(id) generate(control)
dea y1-y2 if inlist(group,0,1) in 2/5, generate(eff) inputs(x1-x2) dmu(id)
assert missing(eff)==missing(control)
assert abs(eff-control)<1e-10 if !missing(eff)
'''
    cases['using_subset_reference']=DATA+'''
dea x1 x2 = y1 y2 if group==1 in 1/5 using "subset weights.dta", dmu(id) gen(eff) ref saver
local refs "`r(refvars)'"
local n=r(N)
matrix L=r(lambda)
use "subset weights.dta", clear
assert _N==`n'
mata: assert(max(abs(st_matrix("L")-st_data(.,tokens(st_local("refs")))))<1e-12)
'''
    cases['savew_minimum']=DATA+'''
dea y1 y2, i(x1 x2) mod(multiplier) dmu(id) savew gen(eff) sav("weights.dta")
use "weights.dta", clear
confirm variable eff _dea_v1 _dea_v2 _dea_u1 _dea_u2 _dea_u0
'''
    cases['using_legacy']=DATA+'''
dea x1 x2 = y1 y2 using "legacy results.dta", dmu(id) generate(eff)
use "legacy results.dta", clear
confirm variable eff
'''
    cases['bootstrap_saveref_saving']=DATA+'''
set seed 12345
dea y1 y2, inputs(x1 x2) dmu(id) rts(vrs) engine(auto) generate(eff) reference saveref bootstrap reps1(20) saving("bootstrap.dta")
local refs "`r(refvars)'"
matrix L=r(lambda)
assert !missing(eff)
use "bootstrap.dta", clear
mata: assert(max(abs(st_matrix("L")-st_data(.,tokens(st_local("refs")))))<1e-12)
'''
    for form in ['equal','inputs']:
        prefix='dea x1 x2 = y1 y2,' if form=='equal' else 'dea y1 y2, inputs(x1 x2)'
        cases['bootstrap_minimum_'+form]=DATA+f'''
set seed 12345
{prefix} dmu(id) generate(eff) boo reps1(20)
assert !missing(eff)
'''
    cases['using_conflict']=DATA+'''
datasignature
local before "`r(datasignature)'"
capture noisily dea x1 x2 = y1 y2 using "legacy.dta", dmu(id) saving("modern.dta")
assert _rc==198
datasignature
assert "`before'"=="`r(datasignature)'"
capture confirm file "legacy.dta"
assert _rc==601
capture confirm file "modern.dta"
assert _rc==601
'''
    invalid=[('conflicting_orientation','ort(i) orientation(input)',198),
             ('missing_reference','saveref',198),
             ('wrong_saveweights_model','saveweights',198),
             ('wrong_engine','engine(unknown)',198),
             ('wrong_rts','rts(unknown)',198),
             ('wrong_stage','stage(3)',198),
             ('wrong_model','model(unknown)',198),
             ('cost_no_prices','model(cost)',198),
             ('bad_no_bad','model(undesirable)',198),
             ('ddfbeta_wrong_model','ddfbeta(beta)',198),
             ('ddfbeta_and_generate','model(ddf) ddfbeta(beta) generate(eff)',198)]
    for form in ['equal','inputs']:
        prefix='dea x1 x2 = y1 y2,' if form=='equal' else 'dea y1 y2, inputs(x1 x2)'
        for name, opts, expected in invalid:
            cases[f'invalid_{form}_{name}']=DATA+'''
datasignature
local before "`r(datasignature)'"
'''+f'''capture noisily {prefix} dmu(id) {opts}
assert _rc=={expected}
datasignature
assert "`before'"=="`r(datasignature)'"
'''

    (out/'help_line_map.json').write_text(json.dumps(example_map, indent=2))
    (out/'manifest.json').write_text(json.dumps({'source_sha256':hashes, 'only':args.only,
        'help_blocks':len(example_map), 'cases_available':len(cases)}, indent=2))
    results=[]
    for name, body in cases.items():
        if args.only and not re.search(args.only, name):
            continue
        work=out/name
        work.mkdir()
        header='version 16\nclear all\nset more off\nset linesize 120\n'
        for folder in ['plus','personal','site','oldplace']:
            (work/folder).mkdir()
            header+=f'sysdir set {folder.upper()} "{work/folder}"\n'
        header+=f'adopath ++ "{package}"\n'
        (work/'test.do').write_text(header+body+'\ndisplay "TEST_PASS"\n')
        start=time.monotonic()
        try:
            proc=subprocess.run([args.stata,'-b','do','test.do'],cwd=work,timeout=180,capture_output=True,text=True)
            log=(work/'test.log').read_text(errors='replace') if (work/'test.log').exists() else proc.stderr
            result={'test':name,'passed':'TEST_PASS' in log.splitlines(),'seconds':round(time.monotonic()-start,3)}
            if not result['passed']:
                result['tail']=log[-2000:]
        except subprocess.TimeoutExpired:
            result={'test':name,'passed':False,'error':'timeout'}
        results.append(result)
        (out/'results.json').write_text(json.dumps(results,indent=2))
        print(json.dumps({k:v for k,v in result.items() if k!='tail'}),flush=True)
    failed=[r['test'] for r in results if not r['passed']]
    print(json.dumps({'total':len(results),'passed':len(results)-len(failed),'failed':failed}))
    raise SystemExit(bool(failed))


if __name__=='__main__':
    main()
