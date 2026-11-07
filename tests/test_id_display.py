#!/usr/bin/env python3
"""DMU labels remain data in core, FDH, interval, and reference displays.

Tests special-character labels against numeric-ID controls, including exact
numerical returns, source data, saved files, caller objects, and literal output.
Uses the Python standard library and licensed Stata 16+.
"""
from _support import package_files
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import re
import shutil
import subprocess
import time

from test_reference_state import DATA, guard, reference_lines

REPO = Path(__file__).resolve().parents[1]
PROFILES = {
    'plain': ('strL', '"unit "+string(obsno)'),
    'quote': ('strL', '"U"+char(34)+"Q"+string(obsno)'),
    'compound_quote': ('strL', '"U"+char(96)+char(34)+"Q"+char(34)+char(39)+string(obsno)'),
    'macro': ('strL', 'char(36)+"F27_GLOBAL "+char(96)+"dmu_text"+char(39)+string(obsno)'),
    'smcl': ('strL', '"{bf:raw}{hline 2}"+string(obsno)'),
    'unicode': ('strL', '"기관 "+char(34)+"가"+char(34)+string(obsno)'),
    'long': ('strL', '"START"'),
    'blank_duplicate': ('strL', 'cond(mod(obsno,3)==0,"",cond(mod(obsno,2)," ","same"))'),
    'numeric': ('double', 'obsno'),
    'fraction': ('double', '-obsno/7'),
    'missing': ('double', 'cond(obsno==2,.a,cond(obsno==4,.z,obsno))'),
    'large': ('double', '1e15+obsno'),
}
SCALARS = ['N', 'N_eff', 'mean_eff', 'min_eff', 'max_eff', 'n_failed',
           'n_radial_failed', 'n_stage2_failed', 'n_crs', 'n_irs', 'n_drs',
           'mean_eff_opt', 'n_failed_primary', 'n_failed_opt']
MATRICES = ['efficiency', 'lambda', 'slack_in', 'slack_out', 'vrsfrontier', 'vrsobs']
QUALIFIERS = [('', '1'), ('if obsno!=2', 'obsno!=2'),
              ('in 2/8', 'inrange(_n,2,8)'),
              ('if obsno!=2 in 1/7', 'obsno!=2 & inrange(_n,1,7)')]


def fixture(profile):
    kind, expression = PROFILES[profile]
    s = DATA + f'gen {kind} id={expression}\ngen long control_id=_n\n'
    if profile == 'long':
        s += 'forvalues j=1/220 {\n quietly replace id=id+"aaaaaaaaaa"\n}\nreplace id=id+char(34)+"TAIL"+string(obsno)\n'
    if kind == 'double':
        s += 'format id %9.2f\ngen strL display_label=strtrim(string(id,"%18.0g"))\n'
    else:
        s += 'gen strL display_label=id\n'
    s += 'global F27_GLOBAL "SHOULD_NEVER_REPLACE_LABEL"\n'
    return s


def command(options, label, prefix, form, qualifier, extra=''):
    head = f'dea x1 x2=y1 y2 {qualifier},' if form == 'A' else f'dea y1 y2 {qualifier}, inputs(x1 x2)'
    return f'{head} {options} dmu({label}) gen({prefix}) {extra}\n'


def paired(profile, options, index, references=False):
    q, selected = QUALIFIERS[index % len(QUALIFIERS)]
    form = 'A' if index % 2 else 'B'
    s = fixture(profile) + f'gen byte selected={selected}\n'
    s += command(options, 'control_id', 'control', form, q, 'reference' if references else '')
    s += 'matrix CStats=(' + ','.join('r('+v+')' for v in SCALARS) + ')\n'
    for j, name in enumerate(MATRICES):
        s += f'capture matrix C{j}=r({name})\nlocal has{j}=(_rc==0)\n'
    s += 'quietly ds\nlocal original_vars "`r(varlist)\'"\ntempfile original\nsave `original\'\n'
    setup, check = guard('mixed')
    s += setup
    if index % 3 == 0:
        s += 'gen str12 res="replace me"\n'
    s += 'di "F27_DISPLAY_BEGIN"\n'
    extra = 'replacevars saving("saved.dta")' + (' reference saveref' if references else '')
    s += command(options, 'id', 'res', form, q, extra)
    s += check
    s += 'matrix RStats=(' + ','.join('r('+v+')' for v in SCALARS) + ')\n'
    s += 'mata: assert(all(st_matrix("CStats"):==st_matrix("RStats")))\n'
    for j, name in enumerate(MATRICES):
        s += f'''capture matrix R{j}=r({name})
assert (_rc==0)==`has{j}'
if `has{j}' {{
    mata: assert(all(st_matrix("C{j}"):==st_matrix("R{j}")))
}}
'''
    s += 'local refs "`r(refvars)\'"\nlocal auxiliaries "`r(idea_vars)\'"\n'
    s += 'assert res==control\nassert missing(res) if !selected\n'
    s += 'capture confirm variable control_opt\nif !_rc assert res_opt==control_opt\n'
    s += 'capture confirm variable control_mid\nif !_rc assert res_mid==control_mid\n'
    s += 'cf `original_vars\' using `original\',all\n'
    s += 'assert "`:char note[audit]\'"=="unchanged"\n'
    s += 'assert "$F27_GLOBAL"=="SHOULD_NEVER_REPLACE_LABEL"\n'
    s += 'di "F27_DISPLAY_END"\n'
    s += '''preserve
keep if selected
export delimited display_label using "labels.csv",replace
keep id x1 x2 y1 y2 res `refs' `auxiliaries'
tempfile expected
save `expected'
use "saved.dta",clear
cf _all using `expected',all
restore
capture ds __*
assert _rc==111
'''
    return s


def literal_output(folder, profile):
    lines = (folder/'test.log').read_text(errors='replace').splitlines()
    start, end = lines.index('F27_DISPLAY_BEGIN'), lines.index('F27_DISPLAY_END')
    output = '\n'.join(lines[start+1:end])
    labels = [row['display_label'] for row in csv.DictReader((folder/'labels.csv').open())]
    if profile == 'long':
        # The 12-character score column may clip a label and Stata wraps long
        # reference lines. Source/saved strL values are checked in full by cf.
        assert 'START' in output
    else:
        for label in labels:
            if label.strip():
                assert label in output, (profile, label, output[-1200:])
    assert 'SHOULD_NEVER_REPLACE_LABEL' not in '\n'.join(
        line for line in output.splitlines() if not line.startswith(('.', '>'))
    ), 'A global macro was expanded inside a label'
    return {'labels_checked': len(labels)}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata', required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--baseline-ado', type=Path)
    ap.add_argument('--only')
    args = ap.parse_args()
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    pkg = out/'current_ado'; pkg.mkdir()
    hashes = {}
    for name in package_files(REPO):
        raw = (REPO/name).read_bytes(); (pkg/Path(name).name).write_bytes(raw)
        hashes[name] = hashlib.sha256(raw).hexdigest()
    # Private failure injection covers FDH's otherwise finite-only display.
    fault = out/'forced_fdh_failure_ado'; shutil.copytree(pkg, fault)
    path = fault/'dea.ado'; source = path.read_text()
    anchor = '    quietly mata: _dea_fdh("`inputs\'", "`outputs\'", "`efficiency\'", "`touse\'", "`orientation\'")\n'
    assert source.count(anchor) == 1
    path.write_text(source.replace(anchor, anchor + '    quietly replace `efficiency\'=.a if `touse\' & _n==3 // TEST ONLY\n'))
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,baseline=str(args.baseline_ado),stata=args.stata),indent=2)+'\n')
    results = []

    def run(name, body, source=pkg, profile=None):
        if args.only and not re.search(args.only, name):
            return
        d = out/name; d.mkdir()
        head = 'version 16\nclear all\nset more off\nset linesize 255\n'
        for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:
            (d/key).mkdir(); head += f'sysdir set {key} "{d/key}"\n'
        head += f'adopath ++ "{source.resolve()}"\n'
        (d/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n')
        row = dict(test=name,passed=False); start = time.monotonic()
        try:
            proc = subprocess.run([args.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=180)
            (d/'process.json').write_text(json.dumps(dict(returncode=proc.returncode,stderr=proc.stderr.decode(errors='replace')),indent=2)+'\n')
            log = (d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else ''
            assert proc.returncode == 0 and 'TEST_PASS' in log.splitlines(), log[-2600:]
            if profile:
                row.update(literal_output(d,profile))
            row['passed'] = True
        except Exception as exc:
            row['error'] = repr(exc)
        row['seconds'] = round(time.monotonic()-start,3); results.append(row)
        (out/'results.json').write_text(json.dumps(results,ensure_ascii=False,indent=2)+'\n')
        if not row['passed']:
            print(json.dumps(dict(test=name,passed=False,error=row['error'][-400:]),ensure_ascii=False),flush=True)

    for i, (profile, engine, rts, ori) in enumerate(itertools.product(PROFILES, ['simplex','ipm','auto'], ['crs','vrs'], ['input','output'])):
        options = f'engine({engine}) rts({rts}) orientation({ori}) stage({1+i//2%2})'
        run(f'core_{profile}_{engine}_{rts}_{ori}', paired(profile, options, i, True), profile=profile)
    for i, (profile, ori, kind) in enumerate(itertools.product(PROFILES,['input','output'],['fdh','interval_point','interval_width','legacy'])):
        options = {'fdh':'model(fdh)', 'interval_point':'xlo(x1 x2) xhi(x1 x2)',
                   'interval_width':'xlo(lo x2) xhi(hi x2)', 'legacy':'xord(x1)'}[kind]
        s = paired(profile, options+f' orientation({ori})', i)
        if kind == 'interval_width':
            s = s.replace('gen long control_id=_n\n', 'gen long control_id=_n\ngen double lo=x1*.8\ngen double hi=x1*1.2\n')
        # Legacy IDEA has an aggregate-only display; still check complete data,
        # result, and saved-file preservation for every label profile.
        run(f'{kind}_{profile}_{ori}',s,profile=profile if kind!='legacy' else None)
    for i, profile in enumerate(['quote','compound_quote','macro','smcl','unicode','numeric']):
        run('fdh_failure_'+profile,paired(profile,'model(fdh)',i),fault,profile)
        s = fixture(profile)+'gen byte selected=1\ngen double score=1\nreplace score=.a in 2\n'
        for j in range(1,9):
            s += f'gen double w{j}=0\nreplace w{j}=.125 in 1\n'
        s += 'replace w2=2e-7 in 2\nreplace w3=.a in 2\nreplace w4=-1 in 2\n'
        s += f'quietly run "{pkg/"dea.ado"}"\n'
        s += 'di "F27_DISPLAY_BEGIN"\n'
        s += '_dea_display, model(te) rts(crs) orientation(input) gen(score) touse(selected) tolerance(1e-7) stage(2) dmu(id) reference lambda_vars("w1 w2 w3 w4 w5 w6 w7 w8")\n'
        s += 'assert r(n_failed)==1\ndi "F27_DISPLAY_END"\nexport delimited display_label using "labels.csv",replace\n'
        run('reference_failure_'+profile,s,profile=profile)
    # Exhaust the printable ASCII range in saved peer identifiers, including
    # backticks retained by Stata 16's strtoname(). Validate all name collisions
    # and metadata without an LP or command-string interpolation of ID values.
    s = f'quietly run "{pkg/"dea.ado"}"\n'
    s += '''set obs 95
gen str1 id=char(_n+31)
gen byte selected=1
mata: _dea_saveref_plan("id","selected")
foreach v of local _ref_names {
    confirm name `v'
    confirm new variable `v'
}
matrix L=I(95)
mata: _dea_saveref_store("id","selected",st_local("_ref_names"),"L")
mata:
names=tokens(st_local("_ref_names"))
assert(cols(names)==95)
assert(all(st_data(.,names):==I(95)))
for (j=1;j<=95;j++) {
    assert(st_global(names[j]+"[dea_peer_id]")==st_sdata(j,"id"))
    assert(st_global(names[j]+"[dea_peer_index]")==strofreal(j))
}
end
'''
    run('saveref_printable_ascii',s)
    if args.baseline_ado:
        for i, (engine,rts,ori) in enumerate(itertools.product(['simplex','ipm','auto'],['crs','vrs'],['input','output'])):
            profile = ['plain','numeric','fraction','missing','large'][i%5]
            name = f'history_{engine}_{rts}_{ori}'
            cmd = command(f'engine({engine}) rts({rts}) orientation({ori})','id','res','A','','reference saveref')
            run(name+'_before',fixture(profile)+cmd+'keep res _ref_*\nsave "expected.dta"\n',args.baseline_ado)
            run(name+'_after',fixture(profile)+cmd+f'cf res _ref_* using "{out/(name+"_before")/"expected.dta"}",all\n')
            if not args.only or re.search(args.only,name):
                before = reference_lines((out/(name+'_before')/'test.log').read_text(errors='replace'))
                after = reference_lines((out/(name+'_after')/'test.log').read_text(errors='replace'))
                assert before == after, (name,before,after)
    summary = dict(total=len(results),passed=sum(r['passed'] for r in results),failed=sum(not r['passed'] for r in results))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary),flush=True)
    raise SystemExit(summary['failed'] > 0)


if __name__ == '__main__':
    main()
