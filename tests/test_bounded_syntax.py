#!/usr/bin/env python3
"""Restore bounded-string syntax and check endpoint parsing and safe rejection.

Requires Python's standard library and licensed Stata 16+. Optional pre-F29
private-worker controls preserve unaffected historical scores. Ordinal-output
controls use the corrected current worker after F30; independent transformation
and permutation oracles are in test_ordinal_order.py. This is not a methodology
certification. All packages, files and instrumentation are isolated copies.
"""
from _support import package_files
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import re
import subprocess
import time

from test_fdh_idea_returns import checks
REPO=Path(__file__).resolve().parents[1]
DATA='''set obs 6
gen long id=_n
gen double x=2+mod(id*7,9)
gen double y=2+mod(id*3,7)
gen double qx=1+mod(id*2,5)
gen double qy=1+mod(id*3,5)
gen str24 xb="["+strtrim(string(x-.5,"%21.17g"))+","+strtrim(string(x+.5,"%21.17g"))+"]"
gen str24 yb="["+strtrim(string(y-.5,"%21.17g"))+","+strtrim(string(y+.5,"%21.17g"))+"]"
gen str20 tag=cond(mod(id,3),"keep = this, yes","omit")
'''
PROFILES={'input':('xb qx','y qy','xb'), 'output':('x qx','yb qy','yb'), 'both':('xb qx','yb qy','xb yb')}
INVALID=['',' ','1','1,2','[1]','[1,2,3]','[[1,2]]','[1,2]]','junk[1,2]','[1,2]junk',
         '(1,2)','[.,2]','[.a,2]','[1,.z]','[NaN,2]','[1,Inf]','[-Inf,2]',
         '[1e309,2e309]','[2,1]','[+,2]','[1 2,3]','[1e,2]','[1,,2]','[1,]','[,2]']
FORMATS={'standard':'"[1,2]"','decimals':'"[1.000,2.]"','signed':'"[+1,+2]"',
         'scientific':'"[.10e1,20E-1]"','leading_zero':'"[01,002]"',
         'inner_space':'"[  1 ,  2 ]"','outer_space':'"   [1,2]   "',
         'tabs':'char(9)+"[1,"+char(9)+"2]"+char(9)',
         'newlines':'char(10)+"[1,2]"+char(13)',
         'unicode_space':'uchar(160)+"["+uchar(160)+"1,2]"+uchar(160)'}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--baseline-ado',type=Path);ap.add_argument('--only')
    args=ap.parse_args()
    if args.baseline_ado and not (args.baseline_ado/'dea.ado').is_file():ap.error('Invalid pre-F29 ado snapshot')
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'current_ado';pkg.mkdir();hashes={}
    for name in package_files(REPO):
        data=(REPO/name).read_bytes();(pkg/Path(name).name).write_bytes(data);hashes[name]=hashlib.sha256(data).hexdigest()
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,stata=args.stata,baseline=str(args.baseline_ado)),indent=2))
    results=[]
    def chosen(name):return not args.only or re.search(args.only,name)
    def run(name,script,source=pkg,force=False):
        if not force and not chosen(name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 120\n'
        for folder in ['PLUS','PERSONAL','SITE','OLDPLACE']:
            (d/folder).mkdir();head+=f'sysdir set {folder} "{d/folder}"\n'
        head+=f'adopath ++ "{source.resolve()}"\n';(d/'test.do').write_text(head+script+'\ndi "TEST_PASS"\n')
        start=time.monotonic();row=dict(test=name,passed=False)
        try:
            p=subprocess.run([args.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=240)
            log=(d/'test.log').read_text(errors='replace')
            if p.returncode or 'TEST_PASS' not in log.splitlines():raise RuntimeError(log[-2300:])
            if (d/'target_before.dta').exists():
                assert (d/'target.dta').read_bytes()==(d/'target_before.dta').read_bytes(),'Rejected call changed the existing file'
            row['passed']=True
        except Exception as e:row['error']=repr(e)
        row['seconds']=round(time.monotonic()-start,3);results.append(row)
        (out/'results.json').write_text(json.dumps(results,indent=2));print(json.dumps(row),flush=True)

    # Public syntax must match the worker. F30 intentionally changes erroneous
    # ordinal-output scores, so only no-output-ordinal cases use pre-F29 scores.
    # Current ordinal-output controls check routing; the separate F30 suite
    # checks transforms independently and tests permutation invariance.
    for i,((kind,(inputs,outputs,bounded)),ordkind,rts,ori,engine) in enumerate(itertools.product(PROFILES.items(),['none','input','output','both'],['crs','vrs','nirs','ndrs'],['input','output'],['simplex','ipm','auto'])):
        name=f'grid_{kind}_{ordkind}_{rts}_{ori}_{engine}'
        if not chosen(name):continue
        ordinal={'none':'','input':'qx','output':'qy','both':'qx qy'}[ordkind]
        setup=DATA
        if i%5==0:setup+='gsort -id\n'
        if i%7==0:setup+='replace qx=.a in 5\n'
        condition=['1','tag=="keep = this, yes"','inrange(_n,2,5)','mod(id,3)!=0 & inrange(_n,2,5)'][i%4]
        qualifier=['','if tag=="keep = this, yes"','in 2/5','if mod(id,3)!=0 in 2/5'][i%4]
        setup+=f'gen byte selected=({condition}) & !missing(qx,qy)\n'
        internal_rts={'nirs':'drs','ndrs':'irs'}.get(rts,rts)
        controlpkg=(args.baseline_ado or pkg) if ordkind in ['none','input'] else pkg
        # engine(auto) resolves to simplex; this fixture has no retry failures.
        worker_engine='simplex' if engine=='auto' else engine
        control=setup+f'quietly run "{controlpkg.resolve()}/dea.ado"\n'
        control+=f'_dea_imprecise_legacy, inputs({inputs}) outputs({outputs}) bounded({bounded}) ordinal({ordinal}) dmu(id) rts({internal_rts}) orientation({ori}) stage({1+i%2}) engine({worker_engine}) touse(selected) generate(res) tolerance(1e-7)\n'
        control+='assert !missing(res) if selected\nkeep res res_opt res_mid\nsave "expected.dta"\n'
        run('control_'+name,control,controlpkg,force=True)
        path=out/('control_'+name)/'expected.dta'
        setup+='quietly ds\nlocal original_vars "`r(varlist)\'"\ntempfile original\nsave `original\'\n'
        setup+='foreach v in CRS_TE VRS_TE SCALE RTS_NUM RTS {\n gen str10 `v\'="untouched"\n}\n'
        replace=i%3==0
        if replace:setup+='foreach v in res res_opt res_mid {\n gen str10 `v\'="old"\n}\n'
        prefix=f'dea {inputs}={outputs} {qualifier},' if i%2 else f'dea {outputs} {qualifier}, inputs({inputs})'
        model=['','model(te)','model(imprecise)'][i%3];requested='imprecise' if i%3==2 else 'te'
        setup+=f'local savevars id {inputs} {outputs}\n'
        setup+=f'{prefix} bounded({bounded}) ordinal({ordinal}) {model} rts({rts}) orientation({ori}) engine({engine}) stage({1+i%2}) dmu(id) gen(res) saving("saved.dta")'+(' replacevars' if replace else '')+'\n'
        setup+=checks('imprecise',requested,'legacy','res res_opt res_mid',internal_rts,saving=True)
        setup+=f'cf res res_opt res_mid using "{path}", all\ncf `original_vars\' using `original\', all\n'
        setup+='foreach v in CRS_TE VRS_TE SCALE RTS_NUM RTS {\n assert `v\'=="untouched"\n}\n'
        run(name,setup)

    def reject(setup,cmd,rc):
        s=setup+'''gen str12 res="old score"
gen str12 res_opt="old opt"
gen str12 res_mid="old mid"
label variable res "original label"
char res[review] "original characteristic"
char _dta[review] "dataset characteristic"
tempfile original
save `original'
save "target.dta"
copy "target.dta" "target_before.dta"
'''
        s+='capture noisily '+cmd+f'\nassert _rc=={rc}\n'
        s+='''cf _all using `original',all
assert "`:variable label res'"=="original label"
assert "`:char res[review]'"=="original characteristic"
assert "`:char _dta[review]'"=="dataset characteristic"
capture confirm variable _dmu,exact
assert _rc==111
'''
        return s

    for i,(value,(kind,(inputs,outputs,bounded)),form,replace) in enumerate(itertools.product(INVALID,PROFILES.items(),['A','B'],[False,True])):
        target='xb' if kind=='input' else 'yb'
        setup=DATA+f'replace {target}="{value}" in 3\n'
        prefix=f'dea {inputs}={outputs} in 2/6,' if form=='A' else f'dea {outputs} in 2/6, inputs({inputs})'
        cmd=f'{prefix} bounded({bounded}) gen(res) saving("target.dta") replace'+(' replacevars' if replace else '')
        run(f'invalid_{i:03d}',reject(setup,cmd,459))

    for role,value,form in itertools.product(['x','y'],['.','.a','.z'],['A','B']):
        setup=DATA+f'replace {role}={value} in 3\n'
        prefix='dea x qx=y qy,' if form=='A' else 'dea y qy, inputs(x qx)'
        run(f'numeric_missing_{role}_{value[1:] or "dot"}_{form}',reject(setup,f'{prefix} bounded({role}) gen(res) replacevars saving("target.dta") replace',459))

    for (kind,(inputs,outputs,bounded)),fmt,form in itertools.product(PROFILES.items(),FORMATS,['A','B']):
        setup=DATA
        for v in bounded.split():setup+=f'replace {v}="[1,2]"\n'
        setup+=f'dea {inputs}={outputs}, bounded({bounded}) dmu(id) gen(control)\n'
        for v in bounded.split():setup+=f'recast strL {v}\nreplace {v}={FORMATS[fmt]}\n'
        prefix=f'dea {inputs}={outputs},' if form=='A' else f'dea {outputs}, inputs({inputs})'
        setup+=f'{prefix} bounded({bounded}) dmu(id) gen(res)\nassert res==control\n'
        run(f'format_{kind}_{fmt}_{form}',setup)

    # Selection precedes bounded validation. Bad bounds outside the effective
    # sample cannot be accidentally interpreted by the per-DMU worker.
    for kind,(inputs,outputs,bounded) in PROFILES.items():
        target=bounded.split()[-1]
        for how in ['if','in','ordinary_missing','none']:
            setup=DATA+f'replace {target}="not an interval" in 3\n'
            q={'if':'if id!=3','in':'in 4/6','ordinary_missing':'','none':'if 0'}[how]
            if how=='ordinary_missing':setup+='replace qx=.a in 3\n'
            cond={'if':'id!=3','in':'_n>=4','ordinary_missing':'!missing(qx)','none':'0'}[how]
            setup+=f'gen byte selected={cond}\nlocal savevars id {inputs} {outputs}\n'
            cmd=f'dea {inputs}={outputs} {q}, bounded({bounded}) dmu(id) gen(res) saving("saved.dta")'
            if how=='none':setup+='capture noisily '+cmd+'\nassert _rc==2000\n'
            else:setup+=cmd+'\n'+checks('imprecise','te','legacy','res res_opt res_mid',saving=True)
            setup+=f'assert {target}=="not an interval" in 3\n'
            run('selection_'+kind+'_'+how,setup)

    # Point-valued numeric and string bounds have identical meanings.
    for form,engine,ori in itertools.product(['A','B'],['simplex','ipm','auto'],['input','output']):
        setup=DATA+'replace xb="["+string(x)+","+string(x)+"]"\nreplace yb="["+string(y)+","+string(y)+"]"\n'
        setup+=f'dea x qx=y qy, bounded(x y) engine({engine}) orientation({ori}) dmu(id) gen(control)\n'
        prefix='dea xb qx=yb qy,' if form=='A' else 'dea yb qy, inputs(xb qx)'
        setup+=f'{prefix} bounded(xb yb) engine({engine}) orientation({ori}) dmu(id) gen(res)\nassert res==control\n'
        run(f'point_bounds_{form}_{engine}_{ori}',setup)

    # Independent CRS oracle for one input and one output, without ordinal
    # transformations: theta = (own y / own x) / max_j(transformed y_j/x_j).
    # The own DMU uses favorable endpoints; reference DMUs use unfavorable ones.
    for kind,form,engine,ori in itertools.product(PROFILES,['A','B'],['simplex','ipm','auto'],['input','output']):
        ins='xb' if kind in ['input','both'] else 'x'
        outs='yb' if kind in ['output','both'] else 'y'
        bounds=' '.join(v for v in [ins,outs] if v.endswith('b'))
        setup=DATA+f'gen double xlo=x'+('-.5' if ins=='xb' else '')+'\n'
        setup+=f'gen double xhi=x'+('+.5' if ins=='xb' else '')+'\n'
        setup+=f'gen double ylo=y'+('-.5' if outs=='yb' else '')+'\n'
        setup+=f'gen double yhi=y'+('+.5' if outs=='yb' else '')+'\n'
        setup+='gen double expected=.\n'
        setup+='''forvalues i=1/6 {
 gen double productivity=ylo/xhi
 replace productivity=yhi/xlo in `i'
 quietly summarize productivity,meanonly
 replace expected=(yhi/xlo)/r(max) in `i'
 drop productivity
}
'''
        prefix=f'dea {ins}={outs},' if form=='A' else f'dea {outs}, inputs({ins})'
        setup+=f'{prefix} bounded({bounds}) rts(crs) orientation({ori}) engine({engine}) dmu(id) gen(res)\n'
        setup+='assert !missing(res) & abs(res-expected)<1e-7\n'
        run(f'analytic_{kind}_{form}_{engine}_{ori}',setup)

    # Parse endpoint values directly as well as through DEA: signed decimal /
    # exponent grammar, equal endpoints, numeric zero and strL whitespace.
    setup=DATA+f'quietly run "{pkg}/dea.ado"\n'
    setup+='''clear
input str40 bound double expected_lo double expected_hi
"[1,2]" 1 2
"[-2,-1]" -2 -1
"[-1,+1]" -1 1
"[0,0]" 0 0
"[.5,2.]" .5 2
"[-.5e+1,2E1]" -5 20
"[1e-20,2e-20]" 1e-20 2e-20
"[1e20,2e20]" 1e20 2e20
end
_dea_bounded_parse bound, lower(lo) upper(hi)
assert lo==expected_lo & hi==expected_hi
'''
    run('endpoint_values',setup)
    setup=DATA+'replace xb="[1,2]"\ndea xb qx=y qy, bounded(xb) dmu(id) gen(control)\n'
    setup+='gen strL pad=" "\nforvalues k=1/14 {\n replace pad=pad+pad\n}\nrecast strL xb\nreplace xb=pad+"["+pad+"1,"+pad+"2]"+pad\n'
    setup+='dea y qy, inputs(xb qx) bounded(xb) dmu(id) gen(res)\nassert res==control\n'
    run('long_string',setup)

    # String model dimensions remain prohibited unless specifically designated
    # bounded in the legacy route. Check both syntaxes and ordinary model types.
    for dim,model,form in itertools.product(['x','y'],['','model(ddf)','fdh','model(sbm)'],['A','B']):
        setup=DATA+f'tostring {dim},replace\n'
        prefix='dea x qx=y qy,' if form=='A' else 'dea y qy, inputs(x qx)'
        run('unmarked_string_'+dim+'_'+(model or 'te').replace('(','').replace(')','')+'_'+form,
            reject(setup,f'{prefix} {model} gen(res) replacevars saving("target.dta") replace',109))
    for kind,(inputs,outputs,bounded) in PROFILES.items():
        for suffix in ['','_opt','_mid']:
            name='res'+suffix
            setup=DATA+f'gen str10 {name}="old"\n'
            setup+=f'capture noisily dea {inputs}={outputs}, bounded({bounded}) gen(res) dmu(id)\nassert _rc==110\nassert {name}=="old"\n'
            run('protect_result_'+kind+(suffix or '_main'),setup)
            source=bounded.split()[-1]
            setup=DATA+f'rename {source} {name}\ntempfile original\nsave `original\'\n'
            ins=inputs.replace(source,name);outs=outputs.replace(source,name);bounds=bounded.replace(source,name)
            setup+=f'capture noisily dea {ins}={outs}, bounded({bounds}) gen(res) dmu(id) replacevars\nassert _rc==110\ncf _all using `original\',all\n'
            run('protect_bounded_'+kind+(suffix or '_main'),setup)

    for kind,(inputs,outputs,bounded) in PROFILES.items():
        for form in ['A','B']:
            name='original_'+kind+'_'+form
            if not chosen(name):continue
            prefix=f'dea {inputs}={outputs},' if form=='A' else f'dea {outputs}, inputs({inputs})'
            cmd=f'{prefix} bounded({bounded}) ordinal(qx) dmu(id) gen(res)'
            if args.baseline_ado and (kind!='input' or form=='B'):
                run('old_'+name,DATA+'capture noisily '+cmd+'\nassert _rc==109\n',args.baseline_ado,force=True)
            run(name,DATA+cmd+'\nassert r(N)==6\nassert "`r(idea_method)\'"=="legacy"\nassert !missing(res)\n')

    for kind,(inputs,outputs,bounded) in PROFILES.items():
        setup=DATA+'gen byte selected=id!=3\n'
        setup+=f'dea {inputs}={outputs} if selected using "saved.dta", bounded({bounded}) gen(res)\n'
        setup+='assert r(generated_dmu)==1\nassert _dmu=="" if !selected\n'
        setup+=f'local savevars _dmu {inputs} {outputs}\n'+checks('imprecise','te','legacy','res res_opt res_mid',saving=True)
        setup+=f'dea {outputs} if selected, inputs({inputs}) bounded({bounded}) dmu(_dmu) gen(res) replacevars saving("saved.dta") replace\n'
        setup+=checks('imprecise','te','legacy','res res_opt res_mid',saving=True)
        run('autoid_using_replace_'+kind,setup)

    summary=dict(total=len(results),passed=sum(r['passed'] for r in results),failed=sum(not r['passed'] for r in results))
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
    raise SystemExit(summary['failed']>0)


if __name__=='__main__':main()
