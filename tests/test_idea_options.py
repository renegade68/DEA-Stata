#!/usr/bin/env python3
"""Exercise IDEA activation, option compatibility and safe preflight rejection.

Uses Python's standard library and licensed Stata 16+. An optional pre-F10
snapshot reproduces ignored ordinal options and repeated transformations.
No installed ado is modified. Existing ordinal/interval formulations are not
certified here; the test compares routing with the established legacy path.
"""
from _support import package_files
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import re
import shutil
import subprocess
import time

from test_saveref import BASIC, instrument
from test_fdh_idea_returns import checks

REPO=Path(__file__).resolve().parents[1]
DATA='''set obs 8
gen long id=_n
gen double x1=1+mod(id*7,9)
gen double x2=1+mod(id*5,11)
gen double y1=1+mod(id*3,7)
gen double y2=1+mod(id*2,5)
gen double z=id+17
gen double w1=1
gen double w2=2
gen double p1=2
gen double p2=1
gen double bad=1+id/10
gen str20 tag=cond(mod(id,3),"keep = this, yes","omit")
'''
VARIANTS={
 'xonly':('xord(x1)','ordinal(x1)'),
 'yonly':('yord(y1)','ordinal(y1)'),
 'both':('xord(x1) yord(y1)','ordinal(x1 y1)'),
 'multix':('xord(x1 x2)','ordinal(x1 x2)'),
 'multiy':('yord(y1 y2)','ordinal(y1 y2)'),
 'mixed_alias':('ordinal(x2 y2) xord(x1) yord(y1)','ordinal(x2 y2 x1 y1)'),
 'duplicates':('ordinal(x1 y1 x1) xord(x1 x1) yord(y1 y1)','ordinal(x1 y1)'),
 'bounded_disjoint':('bounded(x2) xord(x1) yord(y1)','bounded(x2) ordinal(x1 y1)'),
 'wildcards':('xord(x*) yord(y*)','ordinal(x1 x2 y1 y2)'),
}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--baseline-ado',type=Path);ap.add_argument('--only')
    args=ap.parse_args()
    if args.baseline_ado and not (args.baseline_ado/'dea.ado').is_file():ap.error('Invalid pre-F10 ado snapshot')
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'current_ado';pkg.mkdir();hashes={}
    files=package_files(REPO)
    for name in files:
        raw=(REPO/name).read_bytes();(pkg/Path(name).name).write_bytes(raw);hashes[name]=hashlib.sha256(raw).hexdigest()
    retry=out/'retry_ado';instrument(pkg,retry)
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
            if p.returncode or 'TEST_PASS' not in log.splitlines():raise RuntimeError(log[-2500:])
            row['passed']=True
        except Exception as e:row['error']=repr(e)
        row['seconds']=round(time.monotonic()-start,3);results.append(row)
        (out/'results.json').write_text(json.dumps(results,indent=2));print(json.dumps(row),flush=True)

    for i,((kind,(opts,canonical)),rts,ori,engine) in enumerate(itertools.product(VARIANTS.items(),['crs','vrs','nirs','ndrs'],['input','output'],['simplex','ipm','auto'])):
        name=f'grid_{kind}_{rts}_{ori}_{engine}'
        setup=DATA
        if i%5==0:setup+='gsort -id\n'
        if i%7==0:setup+='replace x2=.a in 5\n'
        condition=['1','tag=="keep = this, yes"','inrange(_n,2,7)','mod(id,3)!=0 & inrange(_n,2,7)'][i%4]
        qualifier=['','if tag=="keep = this, yes"','in 2/7','if mod(id,3)!=0 in 2/7'][i%4]
        # Missing numeric bounded data remain a separate F14 concern; the grid
        # puts missing values only in ordinary numeric dimensions.
        if kind=='bounded_disjoint' and i%7==0:setup+='replace x2=3 in 5\n'
        setup+=f'gen byte selected=({condition}) & !missing(x1,x2,y1,y2)\n'
        setup+='local savevars id x1 x2 y1 y2\nquietly ds\nlocal original_vars "`r(varlist)\'"\ntempfile original\nsave `original\'\n'
        base=f'rts({rts}) orientation({ori}) engine({engine}) stage({1+i%2}) dmu(id)'
        setup+=f'dea x1 x2=y1 y2 if selected, {canonical} {base} gen(control)\nmatrix C=r(efficiency)\n'
        setup+='foreach v in CRS_TE VRS_TE SCALE RTS_NUM RTS {\n gen str10 `v\'="untouched"\n}\n'
        replace=i%3==0
        if replace:setup+='foreach v in res res_opt res_mid {\n gen str10 `v\'="old"\n}\n'
        model=['','model(te)','model(imprecise)'][i%3];requested='imprecise' if i%3==2 else 'te'
        prefix=f'dea x1 x2=y1 y2 {qualifier},' if i%2 else f'dea y1 y2 {qualifier}, inputs(x1 x2)'
        setup+=f'{prefix} {opts} {model} {base} gen(res) saving("saved.dta")'+(' replacevars' if replace else '')+'\n'
        setup+=checks('imprecise',requested,'legacy','res res_opt res_mid',{'nirs':'drs','ndrs':'irs'}.get(rts,rts),saving=True)
        setup+='mata: assert(all(st_matrix("E"):==st_matrix("C")))\nassert res==control\ncf `original_vars\' using `original\', all\n'
        setup+='foreach v in CRS_TE VRS_TE SCALE RTS_NUM RTS {\n assert `v\'=="untouched"\n}\n'
        run(name,setup)

    # Assert the full data state, including labels, characteristics, order and
    # an existing save destination, is intact after every unsupported request.
    def reject(opts,form='A',replace=False,autoid=False,setup='',prefix=None,rc=198):
        s=DATA+setup+'''gen str12 res="old score"
gen str12 res_opt="old opt"
gen str12 res_mid="old mid"
label variable res "original label"
char res[review] "original characteristic"
label define groups 0 "no" 1 "yes"
label values w1 groups
char _dta[review] "original dataset characteristic"
tempfile original
save `original'
save "target.dta"
'''
        if prefix is None:prefix='dea x1 x2=y1 y2,' if form=='A' else 'dea y1 y2, inputs(x1 x2)'
        cmd=f'{prefix} {opts} '+('' if autoid else 'dmu(id) ')+('replacevars ' if replace else '')+'gen(res) saving("target.dta") replace'
        s+='capture noisily '+cmd+f'\nassert _rc=={rc}\n'
        s+='''cf _all using `original', all
assert "`:variable label res'"=="original label"
assert "`:char res[review]'"=="original characteristic"
assert "`:char _dta[review]'"=="original dataset characteristic"
assert "`:value label w1'"=="groups"
assert "`:label groups 1'"=="yes"
capture confirm variable _dmu,exact
assert _rc==111
preserve
use "target.dta",clear
cf _all using `original', all
restore
'''
        return s

    for i,(bound,legacy,model,form,replace) in enumerate(itertools.product(
        ['xlo(x1 x2)','xhi(x1 x2)','ylo(y1 y2)','yhi(y1 y2)'],
        ['bounded(x1)','ordinal(x1)','xord(x1)','yord(y1)'],
        ['','model(te)','model(imprecise)'],['A','B'],[False,True])):
        # Full lists containing point placeholders now use the legacy path.
        setup=DATA+'gen byte selected=1\nlocal savevars id x1 x2 y1 y2\n'
        setup+=f'dea x1 x2=y1 y2, {legacy} dmu(id) gen(control)\n'
        if replace:setup+='foreach v in res res_opt res_mid {\n gen str10 `v\'="old"\n}\n'
        prefix='dea x1 x2=y1 y2,' if form=='A' else 'dea y1 y2, inputs(x1 x2)'
        setup+=f'{prefix} {bound} {legacy} {model} dmu(id) gen(res) saving("saved.dta")'+(' replacevars' if replace else '')+'\n'
        setup+=checks('imprecise','imprecise' if model=='model(imprecise)' else 'te','legacy','res res_opt res_mid',saving=True)
        setup+='assert res==control\n'
        run(f'mixed_point_{i:03d}',setup)

    models=['fdh','model(fdh)','model(cost) w(w1 w2)','model(revenue) p(p1 p2)',
            'model(profit) w(w1 w2) p(p1 p2)','model(allocative) w(w1 w2)',
            'model(ddf)','model(sbm)','model(supereff)','model(supersbm)',
            'model(multiplier)','model(additive)','model(congestion)',
            'model(undesirable) bad(bad)','model(sbmundesirable) bad(bad)']
    options=['xlo(x1 x2)','xhi(x1 x2)','ylo(y1 y2)','yhi(y1 y2)','bounded(x1)','ordinal(x1)','xord(x1)','yord(y1)']
    for i,(model,option,form) in enumerate(itertools.product(models,options,['A','B'])):
        run(f'reject_model_{i:03d}',reject(model+' '+option,form,i%2==0,i%3==0))

    invalid=['xord(y1)','yord(x1)','xord(z)','yord(z)','ordinal(z)','bounded(z)',
             'xord(x1 z)','yord(y1 z)','ordinal(x1 z)','bounded(x1 z)',
             'bounded(x1) ordinal(x1)','bounded(x1) xord(x1)',
             'bounded(y1) yord(y1)','ordinal(x1) xord(y1)',
             'model(imprecise)','xlo(x1)','xhi(x1)','ylo(y1)','yhi(y1)']
    for i,(opts,form,replace) in enumerate(itertools.product(invalid,['A','B'],[False,True])):
        run(f'reject_role_{i:03d}',reject(opts,form,replace,i%2==0))
    for opt in ['xord(x1)','yord(x1)','ordinal(x1)','bounded(x1)']:
        for form in ['A','B']:
            prefix='dea x1 x2=x1 y1,' if form=='A' else 'dea x1 y1, inputs(x1 x2)'
            run('reject_ambiguous_'+re.sub('[^a-z]','',opt)+'_'+form,reject(opt,form,True,True,prefix=prefix))
    run('reject_string_ordinal',reject('ordinal(x1)','A',True,True,setup='tostring x1,replace\n'))

    for option,extra in itertools.product(['xord(x1)','yord(y1)','xord(x1) yord(y1)'],['bootstrap','reference saveref','contractcheck','saveweights']):
        run('reject_extra_'+re.sub('[^a-z]','',option+extra),reject(option+' '+extra,replace=True,autoid=True))

    # Existing-output and model-variable protections must include _opt/_mid
    # now that xord/yord alone really dispatch the legacy transformation.
    for kind,opts in [('x','xord(x1)'),('y','yord(y1)'),('both','xord(x1) yord(y1)')]:
        for suffix in ['','_opt','_mid']:
            target='res'+suffix
            setup=DATA+f'gen str10 {target}="untouched"\n'
            setup+=f'capture noisily dea x1 x2=y1 y2, {opts} gen(res) dmu(id)\nassert _rc==110\nassert {target}=="untouched"\n'
            run('reject_existing_'+kind+(suffix or '_main'),setup)
            # Keep the data-role collision valid before checking output protection.
            changed=opts.replace('x1',target)
            setup=DATA+f'rename x1 {target}\ntempfile original\nsave `original\'\n'
            setup+=f'capture noisily dea {target} x2=y1 y2, {changed} gen(res) dmu(id) replacevars\nassert _rc==110\ncf _all using `original\', all\n'
            run('protect_input_'+kind+(suffix or '_main'),setup)

    for opts,canonical in [('xord(x1)','ordinal(x1)'),('yord(y1)','ordinal(y1)'),('xord(x1) yord(y1)','ordinal(x1 y1)')]:
        for ori in ['input','output']:
            name='retry_'+re.sub('[^a-z]','',opts)+'_'+ori
            setup=DATA+'gen byte selected=mod(id,3)!=0\nlocal savevars id x1 x2 y1 y2\n'
            setup+=f'dea x1 x2=y1 y2 if selected, {canonical} engine(ipm) orientation({ori}) gen(control) dmu(id)\n'
            setup+=f'dea x1 x2=y1 y2 if selected, {opts} engine(auto) orientation({ori}) gen(res) dmu(id) saving("saved.dta")\n'
            setup+='assert "`r(engine)\'"=="ipm"\n'+checks('imprecise','te','legacy','res res_opt res_mid',saving=True)+'assert res==control\n'
            run(name,setup,retry)

    # Order information should be invariant to strictly increasing relabelling.
    for opts,variables in [('xord(x1)','x1'),('yord(y1)','y1'),('xord(x1) yord(y1)','x1 y1')]:
        for engine,ori in itertools.product(['simplex','ipm','auto'],['input','output']):
            name='monotone_'+re.sub('[^a-z]','',opts)+'_'+engine+'_'+ori
            setup=DATA+f'dea x1 x2=y1 y2, {opts} engine({engine}) orientation({ori}) gen(control) dmu(id)\n'
            setup+=f'foreach v in {variables} {{\n replace `v\'=19+`v\'^3\n}}\n'
            setup+=f'dea x1 x2=y1 y2, {opts} engine({engine}) orientation({ori}) gen(res) dmu(id)\nassert res==control\n'
            run(name,setup)

    # Numeric/string bounded dimensions can coexist with distinct ordinal ones.
    for engine,ori in itertools.product(['simplex','ipm','auto'],['input','output']):
        name=f'bounded_string_{engine}_{ori}'
        setup=DATA+'tostring x2,replace\nreplace x2="["+x2+","+string(real(x2)*1.2)+"]"\ngen byte selected=mod(id,3)!=0\nlocal savevars id x1 x2 y1 y2\n'
        setup+=f'dea x1 x2=y1 y2 if selected, bounded(x2) ordinal(x1 y1) engine({engine}) orientation({ori}) dmu(id) gen(control)\n'
        setup+=f'dea x1 x2=y1 y2 if selected, bounded(x2) xord(x1) yord(y1) engine({engine}) orientation({ori}) dmu(id) gen(res) saving("saved.dta")\n'
        setup+=checks('imprecise','te','legacy','res res_opt res_mid',saving=True)+'assert res==control\n'
        run(name,setup)

    # Original defect and duplicate alias transformations, using a small fixture
    # where ignoring the ordering or transforming twice changes the raw scores.
    for kind in ['x','y']:
        opt='xord(x)' if kind=='x' else 'yord(y)';var='x' if kind=='x' else 'y'
        common=BASIC+('replace y=5-y\n' if kind=='y' else '')+f'dea x=y, ordinal({var}) dmu(id) gen(control)\n'
        for engine,ori in itertools.product(['simplex','ipm','auto'],['input','output']):
            name=f'original_{kind}_{engine}_{ori}'
            if not chosen(name):continue
            cmd=f'dea x=y, {opt} engine({engine}) orientation({ori}) dmu(id) gen(res)\n'
            # Match the control's orientation/engine, and isolate temporary results.
            base=BASIC+f'dea x=y, ordinal({var}) engine({engine}) orientation({ori}) dmu(id) gen(control)\n'
            if args.baseline_ado:
                old=base+cmd+'assert "`r(model)\'"=="te"\n'
                old+='gen double gap=abs(res-control)\nsummarize gap,meanonly\nassert r(max)>.01\n'
                old+=f'dea x=y, engine({engine}) orientation({ori}) dmu(id) gen(te)\nassert res==te\n'
                run('old_'+name,old,args.baseline_ado,force=True)
            run(name,base+cmd+'assert "`r(model)\'"=="imprecise"\nassert res==control\n')
        name='duplicate_'+kind
        if chosen(name):
            cmd=f'dea x=y, ordinal({var}) {opt} dmu(id) gen(res)\n'
            if args.baseline_ado:
                old=common+cmd+'gen double gap=abs(res-control)\nsummarize gap,meanonly\nassert r(max)>.01\n'
                run('old_'+name,old,args.baseline_ado,force=True)
            run(name,common+cmd+'assert res==control\n')

    # F30 deliberately changes ordinal-output scores. Preserve historical
    # equality only for unaffected input-only transforms; output transforms are
    # checked against independent oracles in test_ordinal_order.py.
    if args.baseline_ado:
        for kind,opts in [('x','xord(x1)'),('bounded_input','bounded(x2) xord(x1)')]:
            for engine,ori,rts in itertools.product(['simplex','ipm','auto'],['input','output'],['crs','vrs']):
                name=f'history_{kind}_{engine}_{ori}_{rts}'
                if not chosen(name):continue
                base=DATA+f'dea x1 x2=y1 y2 if id!=3 in 2/8, model(imprecise) {opts} engine({engine}) orientation({ori}) rts({rts}) dmu(id) gen(res)\n'
                run('old_'+name,base+'keep res res_opt res_mid\nsave "expected.dta"\n',args.baseline_ado,force=True)
                path=out/('old_'+name)/'expected.dta'
                run(name,base+f'cf res res_opt res_mid using "{path}",all\n')

    summary=dict(total=len(results),passed=sum(r['passed'] for r in results),failed=sum(not r['passed'] for r in results))
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
    raise SystemExit(summary['failed']>0)


if __name__=='__main__':main()
