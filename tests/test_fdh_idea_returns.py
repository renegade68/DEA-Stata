#!/usr/bin/env python3
"""FDH/IDEA result contracts, dispatch identity, saving and historical scores.

Requires Python's standard library and licensed Stata 16+. Test-only fault
injection checks missing-result bookkeeping, not failure probabilities. No
production solver or installed ado is changed. Optional pre-F09 scores must
match bit for bit; this does not certify the imprecise-data methodology.
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

REPO = Path(__file__).resolve().parents[1]
DATA = '''set obs 9
gen long id=_n
gen double x1=1+mod(id*7,9)
gen double x2=1+mod(id*5,11)
gen double y1=1+mod(id*3,7)
gen double y2=1+mod(id*2,5)
foreach v of varlist x1 x2 y1 y2 {
 gen double `v'lo=.8*`v'
 gen double `v'hi=1.2*`v'
}
'''
BOUNDS = ['xlo(x1 x2) xhi(x1 x2)', 'xlo(x1lo x2lo)', 'xhi(x1hi x2hi)',
          'ylo(y1lo y2lo)', 'yhi(y1hi y2hi)',
          'xlo(x1lo x2lo) xhi(x1hi x2hi)',
          'ylo(y1lo y2lo) yhi(y1hi y2hi)',
          'xlo(x1lo x2lo) xhi(x1hi x2hi) ylo(y1lo y2lo) yhi(y1hi y2hi)']
LEGACY = {'bounded_numeric': ('bounded(x)', ''),
          'bounded_string': ('bounded(x)', 'tostring x, replace\nreplace x="["+x+","+string(real(x)*1.2)+"]"\n'),
          'ordinal_input': ('ordinal(x)', ''),
          'ordinal_output': ('ordinal(y)', ''),
          'xord_explicit': ('model(imprecise) xord(x)', ''),
          'yord_explicit': ('model(imprecise) yord(y)', '')}


def checks(model, requested, kind='', outputs='res', expected_rts='crs', saving=False):
    s='''matrix E=r(efficiency)
foreach k in N N_eff mean_eff min_eff max_eff n_failed tolerance {
 local got_`k'=r(`k')
}
'''
    s+=f'''assert "`r(model)'"=="{model}"
assert "`r(model_requested)'"=="{requested}"
assert "`r(rts)'"=="{expected_rts}"
assert "`r(idea_method)'"=="{kind}"
'''
    if kind:
        s+=f'''assert "`r(idea_vars)'"=="{outputs}"
foreach k in N_eff_opt mean_opt min_opt max_opt n_failed_primary n_failed_opt {{
 local got_`k'=r(`k')
}}
'''
    else:
        s+='assert "`r(idea_vars)\'"==""\nassert missing(r(mean_opt))\n'
    s+='di "F09_METRICS " `got_N\' " " `got_n_failed\' " " `got_N_eff\'\n'
    if kind:s+='di "F09_AUX " `got_n_failed_primary\' " " `got_n_failed_opt\'\n'
    s+='''capture matrix list r(vrsfrontier)
assert _rc==111
capture matrix list r(vrsobs)
assert _rc==111
count if selected
assert r(N)==`got_N'
assert rowsof(E)==`got_N' & colsof(E)==1
mata: assert(all(st_matrix("E"):==st_data(selectindex(st_data(.,"selected")),"res")))
count if selected & abs(res-1)<=`got_tolerance'
assert r(N)==`got_N_eff'
quietly summarize res if selected, meanonly
assert (missing(r(mean)) & missing(`got_mean_eff')) | abs(r(mean)-`got_mean_eff')<1e-12
assert (missing(r(min)) & missing(`got_min_eff')) | abs(r(min)-`got_min_eff')<1e-12
assert (missing(r(max)) & missing(`got_max_eff')) | abs(r(max)-`got_max_eff')<1e-12
assert missing(res) if !selected
'''
    if kind:
        s+='''count if selected & missing(res,res_opt)
assert r(N)==`got_n_failed'
count if selected & missing(res)
assert r(N)==`got_n_failed_primary'
count if selected & missing(res_opt)
assert r(N)==`got_n_failed_opt'
count if selected & abs(res_opt-1)<=`got_tolerance'
assert r(N)==`got_N_eff_opt'
quietly summarize res_opt if selected, meanonly
assert (missing(r(mean)) & missing(`got_mean_opt')) | abs(r(mean)-`got_mean_opt')<1e-12
assert (missing(r(min)) & missing(`got_min_opt')) | abs(r(min)-`got_min_opt')<1e-12
assert (missing(r(max)) & missing(`got_max_opt')) | abs(r(max)-`got_max_opt')<1e-12
assert missing(res_opt) if !selected
capture confirm variable res_pess, exact
assert _rc==111
'''
        if kind=='legacy':s+='assert res==res_opt & res==res_mid\n'
        else:s+='capture confirm variable res_mid, exact\nassert _rc==111\n'
    else:s+='count if selected & missing(res)\nassert r(N)==`got_n_failed\'\n'
    if saving:
        s+=f'''preserve
keep if selected
keep `savevars' {outputs}
tempfile expected
save `expected'
use "saved.dta",clear
cf _all using `expected', all
restore
'''
    return s


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--baseline-ado',type=Path);ap.add_argument('--only')
    args=ap.parse_args()
    if args.baseline_ado and not (args.baseline_ado/'dea.ado').is_file():ap.error('Invalid pre-F09 snapshot')
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    pkg=out/'current_ado';pkg.mkdir();hashes={}
    files=package_files(REPO)
    for name in files:
        raw=(REPO/name).read_bytes();(pkg/Path(name).name).write_bytes(raw);hashes[name]=hashlib.sha256(raw).hexdigest()
    retry=out/'retry_ado';instrument(pkg,retry)
    fault=out/'fault_ado';shutil.copytree(pkg,fault)
    p=fault/'dea.ado';s=p.read_text()
    for program,anchor,vars_ in [('_dea_fdh_dispatch','    // Copy to final variable',['efficiency']),
                                 ('_dea_imprecise_dispatch','    // Create output variables',['eff_pess','eff_opt']),
                                 ('_dea_imprecise_legacy',"    qui gen double `generate' = `__idea_eff' if `touse'",['__idea_eff'])]:
        a=s.index('program define '+program+', rclass');b=s.index('\nend\n',a)+5;t=s[a:b];assert t.count(anchor)==1
        inj='    // TEST COPY ONLY: deterministic missing results / frontier ties.\n'
        for j,v in enumerate(vars_):
            name='`'+v+"'"
            inj+=f'''    if "$F09_CASE" == "all" quietly replace {name}=.z if `touse'
    if "$F09_CASE" == "mixed" {{
        quietly replace {name}=.a if `touse' & _n=={2+j}
        quietly replace {name}=.z if `touse' & _n==4
    }}
    if "$F09_CASE" == "ties" {{
        quietly replace {name}=1-`tolerance' if `touse' & _n==1
        quietly replace {name}=1+`tolerance' if `touse' & _n==2
        quietly replace {name}=1-2*`tolerance' if `touse' & _n==3
        quietly replace {name}=.a if `touse' & _n==4
    }}
'''
        s=s[:a]+t.replace(anchor,inj+anchor)+s[b:]
    p.write_text(s)
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,stata=args.stata,baseline=str(args.baseline_ado)),indent=2))
    results=[]
    def selected(name):return not args.only or re.search(args.only,name)
    def run(name,script,source=pkg,force=False,oracle=None):
        if not force and not selected(name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 120\n'
        for folder in ['PLUS','PERSONAL','SITE','OLDPLACE']:
            (d/folder).mkdir();head+=f'sysdir set {folder} "{d/folder}"\n'
        head+=f'adopath ++ "{source.resolve()}"\n';(d/'test.do').write_text(head+script+'\ndi "TEST_PASS"\n')
        start=time.monotonic();row=dict(test=name,passed=False)
        try:
            p=subprocess.run([args.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=240)
            log=(d/'test.log').read_text(errors='replace')
            if p.returncode or 'TEST_PASS' not in log.splitlines():raise RuntimeError(log[-2300:])
            metrics=re.findall(r'^F09_METRICS (.*)$',log,re.M)
            aux=re.findall(r'^F09_AUX (.*)$',log,re.M)
            if metrics:row['metrics']=[float(x) for x in metrics[-1].split()]
            if aux:row['missing_primary_opt']=[int(x) for x in aux[-1].split()]
            if oracle:oracle(d)
            row['passed']=True
        except Exception as e:row['error']=repr(e)
        row['seconds']=round(time.monotonic()-start,3);results.append(row)
        (out/'results.json').write_text(json.dumps(results,indent=2));print(json.dumps(row),flush=True)

    # Engines and stages are accepted FDH selectors; dominance scores ignore both.
    for i,(spelling,engine,ori,rts) in enumerate(itertools.product(['fdh','model(fdh)'],['simplex','ipm','auto'],['input','output'],['','crs','vrs'])):
        name=f'fdh_{i:03d}'
        setup=DATA+'gen byte selected=inrange(_n,2,8) & id!=4\nlocal savevars id x1 x2 y1 y2\n'
        prefix='dea y1 y2 if selected, inputs(x1 x2)' if i%2 else 'dea x1 x2=y1 y2 if id!=4 in 2/8,'
        setup+=prefix+f' {spelling} engine({engine}) orientation({ori}) stage({1+i%2}) '+(f'rts({rts}) ' if rts else '')+'dmu(id) gen(res) saving("saved.dta")\n'
        setup+=checks('fdh','te' if spelling=='fdh' else 'fdh',expected_rts='vrs',saving=True)
        # Exhaustive positive-data dominance oracle (independent Stata/Mata code).
        setup+=f'''mata:
X=st_data(.,tokens("x1 x2"),"selected");Y=st_data(.,tokens("y1 y2"),"selected")
E=st_data(.,"res","selected"); answer=J(rows(X),1,.)
for(k=1;k<=rows(X);k++) {{
 if ("{ori}"=="input") {{
  candidates=selectindex(rowmin(Y:-Y[k,.]):>=0)
  answer[k]=min(rowmax(X[candidates,.]:/X[k,.]))
 }}
 else {{
  candidates=selectindex(rowmax(X:-X[k,.]):<=0)
  answer[k]=1/max(rowmin(Y[candidates,.]:/Y[k,.]))
 }}
}}
assert(max(abs(answer-E))<1e-12)
end
'''
        run(name,setup)

    for i,(bounds,rts,ori,engine) in enumerate(itertools.product(BOUNDS,['crs','vrs','nirs','ndrs'],['input','output'],['simplex','ipm','auto'])):
        name=f'interval_{i:03d}'
        explicit=(i%2==0);requested='imprecise' if explicit else 'te'
        setup=DATA+'gen byte selected=inrange(_n,2,8) & id!=4\nlocal savevars id x1 x2 y1 y2\n'
        # Existing decomposition fields belong to the user, even with VRS IDEA.
        setup+='foreach v in CRS_TE VRS_TE SCALE RTS_NUM RTS {\n gen str10 `v\'="untouched"\n}\n'
        if i%5==0:setup+='gen str10 res="old"\ngen str10 res_opt="old"\n'
        prefix='dea y1 y2 if selected, inputs(x1 x2)' if i%2 else 'dea x1 x2=y1 y2 if id!=4 in 2/8,'
        setup+=prefix+f' {bounds} '+('model(imprecise) ' if explicit else '')+f'rts({rts}) orientation({ori}) engine({engine}) stage({1+i%2}) dmu(id) gen(res) saving("saved.dta")'+(' replacevars' if i%5==0 else '')+'\n'
        setup+=checks('imprecise',requested,'interval','res res_opt',{'nirs':'drs','ndrs':'irs'}.get(rts,rts),saving=True)
        setup+='foreach v in CRS_TE VRS_TE SCALE RTS_NUM RTS {\n assert `v\'=="untouched"\n}\n'
        setup+='preserve\nkeep if selected\nformat res res_opt %24.17g\nexport delimited id x1 x2 y1 y2 res res_opt using "scores.csv",replace datafmt\nrestore\n'
        run(name,setup)

    for i,((kind,(opts,change)),rts,ori,engine) in enumerate(itertools.product(LEGACY.items(),['crs','vrs'],['input','output'],['simplex','ipm','auto'])):
        name=f'legacy_{kind}_{rts}_{ori}_{engine}'
        explicit='model(' in opts or i%2==0
        if explicit and 'model(' not in opts:opts='model(imprecise) '+opts
        setup=BASIC+change+'gen byte selected=_n!=2\nlocal savevars id x y\n'
        setup+=f'dea x=y if selected, {opts} rts({rts}) orientation({ori}) engine({engine}) stage({1+i%2}) dmu(id) gen(res) saving("saved.dta")\n'
        setup+=checks('imprecise','imprecise' if explicit else 'te','legacy','res res_opt res_mid',rts,saving=True)
        run(name,setup)

    for kind,case,engine in itertools.product(['fdh','interval','legacy'],['mixed','all','ties'],['simplex','ipm','auto']):
        opts={'fdh':'fdh','interval':'xlo(x) xhi(x)','legacy':'bounded(x)'}[kind]
        setup=BASIC+'gen byte selected=1\nlocal savevars id x y\n'+f'global F09_CASE {case}\n'
        # 1/8 is exactly representable, so equality tests cover the frontier edge.
        setup+=f'dea x=y, {opts} engine({engine}) dmu(id) gen(res) tolerance(.125) saving("saved.dta")\n'
        setup+='assert "`r(engine)\'"=="'+('simplex' if engine=='auto' and kind!='legacy' else engine if engine!='auto' else 'ipm')+'"\n'
        setup+=checks('fdh' if kind=='fdh' else 'imprecise','te','' if kind=='fdh' else kind,
                      {'fdh':'res','interval':'res res_opt','legacy':'res res_opt res_mid'}[kind],
                      'vrs' if kind=='fdh' else 'crs',saving=True)
        if case=='all':setup+='assert `got_N_eff\'==0 & `got_n_failed\'==4\n'
        elif case=='ties':setup+='assert `got_N_eff\'==2 & `got_n_failed\'==1\n'
        else:setup+=f'assert `got_n_failed\'=={3 if kind=="interval" else 2}\n'
        run(f'fault_{kind}_{case}_{engine}',setup,fault)

    for ori in ['input','output']:
        setup=BASIC+'gen byte selected=1\nlocal savevars id x y\n'
        setup+=f'dea x=y, bounded(x) engine(ipm) orientation({ori}) dmu(id) gen(control)\n'
        setup+=f'dea x=y, bounded(x) engine(auto) orientation({ori}) dmu(id) gen(res) saving("saved.dta")\n'
        setup+='assert "`r(engine)\'"=="ipm"\nassert r(n_failed)==0\n'
        setup+=checks('imprecise','te','legacy','res res_opt res_mid',saving=True)+'assert res==control\n'
        run('retry_legacy_'+ori,setup,retry)

    for spelling,rts in itertools.product(['fdh','model(fdh)'],['drs','irs','nirs','ndrs']):
        setup=BASIC+'gen str10 res="original"\n'
        setup+=f'capture noisily dea x=y, {spelling} rts({rts}) gen(res) replacevars\nassert _rc==198\nassert res=="original"\n'
        setup+='capture confirm variable _dmu, exact\nassert _rc==111\n'
        run('reject_fdh_'+str(spelling=='fdh')+'_'+rts,setup)

    for options in ['xlo(x)','xhi(x)','ylo(y)','yhi(y)','ordinal(x)','bounded(x)','model(imprecise)']:
        run('reject_contract_'+re.sub('[^a-z]','',options),BASIC+f'capture noisily dea x=y, {options} contractcheck\nassert _rc==198\n')

    for kind,opts in [('fdh','fdh'),('interval','xlo(x) xhi(x)'),('legacy','ordinal(x)')]:
        setup=BASIC+'gen byte selected=1\n'+f'dea x=y, {opts} dmu(id) gen(first)\n'
        setup+='dea x=y, dmu(id) gen(res)\nassert "`r(model)\'"=="te"\nassert "`r(model_requested)\'"=="te"\n'
        setup+='assert "`r(idea_method)\'`r(idea_vars)\'"==""\nassert missing(r(N_eff_opt),r(mean_opt),r(n_failed_opt))\n'
        run('stale_'+kind,setup)

    for kind,opts in [('fdh','fdh'),('interval','xlo(x) xhi(x)'),('legacy','bounded(x)')]:
        for identity in ['automatic','duplicate','string']:
            setup=BASIC+'gen byte selected=_n!=2\n'
            idopt='dmu(id)'
            if identity=='automatic':
                idopt='';setup+='local savevars _dmu x y\n'
            else:
                setup+='local savevars id x y\n'
                setup+='replace id=1\n' if identity=='duplicate' else 'tostring id,replace\nreplace id="기관 "+id\n'
            outputs={'fdh':'res','interval':'res res_opt','legacy':'res res_opt res_mid'}[kind]
            setup+=f'dea x=y if selected using "saved.dta", {opts} {idopt} gen(res)\n'
            setup+=checks('fdh' if kind=='fdh' else 'imprecise','te','' if kind=='fdh' else kind,outputs,'vrs' if kind=='fdh' else 'crs',saving=True)
            setup+="foreach v in CRS_TE VRS_TE SCALE RTS_NUM RTS {\n capture confirm variable `v',exact\n assert _rc==111\n}\n"
            # Replacing a disk file and replacing in-memory columns are separate.
            setup+=f'dea x=y if selected, {opts} '+('dmu(_dmu)' if identity=='automatic' else idopt)+' gen(res) replacevars saving("saved.dta") replace\n'
            setup+=checks('fdh' if kind=='fdh' else 'imprecise','te','' if kind=='fdh' else kind,outputs,'vrs' if kind=='fdh' else 'crs',saving=True)
            run('identity_'+kind+'_'+identity,setup)

    if args.baseline_ado:
        for engine in ['simplex','auto']:
            # A fresh old session computes IDEA, then crashes in the spurious
            # VRS decomposition because its simplex helper has not been loaded.
            name='cold_interval_vrs_'+engine
            if not selected(name):continue
            base=BASIC+'gen byte selected=1\n'
            cmd=f'dea x=y, xlo(x) xhi(x) rts(vrs) engine({engine}) dmu(id) gen(res)'
            run('old_'+name,base+'capture noisily '+cmd+'\nassert _rc==3499\nassert !missing(res,res_opt)\n',args.baseline_ado,force=True)
            run(name,base+cmd+'\n'+checks('imprecise','te','interval','res res_opt','vrs'))

    # Each endpoint uses a consistent realization, including its own reference
    # row. A singleton therefore has lower=upper=1. Preserve the pre-F09 control.
    for rts,ori in itertools.product(['crs','vrs'],['input','output']):
        name=f'singleton_formulation_{rts}_{ori}'
        if not selected(name):continue
        setup='input id x y\n1 1 1\nend\ngen byte selected=1\ngen double lo=.8\ngen double hi=1.2\n'
        setup+=f'dea x=y, model(imprecise) xlo(lo) xhi(hi) ylo(lo) yhi(hi) rts({rts}) orientation({ori}) dmu(id) gen(res)\n'
        expected=4/9 if rts=='crs' else 2/3
        equations=f'assert abs(res-{expected:.17g})<1e-7\n'
        equations+=('assert abs(res_opt-2.25)<1e-7\n' if rts=='crs' and ori=='input' else 'assert missing(res_opt)\n')
        if args.baseline_ado:run('old_'+name,setup+equations,args.baseline_ado,force=True)
        run(name,setup+checks('imprecise','imprecise','interval','res res_opt',rts)+'assert abs(res-1)<1e-7 & abs(res_opt-1)<1e-7\n')

    profiles={'fdh':('model(fdh)','res'), 'interval':('model(imprecise) xlo(x) xhi(x)','res res_opt'),
              'implicit':('xlo(x) xhi(x)','res res_opt'), 'legacy':('ordinal(x)','res res_opt res_mid')}
    if args.baseline_ado:
        for kind,(opts,cols) in profiles.items():
            for engine,ori,rts in itertools.product(['simplex','ipm','auto'],['input','output'],['crs','vrs']):
                name=f'history_{kind}_{engine}_{ori}_{rts}'
                if not selected(name):continue
                base=BASIC+f'dea x=y, {opts} engine({engine}) orientation({ori}) rts({rts}) dmu(id) gen(res)\n'
                # Historical implicit VRS fails in a fresh session unless the
                # unused simplex decomposition backend was already loaded.
                old='quietly run "'+str(args.baseline_ado.resolve()/'_dea_simplex.ado')+'"\n'+base
                if kind in ['fdh','interval','implicit']:old+='assert missing(r(mean_eff),r(n_failed))\n'
                if kind=='implicit':old+='assert "`r(model)\'"=="te"\n'
                old+=f'keep {cols}\nsave "expected.dta",replace\n'
                run('old_'+name,old,args.baseline_ado,force=True)
                path=out/('old_'+name)/'expected.dta'
                comparison=f'cf {cols} using "{path}", all\n'
                if kind in ['interval','implicit']:
                    # Point-interval equations are unchanged, but R09 tightens
                    # convergence. The independent formulation suite checks
                    # accuracy; this historical control allows solver rounding.
                    comparison=f'''assert !missing(res,res_opt)
preserve
use "{path}",clear
assert !missing(res,res_opt)
mkmat res res_opt,matrix(OLD)
restore
mata: assert(max(abs(st_data(.,tokens("res res_opt"))-st_matrix("OLD")))<1e-6)
'''
                run(name,base+'assert !missing(r(mean_eff),r(n_failed))\n'+comparison)
    summary=dict(total=len(results),passed=sum(r['passed'] for r in results),failed=sum(not r['passed'] for r in results))
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
    raise SystemExit(summary['failed']>0)


if __name__=='__main__':main()
