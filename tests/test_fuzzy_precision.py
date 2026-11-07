#!/usr/bin/env python3
"""Fuzzy LP precision, units, validation, and bounded failure recovery.

Requires NumPy/SciPy and Stata 16+. Instrumentation/faults affect isolated test
copies only. Independent HiGHS matrices are built from the unscaled fixture.
"""
from _support import package_files
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time

import numpy as np
import scipy
from test_fuzzy_rts import DATA, ALPHAS, SPECS, body, read, oracle
from test_output_protection import unchanged

REPO=Path(__file__).resolve().parents[1]
CMD='dea_fuzzy, inputs(x1 x2) outputs(y1 y2) xlo(lo_x1 lo_x2) xhi(hi_x1 hi_x2) ylo(lo_y1 lo_y2) yhi(hi_y1 hi_y2) dmu(id) alpha(0 .5 1) gen(res)'
FIXTURE=DATA+'keep if id>=2 & id<=7\n'
OLD=''.join(f'gen str12 res_{s}_a{a}="old result"\n' for s in ['pess','opt'] for a in [0,50,100])
PROBE_INIT='scalar probe_calls=0\nscalar probe_first=.\nscalar probe_min=.\n'
RECORD='''st_numscalar("probe_calls",st_numscalar("probe_calls")+1)
        if (st_numscalar("probe_calls")==1) st_numscalar("probe_first",use_tol)
        st_numscalar("probe_min",min((st_numscalar("probe_min"),use_tol)))
        value=S.optimize()'''


def compare(d, dim, orientation, mode, threshold=1e-6):
    data=read(d/'results.csv');ref,count=oracle(data,dim,orientation,mode)
    chosen=data['selected']==1;gap=0.;inversion=0.;nesting=0.
    for spec,norm in SPECS.items():
        for sense in ['pess','opt']:
            previous=None
            for alpha in ALPHAS:
                a=data[f'f_{spec}_{sense}_a{int(alpha*100)}'][chosen]
                b=ref[norm,sense,int(alpha*100)]
                assert np.isfinite(a).all()
                gap=max(gap,float(abs(a-b).max()))
                assert np.array_equal(a,data[f'f_{norm}_{sense}_a{int(alpha*100)}'][chosen])
                if previous is not None:
                    nesting=max(nesting,float(np.max(previous-a if sense=='pess' else a-previous)))
                previous=a
        for alpha in ALPHAS:
            inversion=max(inversion,float(np.max(data[f'f_{spec}_pess_a{int(alpha*100)}'][chosen]-data[f'f_{spec}_opt_a{int(alpha*100)}'][chosen])))
    assert gap<threshold,(gap,threshold)
    assert inversion<=1e-8 and nesting<=1e-8,(inversion,nesting)
    return dict(lp_count=count,max_oracle_gap=gap,max_interval_inversion=inversion,max_nesting_violation=nesting,oracle_threshold=threshold)


def grid_script(data, orientation, extra='', scaled=False):
    script=data+'gen byte selected=1\n'
    x='x1 x2';y='y1 y2'
    if scaled:
        for v,factor in [('x1','1e-100'),('x2','1e100'),('y1','1e100'),('y2','1e-100')]:
            for prefix in ['', 'lo_', 'hi_']:
                script+=f'gen double s_{prefix}{v}={prefix}{v}*{factor}\n'
        x='s_x1 s_x2';y='s_y1 s_y2'
    prefix='s_' if scaled else ''
    bounds=' '.join(f'{option}({prefix}{p}_{v}1 {prefix}{p}_{v}2)' for option,p,v in [('xlo','lo','x'),('xhi','hi','x'),('ylo','lo','y'),('yhi','hi','y')])
    for spec in SPECS:
        script+=f'dea_fuzzy, inputs({x}) outputs({y}) {bounds} rts({spec}) ort({orientation}) dmu(id) gen(f_{spec}) {extra}\n'
        script+='assert r(lp_missing)==0 & r(max_primal_violation)<=1e-8\n'
    return script+'export delimited using "results.csv", replace\n'


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--only');args=ap.parse_args()
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    pkg=out/'current_ado';pkg.mkdir();hashes={}
    names=package_files(REPO)
    for name in names:
        data=(REPO/name).read_bytes();(pkg/Path(name).name).write_bytes(data);hashes[name]=hashlib.sha256(data).hexdigest()
    source=(pkg/'dea_fuzzy.ado').read_text()
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,stata=args.stata,numpy=np.__version__,scipy=scipy.__version__),indent=2))
    results=[]

    def run(name,script,edited=None,check=None):
        if args.only and not re.search(args.only,name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 120\n'
        for folder in ['PLUS','PERSONAL','SITE','OLDPLACE']:
            (d/folder).mkdir();head+=f'sysdir set {folder} "{d/folder}"\n'
        head+=f'adopath ++ "{pkg}"\n'
        if edited is not None:
            (d/'injected').mkdir();(d/'injected/dea_fuzzy.ado').write_text(edited)
            head+=f'adopath ++ "{d/"injected"}"\n'
        (d/'test.do').write_text(head+script+'\ndi "TEST_PASS"\n')
        row=dict(test=name,passed=False,instrumented=edited is not None);start=time.monotonic()
        try:
            subprocess.run([args.stata,'-b','do','test.do'],cwd=d,timeout=180,capture_output=True)
            log=(d/'test.log').read_text(errors='replace')
            if 'TEST_PASS' not in log.splitlines():raise RuntimeError(log[-1800:])
            if check:row.update(check(d))
            row['passed']=True
        except Exception as e:row['error']=str(e)
        row['seconds']=round(time.monotonic()-start,3);results.append(row)
        (out/'results.json').write_text(json.dumps(results,indent=2));print(json.dumps(row),flush=True)

    # All four LP construction branches must pass the caller's tolerance.
    instrumented=source.replace('value=S.optimize()',RECORD)
    for tol in [None,'1e-3','1e-7','1e-10','1e-12']:
        option='' if tol is None else f'tolerance({tol})';expected=tol or '1e-10'
        for orientation in ['input','output']:
            run(f'tolerance_wiring_{orientation}_{tol}',FIXTURE+PROBE_INIT+CMD+f' ort({orientation}) {option}\n'+f'''assert r(tolerance)=={expected} & r(check_tolerance)==1e-8
assert probe_first=={expected}
assert probe_calls==r(lp_count)+r(lp_retries)
assert r(lp_missing)==0 & r(max_primal_violation)<=1e-8
''',instrumented)

    # Recheck the two published counterexamples at explicit tolerances.
    for dim,orientation,mode,sample in [(1,'input','outputs','ifin'),(2,'output','inputs','if')]:
        for tol in ['1e-7','1e-10','1e-12']:
            script=body(dim,orientation,mode,sample)
            script='\n'.join(line+f' tolerance({tol})' if line.startswith('dea_fuzzy ') else line for line in script.splitlines())+'\n'
            run(f'counterexample_{dim}_{orientation}_{tol}',script,check=lambda d,dim=dim,orientation=orientation,mode=mode:compare(d,dim,orientation,mode))

    for orientation in ['input','output']:
        run('units_'+orientation,grid_script(FIXTURE,orientation,scaled=True),check=lambda d,o=orientation:compare(d,2,o,'full'))
        for seed in [17,836,2901]:
            random=f'''set obs 20
set seed {seed}
gen id=_n
foreach v in x1 x2 y1 y2 {{
 gen double `v'=exp(3*runiform()-1.5)
 gen double lo_`v'=`v'*(.3+.6*runiform())
 gen double hi_`v'=`v'*(1.1+runiform())
}}
'''
            run(f'random_{orientation}_{seed}',grid_script(random,orientation),check=lambda d,o=orientation:compare(d,2,o,'full'))

    for tol in ['0','-1','.','.a']:
        run('invalid_tolerance_'+tol,FIXTURE+OLD+unchanged(CMD+f' replacevars tolerance({tol})',198))
    for variable,expr in [('lo_x1','-.1'),('lo_y1','-.1'),('lo_x2','x2+.1'),('hi_x2','x2-.1'),('lo_y2','y2+.1'),('hi_y2','y2-.1')]:
        run('invalid_bound_'+variable+'_'+expr.replace('.',''),FIXTURE+f'replace {variable}={expr} in 2\n'+OLD+unchanged(CMD+' replacevars',459))
        run('excluded_bound_'+variable+'_'+expr.replace('.',''),FIXTURE+f'replace {variable}={expr} in 2\n'+CMD.replace('dea_fuzzy,','dea_fuzzy if _n!=2,')+'\nassert r(N)==5\nassert missing(res_pess_a0) in 2\n')

    # Failures are injected after a real LP optimization; count calls to prove
    # bounded retry, and test the complete dataset/metadata on terminal errors.
    injections={
        'nonconvergence':('status=S.errorcode()','status=S.errorcode()\n        {condition} status=3'),
        'primal':('x=S.parameters()','x=S.parameters()\n            {condition} x[1]=-1'),
        'objective':('status=S.errorcode()','status=S.errorcode()\n        {condition} value=value+.1'),
        'missing_parameter':('x=S.parameters()','x=S.parameters()\n            {condition} x[1]=.'),
    }
    for kind,(needle,replacement) in injections.items():
        for persistent in [False,True]:
            text=instrumented.replace(needle,replacement.replace('{condition}','' if persistent else 'if (attempt==1)'))
            if persistent:
                script=FIXTURE+OLD+PROBE_INIT+unchanged(CMD+' replacevars',498)+'assert probe_calls==24\n'
            else:
                script=FIXTURE+PROBE_INIT+CMD+'\nassert r(lp_count)==36 & r(lp_retries)==36 & probe_calls==72\nassert r(lp_missing)==0\nassert probe_min<=1e-12\n'
            run(f'inject_{kind}_{persistent}',script,text)

    for status in [1,2]:
        text=instrumented.replace('status=S.errorcode()',f'status={status}')
        run(f'inject_consistent_status_{status}',FIXTURE+PROBE_INIT+CMD+f'\nassert r(lp_missing)==36 & r(lp_count)==36 & r(lp_retries)==36\nassert r(lp_infeasible)=={36 if status==1 else 0} & r(lp_unbounded)=={36 if status==2 else 0}\nassert missing(res_pess_a0,res_opt_a100)\n',text)
    text=instrumented.replace('status=S.errorcode()','status=attempt')
    run('inject_inconsistent_status',FIXTURE+OLD+PROBE_INIT+unchanged(CMD+' replacevars',498)+'assert probe_calls==24\n',text)

    for kind in ['interval','nesting']:
        for persistent in [False,True]:
            change='eff_pess[1]=eff_opt[1]+.01' if kind=='interval' else 'eff_pess[1]=(alpha==0 ? .9 : .1); eff_opt[1]=1'
            condition='if (1)' if persistent else 'if (tol>=1e-10)'
            text=source.replace('st_matrix(stats_name,lpstats)',f'{condition} {{\n        {change}\n    }}\n    st_matrix(stats_name,lpstats)')
            script=FIXTURE+OLD+unchanged(CMD+' replacevars',498) if persistent else FIXTURE+CMD+'\nassert r(grid_retries)==1 & r(lp_count)==72\nassert r(lp_missing)==0\n'
            run(f'inject_grid_{kind}_{persistent}',script,text)

    for orientation in ['input','output']:
        for rts in ['crs','vrs','drs','irs']:
            run(f'single_{orientation}_{rts}',FIXTURE+'keep in 1\n'+CMD+f' ort({orientation}) rts({rts})\nassert r(lp_missing)==0\nassert abs(res_pess_a0-1)<1e-8 & abs(res_opt_a100-1)<1e-8\n')
            # A redundant zero column is allowed when each row has another
            # positive quantity on that side; an all-zero side is tested below.
            zero=FIXTURE+'replace x2=0\nreplace lo_x2=0\nreplace hi_x2=0\n'
            run(f'zero_column_{orientation}_{rts}',zero+CMD+f' ort({orientation}) rts({rts})\nassert r(lp_missing)==0\n')
    # All-zero output rows are now rejected before publishing results.
    run('zero_output_domain','input id x1 y1\n1 1 0\n2 2 0\nend\n'+unchanged('dea_fuzzy, inputs(x1) outputs(y1) dmu(id) ort(output) alpha(1)',459))
    summary=dict(total=len(results),passed=sum(r['passed'] for r in results),lp_count=sum(r.get('lp_count',0) for r in results),max_oracle_gap=max((r.get('max_oracle_gap',0) for r in results),default=0))
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
    raise SystemExit(not all(r['passed'] for r in results))


if __name__=='__main__':main()
