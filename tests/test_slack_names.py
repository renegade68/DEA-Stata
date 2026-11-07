#!/usr/bin/env python3
"""Verify public TE slack column identities, order, and numerical preservation.

Uses Python's standard library and Stata 16+. Optional pre-R02 sources provide
bit-exact comparisons of every returned numeric matrix and scalar. Instrumented
copies only exercise auto fallback and missing-value metadata, not production.
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

from test_saveref import instrument, HOOK
from test_output_protection import RANDOM as ENV_DATA

REPO = Path(__file__).resolve().parents[1]
DATA = '''input id x1 x2 y1 y2
1 1 1 1 1
2 1 1 1 2
3 1 2 1 1
4 2 2 1 1
5 2 1 2 1
6 3 3 2 3
7 4 1 2 1
8 1 3 1 4
end
gen double x3=x1+2*x2
gen double y3=2*y1+y2
gen str30 note="original"
char note[audit] "preserve"
order id x1 x2 x3 y1 y2 y3
'''
LONG_NAMES = ['input_a_'+'x'*24, 'input_b_'+'x'*24,
              'output_a_'+'y'*23, 'output_b_'+'y'*23]
PROFILES = {
    'basic': ('x1 x2','y1 y2',''),
    'reversed': ('x2 x1','y2 y1',''),
    'asymmetric': ('x3 x1 x2','y2',''),
    'long': (' '.join(LONG_NAMES[:2]),' '.join(LONG_NAMES[2:]),
             ''.join(f'rename {a} {b}\n' for a,b in zip(['x1','x2','y1','y2'],LONG_NAMES))),
}
QUALIFIERS = [('', '1'),('if id!=2','id!=2'),('in 2/8','inrange(_n,2,8)'),
              ('if id!=2 in 1/7','id!=2 & inrange(_n,1,7)')]
EXPORTER = '''mata:
void r02_dump(string scalar path, real matrix A)
{
    real scalar f
    real matrix B
    B=bufio();bufbyteorder(B,2);f=fopen(path,"w")
    fbufput(B,f,"%4b",rows(A));fbufput(B,f,"%4b",cols(A))
    fbufput(B,f,"%8z",A);fclose(f)
}
end
'''


def names_check(x,y,old=False):
    operator = '!=' if old else '=='
    return f'''assert "`:colnames r(slack_in)'"{operator}"{x}"
assert "`:colnames r(slack_out)'"{operator}"{y}"
assert "`:rowfullnames r(slack_in)'"=="`:rowfullnames r(efficiency)'"
assert "`:rowfullnames r(slack_out)'"=="`:rowfullnames r(efficiency)'"
assert colsof(r(slack_in))=={len(x.split())}
assert colsof(r(slack_out))=={len(y.split())}
assert rowsof(r(slack_in))==r(N)
assert rowsof(r(slack_out))==r(N)
'''


def capture_numeric():
    return '''local mats : r(matrices)
foreach name of local mats {
    mata: r02_dump("matrix_`name'.bin",st_matrix("r(`name')"))
}
local scalars : r(scalars)
foreach name of local scalars {
    mata: r02_dump("scalar_`name'.bin",st_numscalar("r(`name')"))
}
'''


def case_body(x,y,engine,rts,ori,form,index,old=False,data=DATA,qual=None,
              extra='',entry='dea',missing=False,repeat=True):
    q, selected = QUALIFIERS[index%4] if qual is None else qual
    s = data+f'gen byte selected=({selected}) & !missing({",".join((x+" "+y).split())})\n'
    s += 'quietly ds\nlocal original_vars "`r(varlist)\'"\ntempfile original\nsave `original\'\n'
    if index%3==0:s+='gen str12 res="replace me"\n'
    refs = index%2==0
    io = f'{x}={y} {q},' if form=='A' else f'{y} {q}, inputs({x})'
    cmd = f'{entry} {io} dmu(id) engine({engine}) rts({rts}) orientation({ori}) stage(2) gen(res) replacevars {extra}'
    cmd += ' reference saveref' if refs else ''
    s += 'set seed 87312\n'+cmd+' saving("saved.dta")\n'
    s += names_check(x,y,old)
    s += capture_numeric()
    s += '''matrix SI=r(slack_in)
matrix SO=r(slack_out)
matrix E=r(efficiency)
local refs "`r(refvars)'"
assert "`r(inputs)'"=="'''+x+'''"
assert "`r(outputs)'"=="'''+y+'''"
'''
    if refs and not missing:
        # The columns must correspond to each stated variable, not just have
        # familiar names. Check the projection equations in the original units.
        s+='matrix L=r(lambda)\n'
        s+=f'''mata:
X=st_data(selectindex(st_data(.,"selected")),tokens("{x}"))
Y=st_data(selectindex(st_data(.,"selected")),tokens("{y}"))
L=st_matrix("L");SI=st_matrix("SI");SO=st_matrix("SO");E=st_matrix("E")
for (k=1;k<=rows(E);k++) {{
    if (!hasmissing((E[k,.],L[k,.],SI[k,.],SO[k,.]))) {{
        targetx={ 'E[k,1]*X[k,.]' if ori=='input' else 'X[k,.]' }
        targety={ 'Y[k,.]' if ori=='input' else 'Y[k,.]/E[k,1]' }
        assert(max(abs(targetx-L[k,.]*X-SI[k,.]):/(1:+abs(targetx)))<1e-6)
        assert(max(abs(L[k,.]*Y-targety-SO[k,.]):/(1:+abs(targety)))<1e-6)
    }}
}}
end
'''
    if missing:
        s+='mata: assert(all(st_matrix("SI"):==.a) & all(st_matrix("SO"):==.z))\n'
    s+='cf `original_vars\' using `original\',all\nassert "`:char note[audit]\'"=="preserve"\n'
    s+='''preserve
keep if selected
keep id '''+x+' '+y+''' res `refs'
tempfile expected
save `expected'
use "saved.dta",clear
'''
    # Bootstrap files contain additional inference columns; compare the common
    # result fields while retaining independent matrix comparisons above.
    s+='cf id '+x+' '+y+' res `refs\' using `expected\',all\nrestore\n'
    if repeat:
        s+='''matrix C=SI
matrix D=SO
local oldrows : rowfullnames SI
tempvar noise1 noise2 noise3
gen double `noise1'=123
gen double `noise2'=456
gen strL `noise3'="allocate a different temporary name"
'''
        # No persistent references on the repeated call: only the request for
        # displaying references remains, so existing saved reference vars stay.
        repeated=cmd.replace(' reference saveref',' reference')
        s+='set seed 87312\n'+repeated+'\n'+names_check(x,y,old)
        s+='mata: assert(all(st_matrix("C"):==st_matrix("r(slack_in)")) & all(st_matrix("D"):==st_matrix("r(slack_out)")))\n'
    return s


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--baseline-ado',type=Path);ap.add_argument('--only')
    a=ap.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    pkg=out/'current_ado';pkg.mkdir();hashes={}
    for name in package_files(REPO):
        raw=(REPO/name).read_bytes();(pkg/Path(name).name).write_bytes(raw);hashes[name]=hashlib.sha256(raw).hexdigest()
    retry=out/'retry_ado';instrument(pkg,retry)
    absent=out/'missing_slack_ado';shutil.copytree(pkg,absent)
    p=absent/'dea.ado';text=p.read_text();assert text.count(HOOK)==1
    inject='''    // TEST ONLY: exercise labels and preservation of extended missings.
    foreach v of local slack_in_vars {
        quietly replace `v'=.a if `touse'
    }
    foreach v of local slack_out_vars {
        quietly replace `v'=.z if `touse'
    }
'''
    p.write_text(text.replace(HOOK,inject+HOOK))
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,stata=a.stata,baseline=str(a.baseline_ado)),indent=2)+'\n')
    results=[]

    def run(name,body,source=pkg,compare=None):
        if a.only and not re.search(a.only,name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 160\n'
        for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:
            (d/key).mkdir();head+=f'sysdir set {key} "{d/key}"\n'
        head+=f'adopath ++ "{source.resolve()}"\n'+EXPORTER
        (d/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n')
        row=dict(test=name,passed=False);start=time.monotonic()
        try:
            proc=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=300)
            (d/'process.json').write_text(json.dumps(dict(returncode=proc.returncode,stderr=proc.stderr.decode(errors='replace')),indent=2)+'\n')
            log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else ''
            assert proc.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-1800:]
            if compare:
                before={p.name:p.read_bytes() for p in compare.glob('*.bin')}
                after={p.name:p.read_bytes() for p in d.glob('*.bin')}
                assert before and before==after, 'Historical numeric matrices/scalars differ'
                row['exact_numeric_files_compared']=len(before)
            row['passed']=True
        except Exception as e:row['error']=repr(e)
        row['seconds']=round(time.monotonic()-start,3);results.append(row)
        (out/'results.json').write_text(json.dumps(results,indent=2)+'\n')
        if not row['passed']:print(json.dumps(dict(test=name,error=row['error'][-500:])),flush=True)

    for i,(profile,engine,rts,ori,form) in enumerate(itertools.product(PROFILES,['simplex','ipm','auto'],['crs','vrs','nirs','ndrs'],['input','output'],['A','B'])):
        x,y,setup=PROFILES[profile];name=f'grid_{profile}_{engine}_{rts}_{ori}_{form}'
        kwargs=dict(x=x,y=y,engine=engine,rts=rts,ori=ori,form=form,index=i,data=DATA+setup)
        compare=None
        if a.baseline_ado:
            run(name+'_before',case_body(**kwargs,old=True),a.baseline_ado)
            compare=out/(name+'_before')
        run(name,case_body(**kwargs),compare=compare)
    for i,(engine,ori) in enumerate(itertools.product(['simplex','ipm','auto'],['input','output'])):
        for token,x,y,form in [('wildcard','x1 x2 x3','y1 y2 y3','A'),('range','x1 x2 x3','y1 y2 y3','B')]:
            # These expanded columns include exact linear combinations. Native
            # IPM's default tolerance gave a 1.039e-6 projection residual both
            # before and after R02. Request tighter solver precision here;
            # retain the independent projection assertion at 1e-6.
            body=case_body(x,y,engine,'vrs',ori,form,i,extra='tolerance(1e-9)')
            if token=='wildcard':body=body.replace('dea '+x+'='+y,'dea x*=y*')
            else:body=body.replace('dea '+y,'dea y1-y3').replace('inputs('+x+')','inputs(x1-x3)')
            run(f'expanded_{token}_{engine}_{ori}',body)
        run(f'alias_{engine}_{ori}',case_body('x2 x1','y2 y1',engine,'vrs',ori,'B',i,entry='dea_lp'))
        run(f'missing_{engine}_{ori}',case_body('x1 x2','y1 y2',engine,'crs',ori,'A',i,missing=True),absent)
    for i,(rts,ori) in enumerate(itertools.product(['crs','vrs','nirs','ndrs'],['input','output'])):
        s=case_body('x2 x1','y2 y1','auto',rts,ori,'B',i)
        s+='assert "`r(engine)\'"=="ipm"\n'
        run(f'retry_{rts}_{ori}',s,retry)
    for i,(engine,rts,ori) in enumerate(itertools.product(['simplex','ipm','auto'],['crs','vrs'],['input','output'])):
        run(f'bootstrap_{engine}_{rts}_{ori}',case_body('x1 x2','y1 y2',engine,rts,ori,'A',i,extra='bootstrap reps1(50)',repeat=False))
    env=ENV_DATA+'char note[audit] "preserve"\n'
    run('environmental_bootstrap',case_body('x','y','simplex','crs','input','B',0,
        data=env,qual=('if id<=75','id<=75'),extra='bootstrap envars(z) reps1(50) reps2(20)',repeat=False))
    for engine,form in itertools.product(['simplex','ipm','auto'],['A','B']):
        s=case_body('x1 x2','y1 y2',engine,'crs','input',form,1,repeat=False)
        s+=f'dea x1 x2=y1 y2, engine({engine}) stage(1) dmu(id) gen(stage1)\n'
        s+='capture matrix list r(slack_in)\nassert _rc==111\ncapture matrix list r(slack_out)\nassert _rc==111\n'
        s+='dea x1 x2=y1 y2, model(fdh) dmu(id) gen(fdh_result)\n'
        s+='capture matrix list r(slack_in)\nassert _rc==111\ncapture matrix list r(slack_out)\nassert _rc==111\n'
        run(f'no_stale_{engine}_{form}',s)
    summary=dict(total=len(results),passed=sum(r['passed'] for r in results),failed=sum(not r['passed'] for r in results),
                 historical_pairs=sum('exact_numeric_files_compared' in r for r in results))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
    raise SystemExit(summary['failed']>0)


if __name__=='__main__':main()
