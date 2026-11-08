#!/usr/bin/env python3
"""Check strategy-specific minimum samples without changing crossgame scores.

Each case uses an isolated Stata session. Optional pre-R04 sources verify exact
numeric preservation for singleton Nash and multi-DMU results. A private counter
checks whether the computation entry point was reached; no engine is modified.
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

from test_sample_selection import DATA
from test_output_protection import unchanged
from test_slack_names import EXPORTER

REPO=Path(__file__).resolve().parents[1]
STRATEGIES=['aggressive','benevolent','nash','all','default']
SINGLETONS={
    'whole':('keep in 4\n','','1'),
    'if':('','if id==4','id==4'),
    'in':('','in 4/4','_n==4'),
    'combined':('','if group==1 in 3/5','id==4'),
    'missing_input':('replace x1=. if id!=4\n','','id==4'),
    'missing_output':('replace y2=.a if id!=4\n','','id==4'),
    'combined_missing':('replace x2=.z if id!=4\n','if id<=6 in 2/8','id==4'),
}
EMPTY={
    'whole':('keep if 0\n',''),
    'if':('','if id<0'),
    'combined':('','if group==0 in 4/4'),
    'missing':('replace x1=.a\n',''),
}
OLD='''gen str12 res_agg="old aggressive"
gen double res_ben=991
gen float res_nash=992
gen str12 res_keep="untouched"
'''

def outputs(strategy):
    return {'aggressive':['res_agg'],'benevolent':['res_ben'],'nash':['res_nash'],
            'all':['res_agg','res_ben','res_nash'],'default':['res_agg','res_ben','res_nash']}[strategy]

def cmd(strategy,form='A',q='',dmu='id',replace=False,inputs='x1 x2',rhs='y1 y2',extra=''):
    io=f'{inputs}={rhs} {q},' if form=='A' else f'{q}, inputs({inputs}) outputs({rhs})'
    return f'dea_crossgame {io} gen(res)'+(' strategy('+strategy+')' if strategy!='default' else '')+(' dmu('+dmu+')' if dmu else '')+(' replacevars' if replace else '')+' '+extra

def successful(setup,selection,command,strategy,n,analytic=False):
    outs=' '.join(outputs(strategy))
    s=setup+f'gen byte selected=({selection})\n'
    s+='quietly ds\nlocal originals "`r(varlist)\'"\ntempfile original\nsave `original\'\n'
    s+=command+f'\nassert r(n)=={n}\nassert "`r(strategy)\'"=="{strategy if strategy!="default" else "all"}"\n'
    s+='mata: r02_dump("return_n.bin",st_numscalar("r(n)"))\n'
    s+='mata: r02_dump("scores.bin",st_data(.,tokens("'+outs+'")))\n'
    if analytic:
        s+='gen double ratio=y1/x1 if selected\nquietly summarize ratio,meanonly\nlocal frontier=r(max)\n'
        s+='foreach v of varlist '+outs+' {\nassert !missing(`v\') if selected\nassert abs(`v\'-ratio/`frontier\')<1e-6 if selected\n}\n'
    elif n==1:
        s+='assert !missing(res_nash) & abs(res_nash-1)<1e-6 if selected\n'
    else:
        s+='foreach v of varlist '+outs+' {\ncount if selected & !missing(`v\')\nassert r(N)>0\n}\n'
    s+='foreach v of varlist '+outs+' {\nassert missing(`v\') if !selected\n}\n'
    s+='local replaced "'+outs+'"\nlocal originals : list originals - replaced\ncf `originals\' using `original\',all\n'
    return s

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--baseline-ado',type=Path);ap.add_argument('--only')
    a=ap.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    pkg=out/'ado';pkg.mkdir();hashes={}
    for name in package_files(REPO):
        raw=(REPO/name).read_bytes();(pkg/Path(name).name).write_bytes(raw);hashes[name]=hashlib.sha256(raw).hexdigest()
    audit=out/'counter_ado';shutil.copytree(pkg,audit)
    p=audit/'_dea_crossgame_core.ado';code=p.read_text();needle='    X = st_data(., tokens(inputs), touse)'
    assert code.count(needle)==1
    p.write_text(code.replace(needle,'    st_numscalar("CG_CALLS",st_numscalar("CG_CALLS")+1) // TEST ONLY\n'+needle))
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,baseline=str(a.baseline_ado),stata=a.stata),indent=2)+'\n')
    rows=[]
    def run(name,body,source=pkg,compare=None):
        if a.only and not re.search(a.only,name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 160\n'
        for sub in ['PLUS','PERSONAL','SITE','OLDPLACE']:
            (d/sub).mkdir();head+=f'sysdir set {sub} "{d/sub}"\n'
        head+=f'adopath ++ "{source.resolve()}"\n'+EXPORTER
        (d/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n')
        start=time.monotonic();row=dict(test=name,passed=False)
        try:
            proc=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=120)
            log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else proc.stderr.decode(errors='replace')
            assert proc.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-1800:]
            if compare:
                old={p.name:p.read_bytes() for p in compare.glob('*.bin')}
                new={p.name:p.read_bytes() for p in d.glob('*.bin')}
                assert old and old==new,'Historical numeric results differ'
                row['exact_numeric_files']=len(old)
            row['passed']=True
        except Exception as e:row['error']=repr(e)
        row['seconds']=round(time.monotonic()-start,3);rows.append(row)
        (out/'results.json').write_text(json.dumps(rows,indent=2)+'\n')
        if not row['passed']:print(json.dumps(row),flush=True)

    # Cardinality is counted AFTER qualifiers and missing-value exclusions.
    for i,(strategy,form,dmu,kind,replace) in enumerate(itertools.product(STRATEGIES,['A','B'],['id','sid',''],SINGLETONS,[False,True])):
        setup,q,selection=SINGLETONS[kind]
        data=DATA+setup+(OLD if replace else 'gen str12 res_keep="untouched"\n')
        extra=['','matrix','detail','debug','matrix detail debug'][i%5]
        command=cmd(strategy,form,q,dmu,replace,extra=extra)
        name=f'one_{strategy}_{form}_{dmu or "auto"}_{kind}_{int(replace)}'
        if strategy=='nash':
            body=successful(data,selection,command,strategy,1)
            prior=None
            if a.baseline_ado:
                run(name+'_before',body,a.baseline_ado);prior=out/(name+'_before')
            run(name,body,compare=prior)
        else:
            run(name,data+unchanged(command,2001))

    for strategy,form,kind,replace in itertools.product(STRATEGIES,['A','B'],EMPTY,[False,True]):
        setup,q=EMPTY[kind]
        run(f'empty_{strategy}_{form}_{kind}_{int(replace)}',DATA+setup+(OLD if replace else '')+unchanged(cmd(strategy,form,q,'',replace),2000))

    for i,(profile,strategy,form) in enumerate(itertools.product(['two','duplicate','multi'],STRATEGIES,['A','B'])):
        dmu=['id','sid',''][i%3]
        setup=DATA
        if profile=='duplicate':
            setup+='replace x1=2 in 4\nreplace x1=2 in 8\nreplace y1=1 in 4\nreplace y1=1 in 8\nreplace id=4 in 8\nreplace sid="same label" in 4\nreplace sid="same label" in 8\n'
        q='if group' if profile=='multi' else 'if _n==4 | _n==8'
        selection='group' if profile=='multi' else '_n==4 | _n==8'
        inputs='x1 x2' if profile=='multi' else 'x1';rhs='y1 y2' if profile=='multi' else 'y1'
        command=cmd(strategy,form,q,dmu,inputs=inputs,rhs=rhs,extra='matrix detail maxiter(30)')
        body=successful(setup,selection,command,strategy,4 if profile=='multi' else 2,analytic=profile!='multi')
        name=f'normal_{profile}_{strategy}_{form}';prior=None
        if a.baseline_ado:
            run(name+'_before',body,a.baseline_ado);prior=out/(name+'_before')
        run(name,body,compare=prior)

    # Instrument only private copies. The invalid sample must never reach Mata.
    for n,strategy,form in itertools.product([0,1,2],STRATEGIES,['A','B']):
        q=['if 0','in 4/4','if _n==4 | _n==8'][n]
        rc=2000 if n==0 else (2001 if n==1 and strategy!='nash' else 0)
        s=DATA+'scalar CG_CALLS=0\n'
        s+='capture noisily '+cmd(strategy,form,q,'',extra='matrix detail debug')+f'\nassert _rc=={rc}\nassert CG_CALLS=={int(rc==0)}\n'
        run(f'counter_{n}_{strategy}_{form}',s,audit)

    for strategy in STRATEGIES:
        # A qualifier may refer to a result that replacevars would otherwise drop.
        setup=DATA+'gen double res_agg=(id==4)\ngen str8 res_ben="previous"\ngen str8 res_nash="previous"\n'
        command=cmd(strategy,q='if res_agg==1',replace=True)
        if strategy=='nash':body=successful(setup,'id==4',command,strategy,1)
        else:body=setup+unchanged(command,2001)
        run('result_qualifier_'+strategy,body)
    run('name_collision_precedes_count',DATA+'keep in 4\n'+OLD+unchanged(cmd('aggressive'),110))
    run('model_protected_precedes_count',DATA+'keep in 4\nrename x1 res_agg\n'+unchanged(cmd('all',replace=True,inputs='res_agg x2'),110))
    summary=dict(total=len(rows),passed=sum(r['passed'] for r in rows),failed=[r['test'] for r in rows if not r['passed']],
                 historical_pairs=sum('exact_numeric_files' in r for r in rows))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
    raise SystemExit(bool(summary['failed']))

if __name__=='__main__':main()
