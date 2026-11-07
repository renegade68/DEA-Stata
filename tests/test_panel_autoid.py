#!/usr/bin/env python3
"""Stable selected-position panel IDs, exact reference controls and old defects.

Python standard library and Stata 16+. Private trace copies capture the actual
ID/sample mapping. The expected ranks come from Python, before Stata runs.
"""
from _support import package_files
import argparse
from collections import Counter
import difflib
import hashlib
import itertools
import json
from pathlib import Path
import random
import re
import shutil
import subprocess
import time
from test_output_protection import unchanged

REPO=Path(__file__).resolve().parents[1]
SCALARS={'window':['n_periods','n_dmus','width','n_windows','mean_eff','sd_eff','n_obs'],
         'adjacent':['n_periods','n_dmus','mean_ec','mean_tc','mean_mpi'],
         'global':['n_periods','n_dmus','mean_ec','mean_tc','mean_mpi']}
KINDS=['adjacent','global','window']


def command(kind,qualifier='',form='A',prefix='res',dmu='',rts='crs',ori='input',dim=2,width=2,extra=''):
    x='x1'+(' x2' if dim==2 else '');y='y1'+(' y2' if dim==2 else '')
    io=f'{x} = {y} {qualifier},' if form=='A' else f'{qualifier}, inputs({x}) outputs({y})'
    model='window' if kind=='window' else 'malmquist'
    spec=f'width({width})' if kind=='window' else 'global' if kind=='global' else ''
    return f'dea_{model} {io} time(period) '+(f'dmu({dmu}) ' if dmu else '')+f'{spec} gen({prefix}) rts({rts}) ort({ori}) {extra}'


def outputs(kind,prefix='res',width=2,nperiods=3):
    if kind!='window':return [prefix+s for s in ['_ec','_tc','_mpi']]
    width=min(3,nperiods) if width==0 else width
    return [prefix+s for s in ['_mean','_sd']]+[f'{prefix}_w{w}_p{p}' for w in range(1,nperiods-width+2) for p in range(1,width+1)]


def fixture(mode,order='period',gaps=False):
    rows=[];periods=[2.5,4.5,9.5] if gaps else [1,2,3]
    for t,period in enumerate(periods):
        for ident in range(1,5):
            x=[1,2,4,5][ident-1];y=[1,3,4,6][ident-1]
            if ident==2:x*=1+.05*t
            if ident==3:x*=1-.04*t
            rows.append(dict(truth=ident,period=period,pick=0 if mode=='unbalanced' and t==1 and ident==2 else 1,x1=x,x2=x*(1+.1*ident),y1=y,y2=y*(1+.05*ident)))
        if mode in ['if','ifin','missing']:
            for j in range(t+1):
                row=dict(truth=0,period=period,pick=0,x1=99.,x2=100.,y1=99.,y2=100.)
                if mode=='missing':
                    row['pick']=1;row[['x1','y1','period'][j%3]]=None
                rows.insert(max(0,len(rows)-2*j),row)
    if order=='interleaved':rows.sort(key=lambda r:(r['truth'],r['period'] if r['period'] is not None else 99))
    if order=='reverse':rows.reverse()
    if order=='shuffle':random.Random(711).shuffle(rows)
    if mode in ['in','ifin']:
        rows=[dict(truth=0,period=periods[0],pick=0,x1=99.,x2=100.,y1=99.,y2=100.)]*2+rows+[dict(truth=0,period=periods[-1],pick=0,x1=99.,x2=100.,y1=99.,y2=100.)]
    # Independent positional oracle: walk original rows; never sort by Stata keys.
    count=Counter();table=[];q={'all':'','if':'if pick==1','in':f'in 3/{len(rows)-1}',
        'ifin':f'if pick==1 in 3/{len(rows)-1}','missing':'','unbalanced':'if pick==1'}[mode]
    for i,row in enumerate(rows,1):
        selected=all(row[v] is not None for v in ['period','x1','x2','y1','y2'])
        if mode in ['if','ifin','unbalanced']:selected=selected and row['pick']==1
        if mode in ['in','ifin']:selected=selected and 3<=i<=len(rows)-1
        expected=None
        if selected:count[row['period']]+=1;expected=count[row['period']]
        values=[i,row['truth'],row['period'],row['pick'],row['x1'],row['x2'],row['y1'],row['y2'],int(selected),expected]
        table.append(' '.join('.' if x is None else str(x) for x in values))
    setup='input double(obsno truth period pick x1 x2 y1 y2 selected expected_id)\n'+'\n'.join(table)+'\nend\ngen str20 note="preserve me"\nchar note[audit] "metadata"\nsort obsno\n'
    return setup,q,count


def instrument(source,dest):
    shutil.copytree(source,dest)
    for key in ['malmquist','window']:
        p=dest/f'dea_{key}.ado';old=p.read_text()
        anchor="    _dea_panel_check, dmu(`dmu') time(`time') touse(`touse')"
        assert old.count(anchor)==1
        insert='''    // TEST COPY ONLY: inspect the selected sample and actual generated IDs.
    if "`_dmu_auto'" != "" {
        mata: st_matrix("AUDIT_PANEL_MAP", (st_data(., "obsno"), st_data(., "`touse'"), st_data(., "`dmu'")))
    }
'''
        new=old.replace(anchor,insert+anchor);p.write_text(new)
        (dest.parent/(dest.name+'_'+key+'.patch')).write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True))))


def mapping_checks():
    return '''mata:
A=st_matrix("AUDIT_PANEL_MAP")
assert(all(A:==st_data(.,tokens("obsno selected expected_id"))))
end
assert "`r(dmu_method)'"=="position"
'''


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stata',required=True);p.add_argument('--output',required=True,type=Path);p.add_argument('--baseline-ado',type=Path);p.add_argument('--only');args=p.parse_args()
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'current_ado';pkg.mkdir();hashes={}
    for name in package_files(REPO):
        data=(REPO/name).read_bytes();(pkg/Path(name).name).write_bytes(data);hashes[name]=hashlib.sha256(data).hexdigest()
    trace=out/'trace_ado';instrument(pkg,trace)
    if args.baseline_ado:instrument(args.baseline_ado.resolve(),out/'trace_before')
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,baseline=str(args.baseline_ado),stata=args.stata),indent=2))
    results=[]
    def run(name,script,source=pkg):
        if args.only and not re.search(args.only,name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 120\n'
        for f in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/f).mkdir();head+=f'sysdir set {f} "{d/f}"\n'
        head+=f'adopath ++ "{source.resolve()}"\n';(d/'test.do').write_text(head+script+'\ndi "TEST_PASS"\n');row=dict(test=name,passed=False);start=time.monotonic()
        try:
            proc=subprocess.run([args.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=180)
            if proc.returncode or not (d/'test.log').exists():(d/'process.json').write_text(json.dumps(dict(returncode=proc.returncode,stderr=proc.stderr.decode(errors='replace')),indent=2))
            log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else proc.stderr.decode(errors='replace')
            assert proc.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2300:]
            row['passed']=True
        except Exception as e:row['error']=str(e)
        row['seconds']=round(time.monotonic()-start,3);results.append(row);(out/'results.json').write_text(json.dumps(results,indent=2));print(json.dumps(row),flush=True)

    # Both syntax forms/orientations/all RTS plus selected sample and source-order variants.
    for kind in KINDS:
        for rts,ori,mode,order in itertools.product(['crs','vrs','drs','irs'] if kind=='window' else ['crs','vrs'],['input','output'],['all','if','in','ifin','missing','unbalanced'],['period','interleaved','reverse','shuffle']):
            setup,q,count=fixture(mode,order,gaps=order=='reverse');nperiods=len(count);width=0 if mode=='all' else 2
            opts=dict(kind=kind,qualifier=q,rts=rts,ori=ori,dim=1 if order=='period' else 2,width=width)
            scalars=','.join('r('+v+')' for v in SCALARS[kind])
            s=setup+'''quietly ds
local original_vars "`r(varlist)'"
quietly _datasignature `original_vars'
local before "`r(datasignature)'"
local sort_before : sortedby
'''
            s+=command(**opts,form='A',prefix='ref',dmu='expected_id')+'\nassert "`r(dmu_method)\'"=="variable"\nmatrix Reference=('+scalars+')\n'
            for seed in [11,732]:
                s+=f'set sortseed {seed}\n'+command(**opts,form='B',extra='replacevars detail debug')+'\nmatrix Got=('+scalars+')\n'+mapping_checks()
                s+='mata: assert(max(abs(st_matrix("Got")-st_matrix("Reference")))<1e-12)\n'
                for a,b in zip(outputs(kind,width=width,nperiods=nperiods),outputs(kind,'ref',width=width,nperiods=nperiods)):
                    s+=f'assert {a}=={b}\nassert missing({a}) if !selected\n'
                expected_dm=max(count.values()) if kind=='window' else min(count[t] for t in sorted(count)[-2:])
                s+=f'assert r(n_dmus)=={expected_dm}\n'
                s+='quietly _datasignature `original_vars\'\nassert "`before\'"=="`r(datasignature)\'"\nassert "`:char note[audit]\'"=="metadata"\n'
            s+='assert "`:sortedby\'"=="`sort_before\'"\ncapture ds __*\nassert _rc==111\n'
            run(f'grid_{kind}_{rts}_{ori}_{mode}_{order}',s,trace)

    # Exact observed defect: automatic IDs previously counted excluded rows.
    defect='''input id period x1 y1
1 1 1 1
2 1 2 3
3 1 4 4
0 2 99 99
1 2 1 1.2
2 2 2 3.6
3 2 4 4.8
end
gen double x2=x1
gen double y2=y1
gen obsno=_n
'''
    for kind in KINDS:
        old=defect+command(kind,'if id>0')+'\nassert r(n_dmus)=='+('4' if kind=='window' else '2')+'\n'
        if args.baseline_ado:run('original_defect_'+kind,old,args.baseline_ado)
        s=defect+command(kind,'if id>0',dmu='id',prefix='ref')+'\nassert r(n_dmus)==3\n'+command(kind,'if id>0')+'\nassert r(n_dmus)==3\nassert "`r(dmu_method)\'"=="position"\n'
        for a,b in zip(outputs(kind,nperiods=2),outputs(kind,'ref',nperiods=2)):s+=f'assert {a}=={b}\n'
        run('original_corrected_'+kind,s)

    # Real-world identities can differ from selected ranks; this remains an explicit assumption.
    for kind in KINDS:
        setup,q,count=fixture('unbalanced')
        s=setup+command(kind,q,dmu='truth',prefix='truthfit')+'\n'+command(kind,q)+'\nassert "`r(dmu_method)\'"=="position"\n'
        # No assertion of identity equivalence: actual map must match the independent positional oracle.
        s+=mapping_checks()
        run('unbalanced_identity_assumption_'+kind,s,trace)
        for qname,qualifier in [('empty','if 0'),('one_row','in 1'),('one_period','if period==1')]:
            s=setup
            for v in outputs(kind):s+=f'gen str12 {v}="old result"\n'
            s+=unchanged(command(kind,qualifier,extra='replacevars'),198)
            run(f'atomic_{kind}_{qname}',s)
        # A qualifier using the old result is evaluated before replacement.
        setup,q,count=fixture('if')
        s=setup+'gen double res_mean=pick\n'
        s+=command(kind,'if res_mean==1',dmu='expected_id',prefix='ref')+'\n'+command(kind,'if res_mean==1',extra='replacevars')+'\n'+mapping_checks()
        for a,b in zip(outputs(kind),outputs(kind,'ref')):s+=f'assert {a}=={b}\n'
        run('replace_sample_'+kind,s,trace)
        # Explicit IDs preserve meaning under reordering and expose their method metadata.
        for label,expr in [('numeric','truth'),('string','"DMU "+string(truth)')]:
            setup,q,count=fixture('all')
            s=setup+('gen double key=' if label=='numeric' else 'gen str20 key=')+expr+'\n'+command(kind,dmu='key')+'\nassert "`r(dmu_method)\'"=="variable"\nrename res_* old_*\ngsort -obsno\n'+command(kind,dmu='key')+'\nassert "`r(dmu_method)\'"=="variable"\n'
            for a,b in zip(outputs(kind),outputs(kind,'old')):s+=f'assert abs({a}-{b})<1e-6 if !missing({b})\nassert missing({a})==missing({b})\n'
            run(f'explicit_{kind}_{label}',s)
    summary=dict(total=len(results),passed=sum(r['passed'] for r in results))
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
    raise SystemExit(summary['total']!=summary['passed'])

if __name__=='__main__':main()
