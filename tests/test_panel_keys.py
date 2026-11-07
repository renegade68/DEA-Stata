#!/usr/bin/env python3
"""Panel-key preflight: duplicates, selected samples, compatibility and LP scores.

Requires NumPy/SciPy and Stata 16+. Historical comparisons use a pre-F13 ado
snapshot. All generated sources/data stay in a new --output directory.
"""
from _support import package_files
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import re
import subprocess
import time

import numpy as np
import scipy
from test_output_protection import unchanged
from verify_bootstrap_lp import solve

REPO=Path(__file__).resolve().parents[1]
DATA='''input id x y
1 1 1
2 2 3
3 4 4
4 5 6
end
expand 3
bysort id: gen period=_n
replace x=x*(1+.05*(period-1)) if id==2
replace x=x*(1-.04*(period-1)) if id==3
gen double x2=x*(1+.1*id)
gen double y2=y*(1+.05*id)
sort period id
gen long rowid=_n
gen str20 note="preserve me"
'''
IDS={
 'numeric':'gen double key=id\n',
 'fractional':'gen double key=id/7\n',
 'large':'gen double key=1e15+id\n',
 'string':'gen str20 key="unit "+string(id)\n',
 'unicode':'gen str30 key="기관 "+string(id)\n',
 'long':'gen strL key="'+('a'*2500)+'"+string(id)\n',
 'space':'gen str20 key=substr("    ",1,id)+"same"\n',
 'quotes':'gen str30 key=char(34)+"group, "+string(id)+char(34)\n',
}
KINDS=['adjacent','global','window']

def cmd(kind,form='A',qualifier='',prefix='res',rts='crs',ori='input',dim=1,extra='',dmu='key',width=2):
    x='x'+(' x2' if dim==2 else '');y='y'+(' y2' if dim==2 else '')
    io=f'{x} = {y} {qualifier},' if form=='A' else f'{qualifier}, inputs({x}) outputs({y})'
    name='window' if kind=='window' else 'malmquist'
    opts=f'width({width})' if kind=='window' else 'global' if kind=='global' else ''
    return f'dea_{name} {io} time(period) '+(f'dmu({dmu}) ' if dmu else '')+f'{opts} gen({prefix}) rts({rts}) ort({ori}) {extra}'


def outputs(kind,prefix='res',width=2,nperiods=3):
    if kind!='window':return [prefix+s for s in ['_ec','_tc','_mpi']]
    return [prefix+s for s in ['_mean','_sd']]+[f'{prefix}_w{w}_p{p}' for w in range(1,nperiods-width+2) for p in range(1,width+1)]


def read(path):
    with path.open() as f:rows=list(csv.DictReader(f))
    return {k:np.array([float(row[k]) if row[k] not in ('','.') else np.nan for row in rows]) for k in rows[0]}


def oracle(d,kind,rts,ori,dim,width):
    data=read(d/'results.csv');selected=data['selected']==1;times=data['period'];ids=data['id']
    x=np.column_stack([data[v] for v in (['x','x2'] if dim==2 else ['x'])]);y=np.column_stack([data[v] for v in (['y','y2'] if dim==2 else ['y'])])
    periods=np.unique(times[selected]);expected={v:np.full(len(ids),np.nan) for v in outputs(kind,width=width,nperiods=len(periods))};count=0
    def dist(k,mask):
        nonlocal count
        count+=1;return solve(x[k],y[k],x[mask],y[mask],ori,rts)
    if kind=='window':
        for w in range(len(periods)-width+1):
            ref=selected&(times>=periods[w])&(times<=periods[w+width-1])
            for k in np.flatnonzero(ref):
                pos=int(np.flatnonzero(periods==times[k])[0])-w+1
                expected[f'res_w{w+1}_p{pos}'][k]=dist(k,ref)
        for k in np.flatnonzero(selected):
            values=[v[k] for name,v in expected.items() if '_w' in name and np.isfinite(v[k])]
            expected['res_mean'][k]=np.mean(values)
            expected['res_sd'][k]=np.std(values,ddof=1) if len(values)>1 else 0
    else:
        for t0,t1 in zip(periods,periods[1:]):
            m0=selected&(times==t0);m1=selected&(times==t1)
            for ident in np.intersect1d(ids[m0],ids[m1]):
                k0=np.flatnonzero(m0&(ids==ident))[0];k1=np.flatnonzero(m1&(ids==ident))[0]
                own0=dist(k0,m0);own1=dist(k1,m1);ec=own1/own0
                if kind=='global':mpi=dist(k1,selected)/dist(k0,selected);tc=mpi/ec
                else:tc=np.sqrt(dist(k1,m0)/own1*own0/dist(k0,m1));mpi=ec*tc
                for suffix,value in [('ec',ec),('tc',tc),('mpi',mpi)]:expected['res_'+suffix][k1]=value
    gap=0.;cells=0
    for name,values in expected.items():
        actual=data[name];assert np.array_equal(np.isnan(actual),np.isnan(values)),name
        valid=np.isfinite(values);cells+=sum(valid)
        delta=float(np.max(abs(actual[valid]-values[valid]))) if valid.any() else 0.
        gap=max(gap,delta);assert delta<1e-6,(name,delta)
    return dict(lp_count=count,oracle_cells=int(cells),max_oracle_gap=gap)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stata',required=True);p.add_argument('--output',required=True,type=Path);p.add_argument('--baseline-ado',type=Path);p.add_argument('--only');args=p.parse_args()
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'current_ado';pkg.mkdir();hashes={}
    for name in package_files(REPO):
        data=(REPO/name).read_bytes();(pkg/Path(name).name).write_bytes(data);hashes[name]=hashlib.sha256(data).hexdigest()
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,numpy=np.__version__,scipy=scipy.__version__,baseline=str(args.baseline_ado)),indent=2))
    rows=[]
    def run(name,script,source=pkg,check=None):
        if args.only and not re.search(args.only,name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 120\n'
        for f in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/f).mkdir();head+=f'sysdir set {f} "{d/f}"\n'
        head+=f'adopath ++ "{source.resolve()}"\n';(d/'test.do').write_text(head+script+'\ndi "TEST_PASS"\n')
        row=dict(test=name,passed=False);start=time.monotonic()
        try:
            proc=subprocess.run([args.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=180)
            if proc.returncode or not (d/'test.log').exists():(d/'process.json').write_text(json.dumps(dict(returncode=proc.returncode,stderr=proc.stderr.decode(errors='replace')),indent=2))
            log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else proc.stderr.decode(errors='replace')
            assert proc.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2100:]
            if check:row.update(check(d))
            row['passed']=True
        except Exception as e:row['error']=str(e)
        row['seconds']=round(time.monotonic()-start,3);rows.append(row);(out/'results.json').write_text(json.dumps(rows,indent=2));print(json.dumps(row),flush=True)

    for kind,identity,form,period,replace in itertools.product(KINDS,IDS,['A','B'],[1,3],[False,True]):
        data=DATA+IDS[identity]+f'expand 2 if id==2 & period=={period}, gen(added)\n'
        if period==3:data+='replace x=x*2 if added\n'
        data+='gsort -rowid\n'
        if replace:
            for v in outputs(kind):data+=f'gen str8 {v}="old"\n'
        run(f'duplicate_{kind}_{identity}_{form}_{period}_{replace}',data+unchanged(cmd(kind,form,extra='replacevars' if replace else ''),459))
    for kind in KINDS:
        for reverse in [False,True]:
            data=DATA+IDS['numeric']+'expand 3 if id==2 | id==3\n'+('gsort -rowid\n' if reverse else '')
            run(f'multiple_groups_{kind}_{reverse}',data+unchanged(cmd(kind),459))
        run(f'string_time_{kind}',DATA+IDS['numeric']+'tostring period, replace\n'+unchanged(cmd(kind,extra='replacevars'),7))
        # Missing keys are excluded before uniqueness validation, for numeric and string IDs.
        for identity in ['numeric','string']:
            setup=DATA+IDS[identity]+'expand 2 if id==2, gen(added)\n'
            for reason in ['if','in','ifin','x','y','id','time']:
                data=setup
                qualifier='';selected='added==0'
                if reason=='if':qualifier='if added==0'
                if reason=='in':data+='sort added rowid\n';qualifier='in 1/12'
                if reason=='ifin':data+='sort added rowid\n';qualifier='if added==0 in 1/15'
                if reason in ['x','y']:data+=f'replace {reason}=. if added\n'
                if reason=='id':data+='replace key='+('""' if identity=='string' else '.a')+' if added\n'
                if reason=='time':data+='replace period=.z if added\n'
                data+='gen byte selected=added==0\n'
                data+='preserve\nkeep if selected\n'+cmd(kind,prefix='control')+'\nkeep rowid control_*\ntempfile reference\nsave `reference\'\nrestore\n'
                data+=cmd(kind,qualifier=qualifier)+'\n'
                data+='preserve\nkeep if selected\n'
                data+='merge 1:1 rowid using `reference\',assert(match) nogen\n'
                for v,c in zip(outputs(kind),outputs(kind,'control')):data+=f'assert {v}=={c}\n'
                data+='restore\n'
                for v in outputs(kind):data+=f'assert missing({v}) if !selected\n'
                run(f'excluded_{kind}_{identity}_{reason}',data)

    # Normal results must retain exact historical values, including unbalanced panels.
    for kind in KINDS:
        for rts,ori,dim,shape in itertools.product(['crs','vrs','drs','irs'] if kind=='window' else ['crs','vrs'],['input','output'],[1,2],['balanced','unbalanced','gaps']):
            name=f'normal_{kind}_{rts}_{ori}_{dim}_{shape}';data=DATA+IDS['string'];qualifier=''
            if shape=='unbalanced':qualifier='if !(id==2 & period==2)'
            if shape=='gaps':data+='replace period=period*2+.25\n'
            width=3 if shape=='gaps' else 2
            data+='gen byte selected=1\n'+('replace selected=0 if id==2 & period==2\n' if shape=='unbalanced' else '')
            command=cmd(kind,'B' if dim==2 else 'A',qualifier,rts=rts,ori=ori,dim=dim,width=width,extra='detail debug')
            if args.baseline_ado:run(name+'_before',data+command+'\nkeep rowid res_*\nsave "expected.dta",replace\n',args.baseline_ado)
            script=data+command+'\nassert note=="preserve me"\n'
            if args.baseline_ado:script+=f'cf res_* using "{out/(name+"_before")/"expected.dta"}",all\n'
            script+='export delimited id period selected x x2 y y2 res_* using "results.csv",replace\n'
            run(name,script,check=lambda d,kind=kind,rts=rts,ori=ori,dim=dim,width=width:oracle(d,kind,rts,ori,dim,width))

    if args.baseline_ado:
        # Both old commands accepted duplicate keys; saved per-row summaries were incomplete.
        for kind in KINDS:
            for orientation in ['input','output']:
                data=DATA+IDS['numeric']+'expand 2 if id==2 & period==3,gen(added)\n'+cmd(kind,ori=orientation)+'\n'
                if kind=='window':data+='count if id==2 & period==3 & missing(res_mean)\nassert r(N)==1\n'
                else:data+='count if id==2 & period==3 & missing(res_mpi)\nassert r(N)==1\n'
                run(f'old_duplicate_{kind}_{orientation}',data,args.baseline_ado)
    # Exact string and double identity, with row order deliberately changed.
    for kind,identity in itertools.product(KINDS,IDS):
        data=DATA+IDS[identity]+cmd(kind)+'\nrename res_* before_*\ngsort -rowid\n'+cmd(kind)+'\n'
        for v,c in zip(outputs(kind),outputs(kind,'before')):data+=f'assert abs({v}-{c})<1e-6 if !missing({c})\nassert missing({v})==missing({c})\n'
        run(f'identity_order_{kind}_{identity}',data)
    # The checker alone supports strL labels beyond the widths used by displays.
    run('helper_exact_long_strings',DATA+IDS['long']+'''gen byte selected=1
_dea_panel_check, dmu(key) time(period) touse(selected)
assert r(duplicates)==0
expand 2 if id==2 & period==3
'''+unchanged('_dea_panel_check, dmu(key) time(period) touse(selected)',459))
    summary=dict(total=len(rows),passed=sum(r['passed'] for r in rows),lp_count=sum(r.get('lp_count',0) for r in rows),oracle_cells=sum(r.get('oracle_cells',0) for r in rows),max_oracle_gap=max((r.get('max_oracle_gap',0) for r in rows),default=0))
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
    raise SystemExit(summary['total']!=summary['passed'])

if __name__=='__main__':main()
