#!/usr/bin/env python3
"""F30 ordinal transformation and observation-order invariance.

Standard-library runner, licensed Stata 16+. Exports from isolated instrumented
copies are compared to an independent Python transform oracle. No production
solver or installed ado is instrumented. A four-DMU analytic oracle covers all
24 permutations and both orientations/four RTS restrictions/three engines.
"""
from _support import package_files
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import random
import re
import shutil
import subprocess
import time
from test_fdh_idea_returns import checks

REPO=Path(__file__).resolve().parents[1]
PROFILES={
 'single':('x','y qy','','qy'),
 'multiple':('x','y qy qz','','qy qz'),
 'both':('x qx','y qy','','qx qy'),
 'mixed':('x qx xb','y qy qz yb','xb yb','qx qy qz'),
 'ordinal_only':('x','qy','','qy'),
 'input_control':('x qx','y','','qx'),
}
BASE=[dict(id=i,x=x,y=y,qx=qx,qy=qy,qz=qz,xblo=x-.5,xbhi=x+.5,yblo=y-.5,ybhi=y+.5)
      for i,x,y,qx,qy,qz in [(1,5,1,2,1,3),(2,4,3,2,3,1),(3,2,2,1,2,2),(4,1,1,4,4,1),(5,6,4,3,3,2),(6,3,2,1,1,3)]]
SMALL='input double(id x y q)\n1 5 1 1\n2 4 1 3\n3 2 1 2\n4 1 1 4\nend\n'


def dataset(rows):
    names=['id','x','y','qx','qy','qz','xblo','xbhi','yblo','ybhi']
    s='input double('+ ' '.join(names)+')\n'
    s+='\n'.join(' '.join(format(row[k],'.17g') for k in names) for row in rows)+'\nend\n'
    s+='gen strL xb="["+strtrim(string(xblo,"%24.17g"))+","+strtrim(string(xbhi,"%24.17g"))+"]"\n'
    s+='gen strL yb="["+strtrim(string(yblo,"%24.17g"))+","+strtrim(string(ybhi,"%24.17g"))+"]"\n'
    return s


def order_code(order):
    s='gen long visit=.\n'
    for k,i in enumerate(order,1):s+=f'replace visit={k} if id=={i}\n'
    return s+'sort visit\n'


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--baseline-ado',type=Path);ap.add_argument('--only')
    args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    pkg=out/'current_ado';pkg.mkdir();hashes={}
    for n in package_files(REPO):
        data=(REPO/n).read_bytes();(pkg/Path(n).name).write_bytes(data);hashes[n]=hashlib.sha256(data).hexdigest()
    trace=out/'trace_ado';shutil.copytree(pkg,trace);p=trace/'dea.ado';s=p.read_text()
    anchor="            foreach _drop_output in `__idea_score' {"
    assert s.count(anchor)==1
    export='''            quietly format id `inputs' `outputs' %24.17g
            quietly export delimited id `inputs' `outputs' using "transform_${F30_capture_tag}_`dmui'.csv", replace datafmt

'''
    p.write_text(s.replace(anchor,export+anchor,1))
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,stata=args.stata,baseline=str(args.baseline_ado),trace='Export transformed data just before each ordinary DEA call; production unchanged'),indent=2))
    results=[]
    def chosen(name):return not args.only or re.search(args.only,name)
    def run(name,code,source=pkg,audit=None,force=False):
        if not force and not chosen(name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 120\n'
        for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:
            (d/key).mkdir();head+=f'sysdir set {key} "{d/key}"\n'
        head+=f'adopath ++ "{source.resolve()}"\n'
        (d/'test.do').write_text(head+code+'\ndi "TEST_PASS"\n')
        start=time.monotonic();row=dict(test=name,passed=False)
        try:
            proc=subprocess.run([args.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=300)
            log=(d/'test.log').read_text(errors='replace')
            if proc.returncode or 'TEST_PASS' not in log.splitlines():raise RuntimeError(log[-2400:])
            if audit:row.update(audit(d))
            row['passed']=True
        except Exception as e:row['error']=repr(e)
        row['seconds']=round(time.monotonic()-start,3);results.append(row)
        (out/'results.json').write_text(json.dumps(results,indent=2));print(json.dumps(row),flush=True)

    # Original defect: known scores differ by 0.8 only due to row order.
    for engine,ori in itertools.product(['simplex','ipm','auto'],['input','output']):
        name=f'original_{engine}_{ori}'
        if not chosen(name):continue
        setup=SMALL+f'dea x=y q, ordinal(q) engine({engine}) orientation({ori}) dmu(id) gen(first)\ngsort -id\n'
        setup+=f'dea x=y q, yord(q) engine({engine}) orientation({ori}) dmu(id) gen(second)\nsort id\n'
        if args.baseline_ado:run('old_'+name,setup+'gen double gap=abs(first-second)\nsummarize gap,meanonly\nassert abs(r(max)-.8)<1e-7\n',args.baseline_ado,force=True)
        run(name,setup+'assert abs(first-second)<1e-7\nassert abs(first-cond(id==2,.5,1))<1e-7\n')

    # Six profiles, four RTS, two orientations, three engines, both stages.
    # Each case retains the same effective sample under identity/reverse/random
    # permutations; in-based samples are permuted within their position range.
    for idx,((kind,(inputs,outputs,bounded,ordinal)),rts,ori,engine,stage) in enumerate(itertools.product(PROFILES.items(),['crs','vrs','nirs','ndrs'],['input','output'],['simplex','ipm','auto'],[1,2])):
        name=f'grid_{kind}_{rts}_{ori}_{engine}_s{stage}'
        if not chosen(name):continue
        rows=[dict(r) for r in BASE]
        if idx%3==1:
            for row in rows:
                row['qy']=1e12+row['qy']*(2**-10)
                row['qx']=-1e12+row['qx']*(2**-10)
        elif idx%3==2:
            for row in rows:row['qy']=1e20+row['qy']*16384;row['qz']=-row['qz']/8
        if idx%11==0:
            for row in rows:row['qy']=2 # All ties, with multiple output qz left variable.
        mode=idx%5
        ids=list(range(1,7));sample={0:ids,1:[1,2,4,5,6],2:[2,3,4,5],3:[2,4,5],4:[1,2,3,4,5]}[mode]
        qualifier={0:'',1:'if id!=3',2:'in 2/5',3:'if id!=3 in 2/5',4:''}[mode]
        mutable=[2,3,4,5] if mode in [2,3] else ids
        reverse=ids[:];random_order=ids[:];shuffled=mutable[:];random.Random(20261101+idx).shuffle(shuffled)
        for at,new,random_id in zip([x-1 for x in mutable],reversed(mutable),shuffled):reverse[at]=new;random_order[at]=random_id
        orders=[ids,reverse,random_order]
        setup=dataset(rows)
        if mode==4:setup+='replace x=.a if id==6\n'
        setup+='gen byte selected='+(' | '.join('id=='+str(i) for i in sample))+'\n'
        setup+='tempfile fixture original expected_order\nsave `fixture\'\n'
        ordins=' '.join(v for v in ordinal.split() if v in inputs.split());ordouts=' '.join(v for v in ordinal.split() if v in outputs.split())
        opts=f'ordinal({ordinal})' if idx%3==0 else (f'xord({ordins}) yord({ordouts})' if idx%3==1 else f'ordinal({ordinal}) xord({ordins}) yord({ordouts})')
        # Same variable repeated across aliases must still transform only once.
        opts+=f' bounded({bounded})'
        requested='imprecise' if idx%3==2 else 'te';model=['','model(te)','model(imprecise)'][idx%3]
        for j,order in enumerate(orders):
            setup+='use `fixture\',clear\n'+order_code(order)
            setup+='quietly ds\nlocal original_vars "`r(varlist)\'"\nsave `original\',replace\n'
            setup+='gen str8 res="old"\ngen str8 res_opt="old"\ngen str8 res_mid="old"\n'
            setup+=f'global F30_capture_tag p{j}\nlocal savevars id {inputs} {outputs}\n'
            prefix=f'dea {inputs}={outputs} {qualifier},' if (idx//2)%2 else f'dea {outputs} {qualifier}, inputs({inputs})'
            setup+=f'{prefix} {opts} {model} rts({rts}) orientation({ori}) engine({engine}) stage({stage}) dmu(id) gen(res) replacevars saving("saved.dta") replace\n'
            metrics='N N_eff mean_eff min_eff max_eff n_failed N_eff_opt mean_opt min_opt max_opt n_failed_primary n_failed_opt'
            if j==0:
                for metric in metrics.split():setup+=f'local baseline_{metric}=r({metric})\n'
            else:
                for metric in metrics.split():setup+=f'assert abs(r({metric})-`baseline_{metric}\')<1e-7\n'
            setup+=checks('imprecise',requested,'legacy','res res_opt res_mid',{'nirs':'drs','ndrs':'irs'}.get(rts,rts),saving=True)
            setup+='assert !missing(res) if selected\ncf `original_vars\' using `original\',all\n'
            if j==0:
                setup+='preserve\nkeep id res\nrename res expected\nsave `expected_order\'\nrestore\n'
            else:
                setup+='merge 1:1 id using `expected_order\',nogen\nassert missing(res)==missing(expected)\nassert abs(res-expected)<1e-7 if selected\n'
        def audit(d,rows=rows,orders=orders,sample=sample,inputs=inputs,outputs=outputs,bounded=bounded,ordinal=ordinal):
            byid={r['id']:r for r in rows};cells=0;matrices=0
            for tag,order in enumerate(orders):
                selected=[i for i in order if i in sample]
                for k,own in enumerate(selected,1):
                    with (d/f'transform_p{tag}_{k}.csv').open() as f:actual=list(csv.DictReader(f))
                    assert [int(r['id']) for r in actual]==selected
                    for got in actual:
                        i=int(got['id']);row=byid[i]
                        for var in (inputs+' '+outputs).split():
                            if var in bounded.split():
                                upper=(i!=own) if var in inputs.split() else (i==own)
                                expected=row[var+('hi' if upper else 'lo')]
                            elif var in ordinal.split():
                                expected=(len(selected) if row[var]<byid[own][var] else 1) if var in inputs.split() else int(row[var]<=byid[own][var])
                            else:expected=row[var]
                            assert float(got[var])==expected,(tag,own,i,var,got[var],expected)
                            cells+=1
                    matrices+=1
            return dict(transformed_matrices=matrices,oracle_cells=cells)
        run(name,setup,trace,audit)

    # Full 4! permutations, tied to analytic scores, not to another DEA call.
    # With constant cardinal output y=1, eligible peers have q_j <= q_i.
    # Input score=min eligible x / own x for every RTS. Output VRS/NIRS is 1;
    # output CRS/NDRS has the same min-x ratio. The four scores are 1,.5,1,1.
    for rts,ori,engine in itertools.product(['crs','vrs','nirs','ndrs'],['input','output'],['simplex','ipm','auto']):
        name=f'permutations_{rts}_{ori}_{engine}'
        if not chosen(name):continue
        setup=SMALL+'tempfile fixture\nsave `fixture\'\n'
        for k,order in enumerate(itertools.permutations(range(1,5))):
            setup+='use `fixture\',clear\n'+order_code(order)
            setup+=f'dea x=y q, ordinal(q) rts({rts}) orientation({ori}) engine({engine}) stage({1+k%2}) dmu(id) gen(res)\n'
            expect='1' if ori=='output' and rts in ['vrs','nirs'] else 'cond(id==2,.5,1)'
            setup+=f'assert !missing(res) & abs(res-({expect}))<1e-7\n'
        run(name,setup)

    # DMU labels must not substitute for observation identity.
    for label,engine,ori in itertools.product(['duplicate','string','missing'],['simplex','ipm','auto'],['input','output']):
        name=f'identity_{label}_{engine}_{ori}'
        setup=SMALL+{'duplicate':'gen label=1\n','string':'gen str30 label="기관 "+string(mod(id,2))\n','missing':'gen double label=cond(id==1,.a,id)\n'}[label]
        setup+=f'dea x=y q, yord(q) dmu(label) engine({engine}) orientation({ori}) gen(first)\ngsort -id\n'
        setup+=f'dea x=y q, ordinal(q) dmu(label) engine({engine}) orientation({ori}) gen(res)\n'
        setup+='assert abs(first-res)<1e-7\nassert abs(res-cond(id==2,.5,1))<1e-7\n'
        run(name,setup)
    # Missing ordinal values must be excluded before any threshold is read.
    for role,value,engine,ori in itertools.product(['input','output'],['.','.a','.z'],['simplex','ipm','auto'],['input','output']):
        name=f'missing_{role}_{value[1:] or "dot"}_{engine}_{ori}'
        formula='x q=y' if role=='input' else 'x=y q'
        setup=SMALL+f'replace q={value} if id==2\n'
        setup+=f'dea {formula} if id!=2, ordinal(q) engine({engine}) orientation({ori}) dmu(id) gen(control)\n'
        setup+='gsort -id\n'
        setup+=f'dea {formula}, ordinal(q) engine({engine}) orientation({ori}) dmu(id) gen(res)\n'
        setup+='assert r(N)==3\nassert missing(res) if id==2\nassert !missing(res) & abs(res-control)<1e-7 if id!=2\n'
        run(name,setup)
    for rts,ori,engine in itertools.product(['crs','vrs','nirs','ndrs'],['input','output'],['simplex','ipm','auto']):
        setup=dataset(BASE)
        setup+=f'dea x qx xb=y qy qz yb if id==3, bounded(xb yb) ordinal(qx qy qz) rts({rts}) orientation({ori}) engine({engine}) dmu(id) gen(res)\n'
        setup+='assert r(N)==1 & r(N_eff)==1 & r(n_failed)==0\nassert abs(res-1)<1e-7 if id==3\nassert missing(res) if id!=3\n'
        run(f'singleton_{rts}_{ori}_{engine}',setup)

    summary=dict(total=len(results),passed=sum(r['passed'] for r in results),failed=sum(not r['passed'] for r in results),transformed_matrices=sum(r.get('transformed_matrices',0) for r in results),oracle_cells=sum(r.get('oracle_cells',0) for r in results))
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
    raise SystemExit(summary['failed']>0)

if __name__=='__main__':main()
