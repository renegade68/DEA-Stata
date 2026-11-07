#!/usr/bin/env python3
"""Numeric bounds + legacy options: equivalence, independent transforms and protection.

Python standard library, licensed Stata 16+. Optional pre-mixed (2.0.39) controls
use string bounded() plus ordinal(), including the F30 order fix. Test copies
export transformed data; production and installed ado files are never modified.
"""
from _support import package_files
import argparse,csv,hashlib,itertools,json,random,re,shutil,subprocess,time
from pathlib import Path
from test_fdh_idea_returns import checks
from test_ordinal_order import BASE
from test_saveref import instrument
REPO=Path(__file__).resolve().parents[1]
CONFIGS={'xboth':['xlo','xhi'],'yboth':['ylo','yhi'],'both':['xlo','xhi','ylo','yhi'],
         'xlo':['xlo'],'xhi':['xhi'],'ylo':['ylo'],'yhi':['yhi'],'cross':['xlo','yhi']}
DATA='input double(id x qx y qy)\n'+'\n'.join(' '.join(str(r[k]) for k in ['id','x','qx','y','qy']) for r in BASE)+'\nend\n'
DATA+='recast int x y\ngen double xl=x-.25\ngen double xh=x+.75\ngen double yl=y-.25\ngen double yh=y+.75\n'
DATA+='gen strL xb="["+strtrim(string(xl,"%24.17g"))+","+strtrim(string(xh,"%24.17g"))+"]"\n'
DATA+='gen strL yb="["+strtrim(string(yl,"%24.17g"))+","+strtrim(string(yh,"%24.17g"))+"]"\n'


def encode(plan):
    # Freeze all textual representations before replacing nominal columns,
    # so controls also support bounds sourced from other model dimensions.
    s=''
    for j,(v,lo,hi) in enumerate(plan):s+=f'gen strL enc{j}="["+strtrim(string({lo},"%24.17g"))+","+strtrim(string({hi},"%24.17g"))+"]"\n'
    for j,(v,lo,hi) in enumerate(plan):s+=f'drop {v}\nrename enc{j} {v}\n'
    return s


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--baseline-ado',type=Path);ap.add_argument('--only')
    args=ap.parse_args()
    if args.baseline_ado and not (args.baseline_ado/'dea.ado').is_file():ap.error('Invalid pre-mixed ado directory')
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'current_ado';pkg.mkdir();hashes={}
    for n in package_files(REPO):
        raw=(REPO/n).read_bytes();(pkg/Path(n).name).write_bytes(raw);hashes[n]=hashlib.sha256(raw).hexdigest()
    def trace_copy(source,target):
        shutil.copytree(source,target);p=target/'dea.ado';s=p.read_text();anchor="            foreach _drop_output in `__idea_score' {";assert s.count(anchor)==1
        s=s.replace(anchor,'''            quietly format id `inputs' `outputs' %24.17g
            quietly export delimited id `inputs' `outputs' using "transform_`dmui'.csv", replace datafmt
'''+anchor);p.write_text(s)
    curtrace=out/'current_trace';oldtrace=out/'legacy_trace';trace_copy(pkg,curtrace);trace_copy(args.baseline_ado or pkg,oldtrace)
    retry=out/'retry_ado';instrument(pkg,retry)
    manifest=dict(source_sha256=hashes,stata=args.stata,baseline=str(args.baseline_ado),legacy_control='String encoding uses fixed endpoint values before renaming; pre-mixed worker includes F30 fix')
    if args.baseline_ado:manifest['baseline_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in args.baseline_ado.iterdir() if p.is_file()}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2));results=[]
    def chosen(n):return not args.only or re.search(args.only,n)
    def run(name,script,source=pkg,audit=None,force=False):
        if not force and not chosen(name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 120\n'
        for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/key).mkdir();head+=f'sysdir set {key} "{d/key}"\n'
        head+=f'adopath ++ "{source.resolve()}"\n';(d/'test.do').write_text(head+script+'\ndi "TEST_PASS"\n');start=time.monotonic();row=dict(test=name,passed=False)
        try:
            proc=subprocess.run([args.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=300)
            if proc.returncode or not (d/'test.log').exists():
                (d/'process.json').write_text(json.dumps(dict(returncode=proc.returncode,stdout=proc.stdout.decode(errors='replace'),stderr=proc.stderr.decode(errors='replace')),indent=2))
            if not (d/'test.log').exists():raise RuntimeError(f'Stata did not create test.log; see process.json (return code {proc.returncode})')
            log=(d/'test.log').read_text(errors='replace')
            if proc.returncode or 'TEST_PASS' not in log.splitlines():raise RuntimeError(log[-2200:])
            if (d/'target_before.dta').exists():assert (d/'target.dta').read_bytes()==(d/'target_before.dta').read_bytes()
            if audit:row.update(audit(d))
            row['passed']=True
        except Exception as e:row['error']=repr(e)
        row['seconds']=round(time.monotonic()-start,3);results.append(row);(out/'results.json').write_text(json.dumps(results,indent=2));print(json.dumps(row),flush=True)

    # Same named dimensions, original observations and effective sample in the
    # two encodings. Controls use the already repaired legacy string route.
    for i,((kind,keys),ordkind,rts,ori,engine) in enumerate(itertools.product(CONFIGS.items(),['input','output','both'],['crs','vrs','nirs','ndrs'],['input','output'],['simplex','ipm','auto'])):
        name=f'grid_{kind}_{ordkind}_{rts}_{ori}_{engine}'
        if not chosen(name):continue
        ordinal={'input':'qx','output':'qy','both':'qx qy'}[ordkind];plan=[];lists={}
        for side,dims in [('x',['x','qx']),('y',['y','qy'])]:
            if not any(k.startswith(side) for k in keys):continue
            for end,suffix in [('lo','l'),('hi','h')]:
                if side+end in keys:lists[side+end]=[side+suffix,dims[1]]
            for j,v in enumerate(dims):
                if v not in ordinal.split():plan.append((v,lists.get(side+'lo',dims)[j],lists.get(side+'hi',dims)[j]))
        opts=' '.join(k+'('+' '.join(v)+')' for k,v in lists.items());numeric=' '.join(v for v,lo,hi in plan)
        setup=DATA;order=list(range(1,7))
        if i%3==0:setup+='gsort -id\n';order.reverse()
        condition=['1','id!=3','inrange(_n,2,6)','mod(id,2) & inrange(_n,2,6)'][i%4]
        qualifier=['','if id!=3','in 2/6','if mod(id,2) in 2/6'][i%4]
        if i%11==0:setup+=f'replace {ordinal.split()[0]}=.a if id==6\n'
        if i%17==0 and 'xlo' in keys and 'xhi' in keys:setup+='replace x=.z if id==2\n'
        setup+=f'gen byte selected=({condition}) & !missing('+','.join(ordinal.split())+')\n'
        selected=[v for pos,v in enumerate(order,1) if [True,v!=3,2<=pos<=6,v%2 and 2<=pos<=6][i%4] and not(i%11==0 and v==6)]
        stage=1+(i//3)%2;controls=setup+encode(plan)
        controls+=f'dea x qx=y qy if selected, bounded({numeric}) ordinal({ordinal}) rts({rts}) orientation({ori}) engine({engine}) stage({stage}) dmu(id) gen(res)\n'
        controls+='assert !missing(res) if selected\nkeep res res_opt res_mid\nsave "expected.dta"\n'
        run('control_'+name,controls,oldtrace,force=True);path=out/('control_'+name)/'expected.dta'
        setup+='quietly ds\nlocal original_vars "`r(varlist)\'"\ntempfile original\nsave `original\'\n'
        if i%2:setup+='foreach v in res res_opt res_mid {\n gen str10 `v\'="old"\n}\n'
        alias=f'ordinal({ordinal})' if i%3==0 else (f'xord(qx)' if ordkind=='input' else 'yord(qy)' if ordkind=='output' else 'xord(qx) yord(qy)')
        if i%3==2:alias+=f' ordinal({ordinal})'
        model=['','model(te)','model(imprecise)'][i%3];requested='imprecise' if i%3==2 else 'te'
        prefix=f'dea x qx=y qy {qualifier},' if (i//3)%2 else f'dea y qy {qualifier}, inputs(x qx)'
        setup+='local savevars id x qx y qy\n'
        setup+=f'{prefix} {opts} {alias} {model} rts({rts}) orientation({ori}) engine({engine}) stage({stage}) dmu(id) gen(res) saving("saved.dta")'+(' replacevars' if i%2 else '')+'\n'
        setup+=f'assert "`r(idea_numeric_vars)\'"=="{numeric}"\n'
        setup+=checks('imprecise',requested,'legacy','res res_opt res_mid',{'nirs':'drs','ndrs':'irs'}.get(rts,rts),saving=True)
        setup+=f'cf res res_opt res_mid using "{path}",all\ncf `original_vars\' using `original\',all\n'
        def audit(d,selected=selected,plan=plan,ordinal=ordinal,name=name):
            byid={r['id']:dict(r,xl=r['x']-.25,xh=r['x']+.75,yl=r['y']-.25,yh=r['y']+.75) for r in BASE};cells=0
            for k,own in enumerate(selected,1):
                with (d/f'transform_{k}.csv').open() as f:got=list(csv.DictReader(f))
                with (out/('control_'+name)/f'transform_{k}.csv').open() as f:old=list(csv.DictReader(f))
                assert got==old,('encoding mismatch',k)
                assert [int(v['id']) for v in got]==selected
                for row in got:
                    who=int(row['id'])
                    for var in ['x','qx','y','qy']:
                        spec=next((p for p in plan if p[0]==var),None)
                        if spec:
                            upper=(who!=own) if var in ['x','qx'] else (who==own)
                            expected=byid[who][spec[2] if upper else spec[1]]
                        elif var in ordinal.split():expected=(len(selected) if byid[who][var]<byid[own][var] else 1) if var=='qx' else int(byid[who][var]<=byid[own][var])
                        else:expected=byid[who][var]
                        assert float(row[var])==expected,(k,who,var,row[var],expected);cells+=1
            return dict(transform_pairs=len(selected),oracle_cells=cells)
        run(name,setup,curtrace,audit)

    # Previously rejected combinations, with unchanged dataset after rejection.
    if args.baseline_ado:
        for kind,keys in CONFIGS.items():
            for form in ['A','B']:
                name=f'old_reject_{kind}_{form}'
                opts=' '.join(k+'('+({'xlo':'xl qx','xhi':'xh qx','ylo':'yl qy','yhi':'yh qy'}[k])+')' for k in keys)
                prefix='dea x qx=y qy,' if form=='A' else 'dea y qy, inputs(x qx)'
                run(name,DATA+'tempfile original\nsave `original\'\n'+f'capture noisily {prefix} {opts} ordinal(qx qy)\nassert _rc==198\ncf _all using `original\',all\n',args.baseline_ado)

    def reject(setup,opts,rc,form,replace=True,formula=('x qx','y qy')):
        inputs,outputs=formula
        s=setup+'''gen str12 res="old"
gen str12 res_opt="old opt"
gen str12 res_mid="old mid"
label variable res "original label"
char res[test] "original characteristic"
char _dta[test] "dataset characteristic"
tempfile original
save `original'
save "target.dta"
copy "target.dta" "target_before.dta"
'''
        prefix=f'dea {inputs}={outputs},' if form=='A' else f'dea {outputs}, inputs({inputs})'
        s+=f'capture noisily {prefix} {opts} gen(res) saving("target.dta") replace'+(' replacevars' if replace else '')+f'\nassert _rc=={rc}\n'
        s+='''cf _all using `original',all
assert "`:variable label res'"=="original label"
assert "`:char res[test]'"=="original characteristic"
assert "`:char _dta[test]'"=="dataset characteristic"
capture confirm variable _dmu,exact
assert _rc==111
'''
        return s
    for endpoint,value,form,replace in itertools.product(['xl','xh','yl','yh'],['.','.a','.z'],['A','B'],[False,True]):
        run(f'invalid_{endpoint}_{value[1:] or "dot"}_{form}_{replace}',reject(DATA+f'replace {endpoint}={value} in 3\n','xlo(xl qx) xhi(xh qx) ylo(yl qy) yhi(yh qy) ordinal(qx qy)',459,form,replace))
    for role,form,replace in itertools.product(['x','y'],['A','B'],[False,True]):
        setup=DATA+f'replace {role}l={role}h+1 in 3\n'
        run(f'reversed_{role}_{form}_{replace}',reject(setup,'xlo(xl qx) xhi(xh qx) ylo(yl qy) yhi(yh qy) ordinal(qx qy)',459,form,replace))
    invalid=[('xlo(xl xl) ordinal(qx)',198,DATA,('x qx','y qy')),
             ('xhi(xh xh) xord(qx)',198,DATA,('x qx','y qy')),
             ('ylo(yl yl) ordinal(qy)',198,DATA,('x qx','y qy')),
             ('yhi(yh yh) yord(qy)',198,DATA,('x qx','y qy')),
             ('xlo(xl qx) bounded(x)',198,DATA,('x qx','y qy')),
             ('xlo(xl xl qx) bounded(xb) ordinal(qx)',198,DATA,('x xb qx','y qy')),
             ('xlo(xl) ordinal(qx)',198,DATA,('x qx','y qy')),
             ('yhi(yh) ordinal(qy)',198,DATA,('x qx','y qy')),
             ('xlo(xb qx) ordinal(qx)',109,DATA,('x qx','y qy')),
             ('xlo(xl qx) xhi(xh qx) ordinal(qx)',109,DATA+'tostring x,replace\n',('x qx','y qy')),
             ('xlo(xl xl qx) ordinal(qx)',198,DATA,('x x qx','y qy')),
             ('xlo(xl qx) ordinal(qx)',198,DATA,('x qx','x qy')),
             ('xlo(xb qx)',109,DATA,('x qx','y qy'))]
    for i,((opts,rc,setup,formula),form,replace) in enumerate(itertools.product(invalid,['A','B'],[False,True])):
        run(f'reject_plan_{i:03d}',reject(setup,opts,rc,form,replace,formula))
    for opt in ['xlo(xl qx)','xhi(xh qx)','ylo(yl qy)','yhi(yh qy)']:
        dim=opt[0];setup=DATA+f'replace {dim}=.a in 3\n'
        for form in ['A','B']:run('missing_fallback_'+opt[:3]+'_'+form,reject(setup,opt+' ordinal(qx qy)',459,form))

    # Mixed representations on the same side, frozen cross-dimension endpoints,
    # and fractional endpoints applied to int nominal variables.
    extras=[('same_side',DATA,'x xb qx','y yb qy','xb yb','qx qy','xlo(xl xb qx) xhi(xh xb qx) ylo(yl yb qy) yhi(yh yb qy)',[('x','xl','xh'),('y','yl','yh')]),
            ('alias_sources',DATA+'gen double cap=10\n','x y qx','qy','','qx qy','xlo(y x qx) xhi(cap cap qx)',[('x','y','cap'),('y','x','cap')]),
            ('precision',DATA+'replace xl=x+1e-10\nreplace xh=x+2e-10\n','x qx','y qy','','qx qy','xlo(xl qx) xhi(xh qx)',[('x','xl','xh')]),
            ('ordinal_source',DATA+'replace xl=.25\n','x qx','y qy','','qx qy','xlo(xl qx) xhi(qx qx)',[('x','xl','qx')]),
            ('bounded_only',DATA,'x xb','y yb','xb yb','','xlo(xl xb) xhi(xh xb) ylo(yl yb) yhi(yh yb)',[('x','xl','xh'),('y','yl','yh')])]
    for (kind,setup,ins,outs,bounded,ordinal,opts,plan),engine,ori,form in itertools.product(extras,['simplex','ipm','auto'],['input','output'],['A','B']):
        name=f'extra_{kind}_{engine}_{ori}_{form}'
        s=setup+'tempfile fixture scores\nsave `fixture\'\n'+encode(plan)
        allbounds=' '.join([bounded]+[v for v,lo,hi in plan])
        s+=f'dea {ins}={outs}, bounded({allbounds}) ordinal({ordinal}) engine({engine}) orientation({ori}) dmu(id) gen(control)\nkeep id control\nsave `scores\'\nuse `fixture\',clear\n'
        prefix=f'dea {ins}={outs},' if form=='A' else f'dea {outs}, inputs({ins})'
        s+=f'{prefix} {opts} bounded({bounded}) ordinal({ordinal}) engine({engine}) orientation({ori}) dmu(id) gen(res)\nmerge 1:1 id using `scores\',nogen\nassert res==control\n'
        run(name,s)
    for how in ['if','in','ordinal_missing','none']:
        s=DATA+'replace xl=.a in 3\n'
        qual={'if':'if id!=3','in':'in 4/6','ordinal_missing':'','none':'if 0'}[how]
        if how=='ordinal_missing':s+='replace qx=.z in 3\n'
        cmd=f'dea x qx=y qy {qual}, xlo(xl qx) xhi(xh qx) ordinal(qx qy) dmu(id) gen(res)'
        s+=('capture noisily '+cmd+'\nassert _rc==2000\n') if how=='none' else (cmd+'\nassert missing(res) if id==3\nassert !missing(res) if '+('id!=3' if how!='in' else 'id>=4')+'\n')
        run('excluded_'+how,s)
    for ori in ['input','output']:
        s=DATA+f'dea x qx=y qy, xlo(xl qx) xhi(xh qx) ordinal(qx qy) engine(ipm) orientation({ori}) dmu(id) gen(control)\n'
        s+=f'dea x qx=y qy, xlo(xl qx) xhi(xh qx) ordinal(qx qy) engine(auto) orientation({ori}) dmu(id) gen(res)\n'
        s+='assert "`r(engine)\'"=="ipm"\nassert res==control\n'
        run('retry_'+ori,s,retry)
    for rts,ori,engine in itertools.product(['crs','vrs','nirs','ndrs'],['input','output'],['simplex','ipm','auto']):
        opts=f'xlo(xl qx) xhi(xh qx) ylo(yl qy) yhi(yh qy) ordinal(qx qy) rts({rts}) orientation({ori}) engine({engine}) dmu(id)'
        s=DATA+f'dea x qx=y qy, {opts} gen(control)\ngsort -id\n'
        s+=f'dea y qy, inputs(x qx) {opts} gen(res)\nassert !missing(res) & abs(res-control)<1e-7\n'
        run(f'permutation_{rts}_{ori}_{engine}',s)
    # New result metadata must not survive an unrelated subsequent call.
    for following in ['','ordinal(qx qy)','xlo(xl qx) xhi(xh qx)']:
        s=DATA+'dea x qx=y qy, xlo(xl qx) xhi(xh qx) ordinal(qx qy) dmu(id) gen(first)\n'
        s+='assert "`r(idea_numeric_vars)\'"=="x"\n'
        s+=f'dea x qx=y qy, {following} dmu(id) gen(res)\nassert "`r(idea_numeric_vars)\'"==""\n'
        run('stale_'+('te' if not following else 'ordinal' if following.startswith('ordinal') else 'interval'),s)
    for suffix in ['', '_opt', '_mid']:
        target='res'+suffix
        s=DATA+f'rename xl {target}\ntempfile original\nsave `original\'\n'
        s+=f'capture noisily dea x qx=y qy, xlo({target} qx) xhi(xh qx) ordinal(qx qy) gen(res) replacevars\nassert _rc==110\ncf _all using `original\',all\n'
        run('protect_endpoint_'+(suffix or 'main'),s)

    for form,engine,ori in itertools.product(['A','B'],['simplex','ipm','auto'],['input','output']):
        s=DATA+'gen byte selected=id!=3\n'
        opts=f'xlo(xl qx) xhi(xh qx) ylo(yl qy) yhi(yh qy) ordinal(qx qy) engine({engine}) orientation({ori}) gen(res)'
        cmd=f'dea x qx=y qy if selected using "saved.dta", {opts}' if form=='A' else f'dea y qy if selected, inputs(x qx) {opts} saving("saved.dta")'
        s+=cmd+'\nassert r(generated_dmu)==1\nassert "`r(idea_numeric_vars)\'"=="x y"\nlocal savevars _dmu x qx y qy\n'
        s+=checks('imprecise','te','legacy','res res_opt res_mid',saving=True)
        s+=f'dea x qx=y qy if selected, {opts} dmu(_dmu) replacevars saving("saved.dta") replace\n'
        s+=checks('imprecise','te','legacy','res res_opt res_mid',saving=True)
        run(f'using_autoid_{form}_{engine}_{ori}',s)

    # Endpoint sources are auxiliary data, not ordinary DEA output columns.
    for low,high,engine,ori,form in itertools.product(['CRS_TE','SCALE'],['VRS_TE','RTS_NUM','RTS'],['simplex','ipm','auto'],['input','output'],['A','B']):
        s=DATA+f'rename xl {low}\nrename xh {high}\n'
        s+='tempfile original\nsave `original\'\n'
        prefix='dea x qx=y qy,' if form=='A' else 'dea y qy, inputs(x qx)'
        s+=f'{prefix} xlo({low} qx) xhi({high} qx) ordinal(qx qy) rts(vrs) engine({engine}) orientation({ori}) dmu(id) gen(res)\n'
        s+=f'assert !missing(res)\ncf id x qx y qy {low} {high} using `original\',all\n'
        run(f'endpoint_name_{low}_{high}_{engine}_{ori}_{form}',s)

    summary=dict(total=len(results),passed=sum(r['passed'] for r in results),failed=sum(not r['passed'] for r in results),transform_pairs=sum(r.get('transform_pairs',0) for r in results),oracle_cells=sum(r.get('oracle_cells',0) for r in results))
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True);raise SystemExit(summary['failed']>0)

if __name__=='__main__':main()
