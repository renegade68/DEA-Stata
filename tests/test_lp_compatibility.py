#!/usr/bin/env python3
"""Standalone DEA alias and table/matrix LP/MILP: isolated, analytic and enumerated checks."""
from _support import package_files
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import random
import re
import subprocess
import time

REPO = Path(__file__).resolve().parents[1]
DEA = '''input double id x y
1 1 1
2 2 3
3 4 4
4 5 3
end
gen str30 tag=cond(id==4,"omit","keep = yes, now")
'''
TABLE = '''input double x1 x2 str2 rel double rhs
3 2 "=" 991
2 1 "<=" 4
1 2 "<=" 4
end
gen str20 rowid="row"+string(_n)
'''
MAT = 'matrix C=(3,2)\nmatrix A=(2,1\\1,2)\nmatrix B=(4\\4)\n'


def protect(command, rc=None):
    s='''scalar caller_scalar=991
matrix caller_matrix=(1,2\\3,4)
datasignature
local before "`r(datasignature)'"
local more_before "`c(more)'"
local rng_before "`c(rngstate)'"
'''
    s += command+'\n' if rc is None else 'capture noisily '+command+f'\nassert _rc=={rc}\n'
    s += '''datasignature
assert "`before'"=="`r(datasignature)'"
assert "`more_before'"=="`c(more)'"
assert "`rng_before'"=="`c(rngstate)'"
assert scalar(caller_scalar)==991
mata: assert(all(st_matrix("caller_matrix"):==(1,2\\3,4)))
'''
    return s


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata',required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--legacy-lp',type=Path)
    ap.add_argument('--legacy-dea',type=Path)
    ap.add_argument('--only')
    a=ap.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    pkg=out/'ado';pkg.mkdir();hashes={}
    for n in package_files(REPO):
        raw=(REPO/n).read_bytes();(pkg/Path(n).name).write_bytes(raw);hashes[n]=hashlib.sha256(raw).hexdigest()
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,stata=a.stata),indent=2)+'\n')
    rows=[];oracles=[]

    def run(name,s,source=pkg,force=False):
        if not force and a.only and not re.search(a.only,name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 120\n'
        for p in ['PLUS','PERSONAL','SITE','OLDPLACE']:
            (d/p).mkdir();head+=f'sysdir set {p} "{d/p}"\n'
        head+=f'adopath ++ "{source.resolve()}"\n'
        (d/'test.do').write_text(head+s+'\ndi "TEST_PASS"\n');row=dict(test=name,passed=False);start=time.monotonic()
        try:
            proc=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=180)
            log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else ''
            if proc.returncode or not log:(d/'process.json').write_text(json.dumps(dict(rc=proc.returncode,stderr=proc.stderr.decode(errors='replace')),indent=2))
            assert proc.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2000:] or proc.stderr.decode(errors='replace')
            if (d/'target_before.dta').exists():assert (d/'target.dta').read_bytes()==(d/'target_before.dta').read_bytes()
            assert not (d/'lp.log').exists() or source!=pkg
            row['passed']=True
        except Exception as e:row['error']=str(e)
        row['seconds']=round(time.monotonic()-start,3);rows.append(row)
        (out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(row),flush=True)

    # Exact parity with current dea given the same explicit options, in fresh
    # and already-loaded sessions. The current wrapper must not be overwritten.
    for i,(load,form,rts,engine,stage) in enumerate(itertools.product(
            ['fresh','dea','run'],['A','B'],['crs','vrs','nirs','ndrs'],['simplex','ipm','auto'],[1,2])):
        q=['','if id!=4','in 2/4','if tag=="keep = yes, now" in 1/4'][i%4]
        spec=f'x=y {q},' if form=='A' else f'y {q}, inputs(x)'
        options=f'dmu(id) rts({rts}) engine({engine}) stage({stage})'
        pre=''
        if load=='dea':pre='dea x=y, dmu(id) gen(initial)\n'
        if load=='run':pre=f'run "{pkg/"dea.ado"}"\n'
        s=DEA+pre+f'dea_lp {spec} {options} gen(res)\nassert "`r(cmd)\'"=="dea_lp"\n'
        s+=f'dea {spec} {options} gen(ctrl) replacevars\nassert res==ctrl\n'
        run(f'alias_{load}_{form}_{rts}_{engine}_{stage}',s)

    for form,solver in itertools.product(['A','B'],['','method(ipm)','met(simplex)','eng(auto)','engine(ipm) method(simplex)']):
        spec='x=y,' if form=='A' else 'y, inputs(x)'
        s=DEA+f'dea_lp {spec} dmu(id) gen(res) {solver}\n'
        s+='assert "`r(rts)\'"=="vrs"\nassert abs(res-cond(id==4,.4,1))<1e-6\n'
        expected='simplex' if solver in ['met(simplex)','eng(auto)'] else 'ipm'
        s+=f'assert "`r(engine)\'"=="{expected}"\n'
        run('alias_default_'+form+'_'+str(['','method(ipm)','met(simplex)','eng(auto)','engine(ipm) method(simplex)'].index(solver)),s)
    for bad,rc in [('tolerance(0)',198),('engine(bad)',198),('rts(bad)',198),('scale(01)',198)]:
        run('alias_invalid_'+bad.split('(')[0],DEA+protect('dea_lp x=y, dmu(id) gen(z) '+bad,rc))
    run('alias_result_collision',DEA+'gen double res=991\n'+protect('dea_lp x=y, dmu(id) gen(res)',110))
    run('alias_result_replace',DEA+'gen double res=991\ndea_lp x=y, dmu(id) gen(res) replacevars\nassert abs(res-cond(id==4,.4,1))<1e-6\n')
    run('alias_repeated_load',DEA+f'run "{pkg/"dea.ado"}"\n'+
        'dea_lp x=y, dmu(id) gen(first)\nprogram drop dea_lp\n'+
        'dea_lp x=y, dmu(id) gen(second) replacevars\nassert first==second\n')

    # Both table and matrix APIs: analytic optima, senses, bounds, empty constraint set.
    for engine,scale,mode in itertools.product(['simplex','ipm','auto'],['none','row'],['lp','mixed','integer']):
        extra={'lp':'','mixed':'intvars(x2)','integer':'intvars(x1 x2)'}[mode]
        answer={'lp':'(4/3,4/3)','mixed':'(1.5,1)','integer':'(2,0)'}[mode]
        obj={'lp':'20/3','mixed':'6.5','integer':'6'}[mode]
        cmd=f'lp x1 x2, max {extra} engine({engine}) scaling({scale}) rownames(rowid)'
        s=TABLE+cmd+f'\nassert r(status)==0\nassert abs(r(objective)-({obj}))<1e-6\n'
        s+='assert r(nvars)==2 & r(nconstraints)==2 & r(nslacks)==2\nassert r(tolerance)==1e-8\n'
        s+='matrix S=r(x)\nmatrix R=r(lprslt)\nmatrix P=r(problem)\n'
        s+=f'mata: assert(max(abs(st_matrix("S"):-{answer}))<1e-6)\n'
        s+='assert P[1,3]==0 & rhs[1]==991\nassert abs(R[1,4]-(4-2*S[1,1]-S[1,2]))<1e-6\n'
        s+='assert abs(R[1,5]-(4-S[1,1]-2*S[1,2]))<1e-6\n'
        s+=protect(cmd)
        run(f'table_{mode}_{engine}_{scale}',s)
        if mode!='lp':
            ints='2' if mode=='mixed' else '1 2'
            run(f'matrix_{mode}_{engine}_{scale}',MAT+
                f'milp_bb, c(C) aub(A) bub(B) integer({ints}) lpengine({engine}) lpscaling({scale})\n'+
                f'assert r(status)==0 & abs(r(objective)-({obj}))<1e-6\n')

    for engine,scale in itertools.product(['simplex','ipm','auto'],['none','row']):
        for relation in ['<=','<','>=','>','=']:
            # One variable: max x with x<=2, min x with x>=2 or x=2.
            sense='max' if relation in ['<=','<'] else 'min'
            body=f'input double x1 str2 rel double rhs\n1 "" .\n1 "{relation}" 2\nend\n'
            run(f'relation_{engine}_{scale}_{relation.replace("<","lt").replace(">","gt").replace("=","eq")}',body+
                f'lp x1, {sense} engine({engine}) scaling({scale})\nassert r(status)==0 & abs(r(objective)-2)<1e-6\n')
        body='input double x1 x2 str2 rel double rhs\n1 2 "" .\n1 1 "=" 3\nend\nmatrix U=(2,3)\n'
        run(f'min_equal_{engine}_{scale}',body+f'lp x1 x2, min ub(U) engine({engine}) scaling({scale})\nassert r(status)==0 & abs(r(objective)-4)<1e-6\n')
        body='input double x1 str2 rel double rhs\n1 "" .\nend\nmatrix L=(-2)\nmatrix U=(3)\n'
        run(f'bound_only_{engine}_{scale}',body+f'lp x1, min lb(L) ub(U) engine({engine}) scaling({scale})\nassert r(status)==0 & abs(r(objective)+2)<1e-6\n')

    for mode,extra in [('lp',''),('mixed','intvars(x2)'),('integer','intvars(x1 x2)')]:
        for q in ['','if include','in 2/4','if include in 1/5']:
            s=TABLE+'gen byte include=1\nset obs 5\nreplace include=0 in 4/5\ngen byte order=_n\nreplace order=0 in 4\nsort order\n'
            # Every qualified variant selects exactly the objective and two rows.
            qual=q or 'if include==1'
            run('selection_'+mode+'_'+str(['','if include','in 2/4','if include in 1/5'].index(q)),s+
                protect(f'lp x1-x2 {qual}, max {extra}'))

    # Using a different dataset must preserve the caller even when variable names
    # exist only in the file; filenames include equals signs, commas and spaces.
    for mode,extra,obj in [('lp','',20/3),('mixed','intvars(x2)',6.5),('integer','intvars(x1 x2)',6)]:
        s=TABLE+'save "input = a,b file.dta"\nclear\nset obs 2\ngen str20 caller="untouched"\n'
        cmd=f'lp x1 x2 using "input = a,b file.dta", max {extra} saving("solution = a,b file")'
        s+=protect(cmd)
        s+='preserve\nuse "solution = a,b file.dta", clear\nassert _N==2\nassert parameter=="x1" in 1\nassert parameter=="x2" in 2\n'
        s+=f'assert status==0 & abs(objective-{obj})<1e-6\nrestore\n'
        run('using_saving_'+mode,s)
    run('existing_file_protected',TABLE+'save "target.dta"\ncopy "target.dta" "target_before.dta"\n'+protect('lp x1 x2, max saving(target)',602))
    run('saving_replace',TABLE+'save "target.dta"\nlp x1 x2, max saving(target) replace\nassert _N==3\npreserve\nuse target, clear\nassert status==0\nrestore\n')
    run('using_bad_restores',TABLE+'replace x1=. in 2\nsave input\nclear\nset obs 2\ngen byte caller=99\n'+protect('lp x1 x2 using input, max',459))
    run('saving_failure_restores',TABLE+protect('lp x1 x2, max saving("nonexistent/result")',603))

    invalid=[('sense_missing','','',198),('sense_both','','min max',198),
        ('bad_engine','','max engine(wrong)',198),('bad_scaling','','max scaling(wrong)',198),
        ('bad_relation','replace rel="??" in 2\n','max',198),
        ('coefficient_missing','replace x1=.a in 2\n','max',459),
        ('objective_missing','replace x1=. in 1\n','max',459),
        ('rhs_missing','replace rhs=.z in 3\n','max',459),
        ('rhs_same','','max rhs(x1)',198),('integer_outside','','max intvars(rhs)',198),
        ('empty_sample','','if 0, max',2000),
        ('tolerance_conflict','','max tol1(1e-14) tol2(1e-8)',198)]
    for name,setup,opts,rc in invalid:
        cmd='lp x1 x2 '+opts if name=='empty_sample' else 'lp x1 x2, '+opts
        run('invalid_'+name,TABLE+setup+protect(cmd,rc))
    for option,bad in itertools.product(['tolerance','tol1','tol2','inttol'],['0','-1','.','.a','.z']):
        run(f'invalid_{option}_{bad.replace(".","dot")}',TABLE+protect(f'lp x1 x2, max {option}({bad})',198))
    for n in ['0','-1','.']:
        run('invalid_nodes_'+n,TABLE+protect(f'lp x1 x2, max maxnodes({n})',198))
    for opt in ['tolerance(1e-8)','tol1(1e-8)','tol2(1e-8)','tol(1e-8)','tol1(1e-8) tol2(1e-8)','tolerance(1e-8) tol1(1e-8) tol2(1e-8)']:
        run('tolerance_alias_'+str(len(opt))+'_'+opt.split('(')[0],TABLE+f'lp x1 x2, max {opt}\nassert r(status)==0 & r(tolerance)==1e-8\n')

    # Failure states of continuous and integer models, without conflating rc=0
    # with an optimal solution. Inputs and caller objects must remain intact.
    for engine in ['simplex','ipm','auto']:
        for mode in ['lp','milp']:
            extra='' if mode=='lp' else 'intvars(x1)'
            body='input double x1 str2 rel double rhs\n1 "" 0\n1 "<=" 0\n1 ">=" 1\nend\n'
            cmd=f'lp x1, max {extra} engine({engine})'
            run(f'infeasible_{mode}_{engine}',body+cmd+'\nassert r(status)==1\nassert missing(r(objective))\n'+protect(cmd))
            body='input double x1 str2 rel double rhs\n1 "" 0\n1 ">=" 0\nend\n'
            # Stata 16.1's direct LinearProgram call reports infeasible for this
            # unbounded one-variable LP. Retain that independent native control;
            # keep the known solver classification separate from successful solves.
            direct=''
            if engine=='ipm':
                direct='''mata:
S=LinearProgram()
S.setCoefficients(1)
S.setMaxOrMin("max")
S.setInequality(-1,0)
S.setBounds(0,.)
S.setTol(1e-8)
S.setTrace("off")
v=S.optimize()
st_numscalar("native_report",S.errorcode())
assert(S.converged()==0)
end
assert inlist(scalar(native_report),1,2)
'''
            expected='2'  # Exact origin and positive coordinate ray override the raw report.
            run(f'unbounded_{mode}_{engine}',body+direct+f'lp x1, max {extra} engine({engine})\nassert r(status)=={expected}\nassert missing(r(objective))\n')
        run('integer_only_infeasible_'+engine,'input double x1 str2 rel double rhs\n1 "" .\n1 "=" .5\nend\n'+
            f'lp x1, min intvars(x1) engine({engine})\nassert r(status)==1\n')
        for api in ['lp','matrix']:
            s=TABLE+MAT
            cmd=(f'lp x1 x2, max intvars(x1 x2) engine({engine}) maxnodes(1)' if api=='lp' else
                 f'milp_bb, c(C) aub(A) bub(B) integer(1 2) lpengine({engine}) maxnodes(1)')
            run(f'node_limit_{api}_{engine}',s+cmd+'\nassert r(status)==4 & r(hasinc)==0\nassert missing(r(objective))\n')
        run('node_limit_with_incumbent_'+engine,TABLE+
            f'lp x1 x2, max intvars(x1 x2) engine({engine}) maxnodes(2)\n'+
            'assert r(status)==3 & r(hasinc)==1\nassert abs(r(objective)-6)<1e-6\n')
        run('forced_fallback_'+engine,MAT+
            ('milp_bb, c(C) aub(A) bub(B) integer(1 2) lpengine(auto) forcelpfallback\nassert r(status)==0 & r(lp_fallbacks)>0\n' if engine=='auto' else
             f'capture noisily milp_bb, c(C) aub(A) bub(B) integer(1 2) lpengine({engine}) forcelpfallback\nassert _rc==198\n'))

    for engine,api,shape in itertools.product(['simplex','ipm','auto'],['lp','milp'],['bounds','bounds_matrix','cost','rhs','missing_cost','missing_lower','extended_upper','inverted']):
        s=MAT+'matrix L=(0,0)\nmatrix U=(5,5)\n'
        change={'bounds':'matrix L=(0,0,0)\n','bounds_matrix':'matrix L=(0,0\0,0)\n','extended_upper':'matrix U=(.a,5)\n','cost':'matrix C=(1,2\\3,4)\n','rhs':'matrix B=(1,2\\3,4)\n','missing_cost':'matrix C=(.,2)\n','missing_lower':'matrix L=(.,0)\n','inverted':'matrix L=(6,0)\n'}[shape]
        cmd=(f'_dea_lpcontract, c(C) aub(A) bub(B) lb(L) ub(U) engine({engine})' if api=='lp' else
             f'milp_bb, c(C) aub(A) bub(B) lb(L) ub(U) integer(1 2) lpengine({engine})')
        run(f'matrix_invalid_{api}_{engine}_{shape}',s+change+cmd+f'\nassert r(status)=={3 if api=="lp" else 5}\nassert missing(r(objective))\n')
    for opt,bad in itertools.product(['lptol','inttol','maxnodes'],['0','-1','.','.a']):
        run('matrix_bad_'+opt+'_'+bad.replace('.','dot'),MAT+f'capture noisily milp_bb, c(C) aub(A) bub(B) integer(1 2) {opt}({bad})\nassert _rc==198\n')
    run('matrix_rhs_rows',MAT+'matrix B=B\'\nmilp_bb, c(C) aub(A) bub(B) integer(1 2)\nassert r(status)==0 & abs(r(objective)-6)<1e-6\n')

    rng=random.Random(1719)
    for i in range(24):
        c=[rng.randint(-5,7) for _ in range(2)]
        A=[[rng.randint(1,4) for _ in range(2)] for _ in range(2)]
        b=[rng.randint(4,17) for _ in range(2)];sense='min' if i%2 else 'max'
        feasible=[(x,y) for x in range(6) for y in range(6) if all(a[0]*x+a[1]*y<=rhs for a,rhs in zip(A,b))]
        obj=(min if sense=='min' else max)(c[0]*x+c[1]*y for x,y in feasible)
        oracles.append(dict(c=c,A=A,b=b,ub=[5,5],sense=sense,optimum=obj))
        body=f'input double x1 x2 str2 rel double rhs\n{c[0]} {c[1]} "" .\n{A[0][0]} {A[0][1]} "<=" {b[0]}\n{A[1][0]} {A[1][1]} "<=" {b[1]}\nend\n'
        body+='matrix U=(5,5)\n'
        for engine in ['simplex','ipm','auto']:
            s=body+f'lp x1 x2, {sense} intvars(x1 x2) ub(U) engine({engine})\n'
            s+=f'assert r(status)==0 & abs(r(objective)-{obj})<1e-6\nmatrix S=r(x)\n'
            s+='mata: assert(max(abs(st_matrix("S"):-round(st_matrix("S"))))<1e-7)\n'
            s+='assert S[1,1]>=0 & S[1,2]>=0 & S[1,1]<=5 & S[1,2]<=5\n'
            for row,bound in zip(A,b):s+=f'assert {row[0]}*S[1,1]+{row[1]}*S[1,2]<={bound}+1e-6\n'
            run(f'enumerated_{i}_{engine}',s)
    (out/'integer_oracles.json').write_text(json.dumps(oracles,indent=2)+'\n')

    if a.legacy_lp:
        for mode,extra in [('lp',''),('mixed','intvars(x2)'),('integer','intvars(x1 x2)')]:
            name='historical_lp_'+mode
            if a.only and not re.search(a.only,name):continue
            cmd='lp x1 x2, max '+extra
            old=TABLE+'mata: mata mlib index\n'+cmd+'\nmatrix E=r(lprslt)\nclear\nsvmat double E\nkeep E1 E2 E3\nsave expected\n'
            run(name+'_before',old,a.legacy_lp,True)
            now=TABLE+cmd+f'\nmatrix E=r(lprslt)\nclear\nsvmat double E\nkeep E1 E2 E3\ncf _all using "{out/(name+"_before")/"expected.dta"}",all\n'
            run(name,now,force=True)
    if a.legacy_dea:
        for form,rts in itertools.product(['A','B'],['vrs','crs']):
            name=f'historical_dea_{form}_{rts}'
            if a.only and not re.search(a.only,name):continue
            spec='x=y,' if form=='A' else 'y, inputs(x)'
            cmd=f'dea_lp {spec} dmu(id) rts({rts}) gen(res)'
            run(name+'_before',DEA+cmd+'\nkeep res\nsave expected\n',a.legacy_dea,True)
            s=DEA+cmd+f'\npreserve\nuse "{out/(name+"_before")/"expected.dta"}",clear\nrename res old\ntempfile ref\nsave `ref\'\nrestore\nmerge 1:1 _n using `ref\', nogen\nassert abs(res-old)<1e-6\n'
            run(name,s,force=True)

    summary=dict(total=len(rows),passed=sum(r['passed'] for r in rows),integer_oracles=len(oracles))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary));raise SystemExit(summary['total']!=summary['passed'])


if __name__=='__main__':main()
