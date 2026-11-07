#!/usr/bin/env python3
"""Positive finite effective tolerance: preflight, preservation and valid controls.

Python standard library and licensed Stata 16+. Copies the package into isolated
ado directories. Optional baseline is the pre-F16 package (core 2.0.43).
"""
from _support import package_files
import argparse,hashlib,itertools,json,re,subprocess,time
from pathlib import Path
from test_output_protection import unchanged
from test_sw2007 import DATA as ENV_DATA
REPO=Path(__file__).resolve().parents[1]
BASE='''input double(id x y q)
1 1 1 1
2 2 3 2
3 4 4 2
4 5 3 3
5 3 2 3
6 6 5 4
end
gen double w=.1
gen double p=10
gen double gx=1
gen double gy=1
gen double bad=1
gen double lo=.8*x
gen double hi=1.2*x
gen str20 xs="["+string(lo)+","+string(hi)+"]"
gen str20 note="keep original"
'''
PROFILES={
 'te':{}, 'cost':dict(model='cost',opts='w(w)'), 'revenue':dict(model='revenue',opts='p(p)'),
 'profit':dict(model='profit',opts='w(w) p(p)',rts='vrs'),
 'allocative_input':dict(model='allocative',opts='w(w)'),
 'allocative_output':dict(model='allocative',opts='p(p)',ori='output'),
 'ddf':dict(model='ddf',opts='gx(gx) gy(gy)'), 'sbm':dict(model='sbm'),
 'additive':dict(model='additive'), 'supereff':dict(model='supereff'),
 'multiplier':dict(model='multiplier'), 'supersbm':dict(model='supersbm'),
 'fdh':dict(model='fdh'), 'fdh_flag':dict(opts='fdh'),
 'interval':dict(model='imprecise',opts='xlo(lo) xhi(hi)'),
 'bounded_numeric':dict(opts='bounded(x)'),
 'bounded_string':dict(opts='bounded(xs)',x='xs'),
 'ordinal_input':dict(opts='ordinal(x)'), 'ordinal_output':dict(opts='ordinal(y)'),
 'mixed':dict(x='x q',opts='xlo(lo q) xhi(hi q) ordinal(q)'),
 'undesirable':dict(model='undesirable',opts='bad(bad)'),
 'congestion':dict(model='congestion',rts='vrs'),
 'sbmundesirable':dict(model='sbmundesirable',opts='bad(bad)'),
 'bootstrap':dict(opts='bootstrap reps1(50)'),
 'environment':dict(opts='bootstrap envars(z) reps1(50) reps2(2)',data=ENV_DATA),
}
INVALID=['0','-1','-1e-12','.','.a','.z']
QUALIFIERS=['','if id!=2','in 2/6','if id!=2 in 1/5','if 0','if runiform()>.25']

def command(profile='te',form='A',engine='auto',tol='1e-7',option='tolerance',q='',extra='',rts=None,ori=None,stage=2,explicit=True):
    d=PROFILES[profile];x=d.get('x','x');y=d.get('y','y');rts=rts or d.get('rts','crs');ori=ori or d.get('ori','input')
    head=f'dea {x}={y} {q},' if form=='A' else f'dea {y} {q}, inputs({x})'
    return f'{head} model({d.get("model","te")}) {d.get("opts","")} rts({rts}) orientation({ori}) stage({stage}) engine({engine}) gen(res) '+('dmu(id) ' if explicit else '')+(f'{option}({tol}) ' if tol is not None else '')+extra


def protected(body,cmd,rc=198,files=True):
    s=body+'''foreach v in res res_opt res_mid res_bc res_bias res_se res_lo res_hi {
 gen str20 `v'="existing result"
}
scalar _dmu_k=991
matrix _lam_val=(1,2\\3,4)
local rng_before "`c(rngstate)'"
'''
    if files:s+='save "target.dta"\ncopy "target.dta" "target_before.dta"\n'
    s+=unchanged(cmd,rc)
    s+='assert "`rng_before\'"=="`c(rngstate)\'"\nassert scalar(_dmu_k)==991\nmata: assert(all(st_matrix("_lam_val"):==(1,2\\3,4)))\n'
    s+='capture confirm variable _dmu,exact\nassert _rc==111\n'
    return s


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',required=True,type=Path);ap.add_argument('--baseline-ado',type=Path);ap.add_argument('--only');a=ap.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'current_ado';pkg.mkdir();hashes={}
    for n in package_files(REPO):
        raw=(REPO/n).read_bytes();(pkg/Path(n).name).write_bytes(raw);hashes[n]=hashlib.sha256(raw).hexdigest()
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,baseline=str(a.baseline_ado),stata=a.stata),indent=2));results=[]
    def run(name,s,source=pkg,message=False,force=False):
        if not force and a.only and not re.search(a.only,name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 120\n'
        for f in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/f).mkdir();head+=f'sysdir set {f} "{d/f}"\n'
        head+=f'adopath ++ "{source.resolve()}"\n';(d/'test.do').write_text(head+s+'\ndi "TEST_PASS"\n');row=dict(test=name,passed=False);start=time.monotonic()
        try:
            proc=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=300)
            if proc.returncode or not (d/'test.log').exists():(d/'process.json').write_text(json.dumps(dict(rc=proc.returncode,stderr=proc.stderr.decode(errors='replace')),indent=2))
            log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else proc.stderr.decode(errors='replace')
            assert proc.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2100:]
            if message:assert 'tolerance() must be positive and finite' in log
            if (d/'target_before.dta').exists():assert (d/'target.dta').read_bytes()==(d/'target_before.dta').read_bytes()
            row['passed']=True
        except Exception as e:row['error']=str(e)
        row['seconds']=round(time.monotonic()-start,3);results.append(row);(out/'results.json').write_text(json.dumps(results,indent=2));print(json.dumps(row),flush=True)
    for i,(profile,form,engine,bad) in enumerate(itertools.product(PROFILES,['A','B'],['simplex','ipm','auto'],INVALID)):
        body=PROFILES[profile].get('data',BASE);q=QUALIFIERS[i%len(QUALIFIERS)]
        cmd=command(profile,form,engine,bad,q=q,explicit=False,stage=1+i%2,extra='replacevars saving("target.dta") replace')
        run(f'reject_{profile}_{form}_{engine}_{bad.replace(".","dot")}',protected(body,cmd),message=True)
    for form,value in itertools.product(['A','B'],['-0','1e309','-1e309','1e-400','abc','1/0']):
        run(f'parser_{form}_{value.replace("/","_")}',protected(BASE,command(form=form,tol=value,explicit=False,extra='replacevars saving("target.dta") replace')))
    for option,bad in itertools.product(['tol','tole','toler','tolera','toleran','toleranc','tolerance'],INVALID):
        run(f'prefix_{option}_{bad.replace(".","dot")}',protected(BASE,command(form='B',option=option,tol=bad,extra='replacevars saving("target.dta") replace')),message=True)
    for engine,bad in itertools.product(['simplex','ipm','auto'],INVALID):
        run(f'legacy_active_{engine}_{bad.replace(".","dot")}',protected(BASE,command(engine=engine,tol=bad,option='tol1',extra='replacevars saving("target.dta") replace')),message=True)
    # Public forwarding commands reach the same core validator.
    for wrapper,form,bad in itertools.product(['dea_congestion','dea_sbmundesirable'],['A','B'],INVALID):
        head=f'{wrapper} x=y,' if form=='A' else f'{wrapper}, inputs(x) outputs(y)'
        cmd=head+(' bad(bad)' if wrapper=='dea_sbmundesirable' else '')+f' tolerance({bad}) gen(res) replacevars'
        run(f'wrapper_{wrapper}_{form}_{bad.replace(".","dot")}',protected(BASE,cmd,files=False),message=True)
    # The contract helpers reject missing values as well as nonpositive values.
    matrices='matrix C=(1,2)\nmatrix AE=(1,1)\nmatrix BE=(3)\nmatrix LB=(0,0)\nmatrix UB=(2,3)\n'
    for helper,engine,bad in itertools.product(['lp','te'],['simplex','ipm','auto'],INVALID):
        cmd=(f'_dea_lpcontract, c(C) aeq(AE) beq(BE) lb(LB) ub(UB) engine({engine}) tolerance({bad})' if helper=='lp' else f'_dea_tecontract, inputs(x) outputs(y) gen(res) engine({engine}) tolerance({bad})')
        s=protected(BASE+matrices,cmd,files=False)+'mata: assert(all(st_matrix("C"):==(1,2)))\nmata: assert(all(st_matrix("UB"):==(2,3)))\n'
        run(f'contract_bad_{helper}_{engine}_{bad.replace(".","dot")}',s,message=True)
    # Numeric parser acceptance has no artificial upper/lower bound on positives.
    for form,value in itertools.product(['A','B'],['1e-300','1e-14','1e-10','1e-7','.125','1','2','1e300']):
        parser='_dea_parse_equal x=y,' if form=='A' else '_dea_parse_standard y, inputs(x)'
        s=BASE+'quietly run "'+str(pkg/'dea.ado')+'"\n'+parser+f' tolerance({value})\nassert !missing(r(tolerance)) & r(tolerance)>0\nassert r(tolerance)=={value}\n'
        run(f'positive_parse_{form}_{value}',s)
    legacy=[('', '1e-7'), ('tol1(1e-14)', '1e-14'), ('tol1(1e-6)', '1e-6'), ('tolerance(1e-7)', '1e-7')]
    for i,(options,value) in enumerate(legacy):
        s=BASE+'quietly run "'+str(pkg/'dea.ado')+'"\n_dea_parse_equal x=y, '+options+f'\nassert r(tolerance)=={value}\n'
        run(f'legacy_precedence_{i}',s)
    for i,options in enumerate(['tolerance(1e-7) tol1(1e-6)', 'tol1(1e-7) tolerance(1e-7)', 'tolerance(1e-8) tol1(0)', 'tol2(0)', 'tol2(.)', 'tol2(1e-8)', 'tol1(1e-6) tol2(1e-8)']):
        run(f'legacy_reject_{i}',protected(BASE,command(tol=None,explicit=False,extra=options+' replacevars saving("target.dta") replace')))
    # Positive changes do not change solver outputs. Same seed for bootstrap.
    def historical(name,setup,cmd):
        if a.only and not re.search(a.only,name):return
        control=a.baseline_ado or pkg
        old=setup+'set seed 883\n'+cmd+'\nscalar returned_tolerance=r(tolerance)\nkeep res*\nsave "expected.dta"\n'
        run(name+'_before',old,control,force=True)
        now=setup+'set seed 883\n'+cmd+f'\nassert !missing(r(tolerance)) & r(tolerance)>0\ncf res* using "{out/(name+"_before")/"expected.dta"}",all\n'
        run(name,now,force=True)
    for i,(engine,rts,ori,stage,form,tol) in enumerate(itertools.product(['simplex','ipm','auto'],['crs','vrs','nirs','ndrs'],['input','output'],[1,2],['A','B'],['1e-4','1e-7','1e-10'])):
        cmd=command(form=form,engine=engine,rts=rts,ori=ori,stage=stage,tol=tol,q=['','if id!=2','in 2/6','if id!=2 in 1/5'][i%4],extra='reference')
        historical(f'healthy_radial_{engine}_{rts}_{ori}_{stage}_{form}_{tol}',BASE,cmd)
    for i,(profile,form,tol) in enumerate(itertools.product(PROFILES,['A','B'],[None,'1e-8'])):
        historical(f'healthy_model_{profile}_{form}_{tol or "default"}',PROFILES[profile].get('data',BASE),command(profile,form,['simplex','ipm','auto'][i%3],tol=tol))
    analytic='''input double x y
1 1
2 3
4 4
5 3
end
'''
    for engine,tol,scaling in itertools.product(['simplex','ipm','auto'],['1e-4','1e-7','1e-10'],['none','row']):
        s=matrices+f'_dea_lpcontract, c(C) aeq(AE) beq(BE) lb(LB) ub(UB) engine({engine}) tolerance({tol}) scaling({scaling})\n'
        s+=f'assert r(status)==0\nassert abs(r(objective)-4)<max(1e-7,10*{tol})\nmatrix S=r(x)\nmata: assert(max(abs(st_matrix("S"):- (2,1)))<max((1e-7,10*{tol})))\n'
        run(f'analytic_lp_{engine}_{tol}_{scaling}',s)
    for engine,tol,stage,ori in itertools.product(['simplex','ipm','auto'],['1e-4','1e-7','1e-10'],[1,2],['input','output']):
        name=f'analytic_te_{engine}_{tol}_{stage}_{ori}'
        s=analytic+f'_dea_tecontract, inputs(x) outputs(y) gen(res) status(status) rts(crs) orientation({ori}) stage({stage}) engine({engine}) tolerance({tol})\n'
        # Bounded native retries recover the former loose/tight stage-2
        # failures. Exact analytic efficiencies still apply to all four DMUs.
        s+=f'assert r(n_failed)==0 & r(stage2_failures)==0\nassert abs(res-y/(1.5*x))<max(1e-7,10*{tol})\nassert status==0\n'
        if engine=='ipm' and stage==2 and (tol=='1e-4' or (tol=='1e-10' and ori=='input')):
            s+='assert r(native_retries)>0 & r(status_conflicts)==0\n'
        run(name,s)
    if a.baseline_ado:
        for form,value in itertools.product(['A','B'],['0','-1','.','.z']):
            run(f'old_core_{form}_{value.replace(".","dot")}',BASE+command(form=form,engine='simplex',tol=value)+'\nassert r(tolerance)=='+value+'\n',a.baseline_ado)
        for helper,value in itertools.product(['lp','te'],['.','.z']):
            cmd=(f'_dea_lpcontract, c(C) aeq(AE) beq(BE) lb(LB) ub(UB) engine(ipm) tolerance({value})' if helper=='lp' else f'_dea_tecontract, inputs(x) outputs(y) gen(res) engine(ipm) tolerance({value})')
            run(f'old_contract_{helper}_{value.replace(".","dot")}',BASE+matrices+cmd+'\n',a.baseline_ado)
    summary=dict(total=len(results),passed=sum(r['passed'] for r in results));(out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True);raise SystemExit(summary['total']!=summary['passed'])

if __name__=='__main__':main()
