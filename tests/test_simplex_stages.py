#!/usr/bin/env python3
"""Stage routing, independent radial/slack LPs, reference feasibility and failures.

Requires NumPy/SciPy and Stata 16+. Counters and faults exist only in private
copies. --baseline-ado must be the pre-F12 package (dea 2.0.40).
"""
from _support import package_files
import argparse
import csv
import difflib
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import time

import numpy as np
import scipy
from scipy.optimize import linprog
from verify_bootstrap_lp import solve
from test_output_protection import unchanged

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
gen str20 sentinel="preserve me"
'''
RANDOM = '''set seed 7351
set obs 12
gen long id=_n
gen double x1=1+10*runiform()
gen double x2=1+10*runiform()
gen double y1=1+5*runiform()
gen double y2=1+5*runiform()
gen str20 sentinel="preserve me"
'''
FIXTURES = {'slack': DATA, 'random': RANDOM, 'scaled': RANDOM+'''replace x1=x1*1e7
replace x2=x2*1e-4
replace y1=y1*1e-3
replace y2=y2*1e6
''', 'duplicate': DATA+'expand 2 if id==2\nreplace id=_n\n'}
RTS = {'crs':'crs', 'ccr':'crs', 'vrs':'vrs', 'bcc':'vrs', 'nirs':'drs', 'drs':'drs', 'ndrs':'irs', 'irs':'irs'}
INIT = '''scalar audit_s1=0
scalar audit_s2=0
scalar audit_dual=0
scalar fault_primal=0
scalar fault_dual=0
scalar fault_stage2=0
'''


def instrument(source, dest):
    shutil.copytree(source, dest)
    p=dest/'_dea_simplex.ado'; old=p.read_text(); s=old
    a=s.index('real colvector _dea_solve_dmu('); b=s.index('// BUILD LP: STAGE 1',a)
    part=s[a:b]
    for target, code, count in [
        ('res1 = _rsm_two_phase(prob1, tol, trace)', 'st_numscalar("audit_s1",st_numscalar("audit_s1")+1);\n    if (st_numscalar("fault_primal")) res1.status=3;',1),
        ('res_m = _rsm_two_phase(prob_m, tol, trace)', 'st_numscalar("audit_dual",st_numscalar("audit_dual")+1);\n        if (st_numscalar("fault_dual")) res_m.status=3;',1),
        ('res2 = _rsm_two_phase(prob2, tol, trace)', 'st_numscalar("audit_s2",st_numscalar("audit_s2")+1);\n    if (st_numscalar("fault_stage2")) res2.status=3;',2),
    ]:
        targets=[target]
        if target=='res2 = _rsm_two_phase(prob2, tol, trace)':
            targets.append('res2 = _rsm_two_phase(prob2, linktol, trace)')
        assert sum(part.count(t) for t in targets)==count,target
        for t in targets: part=part.replace(t,t+'\n    '+code)
    s=s[:a]+part+s[b:];p.write_text(s)
    (dest.parent/(dest.name+'.patch')).write_text(''.join(difflib.unified_diff(old.splitlines(True),s.splitlines(True))))


def command(engine, rts, ori, stage, sample='', extra='', form='A', model='te'):
    io = f'x1 x2 = y1 y2 {sample},' if form=='A' else f'y1 y2 {sample}, inputs(x1 x2)'
    return f'dea {io} dmu(id) engine({engine}) rts({rts}) ort({ori}) model({model}) stage({stage}) gen(e{stage}) {extra}'


def capture(stage, expected=None):
    s=f'''assert r(stage)=={stage}
local engine{stage} "`r(engine)'"
matrix L{stage}=r(lambda)
assert r(n_radial_failed)==0
assert r(n_failed)==r(n_stage2_failed)
matrix SS{stage}=r(solution_source)
matrix ST{stage}=r(stage_status)
mata: assert(all(st_matrix("SS{stage}")[.,1]:==1) | "`engine{stage}'"=="ipm")
mata: assert(all(st_matrix("SS{stage}")[.,2]:==({stage}:*(st_matrix("ST{stage}")[.,{1 if stage==1 else 2}]:==0))))
'''
    if stage==1:
        s+='capture matrix list r(slack_in)\nassert _rc==111\ncapture matrix list r(slack_out)\nassert _rc==111\n'
    else:
        s+='''matrix SI2=r(slack_in)
matrix SO2=r(slack_out)
matrix ST2=r(stage_status)
'''
    if expected:
        second_calls="0" if stage==1 else "`n'"
        s+=f'''if "`engine{stage}'"=="simplex" {{
 assert audit_s1=={expected} & audit_s2=={second_calls} & audit_dual==0
}}
else assert audit_s1==0 & audit_s2==0
'''
    s+=f'''local refs{stage} "`r(refvars)'"
'''
    return s


def export_pair(sample):
    return '''preserve
keep if selected
foreach s in 1 2 {
 svmat double L`s', names(l`s'_)
}
svmat double SI2, names(si2_)
svmat double SO2, names(so2_)
svmat double ST2, names(stage2_)
local statuscols "stage2_*"
gen byte ipm1=("`engine1'"=="ipm")
export delimited id x1 x2 y1 y2 e1 e2 ipm1 l1_* l2_* si2_* so2_* `statuscols' using "scores.csv", replace
restore
assert sentinel=="preserve me"
'''


def read(path):
    with path.open() as f: rows=list(csv.DictReader(f))
    return {k:np.array([float(r[k]) if r[k] not in ('.', '') else np.nan for r in rows]) for k in rows[0]}


def oracle(d, ori, rts):
    v=read(d/'scores.csv');x=np.column_stack([v['x1'],v['x2']]);y=np.column_stack([v['y1'],v['y2']]);n=len(x)
    xs=np.max(abs(x),axis=0);ys=np.max(abs(y),axis=0);x=x/xs;y=y/ys
    gap=primal=slackgap=0.;lpdiff=0.;count=0
    for k in range(n):
        eff=solve(x[k],y[k],x,y,ori,RTS[rts]);count+=1
        radial=eff if ori=='input' else 1/eff
        bx=radial*x[k] if ori=='input' else x[k]
        by=y[k] if ori=='input' else radial*y[k]
        # Independent stage-2 LP; unit-normalized sum matches documented scaling.
        ae=np.block([[x.T,np.eye(2),np.zeros((2,2))],[y.T,np.zeros((2,2)),-np.eye(2)]])
        be=np.r_[bx,by];au=bu=None;row=np.r_[np.ones(n),np.zeros(4)]
        if RTS[rts]=='vrs': ae=np.vstack([ae,row]);be=np.r_[be,1]
        if RTS[rts] in ['drs','irs']:
            sign=1 if RTS[rts]=='drs' else -1;au=(sign*row)[None,:];bu=np.array([sign])
        fit=linprog(np.r_[np.zeros(n),-np.ones(4)],A_eq=ae,b_eq=be,A_ub=au,b_ub=bu,bounds=(0,None),method='highs',options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
        assert fit.success,fit.message;count+=1
        first_slack=np.nan
        for stage in [1,2]:
            actual=v[f'e{stage}'][k];gap=max(gap,abs(actual-eff));assert np.isfinite(actual) and abs(actual-eff)<1e-6
            lam=np.array([v[f'l{stage}_{j}'][k] for j in range(1,n+1)])
            if stage==2 and 'stage2_2' in v and v['stage2_2'][k]!=0:
                assert np.isnan(lam).all()
                assert all(np.isnan(v[f'{side}2_{j}'][k]) for side in ['si','so'] for j in [1,2])
                continue
            if stage==1 and v.get('ipm1',np.zeros(n))[k]==1 and np.isnan(lam).all():
                # Checked dual recovery has a score but no primal weights.
                continue
            assert np.isfinite(lam).all() and min(lam)>-1e-7
            radial_actual=actual if ori=='input' else 1/actual
            si=(radial_actual*x[k] if ori=='input' else x[k])-lam@x
            so=lam@y-(y[k] if ori=='input' else radial_actual*y[k])
            residual=max(0.,float(-min(si)),float(-min(so)))
            if RTS[rts]=='vrs':residual=max(residual,abs(sum(lam)-1))
            if RTS[rts]=='drs':residual=max(residual,sum(lam)-1)
            if RTS[rts]=='irs':residual=max(residual,1-sum(lam))
            primal=max(primal,residual);assert residual<1e-6,(stage,k,residual)
            if stage==1:first_slack=sum(si)+sum(so)
            if stage==2 and 'si2_1' in v:
                observed=np.r_[[v[f'si2_{j}'][k]/xs[j-1] for j in [1,2]],[v[f'so2_{j}'][k]/ys[j-1] for j in [1,2]]]
                assert np.isfinite(observed).all() and np.max(abs(observed-np.r_[si,so]))<1e-6
                delta=abs(sum(observed)+fit.fun);slackgap=max(slackgap,delta);assert delta<1e-6,(k,delta)
                if np.isfinite(first_slack): lpdiff=max(lpdiff,sum(observed)-first_slack)
    if d.name=='grid_simplex_crs_input_slack': assert lpdiff>.1
    return dict(lp_count=count,max_score_gap=float(gap),max_primal_violation=float(primal),max_slack_gap=float(slackgap),max_slack_improvement=float(lpdiff))


def fallback_oracle(d, ori, rts, missing):
    v=read(d/'fallback.csv');x=np.column_stack([v['x1'],v['x2']]);y=np.column_stack([v['y1'],v['y2']]);gap=0.;count=0
    score=v['e1'] if 'e1' in v else v['e2']
    if missing:
        assert np.isnan(score).all()
    else:
        for k in range(len(x)):
            expected=solve(x[k],y[k],x,y,ori,rts);count+=1
            assert np.isfinite(score[k]);gap=max(gap,abs(score[k]-expected))
        assert gap<1e-6,gap
    return dict(lp_count=count,max_score_gap=float(gap))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stata',required=True);p.add_argument('--output',required=True,type=Path);p.add_argument('--baseline-ado',type=Path);p.add_argument('--only');args=p.parse_args()
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'current_ado';pkg.mkdir();hashes={}
    for name in package_files(REPO):
        data=(REPO/name).read_bytes();(pkg/Path(name).name).write_bytes(data);hashes[name]=hashlib.sha256(data).hexdigest()
    audit=out/'instrumented_ado';instrument(pkg,audit)
    before=out/'instrumented_before'
    if args.baseline_ado: instrument(args.baseline_ado.resolve(),before)
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,numpy=np.__version__,scipy=scipy.__version__,baseline=str(args.baseline_ado)),indent=2))
    rows=[]
    def run(name,script,source=pkg,check=None):
        if args.only and not re.search(args.only,name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 120\n'
        for f in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/f).mkdir();head+=f'sysdir set {f} "{d/f}"\n'
        head+=f'adopath ++ "{source.resolve()}"\n'+INIT
        (d/'test.do').write_text(head+script+'\ndi "TEST_PASS"\n');start=time.monotonic();row=dict(test=name,passed=False)
        try:
            proc=subprocess.run([args.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=240)
            log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else proc.stderr.decode(errors='replace')
            assert 'TEST_PASS' in log.splitlines(),log[-2300:]
            if check:row.update(check(d))
            row['passed']=True
        except Exception as e:row['error']=str(e)
        row['seconds']=round(time.monotonic()-start,3);rows.append(row);(out/'results.json').write_text(json.dumps(rows,indent=2));print(json.dumps(row),flush=True)

    for engine in ['simplex','auto','ipm']:
        for rts in RTS:
            for ori in ['input','output']:
                for fixture, data in FIXTURES.items():
                    sample='if id>=2 in 1/8' if fixture=='duplicate' else ''
                    s=data+'gen byte selected=0\nreplace selected=1 '+sample+'\nquietly count if selected\nlocal n=r(N)\n'
                    for stage in [1,2]:
                        s+=INIT+command(engine,rts,ori,stage,sample,extra='reference saveref replacevars',form='A' if stage==1 else 'B')+'\n'
                        expected="4*`n'" if RTS[rts]=='vrs' else "`n'"
                        s+=capture(stage,expected)
                    # Both scores come from the same radial LP with simplex.
                    s+='if "`engine1\'"=="simplex" & "`engine2\'"=="simplex" assert e1==e2\n'
                    # Copy selected rows and preserve unavailable IPM projections.
                    s+=export_pair('ipm' if engine=='ipm' or (engine=='auto' and fixture=='scaled') else '')
                    run(f'grid_{engine}_{rts}_{ori}_{fixture}',s,source=audit,check=lambda d,ori=ori,rts=rts:oracle(d,ori,rts))

    if args.baseline_ado:
        for rts in ['crs','vrs','drs','irs']:
            for ori in ['input','output']:
                for stage in [1,2]:
                    key=f'history_{rts}_{ori}_{stage}'
                    script=DATA+command('simplex',rts,ori,stage,extra='reference saveref')+'\n'
                    old_calls=8*(4 if rts=='vrs' else 1)
                    run(key+'_before',script+f'assert audit_s2=={old_calls}\nkeep e{stage} _ref_*\nsave "old.dta",replace\n',before)
                    # Reference weights are allowed to change for stage 1; stage 2 must remain exact.
                    variables=f'e{stage}'+(' _ref_*' if stage==2 else '')
                    run(key,script+f'cf {variables} using "{out/(key+"_before")/"old.dta"}",all\n')

    for engine in ['simplex','auto']:
        for rts in ['crs','vrs','drs','irs']:
            for ori in ['input','output']:
                for stage in [1,2]:
                    for mode,pf,df,sf in [('primal',1,0,0),('primal_stage2',1,0,1),('stage2',0,0,1),('both',1,1,0)]:
                        s=DATA+f'scalar fault_primal={pf}\nscalar fault_dual={df}\nscalar fault_stage2={sf}\n'
                        s+=command(engine,rts,ori,stage,extra='reference saveref saving("saved.dta")')+'\nlocal refs "`r(refvars)\'"\nmatrix L=r(lambda)\n'
                        nlp=8*(4 if rts=='vrs' and not(mode=='both' and engine=='auto') else 1)
                        expected2=8 if stage==2 and mode!='both' else 0
                        s+=f'assert audit_s1=={nlp} & audit_s2=={expected2} & audit_dual=={nlp*pf}\n'
                        recovered=(mode=='stage2' and stage==1) or (mode=='primal' and stage==2) or (mode=='both' and engine=='auto')
                        if stage==1:
                            s+='capture matrix list r(slack_in)\nassert _rc==111\ncapture matrix list r(slack_out)\nassert _rc==111\n'
                        elif not(mode=='both' and engine=='auto'):
                            s+='matrix SI=r(slack_in)\nmatrix SO=r(slack_out)\n'
                            if mode in ['primal_stage2','both','stage2']:
                                s+='mata: assert(all(st_matrix("SI"):>=.) & all(st_matrix("SO"):>=.))\n'
                            else:
                                s+='mata: assert(!hasmissing(st_matrix("SI")) & !hasmissing(st_matrix("SO")))\n'
                        if mode=='both' and engine=='simplex':s+=f'assert missing(e{stage})\nassert r(n_failed)==8\n'
                        elif mode=='both' and engine=='auto':
                            s+=f'assert !missing(e{stage})\nassert r(n_radial_failed)==0\nassert r(n_failed)==r(n_stage2_failed)\n'
                            if stage==2:
                                s+='matrix ST=r(stage_status)\nmata: assert(st_numscalar("r(n_failed)")==sum(st_matrix("ST")[.,2]:!=0))\n'
                        else:s+=f'assert !missing(e{stage})\nassert r(n_failed)=={8 if stage==2 and mode in ["primal_stage2","stage2"] else 0}\n'
                        s+='matrix ST=r(stage_status)\nmatrix SS=r(solution_source)\n'
                        if not(mode=='both' and engine=='auto'):
                            source=0 if mode=='both' else (2 if pf else 1)
                            refstage=stage if recovered else 0
                            code=3 if stage==2 and sf and mode!='both' else (0 if stage==2 and mode=='primal' else '.')
                            s+=f'mata: assert(all(st_matrix("SS")[.,1]:=={source}) & all(st_matrix("SS")[.,2]:=={refstage}))\n'
                            s+=f'mata: assert(all(st_matrix("ST")[.,1]:=={int(mode=="both")}) & all(st_matrix("ST")[.,2]:=={code}))\n'
                        if recovered:
                            if mode!='both' or engine!='auto': s+='mata: assert(!hasmissing(st_matrix("L")))\n'
                            s+='local final_engine "`r(engine)\'"\n'
                            # Validate fallback references against the actual scores and data.
                            s+=f'''mata:
X=st_data(.,tokens("x1 x2"));Y=st_data(.,tokens("y1 y2"));E=st_data(.,"e{stage}");L=st_matrix("L")
if (st_global("r(engine)")=="ipm") {{
 for (k=1;k<=rows(L);k++) assert(rowmissing(L[k,.])==0 | rowmissing(L[k,.])==cols(L))
 good=selectindex(rowmissing(L):==0)
 Xe=X[good,.];Ye=Y[good,.];E=E[good,.];L=L[good,.]
}}
else {{ Xe=X;Ye=Y; }}
assert(rows(L)>0 & min(L)>=-1e-7)
assert(max(L*X - {'Xe:*E' if ori=='input' else 'Xe'})<1e-5)
assert(max({'Ye' if ori=='input' else 'Ye:/E'}-L*Y)<1e-5)
end
'''
                        else:s+='mata: assert(all(st_matrix("L"):>=.))\n'
                        s+='capture confirm variable _dea_radial_status\nassert _rc==111\npreserve\nuse "saved.dta",clear\nmata: assert(all(st_data(.,tokens(st_local("refs"))):==st_matrix("L")))\nmata: assert(all(st_data(.,tokens("_dea_radial_status _dea_stage2_status _dea_score_source _dea_ref_stage")):==(st_matrix("ST"),st_matrix("SS"))))\nrestore\n'
                        s+=f'export delimited id x1 x2 y1 y2 e{stage} using "fallback.csv",replace\n'
                        run(f'fault_{engine}_{rts}_{ori}_{stage}_{mode}',s,audit,
                            check=lambda d,ori=ori,rts=rts,missing=(mode=='both' and engine=='simplex'):fallback_oracle(d,ori,rts,missing))

    for stage in [1,2]:
        for j,(model,extra) in enumerate([('multiplier',''),('multiplier','saveweights'),('te','ordinal(x1)'),('te','xlo(lo) ordinal(x2)')]):
            s=DATA+'gen double lo=.9*x1\n'
            # For mixed lists the ordinal slot is an explicit placeholder.
            extra=extra.replace('xlo(lo)','xlo(lo x2)')
            s+=command('simplex','crs','input',stage,extra=extra,model=model)+'\n'
            if stage==1:s+='assert audit_s2==0\n'
            elif model=='te' or extra!='saveweights':s+='assert audit_s2>0\n'
            run(f'routing_{stage}_{j}',s,audit)
    for rts in ['crs','vrs','drs','irs']:
        for ori in ['input','output']:
            s=DATA+'gen byte selected=1\nfindfile _dea_simplex.ado\nquietly run "`r(fn)\'"\n'
            for i,tail in enumerate(['',' , ""', ', "", 1', ', "", 2']):
                s+=f'''scalar audit_s2=0
mata: _dea_simplex_full("x1 x2","y1 y2","selected","{ori}","{rts}",1e-7{tail})
matrix R{i}=r(dea_results)
assert audit_s2=={0 if i==2 else 8}
'''
            s+='mata: assert(all(st_matrix("R0"):==st_matrix("R1")) & all(st_matrix("R0"):==st_matrix("R3")))\n'
            for i,tail in enumerate(['', ', ""', ', "", ""', ', "", "", 1', ', "", "", 2']):
                s+=f'''gen double e{i}=.
scalar audit_s2=0
mata: _dea_simplex_envelopment("x1 x2","y1 y2","e{i}","selected","{ori}","{rts}",1e-7,""{tail})
assert audit_s2=={0 if i==3 else 8}
assert e{i}==e0
'''
            run(f'optional_mata_{rts}_{ori}',s,audit)
    for stage in [-1,0,3]:run(f'invalid_stage_{stage}',DATA+unchanged(command('simplex','crs','input',stage),198))
    summary=dict(total=len(rows),passed=sum(r['passed'] for r in rows),lp_count=sum(r.get('lp_count',0) for r in rows))
    for metric in ['max_score_gap','max_primal_violation','max_slack_gap','max_slack_improvement']:summary[metric]=max((r.get(metric,0) for r in rows),default=0)
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
    raise SystemExit(summary['total']!=summary['passed'])

if __name__=='__main__':main()
