#!/usr/bin/env python3
"""IPM TE stage-2 technology, projection identity, status, scope and old scores.

Requires NumPy/SciPy and licensed Stata 16+. Faults/counters exist only in copied
engines. Native numerical failures are checked for truthful incomplete results,
not counted as successful LP solves. --baseline-ado enables a strict historical
comparison; pre-2.0.47 snapshots intentionally differ after the tolerance change.
"""
from _support import package_files
import argparse, csv, hashlib, itertools, json, re, shutil, subprocess, time
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from test_simplex_stages import DATA, RANDOM, RTS, FIXTURES as OLD_FIXTURES
from verify_bootstrap_lp import solve
from test_saveref import instrument as retry_instrument

REPO=Path(__file__).resolve().parents[1]
FIXTURES=dict(OLD_FIXTURES)
FIXTURES.update({
    'p03_crs':'input double id x1 x2 y1 y2\n1 1 1 1 1\n2 2 1 2 2\n3 3 6 1 1\nend\ngen str20 sentinel="preserve me"\n',
    'p03_vrs':'input double id x1 x2 y1 y2\n1 1 1 1 1\n2 1 .5 1 1\n3 3 6 1 1\nend\ngen str20 sentinel="preserve me"\n',
    'p03_irs':'input double id x1 x2 y1 y2\n1 1 1 1 1\n2 2 1 2 2\n3 3 6 1 1\n4 4 2 4 4\n5 .5 .25 .2 .2\n6 1 1 1 1\nend\ngen str20 sentinel="preserve me"\n',
})
INIT='scalar audit_s1=0\nscalar audit_s2=0\nscalar fault_primal=0\nscalar fault_dual=0\nscalar fault_stage2=0\n'

def instrument(src,dst):
    shutil.copytree(src,dst);p=dst/'_dea_ipm.ado';s=p.read_text()
    start=s.index('void _dea_ipm_envelopment(');end=s.index('real scalar _dea_ipm_slacks(',start)
    a=s[start:end]
    anchor='rc = S.optimize()';assert a.count(anchor)==2
    a=a.replace(anchor,anchor+'\n            st_numscalar("audit_s1",st_numscalar("audit_s1")+1)')
    anchor='if (S.converged() && S.errorcode() == 0)';assert a.count(anchor)==2
    a=a.replace(anchor,'if (!st_numscalar("fault_primal") && S.converged() && S.errorcode() == 0)')
    # Persistent primal faults deliberately disable its new recovery path.
    a=a.replace('primal_status=_dea_ipm_finish(', 'if (!st_numscalar("fault_primal")) primal_status=_dea_ipm_finish(')
    anchor='lam_k = J(1, n, .)  // lambda unavailable from dual';assert a.count(anchor)==1
    a=a.replace(anchor,'if (st_numscalar("fault_dual")) eff_score=.\n            '+anchor)
    s=s[:start]+a+s[end:]
    start=s.index('real scalar _dea_ipm_slacks(');end=s.index('real scalar _dea_ipm_dual_input(',start)
    a=s[start:end];anchor='rc = S2.optimize()';assert a.count(anchor)==1
    a=a.replace(anchor,anchor+'\n    st_numscalar("audit_s2",st_numscalar("audit_s2")+1)')
    anchor='if (S2.converged() && S2.errorcode() == 0)';assert a.count(anchor)==1
    a=a.replace(anchor,'if (k!=st_numscalar("fault_stage2") && S2.converged() && S2.errorcode() == 0)')
    # Persistent projection faults must not be repaired by a healthy retry.
    a=a.replace('status=_dea_ipm_finish(', 'if (k!=st_numscalar("fault_stage2")) status=_dea_ipm_finish(')
    s=s[:start]+a+s[end:];p.write_text(s)

def read(path):
    with path.open() as f:rows=list(csv.DictReader(f))
    return {k:np.array([float(r[k]) if r[k] not in ('.','') else np.nan for r in rows]) for k in rows[0]}

def oracle(d,rts,ori,tol=1e-7):
    a=read(d/'scores.csv');n=len(a['id']);x=np.c_[a['x1'],a['x2']];y=np.c_[a['y1'],a['y2']]
    xs=np.max(abs(x),axis=0);ys=np.max(abs(y),axis=0);x=x/xs;y=y/ys
    maxrad=maxeq=maxslack=maxscaled=0.;success=failed=indeterminate=0
    for k in range(n):
        e=solve(x[k],y[k],x,y,ori,RTS[rts]);maxrad=max(maxrad,abs(e-a['e2'][k]))
        assert np.isfinite(a['e2'][k]) and abs(e-a['e2'][k])<1e-6
        assert a['e1'][k]==a['e2'][k]
        lam=np.array([a[f'l2_{j}'][k] for j in range(1,n+1)])
        si=np.array([a[f'si_{j}'][k]/xs[j-1] for j in [1,2]])
        so=np.array([a[f'so_{j}'][k]/ys[j-1] for j in [1,2]])
        if 'status_2' in a and a['status_2'][k]!=0:
            assert np.isnan(lam).all() and np.isnan(si).all() and np.isnan(so).all();failed+=1;continue
        assert np.isfinite(np.r_[lam,si,so]).all()
        radial=a['e2'][k] if ori=='input' else 1/a['e2'][k]
        bx=radial*x[k] if ori=='input' else x[k];by=y[k] if ori=='input' else radial*y[k]
        residual=max(np.max(abs(lam@x+si-bx)),np.max(abs(lam@y-so-by)),max(0.,-np.min(np.r_[lam,si,so])))
        if RTS[rts]=='vrs':residual=max(residual,abs(sum(lam)-1))
        if RTS[rts]=='drs':residual=max(residual,sum(lam)-1)
        if RTS[rts]=='irs':residual=max(residual,1-sum(lam))
        maxeq=max(maxeq,residual);assert residual<1e-6,(k,residual)
        ae=np.block([[x.T,np.eye(2),np.zeros((2,2))],[y.T,np.zeros((2,2)),-np.eye(2)]])
        be=np.r_[bx,by];row=np.r_[np.ones(n),np.zeros(4)];au=bu=None
        if RTS[rts]=='vrs':ae=np.vstack([ae,row]);be=np.r_[be,1]
        if RTS[rts] in ['drs','irs']:
            sign=1 if RTS[rts]=='drs' else -1;au=(sign*row)[None,:];bu=np.array([sign])
        if 'status_2' in a:
            from test_lp_residuals import residual as primal_residual
            z=np.r_[lam,si,so]
            diag=primal_residual(z,ae,be,au if au is not None else np.empty((0,len(z))),
                bu if bu is not None else np.empty(0),np.zeros(len(z)),np.full(len(z),np.inf))
            maxscaled=max(maxscaled,max(diag[3:]))
            assert max(diag[3:])<=float(tol)+4e-16,(k,tol,diag)
        f=linprog(np.r_[np.zeros(n),-np.ones(4)],A_eq=ae,b_eq=be,A_ub=au,b_ub=bu,bounds=(0,None),method='highs-ds',
                  options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
        if f.success:
            gap=abs(sum(si)+sum(so)+f.fun);maxslack=max(maxslack,gap);assert gap<1e-6,(k,gap)
        else:
            # An approximate radial point may be strictly outside the face.
            # This is a remaining numerical-policy issue, not an oracle pass.
            indeterminate+=1
        success+=1
    return dict(dmus=n,projection_successes=success,projection_failures=failed,oracle_nonoptimal=indeterminate,
                max_radial_gap=float(maxrad),max_projection_residual=float(maxeq),max_scaled_projection_residual=float(maxscaled),max_slack_gap=float(maxslack))

def script(data,engine,rts,ori,tol,old=False,qual='',extra=''):
    s=data+('gen byte selected=inrange(_n,2,8) & id>=2\n' if qual else 'gen byte selected=1\n')
    s+='quietly count if selected\nlocal n=r(N)\n'
    for stage in [1,2]:
        io=f'x1 x2=y1 y2 {qual},' if stage==1 else f'y1 y2 {qual}, inputs(x1 x2)'
        s+=f'dea {io} dmu(id) engine({engine}) rts({rts}) ort({ori}) stage({stage}) tolerance({tol}) gen(e{stage}) reference replacevars {extra}\n'
        s+=f'matrix L{stage}=r(lambda)\n'
        if not old:
            s+=f'local eng{stage} "`r(engine)\'"\n'
            if stage==1:
                s+='capture matrix list r(slack_in)\nassert _rc==111\ncapture matrix list r(slack_out)\nassert _rc==111\n'
                s+='''if "`eng1'"=="ipm" {
 matrix S1=r(stage_status)
 assert r(n_radial_failed)==0 & r(n_stage2_failed)==0 & r(n_failed)==0
 mata: assert(all(st_matrix("S1")[.,1]:==0) & all(st_matrix("S1")[.,2]:>=.))
}
'''
            else:
                s+='''matrix SI=r(slack_in)
matrix SO=r(slack_out)
assert rowsof(SI)==`n' & colsof(SI)==2
assert rowsof(SO)==`n' & colsof(SO)==2
matrix ST=r(stage_status)
 scalar nf=r(n_failed)
 scalar ns=r(n_stage2_failed)
 assert r(n_radial_failed)==0
 assert "`:colnames ST'"=="radial stage2"
 mata: assert(st_numscalar("nf")==sum(st_matrix("ST")[.,2]:!=0))
 assert nf==ns
'''
    s+='assert sentinel=="preserve me"\nassert missing(e1,e2) if !selected\nkeep if selected\n'
    # Exact double comparison, independent of CSV formatting/round-tripping.
    s+='''mata:
bh=bufio();bufbyteorder(bh,2)
fh=fopen("history_scores.bin","w")
fbufput(bh,fh,"%8z",st_data(.,tokens("e1 e2")));fclose(fh)
fh=fopen("history_lambda1.bin","w")
fbufput(bh,fh,"%8z",st_matrix("L1"));fclose(fh)
end
'''
    s+='svmat double L1, names(l1_)\nsvmat double L2, names(l2_)\n'
    columns='id x1 x2 y1 y2 e1 e2 l1_* l2_*'
    if not old:
        s+='svmat double SI, names(si_)\nsvmat double SO, names(so_)\n'
        s+='if "`eng2\'"=="ipm" svmat double ST, names(status_)\n'
        s+='if "`eng2\'"=="ipm" local statuscols "status_*"\n';columns+=" si_* so_* `statuscols'"
    return s+'recast double x1 x2 y1 y2\n'+f'format {columns} %24.17g\nexport delimited {columns} using "scores.csv", replace datafmt\n'

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stata',required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--baseline-ado',type=Path);p.add_argument('--only');args=p.parse_args()
    out=args.output.resolve();out.mkdir(parents=True);pkg=out/'current_ado';pkg.mkdir();hashes={}
    for f in package_files(REPO):
        a=(REPO/f).read_bytes();(pkg/Path(f).name).write_bytes(a);hashes[f]=hashlib.sha256(a).hexdigest()
    audit=out/'instrumented_ado';instrument(pkg,audit)
    retry=out/'forced_retry_ado';retry_instrument(pkg,retry)
    (out/'manifest.json').write_text(json.dumps(hashes,indent=2)+'\n');results=[]
    def run(name,s,source=pkg,check=None):
        if args.only and not re.search(args.only,name):return None
        d=out/name;d.mkdir();h='version 16\nclear all\nset more off\nset linesize 120\n'
        for f in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/f).mkdir();h+=f'sysdir set {f} "{d/f}"\n'
        h+=f'adopath ++ "{source.resolve()}"\n'+INIT;(d/'test.do').write_text(h+s+'\ndi "TEST_PASS"\n')
        row=dict(test=name,passed=False);start=time.monotonic()
        try:
            subprocess.run([args.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=240)
            log=(d/'test.log').read_text(errors='replace');assert 'TEST_PASS' in log.splitlines(),log[-2400:]
            if check:row.update(check(d))
            row['passed']=True
        except Exception as e:row['error']=repr(e)
        row['seconds']=round(time.monotonic()-start,3);results.append(row)
        (out/'results.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(row),flush=True)
        return d
    for fixture,engine,rts,ori,tol in itertools.product(FIXTURES,['ipm','simplex','auto'],['crs','vrs','drs','irs'],['input','output'],['1e-7','1e-10']):
        name=f'grid_{fixture}_{engine}_{rts}_{ori}_{tol}';q='if id>=2 in 2/8' if fixture=='duplicate' else ''
        d=run(name,script(FIXTURES[fixture],engine,rts,ori,tol,qual=q),check=lambda d,rts=rts,ori=ori,tol=tol:oracle(d,rts,ori,tol))
        if d and args.baseline_ado and engine=='ipm' and tol=='1e-7':
            def history(b,d=d):
                before=read(b/'scores.csv');after=read(d/'scores.csv')
                cols=[k for k in before if k in ['e1','e2'] or k.startswith('l1_')]
                assert all(np.array_equal(before[k],after[k],equal_nan=True) for k in cols)
                for f in ['history_scores.bin','history_lambda1.bin']:
                    assert (b/f).read_bytes()==(d/f).read_bytes(),f
                return dict(unchanged_radial_and_stage1_references=True,binary_double_equality=True)
            run(name+'_before',script(FIXTURES[fixture],engine,rts,ori,tol,old=True,qual=q),args.baseline_ado,history)
    for rts,ori,engine in itertools.product(['ccr','bcc','nirs','ndrs'],['input','output'],['ipm','simplex','auto']):
        run(f'alias_{engine}_{rts}_{ori}',script(DATA,engine,rts,ori,'1e-7'),check=lambda d,rts=rts,ori=ori:oracle(d,rts,ori))
    for rts,ori,stage,mode in itertools.product(['crs','vrs','drs','irs'],['input','output'],[1,2],['none','stage2','primal','both']):
        pf=int(mode in ['primal','both']);df=int(mode=='both');sf=2 if mode=='stage2' else 0
        s=FIXTURES['p03_vrs']+f'scalar fault_primal={pf}\nscalar fault_dual={df}\nscalar fault_stage2={sf}\n'
        s+=f'dea x1 x2=y1 y2, dmu(id) rts({rts}) ort({ori}) stage({stage}) engine(ipm) gen(e) reference saveref saving("saved.dta")\n'
        s+='matrix L=r(lambda)\nmatrix ST=r(stage_status)\nlocal refs "`r(refvars)\'"\n'
        s+='scalar nf=r(n_failed)\nscalar nr=r(n_radial_failed)\nscalar ns=r(n_stage2_failed)\n'
        s+=f'assert audit_s1==3\nassert audit_s2=={3 if stage==2 and not pf else 0}\n'
        s+=f'assert nr=={3 if df else 0}\n'
        s+='matrix SS=r(solution_source)\n'
        s+=f'mata: assert(all(st_matrix("SS")[.,1]:=={0 if df else (2 if pf else 1)}))\n'
        ref_expr='(st_matrix("SS")[.,1]:==1)' if stage==1 else '(2:*(st_matrix("ST")[.,2]:==0))'
        s+=f'mata: assert(all(st_matrix("SS")[.,2]:=={ref_expr}))\n'
        if stage==1:
            s+='assert ns==0 & nf==nr\nmata: assert(all(st_matrix("ST")[.,2]:>=.))\n'
            s+=('mata: assert(all(st_matrix("L"):>=.))\n' if pf else 'mata: assert(!hasmissing(st_matrix("L")))\n')
        else:
            s+='matrix SI=r(slack_in)\nmatrix SO=r(slack_out)\n'
            if pf:s+='assert ns==3 & nf==3\nmata: assert(all(st_matrix("ST")[.,2]:>=.))\n'
            if sf:s+='assert ST[2,2]==3\nassert ns>=1 & nf>=1\n'
            s+='''mata:
st=st_matrix("ST");L=st_matrix("L");SI=st_matrix("SI");SO=st_matrix("SO")
assert(st_numscalar("ns")==sum(st[.,2]:!=0))
assert(st_numscalar("nf")==sum((st[.,1]:!=0):|(st[.,2]:!=0)))
for (k=1;k<=rows(st);k++) {
 if (st[k,2]!=0) assert(all(L[k,.]:>=.) & all(SI[k,.]:>=.) & all(SO[k,.]:>=.))
 else assert(!hasmissing((L[k,.],SI[k,.],SO[k,.])))
}
end
'''
        s+='use "saved.dta",clear\nmata: assert(all(st_data(.,tokens(st_local("refs"))):==st_matrix("L")))\n'
        s+='mata: assert(all(st_data(.,tokens("_dea_radial_status _dea_stage2_status _dea_score_source _dea_ref_stage")):==(st_matrix("ST"),st_matrix("SS"))))\n'
        run(f'fault_{rts}_{ori}_{stage}_{mode}',s,audit)
    for rts,ori in itertools.product(['crs','vrs','drs','irs'],['input','output']):
        run(f'retry_{rts}_{ori}',script(DATA,'auto',rts,ori,'1e-7'),retry,check=lambda d,rts=rts,ori=ori:oracle(d,rts,ori))
    for rts,ori in itertools.product(['crs','vrs'],['input','output']):
        data=DATA+'replace x1=. in 2\nreplace y2=.a in 5\n'
        s=script(data,'ipm',rts,ori,'1e-7').replace('gen byte selected=1','gen byte selected=!missing(x1,x2,y1,y2)')
        run(f'missing_sample_{rts}_{ori}',s,check=lambda d,rts=rts,ori=ori:oracle(d,rts,ori))
    for stage,reference in itertools.product([1,2],[False,True]):
        s=FIXTURES['p03_vrs']+f'dea x1 x2=y1 y2, dmu(id) rts(vrs) stage({stage}) engine(ipm) gen(e) '+('reference' if reference else '')+'\n'
        s+='matrix ST=r(stage_status)\nassert rowsof(ST)==3\n'
        s+='capture matrix list r(lambda)\n'+f'assert _rc=={0 if reference else 111}\n'
        s+='capture matrix list r(slack_in)\n'+f'assert _rc=={0 if stage==2 else 111}\n'
        run(f'optional_reference_{stage}_{reference}',s)
    for rts,stage in itertools.product(['crs','vrs','drs','irs'],[1,2]):
        s=FIXTURES['p03_vrs']+f'dea x1 x2=y1 y2, dmu(id) rts({rts}) stage({stage}) engine(ipm) gen(e0) reference\n'
        s+='matrix ST0=r(stage_status)\nmatrix L0=r(lambda)\nscalar NF0=r(n_failed)\n'
        if stage==2:s+='matrix SI0=r(slack_in)\nmatrix SO0=r(slack_out)\n'
        s+='set seed 7123\n'+f'dea x1 x2=y1 y2, dmu(id) rts({rts}) stage({stage}) engine(ipm) gen(e) reference bootstrap reps1(50) replacevars\n'
        s+='assert e==e0\nassert r(n_failed)==NF0\nassert colsof(r(bootstrap))==5\n'
        s+='mata: assert(all(st_matrix("ST0"):==st_matrix("r(stage_status)")))\n'
        s+='mata: assert(all(st_matrix("L0"):==st_matrix("r(lambda)")))\n'
        if stage==2:
            s+='mata: assert(all(st_matrix("SI0"):==st_matrix("r(slack_in)")) & all(st_matrix("SO0"):==st_matrix("r(slack_out)")))\n'
        run(f'bootstrap_preservation_{rts}_{stage}',s)
    # Preserve the pre-existing direct Mata call with eleven arguments.
    s=FIXTURES['p03_vrs']+'gen byte chosen=1\ngen double e=.\ngen double si1=.\ngen double si2=.\ngen double so1=.\ngen double so2=.\n'
    s+='findfile _dea_ipm.ado\nquietly run "`r(fn)\'"\n'
    s+='mata: _dea_ipm_envelopment("x1 x2","y1 y2","e","","si1 si2","so1 so2","chosen","input","vrs",1,1e-7)\n'
    s+='assert abs(e[3]-1/3)<1e-6\nassert abs(si2[3]-1.5)<1e-6\n'
    run('mata_legacy_arguments',s)
    assert all(r['passed'] for r in results),f'{sum(not r["passed"] for r in results)} failed; see results.json'
    print(f'All {len(results)} command/contract cases passed.')

if __name__=='__main__':main()
