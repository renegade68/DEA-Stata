#!/usr/bin/env python3
"""Bounded exact game outputs, small gains, limits, faults and loader contracts.

Synthetic score-rule tests are explicitly distinct from feasible DEA datasets.
Needs Stata 16+, NumPy/SciPy. No Git, R, network or private research data needed.
"""
from _support import package_files
from pathlib import Path
import argparse, hashlib, json, re, shutil, subprocess, time
import numpy as np
from games_witness import F, check_derived, check_scores, capture_results, input_data
from test_core_retries import WRITER, readbin
from test_cooperative_nash import SUFFIX
from test_output_protection import unchanged
REPO=Path(__file__).resolve().parents[1]

def body_for(a,mode='nash',rts='crs',tol=1e-10,m=1,s=1,detail=False):
    xs=[f'x{i+1}' for i in range(m)];ys=[f'y{i+1}' for i in range(s)]
    b=input_data(a,['id']+xs+ys)+f'dea_{mode} '+ ' '.join(xs)+'='+' '.join(ys)+f',dmu(id) rts({rts}) tol({tol:.17g}) gen(res) matrix'+(' detail' if detail else '')+'\n'
    b+=capture_results(mode)+'mata: savebin("data.bin",st_data(.,tokens("id '+' '.join(xs+ys)+'")))\n'
    return b

def integer(z):
    sign=(z>0)-(z<0);z=abs(z);a=[]
    while z:a.append(z%2**26);z//=2**26
    return '('+','.join(map(str,[sign]+(a or [0])))+')'

def rational(q,name='q'):
    q=F(q);return f'{name}.N={integer(q.numerator)};{name}.D={integer(q.denominator)}\n'

def synthetic(values,tol):
    b='quietly _dea_games_load\nmata:\nvoid rulecheck() {\nstruct _T128_Scores scalar G\nstruct _G129_Derived scalar B\nstruct _C84_Q scalar q\nreal rowvector meta\nreal matrix NR\n'
    n=len(values);b+=f'G.D=J({n},4,0);G.Q=J({n},36,0);G.Q[.,7]=J({n},1,2);G.exact=asarray_create("real",1)\n'
    for i,q in enumerate(values,1):b+=rational(q)+f'asarray(G.exact,{i},q);G.D[{i},2]=_c84_value(q)\n'
    b+=f'NR=_dea_games_allocate(G.D[.,2],{tol:.17g},meta);B=_g129_nash(G,{tol:.17g},NR,meta)\n'
    for file,expr in [('results','B.R'),('quality','B.Q'),('reason','B.why'),('aggregate','B.A'),('meta','(J(1,6,0),B.meta)')]:b+=f'savebin("{file}.bin",{expr})\n'
    return b+'}\nrulecheck()\nend\n'

def truth1(a):
    ratios=[F(float(y))/F(float(x)) for _,x,y in a]
    return [q/max(ratios) for q in ratios]

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--only');a=ap.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir();hashes={}
    for n in package_files(REPO):
        raw=(REPO/n).read_bytes();(pkg/Path(n).name).write_bytes(raw);hashes[n]=hashlib.sha256(raw).hexdigest()
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,stata=a.stata),indent=2)+'\n');rows=[]
    def selected(name):return not a.only or re.search(a.only,name)
    def run(name,body,source=pkg,check=None):
        if not selected(name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset processors 1\nset type double\nset linesize 200\n'
        for k in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/k).mkdir();head+=f'sysdir set {k} "{d/k}"\n'
        head+=f'adopath ++ "{source}"\n'+WRITER;(d/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n');t=time.monotonic();record=dict(test=name,passed=False)
        try:
            p=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=900)
            log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else ''
            if p.returncode or not log:(d/'process.json').write_text(json.dumps(dict(returncode=p.returncode,stdout=p.stdout.decode(errors='replace'),stderr=p.stderr.decode(errors='replace'))))
            assert p.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2500:]
            if check:record.update(check(d) or {})
            record['passed']=True
        except Exception as error:record['error']=repr(error)
        record['seconds']=round(time.monotonic()-t,3);rows.append(record);(out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(record),flush=True)
    def exact(a,mode,tol=1e-10,detail=False):
        truth=truth1(a)
        return lambda d:dict(**check_scores(d,truth,1,1,'crs',tol),**check_derived(d,truth,mode,tol,detail))
    # Analytic 1x1 CRS oracle: max(y/x) defines the frontier exactly.
    for exponent in [20,30,40,52]:
        u=2.**-exponent;original=np.array([[1,1+2*u,1+u],[2,1+u,1]],float)
        for variant in ['base','reverse','dyadic']:
            data=original[::-1].copy() if variant=='reverse' else original.copy()
            if variant=='dyadic':data[:,1]*=2.**20;data[:,2]*=2.**-20
            for mode in SUFFIX:run(f'small_e{exponent}_{variant}_{mode}',body_for(data,mode,detail=True),check=exact(data,mode,detail=True))
    for n in [8,12,16,24,32]:
        rng=np.random.default_rng(129+n);original=np.c_[np.arange(1,n+1),rng.uniform(.2,7,(n,2))]
        for variant in ['base','reverse','dyadic']:
            data=original[::-1].copy() if variant=='reverse' else original.copy()
            if variant=='dyadic':data[:,1]*=2.**20;data[:,2]*=2.**-20
            run(f'growth_n{n}_{variant}',body_for(data),check=exact(data,'nash'))
    u=2.**-52;data=np.array([[1,1,2],[2,1+u,1],[3,1+2*u,1+u]])
    for mode in SUFFIX:run('rounded_collision_'+mode,body_for(data,mode,detail=True),check=exact(data,mode,detail=True))
    # Synthetic fractions exercise formula boundaries, not an all-inefficient
    # observed DEA frontier. Exact 1/100 and 1/5 constants belong to the rule.
    cases=[('zero',[F(0)]*3,1e-10),('one',[F(1)]*3,1e-10),('capped',[F(1),F(1,2),F(1,4)],1e-10),('positive',[F(1,2)]*3,1e-10),('floor',[F(99,100)-F(1,10**18),F(99,100),F(99,100)+F(1,10**18)],1e-10),('threshold',[1-F(1e-7)-F(1,10**22),1-F(1e-7),1-F(1e-7)+F(1,10**22)],1e-7),('strict',[F(1,2)]*2,1e-18),('subnormal_exact',[1-F(1,2**40)]*26,1e-10),('subnormal_inexact',[1-F(9,10*2**40)]*26,1e-10),('product_underflow',[1-F(1,2**40)]*27,1e-10),('gain_underflow',[1-F(1,2**1100)]*2,1e-10)]
    for tag,values,tol in cases:run('rule_'+tag,synthetic(values,tol),check=lambda d,v=values,t=tol:check_derived(d,v,'nash',t))
    # n and multiplier-column scope boundaries. Outside scope the native
    # estimate is retained, explicitly labelled unverified, never certified.
    for mode in SUFFIX:
        for rts in ['crs','vrs']:
            for n,p in [(32,4),(33,4),(2,16),(2,17)]:
                m=(p-2*(rts=='vrs'))//2;s=p-2*(rts=='vrs')-m
                data=[[i+1]+[1.+i%2]*m+[1.]*s for i in range(n)]
                def scope(d,n=n,p=p):
                    Q=readbin(d/'score_quality.bin');D=readbin(d/'lp_diagnostics.bin');N=readbin(d/'native_lp.bin')
                    if n<=32 and p<=16:assert np.all(Q[:,6]==2)
                    else:
                        assert np.all(Q[:,7]==(201 if n>32 else 202)) and np.array_equal(D,N,equal_nan=True) and np.all(Q[:,6]==N[:,2])
                        assert 'Unverified numerical DEA scores' in (d/'test.log').read_text()
                    return dict(rows=n,certified=int(sum(Q[:,6]==2)))
                run(f'scope_{mode}_{rts}_n{n}_p{p}',body_for(data,mode,rts,m=m,s=s),check=scope)
    def changed(tag,file,old,new):
        target=out/('fault_ado_'+tag);shutil.copytree(pkg,target);p=target/file;text=p.read_text();assert text.count(old)==1,(tag,text.count(old));p.write_text(text.replace(old,new));return target
    data=[[1,3,1],[2,2,1],[3,1,1]]
    faults=[('total_budget','W=_c84_work();total=','W=_c84_work();W.limit=1;total='),('weights_budget','W=_c84_work();wsum=','W=_c84_work();W.limit=1;wsum='),('ratio_budget','W=_c84_work();unit=','W=_c84_work();W.limit=1;unit='),('row_budget','share=_g129_mul(unit,w,W);','if(i==1) W.limit=1\n  share=_g129_mul(unit,w,W);'),('summary_budget','W=_c84_work();allocated=','W=_c84_work();W.limit=1;allocated='),('mean_budget','meanq=_g129_div(allocated,','W.limit=1;meanq=_g129_div(allocated,'),('scenario_width','W=_c84_work();one=_c84_number(1,W);coef=','W=_c84_work();W.limbs=1;one=_c84_number(1,W);coef='),('bad_export','value=_c84_value(q)','value=_c84_value(q)*.5'),('wrong_division','q=_c58b_div(a,b,W)','q=_c58b_div(a,b,W);q[2]=q[2]+1')]
    for tag,old,new in faults:
        if not selected('fault_'+tag):continue
        target=changed(tag,'_dea_games_quality.ado',old,new)
        def faultcheck(d,tag=tag):
            Q=readbin(d/'quality.bin');A=readbin(d/'aggregate.bin');R=readbin(d/'results.bin');M=readbin(d/'meta.bin').ravel();assert np.all(Q[:,0]==2)
            if tag=='total_budget':assert np.all(Q[:,:2]==2) and np.all(Q[:,2:]==3) and np.all(A[:,0]==3)
            elif tag in ['weights_budget','ratio_budget']:assert np.all(Q[:,:2]==2) and np.all(Q[:,2:]==3) and A[0,0]==2 and np.all(A[1:,0]==3)
            elif tag=='row_budget':assert np.all(Q[:,:2]==2) and np.all(Q[0,2:]==3) and np.all(Q[1:,2:]==2) and A[0,0]==2 and np.all(A[1:,0]==3) and np.isnan(M[12:14]).all() and M[15]==0
            elif tag=='summary_budget':assert np.all(Q==2) and list(A[:5,0])==[2,3,2,3,2]
            elif tag=='mean_budget':assert np.all(Q==2) and list(A[:5,0])==[2,2,2,3,2]
            elif tag=='scenario_width':assert np.all(Q[:,1:]==3) and np.all(A[:,0]==3)
            elif tag=='bad_export':assert np.all(Q[:,1:3]==3) and Q[2,3]==2 and R[2,3]==0 and np.all(A[:4,0]==3) and A[4,0]==2
            elif tag=='wrong_division':assert np.any(Q[:,1:]==3) and np.all(A[:,0]==3)
            assert 'Warning:' in (d/'test.log').read_text()
        run('fault_'+tag,body_for(data,rts='vrs'),target,faultcheck)
    for mode in SUFFIX:
        for helper in ['_dea_games_load.ado','_dea_games_exact.ado','_dea_games_quality.ado']:
            name='missing_'+helper+'_'+mode
            if not selected(name):continue
            target=out/('fault_ado_'+name);shutil.copytree(pkg,target);(target/helper).unlink()
            run(name,input_data(data,['id','x','y'])+f'gen str8 res_{SUFFIX[mode].split()[0]}="keep"\n'+unchanged(f'dea_{mode} x=y,gen(res) replacevars',601),target)
        body=input_data(data,['id','x','y'])+'scalar user_scalar=123\nmatrix user_matrix=(9,8)\n'
        body+=f'dea_{mode} x=y,gen(first) matrix\n'
        for i,fn in enumerate(['_g130_sum','_g129_nash','_t128_certify','_dea_games_compute'],1):
            body+=f'mata: mata drop {fn}()\ndea_{mode} x=y,gen(next{i}) matrix\nassert r(n_certified)==3\nassert scalar(user_scalar)==123\nmata: asserteq(st_matrix("user_matrix"),(9,8))\n'
            for suffix in SUFFIX[mode].split():body+=f'assert next{i}_{suffix}==first_{suffix}\n'
        run('reload_'+mode,body)
    # Both proof and native pathways unavailable: keep explicit missing values.
    if any(selected('both_unavailable_'+mode) for mode in SUFFIX):
        target=changed('both','_dea_independent_exact.ado','if(n>32)','if(n>0)');p=target/'_dea_games.ado';text=p.read_text();old='nrc=_dea_nonte_optimize(S,eps)';assert text.count(old)==1;p.write_text(text.replace(old,old+';nrc=1'))
        for mode in SUFFIX:
            def absent(d):
                assert np.all(readbin(d/'score_quality.bin')[:,6]==0) and np.isnan(readbin(d/'results.bin')).all()
                assert readbin(d/'counts.bin').ravel()[7]==0 and 'Unavailable DEA scores: 3.' in (d/'test.log').read_text()
            run('both_unavailable_'+mode,body_for(data,mode),target,absent)
    # Regression inputs that previously exhausted one monolithic work budget.
    # Keep RNG state across sizes to reproduce that specific positive dataset.
    rng=np.random.default_rng(129)
    for n in [8,12,16,24,32]:
        original=np.c_[np.arange(1,n+1),rng.uniform(.2,4,n),rng.uniform(.2,4,n)]
        for variant in ['base','reverse','dyadic']:
            data=original[::-1].copy() if variant=='reverse' else original.copy()
            if variant=='dyadic':data[:,1]*=2.**20;data[:,2]*=2.**-20
            run(f'prior_growth_n{n}_{variant}',body_for(data),check=exact(data,'nash'))
    # Aggregate-only failures must not erase certified rows or unrelated totals.
    for tag,old,new in [('loss_budget','W=_c84_work();unallocated=','W=_c84_work();W.limit=1;unallocated='),('product_budget','WP=_c84_work();product=','WP=_c84_work();WP.limit=1;product='),('log_budget','W=_c84_work();q=source','W=_c84_work();W.limit=1;q=source')]:
        if not selected('fault_'+tag):continue
        target=changed(tag,'_dea_games_quality.ado',old,new);values=[F(1),F(1),F(1,2)] if tag=='loss_budget' else [F(1,2)]*3
        def aggregate_fault(d,tag=tag):
            Q=readbin(d/'quality.bin');A=readbin(d/'aggregate.bin');M=readbin(d/'meta.bin').ravel();assert np.all(Q==2)
            if tag=='loss_budget':assert A[0,0]==2 and np.all(A[1:4,0]==3) and A[4,0]==2 and A[4,2]==0 and M[12]==M[13]==2 and M[15]==0
            elif tag=='product_budget':assert np.all(A[:4,0]==2) and A[4,0]==3 and A[5,0]==1 and M[15]==1
            else:assert np.all(A[:5,0]==2) and A[5,0]==3 and M[15]==1
        run('fault_'+tag,synthetic(values,1e-10),target,aggregate_fault)
    if selected('fault_coalition_budget'):
        target=changed('coalition_budget','_dea_games.ado','shown=1;WC=_c84_work()','shown=1;WC=_c84_work();WC.limit=1')
        def coalition_fault(d):
            assert np.all(readbin(d/'quality.bin')==2) and readbin(d/'aggregate.bin')[0,0]==2
            assert '7 output(s) withheld' in (d/'test.log').read_text()
        run('fault_coalition_budget',body_for([[1,3,1],[2,2,1],[3,1,1]],'cooperative',detail=True),target,coalition_fault)
    for mode in SUFFIX:
        for tag,data,tol in [('strict_fraction',[[1,3,1],[2,1,1]],1e-18),('strict_dyadic',[[1,2,1],[2,1,1]],1e-18),('score_underflow',[[1,2.**600,1],[2,2.**-500,1]],1e-10)]:
            def score_boundary(d,tag=tag):
                Q=readbin(d/'score_quality.bin')
                if tag=='strict_dyadic':assert np.all(Q[:,6]==2)
                else:assert Q[0,6]!=2 and Q[0,7]!=0 and any(message in (d/'test.log').read_text() for message in ['Unavailable DEA scores:', 'Unverified numerical DEA scores:'])
            run(tag+'_'+mode,body_for(data,mode,tol=tol),check=score_boundary)
    summary=dict(total=len(rows),passed=sum(r['passed'] for r in rows),failed=sum(not r['passed'] for r in rows));(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary);return bool(summary['failed'] or not rows)
if __name__=='__main__':raise SystemExit(main())
