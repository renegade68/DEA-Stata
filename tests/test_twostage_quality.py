#!/usr/bin/env python3
"""Two-stage exact aggregate export, provenance, bounds and reload contracts.

Self-contained synthetic data; no Git, private folders, R or downloaded oracles.
Requires licensed Stata16+ and NumPy. Raw outputs are local diagnostic evidence.
"""
from _support import package_files,executable
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,itertools,json,re,shutil,subprocess,time
import numpy as np
from test_core_retries import WRITER,readbin

REPO=Path(__file__).resolve().parents[1]
MATS=['efficiencies','status','score_quality','score_reason','native_efficiencies',
      'independent_efficiencies','independent_quality','independent_reason',
      'native_independent','lambda_stage1','lambda_stage2',
      'native_lambda_stage1','native_lambda_stage2']

def body_for(a,shape,rts='crs',ort='input',tol=1e-10,method='independent',stack=False,leader='stage1'):
    m,d,s=shape;xs=[f'x{i+1}' for i in range(m)];zs=[f'z{i+1}' for i in range(d)];ys=[f'y{i+1}' for i in range(s)];names=['id']+xs+zs+ys
    body='input double('+' '.join(names)+')\n'+'\n'.join(' '.join(format(v,'.17g') for v in row) for row in a)+'\nend\n'
    body+=f'{"dea_stackelberg" if stack else "dea_network"}, inputs({" ".join(xs)}) intermediates({" ".join(zs)}) outputs({" ".join(ys)}) method({method}) rts({rts}) '+(f'leader({leader})' if stack else f'ort({ort})')+f' gen(res) tol({tol:.17g}) matrix detail\n'
    for name in MATS+([] if method=='noncooperative' else ['aggregate_diagnostics']):
        body+=f'mata: savebin("{name}.bin",st_matrix("r({name})"))\n'
    body+='mata: savebin("counts.bin",(st_numscalar("r(n_certified_overall)"),st_numscalar("r(n_withheld_overall)"),st_numscalar("r(n_unverified_overall)")))\n'
    return body

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--only');a=ap.parse_args()
    stata=executable('stata',a.stata);out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir();manifest={};rows=[]
    for name in package_files(REPO):
        raw=(REPO/name).read_bytes();(pkg/Path(name).name).write_bytes(raw);manifest[name]=hashlib.sha256(raw).hexdigest()
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=manifest),indent=2)+'\n')
    def run(name,body,check=None,source=pkg):
        if a.only and not re.search(a.only,name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset processors 1\nset type double\nset linesize 180\n'
        for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/key).mkdir();head+=f'sysdir set {key} "{d/key}"\n'
        (d/'test.do').write_text(head+f'adopath ++ "{source}"\n'+WRITER+body+'\ndi "TEST_PASS"\n');row=dict(test=name,passed=False);t=time.monotonic()
        try:
            p=subprocess.run([stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=900)
            (d/'stdout.txt').write_bytes(p.stdout);(d/'stderr.txt').write_bytes(p.stderr)
            log=(d/'test.log').read_text(errors='replace');assert p.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2000:]
            if check:row.update(check(d) or {})
            row['passed']=True
        except Exception as e:row['error']=repr(e)
        row['seconds']=round(time.monotonic()-t,3);rows.append(row);(out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(row),flush=True)
    cases=[]
    for e,rts,ort in [(e,'crs',o) for e in [500,537,538,600] for o in ['input','output']]+[(600,'vrs',o) for o in ['input','output']]:
        K=2.**e;cases.append((f'product_e{e}_{rts}_{ort}',[[1,K,1,1],[2,1,1,K]],rts,ort,1e-10,'independent',F(2)**(-2*e if rts=='crs' else -e)))
    for ort in ['input','output']:
        cases.append((f'nonzero_subnormal_{ort}',[[1,2.**538,1,1],[2,3,1,2.**538]],'crs',ort,1e-10,'independent',3*F(2)**-1076))
        cases.append((f'half_subnormal_{ort}',[[1,2.**537,1,1],[2,1,1,2.**538]],'crs',ort,1e-10,'independent',F(2)**-1075))
        for tol in [1e-10,1e-18]:cases.append((f'mean_precision_{ort}_{tol}',[[1,1,1,1],[2,1,1,2.**53]],'crs',ort,tol,'additive',(1+F(2)**-53)/2))
        for method in ['independent','additive']:cases.append((f'fraction_{method}_{ort}',[[1,7,1,1],[2,3,1,11]],'crs',ort,1e-10,method,F(3,77) if method=='independent' else (F(3,7)+F(1,11))/2))
    def boundary(d,truth,tol):
        R=readbin(d/'efficiencies.bin');Q=readbin(d/'score_quality.bin');S=readbin(d/'score_reason.bin');A=readbin(d/'aggregate_diagnostics.bin');status=readbin(d/'status.bin')
        value=float(truth);error=abs(F(value)-truth)/truth;reason=1101 if value==0 else 1103 if error>F(tol) else 0
        assert np.all(Q[:,1:]==2) and np.all(readbin(d/'independent_quality.bin')==2)
        assert Q[0,0]==(3 if reason else 2) and S[0,0]==reason and status[0,0]==(not reason)
        assert np.isnan(R[0,0]) if reason else R[0,0]==value
        if value>0:assert A[0,3]==float(error)
        if reason:assert 'withheld after exact aggregation/export checks' in (d/'test.log').read_text()
        assert np.all(R[:,1:]>0)
        return dict(reason=reason,relative_export_error=float(error))
    for name,data,rts,ort,tol,method,truth in cases:run(name,body_for(data,[1,1,1],rts,ort,tol,method),lambda d,t=truth,v=tol:boundary(d,t,v))
    for rts,ort in itertools.product(['crs','vrs'],['input','output']):
        for n,p in [(32,4),(33,4),(2,16),(2,17)]:
            m=2 if p==4 else (p-2*(rts=='vrs'))//2;z=2 if p==4 else p-2*(rts=='vrs')-m
            data=[[i+1]+[1.+i%2]*m+[1.]*z+[2.]*m for i in range(n)]
            def bounds(d,n=n,p=p):
                Q=readbin(d/'score_quality.bin');S=readbin(d/'score_reason.bin');R=readbin(d/'efficiencies.bin')
                if n<=32 and p<=16:assert np.all(Q==2)
                else:
                    assert np.all(S[:,0]==1001) and np.all(S[:,1:]==(201 if n>32 else 202))
                    assert np.array_equal(R,readbin(d/'native_efficiencies.bin'),equal_nan=True) and np.array_equal(Q,np.isfinite(R))
                    for side in [1,2]:assert np.array_equal(readbin(d/f'lambda_stage{side}.bin'),readbin(d/f'native_lambda_stage{side}.bin'),equal_nan=True)
            run(f'bound_{rts}_{ort}_n{n}_p{p}',body_for(data,[m,z,m],rts,ort),bounds)
    data='input double(id x z y)\n1 3 1 1\n2 1 1 1\nend\n';command='dea_network x=y, intermediates(z) method(independent) rts(vrs)'
    context='''mata:
_rsm_success_policy=2;_rsm_success_recording=1;_dea_link_record=1
saved=_c57_context_save()
end
'''
    check_context='''mata:
after=_c57_context_save();names=_c57_context_names()
for(j=1;j<=cols(names);j++) {
 assert(asarray(saved,names[j]+":existed")==asarray(after,names[j]+":existed"))
 if(asarray_contains(saved,names[j])) assert(asarray(saved,names[j])==asarray(after,names[j]))
}
end
'''
    body=data+'quietly _dea_crossgame_load\n'+context+command+' gen(first) matrix\n'+check_context
    for i,fn in enumerate(['_t123_aggregate','_t123_stageq','_t122_stage','_t120_eval','_c84_mul','_dea_ts_compute','_t121_nc'],1):
        body+=f'mata: mata drop {fn}()\n'+command+f' gen(next{i}) matrix\nassert next{i}_overall==first_overall & next{i}_stage1==first_stage1 & next{i}_stage2==first_stage2\nassert r(n_certified_overall)==2\n'+check_context
    run('context_partial_reload',body)
    for stack in [False,True]:
        cmd=('dea_stackelberg if id<=2 in 1/3, inputs(x) intermediates(z) outputs(y)' if stack else 'dea_network x=y if id<=2 in 1/3, intermediates(z) method(independent)')
        body=data+'set obs 3\nreplace id=3 in 3\n'+cmd+' gen(nomat)\nassert r(n)==2 & r(n_certified_overall)==2 & r(n_withheld_overall)==0\n'
        body+='mata: assert(rows(st_matrix("r(score_quality)"))==0)\nassert missing(nomat_overall) in 3\n'
        body+=cmd+' gen(nomat) replacevars matrix\nassert r(n_certified_stage1)==2 & r(n_certified_stage2)==2\nmata: assert(all(st_matrix("r(score_quality)"):==2))\n'
        run('no_matrix_filter_replace_'+str(stack),body)
    for lead in [1,2]:
        a0=[[1,1,1,0,1],[2,1,2,1,1],[3,1,0,1,1]]
        def zero(d,lead=lead):
            R=readbin(d/'efficiencies.bin');Q=readbin(d/'score_quality.bin');assert np.all(Q==2)
            wanted=[[0,.5,0],[1/3,1,1/3],[0,1,0]] if lead==1 else [[1/3,1/3,1],[1/3,1,1/3],[1/3,1/3,1]]
            assert np.array_equal(R,wanted)
            assert np.all(readbin(d/'independent_quality.bin')==2)
        run(f'exact_zero_follower_leader{lead}',body_for(a0,[1,2,1],method='noncooperative',stack=True,leader=f'stage{lead}'),zero)
    for name,needle,repl,reason in [('wrong_value','value=_c84_value(q)','value=_c84_value(q)*.5',1103),('zero_value','value=_c84_value(q)','value=0',1101),('work_limit','W=_c84_work();zero=','W=_c84_work();W.limit=1;zero=',None),('runtime_error','W=_c84_work();zero=','_error(3498);W=_c84_work();zero=',3498)]:
        target=out/('fault_'+name);shutil.copytree(pkg,target);p=target/'_dea_twostage_quality.ado';s=p.read_text();assert needle in s;p.write_text(s.replace(needle,repl))
        if name=='runtime_error':run(name,data+'quietly _dea_crossgame_load\n'+context+'capture noisily '+command+' gen(res) matrix\nscalar caught=_rc\nassert scalar(caught)==3498\n'+check_context,source=target)
        else:
            def fault(d,reason=reason):
                Q=readbin(d/'score_quality.bin');S=readbin(d/'score_reason.bin');assert np.all(Q[:,0]==3) and np.all(Q[:,1:]==2)
                assert np.all(S[:,0]==reason) if reason else np.all(S[:,0]>1200)
            run(name,body_for([[1,3,1,1],[2,1,1,1]],[1,1,1]),fault,source=target)
    for name,cmd in [('relational','dea_network, inputs(x) intermediates(z) outputs(y)'),('noncooperative','dea_stackelberg, inputs(x) intermediates(z) outputs(y) method(noncooperative)')]:
        for n in [2,33]:
            body=f'set obs {n}\ngen double x=1+mod(_n,2)\ngen double z=1\ngen double y=1\n'+cmd+' gen(res) matrix\n'
            body+=f'assert r(n_certified_overall)=={n if n==2 else 0}\nassert r(n_withheld_overall)==0\n'
            if n==33:body+='mata: assert(all(st_matrix("r(score_reason)"):==201))\nmata: assert(all(st_matrix("r(efficiencies)"):==st_matrix("r(native_efficiencies)")))\n'
            body+=cmd+f' gen(second)\nassert r(n_certified_overall)=={n if n==2 else 0}\nmata: assert(rows(st_matrix("r(score_quality)"))==0)\n'
            run(f'linked_retention_{name}_n{n}',body)
    def partial(d):
        Q=readbin(d/'score_quality.bin');S=readbin(d/'score_reason.bin');R=readbin(d/'efficiencies.bin')
        assert np.all(Q[:,1]==2) and np.all(Q[:,2]!=2) and np.all(S[:,0]==1001) and np.all(S[:,2]==202)
        assert np.array_equal(R[:,0],R[:,1]*R[:,2],equal_nan=True)
    run('partial_stage_certification',body_for([[1,3,1]+[1]*16,[2,1,1]+[2]*16],[1,1,16]),partial)
    summary=dict(total=len(rows),passed=sum(q['passed'] for q in rows),failed=sum(not q['passed'] for q in rows))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary);return bool(summary['failed'])

if __name__=='__main__':raise SystemExit(main())
