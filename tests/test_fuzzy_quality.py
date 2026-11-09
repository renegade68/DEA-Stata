#!/usr/bin/env python3
"""Fuzzy alpha=1 exact certificates, recovery, scope and output contracts.

Requires Stata16+, NumPy/SciPy. No Git, R, network or private research inputs.
"""
from _support import package_files
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,itertools,json,re,shutil,subprocess,time
import numpy as np
from test_core_retries import WRITER,readbin
from test_output_protection import unchanged
from fuzzy_witness import grid_body,check_grid
REPO=Path(__file__).resolve().parents[1]

def faults(run,clean,out):
    DATA='input double id x y\n1 1 1\n2 2 1\n3 4 1\nend\n'
    CMD='dea_fuzzy, inputs(x) outputs(y) dmu(id) alpha(1) gen(res) matrix'
    OLD='gen str12 res_pess_a100="old pess"\ngen str12 res_opt_a100="old opt"\n'
    def case(name,body,source=clean,verify=None):
     def check(d):
      text=(d/'test.log').read_text()
      if name.startswith('retain_') or name in ['budget_retain','dimension_guard','sample_guard']:
       assert 'retain numerical scores without exact certification' in text
      if name.startswith('missing_') and name!='missing_helper_protection' or name=='budget_missing':
       assert 'remain unavailable; certification limits do not prove infeasibility' in text
      if name.startswith('recover_'):
       assert 'scores recovered for' in text and 'native failure reports remain in center_diagnostics' in text
      return verify(d) if verify else {}
     return run(name,body,source=source,check=check)
    def package(name,native=None,reason=None,point=None,budget=False):
     d=out/'fault_ado'/name;shutil.copytree(clean,d)
     if native:
      p=d/'dea_fuzzy.ado';s=p.read_text();needle='status=S.errorcode()';assert s.count(needle)==1;s=s.replace(needle,'status='+native);p.write_text(s)
     if reason:
      p=d/'_dea_fuzzy_exact.ado';s=p.read_text();needle='R.diag[2]=2;R.T=_c58_map';assert needle in s;s=s.replace(needle,f'return(_f135_stop(R,{reason}));\n '+needle);p.write_text(s)
     if point:
      p=d/'_dea_fuzzy_exact.ado';s=p.read_text();needle='R.weights=_c58_point(R.T,R.primary,1)';assert needle in s;s=s.replace(needle,needle+';'+point);p.write_text(s)
     if budget:
      p=d/'_dea_crosseff_cert.ado';s=p.read_text().replace('W.limit=budget','W.limit=0');p.write_text(s)
     return d
    for status in ['1','2','3','attempt']:
     pkg=package('native_'+status,native=status)
     for ori in ['input','output']:
      check=DATA+CMD+f' ort({ori})\nassert r(n_center_certified)==3 & r(n_center_recovered)==3 & r(lp_missing)==0\nassert r(lp_count)==6 & r(lp_retries)==6\nassert res_pess_a100==res_opt_a100\nassert abs(res_pess_a100-1/x)<1e-15\n'
      check+=f'assert r(lp_infeasible)=={6 if status=="1" else 0} & r(lp_unbounded)=={6 if status=="2" else 0}\nmatrix D=r(center_diagnostics)\n'
      check+='mata: assert(all(st_matrix("D")[,(5,11)]:>=.) & all(st_matrix("D")[,(6,12)]:==1))\n'
      flag=7 if status=='1' else 8 if status=='2' else 9
      check+=f'mata: assert(all(st_matrix("D")[,({flag},{flag+6})]:==1))\n'
      case('recover_'+status+'_'+ori,check,pkg)
    for name,opts,reason in [('scope',dict(reason=201),201),('export_score',dict(point='R.score=0'),701),('export_weight',dict(point='R.weights[1]=.'),702)]:
     pkg=package(name,**opts)
     case('retain_'+name,DATA+CMD+f'\nassert r(n_center_certified)==0 & r(n_center_unverified)==3\nmatrix D=r(center_diagnostics)\nassert D[1,3]=={reason}\nassert res_pess_a100>0 & res_opt_a100>0\n',pkg)
     pkg=package(name+'_missing',native='1',**opts)
     case('missing_'+name,DATA+CMD+f'\nassert r(n_center_certified)==0 & r(n_center_unavailable)==3 & r(lp_missing)==6\nmatrix D=r(center_diagnostics)\nassert D[1,3]=={reason}\nassert missing(res_pess_a100,res_opt_a100)\n',pkg)
     pkg=package(name+'_abort',native='3',**opts)
     case('abort_'+name,DATA+OLD+unchanged(CMD+' replacevars',498),pkg)
    pkg=package('budget',budget=True)
    case('budget_retain',DATA+CMD+'\nassert r(n_center_certified)==0 & r(n_center_unverified)==3\nassert res_pess_a100>0\n',pkg)
    pkg=package('budget_missing',native='1',budget=True)
    case('budget_missing',DATA+CMD+'\nassert r(n_center_unavailable)==3 & r(lp_missing)==6\nassert missing(res_pess_a100,res_opt_a100)\n',pkg)
    # Original non-center numerical abort semantics still protect permanent results.
    for status in ['1','2','3','attempt']:
     pkg=package('other_alpha_'+status,native=status)
     cmd=CMD.replace('alpha(1)','alpha(.5)')
     if status in ['1','2']:body=DATA+cmd+'\nassert r(lp_missing)==6 & r(center_attempts)==0\nassert r(n_center_certified)==0\n'
     else:body=DATA+'gen str12 res_pess_a50="old"\ngen str12 res_opt_a50="old"\n'+unchanged(cmd+' replacevars',498)
     case('other_alpha_'+status,body,pkg)
    # No center cut means no new-helper dependency; a missing helper at alpha1 must
    # fail before replacement. Existing installation is never changed.
    pkg=out/'fault_ado/missing_helper';shutil.copytree(clean,pkg);(pkg/'_dea_fuzzy_exact.ado').unlink()
    case('no_center_without_helper',DATA+CMD.replace('alpha(1)','alpha(.5)')+'\nassert r(center_attempts)==0\n',pkg)
    case('missing_helper_protection',DATA+OLD+unchanged(CMD+' replacevars',601),pkg)
    # Scope guards: explicit diagnostic reasons and preserved numerical outcomes.
    case('dimension_guard', 'set obs 1\ngen id=1\n'+''.join(f'gen double x{i}=1\n' for i in range(1,17))+'gen double y=1\n'+CMD.replace('inputs(x)','inputs('+ ' '.join(f'x{i}' for i in range(1,17))+')')+'\nassert r(n_center_unverified)==1\nmatrix D=r(center_diagnostics)\nassert D[1,3]==202\n')
    case('sample_guard','set obs 33\ngen id=_n\ngen double x=_n\ngen double y=1\n'+CMD+'\nassert r(n_center_certified)==0\nassert r(n_center_unverified)+r(n_center_unavailable)==33\nmatrix D=r(center_diagnostics)\nassert D[1,3]==201\n')
    # Save/restore all shared solver state for both successful and exceptional exits.
    pre='quietly _dea_fuzzy_load\nmata: before135=_c57_context_save()\n'
    post='''mata:
    after135=_c57_context_save();nn135=_c57_context_names()
    for(i135=1;i135<=cols(nn135);i135++) {
     assert(asarray(before135,nn135[i135]+":existed")==asarray(after135,nn135[i135]+":existed"))
     if(asarray(before135,nn135[i135]+":existed")) asserteq(asarray(before135,nn135[i135]),asarray(after135,nn135[i135]))
    }
    end
    '''
    case('context_success',DATA+pre+CMD+'\n'+post)
    pkg=out/'fault_ado/throw';shutil.copytree(clean,pkg);p=pkg/'_dea_fuzzy_exact.ado';p.write_text(p.read_text().replace('R=_f135_stage(X,Y,k,rts,ort,tol)','_error(3498);R=_f135_stage(X,Y,k,rts,ort,tol)'))
    case('context_exception',DATA+OLD+pre+unchanged(CMD+' replacevars',3498)+post,pkg)
    case('reload_missing_function',DATA+CMD+'\nmata: mata drop _f135_stage()\n'+CMD+' replacevars\nassert r(n_center_certified)==3\ndiscard\n'+CMD+' replacevars\nassert r(n_center_certified)==3\n')


def contracts(run):
    case=run
    for rts,ori in itertools.product(['crs','vrs','drs','irs'],['input','output']):
     data='''input double x y
    -1 -1
    1 1
    2 1
    4 1
    -1 -1
    end
    gen str8 id="unit"+string(_n)
    gen double xl=x/2
    gen double xu=x*2
    gen double yl=y/2
    gen double yu=y*2
    '''
     cmd=f'dea_fuzzy if _n!=5 in 2/5, inputs(x) outputs(y) dmu(id) rts({rts}) ort({ori})'
     body=data+cmd+' alpha(1) gen(a)\n'+'''assert r(N)==3 & r(n_center_certified)==3 & r(center_attempts)==3
    assert r(n_center_unverified)==0 & r(n_center_unavailable)==0
    assert "`r(center_scope)'"=="alpha1_original_data"
    matrix D=r(center_diagnostics)
    assert rowsof(D)==3 & colsof(D)==43
    assert D[1,1]==2 & D[2,1]==3 & D[3,1]==4
    assert D[1,2]==2 & D[2,2]==2 & D[3,2]==2
    assert missing(a_pess_a100,a_opt_a100) if _n==1 | _n==5
    '''
     body+=cmd+' alpha(1 .5 0) xlo(xl) xhi(xu) ylo(yl) yhi(yu) gen(b) matrix\n'
     body+='''assert r(n_center_certified)==3 & r(center_attempts)==3
    assert a_pess_a100==b_pess_a100 & a_opt_a100==b_opt_a100
    assert b_pess_a100==b_opt_a100
    replace xl=x
    replace xu=x
    replace yl=y
    replace yu=y
    '''
     body+=cmd+' alpha(1) xlo(xl) xhi(xu) ylo(yl) yhi(yu) gen(b) replacevars\n'
     body+='assert a_pess_a100==b_pess_a100 & a_opt_a100==b_opt_a100\n'
     case('sample_bounds_'+rts+'_'+ori,body)
     # Maximum supported dimensions, including redundant zero columns.
     m=15-(2 if rts=='vrs' else 0 if rts=='crs' else 1)
     body='set obs 32\ngen id=_n\ngen double x1=1+mod(_n-1,4)\ngen double y=1\n'
     body+=''.join(f'gen double x{j}=0\n' for j in range(2,m+1))
     body+=f'dea_fuzzy, inputs({" ".join("x"+str(j) for j in range(1,m+1))}) outputs(y) dmu(id) rts({rts}) ort({ori}) alpha(1) gen(a)\n'
     body+='assert r(N)==32 & r(n_center_certified)==32 & r(center_attempts)==32\nassert r(lp_missing)==0\nassert a_pess_a100==a_opt_a100\n'
     expected='1' if ori=='output' and rts in ['vrs','drs'] else '1/x1'
     body+=f'assert abs(a_pess_a100-({expected}))<=1e-15\n'
     case('scope_boundary_'+rts+'_'+ori,body)


def interleave(run):
    groups={
     'cooperative':('dea_cooperative, inputs(x) outputs(y) dmu(id)', ['results','lp_diagnostics','score_diagnostics','native_lp_diagnostics','result_quality','result_reason','aggregate_diagnostics','numerical_results','numerical_aggregates','lambda','weights']),
     'crossgame':('dea_crossgame, inputs(x) outputs(y) dmu(id) strategy(all)', ['self_eff','weights_input','weights_output','ce_initial','ce_agg','ce_ben','ce_nash','status_agg','status_ben']),
     'network':('dea_network, inputs(x) intermediates(z) outputs(y) dmu(id) method(independent)', ['efficiencies','status','score_quality','score_reason','native_efficiencies','independent_efficiencies','independent_quality','independent_reason','native_independent','native_lambda_stage1','native_lambda_stage2','lambda_stage1','lambda_stage2'])
    }
    body='input double id x z y\n1 1 1 1\n2 2 1 1\n3 4 1 1\nend\n'
    def block(suffix):
     text=''
     for name,(cmd,matrices) in groups.items():
      text+=cmd+f' gen({name}{suffix}) matrix\n'
      for key in matrices:text+=f'mata: savebin("{name}_{key}_{suffix}.bin",st_matrix("r({key})"))\n'
     return text
    body+=block('before')
    body+='dea_fuzzy, inputs(x) outputs(y) dmu(id) alpha(0 .5 1) gen(f) matrix\n'
    body+=block('after')
    body+='dea_fuzzy, inputs(x) outputs(y) dmu(id) alpha(0 .5 1) gen(f2) matrix\nassert f_pess_a100==f2_pess_a100 & f_opt_a100==f2_opt_a100\n'
    def verify(d):
     count=0
     for name,(_,matrices) in groups.items():
      for key in matrices:
       before=readbin(d/f'{name}_{key}_before.bin');after=readbin(d/f'{name}_{key}_after.bin')
       assert before.size>0 and np.array_equal(before,after,equal_nan=True),(name,key)
       count+=1
     return dict(matrices=count)
    run('shared_loader_interleave',body,check=verify)

def near_cases(run):
    records=json.loads((REPO/'tests/fixtures/fuzzy_center.json').read_text())['profiles']
    for row in records:
        names=row['columns'];data=np.array(row['data']);m=row['m'];s=row['s']
        body='input double ('+' '.join(names)+')\n'+'\n'.join(' '.join(format(v,'.17g') for v in a) for a in data)+'\nend\n'
        xs=' '.join(f'x{i+1}' for i in range(m));ys=' '.join(f'y{i+1}' for i in range(s))
        bounds=' '.join(f'{side}{end}('+ ' '.join(f'{pre}_{v}' for v in (xs if side=='x' else ys).split())+')' for side,end,pre in [('x','lo','lo'),('x','hi','hi'),('y','lo','lo'),('y','hi','hi')])
        body+=f'dea_fuzzy, inputs({xs}) outputs({ys}) {bounds} dmu(id) rts({row["rts"]}) ort({row["orientation"]}) alpha(0 .001 .004 .37 .8 1) gen(res) matrix\n'
        body+='assert r(n_center_certified)==5 & r(n_center_unverified)==0 & r(lp_missing)==0\n'
        body+='matrix D=r(center_diagnostics)\nmatrix M=(r(n_center_recovered),r(center_attempts),r(lp_infeasible),r(lp_unbounded))\n'
        body+='mata: savebin("diag.bin",st_matrix("D"));savebin("meta.bin",st_matrix("M"));savebin("pess.bin",st_matrix("r(eff_pess)"));savebin("opt.bin",st_matrix("r(eff_opt)"))\n'
        body+='mata: savebin("data.bin",st_data(.,tokens("'+' '.join(names)+'")))\n'
        def check(d,row=row,data=data):
            assert np.array_equal(readbin(d/'data.bin'),data)
            P=readbin(d/'pess.bin');O=readbin(d/'opt.bin');D=readbin(d/'diag.bin');M=readbin(d/'meta.bin').ravel()
            assert np.isfinite(P).all() and np.isfinite(O).all() and np.array_equal(P[:,-1],O[:,-1])
            assert np.all(D[:,1]==2) and np.all(D[:,2]==0) and M[1]==5
            missing=~np.isfinite(D[:,4])|~np.isfinite(D[:,10]);assert np.array_equal(D[:,3],missing.astype(float)) and M[0]==sum(missing)
            gaps=[float(abs(F(float(v))-F(q))/F(q)) for v,q in zip(P[:,-1],row['center_scores'])];assert max(gaps)<2e-16
            return dict(certified_rows=5,recovered_rows=int(M[0]),max_relative_export_error=max(gaps))
        run('near_'+row['name'],body,check=check)

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--only');a=ap.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir();hashes={};rows=[]
    for n in package_files(REPO):
        raw=(REPO/n).read_bytes();(pkg/Path(n).name).write_bytes(raw);hashes[n]=hashlib.sha256(raw).hexdigest()
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,stata=a.stata),indent=2)+'\n')
    def run(name,body,source=pkg,check=None):
        if a.only and not re.search(a.only,name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset processors 1\nset type double\nset linesize 200\n'
        for k in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/k).mkdir();head+=f'sysdir set {k} "{d/k}"\n'
        head+=f'adopath ++ "{source}"\n'+WRITER;(d/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n');start=time.monotonic();record=dict(test=name,passed=False)
        try:
            p=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=900)
            log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else ''
            if p.returncode or not log:(d/'process.json').write_text(json.dumps(dict(returncode=p.returncode,stdout=p.stdout.decode(errors='replace'),stderr=p.stderr.decode(errors='replace'))))
            assert p.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2500:]
            if check:record.update(check(d) or {})
            record['passed']=True
        except Exception as e:record['error']=repr(e)
        record['seconds']=round(time.monotonic()-start,3);rows.append(record);(out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(record),flush=True)
    if not a.only or re.search(a.only,'exact_grid'):
        body,jobs=grid_body();run('exact_grid',body,check=lambda d:check_grid(d,jobs))
    near_cases(run);faults(run,pkg,out);contracts(run);interleave(run)
    summary=dict(total=len(rows),passed=sum(q['passed'] for q in rows),failed=sum(not q['passed'] for q in rows))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)
    return 0 if rows and not summary['failed'] else 1

if __name__=='__main__':raise SystemExit(main())
