#!/usr/bin/env python3
"""Cooper source-example components: exact faces and current public results.

Requires Stata16+ and NumPy. Optional R uses lpSolve on the same LPs. No Git,
downloads, private inputs or general-purpose exact solver are required.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,csv,hashlib,json,shutil,struct,subprocess,time
import numpy as np
from _support import package_files
from test_core_retries import WRITER,readbin
from cooper_certified.table import evaluate,case_data

REPO=Path(__file__).resolve().parents[1]

def identical(value,expected):
    assert struct.pack('<d',float(value))==struct.pack('<d',float(expected))

def check_main(d,industry,fixture,records):
    C=readbin(d/'congestion.bin');status=readbin(d/'status.bin').ravel()
    years=readbin(d/'years.bin').ravel()
    assert C.shape==(17,2) and np.isfinite(C).all()
    assert np.array_equal(years,np.arange(1981,1998)) and np.array_equal(status,np.zeros(17))
    table=np.array(fixture['datasets'][industry]['manuscript_table'],dtype=float)
    actual={(1981+i,j) for i,j in zip(*np.where(abs(C-table)>.00500001))}
    expected={(c['year'],['labor','capital'].index(c['input'])) for c in fixture['cases'] if c['industry']==industry}
    assert actual==expected
    for case in fixture['cases']:
        if case['industry']!=industry:continue
        record=next(r for r in records if r['id']==case['id'] and r['representation']=='normalized_binary64')
        identical(C[case['year']-1981,['labor','capital'].index(case['input'])],F(record['component']))
    return dict(main_DMU=17,manuscript_components=34,printed_interval_differences=len(actual))

def compare_R(executable,out,models):
    d=out/'R';d.mkdir();expected={}
    for p in models:
        A=np.array(p['A'],dtype=float);b=np.array(p['b'],dtype=float);c=np.array(p['c'],dtype=float)
        rowscale=np.maximum(np.max(abs(A),axis=1),abs(b));rowscale[rowscale==0]=1
        cscale=max(1.,max(abs(c)));A/=rowscale[:,None];b/=rowscale;c/=cscale
        expected[p['id']]=float(p['objective'])/cscale
        for name,value in [('A',A),('b',b),('c',c)]:np.savetxt(d/(p['id']+'_'+name+'.csv'),value,delimiter=',',fmt='%.17g')
    (d/'ids.txt').write_text('\n'.join(p['id'] for p in models)+'\n')
    code='''suppressPackageStartupMessages(library(lpSolve))
a<-commandArgs(TRUE);d<-a[1];ids<-readLines(file.path(d,"ids.txt"));results<-list()
for(id in ids){
 A<-as.matrix(read.csv(file.path(d,paste0(id,"_A.csv")),header=FALSE))
 b<-scan(file.path(d,paste0(id,"_b.csv")),quiet=TRUE,sep=",")
 c<-scan(file.path(d,paste0(id,"_c.csv")),quiet=TRUE,sep=",")
 f<-lp("min",c,A,rep("=",length(b)),b,transpose.constraints=TRUE)
 stopifnot(f$status==0)
 results[[length(results)+1]]<-data.frame(id=id,objective=f$objval,
 residual=max(abs(A%*%f$solution-b)),min_x=min(f$solution))
}
write.csv(do.call(rbind,results),file.path(d,"results.csv"),row.names=FALSE)
cat("R_TABLE_PASS",length(ids),"lpSolve",as.character(packageVersion("lpSolve")),"\\n")
'''
    (d/'compare.R').write_text(code)
    p=subprocess.run([executable,str(d/'compare.R'),str(d)],capture_output=True,text=True,timeout=120)
    (d/'console.txt').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr
    with (d/'results.csv').open() as f:rows=list(csv.DictReader(f))
    assert len(rows)==len(models)==132 and {r['id'] for r in rows}==set(expected)
    gap=max(abs(float(r['objective'])-expected[r['id']])/(1+abs(expected[r['id']])) for r in rows)
    residual=max(float(r['residual']) for r in rows);minimum=min(float(r['min_x']) for r in rows)
    assert gap<1e-7 and residual<1e-7 and minimum>=-1e-9
    return dict(R_LPs=len(rows),R_max_relative_objective_gap=gap,R_max_scaled_primal_residual=residual)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata',required=True);ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--rscript');a=ap.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir()
    hashes={}
    for name in package_files(REPO):
        raw=(REPO/name).read_bytes();(pkg/Path(name).name).write_bytes(raw);hashes[name]=hashlib.sha256(raw).hexdigest()
    (out/'manifest.json').write_text(json.dumps(hashes,indent=2)+'\n')
    fixture=json.loads((REPO/'tests/fixtures/cooper_table.json').read_text())
    records,models,rejected=evaluate(fixture)
    (out/'exact-results.json').write_text(json.dumps(records,indent=2)+'\n')
    runs=[]
    def run(name,body,check):
        d=out/name;d.mkdir();header='version 16\nclear all\nset more off\nset type double\nset processors 1\nset linesize 180\n'
        for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:
            (d/key).mkdir();header+=f'sysdir set {key} "{d/key}"\n'
        header+=f'adopath ++ "{pkg}"\n'+WRITER
        (d/'test.do').write_text(header+body+'\ndi "COOPER_TABLE_PASS"\nexit,clear\n')
        started=time.monotonic();p=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=1200)
        log=(d/'test.log').read_text(errors='replace')
        assert p.returncode==0 and 'COOPER_TABLE_PASS' in log.splitlines(),log[-3000:]
        row=dict(name=name,passed=True,seconds=round(time.monotonic()-started,3),**check(d))
        runs.append(row);(out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n');print(json.dumps(row),flush=True)
    def data_body(industry):
        ds=fixture['datasets'][industry]
        return 'input double(year labor capital output)\n'+'\n'.join(str(y)+' '+' '.join(row) for y,row in zip(ds['years'],ds['data']))+'\nend\n'
    for industry in fixture['datasets']:
        command=('dea labor capital=output,model(congestion)' if industry=='textile' else 'dea_congestion,inputs(labor capital) outputs(output)')
        body=data_body(industry)+command+' congmethod(cooper) congweight(raw) engine(simplex) congmatrix dmu(year) generate(cp)\n'
        body+='assert r(N)==17 & r(n_failed)==0\nassert cp_status==0\nassert "`r(actual_engine)\'"=="simplex"\n'
        body+='mata: savebin("congestion.bin",st_matrix("r(input_congestion)")); savebin("status.bin",st_data(.,"cp_status")); savebin("years.bin",st_data(.,"year"))\n'
        run('main_'+industry,body,lambda d,industry=industry:check_main(d,industry,fixture,records))
    for case in fixture['cases']:
        # Target only the evaluated row, retaining all17reference rows. This
        # calls the unchanged internal pipeline, not a restricted public sample.
        k=case['year']-1980;j=['labor','capital'].index(case['input'])
        body=data_body(case['industry'])+f'run "{pkg/"_dea_cooper_data.ado"}"\n'
        body+='mata: _dea_rec_abort(); _rsm_audit_reset(0); _rsm_success_reset(1,0); _dea_link_reset(2,0)\n'
        body+='mata:\nX=st_data(.,("labor","capital"));Y=st_data(.,"output")\n'
        body+='c=cs=ct=sx=sy=lo=hi=J(1,0,.);E=U=J(0,0,.);b=ub=J(0,1,.)\n'
        body+=f'_c53_model(X,Y,{k},"raw",c,E,b,U,ub,lo,hi,cs,ct,sx,sy)\n'
        for variable in ['E','b','cs','ct','lo','hi']:
            body+=f'savebin("{variable}.bin",{variable})\n'
        body+=f'L=_c53_pipeline(c,E,b,U,ub,lo,hi,cs,ct,17,2,1,{k},J(1,17,.),J(1,17,.),"stage1",16,sx,sy,1)\n'
        body+='assert(L.point.usable==1 & L.complete==1 & L.conversion==0)\n'
        body+='savebin("state.bin",L.state); savebin("pmin.bin",L.pmin); savebin("pmax.bin",L.pmax)\n'
        body+='savebin("cmin.bin",L.cmin); savebin("cmax.bin",L.cmax); savebin("congestion.bin",L.point.congestion)\nend\n'
        def check(d,case=case,j=j):
            data,units,weights=case_data(fixture,case,'normalized_binary64')
            base=next(p for p in models if p['id']==case['id']+'_normalized_binary64_radial')
            assert np.array_equal(readbin(d/'E.bin'),np.array(base['A'],dtype=float))
            assert np.array_equal(readbin(d/'b.bin').ravel(),np.array(base['b'],dtype=float)+np.array(base['A'],dtype=float)[:,-1])
            assert np.array_equal(readbin(d/'cs.bin').ravel(),np.array([0]*17+weights+[0],dtype=float))
            assert np.array_equal(readbin(d/'ct.bin').ravel(),np.array([float.fromhex(x) for x in case['normalized_technical_objective']]))
            assert np.array_equal(readbin(d/'lo.bin').ravel(),np.array([0]*20+[1]))
            hi=readbin(d/'hi.bin').ravel();assert np.array_equal(hi[:19],np.array([1]*17+data[case['year']-1981][:2],dtype=float)) and np.isnan(hi[19:]).all()
            assert np.array_equal(readbin(d/'state.bin').ravel(),[1,1,2,1,6,4])
            ref=next(r for r in records if r['id']==case['id'] and r['representation']=='normalized_binary64')
            for key in ['pmin','pmax']:
                for got,want in zip(readbin(d/(key+'.bin')).ravel(),ref['projection']):identical(got,F(want))
            for key in ['cmin','cmax','congestion']:identical(readbin(d/(key+'.bin')).ravel()[j],F(ref['component']))
            return dict(reference_rows=17,target_rows=1,projection_endpoints=6,congestion_endpoints=4,checked_component_endpoints=2)
        run('ranges_'+case['id'],body,check)
    rcheck=compare_R(a.rscript,out,models) if a.rscript else dict(R_requested=False,R_LPs=0)
    summary=dict(passed=True,Stata_batches=len(runs),public_main_commands=2,public_main_DMU=34,targeted_internal_pipelines=4,
        exact_model_sequences=12,exact_certificates=len(models),tampering_rejected=rejected,source_table_components=68,
        source_table_differences=4,full_reference_rows_per_industry=17,full_all_DMU_ranges_run=False,
        final_publication_verified=False,production_solver_changed=False,**rcheck)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)

if __name__=='__main__':main()
