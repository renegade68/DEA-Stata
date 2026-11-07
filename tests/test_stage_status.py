#!/usr/bin/env python3
"""R05 TE status/provenance, save-only diagnostics, and unchanged numerical results.
Requires Stata 16+, NumPy/SciPy, and the immediately preceding dea 2.0.51 ado snapshot.
Fault injection is confined to copied packages; solver tests also live in
 test_simplex_stages.py and test_ipm_stages.py.
"""
from _support import package_files
from pathlib import Path
import argparse, hashlib, itertools, json, re, shutil, subprocess, time
import numpy as np
from test_simplex_stages import DATA, FIXTURES, INIT, instrument
from test_core_retries import WRITER, readbin

REPO=Path(__file__).resolve().parents[1]
COLS='_dea_radial_status _dea_stage2_status _dea_score_source _dea_ref_stage'


def diagnostics(stage, saved=True, sample=''):
    s=f'''matrix ST=r(stage_status)
matrix SS=r(solution_source)
matrix E=r(efficiency)
scalar nf=r(n_failed)
scalar nr=r(n_radial_failed)
scalar ns=r(n_stage2_failed)
local final_engine "`r(engine)'"
local refs "`r(refvars)'"
local sc : colnames ST
assert "`sc'"=="radial stage2"
local sc : colnames SS
assert "`sc'"=="radial reference"
mata:
S=st_matrix("ST");P=st_matrix("SS");E=st_matrix("E")
assert(cols(S)==2 & cols(P)==2 & rows(S)==rows(E) & rows(P)==rows(E))
assert(all(S[.,1]:==(E:>=.)))
assert(all((P[.,1]:==0):==(E:>=.)))
assert(all((P[.,1]:==0):|(P[.,1]:==1):|(P[.,1]:==2)))
assert(st_numscalar("nr")==sum(E:>=.))
if ({stage}==1) {{
 assert(all(S[.,2]:>=.) & st_numscalar("ns")==0)
 assert(all(P[.,2]:==(P[.,1]:==1)))
 assert(st_numscalar("nf")==sum(E:>=.))
}}
else {{
 assert(all(P[.,2]:==(2:*(S[.,2]:==0))))
 assert(st_numscalar("ns")==sum(S[.,2]:!=0))
 assert(st_numscalar("nf")==sum((E:>=.):|(S[.,2]:!=0)))
}}
end
capture matrix list r(lambda)
if !_rc {{
 matrix L=r(lambda)
 mata: assert(all((rowmissing(st_matrix("L")):==cols(st_matrix("L"))):==(st_matrix("SS")[.,2]:==0)))
}}
foreach v in {COLS} {{
 capture confirm variable `v', exact
 assert _rc==111
}}
'''
    if saved:
        s+=f'''preserve
use "saved.dta",clear
assert _N==rowsof(ST)
mata: assert(all(st_data(.,tokens("{COLS}")):==(st_matrix("ST"),st_matrix("SS"))))
if "`refs'"!="" {{
 mata: assert(all(st_data(.,tokens(st_local("refs"))):==st_matrix("L")))
}}
restore
'''
    return s


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stata',required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--baseline-ado',type=Path,required=True);p.add_argument('--only')
    args=p.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    pkg=out/'ado';pkg.mkdir();hashes={}
    for n in package_files(REPO):
        raw=(REPO/n).read_bytes();(pkg/Path(n).name).write_bytes(raw);hashes[n]=hashlib.sha256(raw).hexdigest()
    baseline=args.baseline_ado.resolve();audit=out/'instrumented_ado';instrument(pkg,audit)
    oldaudit=out/'instrumented_before';instrument(baseline,oldaudit)
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,baseline_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in baseline.iterdir() if p.is_file()}),indent=2)+'\n')
    results=[]
    def run(name,body,source=pkg,check=None):
        if args.only and not re.search(args.only,name):return None
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 160\n'
        for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:
            (d/key).mkdir();head+=f'sysdir set {key} "{d/key}"\n'
        head+=f'adopath ++ "{source}"\n'+INIT+WRITER
        (d/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n');start=time.monotonic();row=dict(test=name,passed=False)
        try:
            r=subprocess.run([args.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=240)
            log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else r.stderr.decode(errors='replace')
            assert r.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-3500:]
            if check:row.update(check(d))
            row['passed']=True
        except Exception as e:row['error']=repr(e)
        row['seconds']=round(time.monotonic()-start,3);results.append(row)
        (out/'results.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(row),flush=True)
        return d
    def export(stage):
        mats=['efficiency','lambda']+(['slack_in','slack_out'] if stage==2 else [])
        return ''.join(f'mata: savebin("{m}.bin",st_matrix("r({m})"))\n' for m in mats)
    # The old defect is reproduced separately, not mistaken for naturally failing data.
    for ori in ['input','output']:
        s=DATA+'scalar fault_stage2=1\n'+f'dea x1 x2=y1 y2,dmu(id) engine(simplex) ort({ori}) reference\n'
        s+='assert r(n_failed)==0\nmata: assert(!hasmissing(st_matrix("r(lambda)")) & !hasmissing(st_matrix("r(slack_in)")))\n'
        s+='capture matrix list r(stage_status)\nassert _rc==111\n'
        run(f'reproduce_before_{ori}',s,oldaudit)
    for fixture,engine,rts,ori,stage in itertools.product(['slack','scaled'],['simplex','ipm','auto'],['crs','vrs','drs','irs'],['input','output'],[1,2]):
        name=f'history_{fixture}_{engine}_{rts}_{ori}_{stage}'
        if args.only and not re.search(args.only,name):continue
        body=FIXTURES[fixture]+f'dea x1 x2=y1 y2,dmu(id) engine({engine}) rts({rts}) ort({ori}) stage({stage}) reference gen(e)\n'
        old=run(name+'_before',body+export(stage),baseline)
        def check(d,old=old,stage=stage):
            assert (d/'efficiency.bin').read_bytes()==(old/'efficiency.bin').read_bytes()
            st=readbin(d/'status.bin');good=st[:,1]==0 if stage==2 else np.ones(len(st),dtype=bool)
            exact=1
            for m in ['lambda']+(['slack_in','slack_out'] if stage==2 else []):
                a=readbin(d/f'{m}.bin');b=readbin(old/f'{m}.bin')
                np.testing.assert_array_equal(a[good],b[good])
                if stage==2:assert np.isnan(a[~good]).all()
                if (d/f'{m}.bin').read_bytes()==(old/f'{m}.bin').read_bytes():exact+=1
            return dict(exact_numeric_files=exact,projection_rows_compared=int(good.sum()),unavailable_rows=int((~good).sum()))
        run(name+'_after',body+export(stage)+'mata: savebin("status.bin",st_matrix("r(stage_status)"))\n'+diagnostics(stage,False),check=check)
    # Both syntaxes, aliases, selected rows, duplicate/string IDs, and save/no-ref paths.
    for form,engine,rts,ori,stage,ref in itertools.product(['equal','standard'],['simplex','ipm','auto'],['crs','vrs','drs','irs','nirs','ndrs'],['input','output'],[1,2],['','reference','saveref']):
        sample='if id>=2 in 2/7';io=f'x1 x2=y1 y2 {sample},' if form=='equal' else f'y1 y2 {sample}, inputs(x1 x2)'
        refopt="reference saveref" if ref=="saveref" else ref
        body=DATA+'gen str8 labelid="same"\n'
        body+=f'dea {io} dmu(labelid) engine({engine}) rts({rts}) ort({ori}) stage({stage}) gen(e) {refopt} saving("saved.dta")\n'
        body+=diagnostics(stage)+'assert rowsof(ST)==6\nassert missing(e) if _n==1 | _n==8\nassert sentinel=="preserve me"\n'
        run(f'grid_{form}_{engine}_{rts}_{ori}_{stage}_{ref or "none"}',body)
    # Mixed outcomes in one sample: normal, dual recovery, stage-2 failure,
    # radial failure, and a second stage-2 failure. No double-counting of rows.
    partial=out/'mixed_fault_ado';shutil.copytree(audit,partial)
    f=partial/'_dea_simplex.ado';t=f.read_text()
    for flag,condition in [('primal','(k==2 | k==4)'),('dual','k==4'),('stage2','(k==3 | k==5)')]:
        t=t.replace(f'st_numscalar("fault_{flag}")',f'(st_numscalar("fault_{flag}") & {condition})')
    f.write_text(t)
    for engine,rts,ori,stage in itertools.product(['simplex','auto'],['crs','vrs','drs','irs'],['input','output'],[1,2]):
        body=DATA
        cmd=f'dea x1 x2=y1 y2,dmu(id) rts({rts}) ort({ori}) stage({stage}) reference'
        if engine=='auto':
            body+=cmd+' engine(ipm) gen(control)\n'
            body+='matrix CS=r(stage_status)\nmatrix CP=r(solution_source)\nmatrix CL=r(lambda)\n'
        body+='scalar fault_primal=1\nscalar fault_dual=1\nscalar fault_stage2=1\n'
        body+=cmd+f' engine({engine}) gen(e) saveref saving("saved.dta") replacevars\n'+diagnostics(stage)
        if engine=='simplex':
            body+='mata: assert(all(st_matrix("SS")[.,1]:==(1\\2\\1\\0\\1\\1\\1\\1)))\n'
            body+=f'assert nr==1 & ns=={0 if stage==1 else 3} & nf=={1 if stage==1 else 3}\n'
            refstage='1\\0\\1\\0\\1\\1\\1\\1' if stage==1 else '2\\2\\0\\0\\0\\2\\2\\2'
            body+=f'mata: assert(all(st_matrix("SS")[.,2]:==({refstage})))\n'
        else:
            body+='assert "`final_engine\'"=="ipm"\nassert e==control\n'
            body+='mata: assert(all(st_matrix("ST"):==st_matrix("CS")) & all(st_matrix("SS"):==st_matrix("CP")) & all(st_matrix("L"):==st_matrix("CL")))\n'
        run(f'mixed_{engine}_{rts}_{ori}_{stage}',body,partial)
    # Small positive radial scores must survive a stage-2 failure (R07).
    for ori in ['input','output']:
        values='1 1 1\n2 1e9 1' if ori=='input' else '1 1 1e9\n2 1 1'
        body='input double id x y\n'+values+'\nend\n'
        for stage in [1,2]:
            body+=f'scalar fault_stage2={stage-1}\n'
            body+=f'dea x=y,dmu(id) engine(simplex) ort({ori}) stage({stage}) tolerance(1e-12) gen(e{stage}) reference\n'
        body+='assert e1==e2\nassert e1[2]>0 & abs(e1[2]-1e-9)<1e-19\nassert r(n_stage2_failed)==2\n'
        run(f'near_zero_{ori}',body,audit)
    # Failure metadata must survive the radial bootstrap, without changing its purpose.
    for engine,ori,pf in itertools.product(['simplex','auto'],['input','output'],[0,1]):
        s=DATA+f'scalar fault_primal={pf}\nscalar fault_stage2=1\nset seed 5297\n'
        s+=f'dea x1 x2=y1 y2,dmu(id) engine({engine}) ort({ori}) stage(2) gen(e) reference saveref bootstrap reps1(60) saving("saved.dta")\n'
        s+='assert r(n_failed)==8 & r(n_radial_failed)==0 & r(n_stage2_failed)==8\nassert r(bootstrap_reps)==60\nassert !missing(e_bc,e_se,e_lo,e_hi)\n'
        s+=diagnostics(2)
        run(f'bootstrap_{engine}_{ori}_{pf}',s,audit)
    # A saved diagnostic name may not overwrite selected scientific data or the score.
    for form,name,role in itertools.product(['equal','standard'],COLS.split(),['input','output','id','score']):
        body=DATA;inputs='x1 x2';outputs='y1 y2';idvar='id';score='keepme'
        if role=='input':body+=f'rename x1 {name}\n';inputs=f'{name} x2'
        elif role=='output':body+=f'rename y1 {name}\n';outputs=f'{name} y2'
        elif role=='id':body+=f'rename id {name}\n';idvar=name
        else:score=name
        body+=f'gen double {score}=777\nquietly datasignature\nlocal sig "`r(datasignature)\'"\n'
        io=f'{inputs}={outputs},' if form=='equal' else f'{outputs},inputs({inputs})'
        body+=f'capture noisily dea {io} dmu({idvar}) engine(simplex) gen({score}) saving("saved.dta") replacevars\nassert _rc==110\n'
        body+='quietly datasignature\nassert "`r(datasignature)\'"=="`sig\'"\ncapture confirm file "saved.dta"\nassert _rc==601\n'
        run(f'collision_{form}_{name}_{role}',body)
    # Unretained source columns are left intact in memory; saved diagnostics are fresh.
    s=DATA+''.join(f'gen double {v}=77\n' for v in COLS.split())
    s+='dea x1 x2=y1 y2,dmu(id) engine(simplex) gen(e) saving("saved.dta")\n'
    s+='matrix ST=r(stage_status)\nmatrix SS=r(solution_source)\n'
    s+=''.join(f'assert {v}==77\n' for v in COLS.split())
    s+=f'use "saved.dta",clear\nmata: assert(all(st_data(.,tokens("{COLS}")):==(st_matrix("ST"),st_matrix("SS"))))\n'
    run('unretained_source_columns',s)
    # Non-TE calls must not inherit metadata from a preceding TE call.
    for model,extra in [('supereff',''),('multiplier',''),('te','ordinal(x1)'),('te','fdh')]:
        s=DATA+'dea x1 x2=y1 y2,dmu(id) engine(simplex) gen(e) reference\n'
        s+=f'dea x1 x2=y1 y2,dmu(id) engine(simplex) gen(f) model({model}) {extra} saving("saved.dta")\n'
        s+='capture matrix list r(stage_status)\nassert _rc==111\ncapture matrix list r(solution_source)\nassert _rc==111\nuse "saved.dta",clear\n'
        s+=''.join(f'capture confirm variable {v}\nassert _rc==111\n' for v in COLS.split())
        run(f'no_stale_{model}_{extra.replace("(","_").replace(")","")}',s)
    summary=dict(total=len(results),passed=sum(r['passed'] for r in results),historical_pairs=sum('exact_numeric_files' in r for r in results),exact_numeric_files=sum(r.get('exact_numeric_files',0) for r in results))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
    raise SystemExit(summary['total']!=summary['passed'])

if __name__=='__main__':main()
