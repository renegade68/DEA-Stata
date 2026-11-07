#!/usr/bin/env python3
"""Verify VRS matrix observation identity, duplicate labels and source-row maps.

Uses Python's standard library and licensed Stata 16+. This tests matrix
assembly and result preservation; it does not change or certify the LP kernels.
An optional pre-F08 snapshot reproduces the historical aggregation/order bugs.
"""
from _support import package_files
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import re
import shutil
import subprocess
import time

from test_saveref import BASIC, instrument
from test_output_protection import RANDOM

REPO = Path(__file__).resolve().parents[1]
DATA = '''set obs 12
gen long serial=_n
gen double x1=1+mod(serial*7,9)
gen double x2=1+mod(serial*5,11)
gen double y1=1+mod(serial*3,7)
gen double y2=1+mod(serial*2,5)
gen str20 tag=cond(mod(serial,3),"keep = this, yes","omit")
gen double id=serial
'''
PATTERNS = {
    'numeric_unique': '',
    'numeric_unsorted': 'replace id=mod(serial*7,13)+11\n',
    'numeric_duplicate': 'replace id=ceil(serial/3)\n',
    'numeric_same': 'replace id=1\n',
    'numeric_missing': 'replace id=. in 2\nreplace id=.a in 3\nreplace id=.z in 5\n',
    'numeric_fractional': 'replace id=serial/7\n',
    'numeric_large': 'replace id=1e15+serial\n',
    'numeric_value_labels': 'replace id=mod(serial,3)\nlabel define ids 0 "same" 1 "same" 2 "third"\nlabel values id ids\n',
    'string_unique': 'drop id\ngen str20 id="unit"+string(serial)\n',
    'string_duplicate': 'drop id\ngen str30 id=cond(mod(serial,2),"same unit, A","other: unit = B")\n',
    'string_missing': 'drop id\ngen str20 id="unit"\nreplace id="" in 2/4\nreplace id=" " in 5\n',
    'string_long': 'drop id\ngen str100 id="'+'a'*80+'"+string(serial)\n',
    'string_unicode': 'drop id\ngen str30 id=cond(mod(serial,2),"기관 가","기관 나")\n',
    'string_strL': 'drop id\ngen strL id="long unit "+string(mod(serial,2))\n',
    'automatic': '',
}

CAPTURE = '''matrix V=r(vrsfrontier)
matrix O=r(vrsobs)
matrix E=r(efficiency)
matrix ST=r(stage_status)
matrix SS=r(solution_source)
local n=r(N)
local ncrs=r(n_crs)
local ndrs=r(n_drs)
local nirs=r(n_irs)
local refs "`r(refvars)'"
local eng "`r(engine)'"
'''
CHECK = '''assert rowsof(V)==`n' & colsof(V)==5
assert rowsof(O)==`n' & colsof(O)==1
assert "`:colnames V'"=="CRS_TE VRS_TE DRS_TE SCALE RTS"
assert "`:colnames O'"=="obsno"
mata:
chosen=selectindex(st_data(.,"selected"))
assert(all(st_matrix("O"):==chosen))
assert(all(st_matrix("V")[.,(1,2,4,5)]:==st_data(chosen,tokens("CRS_TE VRS_TE SCALE RTS_NUM"))))
assert(max(abs(st_matrix("V")[.,2]-st_matrix("E")))<1e-7)
assert(all(st_matrixrowstripe("V"):==st_matrixrowstripe("E")))
assert(all(st_matrixrowstripe("O"):==st_matrixrowstripe("E")))
end
count if selected & RTS_NUM==0
assert r(N)==`ncrs'
count if selected & RTS_NUM==-1
assert r(N)==`ndrs'
count if selected & RTS_NUM==1
assert r(N)==`nirs'
assert missing(CRS_TE,VRS_TE,SCALE,RTS_NUM) if !selected
assert RTS=="" if !selected
'''


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--baseline-ado',type=Path);ap.add_argument('--only')
    args=ap.parse_args()
    if args.baseline_ado and not (args.baseline_ado/'dea.ado').is_file():ap.error('Invalid pre-F08 ado snapshot')
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'current_ado';pkg.mkdir();hashes={}
    files=package_files(REPO)
    for name in files:
        data=(REPO/name).read_bytes();(pkg/Path(name).name).write_bytes(data);hashes[name]=hashlib.sha256(data).hexdigest()
    retry=out/'forced_retry_ado';instrument(pkg,retry)
    fault=out/'missing_auxiliary_ado';shutil.copytree(pkg,fault)
    p=fault/'dea.ado';s=p.read_text();start=s.index('program define _dea_vrsfrontier, rclass');end=s.index('\nend\n',start)+5
    fragment=s[start:end];anchor='    // Count DMUs\n';assert fragment.count(anchor)==1
    injection='''    // TEST COPY ONLY: preserve distinct extended missings for one duplicate ID.
    foreach v in `crs_te' `vrs_te' `drs_te' `scale_eff' `rts_type' {
        quietly replace `v'=.a if `touse' & _n==2
        quietly replace `v'=.z if `touse' & _n==4
    }
'''
    p.write_text(s[:start]+fragment.replace(anchor,injection+anchor)+s[end:])
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,stata=args.stata,baseline=str(args.baseline_ado)),indent=2))
    results=[]

    def run(name,script,source=pkg,force=False):
        if not force and args.only and not re.search(args.only,name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 120\n'
        for folder in ['PLUS','PERSONAL','SITE','OLDPLACE']:
            (d/folder).mkdir();head+=f'sysdir set {folder} "{d/folder}"\n'
        head+=f'adopath ++ "{source.resolve()}"\n';(d/'test.do').write_text(head+script+'\ndi "TEST_PASS"\n')
        start=time.monotonic();row=dict(test=name,passed=False)
        try:
            subprocess.run([args.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=240)
            log=(d/'test.log').read_text(errors='replace')
            if 'TEST_PASS' not in log.splitlines():raise RuntimeError(log[-2200:])
            row['passed']=True
        except Exception as e:row['error']=repr(e)
        row['seconds']=round(time.monotonic()-start,3);results.append(row)
        (out/'results.json').write_text(json.dumps(results,indent=2));print(json.dumps(row),flush=True)

    for i,(pattern,engine,ori,stage) in enumerate(itertools.product(PATTERNS,['simplex','ipm','auto'],['input','output'],[1,2])):
        name=f'grid_{pattern}_{engine}_{ori}_{stage}'
        setup=DATA+PATTERNS[pattern]
        if i%3==0:setup+='gsort -serial\n'
        if i%4==3:setup+='replace x1=. in 3\nreplace y2=.a in 8\n'
        q=['','if tag=="keep = this, yes"','in 2/10','if mod(serial,3)!=0 in 2/10'][i%4]
        cond=['1','tag=="keep = this, yes"','inrange(_n,2,10)','mod(serial,3)!=0 & inrange(_n,2,10)'][i%4]
        setup+=f'gen byte selected=({cond}) & !missing(x1,x2,y1,y2)\ngen long control_id=_n\n'
        prefix=f'dea x1 x2=y1 y2 {q},' if i%2 else f'dea y1 y2 {q}, inputs(x1 x2)'
        options=f'rts(vrs) engine({engine}) orientation({ori}) stage({stage})'
        refs=i%3==0;refopts=' reference saveref' if refs else ''
        setup+='quietly ds\nlocal original_vars "`r(varlist)\'"\ntempfile original\nsave `original\'\npreserve\n'
        setup+=prefix+f' {options} dmu(control_id) gen(control)'+(' reference' if refs else '')+'\n'
        setup+='matrix C=r(vrsfrontier)\nmatrix CE=r(efficiency)\nlocal ccrs=r(n_crs)\nlocal cdrs=r(n_drs)\nlocal cirs=r(n_irs)\n'
        if refs:setup+='matrix CL=r(lambda)\n'
        setup+='restore\n'
        if i%5==0:
            setup+='foreach v in res CRS_TE VRS_TE SCALE RTS_NUM RTS {\n gen str20 `v\'="old output"\n}\n'
        idopt='' if pattern=='automatic' else ' dmu(id)'
        setup+=prefix+f' {options}{idopt} gen(res){refopts} saving("saved.dta")'+(' replacevars' if i%5==0 else '')+'\n'
        setup+=CAPTURE
        if refs:setup+='matrix L=r(lambda)\nmata: assert(all(st_matrix("L"):==st_matrix("CL")))\n'
        setup+=CHECK+'mata: assert(all(st_matrix("V"):==st_matrix("C")))\nmata: assert(all(st_matrix("E"):==st_matrix("CE")))\n'
        setup+='assert `ncrs\'==`ccrs\' & `ndrs\'==`cdrs\' & `nirs\'==`cirs\'\ncf `original_vars\' using `original\', all\n'
        if refs:
            setup+='mata: assert(all(st_matrix("L"):==st_data(selectindex(st_data(.,"selected")),tokens(st_local("refs")))))\n'
            setup+='local j=0\nforeach v of local refs {\n local ++j\n assert "`:char `v\'[dea_peer_index]\'"=="`j\'"\n}\n'
        saveid='_dmu' if pattern=='automatic' else 'id'
        setup+=f'preserve\nkeep if selected\nkeep {saveid} x1 x2 y1 y2 res `refs\'\ntempfile expected\nsave `expected\'\nuse "saved.dta",clear\nmata: assert(all(st_data(.,tokens("_dea_radial_status _dea_stage2_status _dea_score_source _dea_ref_stage")):==(st_matrix("ST"),st_matrix("SS"))))\ndrop _dea_radial_status _dea_stage2_status _dea_score_source _dea_ref_stage\ncf _all using `expected\', all\nrestore\n'
        setup+='preserve\nclear\nsvmat double V,names(col)\nsvmat double O,names(col)\nexport delimited using "matrix.csv",replace\nrestore\n'
        run(name,setup)

    for engine,ori in itertools.product(['simplex','ipm','auto'],['input','output']):
        for kind,change in [('duplicate','replace id=1 in 2\n'),('unsorted','replace id=50 in 1\nreplace id=10 in 2\nreplace id=30 in 3\nreplace id=20 in 4\n')]:
            name=f'original_{kind}_{engine}_{ori}'
            cmd=f'dea x=y, rts(vrs) engine({engine}) orientation({ori}) dmu(id) gen(res)\n'
            run(name,BASIC+change+cmd+'gen byte selected=1\n'+CAPTURE+CHECK)
            if args.baseline_ado:
                old=BASIC+change+cmd+'matrix V=r(vrsfrontier)\nmatrix E=r(efficiency)\nassert !missing(res)\n'
                old+='mata: assert(all(st_matrix("V")[4,.]:==0))\n' if kind=='duplicate' else 'mata: assert(max(abs(st_matrix("V")[.,2]-st_matrix("E")))>.1)\n'
                run('old_'+name,old,args.baseline_ado)
        s=BASIC+'replace id=1 in 2\n'+f'dea x=y, rts(vrs) engine({engine}) orientation({ori}) dmu(id) gen(res)\n'
        s+='matrix V=r(vrsfrontier)\nmatrix O=r(vrsobs)\nmata: assert(all(st_matrix("V")[2,.]:==.a))\nmata: assert(all(st_matrix("V")[4,.]:==.z))\nmata: assert(all(st_matrix("O"):==(1::4)))\nassert !missing(res)\n'
        run(f'auxiliary_missing_{engine}_{ori}',s,fault)

    for ori,reference in itertools.product(['input','output'],[False,True]):
        opts=f'rts(vrs) dmu(id) orientation({ori})'
        s=BASIC+'replace id=1 in 2\n'+f'dea x=y, {opts} engine(ipm) gen(control)\nmatrix C=r(vrsfrontier)\nmatrix CO=r(vrsobs)\nlocal failed0=r(n_failed)\nmatrix CST=r(stage_status)\n'
        s+=f'dea x=y, {opts} engine(auto) gen(res) replacevars '+('reference saveref ' if reference else '')+'saving("saved.dta")\n'
        s+='assert "`r(engine)\'"=="ipm"\nassert r(n_failed)==`failed0\'\nmatrix ST=r(stage_status)\nmata: assert(all(st_matrix("ST"):==st_matrix("CST")))\nmatrix V=r(vrsfrontier)\nmatrix O=r(vrsobs)\nmata: assert(all(st_matrix("V"):==st_matrix("C")))\nmata: assert(all(st_matrix("O"):==st_matrix("CO")))\nassert res==control\n'
        run(f'retry_{ori}_{reference}',s,retry)

    for engine,ori in itertools.product(['simplex','ipm','auto'],['input','output']):
        # Verify r-class preservation across bootstrap and saving after VRS.
        s=BASIC+'replace id=1 in 2/3\ngen byte selected=_n>1\nset seed 551\n'
        s+=f'dea x=y if selected, rts(vrs) engine({engine}) orientation({ori}) dmu(id) gen(res) bootstrap reps1(50) saving("saved.dta")\n'
        s+='matrix B=r(bootstrap)\n'+CAPTURE+CHECK+'assert rowsof(B)==`n\'\n'
        run(f'bootstrap_{engine}_{ori}',s)
    s=RANDOM+'gen byte selected=id<=70\nset seed 551\n'
    s+='dea y if selected, inputs(x) envars(z) bootstrap reps1(50) reps2(20) dmu(id) rts(vrs) engine(auto) reference saveref gen(res) saving("saved.dta")\n'
    run('environmental_bootstrap',s+CAPTURE+CHECK)
    s=BASIC+'dea x=y, dmu(id) rts(vrs) gen(first)\nmatrix O=r(vrsobs)\ndea x=y, dmu(id) rts(crs) gen(second)\ncapture matrix list r(vrsobs)\nassert _rc==111\n'
    run('no_stale_map_after_crs',s)
    if args.baseline_ado:
        # Bit-exact preservation of scores/generated decomposition variables
        # against the pre-F08 package, independent of the new row-label contract.
        for engine,ori,stage in itertools.product(['simplex','ipm','auto'],['input','output'],[1,2]):
            name=f'history_{engine}_{ori}_{stage}'
            cmd=f'dea x1 x2=y1 y2 if serial!=3 in 2/11, dmu(id) rts(vrs) engine({engine}) orientation({ori}) stage({stage}) gen(res)\n'
            variables='res CRS_TE VRS_TE SCALE RTS_NUM RTS'
            old=DATA+'replace id=mod(serial*7,13)+11\n'+cmd+'keep '+variables+'\nsave "expected.dta",replace\n'
            run('old_'+name,old,args.baseline_ado,force=not args.only or bool(re.search(args.only,name)))
            path=out/('old_'+name)/'expected.dta'
            s=DATA+'replace id=mod(serial*7,13)+11\n'+cmd+f'cf {variables} using "{path}", all\n'
            run(name,s)
    summary=dict(total=len(results),passed=sum(r['passed'] for r in results),failed=sum(not r['passed'] for r in results))
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
    raise SystemExit(summary['failed']>0)


if __name__=='__main__':main()
