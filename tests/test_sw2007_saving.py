#!/usr/bin/env python3
"""Verify SW2007 observation files, estimation invariance and write failures.

Uses the real installed simarwilson for --external; --intel runs its published
x86_64 plugin on Apple Silicon. Private fault injection never edits installed code.
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

from test_sw2007 import DATA
from test_output_protection import unchanged
from test_slack_names import EXPORTER

REPO=Path(__file__).resolve().parents[1]
QUALS=[('', '1'),('if id<=75','id<=75'),('in 3/76','inrange(_n,3,76)'),
       ('if id!=5 in 2/75','id!=5 & inrange(_n,2,75)')]
FIXTURE=DATA+'''gen double theta=1/delta
gen str20 sid="Unit "+string(id)
char unrelated[audit] "caller metadata"
label variable x "Original input"
char x[audit] "input metadata"
label data "Caller data"
'''

CAPTURE='''local scalars : e(scalars)
local matrices : e(matrices)
local macros : e(macros)
local scalars : list sort scalars
local matrices : list sort matrices
local macros : list sort macros
local rngstate `"`c(rngstate)'"'
foreach s of local scalars {
    tempname exs_`s'
    scalar `exs_`s''=e(`s')
}
foreach s of local matrices {
    tempname exm_`s'
    matrix `exm_`s''=e(`s')
    local rn_`s'=ustrregexra("`:rowfullnames e(`s')'","__[0-9A-Za-z]{6}","TEMP")
    local cn_`s'=ustrregexra("`:colfullnames e(`s')'","__[0-9A-Za-z]{6}","TEMP")
}
foreach s of local macros {
    local ex_`s'=ustrregexra(`"`e(`s')'"',"__[0-9A-Za-z]{6}","TEMP")
}
tempvar esample_before
gen byte `esample_before'=e(sample)
'''
VERIFY='''local got : e(scalars)
local got : list sort got
assert "`got'"=="`scalars'"
local got : e(matrices)
local got : list sort got
assert "`got'"=="`matrices'"
local got : e(macros)
local got : list sort got
assert "`got'"=="`macros'"
assert `"`c(rngstate)'"'==`"`rngstate'"'
assert e(sample)==`esample_before'
foreach s of local scalars {
    assert e(`s')==scalar(`exs_`s'')
}
foreach s of local matrices {
    mata: assert(all(st_matrix("e(`s')"):==st_matrix("`exm_`s''")))
    assert ustrregexra("`:rowfullnames e(`s')'","__[0-9A-Za-z]{6}","TEMP")=="`rn_`s''"
    assert ustrregexra("`:colfullnames e(`s')'","__[0-9A-Za-z]{6}","TEMP")=="`cn_`s''"
}
foreach s of local macros {
    assert ustrregexra(`"`e(`s')'"',"__[0-9A-Za-z]{6}","TEMP")==`"`ex_`s''"'
}
'''

def dump():
    return '''local mats : e(matrices)
foreach m of local mats {
    mata: r02_dump("matrix_`m'.bin",st_matrix("e(`m')"))
}
local scalars : e(scalars)
foreach m of local scalars {
    mata: r02_dump("scalar_`m'.bin",st_numscalar("e(`m')"))
}
'''

def command(route='sequential',rts='crs',ori='output',form='A',q='',nobc=False,panel=False,prefix='a'):
    if form=='A':io=f'x=y {q},'
    elif form=='B':io=f'y {q}, inputs(x)'
    else:io=f'{q}, inputs(x) outputs(y) scores(theta)'
    if route=='scoreonly':io=f'{q}, scores(theta)'
    backend='simarwilson' if route in ['external','scoreonly'] else 'internal'
    proc=route if backend=='internal' else 'sequential'
    return f'dea_sw2007 {io} envars(z) dmu(sid) rts({rts}) orientation({ori}) reps(12) bcreps(15) seed(127) gen({prefix}) backend({backend}) procedure({proc})'+(' nobcorrected' if nobc else '')+(' panel(id) time(year)' if panel else '')

def case(route,rts,ori,form,index,nobc=False,panel=False,old=False):
    q,selected=QUALS[index%4];external=route in ['external','scoreonly']
    data=FIXTURE
    if index%7==0:data+='gsort -id\n'
    if nobc:
        stale='a_delta_bc' if external else 'a_eff_bc'
        data+=f'gen str6 {stale}="stale"\n'
    if panel:
        data+='expand 2\nbysort id: gen year=_n\nreplace y=1.1*y if year==2\n'
        q='if id<=75';selected='id<=75'
    # Unrelated missing values do not reduce the estimation or saved sample.
    data+='gen double unused=.\n'
    if index%5==0:
        data+='replace z=. in 4\n'
    selected=f'({selected}) & !missing(z)'
    s=data+f'gen byte selected=({selected})\n'
    s+='quietly ds\nlocal originals "`r(varlist)\'"\ntempfile original\nsave `original\'\n'
    cmd=command(route,rts,ori,form,q,nobc,panel)
    filename='results = selected, 한글.dta'
    internal=not external
    outs='a_eff'+(' a_eff_bc' if not nobc else '') if internal else 'a_delta'+(' a_delta_bc' if route=='external' and form!='C' and not nobc else '')
    keep='sid '+('id year ' if panel else '')+('' if route=='scoreonly' else 'x y ')+('theta ' if route=='scoreonly' or form=='C' else '')+'z'
    if old:
        s+=cmd+f' saving("{filename}")\n'+dump()+f'capture confirm file "{filename}"\nassert _rc==601\n'
        if internal:s+='mata: r02_dump("data_outputs.bin",st_data(.,tokens("'+outs+'")))\n'
        return s
    s+=cmd+'\n'+CAPTURE
    if internal:
        s+='tempfile expected_outputs\nquietly save `expected_outputs\'\n'
    s+=cmd+f' replacevars saving("{filename}")\n'+VERIFY+dump()
    if internal:
        s+='cf '+outs+' using `expected_outputs\',all\n'
        s+='mata: r02_dump("data_outputs.bin",st_data(.,tokens("'+outs+'")))\n'
    s+='cf `originals\' using `original\',all\nassert "`:char unrelated[audit]\'"=="caller metadata"\n'
    s+='assert "`:char _dta[dea_sw2007_backend]\'"==""\n'
    if external:
        s+='capture confirm variable a_delta, exact\nassert _rc==111\n'
        # Compare the exported scores to the real public command, including
        # its boundary/quantile conventions. No reference source is patched.
        if route=='scoreonly' or form=='C':
            s+='gen double a_delta=1/theta if selected\n'
        else:
            inv='invert' if ori=='input' else ''
            s+=f'tempvar oracle_bc\nset seed 127\nsimarwilson (y=x) z if selected, algorithm(2) rts({rts}) base({ori}) {inv} reps(12) bcreps(15) level(95) noprint tename(a_delta) tebc(`oracle_bc\')\n'
            if not nobc:s+='rename `oracle_bc\' a_delta_bc\n'
    s+='preserve\nkeep if selected\nkeep '+keep+' '+outs+'\norder '+keep+' '+outs+'\ntempfile expected\nsave `expected\'\n'
    s+='foreach v of varlist '+keep+''' {
 local type_`v' : type `v'
 local format_`v' : format `v'
 local label_`v' : variable label `v'
}
'''
    s+=f'use "{filename}",clear\ncf _all using `expected\',all\nquietly ds\nassert "`r(varlist)\'"=="{keep} {outs}"\n'
    s+='foreach v of varlist '+keep+''' {
 assert "`:type `v''"=="`type_`v''"
 assert "`:format `v''"=="`format_`v''"
 assert "`:variable label `v''"=="`label_`v''"
}
'''
    backend='internal' if internal else 'simarwilson';proc=route if internal else 'sequential'
    s+=f'assert "`:char _dta[dea_sw2007_backend]\'"=="{backend}"\nassert "`:char _dta[dea_sw2007_procedure]\'"=="{proc}"\n'
    scale='efficiency' if internal else 'inverse_efficiency'
    s+=f'assert "`:char _dta[dea_sw2007_score_scale]\'"=="{scale}"\nrestore\n'
    return s

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata',required=True);ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--baseline-ado',type=Path);ap.add_argument('--external',action='store_true')
    ap.add_argument('--intel',action='store_true');ap.add_argument('--only')
    a=ap.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    pkg=out/'ado';pkg.mkdir();hashes={}
    for name in package_files(REPO):
        raw=(REPO/name).read_bytes();(pkg/Path(name).name).write_bytes(raw);hashes[name]=hashlib.sha256(raw).hexdigest()
    broken=out/'write_failure_ado';shutil.copytree(pkg,broken)
    p=broken/'dea_sw2007.ado';code=p.read_text();needle='    quietly save `"`saving\'"\', `replace\''
    assert code.count(needle)==1;p.write_text(code.replace(needle,'    exit 603 // TEST ONLY: a final disk-write failure'))
    reference={}
    if a.external:
        plus=Path.home()/'Library/Application Support/Stata/ado/plus'
        for rel in ['s/simarwilson.ado','s/simarwilson.sthlp','t/teradial.ado','t/te_radial.plugin']:
            f=plus/rel;reference[str(f)]=hashlib.sha256(f.read_bytes()).hexdigest()
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,reference_sha256=reference,
        baseline=str(a.baseline_ado),intel=a.intel,stata=a.stata),indent=2)+'\n')
    rows=[]
    def run(name,body,source=pkg,external=False,compare=None):
        if a.only and not re.search(a.only,name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 160\n'
        if not external:
            for sub in ['PLUS','PERSONAL','SITE','OLDPLACE']:
                (d/sub).mkdir();head+=f'sysdir set {sub} "{d/sub}"\n'
        head+=f'adopath ++ "{source.resolve()}"\n'+EXPORTER
        (d/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n')
        cmd=(['arch','-x86_64'] if external and a.intel else [])+[a.stata,'-b','do','test.do']
        row=dict(test=name,passed=False);start=time.monotonic()
        try:
            proc=subprocess.run(cmd,cwd=d,capture_output=True,timeout=240)
            log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else proc.stderr.decode(errors='replace')
            assert proc.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2000:]
            if compare:
                previous={p.name:p.read_bytes() for p in compare.glob('*.bin')}
                current={p.name:p.read_bytes() for p in d.glob('*.bin')}
                assert previous and current==previous,'Historical numeric results differ'
                row['exact_numeric_files']=len(previous)
            if (d/'protected_copy.dta').exists():
                assert (d/'protected_copy.dta').read_bytes()==(d/'existing.dta').read_bytes(),'Protected file changed'
                row['protected_file_bytes_unchanged']=True
            row['passed']=True
        except Exception as e:row['error']=repr(e)
        row['seconds']=round(time.monotonic()-start,3);rows.append(row)
        (out/'results.json').write_text(json.dumps(rows,indent=2)+'\n')
        if not row['passed']:print(json.dumps(row),flush=True)
    plans=list(itertools.product(['sequential','joint'],['crs','vrs','drs','irs'],['input','output'],['A','B','C'],[False,True],[False]))
    plans+=list(itertools.product(['sequential','joint'],['crs','vrs'],['input','output'],['A'],[False,True],[True]))
    if a.external:
        plans+=list(itertools.product(['external'],['crs','vrs','nirs'],['input','output'],['A','B','C'],[False,True],[False]))
        plans+=list(itertools.product(['scoreonly'],['crs'],['input','output'],['C'],[False,True],[False]))
    for i,(route,rts,ori,form,nobc,panel) in enumerate(plans):
        name=f'grid_{route}_{rts}_{ori}_{form}_{int(nobc)}_{int(panel)}';external=route in ['external','scoreonly']
        kwargs=dict(route=route,rts=rts,ori=ori,form=form,index=i,nobc=nobc,panel=panel)
        prior=None
        if a.baseline_ado:
            prior=out/(name+'_before');run(name+'_before',case(**kwargs,old=True),a.baseline_ado,external)
        run(name,case(**kwargs),external=external,compare=prior)
    routes=['sequential','joint']+(['external','scoreonly'] if a.external else [])
    for route in routes:
        external=route in ['external','scoreonly'];cmd=command(route)
        # Failures after file preflight must neither publish a file nor replace
        # an existing one, even when replace was explicitly authorized.
        failures=[('empty_sample','',command(route,q='if id<0'),2000),
                  ('invalid_score','replace theta=1.2 in 1\n',command(route,form='C'),459)]
        if external:
            failures.append(('unsupported_rts','',command(route,rts='irs'),198))
        else:
            failures.append(('short_sample','keep in 1/2\n',cmd,2001))
        for label,setup,bad,rc in failures:
            s=FIXTURE+setup+'save "existing.dta"\ncopy "existing.dta" "protected_copy.dta"\n'
            s+=unchanged(bad+' saving("existing.dta") replace',rc)
            run(f'abort_{route}_{label}',s,external=external)
        for suffix in ['.dta','']:
            # Existing files are checked before seed(), estimation or data edits.
            s=FIXTURE+'save "existing.dta"\nregress y z\n'+CAPTURE
            s+=unchanged(cmd+f' saving("existing{suffix}") replacevars',602)+VERIFY
            s+='drop `esample_before\'\ncf _all using "existing.dta",all\n'
            run(f'protect_file_{route}_{len(suffix)}',s,external=external)
        for filename in ['fresh results','fresh.xyz','dir.with.dot/한글 result']:
            path=filename+'.dta' if '.' not in Path(filename).name else filename
            s=FIXTURE+('mkdir "dir.with.dot"\n' if '/' in filename else '')
            s+=cmd+f' saving("{filename}")\nconfirm file "{path}"\n'
            s+='preserve\nuse "'+path+'",clear\nassert _N==80\nrestore\n'
            run(f'path_{route}_{filename.split("/")[-1].split(" ")[0]}',s,external=external)
        s=FIXTURE+unchanged(cmd+' saving("missing parent/results.dta")',603)
        run(f'bad_parent_{route}',s,external=external)
        s=FIXTURE+unchanged(cmd+' saving("   ")',198)
        run(f'empty_path_{route}',s,external=external)
        s=FIXTURE+'save "existing.dta"\n'
        if not external:s+='gen str8 a_eff="old"\ngen str8 a_eff_bc="old"\n'
        s+=cmd+' saving("existing") replace\nconfirm file "existing.dta"\n'
        s+='assert _N==80\nassert unrelated=="unchanged"\n'
        run(f'replace_file_{route}',s,external=external)
        s=FIXTURE+'save "existing.dta"\n'+unchanged(cmd+' saving("existing.dta") replace',603)
        s+='cf _all using "existing.dta",all\nassert "`e(cmd)\'"==""\n'
        run(f'write_failure_{route}',s,broken,external)
        if external:
            s=FIXTURE+'gen str20 a_delta="caller column"\n'+cmd+' saving("external.dta")\nassert a_delta=="caller column"\n'
            s+='preserve\nuse "external.dta",clear\nconfirm numeric variable a_delta\nrestore\n'
            run(f'external_unrelated_name_{route}',s,external=True)
        else:
            # Identifier roles may name the same variable; save it once.
            s=FIXTURE+'expand 2\nbysort id: gen year=_n\n'
            s+=command(route,panel=True).replace('dmu(sid)','dmu(id)')+' saving("panel.dta")\n'
            s+='preserve\nuse "panel.dta",clear\nquietly ds\nassert "`r(varlist)\'"=="id year x y z a_eff a_eff_bc"\nrestore\n'
            run(f'duplicate_id_roles_{route}',s)
        for prefix in ['abcdefghijklmnopqrstuvwxyz12345','sid']:
            if prefix=='sid':
                target='sid_delta' if external else 'sid_eff'
                data=FIXTURE+f'rename x {target}\n'
                bad=command(route,prefix=prefix).replace('x=y',target+'=y')
                if route=='scoreonly':bad=bad.replace('scores(theta)',f'scores(theta) inputs({target})')
                rc=110
            else:data=FIXTURE;bad=command(route,prefix=prefix);rc=7
            run(f'name_guard_{route}_{prefix}',data+unchanged(bad+' saving("guard.dta") replace',rc)+'capture confirm file "guard.dta"\nassert _rc==601\n',external=external)
    summary=dict(total=len(rows),passed=sum(r['passed'] for r in rows),failed=[r['test'] for r in rows if not r['passed']],
        historical_pairs=sum('exact_numeric_files' in r for r in rows))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
    raise SystemExit(bool(summary['failed']))

if __name__=='__main__':main()
