#!/usr/bin/env python3
"""Sequential SW2007 invariants, legacy reproduction, syntax and public backend.

Optional public-reference tests require installed simarwilson/teradial. On Apple
Silicon use --intel for the published x86_64 plugin. No installed files are edited.
"""
from _support import package_files
import argparse, hashlib, json, subprocess, time, re
from pathlib import Path
from test_output_protection import unchanged

REPO=Path(__file__).resolve().parents[1]
DATA='''set type double
set seed 7623
set obs 80
gen id=_n
gen x=1+runiform()*4
gen z=runiform()*2-1
gen delta=1.5+.2*z+.3*invnormal(normal((1-1.5-.2*z)/.3)+runiform()*(1-normal((1-1.5-.2*z)/.3)))
gen y=x/delta
gen x2=1+runiform()*2
gen y2=y*(.8+.4*runiform())
gen str12 unrelated="unchanged"
'''


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata',required=True); ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--baseline-ado',type=Path);ap.add_argument('--external',action='store_true');ap.add_argument('--align-boundary',action='store_true');ap.add_argument('--intel',action='store_true')
    ap.add_argument('--baseline-procedure',choices=['default','joint'],default='default',
                    help='Use joint when the frozen baseline already defaults to sequential')
    ap.add_argument('--only',help='Run case names matching this regular expression')
    args=ap.parse_args(); out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    pkg=out/'ado';pkg.mkdir()
    names=package_files(REPO)
    hashes={}
    for name in names:
        data=(REPO/name).read_bytes();(pkg/Path(name).name).write_bytes(data);hashes[name]=hashlib.sha256(data).hexdigest()
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,baseline=str(args.baseline_ado),baseline_procedure=args.baseline_procedure,intel=args.intel),indent=2))
    reference_dir=None
    if args.external:
        public=Path.home()/'Library/Application Support/Stata/ado/plus/s/simarwilson.ado'
        original=public.read_text()
        (out/'reference_sha256.txt').write_text(hashlib.sha256(public.read_bytes()).hexdigest()+'\n'+str(public)+'\n')
        if args.align_boundary:
            import difflib
            reference_dir=out/'aligned_reference';reference_dir.mkdir()
            needle="replace `__deascore' = (`__deascore')^(`expinvert')"
            assert original.count(needle)==1
            modified=original.replace(needle,needle+"\n                replace `__deascore' = 1 if abs(`__deascore'-1)<=1e-7")
            (reference_dir/'simarwilson.ado').write_text(modified)
            (out/'reference_boundary_alignment.patch').write_text(''.join(difflib.unified_diff(original.splitlines(True),modified.splitlines(True))))
    results=[]
    def run(name,body,source=pkg,external=False):
        if args.only and not re.search(args.only,name) and not (name.startswith('legacy_control_') and re.search(args.only,name.replace('legacy_control_','legacy_'))):return out/name
        work=out/name;work.mkdir();head='version 16\nclear all\nset more off\nset linesize 120\n'
        if not external:
            for p in ['PLUS','PERSONAL','SITE','OLDPLACE']:
                d=work/p;d.mkdir();head+=f'sysdir set {p} "{d}"\n'
        head+=f'adopath ++ "{source}"\n'
        if external and reference_dir:head+=f'adopath ++ "{reference_dir}"\n'
        (work/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n')
        cmd=(['arch','-x86_64'] if args.intel else [])+[args.stata,'-b','do','test.do']
        start=time.monotonic()
        try:
            subprocess.run(cmd,cwd=work,timeout=240,capture_output=True)
            log=(work/'test.log').read_text(errors='replace')
            row=dict(test=name,passed='TEST_PASS' in log.splitlines(),seconds=round(time.monotonic()-start,2))
            if not row['passed']:row['tail']=log[-2500:]
        except Exception as exc:row=dict(test=name,passed=False,error=str(exc))
        results.append(row);(out/'results.json').write_text(json.dumps(results,indent=2))
        print(json.dumps(row),flush=True);return work

    if not args.external:
        for rts in ['crs','vrs','drs','irs']:
            for ori in ['input','output']:
                cmd=f'dea_sw2007 x=y, envars(z) rts({rts}) ort({ori}) reps(30) bcreps(40) seed(127) gen(a)'
                body=DATA+cmd+'''
assert e(N)==80 & e(N_trunc)==80 & e(bc_success)==40 & e(boot_success)==30
assert e(dea_reps)==40+e(bc_failed)
assert "`e(procedure)'"=="sequential"
assert e(sample)
assert a_eff_bc>0 & a_eff_bc<1
matrix point=e(b)
matrix initial=e(b_initial)
matrix B=e(b_reps)
matrix V=e(V)
matrix C=e(ci_percentile)
mata: assert(max(abs(st_matrix("V")-variance(st_matrix("B"))))<1e-12)
mata: assert(all(st_matrix("C")[1,.]:<=st_matrix("C")[2,.]))
mata: assert(max(abs(st_matrix("point")-st_matrix("initial")))>1e-5)
gen invbc=1/a_eff_bc
qui truncreg invbc z, ll(1) iterate(200)
mata: assert(max(abs(st_matrix("point")-st_matrix("e(b)")))<1e-6)
'''+cmd.replace('gen(a)','gen(b) procedure(sequential)')+'''
assert a_eff==b_eff & a_eff_bc==b_eff_bc
mata: assert(max(abs(st_matrix("B")-st_matrix("e(b_reps)")))==0)
'''
                run(f'stages_{rts}_{ori}',body)
        for form in ['x=y if inlist(id,1,3,4) | id>=10,','if id!=2 & id<=70, inputs(x) outputs(y)', 'y in 4/75, inputs(x)']:
            qualifier=form[form.find('if '):form.rindex(',')] if 'if ' in form else 'in 4/75'
            run('subset_'+str(len(results)),DATA+f'dea_sw2007 {form} envars(z) rts(crs) reps(20) bcreps(30) seed(712) gen(a)\n'+'''
matrix B=e(b_reps)
keep '''+qualifier+'''
dea_sw2007 x=y, envars(z) rts(crs) reps(20) bcreps(30) seed(712) gen(b)
assert a_eff==b_eff & a_eff_bc==b_eff_bc
mata: assert(max(abs(st_matrix("B")-st_matrix("e(b_reps)")))==0)
''')
        run('panel',DATA+'''expand 2
bysort id: gen year=_n
replace y=y*1.1 if year==2
dea_sw2007 x=y, envars(z) panel(id) time(year) rts(crs) reps(20) bcreps(30) seed(127) gen(a)
assert e(N)==160 & e(N_trunc)==160
assert !missing(a_eff_bc)
''')
        run('suppress_output_only',DATA+'''dea_sw2007 x=y, envars(z) rts(crs) reps(20) bcreps(30) seed(127) gen(a)
matrix B=e(b_reps)
gen b_eff_bc=991
dea_sw2007 x=y, envars(z) rts(crs) reps(20) bcreps(30) seed(127) gen(b) nobcorrected
assert b_eff_bc==991
assert a_eff==b_eff
mata: assert(max(abs(st_matrix("B")-st_matrix("e(b_reps)")))==0)
''')
        run('tau_and_precomputed',DATA+'''dea_sw2007 x=y, envars(z) rts(crs) reps(20) bcreps(30) seed(127) gen(a)
matrix B=e(b_reps)
dea_sw2007 y, inputs(x) scores(a_eff) envars(z) rts(crs) reps(20) bcreps(30) seed(127) gen(b)
assert a_eff_bc==b_eff_bc
mata: assert(max(abs(st_matrix("B")-st_matrix("e(b_reps)")))==0)
dea_sw2007 x=y, envars(z) rts(crs) tau(2) reps(20) bcreps(30) seed(127) gen(c)
assert abs(a_eff_bc-c_eff_bc)<1e-5
mata: assert(max(abs(2*st_matrix("B")-st_matrix("e(b_reps)")))<1e-5)
''')
        for ori in ['input','output']:
            for eng in ['simplex','ipm']:
                run(f'core_{ori}_{eng}',DATA+f'''replace z=. in 2
set seed 127
dea x=y if id<=75, envars(z) bootstrap reps1(50) reps2(30) engine({eng}) rts(crs) ort({ori}) dmu(id) gen(a) saving("saved.dta")
assert r(N)==74 & r(bootstrap_reps2)==30
assert "`r(bootstrap_procedure)'"=="sw2007_sequential"
matrix B=e(b_reps)
assert missing(a_bc) if id==2 | id>75
assert !missing(a_bc) if id!=2 & id<=75
set seed 127
dea_sw2007 x=y if id<=75, scores(a) envars(z) reps(30) bcreps(50) rts(crs) ort({ori}) gen(b)
assert abs(a_bc-b_eff_bc)<1e-12 if id!=2 & id<=75
mata: assert(max(abs(st_matrix("B")-st_matrix("e(b_reps)")))<1e-12)
use "saved.dta",clear
confirm variable a_bc
''')
        run('independent_counts_and_level',DATA+'''dea_sw2007 x=y, envars(z) rts(crs) reps(20) bcreps(30) seed(127) gen(a)
matrix B=e(b_reps)
matrix P=e(b)
matrix C=e(ci_percentile)
dea_sw2007 x=y, envars(z) rts(crs) reps(40) bcreps(30) seed(127) gen(b)
assert a_eff_bc==b_eff_bc
mata: assert(max(abs(st_matrix("P")-st_matrix("e(b)")))==0)
mata: assert(max(abs(st_matrix("B")-st_matrix("e(b_reps)")[1..20,.]))==0)
dea_sw2007 x=y, envars(z) rts(crs) reps(20) bcreps(30) seed(127) level(80) gen(c)
mata: assert(max(abs(st_matrix("B")-st_matrix("e(b_reps)")))==0)
mata: assert(all(st_matrix("C")[1,.]:<=st_matrix("e(ci_percentile)")[1,.]))
mata: assert(all(st_matrix("C")[2,.]:>=st_matrix("e(ci_percentile)")[2,.]))
''')
        run('default_replications',DATA+'''dea_sw2007 x=y, envars(z) rts(crs) seed(127) gen(a)
assert e(reps)==2000 & e(bcreps)==1000
assert e(boot_success)==2000 & e(bc_success)==1000
''')
        run('core_reference_saving_replace',DATA+'''set seed 127
dea y if id<=70, inputs(x) envars(z) bootstrap reps1(50) reps2(20) dmu(id) rts(vrs) engine(auto) reference saveref gen(a) saving("result.dta")
assert r(N)==70 & r(bootstrap_reps2)==20
local refvars "`r(refvars)'"
matrix L=r(lambda)
matrix B=e(b_reps)
mkmat `refvars' if id<=70,matrix(R)
mata: assert(max(abs(st_matrix("L")-st_matrix("R")))<1e-12)
gen old_bc=a_bc
set seed 127
dea x=y if id<=70, envars(z) bootstrap reps1(50) reps2(20) dmu(id) rts(vrs) engine(auto) reference saveref gen(a) replacevars
assert old_bc==a_bc
mata: assert(max(abs(st_matrix("B")-st_matrix("e(b_reps)")))==0)
use "result.dta",clear
assert _N==70
confirm variable a_bc
''')
        for opt in ['procedure(bad)','reps(.)','bcreps(1)','tau(.)','tolerance(-1)','level(100)','backend(simarwilson) procedure(joint)','backend(simarwilson) tau(2)']:
            run('invalid_'+str(len(results)),DATA+unchanged('dea_sw2007 x=y, envars(z) '+opt,198))
        for change,rc in [('keep in 1/2',2001),('gen year=cond(id<=2,1,2)',2001),('replace x=0 in 1',459),('gen theta=1.2',459)]:
            opts='panel(id) time(year)' if 'year=' in change else ('scores(theta)' if 'theta=' in change else '')
            run('bad_data_'+str(len(results)),DATA+change+'\n'+unchanged('dea_sw2007 x=y, envars(z) '+opts,rc))
        run('tiny_level_and_tail',DATA+'''dea_sw2007 x=y, envars(z) rts(crs) reps(3) bcreps(10) seed(127) level(.1)
mata: assert(all(st_matrix("e(ci_percentile)")[1,.]:<=st_matrix("e(ci_percentile)")[2,.]))
capture program drop _dea_sw2007_seq
quietly findfile _dea_sw2007_seq.ado
quietly run `"`r(fn)'"'
gen byte sample=1
gen mu=-40
gen draw=.
mata: _sw07_seq_draw("mu",1,1,"sample","draw")
assert draw>1 & draw<.
''')
        if args.baseline_ado:
            for ori in ['input','output']:
                for reps,bcreps in [(20,30),(30,20)]:
                    cmd=f'dea_sw2007 x=y, envars(z) rts(crs) ort({ori}) reps({reps}) bcreps({bcreps}) seed(127) gen(a)'
                    baseline_cmd=cmd+(' procedure(joint)' if args.baseline_procedure=='joint' else '')
                    old=run(f'legacy_control_{ori}_{reps}',DATA+baseline_cmd+'\nestimates save expected,replace\nkeep a_eff a_eff_bc\nsave expected,replace\n',args.baseline_ado.resolve())
                    run(f'legacy_{ori}_{reps}',DATA+cmd+' procedure(joint)\n'+f'''cf a_eff a_eff_bc using "{old}/expected.dta", all
matrix point=e(b)
matrix B=e(b_boot)
matrix C=e(ci_lo)
estimates use "{old}/expected.ster"
mata: assert(max(abs(st_matrix("point")-st_matrix("e(b)")))==0)
mata: assert(max(abs(st_matrix("B")-st_matrix("e(b_boot)")))==0)
mata: assert(max(abs(st_matrix("C")-st_matrix("e(ci_lo)")))==0)
''')
    else:
        for rts in ['crs','vrs','nirs']:
            for ori in ['input','output']:
                inv='invert' if ori=='input' else ''
                compare='''assert gap<1e-5
mata: assert(max(abs(st_matrix("e(b)")-st_matrix("extb")))<1e-5)
mata: assert(max(abs(st_matrix("e(V)")-st_matrix("extv")))<1e-5)
mata: assert(max(abs(sort(st_matrix("e(b_reps)"),1)-sort(extB,1)))<1e-4)
'''
                if rts=='vrs' and ori=='input' and not args.align_boundary:
                    compare='''assert abs(1/a_eff-ext)<1e-12
count if a_eff==1 & ext>1
assert r(N)>0
assert gap<.03
di "BOUNDARY_SAMPLE_DIFFERENCE_CONFIRMED"
'''
                run(f'public_{rts}_{ori}',DATA+f'''set seed 127
simarwilson (y=x) z, algorithm(2) rts({rts}) base({ori}) {inv} reps(100) bcreps(100) noprint te(ext) tebc(extbc) saveall(extB)
matrix extb=e(b)
matrix extv=e(V)
dea_sw2007 x=y, envars(z) rts({rts}) ort({ori}) reps(100) bcreps(100) seed(127) gen(a)
gen gap=abs(1/a_eff_bc-extbc)
quietly sum gap
di "SCORE_GAP=" %21.14g r(max)
mata: printf("POINT_GAP=%g\\n", max(abs(st_matrix("e(b)")-st_matrix("extb"))))
{compare}
set seed 127
dea_sw2007 x=y, envars(z) rts({rts}) ort({ori}) reps(100) bcreps(100) backend(simarwilson)
mata: assert(max(abs(st_matrix("e(b)")-st_matrix("extb")))<1e-10)
''',external=True)
        run('external_scores',DATA+'''gen theta=1/delta
set seed 127
simarwilson delta z, algorithm(1) nounit reps(50) noprint
matrix expected=e(b)
matrix V=e(V)
dea_sw2007, scores(theta) envars(z) backend(simarwilson) reps(50) seed(127)
mata: assert(max(abs(st_matrix("expected")-st_matrix("e(b)")))<1e-6)
mata: assert(max(abs(st_matrix("V")-st_matrix("e(V)")))<1e-6)
''',external=True)
    print(json.dumps(dict(total=len(results),passed=sum(r['passed'] for r in results))))
    raise SystemExit(any(not r['passed'] for r in results))

if __name__=='__main__':main()
