#!/usr/bin/env python3
"""Exercise output-name preflight, explicit replacement, and unchanged normal results.

Each case uses a fresh Stata process and isolated ado paths. A successful rejection
case means the documented error occurred AND the entire dataset stayed unchanged.
"""
from _support import package_files, require_git_history
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time

REPO=Path(__file__).resolve().parents[1]
DATA='''input id firm t x y
1 1 1 1 1
2 2 1 2 3
3 3 1 4 4
4 1 2 1 1.2
5 2 2 2 3.6
6 3 2 4 4.8
end
gen double z=(x+y)/2
gen double w=1
gen double p=10
gen double bad=1
gen double lo=.9*x
gen double hi=1.1*x
gen str16 note="keep me"
'''
RANDOM='''set obs 80
set seed 27811
gen id=_n
gen double z=rnormal()
gen double x=1+5*runiform()
gen double y=x/(2+exp(.3*z)+abs(rnormal()))
gen str16 note="keep me"
'''


def unchanged(command, rc=110):
    return '''
quietly ds
local original_vars "`r(varlist)'"
local i=0
foreach v of local original_vars {
    local ++i
    label variable `v' "Protected: `v'"
    char `v'[review] "custom metadata"
    local type_`i' : type `v'
    local format_`i' : format `v'
    local label_`i' : variable label `v'
    local vallabel_`i' : value label `v'
}
quietly datasignature
local sig_before "`r(datasignature)'"
'''+f'capture noisily {command}\nlocal actual_rc=_rc\nassert `actual_rc\'=={rc}\n'+'''
quietly datasignature
assert "`sig_before'"=="`r(datasignature)'"
quietly ds
assert "`original_vars'"=="`r(varlist)'"
local i=0
foreach v of local original_vars {
    local ++i
    assert "`:type `v''"=="`type_`i''"
    assert "`:format `v''"=="`format_`i''"
    assert "`:variable label `v''"=="`label_`i''"
    assert "`:value label `v''"=="`vallabel_`i''"
    assert "`:char `v'[review]'"=="custom metadata"
}
'''


def fuzzy_historical_check(outputs, path):
    """F23 changes precision; use bounded historical differences, not bit equality."""
    return f'''tempname observed historical
mkmat {' '.join(outputs)}, matrix(`observed')
preserve
use "{path}", clear
mkmat {' '.join(outputs)}, matrix(`historical')
restore
mata: assert(!hasmissing(st_matrix("`observed'")) & !hasmissing(st_matrix("`historical'")))
mata: assert(max(abs(st_matrix("`observed'")-st_matrix("`historical'")))<1e-6)
'''


def radial_historical_check(key, outputs, path):
    """Retain an analytic oracle and exact replacement checks while allowing
    final-bit changes from basis refactorization against the old Git version.
    """
    numeric=[v for v in outputs if v!='RTS']
    s=f'''tempname observed historical
mkmat {' '.join(numeric)}, matrix(`observed')
preserve
use "{path}", clear
mkmat {' '.join(numeric)}, matrix(`historical')
restore
mata: assert(!hasmissing(st_matrix("`observed'")) & !hasmissing(st_matrix("`historical'")))
mata: assert(max(abs(st_matrix("`observed'")-st_matrix("`historical'")))<1e-12)
'''
    if key=='core_vrs':
        # Lower convex hull of the fixed one-input/one-output DATA fixture.
        s+=f'cf RTS RTS_NUM using "{path}", all\n'
        s+='''assert abs(res-1)<1e-12 if inlist(id,1,4,5,6)
assert abs(res-((y[2]-y[4])/(y[5]-y[4])*(x[5]-x[4])+x[4])/x[2])<1e-12 if id==2
assert abs(res-((y[3]-y[5])/(y[6]-y[5])*(x[6]-x[5])+x[5])/x[3])<1e-12 if id==3
'''
    else:
        # CRS with one input/output: the maximum y/x defines the frontier.
        s+='''tempvar ratio
gen double `ratio'=y/x
quietly summarize `ratio', meanonly
local frontier=r(max)
assert abs(res_eff-`ratio'/`frontier')<1e-12
drop `ratio'
'''
    return s+'assert note=="keep me"\n'


def window_historical_check(outputs, path):
    """R11-C column scaling can change final bits. Check availability, a tight
    historical bound, and the analytic pooled CRS frontier independently.
    Replacement results are also compared exactly with a current clean run.
    """
    return f'''tempname observed historical
mkmat {' '.join(outputs)}, matrix(`observed')
preserve
use "{path}", clear
mkmat {' '.join(outputs)}, matrix(`historical')
restore
mata: assert(all((st_matrix("`observed'"):>=.):==(st_matrix("`historical'"):>=.)))
mata: assert(max(abs(st_matrix("`observed'")-st_matrix("`historical'")))<1e-12)
tempvar ratio
gen double `ratio'=y/x
quietly summarize `ratio', meanonly
local frontier=r(max)
assert abs(res_mean-`ratio'/`frontier')<1e-8
assert res_sd==0
assert abs(res_w1_p1-res_mean)<1e-12 if t==1
assert abs(res_w1_p2-res_mean)<1e-12 if t==2
assert missing(res_w1_p1) if t==2
assert missing(res_w1_p2) if t==1
drop `ratio'
assert note=="keep me"
'''


def crosseff_historical_check(outputs, path):
    """The fixed one-input/one-output CRS fixture has a unique peer ratio.
    Validate the analytic answer and exact current replacements, while allowing
    numerical improvement after R11-D1 removes the weight floor and scales LPs.
    """
    return f'''tempname observed historical
mkmat {' '.join(outputs)}, matrix(`observed')
preserve
use "{path}", clear
mkmat {' '.join(outputs)}, matrix(`historical')
restore
mata: assert(!hasmissing(st_matrix("`observed'")) & !hasmissing(st_matrix("`historical'")))
mata: assert(max(abs(st_matrix("`observed'")-st_matrix("`historical'")))<1e-9)
tempvar ratio
gen double `ratio'=y/x
quietly summarize `ratio', meanonly
local frontier=r(max)
assert abs(res-`ratio'/`frontier')<1e-10
assert abs(res_self-`ratio'/`frontier')<1e-10
assert abs(res_maverick-1)<1e-10
drop `ratio'
assert note=="keep me"
'''



def crossgame_historical_check(outputs, path):
    """Unique 1x1 CRS ratios are analytic; allow the old native solve's error."""
    s=crosseff_historical_check(outputs,path)
    s=s.replace("assert abs(res-`ratio'/`frontier')<1e-10", "assert abs(res_agg-`ratio'/`frontier')<1e-10")
    s=s.replace("assert abs(res_self-`ratio'/`frontier')<1e-10", "assert abs(res_ben-`ratio'/`frontier')<1e-10")
    s=s.replace("assert abs(res_maverick-1)<1e-10", "assert abs(res_nash-`ratio'/`frontier')<1e-10")
    return s.replace(')))<1e-9)', ')))<1e-6)')


def twostage_historical_check(key, outputs, path):
    # This fixture has one intermediate: both relational and independent CRS
    # decompositions reduce to the two analytic radial ratios.
    s=fuzzy_historical_check(outputs,path)
    s+='tempvar ratio1 ratio2 expected1 expected2\ngenerate double `ratio1\'=z/x\ngenerate double `ratio2\'=y/z\nquietly summarize `ratio1\'\ngenerate double `expected1\'=`ratio1\'/r(max)\nquietly summarize `ratio2\'\ngenerate double `expected2\'=`ratio2\'/r(max)\n'
    if key=='network':s+='assert abs(res_stage1-`expected1\')<1e-7\nassert abs(res_stage2-`expected2\')<1e-7\n'
    else:s+='assert abs(res_leader-`expected1\')<1e-7\nassert abs(res_follower-`expected2\')<1e-7\nassert res_leader==res_indep1 & res_follower==res_indep2\n'
    s+='assert abs(res_overall-`expected1\'*`expected2\')<1e-7\nassert note=="keep me"\n'
    return s


def games_historical_check(key, outputs, path):
    # Same characteristic/allocation rule; normalization and explicit LP bounds
    # change numerical last bits, and additive Shapley removes summation noise.
    names=[v for v in outputs if not v.endswith('_rank')]
    text='mata: game_current=st_data(.,tokens("'+' '.join(names)+'"))\npreserve\n'
    text+=f'use "{path}",clear\n'
    text+='mata: assert(!hasmissing(game_current) & max(abs(game_current-st_data(.,tokens("'+' '.join(names)+'"))))<1e-6)\nrestore\n'
    if key=='cooperative':
        text+='assert res_indiv==res_shapley & res_contrib==0\n'
        text+='mata: v=st_data(.,"res_shapley"); assert(all(st_data(.,"res_rank"):==_dea_games_rank(v)))\n'
    else:text+='assert res_gain==res_bargain-res_disagree & res_gain>=0\n'
    return text+'assert note=="keep me"\n'


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata',required=True)
    ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--baseline-ref',default='97cbb4c6c35be93a425934c05649f5682a99bfdd')
    ap.add_argument('--baseline-ado',type=Path,help='Compare with a reviewed source snapshot instead of the Git baseline')
    ap.add_argument('--only')
    args = ap.parse_args()
    if not args.baseline_ado:
        try:
            require_git_history(REPO, args.baseline_ref)
        except ValueError as error:
            ap.error(str(error))
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    pkg=out/'current_ado';pkg.mkdir()
    names=package_files(REPO)
    hashes={}
    for name in names:
        data=(REPO/name).read_bytes();(pkg/Path(name).name).write_bytes(data)
        hashes[name]=hashlib.sha256(data).hexdigest()
    baseline=out/'baseline_ado';baseline.mkdir()
    if args.baseline_ado:
        if not (args.baseline_ado/'dea.ado').is_file():ap.error('Invalid baseline ado directory')
        for src in args.baseline_ado.iterdir():
            if src.is_file():(baseline/src.name).write_bytes(src.read_bytes())
        baseline_label=str(args.baseline_ado.resolve())
    else:
        names=subprocess.check_output(['git','-C',str(REPO),'ls-tree','-r','--name-only',args.baseline_ref,'ado/'],text=True).splitlines()
        for name in names:
            data=subprocess.check_output(['git','-C',str(REPO),'show',args.baseline_ref+':'+name])
            (baseline/Path(name).name).write_bytes(data)
        baseline_label=args.baseline_ref
    baseline_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in baseline.iterdir()}
    (out/'manifest.json').write_text(json.dumps(dict(baseline_ref=baseline_label,baseline_sha256=baseline_hashes,source_sha256=hashes,stata=args.stata),indent=2))
    results=[]
    def run(name,body,source=pkg,force=False):
        if not force and args.only and not re.search(args.only,name):return None
        work=out/name;work.mkdir()
        head='version 16\nclear all\nset more off\nset linesize 120\n'
        for folder in ['plus','personal','site','oldplace']:
            (work/folder).mkdir();head+=f'sysdir set {folder.upper()} "{work/folder}"\n'
        head+=f'adopath ++ "{source}"\n'
        (work/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n')
        start=time.monotonic()
        try:
            proc=subprocess.run([args.stata,'-b','do','test.do'],cwd=work,capture_output=True,timeout=180)
            log=(work/'test.log').read_text(errors='replace') if (work/'test.log').exists() else proc.stderr.decode(errors='replace')
            row=dict(test=name,passed='TEST_PASS' in log.splitlines(),seconds=round(time.monotonic()-start,3))
            if not row['passed']:row['tail']=log[-2200:]
        except subprocess.TimeoutExpired:row=dict(test=name,passed=False,error='timeout')
        results.append(row);(out/'results.json').write_text(json.dumps(results,indent=2))
        print(json.dumps(row),flush=True)
        return work

    # Profile: full public output plan and a reproducible normal command.
    profiles={
      'core_te':(DATA,'dea x = y, dmu(id) gen(res)',['res']),
      'core_vrs':(DATA,'dea x = y, dmu(id) rts(vrs) gen(res)',['res','CRS_TE','VRS_TE','SCALE','RTS_NUM','RTS']),
      'core_fdh':(DATA,'dea x = y, dmu(id) model(fdh) gen(res)',['res']),
      'core_interval':(DATA,'dea x = y, dmu(id) model(imprecise) xlo(lo) xhi(hi) gen(res)',['res','res_opt']),
      'core_ordinal':(DATA,'dea x = y, dmu(id) ordinal(x) gen(res)',['res','res_opt','res_mid']),
      'core_weights':(DATA,'dea x = y, dmu(id) model(multiplier) saveweights gen(res)',['res','_dea_v1','_dea_u1','_dea_u0']),
      'core_refs':(DATA,'dea x = y, dmu(id) reference saveref gen(res)',['res']+[f'_ref_{i}' for i in range(1,7)]),
      'core_bootstrap':(DATA,'dea x = y, dmu(id) bootstrap reps1(50) gen(res)',['res','res_bc','res_bias','res_se','res_lo','res_hi']),
      'core_congestion':(DATA,'dea x = y, dmu(id) model(congestion) gen(res)',['res','res_crs','res_vrs','res_cong_eff','res_congestion','res_se']),
      'bootstrap':(DATA,'dea_bootstrap x = y, dmu(id) reps(50) seed(12) gen(res)',['res','res_bc','res_se','res_lo','res_hi']),
      'crosseff':(DATA,'dea_crosseff, inputs(x) outputs(y) dmu(id) gen(res)',['res','res_self','res_maverick']),
      'cooperative':(DATA,'dea_cooperative, inputs(x) outputs(y) dmu(id) gen(res)',['res_indiv','res_shapley','res_contrib','res_rank']),
      'crossgame':(DATA,'dea_crossgame, inputs(x) outputs(y) dmu(id) gen(res)',['res_agg','res_ben','res_nash']),
      'nash':(DATA,'dea_nash, inputs(x) outputs(y) dmu(id) gen(res)',['res_disagree','res_coop','res_bargain','res_gain']),
      'network':(DATA,'dea_network y, inputs(x) intermediates(z) dmu(id) gen(res)',['res_overall','res_stage1','res_stage2']),
      'stackelberg':(DATA,'dea_stackelberg, inputs(x) intermediates(z) outputs(y) dmu(id) gen(res)',['res_leader','res_follower','res_overall','res_indep1','res_indep2']),
      'malmquist':(DATA,'dea_malmquist y, inputs(x) dmu(firm) time(t) gen(res)',['res_ec','res_tc','res_mpi']),
      'window':(DATA,'dea_window y, inputs(x) dmu(firm) time(t) width(2) gen(res)',['res_mean','res_sd','res_w1_p1','res_w1_p2']),
      'fuzzy':(DATA,'dea_fuzzy, inputs(x) outputs(y) dmu(id) alpha(0 1) gen(res)',['res_pess_a0','res_opt_a0','res_pess_a100','res_opt_a100']),
      'congestion':(DATA,'dea_congestion, inputs(x) outputs(y) dmu(id) gen(res)',['res','res_crs','res_vrs','res_cong_eff','res_congestion','res_se']),
      'sbmundesirable':(DATA,'dea_sbmundesirable, inputs(x) outputs(y) bad(bad) dmu(id) gen(res)',['res']),
      'sw2007':(RANDOM,'dea_sw2007 y, inputs(x) envars(z) dmu(id) rts(crs) reps(10) bcreps(10) seed(553) gen(res)',['res_eff','res_eff_bc']),
    }
    for key,(fixture,cmd,outs) in profiles.items():
        selected=[outs[0]] if len(outs)==1 else [outs[0],outs[-1]]
        for target in selected:
            run(f'reject_{key}_{target}',fixture+f'gen str12 {target}="old result"\n'+unchanged(cmd))
        prep=''.join(f'gen str12 {v}="old result"\n' for v in outs)
        run(f'reject_all_{key}',fixture+prep+unchanged(cmd))
        # Renaming a model input to a planned output must be rejected even with replacement.
        target=outs[-1]
        colliding=re.sub(r'\bx\b',target,cmd)
        for flag in ['', ' replacevars']:
            run(f'protect_input_{key}_{bool(flag)}',fixture+f'rename x {target}\n'+unchanged(colliding+flag))
        idvar=re.search(r'\bdmu\(([^)]+)\)',cmd)[1]
        id_collision=re.sub(r'\bdmu\([^)]+\)',f'dmu({target})',cmd)
        run(f'protect_string_id_{key}',fixture+f'tostring {idvar}, replace\nrename {idvar} {target}\n'+unchanged(id_collision+' replacevars'))
        # Mark the sample before removing an explicitly replaced output used in if.
        qualifier=cmd.replace(',',f' if {outs[0]}>0,',1)
        run(f'replace_qualifier_{key}',fixture+f'gen {outs[0]}=1\n'+qualifier+' replacevars\n')
        # New output names may be prefixes of unrelated existing variables.
        extras=[v+'_extra' for v in outs if len(v+'_extra')<=32]
        run(f'exact_names_{key}',fixture+''.join(f'gen double {v}=991\n' for v in extras)+'set seed 883\n'+cmd+'\n'+''.join(f'assert {v}==991\n' for v in extras))
        # Normal results and replacements match the selected baseline.
        comparison_names=[f'normal_{key}',f'replace_{key}_partial',f'replace_{key}_last',f'replace_{key}_all']
        needs_control=bool(args.only and any(re.search(args.only,n) for n in comparison_names))
        # F03 replaces the broken core bootstrap corrections and repairs the
        # dedicated bootstrap's smoother and CI quantiles. Their old values are not an
        # oracle; test_core_bootstrap.py checks the new inference separately.
        compare_outs=['res'] if key in ['core_bootstrap','bootstrap'] else (['res_eff'] if key == 'sw2007' else outs)
        control=run(f'control_{key}',fixture+'set seed 883\n'+cmd+'\nkeep '+ ' '.join(compare_outs)+'\nsave "expected.dta", replace\n',baseline,force=needs_control)
        if control:
            check=f'cf {" ".join(compare_outs)} using "{control}/expected.dta", all\nassert note=="keep me"\n'
            if key == 'fuzzy':
                check=fuzzy_historical_check(outs,control/'expected.dta')+'assert note=="keep me"\n'
            if key in ['core_vrs','sw2007']:
                check=radial_historical_check(key,compare_outs,control/'expected.dta')
            if key == 'window':
                check=window_historical_check(outs,control/'expected.dta')
            if key == 'crosseff':
                check=crosseff_historical_check(outs,control/'expected.dta')
            if key in ['network','stackelberg']:
                check=twostage_historical_check(key,outs,control/'expected.dta')
            if key in ['cooperative','nash']:
                check=games_historical_check(key,outs,control/'expected.dta')
            if key == 'crossgame':
                check=crossgame_historical_check(outs,control/'expected.dta')
            if key in ['core_congestion','congestion']:
                # P03 delivers the requested native TE convergence tolerance;
                # historical last-bit equality is no longer the contract.
                check='mata: current=st_data(.,tokens("'+' '.join(outs)+'"))\npreserve\n'
                check+=f'use "{control}/expected.dta", clear\n'
                check+='mata: assert(!hasmissing(current) & max(abs(current-st_data(.,tokens("'+' '.join(outs)+'"))))<1e-6)\nrestore\nassert note=="keep me"\n'
            if key in ['core_vrs','core_bootstrap','bootstrap','sw2007','fuzzy','core_congestion','congestion','window','crosseff','crossgame','network','stackelberg','cooperative','nash']:
                current=run(f'current_inference_{key}',fixture+'set seed 883\n'+cmd+'\nkeep '+' '.join(outs)+'\nsave "expected.dta", replace\n',force=needs_control)
                if current:
                    check+=f'cf {" ".join(outs)} using "{current}/expected.dta", all\n'
            run(f'normal_{key}',fixture+'set seed 883\n'+cmd+'\n'+check)
            for mode,old in [('partial',selected[:1]),('last',[outs[-1]]),('all',outs)]:
                run(f'replace_{key}_{mode}',fixture+''.join(f'gen str12 {v}="old result"\n' for v in old)+'set seed 883\n'+cmd+' replacevars\n'+check)

    # Core parsers, engines and primary model output protection, with and without the opt-in.
    models={'te':'','fdh':'','ddf':'gx(x) gy(y)','undesirable':'bad(bad)',
            'cost':'w(w)','revenue':'p(p)','profit':'w(w) p(p) rts(vrs)',
            'allocative':'w(w)','sbm':'','supereff':'','supersbm':'','additive':'',
            'multiplier':'','imprecise':'xlo(lo) xhi(hi)','congestion':'','sbmundesirable':'bad(bad)'}
    for model,options in models.items():
        for engine in ['simplex','ipm']:
            for form in ['A','B']:
                body='dea x = y,' if form=='A' else 'dea y, inputs(x)'
                cmd=f'{body} dmu(id) model({model}) engine({engine}) {options} gen(x)'
                for flag in ['', ' replacevars']:
                    run(f'core_model_input_{model}_{engine}_{form}_{bool(flag)}',DATA+unchanged(cmd+flag))
    # Protect every model-data role, not only primary inputs.
    roles=[('id',''),('y',''),('w','model(cost) w(w)'),('p','model(revenue) p(p)'),
           ('z','bootstrap envars(z) reps1(50)'),('bad','model(sbmundesirable) bad(bad)'),
           ('lo','model(imprecise) xlo(lo) xhi(hi)'),('hi','model(imprecise) xlo(lo) xhi(hi)'),
           ('z','model(ddf) gx(z) gy(y)')]
    for i,(target,opt) in enumerate(roles):
        # The incomplete core envars path now fails option validation before
        # name preflight. Both rejection paths must preserve the full dataset.
        expected_rc=110
        run(f'protected_role_{i}',DATA+unchanged(f'dea x = y, dmu(id) gen({target}) {opt} replacevars',expected_rc))
    run('protected_environment_sw2007',RANDOM+'rename z res_eff\n'+unchanged('dea_sw2007 y, inputs(x) envars(res_eff) dmu(id) gen(res) replacevars reps(10) bcreps(10)'))
    for key,cmd in [('dea','dea x = y, gen(res)'),('fuzzy','dea_fuzzy, inputs(x) outputs(y) alpha(1) gen(res)')]:
        for flag in ['', ' replacevars']:
            run(f'auto_id_existing_{key}_{bool(flag)}',DATA+'gen str20 _dmu="user ID"\n'+unchanged(cmd+flag))
        run(f'auto_id_success_{key}',DATA+cmd+'\nconfirm variable _dmu\nassert _dmu=="dmu"+string(_n)\n')
    run('vrs_internal_duplicate_name',DATA+unchanged('dea x = y, dmu(id) rts(vrs) gen(CRS_TE) replacevars'))
    for flag in ['', ' replacevars']:
        run(f'reserved_all_{bool(flag)}',DATA+unchanged('dea x = y, dmu(id) gen(_all)'+flag,7))
    run('fuzzy_alpha_duplicate',DATA+unchanged('dea_fuzzy, inputs(x) outputs(y) dmu(id) alpha(.001 .001) gen(res) replacevars',198))
    run('long_suffix_atomic',DATA+'gen abcdefghijklmnopqrstuvwxyz12345=991\n'+unchanged('dea x = y, dmu(id) bootstrap gen(abcdefghijklmnopqrstuvwxyz12345) replacevars',7))
    run('fuzzy_partial_collision_atomic',DATA+'gen res_opt_a100=991\n'+unchanged('dea_fuzzy, inputs(x) outputs(y) dmu(id) alpha(0 1) gen(res)'))
    run('file_replace_not_variable_replace',DATA+'gen res=991\n'+unchanged('dea x = y, dmu(id) gen(res) saving("output.dta") replace')+'capture confirm file "output.dta"\nassert _rc==601\n')
    run('variable_replace_not_file_replace',DATA+'''gen res=991
save "output.dta", replace
capture noisily dea x = y, dmu(id) gen(res) saving("output.dta") replacevars
assert _rc==602
use "output.dta", clear
assert res==991
''')
    run('legacy_internal_names_preserved',DATA+'''gen __idea_eff=991
gen __idea_score=992
gen x_2=993
dea x = y, dmu(id) ordinal(x) gen(res)
assert __idea_eff==991 & __idea_score==992 & x_2==993
''')
    run('sw2007_replace_alias_protects_input',RANDOM+'rename x res_eff\n'+unchanged('dea_sw2007 y, inputs(res_eff) envars(z) dmu(id) gen(res) replace reps(10) bcreps(10)'))
    for name in ['malmquist','window']:
        fixture,cmd,outs=profiles[name]
        for role in ['firm','t']:
            target=outs[-1]
            run(f'protect_panel_{name}_{role}',fixture+f'rename {role} {target}\n'+unchanged(re.sub(r'\b'+role+r'\b',target,cmd)+' replacevars'))
    run('inactive_weights_preserved',DATA+'''gen _dea_v1=991
gen _dea_u1=992
gen _dea_u0=993
dea x = y, dmu(id) model(multiplier) stage(1) saveweights gen(res) replacevars
assert _dea_v1==991 & _dea_u1==992 & _dea_u0==993
''')
    run('inactive_strategy_preserved',DATA+'''gen res_ben=991
gen res_nash=992
dea_crossgame, inputs(x) outputs(y) dmu(id) strategy(aggressive) gen(res) replacevars
assert res_ben==991 & res_nash==992
''')
    run('inactive_sw_bias_preserved',RANDOM+'''gen res_eff_bc=991
dea_sw2007 y, inputs(x) envars(z) dmu(id) rts(crs) reps(10) bcreps(10) seed(553) gen(res) nobcorrected replacevars
assert res_eff_bc==991
''')
    run('external_backend_preserves_columns',RANDOM+'''gen res_eff=991
gen res_eff_bc=992
'''+unchanged('dea_sw2007 y, inputs(x) envars(z) dmu(id) rts(crs) reps(10) bcreps(10) gen(res) backend(simarwilson) replacevars',199))
    run('window_auto_id_rejection_order',DATA+'''gsort -id
gen res_mean=991
'''+unchanged('dea_window y, inputs(x) time(t) width(2) gen(res)'))
    run('window_auto_id_success_order',DATA+'''gsort -id
quietly datasignature
local before "`r(datasignature)'"
dea_window y, inputs(x) time(t) width(2) gen(res)
keep id firm t x y z w p bad lo hi note
quietly datasignature
assert "`before'"=="`r(datasignature)'"
''')
    summary=dict(total=len(results),passed=sum(r['passed'] for r in results),failed=sum(not r['passed'] for r in results))
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
    raise SystemExit(bool(summary['failed']) or not results)

if __name__=='__main__':main()
