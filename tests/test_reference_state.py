#!/usr/bin/env python3
"""Reference display preserves caller objects and renders the returned lambda.

Stdlib + Stata 16+. Includes data-name shadowing, fabricated display weights,
private forced retry, historical controls and success/error cleanup checks.
"""
from _support import package_files
import argparse,csv,hashlib,itertools,json,math,re,subprocess,time
from pathlib import Path
from test_saveref import instrument

REPO=Path(__file__).resolve().parents[1]
NAMES=['_dmu_k','_dmu_j','_lam_val']
DATA='''input double(obsno x1 x2 y1 y2)
1 1 2 1 2
2 2 1 3 3
3 3 4 2 5
4 4 3 5 4
5 6 5 7 6
6 7 4 5 8
7 9 6 8 7
8 10 7 11 9
end
gen str20 note="original data"
char note[audit] "unchanged"
'''
QUAL=[('','1'),('if obsno!=2','obsno!=2'),('in 2/8','inrange(_n,2,8)'),('if obsno!=2 in 1/7','obsno!=2 & inrange(_n,1,7)')]
SCALARS=['N','N_eff','mean_eff','min_eff','max_eff','n_failed']


def guard(state):
    setup='scalar user_number=987.25\nscalar user_string="keep user text 한글"\nmatrix user_matrix=(1,2\\3,4)\n'
    expected={'user_number':('number','987.25'),'user_string':('string','keep user text 한글'),'user_matrix':('matrix',1)}
    for i,n in enumerate(NAMES):
        kind={'absent':'none','numeric':'number','string':'string','matrix':'matrix','mixed':['number','string','matrix'][i]}[state]
        if kind=='none':continue
        val='.a' if state=='mixed' and i==0 else str(991+i) if kind=='number' else f'keep {n} 한글' if kind=='string' else 991+i
        expected[n]=(kind,val)
        if kind=='number':setup+=f'scalar {n}={val}\n'
        if kind=='string':setup+=f'scalar {n}="{val}"\n'
        if kind=='matrix':setup+=f'matrix {n}=({val},2\\3,4)\n'
    for n,(kind,val) in expected.items():
        if kind=='matrix':setup+=f'matrix rownames {n}=keptR1 keptR2\nmatrix colnames {n}=keptC1 keptC2\n'
    setup+='local before_s : all scalars\nlocal before_s : list sort before_s\nlocal before_m : all matrices\nlocal before_m : list sort before_m\n'
    check='local after_s : all scalars\nlocal after_s : list sort after_s\nlocal after_m : all matrices\nlocal after_m : list sort after_m\nassert "`before_s\'"=="`after_s\'"\nassert "`before_m\'"=="`after_m\'"\n'
    for n,(kind,val) in expected.items():
        if kind=='number':check+=f'assert scalar({n})=={val}\n'
        if kind=='string':check+=f'assert scalar({n})=="{val}"\n'
        if kind=='matrix':check+=f'mata: assert(all(st_matrix("{n}"):==({val},2\\3,4)))\nassert "`:rownames {n}\'"=="keptR1 keptR2"\nassert "`:colnames {n}\'"=="keptC1 keptC2"\n'
    return setup,check


def fixture(profile,i=0):
    s=DATA
    if i%5==0:s+='gsort -obsno\n'
    if profile=='auto':s+='gen str20 display_label="dmu"+string(_n)\n';dmu=''
    elif profile in ['string','unicode','duplicate']:
        expr={'string':'"unit "+string(obsno)','unicode':'"연구소 "+string(obsno)','duplicate':'"group "+string(mod(obsno,3))'}[profile]
        s+=f'gen str30 id={expr}\ngen str30 display_label=id\n';dmu='dmu(id)'
    else:
        expr={'numeric':'obsno','fraction':'obsno+.25','negative':'-obsno','missing':'cond(obsno==4,.a,obsno)'}.get(profile,'obsno')
        v={'scratch_k':'_dmu_k','scratch_j':'_dmu_j','scratch_lam':'_lam_val'}.get(profile,'id')
        s+=f'gen double {v}={expr}\ngen str30 display_label=strtrim(string({v},"%18.0g"))\n';dmu=f'dmu({v})'
    return s,dmu


def command(form='A',dmu='dmu(id)',engine='auto',rts='crs',ori='input',stage=2,q='',model='te',prefix='res',extra='reference'):
    head=f'dea x1 x2=y1 y2 {q},' if form=='A' else f'dea y1 y2 {q}, inputs(x1 x2)'
    return f'{head} {dmu} model({model}) engine({engine}) rts({rts}) orientation({ori}) stage({stage}) gen({prefix}) {extra}'


def export_reference():
    return '''preserve
keep if selected
keep display_label
export delimited using "labels.csv",replace
restore
preserve
clear
svmat double AuditLambda, names(c)
format c* %24.17g
export delimited using "lambda.csv",replace datafmt
restore
'''


def reference_lines(log):
    blocks=[];lines=log.splitlines()
    for i,line in enumerate(lines):
        if line.startswith('Reference Sets (lambda > '):
            block=[]
            for row in lines[i+2:]:
                if re.fullmatch('-+',row):break
                if row.strip():block.append(re.sub(r'\s+',' ',row).strip())
            blocks.append(block)
    return blocks


def verify_display(d,expected=None):
    if expected is None:
        with (d/'labels.csv').open() as f:labels=[r['display_label'] for r in csv.DictReader(f)]
        with (d/'lambda.csv').open() as f:rows=list(csv.DictReader(f))
        assert len(rows)==len(labels),(len(rows),len(labels))
        expected=[]
        for label,row in zip(labels,rows):
            peers=[]
            for peer,raw in zip(labels,row.values()):
                value=float(raw) if raw and not raw.startswith('.') else (float(raw) if re.match(r'^\.\d',raw) else math.nan)
                if math.isfinite(value) and value>1e-7:peers.append(f'{peer}({value:6.4f})')
            if peers:
                rendered=label+': '+' '.join(peers[:6])+ (f' ... (+{len(peers)-6} more)' if len(peers)>6 else '')
                expected.append(re.sub(r'\s+',' ',rendered).strip())
    blocks=reference_lines((d/'test.log').read_text(errors='replace'));assert blocks,'No reference block'
    assert blocks[-1]==expected,(blocks[-1],expected)
    return dict(reference_rows=len(expected))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stata',required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--baseline-ado',type=Path);p.add_argument('--only');a=p.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'current_ado';pkg.mkdir();hashes={}
    for n in package_files(REPO):
        raw=(REPO/n).read_bytes();(pkg/Path(n).name).write_bytes(raw);hashes[n]=hashlib.sha256(raw).hexdigest()
    retry=out/'retry_ado';instrument(pkg,retry)
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,baseline=str(a.baseline_ado),stata=a.stata),indent=2));results=[]
    def run(name,script,source=pkg,check=None):
        if a.only and not re.search(a.only,name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 255\n'
        for f in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/f).mkdir();head+=f'sysdir set {f} "{d/f}"\n'
        head+=f'adopath ++ "{source.resolve()}"\n';(d/'test.do').write_text(head+script+'\ndi "TEST_PASS"\n');r=dict(test=name,passed=False);start=time.monotonic()
        try:
            proc=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=180)
            if proc.returncode or not (d/'test.log').exists():(d/'process.json').write_text(json.dumps(dict(rc=proc.returncode,stderr=proc.stderr.decode(errors='replace')),indent=2))
            log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else proc.stderr.decode(errors='replace')
            assert proc.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-1800:]
            if check:r.update(check(d) or {})
            r['passed']=True
        except Exception as e:r['error']=str(e)
        r['seconds']=round(time.monotonic()-start,3);results.append(r);(out/'results.json').write_text(json.dumps(results,indent=2));print(json.dumps(r),flush=True)
    def paired(profile,form,engine,rts,ori,stage,state,i,model='te',source=pkg):
        s,dmu=fixture(profile,i);q,cond=QUAL[i%4];s+=f'gen byte selected={cond}\n'
        opts=dict(form=form,dmu=dmu,engine=engine,rts=rts,ori=ori,stage=stage,q=q,model=model)
        control_opts=dict(opts)
        if source==retry:control_opts['engine']='ipm'
        s+=command(**control_opts,prefix='control')+'\nmatrix C=r(efficiency)\nmatrix L=r(lambda)\nmatrix Stats=('+','.join('r('+x+')' for x in SCALARS)+')\n'
        if profile=='auto':s+='drop _dmu\n'
        # Hostile data names used to shadow the public scalar reads.
        for n in NAMES:s+=f'capture confirm variable {n},exact\nif _rc gen double {n}=0\n'
        s+='quietly ds\nlocal original_vars "`r(varlist)\'"\ntempfile original\nsave `original\'\n'
        setup,check=guard(state);s+=setup
        if i%2:s+='gen str12 res="replace me"\n'
        s+=command(**opts,extra='reference saveref replacevars saving("saved.dta")')+'\n'+check
        s+='mata: assert(all(st_matrix("C"):==st_matrix("r(efficiency)")))\nmata: assert(all(st_matrix("L"):==st_matrix("r(lambda)")))\n'
        for j,k in enumerate(SCALARS,1):s+=f'assert r({k})==Stats[1,{j}]\n'
        s+='assert res==control\nlocal refs "`r(refvars)\'"\n'
        s+='mata: assert(all(st_matrix("r(lambda)"):==st_data(selectindex(st_data(.,"selected")),tokens(st_local("refs")))))\n'
        s+='matrix AuditLambda=r(lambda)\ncf `original_vars\' using `original\',all\nassert "`:char note[audit]\'"=="unchanged"\n'
        # Commands after r() capture do not destroy r(lambda); preserve/export are checked by actual output.
        s+=export_reference()
        s+='preserve\nuse "saved.dta",clear\nassert _N==rowsof(C)\nmata: assert(all(st_matrix("L"):==st_data(.,tokens(st_local("refs")))))\nrestore\ncapture ds __*\nassert _rc==111\n'
        return s
    profiles=['numeric','string','unicode','duplicate','fraction','negative','missing','auto','scratch_k','scratch_j','scratch_lam']
    for i,(engine,rts,ori,stage,form,state) in enumerate(itertools.product(['simplex','ipm','auto'],['crs','vrs','nirs','ndrs'],['input','output'],[1,2],['A','B'],['absent','numeric','string','matrix','mixed'])):
        run(f'grid_{engine}_{rts}_{ori}_{stage}_{form}_{state}',paired(profiles[(i//5)%len(profiles)],form,engine,rts,ori,stage,state,i),check=verify_display)
    for i,(engine,ori,stage,state) in enumerate(itertools.product(['simplex','ipm','auto'],['input','output'],[1,2],['numeric','matrix'])):
        run(f'super_{engine}_{ori}_{stage}_{state}',paired('numeric','A',engine,'vrs',ori,stage,state,i,model='supereff'),check=verify_display)
    for i,(rts,ori,state) in enumerate(itertools.product(['crs','vrs','nirs','ndrs'],['input','output'],['string','matrix'])):
        run(f'retry_{rts}_{ori}_{state}',paired('string','B','auto',rts,ori,2,state,i,source=retry),retry,verify_display)
    for state,kind in itertools.product(['absent','numeric','string','matrix','mixed'],['collision','save_error','quoted_label','no_reference']):
        s,dmu=fixture('numeric');setup,check=guard(state);s+=setup
        extra='reference';prefix='res';expected=0
        if kind=='collision':s+='gen res=991\n';expected=110
        if kind=='save_error':extra+=' saving("missing_folder/results.dta")';expected=603
        if kind=='quoted_label':s+='tostring id,replace\nreplace id=char(34)+id in 1\n'
        if kind=='no_reference':extra=''
        s+='capture noisily '+command(dmu=dmu,extra=extra)+f'\nassert _rc=={expected}\n'+check
        run(f'cleanup_{state}_{kind}',s)
    for target,state in itertools.product(NAMES,['absent','numeric','matrix']):
        s,dmu=fixture('numeric',1);s+='gen byte selected=1\n';setup,check=guard(state)
        s+=setup+command(dmu=dmu,prefix=target,extra='reference saveref')+'\n'+check
        s+=f'mata: assert(all(st_matrix("r(efficiency)"):==st_data(.,"{target}")))\n'
        s+='matrix AuditLambda=r(lambda)\n'+export_reference()
        run(f'output_name_{target}_{state}',s,check=verify_display)
    # A known reference matrix isolates display rules from LP nonuniqueness.
    for state in ['absent','numeric','string','matrix','mixed']:
        s='set obs 8\ngen double id=_n\ngen double score=1\ngen byte selected=1\n'
        for j in range(1,9):s+=f'gen double w{j}=0\nreplace w{j}=.125 in 1\n'
        s+='replace w1=1e-7 in 2\nreplace w2=2e-7 in 2\nreplace w3=.a in 2\nreplace w4=-1 in 2\n'
        for n in NAMES:s+=f'gen double {n}=0\n'
        setup,check=guard(state);s+=setup+'quietly run "'+str(pkg/'dea.ado')+'"\n'
        s+='_dea_display, model(te) rts(crs) orientation(input) gen(score) touse(selected) tolerance(1e-7) stage(2) dmu(id) reference lambda_vars("w1 w2 w3 w4 w5 w6 w7 w8")\n'+check
        expected=['1: '+' '.join(f'{j}(0.1250)' for j in range(1,7))+' ... (+2 more)','2: 2(0.0000)']
        run('fabricated_'+state,s,check=lambda d,e=expected:verify_display(d,e))
    if a.baseline_ado:
        for state in ['numeric','string','matrix']:
            s,dmu=fixture('numeric');setup,_=guard(state);s+=setup+command(dmu=dmu)+'\n'
            s+='assert scalar(_dmu_k)==id[_N]\nassert scalar(_dmu_j)<.\nassert scalar(_lam_val)<.\n'
            run('old_overwrite_'+state,s,a.baseline_ado)
        s,dmu=fixture('numeric');s+='gen _dmu_k=0\ngen _dmu_j=0\ngen _lam_val=0\n'+command(dmu=dmu)+'\n'
        run('old_data_shadow',s,a.baseline_ado,lambda d:verify_display(d,[]))
        for i,(engine,rts,ori,stage) in enumerate(itertools.product(['simplex','ipm','auto'],['crs','vrs','nirs','ndrs'],['input','output'],[1,2])):
            name=f'healthy_{engine}_{rts}_{ori}_{stage}';s,dmu=fixture(profiles[i%7]);cmd=command(dmu=dmu,engine=engine,rts=rts,ori=ori,stage=stage,extra='reference saveref')
            old=s+cmd+'\nlocal refs "`r(refvars)\'"\nkeep res `refs\'\nsave "expected.dta"\n'
            run(name+'_before',old,a.baseline_ado)
            now=s+cmd+f'\ncf res _ref_* using "{out/(name+"_before")/"expected.dta"}",all\n'
            def compare(d,name=name):
                old=reference_lines((out/(name+'_before')/'test.log').read_text(errors='replace'))
                new=reference_lines((d/'test.log').read_text(errors='replace'));assert old==new,(old,new)
            run(name,now,check=compare)
    summary=dict(total=len(results),passed=sum(r['passed'] for r in results),reference_rows=sum(r.get('reference_rows',0) for r in results))
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True);raise SystemExit(summary['total']!=summary['passed'])

if __name__=='__main__':main()
