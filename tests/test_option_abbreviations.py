#!/usr/bin/env python3
"""Audit extension option abbreviations against help and real Stata syntax.

Uses the standard library only. Parser probes compile the exact public syntax
block from each snapshotted ado; separate complete-command runs check outputs,
returns, alias conflicts and preservation of the user's data.
"""
from _support import package_files
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time

from test_io_syntax import CAPTURE, VERIFY, OUTS, fixture, command
from test_output_protection import unchanged, RANDOM

REPO=Path(__file__).resolve().parents[1]
VALUES=dict(inputs='x1 x2',outputs='y1 y2',intermediates='z1 z2',dmu='obsid',rts='vrs',
 orientation='output',ort='output',method='independent',vrscross='legacy',generate='res',tolerance='1e-8',
 reps='50',bcreps='10',level='90',seed='551',alpha='0 .5 1',xlo='lo_x1 lo_x2',
 xhi='hi_x1 hi_x2',ylo='lo_y1 lo_y2',yhi='hi_y1 hi_y2',time='period',width='2',
 bad='bad1 bad2',leader='stage2',strategy='all',maxiter='100',envars='env',scores='score',
 panel='firm',tau='1',saving='"trial = data, file.dta"',procedure='sequential',backend='internal',engine='ipm',congmethod='cooper',congweight='raw',
 congpolicy='stage1',conganchor='A',congreference='M',congselectlimit='16')
REQUIRED=dict(network=['intermediates'],stackelberg=['intermediates'],malmquist=['time'],
 window=['time'],sbmundesirable=['bad'],sw2007=['envars'])


def inventory(ado):
    code=ado.read_text()
    # Read the continued syntax command; blank lines are formatting, not syntax.
    start=re.search(r'^    syntax ',code,re.M)
    assert start is not None, ado.name
    lines=code[start.start():].splitlines();parts=[]
    for line in lines:
        parts.append(line)
        if '///' not in line:break
    block='\n'.join(parts)
    options={}
    for m in re.finditer(r'\b([A-Za-z][A-Za-z0-9]*)(?:\(([^()]*)\))?',block.split(',',1)[1].replace('///','')):
        declaration,arg=m.groups();full=declaration.lower()
        options[full]=dict(declaration=declaration,argument=arg is not None,
                          minimum=sum(c.isupper() for c in declaration) or len(full))
    docs={}
    helptext=ado.with_suffix('.sthlp').read_text()
    for short,tail in re.findall(r'\{cmdab:([a-z]+):([^}]+)\}',helptext,re.I):
        full=short+re.match(r'[A-Za-z]*',tail).group();docs[full.lower()]=short.lower()
    for short,tail in re.findall(r'\{opth?:([a-z]+):([a-z]+)',helptext,re.I):
        docs[(short+tail).lower()]=short.lower()
    assert set(docs)<=set(options),(ado.name,set(docs)-set(options))
    return dict(syntax=block,options=options,documented_minimum=docs,
                shared_io='_dea_parse_io `0\'' in code,
                forwards_options=bool(re.search(r'\*\s*\]',block)))


def token(full,info,spelling=None):
    spelling=spelling or full
    return spelling+'('+VALUES[full]+')' if info['argument'] else spelling


def parser_program(info):
    s='program define syntax_probe, rclass\nversion 16\n'+info['syntax']+'\n'
    for full in info['options']:s+=f'return local {full} `"`{full}\'"\'\n'
    if info['forwards_options']:s+='return local passthrough `"`options\'"\'\n'
    return s+'end\n'


def parser_args(key,info,target=None,spelling=None):
    names=['inputs','outputs']+REQUIRED.get(key,[])
    names=[n for n in names if n in info['options']]
    if target and target not in names:names.append(target)
    positional='x1 x2 = y1 y2' if key=='lp' else ''
    return positional+', '+' '.join(token(n,info['options'][n],spelling if n==target else None) for n in names)


def setup(key):
    if key=='sw2007':
        return RANDOM+'gen long obsid=id\ngen byte group=id<=75\ngen double score=y/x\ngen str30 tag=cond(group,"keep = this, yes","omit")\n'
    s=fixture(key)
    if key=='fuzzy':
        s+='foreach v in x1 x2 y1 y2 {\n gen double lo_`v\'=.9*`v\'\n gen double hi_`v\'=1.1*`v\'\n}\n'
    return s


def complete(key):
    """Full names with nondefault options, before abbreviating any token."""
    if key=='lp':
        # Stage-2 forwarding must preserve matrix names as well as values,
        # independent of the temporary variables allocated by each spelling.
        return 'dea_lp y1 y2 if group==1, inputs(x1 x2) dmu(obsid) generate(res) rts(vrs) engine(simplex) method(simplex) stage(2)',['res'],'r'
    if key in OUTS:
        extra=''
        if key in ['crosseff','network','malmquist','window']:extra='rts(vrs) orientation(output) tolerance(1e-8)'
        if key=='crosseff':extra+=' vrscross(legacy)'
        if key=='sbmundesirable':extra='rts(vrs) tolerance(1e-8)'
        return command(key,'B',q='if group==1',extra=extra),OUTS[key],'r'
    base=f'dea_{key} if group==1, inputs(x1 x2) outputs(y1 y2) dmu(obsid) generate(res)'
    if key=='cooperative':return base+' rts(vrs)',['res_indiv','res_shapley','res_contrib','res_rank'],'r'
    if key=='crossgame':return base+' strategy(all) matrix',['res_agg','res_ben','res_nash'],'r'
    if key=='stackelberg':return base+' intermediates(z1 z2) rts(vrs) leader(stage2)',['res_'+s for s in ['leader','follower','overall','indep1','indep2']],'r'
    if key=='bootstrap':return base+' rts(vrs) orientation(output) reps(50) level(90) seed(551)',['res','res_bc','res_se','res_lo','res_hi'],'r'
    if key=='fuzzy':return base+' xlo(lo_x1 lo_x2) xhi(hi_x1 hi_x2) ylo(lo_y1 lo_y2) yhi(hi_y1 hi_y2) alpha(0 .5 1) rts(ndrs) orientation(output) tolerance(1e-10)',[f'res_{s}_a{a}' for a in [0,50,100] for s in ['pess','opt']],'r'
    if key=='sw2007':return 'dea_sw2007 if group==1, inputs(x) outputs(y) envars(z) dmu(id) generate(res) rts(crs) orientation(input) reps(10) bcreps(10) seed(551) procedure(sequential) backend(internal)',['res_eff','res_eff_bc'],'e'
    raise ValueError(key)


def abbreviate(cmd,mapping):
    for full,short in mapping.items():
        cmd=re.sub(r'(?<![A-Za-z0-9_])'+full+r'(?=\(|\s|$)',short,cmd)
    return cmd


def to_form_a(cmd):
    x=re.search(r' inputs\(([^)]+)\)',cmd).group(1);y=re.search(r' outputs\(([^)]+)\)',cmd).group(1)
    cmd=re.sub(r' (inputs|outputs)\([^)]+\)','',cmd)
    name,rest=cmd.split(' ',1)
    return name+' '+x+' = '+y+' '+rest


def equivalent(key,full,short,outputs,returns):
    # Replace result-class references without changing scalar() to scalae().
    capture=re.sub(r'(?<![A-Za-z0-9_])r\(',returns+'(',CAPTURE)
    verify=re.sub(r'(?<![A-Za-z0-9_])r\(',returns+'(',VERIFY)
    return setup(key)+'''tempfile original expected
quietly ds
local original_vars "`r(varlist)'"
save `original'
preserve
set seed 883
'''+full+'\n'+capture+'keep '+' '.join(outputs)+'\nsave `expected\'\nrestore\nset seed 883\n'+short+'\n'+verify+f'count if !missing({outputs[0]})\nassert r(N)>0\ncf {" ".join(outputs)} using `expected\', all\ncf `original_vars\' using `original\', all\n'


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata',required=True);ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--baseline-ado',type=Path);ap.add_argument('--only');args=ap.parse_args()
    if args.baseline_ado and not (args.baseline_ado/'dea_fuzzy.ado').is_file():ap.error('Invalid pre-F24 ado snapshot')
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'current_ado';pkg.mkdir();hashes={}
    names=package_files(REPO)
    for name in names:
        b=(REPO/name).read_bytes();(pkg/Path(name).name).write_bytes(b);hashes[name]=hashlib.sha256(b).hexdigest()
    inventories={p.stem[4:]:inventory(p) for p in sorted(pkg.glob('dea_*.ado')) if p.with_suffix('.sthlp').exists()}
    missing={key:sorted(opt for opt,decl in info['options'].items() if decl['argument'] and opt not in VALUES) for key,info in inventories.items()}
    assert not any(missing.values()),('Add test values for new argument-bearing options',missing)
    (out/'inventory.json').write_text(json.dumps(inventories,indent=2))
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,baseline=str(args.baseline_ado),stata=args.stata),indent=2));rows=[]

    def run(name,script,source=pkg,**metadata):
        if args.only and not re.search(args.only,name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 120\n'
        for f in ['PLUS','PERSONAL','SITE','OLDPLACE']:
            (d/f).mkdir();head+=f'sysdir set {f} "{d/f}"\n'
        head+=f'adopath ++ "{source.resolve()}"\n';(d/'test.do').write_text(head+script+'\ndi "TEST_PASS"\n')
        start=time.monotonic();row=dict(test=name,passed=False,**metadata)
        try:
            subprocess.run([args.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=240)
            log=(d/'test.log').read_text(errors='replace')
            if 'TEST_PASS' not in log.splitlines():raise RuntimeError(log[-2400:])
            row['passed']=True
        except Exception as e:row['error']=str(e)
        row['seconds']=round(time.monotonic()-start,3);rows.append(row)
        (out/'results.json').write_text(json.dumps(rows,indent=2));print(json.dumps(row),flush=True)

    for key,info in inventories.items():
        probe_data=fixture('crosseff')+'gen double env=firm/5\ngen double score=2+firm/10\n'
        probe_data+='foreach v in x1 x2 y1 y2 {\n gen double lo_`v\'=.9*`v\'\n gen double hi_`v\'=1.1*`v\'\n}\n'
        prelude=probe_data+parser_program(info)
        full,outputs,returns=complete(key)
        oldpath=args.baseline_ado/f'dea_{key}.ado' if args.baseline_ado else None
        old=inventory(oldpath) if oldpath and oldpath.exists() else None
        if old:
            assert info['syntax'].lower()==old['syntax'].lower(),key
            assert all(decl['minimum']<=old['options'][opt]['minimum'] for opt,decl in info['options'].items()),key
        for opt,decl in info['options'].items():
            doc_min=len(info['documented_minimum'].get(opt,opt))
            assert decl['minimum']<=doc_min,(key,opt,'help promises an unsupported abbreviation')
            s=prelude+'syntax_probe '+parser_args(key,info,opt)+'\nlocal expected `"`r('+opt+')\'"\'\n'
            spellings=[opt[:n] for n in range(decl['minimum'],len(opt)+1)]
            for spelling in spellings:
                s+='syntax_probe '+parser_args(key,info,opt,spelling)+f'\nassert `"`r({opt})\'"\'==`"`expected\'"\'\n'
            run(f'parser_{key}_{opt}',s,prefixes_checked=len(spellings))
            # Repeated argument-bearing options, including required inputs,
            # must not silently choose one of the spellings or values.
            if decl['argument']:
                duplicate='syntax_probe '+parser_args(key,info,opt,spellings[0])+' '+token(opt,decl)
                if info['forwards_options']:
                    # syntax ..., * retains a duplicate for the delegated
                    # command; public-path duplicate rejection is tested below.
                    script=prelude+duplicate+f'\nassert `"`r({opt})\'"\'!=""\n'
                    script+=f'assert strpos(`"`r(passthrough)\'"\',"{opt}(")>0\n'
                    run(f'parser_duplicate_{key}_{opt}',script)
                else:
                    run(f'parser_duplicate_{key}_{opt}',prelude+unchanged(duplicate,198))
            if decl['minimum']>1:
                shorter=opt[:decl['minimum']-1]
                resolves_elsewhere=any(other.startswith(shorter) and d['minimum']<=len(shorter) for other,d in info['options'].items())
                if not resolves_elsewhere:
                    if info['forwards_options']:
                        # A forwarding wrapper retains unknown options in *;
                        # they must not populate a named wrapper option.
                        script=prelude+'syntax_probe '+parser_args(key,info,opt,shorter)+'\n'
                        script+=f'assert `"`r({opt})\'"\'==""\nassert strpos(`"`r(passthrough)\'"\',"{shorter}(")>0\n'
                        run(f'parser_too_short_{key}_{opt}',script)
                    else:
                        run(f'parser_too_short_{key}_{opt}',prelude+unchanged('syntax_probe '+parser_args(key,info,opt,shorter),198))
            # Prove each advertised, formerly rejected minimum on the old parser.
            if old and opt in info['documented_minimum'] and old['options'][opt]['minimum']>doc_min:
                old_prelude=probe_data+parser_program(old)
                run(f'old_rejected_{key}_{opt}',old_prelude+unchanged('syntax_probe '+parser_args(key,old,opt,opt[:doc_min]),198),args.baseline_ado)
        mapping={k:v for k,v in info['documented_minimum'].items()}
        if info['forwards_options']:
            mapping.update({k:k[:v['minimum']] for k,v in info['options'].items()})
        forms=['B','A'] if info['shared_io'] else ['B']
        for form in forms:
            reference=full if form=='B' else to_form_a(full)
            short=abbreviate(reference,mapping)
            run(f'complete_{key}_{form}',equivalent(key,reference,short,outputs,returns))
            complex_condition='if inrange(group,1,1) & tag=="keep = this, yes" in 1/l'
            complex_full=reference.replace('if group==1',complex_condition)
            run(f'complex_condition_{key}_{form}',equivalent(key,complex_full,abbreviate(complex_full,mapping),outputs,returns))
            # Reject duplicates on the public path before any existing result
            # columns or model data can change, even with replacevars.
            sentinel=''.join(f'gen str12 {v}="old result"\n' for v in outputs)
            for opt in ['generate','rts','orientation','intermediates','engine','method']:
                if opt not in info['options'] or f'{opt}(' not in reference:continue
                current_value=re.search(opt+r'\(([^)]+)\)',reference).group(1)
                alias=mapping.get(opt,opt)
                run(f'complete_duplicate_{key}_{form}_{opt}',setup(key)+sentinel+unchanged(short+' '+opt+'('+current_value+') replacevars',198))
            if 'orientation' in info['options'] and 'orientation(' in reference:
                for value in ['input','output']:
                    run(f'alias_conflict_{key}_{form}_{value}',setup(key)+sentinel+unchanged(short+f' ort({value}) replacevars',198))
                    first=short.replace(',',f', ort({value})',1)
                    run(f'alias_conflict_first_{key}_{form}_{value}',setup(key)+sentinel+unchanged(first+' replacevars',198))
                legacy=abbreviate(reference,{'orientation':'ort'})
                run(f'ort_alias_{key}_{form}',equivalent(key,reference,legacy,outputs,returns))
            # Input/intermediate/output options must resolve independently.
            if key=='network':
                for inp in (['i','in','inp','inputs'] if form=='B' else ['inputs']):
                    for inter in ['int','inte','inter','intermediates']:
                        short=abbreviate(reference,{'inputs':inp,'intermediates':inter,'orientation':'or','rts':'r'})
                        run(f'input_intermediate_{form}_{inp}_{inter}',equivalent(key,reference,short,outputs,returns))
        if old:
            # Full-name historical results remain exactly the same.
            hist='historical_'+key
            d=out/(hist+'_control')
            if not args.only or re.search(args.only,hist):
                run(hist+'_control',setup(key)+'set seed 883\n'+full+'\nkeep '+' '.join(outputs)+'\nsave "expected.dta"\n',args.baseline_ado)
                run(hist,setup(key)+'set seed 883\n'+full+f'\ncf {" ".join(outputs)} using "{d/"expected.dta"}", all\n')
            if key=='fuzzy':
                run('old_fuzzy_ignored_ort',equivalent(key,full,full+' ort(input)',outputs,returns),args.baseline_ado)
    summary=dict(total=len(rows),passed=sum(r['passed'] for r in rows),parser_prefixes=sum(r.get('prefixes_checked',0) for r in rows),failed=[r['test'] for r in rows if not r['passed']])
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
    raise SystemExit(bool(summary['failed']))


if __name__=='__main__':main()
