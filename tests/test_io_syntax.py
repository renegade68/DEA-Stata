#!/usr/bin/env python3
"""Exercise shared DEA IO parsing against physically selected sample controls.

Python standard library and Stata 16+ only. Historical controls require an
explicit pre-F06 ado directory. No installed ado paths are changed.
"""
from _support import package_files
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import re
import subprocess
import time

from test_output_protection import unchanged

REPO=Path(__file__).resolve().parents[1]
KEYS=['crossgame','crosseff','nash','network','malmquist','window','congestion','sbmundesirable']
POSITIONAL={'network','malmquist','window'}
PANEL={'malmquist','window'}
DATA='''input firm x1 x2 y1 y2
1 .1 .2 12 10
2 2 3 2 1
3 3 2 3 2
4 4 5 3 4
5 2 4 1 2
6 5 3 4 3
end
'''
OUTS={
 'crossgame':['res_agg','res_ben','res_nash'],
 'crosseff':['res','res_self','res_maverick'],
 'nash':['res_'+s for s in ['disagree','coop','bargain','gain']],
 'network':['res_'+s for s in ['overall','stage1','stage2']],
 'malmquist':['res_'+s for s in ['ec','tc','mpi']],
 'window':['res_mean','res_sd','res_w1_p1','res_w1_p2','res_w2_p1','res_w2_p2'],
 'congestion':['res','res_crs','res_vrs','res_cong_eff','res_congestion','res_se'],
 'sbmundesirable':['res'],
}
QUALIFIERS={
 'boolean':'if group', 'equality':'if group==1',
 'inequality':'if firm!=1 & firm<=5', 'tilde':'if firm~=1 & firm<=5',
 'string':'if tag=="keep = this, yes"',
 'function':'if inlist(firm,2,3,4,5)',
 'nested':'if inrange(firm,2,5) & (group==1 | firm==3)',
 'paren':'if(group==1)', 'in':'in 2/l',
 'if_in':'if group==1 in 2/l', 'negative_in':'in -4/l',
 'all':'if firm>=1 in 1/l',
}


def fixture(key,clean=False):
    s=DATA
    if key in PANEL:
        s+='expand 3\nbysort firm: gen period=_n\nreplace y1=y1*(1+.1*(period-1))\nreplace y2=y2*(1+.1*(period-1))\nsort period firm\n'
    else:s+='gen period=1\n'
    s+='''gen long obsid=_n
gen byte group=inrange(firm,2,5)
gen double z1=(x1+y1)/2
gen double z2=(x2+y2)/2
gen double bad1=1+.1*firm
gen double bad2=2+.2*firm
'''
    if not clean:s+='gen str30 tag=cond(group,"keep = this, yes","omit")\ngen str12 sid="Firm "+string(firm)\ngen double unused=.\n'
    return s


def command(key,form='B',q='',extra='',dmu=None,inputs='x1 x2',outputs='y1 y2',prefix='res'):
    if dmu is None:dmu='firm' if key in PANEL else 'obsid'
    if form=='A':io=f'{inputs} = {outputs} {q},'
    elif form=='P':io=f'{outputs} {q}, inputs({inputs})'
    else:io=f'{q}, inputs({inputs}) outputs({outputs})'
    req={'network':'intermediates(z1 z2) matrix','malmquist':'time(period)','window':'time(period) width(2)','sbmundesirable':'bad(bad1 bad2)','crossgame':'matrix','crosseff':'matrix'}.get(key,'')
    return f'dea_{key} {io} {req} {"dmu("+dmu+")" if dmu else ""} generate({prefix}) {extra}'


CAPTURE='''local scalar_names : r(scalars)
local macro_names : r(macros)
local matrix_names : r(matrices)
local sidx=0
foreach s of local scalar_names {
 local ++sidx
 tempname exs`sidx'
 scalar `exs`sidx''=r(`s')
}
local midx=0
foreach s of local macro_names {
 local ++midx
 local exm`midx' `"`r(`s')'"'
}
local aidx=0
foreach s of local matrix_names {
 local ++aidx
 tempname exa`aidx'
 matrix `exa`aidx''=r(`s')
 local rn`aidx' : rowfullnames r(`s')
 local cn`aidx' : colfullnames r(`s')
}
'''
VERIFY='''assert "`:r(scalars)'"=="`scalar_names'"
assert "`:r(macros)'"=="`macro_names'"
assert "`:r(matrices)'"=="`matrix_names'"
local sidx=0
foreach s of local scalar_names {
 local ++sidx
 assert r(`s')==scalar(`exs`sidx'')
}
local midx=0
foreach s of local macro_names {
 local ++midx
 // The delegated core command records a call containing its temporary
 // sample-marker name. That call is intentionally different in each run.
 if "`s'"!="cmdline" assert `"`r(`s')'"'==`"`exm`midx''"'
}
local aidx=0
foreach s of local matrix_names {
 local ++aidx
 tempname got
 matrix `got'=r(`s')
 mata: assert(st_matrix("`got'")==st_matrix("`exa`aidx''"))
 assert "`:rowfullnames r(`s')'"=="`rn`aidx''"
 assert "`:colfullnames r(`s')'"=="`cn`aidx''"
}
'''


def pair(key,form='B',q='if group==1',extra='',setup='',dmu=None,inputs='x1 x2',outputs='y1 y2'):
    # Select using Stata itself so commas, string literals and relative in
    # ranges have exactly the public language's semantics.
    outs=OUTS[key]
    if key=='crossgame':
        strategy=re.search(r'strategy\(([^)]+)\)',extra)
        if strategy and strategy[1]!='all':outs=['res_'+dict(aggressive='agg',benevolent='ben',nash='nash')[strategy[1]]]
    if key=='window' and q=='in -12/l':outs=outs[:4]
    s=fixture(key)+setup+f'gen byte selected=0\nreplace selected=1 {q}\n'
    required='x1 x2 y1 y2'+(' z1 z2' if key=='network' else '')+(' bad1 bad2' if key=='sbmundesirable' else '')+(' period' if key in PANEL else '')
    for v in required.split():s+=f'replace selected=0 if missing({v})\n'
    if dmu:s+=f'replace selected=0 if missing({dmu})\n'
    s+='quietly ds\nlocal original_vars "`r(varlist)\'"\ntempfile original expected\nsave `original\'\npreserve\nkeep if selected\n'
    # The control always uses the canonical, fully specified syntax.
    s+=command(key,'B',extra=extra,dmu=dmu)+'\n'+CAPTURE
    s+='keep obsid '+' '.join(outs)+'\nsave `expected\'\nrestore\n'
    s+=command(key,form,q,extra,dmu,inputs,outputs)+'\n'+VERIFY
    for v in outs:s+=f'assert missing({v}) if !selected\n'
    # Do not let equality of two entirely missing runs count as success.
    s+=f'count if selected & !missing({outs[0]})\nassert r(N)>0\n'
    s+='preserve\nkeep if selected\ncf obsid '+' '.join(outs)+' using `expected\', all\nrestore\n'
    s+='local replaced "'+' '.join(outs)+'"\nlocal original_vars : list original_vars - replaced\ncf `original_vars\' using `original\', all\n'
    return s


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata',required=True);ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--baseline-ado',type=Path);ap.add_argument('--only');args=ap.parse_args()
    if args.baseline_ado and not (args.baseline_ado/'dea_network.ado').is_file():ap.error('--baseline-ado must contain the pre-F06 ado files')
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'current_ado';pkg.mkdir()
    names=package_files(REPO);hashes={}
    for name in names:
        data=(REPO/name).read_bytes();(pkg/Path(name).name).write_bytes(data);hashes[name]=hashlib.sha256(data).hexdigest()
    manifest=dict(source_sha256=hashes,stata=args.stata,baseline=str(args.baseline_ado))
    if args.baseline_ado:manifest['baseline_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in args.baseline_ado.iterdir() if p.is_file()}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2));results=[]

    def run(name,script,source=pkg,force=False):
        if not force and args.only and not re.search(args.only,name):return None
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 120\n'
        for folder in ['PLUS','PERSONAL','SITE','OLDPLACE']:
            (d/folder).mkdir();head+=f'sysdir set {folder} "{d/folder}"\n'
        head+=f'adopath ++ "{source.resolve()}"\n';(d/'test.do').write_text(head+script+'\ndi "TEST_PASS"\n')
        start=time.monotonic();row=dict(test=name,passed=False)
        try:
            subprocess.run([args.stata,'-b','do','test.do'],cwd=d,timeout=180,capture_output=True)
            log=(d/'test.log').read_text(errors='replace')
            if 'TEST_PASS' not in log.splitlines():raise RuntimeError(log[-2000:])
            row['passed']=True
        except Exception as e:row['error']=str(e)
        row['seconds']=round(time.monotonic()-start,3);results.append(row)
        (out/'results.json').write_text(json.dumps(results,indent=2));print(json.dumps(row),flush=True)
        return d if row['passed'] else None

    for key in KEYS:
        forms=['A','B']+(['P'] if key in POSITIONAL else [])
        for form in forms:
            for label,q in QUALIFIERS.items():
                # Last four rows would retain only one period in the panel.
                if key in PANEL and label=='negative_in':q='in -12/l'
                run(f'qual_{key}_{form}_{label}',pair(key,form,q))
            run(f'string_id_{key}_{form}',pair(key,form,dmu='sid'))
            for label,x,y in [('range','x1-x2','y1-y2'),('wildcard','x?','y*')]:
                run(f'varlist_{key}_{form}_{label}',pair(key,form,inputs=x,outputs=y))
            for var in ['x1','y2']:
                run(f'missing_{key}_{form}_{var}',pair(key,form,setup=f'replace {var}=. in 4\n'))
            first=OUTS[key][0]
            run(f'replace_{key}_{form}',pair(key,form,q=f'if {first}==1',setup=f'gen {first}=group\n',extra='replacevars'))
            errors=[('duplicate','inputs(x1)',198),('unknown','nonsense',198)]
            # Duplicate inputs is already present in all normalized forms.
            for label,opt,rc in errors:
                run(f'invalid_{key}_{form}_{label}',fixture(key)+unchanged(command(key,form,q='if group==1',extra=opt),rc))
            run(f'protect_{key}_{form}',fixture(key)+f'rename x1 {first}\n'+unchanged(command(key,form,q='if group==1',inputs=f'{first} x2',extra='replacevars')))
        for form in ['A','B']:
            run(f'missing_outputs_{key}_{form}',fixture(key)+unchanged(command(key,form,outputs=''),100 if form=='A' else 198))
            run(f'string_output_{key}_{form}',fixture(key)+unchanged(command(key,form,outputs='tag'),109))

        variants={
            'crossgame':[f'strategy({st}) maxiter({lim}) detail debug tolerance(1e-8)' for st,lim in itertools.product(['all','aggressive','benevolent','nash'],[1,30])],
            'crosseff':[f'rts({r}) ort({o}) method({m}) detail debug tolerance(1e-8)' for r,o,m in itertools.product(['crs','vrs'],['input','output'],['arbitrary','aggressive','benevolent'])],
            'nash':[f'rts({r}) detail debug tolerance(1e-8)' for r in ['crs','vrs']],
            'network':[f'rts({r}) orientation({o}) method({m}) detail debug tolerance(1e-8)' for r,o,m in itertools.product(['crs','vrs'],['input','output'],['multiplicative','additive','independent'])],
            'malmquist':[f'rts({r}) ort({o}) {g} detail debug tolerance(1e-8)' for r,o,g in itertools.product(['crs','vrs'],['input','output'],['','global'])],
            'window':[f'rts({r}) orientation({o}) detail debug tolerance(1e-8)' for r,o in itertools.product(['crs','vrs','drs','irs'],['input','output'])],
            'congestion':[f'engine({e}) tolerance(1e-8)' for e in ['simplex','ipm']],
            'sbmundesirable':[f'rts({r}) tolerance(1e-8)' for r in ['crs','vrs']],
        }[key]
        for i,extra in enumerate(variants):
            for form in ['A','B']:
                run(f'options_{key}_{form}_{i}',pair(key,form,q='if inlist(firm,2,3,4,5)',extra=extra))

        if args.baseline_ado:
            # Old = anywhere bug: keep data numeric so failure has a single cause.
            run(f'old_equality_{key}',fixture(key,clean=True)+unchanged(command(key,'B',q='if group==1'),198),args.baseline_ado)
            if key in ['crosseff','nash','congestion','sbmundesirable']:
                run(f'old_form_a_if_{key}',fixture(key,clean=True)+unchanged(command(key,'A',q='if group'),198),args.baseline_ado)
            name='historical_'+key
            wanted=not args.only or re.search(args.only,name)
            if wanted:
                c=run('before_'+key,fixture(key)+command(key,'A')+'\nkeep obsid '+' '.join(OUTS[key])+'\nsave "expected.dta"\n',args.baseline_ado,force=True)
                if c:run(name,fixture(key)+command(key,'A')+f'\ncf obsid {" ".join(OUTS[key])} using "{c/"expected.dta"}", all\n')
    summary=dict(total=len(results),passed=sum(r['passed'] for r in results),failed=[r['test'] for r in results if not r['passed']])
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
    raise SystemExit(bool(summary['failed']))


if __name__=='__main__':main()
