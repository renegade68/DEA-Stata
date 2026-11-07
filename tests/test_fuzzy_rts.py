#!/usr/bin/env python3
"""Fuzzy RTS aliases: independent interval LPs, crisp limits and historical controls.

Requires NumPy/SciPy and Stata 16+. Tests isolate package sources and ado paths.
All interval LPs are independently constructed from the fixture and solved with
SciPy HiGHS; no production matrices or normalized RTS values serve as the oracle.
"""
from _support import package_files
import argparse,csv,hashlib,json,re,subprocess,time
from pathlib import Path
import numpy as np
import scipy
from verify_bootstrap_lp import solve
from test_output_protection import unchanged

REPO=Path(__file__).resolve().parents[1]
ALPHAS=[0,.25,.5,.75,1]
SPECS={'crs':'crs','ccr':'crs','vrs':'vrs','bcc':'vrs','drs':'drs','nirs':'drs','irs':'irs','ndrs':'irs'}
DATA='''input id x1 y1
1 .1 80
2 1 1
3 2 3
4 4 4
5 5 3
6 3 2
7 2 2
8 6 4
end
gen double x2=1+mod(id*7,5)
gen double y2=1+mod(id*3,4)
gen str8 sid="unit"+string(id)
gen str20 unrelated="preserve me"
foreach v in x1 x2 y1 y2 {
 gen double lo_`v'=`v'*(.70+.02*mod(id,4))
 gen double hi_`v'=`v'*(1.12+.03*mod(id,3))
}
'''
SAMPLES={'full':'','if':'if id>=2 & id<=7','ifin':'if sid!="unit1" in 1/6'}


def body(dim,orientation,mode,sample,old=False):
    x=' '.join(f'x{i}' for i in range(1,dim+1));y=' '.join(f'y{i}' for i in range(1,dim+1))
    lo_x=' '.join('lo_'+v for v in x.split());hi_x=' '.join('hi_'+v for v in x.split())
    lo_y=' '.join('lo_'+v for v in y.split());hi_y=' '.join('hi_'+v for v in y.split())
    bounds=''
    if mode in ['full','inputs','lower']:bounds+=f' xlo({lo_x})'
    if mode in ['full','inputs','upper']:bounds+=f' xhi({hi_x})'
    if mode in ['full','outputs','lower']:bounds+=f' ylo({lo_y})'
    if mode in ['full','outputs','upper']:bounds+=f' yhi({hi_y})'
    q=SAMPLES[sample];dmu='id' if orientation=='input' else 'sid'
    s=DATA+f'gen byte selected=0\nreplace selected=1 {q}\ncount if selected\nlocal n=r(N)\n'
    for spec,canonical in SPECS.items():
        normalized={'nirs':'irs','ndrs':'drs'}.get(spec,canonical) if old else canonical
        alpha='' if sample=='full' else 'alpha(.75 0 1 .25 .5)'
        s+=f'''dea_fuzzy {q}, inputs({x}) outputs({y}) {bounds} rts({spec}) orientation({orientation}) dmu({dmu}) gen(f_{spec}) {alpha}
assert "`r(rts)'"=="{normalized}" & "`r(orientation)'"=="{orientation}"
assert r(N)==`n' & r(m_in)=={dim} & r(m_out)=={dim}
foreach a in 0 25 50 75 100 {{
 assert missing(f_{spec}_pess_a`a') & missing(f_{spec}_opt_a`a') if !selected
 assert !missing(f_{spec}_pess_a`a',f_{spec}_opt_a`a') if selected
 assert f_{spec}_pess_a`a'<=f_{spec}_opt_a`a'+1e-6 if selected
}}
'''
        for prev,nxt in zip([0,25,50,75],[25,50,75,100]):
            s+=f'assert f_{spec}_pess_a{prev}<=f_{spec}_pess_a{nxt}+1e-6 & f_{spec}_opt_a{prev}+1e-6>=f_{spec}_opt_a{nxt} if selected\n'
    for canonical in ['crs','vrs','drs','irs']:
        for engine in ['simplex','ipm']:
            s+=f'''dea {x}={y} {q}, dmu({dmu}) rts({canonical}) ort({orientation}) engine({engine}) gen(c_{canonical}_{engine}) replacevars
assert abs(c_{canonical}_{engine}-f_{canonical}_pess_a100)<1e-6 & abs(c_{canonical}_{engine}-f_{canonical}_opt_a100)<1e-6 if selected
'''
    for alias,canonical in [('NIRS','drs'),('NdRs','irs')]:
        norm={'drs':'irs','irs':'drs'}[canonical] if old else canonical
        tag=alias.lower()
        s+=f'''dea_fuzzy {q}, inputs({x}) outputs({y}) {bounds} rts({alias}) ort({orientation[0]}) dmu({dmu}) gen(u_{tag}) alpha(0 .5 1)
assert "`r(rts)'"=="{norm}"
foreach a in 0 50 100 {{
 assert u_{tag}_pess_a`a'==f_{tag}_pess_a`a' & u_{tag}_opt_a`a'==f_{tag}_opt_a`a'
}}
'''
    return s+'assert unrelated=="preserve me"\nexport delimited using "results.csv",replace\n'


def read(path):
    with path.open() as f:rows=list(csv.DictReader(f))
    return {k:np.array([float(row[k]) if row[k] not in ['','.'] else np.nan for row in rows]) for k in rows[0] if k.lower() not in ['sid','unrelated','rts']}


def oracle(data,dim,orientation,mode):
    chosen=data['selected']==1
    x=np.column_stack([data[f'x{i}'][chosen] for i in range(1,dim+1)])
    y=np.column_stack([data[f'y{i}'][chosen] for i in range(1,dim+1)])
    xl=np.column_stack([data[f'lo_x{i}'][chosen] for i in range(1,dim+1)]) if mode in ['full','inputs','lower'] else x
    xu=np.column_stack([data[f'hi_x{i}'][chosen] for i in range(1,dim+1)]) if mode in ['full','inputs','upper'] else x
    yl=np.column_stack([data[f'lo_y{i}'][chosen] for i in range(1,dim+1)]) if mode in ['full','outputs','lower'] else y
    yu=np.column_stack([data[f'hi_y{i}'][chosen] for i in range(1,dim+1)]) if mode in ['full','outputs','upper'] else y
    answer={};n=len(x)
    for rts in ['crs','vrs','drs','irs']:
        for alpha in ALPHAS:
            lowx=xl+alpha*(x-xl);highx=xu-alpha*(xu-x)
            lowy=yl+alpha*(y-yl);highy=yu-alpha*(yu-y)
            for sense in ['pess','opt']:
                vals=[]
                for k in range(n):
                    if sense=='pess':
                        x0,y0=highx[k],lowy[k];xx,yy=lowx.copy(),highy.copy()
                    else:
                        x0,y0=lowx[k],highy[k];xx,yy=highx.copy(),lowy.copy()
                    # Evaluated DMU remains in the reference set at its own realization.
                    xx[k]=x0;yy[k]=y0
                    vals.append(solve(x0,y0,xx,yy,orientation,rts))
                answer[(rts,sense,int(100*alpha))]=np.array(vals)
    return answer,n*4*len(ALPHAS)*2


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--baseline-ado',type=Path);ap.add_argument('--only');args=ap.parse_args()
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'current_ado';pkg.mkdir()
    names=package_files(REPO);hashes={}
    for name in names:
        dat=(REPO/name).read_bytes();(pkg/Path(name).name).write_bytes(dat);hashes[name]=hashlib.sha256(dat).hexdigest()
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,numpy=np.__version__,scipy=scipy.__version__,baseline=str(args.baseline_ado)),indent=2))
    results=[]
    def run(name,script,source=pkg,check=None):
        if args.only and not re.search(args.only,name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 120\n'
        for f in ['PLUS','PERSONAL','SITE','OLDPLACE']:
            (d/f).mkdir();head+=f'sysdir set {f} "{d/f}"\n'
        head+=f'adopath ++ "{source}"\n';(d/'test.do').write_text(head+script+'\ndi "TEST_PASS"\n')
        row=dict(test=name,passed=False);start=time.monotonic()
        try:
            subprocess.run([args.stata,'-b','do','test.do'],cwd=d,timeout=180,capture_output=True)
            log=(d/'test.log').read_text(errors='replace')
            if 'TEST_PASS' not in log.splitlines():raise RuntimeError(log[-1800:])
            if check:row.update(check(d))
            row['passed']=True
        except Exception as e:row['error']=str(e)
        row['seconds']=round(time.monotonic()-start,3);results.append(row);(out/'results.json').write_text(json.dumps(results,indent=2));print(json.dumps(row),flush=True)

    for dim in [1,2]:
        for ori in ['input','output']:
            for mode in ['crisp','full','inputs','outputs','lower','upper']:
                for sample in SAMPLES:
                    name=f'{dim}d_{ori}_{mode}_{sample}'
                    def check(d,dim=dim,ori=ori,mode=mode,sample=sample):
                        data=read(d/'results.csv');ref,count=oracle(data,dim,ori,mode);chosen=data['selected']==1;gap=0.
                        for spec,norm in SPECS.items():
                            for alpha in ALPHAS:
                                for sense in ['pess','opt']:
                                    a=data[f'f_{spec}_{sense}_a{int(alpha*100)}'][chosen];b=ref[norm,sense,int(alpha*100)]
                                    gap=max(gap,float(abs(a-b).max()));assert np.isfinite(a).all() and np.max(abs(a-b))<1e-6,(spec,alpha,sense,a,b)
                                    assert np.array_equal(data[f'f_{spec}_{sense}_a{int(alpha*100)}'],data[f'f_{norm}_{sense}_a{int(alpha*100)}'],equal_nan=True)
                        historical=0.
                        if args.baseline_ado and mode in ['crisp','full'] and sample=='full':
                            old=d/'historical';old.mkdir();(old/'test.do').write_text('clear all\nset more off\n'+f'adopath ++ "{args.baseline_ado.resolve()}"\n'+body(dim,ori,mode,sample,old=True)+'\ndi "TEST_PASS"\n')
                            subprocess.run([args.stata,'-b','do','test.do'],cwd=old,timeout=180,capture_output=True)
                            assert 'TEST_PASS' in (old/'test.log').read_text(errors='replace').splitlines()
                            before=read(old/'results.csv')
                            for key,value in data.items():
                                if any(key.startswith('f_'+spec+'_') for spec in ['crs','ccr','vrs','bcc','drs','irs']):
                                    historical=max(historical,float(np.nanmax(abs(value-before[key]))))
                                    # F23 deliberately improves solver precision and scaling.
                                    # The independent oracle above remains the accuracy check.
                                    assert np.allclose(value,before[key],rtol=0,atol=1e-6,equal_nan=True),key
                            changed=max(float(np.nanmax(abs(data[k]-before[k]))) for k in data if k.startswith(('f_nirs_','f_ndrs_')))
                            assert changed>1e-3
                        return dict(lp_count=count,max_oracle_gap=gap,max_historical_canonical_gap=historical)
                    run(name,body(dim,ori,mode,sample),check=check)
    for spec in ['nirs','ndrs']:
        run('replace_'+spec,DATA+f'''dea_fuzzy if id>=2, inputs(x1) outputs(y1) xlo(lo_x1) yhi(hi_y1) rts({spec}) dmu(id) gen(a)
foreach k in pess opt {{
 foreach alpha in 0 25 50 75 100 {{
  clonevar expected_`k'_`alpha'=a_`k'_a`alpha'
 }}
}}
dea_fuzzy if id>=2, inputs(x1) outputs(y1) xlo(lo_x1) yhi(hi_y1) rts({spec}) dmu(id) gen(a) replacevars
foreach k in pess opt {{
 foreach alpha in 0 25 50 75 100 {{
  assert expected_`k'_`alpha'==a_`k'_a`alpha'
 }}
}}
''')
        run('protected_'+spec,DATA+'rename x1 a_pess_a0\n'+unchanged(f'dea_fuzzy, inputs(a_pess_a0) outputs(y1) rts({spec}) dmu(id) gen(a) replacevars'))
    run('bad_rts',DATA+unchanged('dea_fuzzy, inputs(x1) outputs(y1) rts(nonsense) dmu(id)',198))
    summary=dict(total=len(results),passed=sum(x['passed'] for x in results),lp_count=sum(x.get('lp_count',0) for x in results),max_oracle_gap=max((x.get('max_oracle_gap',0) for x in results),default=0))
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
    raise SystemExit(not all(x['passed'] for x in results))

if __name__=='__main__':main()
