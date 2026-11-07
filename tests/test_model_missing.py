#!/usr/bin/env python3
"""Model-specific complete cases and IDEA endpoint preflight, with HiGHS controls.

Requires NumPy, SciPy and licensed Stata 16+. Every case has a fresh Stata
process and isolated ado paths. Historical defects are separate expectations.
"""
from _support import package_files
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import re
import subprocess
import time
import numpy as np
from scipy.optimize import linprog
from test_output_protection import unchanged

REPO=Path(__file__).resolve().parents[1]
DATA='''input double(id x1 x2 y1 y2 b1 b2)
1 1 2 1 2 1 2
2 2 1 3 3 2 3
3 3 4 2 5 3 4
4 4 3 5 4 1 3
5 6 5 7 6 2 4
6 7 4 5 8 3 2
7 9 6 8 7 1 4
8 10 7 11 9 2 3
end
foreach v of varlist x1 x2 y1 y2 {
 gen double `v'l=.8*`v'
 gen double `v'h=1.2*`v'
}
gen double gx1=1
gen double gx2=.5
gen double gy1=1
gen double gy2=.8
gen str20 note="original metadata"
char note[audit] "retained"
'''
QUALIFIERS=[('','1'),('if id!=2','id!=2'),('in 2/8','inrange(_n,2,8)'),('if id!=2 in 1/7','id!=2 & inrange(_n,1,7)')]
SCALARS=['N','N_eff','mean_eff','min_eff','max_eff','n_failed']
BOUNDS={'xboth':'xlo(x1l x2l) xhi(x1h x2h)', 'yboth':'ylo(y1l y2l) yhi(y1h y2h)',
        'both':'xlo(x1l x2l) xhi(x1h x2h) ylo(y1l y2l) yhi(y1h y2h)',
        'xlo':'xlo(x1l x2l)', 'xhi':'xhi(x1h x2h)', 'ylo':'ylo(y1l y2l)', 'yhi':'yhi(y1h y2h)',
        'cross':'xlo(x1l x2l) yhi(y1h y2h)'}


def command(model,form='A',dim=2,q='',opts='',prefix='res',engine='auto',rts='crs',ori='input',dmu='dmu(id)',alias=False):
    x='x1'+(' x2' if dim==2 else '');y='y1'+(' y2' if dim==2 else '')
    call=f'dea {x}={y} {q},' if form=='A' else f'dea {y} {q}, inputs({x})'
    modelopt='' if model=='interval' else f'model({model})'
    return f'{call} {modelopt} {opts} engine({engine}) rts({rts}) orientation({ori}) {dmu} '+('ddfbeta' if alias else 'gen')+f'({prefix})'


def oracle(path,model,dim,rts,directions):
    with path.open() as f: rows=list(csv.DictReader(f))
    rows=[r for r in rows if r['selected']=='1']
    def values(prefix): return np.array([[float(r[prefix+str(i)]) for i in range(1,dim+1)] for r in rows])
    X,Y,B=values('x'),values('y'),values('b');GX=values('gx') if 'gx' in directions else X;GY=values('gy') if 'gy' in directions else Y
    n=len(rows);pred=[]
    for k in range(n):
        if model=='ddf':
            c=np.r_[np.zeros(n),-1.];A=np.vstack([np.c_[X.T,GX[k]],np.c_[-Y.T,GY[k]]]);rhs=np.r_[X[k],-Y[k]];eq=None;erhs=None
            if rts=='vrs':eq=np.array([np.r_[np.ones(n),0.]]);erhs=[1.]
            if rts in ['nirs','ndrs']:
                sign=1 if rts=='nirs' else -1;A=np.vstack([A,np.r_[np.full(n,sign),0.]]);rhs=np.r_[rhs,sign]
            fit=linprog(c,A_ub=A,b_ub=rhs,A_eq=eq,b_eq=erhs,bounds=(0,None),method='highs')
            assert fit.success,fit.message
            pred.append(-fit.fun)
        else:
            # Variables: peer weights z, input/good/bad slacks, scaling t.
            size=n+3*dim+1;eq=[];rhs=[];c=np.zeros(size);c[-1]=1;c[n:n+dim]=-1/(dim*X[k])
            for block,values_,sign in [(0,X,1),(1,Y,-1),(2,B,1)]:
                for j in range(dim):
                    row=np.zeros(size);row[:n]=values_[:,j];row[n+block*dim+j]=sign;row[-1]=-values_[k,j];eq.append(row);rhs.append(0.)
            row=np.zeros(size);row[-1]=1;row[n+dim:n+2*dim]=1/(2*dim*Y[k]);row[n+2*dim:n+3*dim]=1/(2*dim*B[k]);eq.append(row);rhs.append(1.)
            if rts=='vrs':row=np.zeros(size);row[:n]=1;row[-1]=-1;eq.append(row);rhs.append(0.)
            fit=linprog(c,A_eq=np.array(eq),b_eq=rhs,bounds=(0,None),method='highs');assert fit.success,fit.message;pred.append(fit.fun)
    actual=np.array([float(r['res']) for r in rows]);gap=float(np.max(np.abs(actual-pred)))
    assert np.isfinite(actual).all() and gap<1e-6,(gap,actual,pred)
    return dict(lp_count=n,max_oracle_gap=gap)


def paired(setup,model,dim=2,q='',condition='1',required='',form='A',opts='',engine='auto',rts='crs',ori='input',alias=False,interval=False):
    names=['res','res_opt'] if interval else ['res'];scalars=SCALARS+(['mean_opt','n_failed_primary','n_failed_opt'] if interval else [])
    s=DATA+setup+f'gen byte selected=({condition})\n'
    for v in required.split():s+=f'replace selected=0 if missing({v})\n'
    s+='quietly ds\nlocal original_vars "`r(varlist)\'"\ntempfile original\nsave `original\'\npreserve\nkeep if selected\n'
    common=dict(model=model,dim=dim,opts=opts,engine=engine,rts=rts,ori=ori,alias=alias)
    s+=command(**common,form='A')+'\nmatrix Expected=('+','.join('r('+v+')' for v in scalars)+')\n'
    s+='mkmat '+' '.join(names)+', matrix(ExpectedScores)\nrestore\n'
    for v in names:s+=f'gen str20 {v}="old result"\n'
    s+=command(**common,form=form,q=q)+' replacevars saving("saved.dta")\nmatrix Actual=('+','.join('r('+v+')' for v in scalars)+')\n'
    s+='matrix E=r(efficiency)\nassert r(N)==rowsof(ExpectedScores)\n'
    s+='mata: assert(all(st_matrix("Actual"):==st_matrix("Expected")))\n'
    s+='mata: assert(all(st_matrix("ExpectedScores"):==st_data(selectindex(st_data(.,"selected")),tokens("'+' '.join(names)+'"))))\n'
    s+='mata: assert(all(st_matrix("E"):==st_data(selectindex(st_data(.,"selected")),"res")))\n'
    for v in names:s+=f'assert missing({v}) if !selected\n'
    if not interval:s+='assert r(n_failed)==0\nassert !missing(res) if selected\n'
    s+='cf `original_vars\' using `original\',all\nassert "`:char note[audit]\'"=="retained"\n'
    s+='export delimited id selected x1 x2 y1 y2 b1 b2 gx1 gx2 gy1 gy2 res using "scores.csv",replace\n'
    s+='preserve\nkeep if selected\nkeep id '+' '.join(['x1','x2'][:dim]+['y1','y2'][:dim]+(['b1','b2'][:dim] if model=='sbmundesirable' else [])+names)+'\ntempfile expectedsave\nsave `expectedsave\'\nuse "saved.dta",clear\ncf _all using `expectedsave\',all\nrestore\ncapture ds __*\nassert _rc==111\n'
    return s


def reject(setup,cmd,rc):
    return DATA+setup+'''gen str20 res="old main"
gen str20 res_opt="old optimum"
gen str20 res_mid="old midpoint"
save "target.dta"
copy "target.dta" "target_before.dta"
'''+unchanged(cmd+' replacevars saving("target.dta") replace',rc)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stata',required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--baseline-ado',type=Path);p.add_argument('--only');a=p.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'current_ado';pkg.mkdir();hashes={}
    for n in package_files(REPO):
        raw=(REPO/n).read_bytes();(pkg/Path(n).name).write_bytes(raw);hashes[n]=hashlib.sha256(raw).hexdigest()
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,baseline=str(a.baseline_ado),stata=a.stata,numpy=np.__version__,score_tolerance=1e-6),indent=2));rows=[]
    def run(name,script,source=pkg,check=None):
        if a.only and not re.search(a.only,name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 120\n'
        for f in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/f).mkdir();head+=f'sysdir set {f} "{d/f}"\n'
        head+=f'adopath ++ "{source.resolve()}"\n';(d/'test.do').write_text(head+script+'\ndi "TEST_PASS"\n');r=dict(test=name,passed=False);start=time.monotonic()
        try:
            proc=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=240)
            if proc.returncode or not (d/'test.log').exists():(d/'process.json').write_text(json.dumps(dict(rc=proc.returncode,stderr=proc.stderr.decode(errors='replace')),indent=2))
            log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else proc.stderr.decode(errors='replace')
            assert proc.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2200:]
            if (d/'target_before.dta').exists():assert (d/'target.dta').read_bytes()==(d/'target_before.dta').read_bytes()
            if check:r.update(check(d))
            r['passed']=True
        except Exception as e:r['error']=str(e)
        r['seconds']=round(time.monotonic()-start,3);rows.append(r);(out/'results.json').write_text(json.dumps(rows,indent=2));print(json.dumps(r),flush=True)

    for i,(engine,rts,form,directions,missing) in enumerate(itertools.product(['simplex','ipm','auto'],['crs','vrs','nirs','ndrs'],['A','B'],['gx','gy','gx gy'],['.','.a','.z'])):
        dim=1+i%2;d=dim;q,cond=QUALIFIERS[i%4];req=' '.join(f'{s}{j}' for s in ['x','y']+directions.split() for j in range(1,d+1));badvar=directions.split()[-1]+str(d)
        opts=' '.join(g+'('+' '.join(g+str(j) for j in range(1,d+1))+')' for g in directions.split())
        setup=f'replace {badvar}={missing} if id==4\n'+('gsort -id\n' if i%5==0 else '')
        script=paired(setup,'ddf',dim,q,cond,req,form,opts,engine,rts,ori=['input','output'][i%2],alias=i%7==0)
        run(f'ddf_{engine}_{rts}_{form}_{directions.replace(" ","_")}_{missing[1:] or "dot"}',script,check=lambda p,d=d,rts=rts,g=directions:oracle(p/'scores.csv','ddf',d,rts,g))
    for i,(engine,rts,form,missing) in enumerate(itertools.product(['simplex','ipm','auto'],['crs','vrs'],['A','B'],['.','.a','.z'])):
        d=1+i%2;q,cond=QUALIFIERS[i%4];req=' '.join(f'{s}{j}' for s in ['x','y','b'] for j in range(1,d+1));opts='bad('+' '.join('b'+str(j) for j in range(1,d+1))+')'
        script=paired(f'replace b{d}={missing} if id==4\n','sbmundesirable',d,q,cond,req,form,opts,engine,rts)
        run(f'sbmu_{engine}_{rts}_{form}_{missing[1:] or "dot"}',script,check=lambda p,d=d,rts=rts:oracle(p/'scores.csv','sbm',d,rts,''))

    for i,(kind,rts,ori,form) in enumerate(itertools.product(BOUNDS,['crs','vrs','nirs','ndrs'],['input','output'],['A','B'])):
        q,cond=QUALIFIERS[i%4];setup='';req='x1 x2 y1 y2'
        if kind in ['xboth','both']:setup+='replace x1=.z if id==4\n';req=req.replace('x1','')
        if kind in ['yboth','both']:setup+='replace y2=.a if id==5\n';req=req.replace('y2','')
        run(f'interval_{kind}_{rts}_{ori}_{form}',paired(setup,'interval',2,q,cond,req,form,BOUNDS[kind],['simplex','ipm','auto'][i%3],rts,ori,interval=True))
    for var,missing,form,engine in itertools.product(['x1l','x2h','y1l','y2h'],['.','.a','.z'],['A','B'],['simplex','ipm','auto']):
        run(f'invalid_{var}_{missing[1:] or "dot"}_{form}_{engine}',reject(f'replace {var}={missing} in 4\n',command('interval',form,opts=BOUNDS['both'],engine=engine,dmu=''),459))
    for side,form,engine in itertools.product(['x','y'],['A','B'],['simplex','ipm','auto']):
        run(f'reversed_{side}_{form}_{engine}',reject(f'replace {side}1l={side}1h+1 in 4\n',command('interval',form,opts=BOUNDS['both'],engine=engine,dmu=''),459))
    for kind,missing,form in itertools.product(['xlo','xhi','ylo','yhi'],['.','.a','.z'],['A','B']):
        run(f'fallback_{kind}_{missing[1:] or "dot"}_{form}',reject(f'replace {kind[0]}1={missing} in 4\n',command('interval',form,opts=BOUNDS[kind],dmu=''),459))
    for kind,form in itertools.product(['xlo','xhi','ylo','yhi'],['A','B']):
        side=kind[0];other='y' if side=='x' else 'x';v=side+'1'+('l' if kind.endswith('lo') else 'h')
        for mode,q,cond,extra in [('if','if id!=4','id!=4',''),('in','in 5/8','inrange(_n,5,8)',''),('ordinary','','1',f'replace {other}1=.a if id==4\n')]:
            setup=f'replace {v}=.z if id==4\n'+extra;req='x1 x2 y1 y2'
            run(f'excluded_{kind}_{form}_{mode}',paired(setup,'interval',q=q,condition=cond,required=req,form=form,opts=BOUNDS[kind],interval=True))
    for model,form,engine in itertools.product(['ddf','sbmundesirable'],['A','B'],['simplex','ipm','auto']):
        var='gx1' if model=='ddf' else 'b1';opts='gx(gx1 gx2)' if model=='ddf' else 'bad(b1 b2)'
        run(f'empty_{model}_{form}_{engine}',reject(f'replace {var}=.z\n',command(model,form,opts=opts,engine=engine,dmu=''),2000))
    for form in ['A','B']:
        run('empty_interval_'+form,reject('replace y1=.\nreplace x1l=.\n',command('interval',form,opts=BOUNDS['xboth'],dmu=''),2000))
        for side in ['x','y']:
            # One nominal column has both a bounded role and an ordinary role.
            call=f'dea x1=x1,' if form=='A' else 'dea x1, inputs(x1)'
            s=DATA+'replace x1=. in 4\n'+call+f' {side}lo(x1l) {side}hi(x1h) dmu(id) gen(res)\nassert r(N)==7\nassert missing(res) if id==4\n'
            run(f'shared_point_role_{side}_{form}',s)
        s=DATA+'replace b1=. in 4\n'+command('sbmundesirable',form,opts='bad(b1 b2)')+'\n'+('dea_sbmundesirable x1 x2=y1 y2,' if form=='A' else 'dea_sbmundesirable, inputs(x1 x2) outputs(y1 y2)')+' bad(b1 b2) dmu(id) gen(wrapped)\nassert res==wrapped\n'
        run('wrapper_parity_'+form,s)
        run('unused_options_'+form,DATA+'replace gx1=.\nreplace b1=.\n'+command('te',form,opts='gx(gx1 gx2) bad(b1 b2)')+'\nassert r(N)==8\n')
    if a.baseline_ado:
        for name,setup,opts,model,engine,rc,N in [
            ('bad','replace b1=. in 4','bad(b1 b2)','sbmundesirable','auto',3351,None),
            ('bound','replace x1l=. in 4',BOUNDS['xboth'],'interval','auto',3351,None),
            ('reversed','replace x1l=x1h+1 in 4',BOUNDS['xboth'],'interval','auto',0,8),
            ('nominal','replace x1=. in 4',BOUNDS['xboth'],'interval','auto',0,7)]+[
            (g+'_'+e,f'replace {g}1=. in 4',f'{g}({g}1 {g}2)','ddf',e,3351 if e=='ipm' else 0,None if e=='ipm' else 8) for g in ['gx','gy'] for e in ['simplex','ipm','auto']]:
            s=DATA+setup+'\ncapture noisily '+command(model,opts=opts,engine=engine)+f'\nassert _rc=={rc}\n'
            if N is not None:s+=f'assert r(N)=={N}\n'
            run('old_defect_'+name,s,a.baseline_ado)
        for model,engine,rts in itertools.product(['ddf','sbmundesirable','interval'],['simplex','ipm','auto'],['crs','vrs']):
            name=f'healthy_{model}_{engine}_{rts}';opts={'ddf':'gx(gx1 gx2) gy(gy1 gy2)','sbmundesirable':'bad(b1 b2)','interval':BOUNDS['both']}[model];cmd=command(model,opts=opts,engine=engine,rts=rts)
            run(name+'_before',DATA+cmd+'\nkeep res*\nsave "expected.dta"\n',a.baseline_ado)
            run(name,DATA+cmd+f'\ncf res* using "{out/(name+"_before")/"expected.dta"}",all\n')
    summary=dict(total=len(rows),passed=sum(r['passed'] for r in rows),lp_count=sum(r.get('lp_count',0) for r in rows),max_oracle_gap=max((r.get('max_oracle_gap',0) for r in rows),default=0))
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True);raise SystemExit(summary['total']!=summary['passed'])

if __name__=='__main__':main()
