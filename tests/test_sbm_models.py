#!/usr/bin/env python3
"""SBM/super-SBM: independent projection LPs and unmodified deaR comparison.

The ordinary oracle eliminates slacks; super-SBM uses dimensionless projection
ratios. Zero outputs follow the documented structural-zero convention. deaR
perturbs zeros, so only strictly positive data are compared with that package.
"""
from _support import package_files
from pathlib import Path
import argparse, csv, hashlib, itertools, json, os, re, subprocess, time
import numpy as np
import scipy
from scipy.optimize import linprog

REPO = Path(__file__).resolve().parents[1]
RTS = ['crs', 'vrs', 'drs', 'irs']
ORT = ['', 'input', 'output']
ENG = ['simplex', 'ipm', 'auto']

def lp(c, u, b, e, d):
    fit = linprog(c, A_ub=np.array(u) if len(u) else None,
                  b_ub=np.array(b) if len(u) else None,
                  A_eq=np.array(e) if len(e) else None,
                  b_eq=np.array(d) if len(e) else None, bounds=(0, None),
                  method='highs', options=dict(primal_feasibility_tolerance=1e-9,
                                              dual_feasibility_tolerance=1e-9))
    if fit.status == 2:
        return np.nan
    assert fit.success, fit.message
    if len(u): assert np.max(np.array(u) @ fit.x-b) < 2e-7
    if len(e): assert np.max(abs(np.array(e) @ fit.x-d)) < 2e-7
    return fit.fun

def ordinary(x, y, k, rts, ort):
    # Only lambda (and t for the fractional model): no production slack LP.
    x=x/np.where(x.max(axis=0)>0,x.max(axis=0),1)
    y=y/np.where(y.max(axis=0)>0,y.max(axis=0),1)
    n=len(x); xp=x[k]>0; yp=y[k]>0
    cx=np.mean(x[:,xp]/x[k,xp],axis=1)
    cy=np.mean(y[:,yp]/y[k,yp],axis=1)
    if ort:
        c=cx if ort=='input' else -cy
        u=list(np.vstack([x.T,-y.T]));b=list(np.r_[x[k],-y[k]])
        rt=np.ones(n);rb=1.;e=[];d=[]
    else:
        c=np.r_[cx,0]
        u=list(np.vstack([np.c_[x.T,-x[k]],np.c_[-y.T,y[k]]]))
        b=[0.]*(x.shape[1]+y.shape[1]);e=[np.r_[cy,0]];d=[1.]
        rt=np.r_[np.ones(n),-1];rb=0.
    if rts=='vrs':e.append(rt);d.append(rb)
    if rts in ['drs','irs']:
        sign=1 if rts=='drs' else -1;u.append(sign*rt);b.append(sign*rb)
    val=lp(c,u,b,e,d)
    return -1/val if ort=='output' else val

def super_score(x, y, k, rts, ort):
    # Relative projected inputs a and outputs b, independent of physical units.
    xr=np.delete(x,k,axis=0)/x[k];yr=np.delete(y,k,axis=0)/y[k]
    n=len(xr);m=x.shape[1];s=y.shape[1];nv=n+m+s+1
    if not n:return np.nan
    def row():return np.zeros(nv)
    c=row();u=[];b=[];e=[];d=[]
    if ort=='output':c[n+m:n+m+s]=-1/s
    else:c[n:n+m]=1/m
    for j in range(m):
        v=row();v[:n]=xr[:,j];v[n+j]=-1;u.append(v);b.append(0)
        v=row();v[n+j]=-1;v[-1]=1;u.append(v);b.append(0)
    for j in range(s):
        v=row();v[:n]=-yr[:,j];v[n+m+j]=1;u.append(v);b.append(0)
        v=row();v[n+m+j]=1;v[-1]=-1;u.append(v);b.append(0)
    if ort:
        v=row();v[-1]=1;e.append(v);d.append(1)
        for j in range(s if ort=='input' else m):
            v=row();v[n+m+j if ort=='input' else n+j]=1;e.append(v);d.append(1)
    else:
        v=row();v[n+m:n+m+s]=1/s;e.append(v);d.append(1)
    v=row();v[:n]=1;v[-1]=-1
    if rts=='vrs':e.append(v);d.append(0)
    if rts in ['drs','irs']:u.append(v if rts=='drs' else -v);b.append(0)
    val=lp(c,u,b,e,d)
    # Zero optimum in output super-SBM has no finite reciprocal score.
    return (-1/val if val<0 else np.nan) if ort=='output' else val

def oracle(x,y,rts,ort,model,tol=1e-7):
    out=[]
    for k in range(len(x)):
        v=ordinary(x,y,k,rts,ort)
        if model=='supersbm' and np.isfinite(v) and v>=1-tol:
            v=super_score(x,y,k,rts,ort)
        out.append(v)
    return np.array(out)

def fixture(m,s,variant):
    x=np.array([[4,3],[6,3],[8,1],[8,1],[2,4]],float)[:,:m]
    y=np.array([[2,3],[2,3],[6,2],[6,1],[1,4]],float)[:,:s]
    if variant=='random':
        rng=np.random.default_rng(20413+10*m+s)
        x=rng.uniform(.5,9,(9,m));y=rng.uniform(.5,9,(9,s))
    if variant=='units':x*=np.array([1e-9,1e9])[:m];y*=np.array([1e8,1e-8])[:s]
    if variant=='tiny':x*=1e-12;y*=1e-12
    if variant=='zero':
        if m>1:x[:2,0]=0
        if s>1:y[:2,0]=0
    if variant=='allzero_column' and m>1:x[:,0]=0
    if variant=='allzero_output' and s>1:y[:,0]=0
    if variant=='singleton':x=x[:1];y=y[:1]
    z=np.c_[np.arange(1,len(x)+1),x,y]
    if variant=='permuted':z=z[[4,1,3,0,2]]
    if variant=='duplicate':z=np.vstack([z,z[2]]);z[-1,0]=6
    return z,['id']+[f'x{j+1}' for j in range(m)]+[f'y{j+1}' for j in range(s)]

def read(p):
    with p.open() as f: rows=list(csv.DictReader(f))
    return {k:np.array([float(r[k]) if r[k] and not r[k].startswith('.') or
                       (r[k].startswith('.') and len(r[k])>1 and r[k][1].isdigit()) else np.nan
                       for r in rows]) for k in rows[0]}

R_SCRIPT=r'''
suppressPackageStartupMessages(library(deaR))
a=commandArgs(TRUE);jobs=read.csv(a[1]);ans=list();ii=0
for(j in seq_len(nrow(jobs))) {
 z=jobs[j,];df=read.csv(z$file);df=df[df$selected==1,]
 quantities=c(paste0('x',1:z$m),paste0('y',1:z$s))
 df[quantities]=lapply(df[quantities],as.numeric)
 data=make_deadata(df,dmus=1,inputs=match(paste0('x',1:z$m),names(df)),outputs=match(paste0('y',1:z$s),names(df)))
 for(rt in c('crs','vrs','nirs','ndrs')) for(ori in c('no','io','oo')) {
  reg=model_sbmeff(data,orientation=ori,rts=rt,silent=TRUE,compute_target=FALSE)
  v=efficiencies(reg)
  raw=suppressWarnings(model_sbmsupereff(data,orientation=ori,rts=rt,silent=TRUE,compute_target=FALSE))
  su=efficiencies(raw);combined=ifelse(v<1-1e-7,v,su)
  ii=ii+1;ans[[ii]]=data.frame(test=z$test,rts=rt,ort=ori,id=df$id,sbm=v,supersbm=combined,raw_super=su)
 }
}
write.csv(do.call(rbind,ans),a[2],row.names=FALSE,na='NA')
writeLines(c(paste('deaR',packageVersion('deaR')),capture.output(sessionInfo())),a[3])
'''

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata',required=True);ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--source-ado',type=Path);ap.add_argument('--only');ap.add_argument('--rscript')
    ap.add_argument('--baseline-ado',type=Path,help='Old engines for an in-session reload check')
    a=ap.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    pkg=out/'ado';pkg.mkdir();src=a.source_ado or REPO/'ado';hashes={}
    if a.source_ado:
        paths=list(src.iterdir())
    else:
        names=package_files(REPO)
        paths=[REPO/name for name in names]
    for p in paths:
        if p.is_file() and p.suffix in ['.ado','.do','.sthlp']:
            raw=p.read_bytes();(pkg/p.name).write_bytes(raw);hashes[p.name]=hashlib.sha256(raw).hexdigest()
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,numpy=np.__version__,scipy=scipy.__version__,stata=a.stata),indent=2)+'\n')
    results=[];jobs=[];expected={}
    def chosen(name):return not a.only or re.search(a.only,name)
    def run(name,body,check=None):
        if not chosen(name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset processors 1\nset type double\nset linesize 160\n'
        for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:
            (d/key).mkdir();head+=f'sysdir set {key} "{d/key}"\n'
        head+=f'adopath ++ "{pkg}"\n';(d/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n')
        row=dict(test=name,passed=False);t=time.monotonic()
        try:
            p=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=240)
            log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else p.stderr.decode(errors='replace')
            assert p.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2000:]
            if check:row.update(check(d))
            row['passed']=True
        except Exception as ex:row['error']=repr(ex)
        row['seconds']=round(time.monotonic()-t,3);results.append(row)
        (out/'results.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(row),flush=True)
    variants=['base','random','permuted','duplicate','units','tiny','sample','zero','allzero_column','allzero_output','singleton']
    for m,s,variant in itertools.product([1,2],[1,2],variants):
        name=f'grid_{m}x{s}_{variant}'
        if not chosen(name):continue
        z,cols=fixture(m,s,variant);f=out/(name+'.csv')
        np.savetxt(f,z,delimiter=',',header=','.join(cols),comments='',fmt='%.17g')
        xs=' '.join(f'x{j+1}' for j in range(m));ys=' '.join(f'y{j+1}' for j in range(s))
        q='if id!=2 in 1/5' if variant=='sample' else ''
        body=f'import delimited using "{f}",clear varnames(1) asdouble\ngen byte selected=0\nreplace selected=1 {q}\n'
        if variant=='sample':body+='replace x1=.z if id==4\nreplace selected=0 if id==4\n'
        models=['sbm'] if variant in ['zero','allzero_column','allzero_output'] else ['sbm','supersbm']
        for ei,ri,oi,mi in itertools.product(range(3),range(4),range(3),range(len(models))):
            model=models[mi];ort=ORT[oi];opt=f'ort({ort})' if ort else ''
            rts=RTS[ri] if ei%2 else ['CCR','BCC','NIRS','NDRS'][ri]
            spec=f'dea {xs}={ys} {q},' if (ei+oi)%2 else f'dea {ys} {q}, inputs({xs})'
            v=f'e{ei}_{ri}_{oi}_{mi}';tol='1e-7' if ei%2 else '1e-9'
            body+=f'{spec} model({model}) rts({rts}) {opt} eng({ENG[ei]}) tol({tol}) stage({1+(ei+oi)%2}) dmu(id) gen({v})\n'
            body+=f'local failures=r(n_failed)\nlocal neff=r(N_eff)\n'
            body+=f'assert "`r(orientation)\'"=="{ort or "non-oriented"}"\n'
            body+=f'count if selected & missing({v})\nassert r(N)==`failures\'\n'
            criterion=f'{v}>=1-{tol}' if model=='supersbm' else f'abs({v}-1)<={tol}'
            body+=f'count if selected & !missing({v}) & {criterion}\nassert r(N)==`neff\'\nassert missing({v}) if !selected\n'
        body+='format _all %24.17g\nexport delimited using "scores.csv",replace datafmt\n'
        def check(d,m=m,s=s,name=name,models=models,variant=variant):
            got=read(d/'scores.csv');sel=got['selected']==1
            x=np.column_stack([got[f'x{j+1}'][sel] for j in range(m)])
            y=np.column_stack([got[f'y{j+1}'][sel] for j in range(s)])
            ex={};gaps=[];errors=[];missing=0;checks=0
            for ri,oi,mi in itertools.product(range(4),range(3),range(len(models))):
                v=oracle(x,y,RTS[ri],ORT[oi],models[mi]);ex[(ri,oi,mi)]=v
                for ei in range(3):
                    val=got[f'e{ei}_{ri}_{oi}_{mi}'][sel];checks+=len(v)
                    if not np.array_equal(np.isfinite(val),np.isfinite(v)):
                        errors.append([ei,ri,oi,mi,'availability',val.tolist(),v.tolist()]);continue
                    finite=np.isfinite(v);missing+=int((~finite).sum())
                    gap=float(np.max(abs(val[finite]-v[finite]))) if finite.any() else 0.
                    gaps.append(gap)
                    if gap>3e-6:errors.append([ei,ri,oi,mi,'score',gap,val.tolist(),v.tolist()])
            (d/'comparison.json').write_text(json.dumps(dict(errors=errors,max_gap=max(gaps,default=0),score_checks=checks,expected_missing=missing),indent=2)+'\n')
            expected[name]=dict(ids=got['id'][sel],ex=ex)
            if (x>0).all() and (y>0).all() and len(x)>1 and len(models)==2 and variant not in ['units','tiny']:
                jobs.append(dict(test=name,file=str(d/'scores.csv'),m=m,s=s))
            assert not errors,f'{len(errors)} comparisons failed; see comparison.json: {errors[:1]}'
            return dict(commands=36*len(models),score_checks=checks,expected_missing=missing,max_score_gap=max(gaps,default=0))
        run(name,body,check)
    # Exact counterexamples for zero-term denominators and reference constraints.
    for eng in ENG:
        body='''input double(id x1 x2 y1 y2)
1 0 2 0 1
2 0 1 0 2
3 1 .1 9 9
end
'''
        for ori,score in [('',.25),('input',.5),('output',.5)]:
            opt=f'ort({ori})' if ori else ''
            body+=f'dea x1 x2=y1 y2, model(sbm) rts(vrs) {opt} eng({eng}) gen(e) dmu(id) replacevars\nassert abs(e[1]-{score})<2e-7\n'
        run('analytic_zero_'+eng,body)
        body='''input double(id x y)
1 1 2
2 2 1
end
'''
        for ori in ORT:
            opt=f'ort({ori})' if ori else ''
            body+=f'dea x=y,model(supersbm) rts(vrs) {opt} eng({eng}) gen(e) dmu(id) replacevars\n'
            body+=('assert abs(e[1]-4)<2e-7\n' if not ori else 'assert missing(e[1])\n')
        run('analytic_super_'+eng,body)
        # Invalid rows outside if/in or excluded for another missing quantity
        # must not invalidate the complete selected sample.
        body='''input double(id x1 x2 y1 y2)
1 1 2 1 1
2 2 1 2 1
3 -2 0 -1 0
4 0 .z 0 0
end
'''
        for model in ['sbm','supersbm']:
            body+=f'dea x1 x2=y1 y2 if id!=3 in 1/4,model({model}) engine({eng}) dmu(id) gen(e) replacevars\nassert r(N)==2\nassert missing(e) if id>=3\n'
        run('sample_guard_'+eng,body)
        # stage() changes neither SBM objective nor the super-SBM screening.
        body='''input double(id x1 x2 y1 y2)
1 4 3 2 3
2 6 3 2 3
3 8 1 6 2
4 8 1 6 1
5 2 4 1 4
end
'''
        for model,ort in itertools.product(['sbm','supersbm'],ORT):
            opt=f'ort({ort})' if ort else ''
            for stage in [1,2]:
                body+=f'dea x1 x2=y1 y2,model({model}) {opt} rts(vrs) engine({eng}) stage({stage}) dmu(id) gen(e{stage}) replacevars\n'
            body+='assert e1==e2\n'
        run('stage_parity_'+eng,body)
        if eng!='auto':
            # The optional internal slack interface must undo column scaling.
            body='''input double(id x1 x2 y1 y2)
1 0 2e9 0 1e-9
2 0 1e9 0 2e-9
3 1 .1e9 9e-9 9e-9
end
gen byte selected=1
foreach v in e si1 si2 so1 so2 {
 gen double `v'=.
}
'''
            body+=f'quietly run "{pkg}/_dea_{eng}.ado"\n'
            orientations=ORT if eng=='simplex' else ['']
            for ort in orientations:
                args='"x1 x2", "y1 y2", "e", "si1 si2", "so1 so2", "selected", "vrs"'
                if eng=='simplex':args+=f', "{ort}"'
                body+=f'mata: _dea_{eng}_sbm({args}, 1e-9)\n'
                score=.25 if not ort else .5
                body+=f'assert abs(e[1]-{score})<2e-7\nassert abs(si1[1])<1e-8 & abs(si2[1]/1e9-1)<2e-7\nassert abs(so1[1]/1e-9)<2e-7 & abs(so2[1]/1e-9-1)<2e-7\n'
            run('slack_units_'+eng,body)
        if a.baseline_ado and eng!='auto':
            old=a.baseline_ado.resolve()
            body=f'quietly run "{old}/_dea_{eng}.ado"\n'
            body+='input double(id x1 x2 y)\n1 0 2 1\n2 0 1 1\nend\n'
            marker=re.search(r'real scalar _dea_'+eng+r'_audit_ready\(\) return\((\d+)\)', (pkg/('_dea_'+eng+'.ado')).read_text())
            assert marker, 'Engine readiness marker is missing'
            expected_ready=int(marker.group(1))
            body+=f'dea x1 x2=y,model(sbm) rts(vrs) engine({eng}) dmu(id) gen(e)\nassert abs(e[1]-.5)<1e-7\nmata: assert(_dea_{eng}_audit_ready()=={expected_ready})\n'
            run('cached_engine_'+eng,body)
    # Early data validation must protect existing results and files.
    for model,eng,kind in itertools.product(['sbm','supersbm'],ENG,['negative_x','negative_y','zero_x','zero_y','zero_element']):
        if model=='sbm' and kind=='zero_element':continue
        changes={'negative_x':'replace x1=-1 in 1','negative_y':'replace y1=-1 in 1','zero_x':'replace x1=0 in 1\nreplace x2=0 in 1','zero_y':'replace y1=0 in 1\nreplace y2=0 in 1','zero_element':'replace x1=0 in 1'}
        body='clear\nset obs 3\ngen double x1=_n\ngen double x2=4-_n\ngen double y1=2\ngen double y2=3\ngen double e=123\n'+changes[kind]+'\n'
        body+=f'capture noisily dea x1 x2=y1 y2, model({model}) engine({eng}) gen(e) replacevars saving("protected.dta")\nassert _rc==459\nassert e==123\nconfirm new variable _dmu\ncapture confirm file "protected.dta"\nassert _rc!=0\n'
        run(f'guard_{model}_{eng}_{kind}',body)
    if a.rscript and jobs:
        with (out/'r_jobs.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=['test','file','m','s']);w.writeheader();w.writerows(jobs)
        (out/'compare.R').write_text(R_SCRIPT)
        env=dict(os.environ,RGL_USE_NULL='TRUE')
        p=subprocess.run([a.rscript,str(out/'compare.R'),str(out/'r_jobs.csv'),str(out/'r_results.csv'),str(out/'R_session.txt')],capture_output=True,env=env,timeout=240)
        (out/'R.log').write_bytes(p.stdout+p.stderr);row=dict(test='R_deaR',passed=False)
        try:
            assert p.returncode==0,(p.stdout+p.stderr).decode()[-2500:]
            count=0;gaps=[]
            with (out/'r_results.csv').open() as f:
                for r in csv.DictReader(f):
                    ex=expected[r['test']];idx=list(ex['ids']).index(float(r['id']))
                    ri=['crs','vrs','nirs','ndrs'].index(r['rts']);oi=['no','io','oo'].index(r['ort'])
                    for mi,model in enumerate(['sbm','supersbm']):
                        v=ex['ex'][(ri,oi,mi)][idx];val=float(r[model]) if r[model]!='NA' else np.nan
                        assert np.isfinite(val)==np.isfinite(v),(r,v)
                        if np.isfinite(v):gaps.append(abs(val-v));assert gaps[-1]<3e-6,(r,v)
                        count+=1
            row.update(passed=True,score_checks=count,max_score_gap=max(gaps))
        except Exception as ex:row['error']=repr(ex)
        results.append(row);print(json.dumps(row),flush=True)
    (out/'results.json').write_text(json.dumps(results,indent=2)+'\n')
    summary=dict(total=len(results),passed=sum(r['passed'] for r in results),failed=sum(not r['passed'] for r in results))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)
    return bool(summary['failed'])

if __name__=='__main__':raise SystemExit(main())
