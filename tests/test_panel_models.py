#!/usr/bin/env python3
"""Malmquist/window: independent radial LPs, panel matching and summary contracts.

Requires NumPy/SciPy, Stata 16+; optional unmodified R Benchmarking comparison.
"""
from _support import package_files
from pathlib import Path
import argparse,csv,hashlib,itertools,json,re,shutil,subprocess,time
import numpy as np
import scipy
from scipy.optimize import linprog
from test_output_protection import unchanged
REPO=Path(__file__).resolve().parents[1]

def distance(x0,y0,x,y,ori,rts):
    # Radial variable FIRST (production puts it last), balanced physical rows.
    sx=np.maximum(x.max(0),x0);sy=np.maximum(y.max(0),y0)
    sx=np.where(sx>0,sx,1);sy=np.where(sy>0,sy,1)
    x=x/sx;x0=x0/sx;y=y/sy;y0=y0/sy;n=len(x)
    c=np.r_[1 if ori=='input' else -1,np.zeros(n)]
    a=np.vstack([np.c_[-x0 if ori=='input' else np.zeros(len(x0)),x.T],np.c_[np.zeros(len(y0)) if ori=='input' else y0,-y.T]])
    b=np.r_[np.zeros(len(x0)) if ori=='input' else x0,-y0 if ori=='input' else np.zeros(len(y0))]
    v=np.r_[0,np.ones(n)];e=d=None
    if rts=='vrs':e=v[None,:];d=[1]
    if rts in ['drs','irs']:
        sign=1 if rts=='drs' else -1;a=np.vstack([a,sign*v]);b=np.r_[b,sign]
    f=linprog(c,A_ub=a,b_ub=b,A_eq=e,b_eq=d,bounds=(0,None),method='highs',options=dict(primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9))
    if f.status==2:return np.nan
    assert f.success,f.message
    assert (a@f.x-b).max()<2e-7
    if e is not None:assert abs((e@f.x)[0]-1)<2e-7
    return f.x[0] if ori=='input' else 1/f.x[0]

def ratio(a,b):return a/b if np.isfinite(a) and np.isfinite(b) and a>0 and b>0 else np.nan

def avg(v):
    v=np.asarray(v);v=v[np.isfinite(v)];return float(v.mean()) if len(v) else np.nan

def fixture(m,s,variant):
    rng=np.random.default_rng(55324+10*m+s);x=rng.uniform(1,8,(5,m));y=rng.uniform(1,8,(5,s));rows=[]
    for ti,t in enumerate([1,3,7,9]):
        for i in range(5):rows.append([len(rows)+1,i+1,t,*x[i],*(y[i]*(1+.15*ti)),1])
    z=np.array(rows,float)
    if variant=='random':z[:,3:3+m+s]=rng.uniform(1,8,(20,m+s))
    if variant=='unbalanced':z=z[~((z[:,1]==2)&(z[:,2]==3)) & ~((z[:,1]==4)&(z[:,2]==7))]
    if variant=='disjoint':z[z[:,2]==9,1]+=10
    if variant=='gaps':z[:,2]=np.select([z[:,2]==1,z[:,2]==3,z[:,2]==7,z[:,2]==9],[1.25,10.5,101.125,1000.75])
    if variant=='precise_time':z[:,2]=1e9+z[:,2]*.00001
    if variant=='units':z[:,3:3+m]*=np.array([1e-9,1e9])[:m];z[:,3+m:3+m+s]*=np.array([1e8,1e-8])[:s]
    if variant=='tiny':z[:,3:3+m+s]*=1e-12
    if variant=='zero':
        if m>1:z[z[:,1]==1,3]=0
        if s>1:z[z[:,1]==2,3+m]=0
    if variant=='sample':
        z[z[:,1]==3,-1]=0;z[0,-1]=0;z[-1,-1]=0
        z[(z[:,1]==2)&(z[:,2]==3),3]=np.nan;z[np.isnan(z[:,3]),-1]=0
        z[z[:,-1]==0,3+m]=-.5 # invalid excluded data must not cause rejection
    if variant=='permuted':z=z[rng.permutation(len(z))]
    return z,['obsno','id','period']+[f'x{i+1}' for i in range(m)]+[f'y{i+1}' for i in range(s)]+['selected']

def oracle(z,m,s,kind,rts,ori,width):
    sel=z[:,-1]==1;ids=z[:,1];ts=z[:,2];x=z[:,3:3+m];y=z[:,3+m:3+m+s];periods=np.unique(ts[sel]);out={};pairs=[];calls=0
    def dist(k,mask):
        nonlocal calls
        calls+=1;return distance(x[k],y[k],x[mask],y[mask],ori,rts)
    if kind=='window':
        allv=[];attempts=0
        for w in range(len(periods)-width+1):
            ref=sel&np.isin(ts,periods[w:w+width])
            for p in range(width):out[f'res_w{w+1}_p{p+1}']=np.full(len(z),np.nan)
            for k in np.flatnonzero(ref):
                p=list(periods[w:w+width]).index(ts[k])+1;v=dist(k,ref);out[f'res_w{w+1}_p{p}'][k]=v;allv.append(v);attempts+=1
        vals=np.column_stack(list(out.values()));out['res_mean']=np.array([avg(v) for v in vals]);out['res_sd']=np.array([np.nan if not np.isfinite(v).any() else (np.nanstd(v,ddof=1) if np.isfinite(v).sum()>1 else 0.) for v in vals])
        good=np.array(allv)[np.isfinite(allv)]
        meta=[len(periods),len(np.unique(ids[sel])),len(periods)-width+1,width,len(good),avg(good),np.std(good,ddof=1) if len(good)>1 else np.nan,attempts-len(good),attempts,int(np.sum(sel&~np.isfinite(out['res_mean'])))]
    else:
        for name in ['ec','tc','mpi']:out['res_'+name]=np.full(len(z),np.nan)
        for t0,t1 in zip(periods,periods[1:]):
            r0=sel&(ts==t0);r1=sel&(ts==t1);common=np.intersect1d(ids[r0],ids[r1]);vv=[]
            for ident in common:
                k0=np.flatnonzero(r0&(ids==ident))[0];k1=np.flatnonzero(r1&(ids==ident))[0]
                d0=dist(k0,r0);d1=dist(k1,r1);ec=ratio(d1,d0)
                if kind=='global':mpi=ratio(dist(k1,sel),dist(k0,sel));tc=ratio(mpi,ec)
                else:
                    a=ratio(dist(k1,r0),d1);b=ratio(d0,dist(k0,r1));tc=np.sqrt(a*b);mpi=ec*tc
                vals=[ec,tc,mpi];vv.append(vals)
                for name,v in zip(['ec','tc','mpi'],vals):out['res_'+name][k1]=v
            vv=np.array(vv).reshape(-1,3);missing=~np.isfinite(vv);pairs.append([t0,t1,len(common),int(missing.any(1).sum()),*missing.sum(0),*[avg(vv[:,c]) for c in range(3)]])
        last=pairs[-1];meta=[len(periods),last[2],*last[7:10],*last[3:7],sum(p[3] for p in pairs)]
    return out,np.array(meta),np.array(pairs),calls

def read(p):
    with p.open() as f:rows=list(csv.DictReader(f))
    def number(v):
        try:return float(v)
        except ValueError:return np.nan
    return {k:np.array([number(r[k]) for r in rows]) for k in rows[0]}

R_SCRIPT=r'''
suppressPackageStartupMessages(library(Benchmarking))
a=commandArgs(TRUE);jobs=read.csv(a[1]);ans=list();ii=0
put=function(test,z,component,values) {ii<<-ii+1;ans[[ii]]<<-data.frame(test=test,obsno=z$obsno,component=component,value=values)}
for(j in seq_len(nrow(jobs))) {
 q=jobs[j,];z=read.csv(q$file);z=z[z$selected==1,];x=as.matrix(z[paste0('x',1:q$m)]);y=as.matrix(z[paste0('y',1:q$s)])
 storage.mode(x)='double';storage.mode(y)='double';x=sweep(x,2,apply(x,2,max),'/');y=sweep(y,2,apply(y,2,max),'/')
 ori=if(q$ori=='input') 'in' else 'out';times=sort(unique(z$period))
 eff=function(k,ref){v=dea(x[k,,drop=FALSE],y[k,,drop=FALSE],XREF=x[ref,,drop=FALSE],YREF=y[ref,,drop=FALSE],RTS=q$rts,ORIENTATION=ori,FAST=TRUE);if(ori=='out') v=1/v;v[!is.finite(v)|v<=0]=NA;v}
 if(q$kind=='window') {
  for(w in seq_len(length(times)-q$width+1)) {
   mask=z$period %in% times[w:(w+q$width-1)];k=which(mask);v=eff(k,mask)
   for(p in seq_len(q$width)) {at=z$period[k]==times[w+p-1];put(q$test,z[k[at],],paste0('res_w',w,'_p',p),v[at])}
  }
 } else {
  for(t in 2:length(times)) {
   m0=z$period==times[t-1];m1=z$period==times[t];common=intersect(z$id[m0],z$id[m1]);if(!length(common)) next
   k0=match(common,z$id[m0]);k1=match(common,z$id[m1]);i0=which(m0)[k0];i1=which(m1)[k1]
   if(q$kind=='adjacent') {
    fit=malmq(x[m0,,drop=FALSE],y[m0,,drop=FALSE],z$id[m0],x[m1,,drop=FALSE],y[m1,,drop=FALSE],z$id[m1],RTS=q$rts,ORIENTATION=ori,SAMEREF=FALSE)
    ec=fit$ec;tc=fit$tc;mpi=fit$mq
    if(ori=='out'){ec=1/ec;tc=1/tc;mpi=1/mpi}
    # A zero returned by reciprocal conversion of Inf means unavailable here.
    ec[!is.finite(ec)|ec<=0]=NA;tc[!is.finite(tc)|tc<=0]=NA;mpi[!is.finite(mpi)|mpi<=0]=NA
    i1=which(m1)[match(fit$id,z$id[m1])]
   } else {ec=eff(i1,m1)/eff(i0,m0);mpi=eff(i1,rep(TRUE,nrow(z)))/eff(i0,rep(TRUE,nrow(z)));tc=mpi/ec}
   put(q$test,z[i1,],'res_ec',ec);put(q$test,z[i1,],'res_tc',tc);put(q$test,z[i1,],'res_mpi',mpi)
  }
 }
}
write.csv(do.call(rbind,ans),a[2],row.names=FALSE,na='NA')
writeLines(c(paste('Benchmarking',packageVersion('Benchmarking')),capture.output(sessionInfo())),a[3])
'''

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--source-ado',type=Path);ap.add_argument('--only');ap.add_argument('--rscript');a=ap.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir();h={}
    paths=list(a.source_ado.iterdir()) if a.source_ado else [REPO/n for n in package_files(REPO)]
    for p in paths:
        if p.suffix in ['.ado','.do','.sthlp']:raw=p.read_bytes();(pkg/p.name).write_bytes(raw);h[p.name]=hashlib.sha256(raw).hexdigest()
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=h,stata=a.stata,numpy=np.__version__,scipy=scipy.__version__),indent=2)+'\n')
    results=[];jobs=[];expect={}
    def chosen(n):return not a.only or re.search(a.only,n)
    def run(name,body,check=None):
        if not chosen(name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset processors 1\nset type double\nset linesize 160\n'
        for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/key).mkdir();head+=f'sysdir set {key} "{d/key}"\n'
        head+=f'adopath ++ "{pkg}"\n';(d/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n');r=dict(test=name,passed=False);t=time.monotonic()
        try:
            p=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=180)
            log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else p.stderr.decode(errors='replace')
            assert p.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-1800:]
            if check:r.update(check(d))
            r['passed']=True
        except Exception as e:r['error']=repr(e)
        r['seconds']=round(time.monotonic()-t,3);results.append(r);(out/'results.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(r),flush=True)
    variants=['base','random','unbalanced','disjoint','gaps','precise_time','units','tiny','zero','sample','permuted']
    for m,s,variant in itertools.product([1,2],[1,2],variants):
        z,cols=fixture(m,s,variant);file=out/f'data_{m}x{s}_{variant}.csv'
        with file.open('w') as f:
            writer=csv.writer(f);writer.writerow(cols);writer.writerows(['' if np.isnan(v) else format(v,'.17g') for v in row] for row in z)
        xs=' '.join(f'x{i+1}' for i in range(m));ys=' '.join(f'y{i+1}' for i in range(s))
        configs=[(k,r,o,0) for k,r,o in itertools.product(['adjacent','global'],['crs','vrs'],['input','output'])]+[('window',r,o,w) for r,o,w in itertools.product(['crs','vrs','drs','irs'],['input','output'],[2,3,4])]
        for ci,(kind,rts,ori,width) in enumerate(configs):
            name=f'{kind}_{m}x{s}_{variant}_{rts}_{ori}_{width}'
            if not chosen(name):continue
            cmd='dea_window' if kind=='window' else 'dea_malmquist';opt=f'width({width})' if kind=='window' else ('global' if kind=='global' else '')
            head=f'{cmd} {xs}={ys} if selected,' if ci%2 else f'{cmd} if selected,inputs({xs}) outputs({ys})'
            tol='1e-9' if ci%2 else '1e-7'
            body=f'import delimited using "{file}", clear varnames(1) asdouble\n'+f'{head} time(period) dmu(id) rts({rts}) ort({ori}) {opt} tol({tol}) gen(res)\n'
            meta=['n_periods','n_dmus','n_windows','width','n_obs','mean_eff','sd_eff','n_failed','n_attempted','n_failed_obs'] if kind=='window' else ['n_periods','n_dmus','mean_ec','mean_tc','mean_mpi','n_failed','n_failed_ec','n_failed_tc','n_failed_mpi','n_failed_total']
            body+='matrix meta=('+','.join(f'r({v})' for v in meta)+')\n'
            if kind!='window':body+='capture matrix pairs=r(pair_summary)\nif !_rc svmat double pairs\n'
            body+='svmat double meta\nformat _all %24.17g\nexport delimited using "scores.csv",replace datafmt\n'
            def check(d,z=z,m=m,s=s,kind=kind,rts=rts,ori=ori,width=width,name=name,file=file):
                got=read(d/'scores.csv');exp,meta,pairs,lp_count=oracle(z,m,s,kind,rts,ori,width);errors=[];gaps=[];checks=0;miss=0
                for v,e in exp.items():
                    actual=got[v];checks+=len(e);miss+=int((~np.isfinite(e)).sum())
                    if not np.array_equal(np.isfinite(actual),np.isfinite(e)):errors.append([v,'availability',actual.tolist(),e.tolist()]);continue
                    mask=np.isfinite(e);delta=float(np.max(abs(actual[mask]-e[mask]))) if mask.any() else 0;gaps.append(delta)
                    if delta>3e-6:errors.append([v,'value',delta])
                actual=np.array([got[f'meta{i+1}'][0] for i in range(10)])
                if not np.allclose(actual,meta,atol=3e-6,rtol=1e-8,equal_nan=True):errors.append(['summary',actual.tolist(),meta.tolist()])
                if kind!='window':
                    if 'pairs10' not in got:errors.append(['pair_summary','absent'])
                    else:
                        actual_pairs=np.column_stack([got[f'pairs{i+1}'][:len(pairs)] for i in range(10)])
                        if not np.allclose(actual_pairs,pairs,atol=3e-6,rtol=1e-12,equal_nan=True):errors.append(['pair_summary',actual_pairs.tolist(),pairs.tolist()])
                (d/'expected.json').write_text(json.dumps(dict(values={k:v.tolist() for k,v in exp.items()},meta=meta.tolist(),pairs=pairs.tolist(),errors=errors),indent=2)+'\n')
                expect[name]=(z[:,0],exp)
                if variant not in ['precise_time']:
                    jobs.append(dict(test=name,file=str(file),m=m,s=s,kind=kind,rts=rts,ori=ori,width=width))
                assert not errors,str(errors[:2])
                return dict(independent_LPs=lp_count,cells=checks,expected_missing=miss,max_score_gap=max(gaps,default=0))
            run(name,body,check)
    for kind,bad in itertools.product(['malmquist','window'],['0','-1','.','.a','.z']):
        body='input double(id period x y)\n1 1 1 1\n2 1 2 1\n1 2 1 2\n2 2 2 2\nend\n'
        for v in (['res_ec','res_tc','res_mpi'] if kind=='malmquist' else ['res_w1_p1','res_w1_p2','res_mean','res_sd']):body+=f'gen double {v}=123\n'
        body+=unchanged(f'dea_{kind} x=y, dmu(id) time(period) gen(res) replacevars tol({bad})',198)
        run(f'invalid_tol_{kind}_{bad.replace(".","dot")}',body)
    for kind,bad in itertools.product(['malmquist','window'],['negative_x','negative_y','zero_x','zero_y']):
        x,y={'negative_x':('-1','1'),'negative_y':('1','-1'),'zero_x':('0','1'),'zero_y':('1','0')}[bad]
        body=f'input double(id period x y)\n1 1 {x} {y}\n2 1 2 1\n1 2 1 2\n2 2 2 2\nend\n'+unchanged(f'dea_{kind} x=y,dmu(id) time(period) gen(res)',459)
        run(f'invalid_data_{kind}_{bad}',body)
    if a.rscript and jobs:
        with (out/'r_jobs.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(jobs[0]));w.writeheader();w.writerows(jobs)
        (out/'compare.R').write_text(R_SCRIPT)
        p=subprocess.run([a.rscript,str(out/'compare.R'),str(out/'r_jobs.csv'),str(out/'r_scores.csv'),str(out/'R_session.txt')],capture_output=True,timeout=600);(out/'R.log').write_bytes(p.stdout+p.stderr)
        row=dict(test='R_Benchmarking',passed=False)
        try:
            assert p.returncode==0,(p.stdout+p.stderr).decode()[-1800:]
            cells=0;gaps=[]
            with (out/'r_scores.csv').open() as f:
                for rr in csv.DictReader(f):
                    ids,ex=expect[rr['test']];k=list(ids).index(float(rr['obsno']));v=ex[rr['component']][k];actual=float(rr['value']) if rr['value']!='NA' else np.nan
                    assert np.isfinite(v)==np.isfinite(actual),(rr,v)
                    if np.isfinite(v):gaps.append(abs(actual-v));assert gaps[-1]<3e-6,(rr,v)
                    cells+=1
            row.update(passed=True,cells=cells,max_score_gap=max(gaps,default=0))
        except Exception as e:row['error']=repr(e)
        results.append(row);print(json.dumps(row),flush=True)
    (out/'results.json').write_text(json.dumps(results,indent=2)+'\n');summary=dict(total=len(results),passed=sum(r['passed'] for r in results),failed=sum(not r['passed'] for r in results))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary);return bool(summary['failed'])
if __name__=='__main__':raise SystemExit(main())
