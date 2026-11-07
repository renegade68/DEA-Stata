#!/usr/bin/env python3
"""Convex weak-disposal DDF: independent multiplier dual and mu<=lambda primal,
closure/RTS relations, domains, option combinations and atomic output contracts.
Explicit simplex/native routes and auto fallback; requires NumPy/SciPy; optional R.
"""
from _support import package_files
from pathlib import Path
import argparse, hashlib, itertools, json, re, shutil, subprocess, time
import numpy as np
from scipy.optimize import linprog
from test_core_retries import WRITER, readbin
from test_output_protection import unchanged
from test_undesirable_ddf import oracle as equality_oracle

REPO = Path(__file__).resolve().parents[1]
MATS = 'active idle intensity target_inputs target_outputs target_bad input_slack output_slack diagnostics obs_index'.split()
SMALL = 'input double(id x y b)\n1 1 1 1\n2 1 4 2\n3 2 1 2\nend\n'
BASE = 'model(undesirable) bad(b) disposability(weak)'
RTS = dict(crs='crs', ccr='crs', vrs='vrs', bcc='vrs', drs='drs', nirs='drs', irs='irs', ndrs='irs')


ORACLE_ATTEMPTS = 0
ORACLE_RETRIES = 0

def lp(c, U=None, u=None, E=None, e=None, bounds=(0,None)):
    global ORACLE_ATTEMPTS, ORACLE_RETRIES
    bb = [bounds]*len(c) if isinstance(bounds,tuple) else bounds
    history=[]
    for attempt,(method,pre) in enumerate([('highs-ds',False),('highs-ipm',False),('highs-ds',True),('highs-ipm',True)]):
        ORACLE_ATTEMPTS+=1
        if attempt: ORACLE_RETRIES+=1
        f=linprog(c,A_ub=U,b_ub=u,A_eq=E,b_eq=e,bounds=bounds,method=method,
                  options=dict(presolve=pre,primal_feasibility_tolerance=1e-10,dual_feasibility_tolerance=1e-10))
        valid=f.success and np.isfinite(f.x).all()
        if valid:
            if U is not None:valid &= np.max(U@f.x-u,initial=0)<1e-8
            if E is not None:valid &= np.max(abs(E@f.x-e),initial=0)<1e-8
            for j,(lo,hi) in enumerate(bb):
                if lo is not None:valid &= f.x[j]>=lo-1e-9
                if hi is not None:valid &= f.x[j]<=hi+1e-9
        history.append((method,pre,f.message))
        if valid:return f.fun
    raise AssertionError(history)


def independent(x,y,b,k,direction,rts):
    """Dual variables v,u,w,kappa; separately total lambda and active mu primal."""
    n,m=x.shape; s=y.shape[1]; q=b.shape[1]
    gx=x[k] if direction!='yb' else np.zeros(m)
    gy=y[k] if direction!='x' else np.zeros(s)
    gb=b[k] if direction!='x' else np.zeros(q)
    # Dual: active and idle peer inequalities; free bad prices.
    obj=np.r_[x[k],-y[k],b[k]]
    U=-np.c_[x,-y,b]; idle=-np.c_[x,np.zeros((n,s+q))]
    U=np.r_[U,idle]; vdir=-np.r_[gx,gy,gb]
    bounds=[(0,None)]*(m+s)+[(None,None)]*q
    if rts!='crs':
        U=np.c_[U,-np.ones(2*n)];obj=np.r_[obj,1.];vdir=np.r_[vdir,0]
        bounds += [dict(vrs=(None,None),drs=(0,None),irs=(None,0))[rts]]
    U=np.r_[U,vdir[None,:]];u=np.r_[np.zeros(2*n),-1.]
    dual=lp(obj,U,u,bounds=bounds)
    # Primal with beta first, lambda total, mu active and explicit mu<=lambda.
    U=np.r_[np.c_[gx,x.T,np.zeros((m,n))],np.c_[gy,np.zeros((s,n)),-y.T],
            np.c_[np.zeros(n),-np.eye(n),np.eye(n)]]
    u=np.r_[x[k],-y[k],np.zeros(n)]
    E=np.c_[gb,np.zeros((q,n)),b.T];e=b[k]
    scale=np.r_[0.,np.ones(n),np.zeros(n)]
    if rts=='vrs': E=np.r_[E,scale[None,:]];e=np.r_[e,1.]
    if rts in ['drs','irs']:
        sign=1 if rts=='drs' else -1
        U=np.r_[U,(sign*scale)[None,:]];u=np.r_[u,sign]
    primal=-lp(np.r_[-1.,np.zeros(2*n)],U,u,E,e)
    assert abs(primal-dual)<2e-8*(1+abs(dual)), (k,direction,rts,primal,dual)
    old=equality_oracle(x,y,b,k,direction,rts)
    assert dual>=old-2e-8
    if rts in ['crs','drs']: assert abs(dual-old)<2e-8*(1+abs(dual))
    return dual,old


def audit(d,z,m,s,direction,rts,tol=1e-7):
    attempts_before=ORACLE_ATTEMPTS; retries_before=ORACLE_RETRIES
    a=readbin(d/'scores.bin');mat={key:readbin(d/(key+'.bin')) for key in MATS}
    meta=readbin(d/'meta.bin').ravel();sc=z.max(0);sc[sc==0]=1
    x=z[:,:m]/sc[:m];y=z[:,m:m+s]/sc[m:m+s];b=z[:,m+s:]/sc[m+s:]
    n=len(z);good=np.isfinite(a[:,0]);gap=resid=0.;oracle=[];old=[]
    assert meta[0]==n and meta[1]==sum(~good) and meta[2]==sum(good & (abs(a[:,0])<=tol))
    assert np.array_equal(a[:,1],np.where(good,0,1))
    assert np.array_equal(mat['diagnostics'][:,0],a[:,1])
    assert meta[3]==sum(mat['diagnostics'][:,3]>0) and meta[4]==sum(mat['diagnostics'][:,3]==2)
    assert meta[5]==sum(good & (mat['diagnostics'][:,1]==1)) and meta[6]==sum(good & (mat['diagnostics'][:,1]==2))
    assert meta[7]==sum(mat['diagnostics'][:,2]>0) and meta[8]==sum(mat['diagnostics'][:,2]==2)
    assert meta[9]==sum(mat['diagnostics'][:,1]==2)
    for k in range(n):
        target,previous=independent(x,y,b,k,direction,rts);oracle.append(target);old.append(previous)
        if not good[k]:
            for key in MATS[:8]: assert np.isnan(mat[key][k]).all(),key
            assert np.isnan(mat['diagnostics'][k,4:]).all()
            continue
        beta=a[k,0];g=abs(beta-target)/(1+abs(target));gap=max(gap,g)
        assert g<2e-7,(d.name,k,beta,target,g)
        active=mat['active'][k];idle=mat['idle'][k];lam=mat['intensity'][k]
        assert np.max(abs(lam-active-idle))<1e-12
        gx=x[k] if direction!='yb' else np.zeros(m)
        gy=y[k] if direction!='x' else np.zeros(s)
        gb=b[k] if direction!='x' else np.zeros(b.shape[1])
        tx=x[k]-beta*gx;ty=y[k]+beta*gy;tb=b[k]-beta*gb
        err=max(0.,-beta,-active.min(),-idle.min(),np.max(lam@x-tx),np.max(ty-active@y),np.max(abs(active@b-tb)))
        if rts=='vrs':err=max(err,abs(lam.sum()-1))
        if rts=='drs':err=max(err,lam.sum()-1)
        if rts=='irs':err=max(err,1-lam.sum())
        resid=max(resid,err);assert err<5e-8,(d.name,k,err)
        for key,want,scale in [('target_inputs',tx,sc[:m]),('target_outputs',ty,sc[m:m+s]),('target_bad',tb,sc[m+s:]),('input_slack',tx-lam@x,sc[:m]),('output_slack',active@y-ty,sc[m:m+s])]:
            assert np.max(abs(mat[key][k]/scale-want))<1e-12,key
        assert np.max(mat['diagnostics'][k,7:])<=min(tol,1e-9)
        # Closure witness at the same inputs: active'=alpha active,
        # idle'=idle+(1-alpha)active. No optimization or new data substitution.
        for alpha in [0.,.1,.5,1.]:
            za=alpha*active;wa=idle+(1-alpha)*active
            assert np.max(abs(za+wa-lam))<1e-12
            assert np.max(alpha*ty-za@y)<5e-8
            assert np.max(abs(za@b-alpha*tb))<5e-8
    np.savetxt(d/'independent.csv',np.c_[oracle,old],delimiter=',',header='weak,equality',comments='')
    return dict(independent_LPs=3*n,dual_active_LP_attempts=ORACLE_ATTEMPTS-attempts_before,oracle_retries=ORACLE_RETRIES-retries_before,DMUs=n,missing=int(sum(~good)),max_score_gap=gap,max_normalized_residual=resid,
                changed_from_equality=int(sum(np.array(oracle)-np.array(old)>1e-6)))


R_CODE = r'''
suppressPackageStartupMessages(library(lpSolve))
a<-commandArgs(TRUE);d<-read.csv(a[1]);m<-as.integer(a[2]);s<-as.integer(a[3]);z<-as.matrix(d[,-1]);sc<-apply(z,2,max);sc[sc==0]<-1;z<-sweep(z,2,sc,'/');n<-nrow(z);x<-z[,1:m,drop=FALSE];y<-z[,(m+1):(m+s),drop=FALSE];b<-z[,(m+s+1):ncol(z),drop=FALSE];q<-ncol(b);ans<-c()
for(k in 1:n){
 gx<-if(a[4]=='yb')rep(0,m) else x[k,];gy<-if(a[4]=='x')rep(0,s) else y[k,];gb<-if(a[4]=='x')rep(0,q) else b[k,]
 A<-rbind(cbind(gx,t(x),matrix(0,m,n)),cbind(gy,matrix(0,s,n),-t(y)),cbind(gb,matrix(0,q,n),t(b)),cbind(rep(0,n),-diag(n),diag(n)))
 rhs<-c(x[k,],-y[k,],b[k,],rep(0,n));dir<-c(rep('<=',m+s),rep('=',q),rep('<=',n))
 if(a[5]!='crs'){A<-rbind(A,c(0,rep(1,n),rep(0,n)));rhs<-c(rhs,1);dir<-c(dir,switch(a[5],vrs='=',drs='<=',irs='>='))}
 f<-lp('max',c(1,rep(0,2*n)),A,dir,rhs);stopifnot(f$status==0);ans<-c(ans,f$objval)
}
write.table(ans,a[6],row.names=FALSE,col.names=FALSE,sep=',');cat('Independent transcription via lpSolve',as.character(packageVersion('lpSolve')),'\n')
'''


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--only');ap.add_argument('--rscript');a=ap.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir();hashes={};results=[]
    for n in package_files(REPO):
        data=(REPO/n).read_bytes();(pkg/Path(n).name).write_bytes(data);hashes[n]=hashlib.sha256(data).hexdigest()
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes),indent=2)+'\n')
    def run(name,body,check=None,source=pkg):
        if a.only and not re.search(a.only,name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset type double\nset processors 1\nset linesize 180\n'
        for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/key).mkdir();head+=f'sysdir set {key} "{d/key}"\n'
        (d/'test.do').write_text(head+f'adopath ++ "{source}"\n'+WRITER+body+'\ndi "TEST_PASS"\n');r=dict(test=name,passed=False);t=time.monotonic()
        try:
            p=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=180);log=(d/'test.log').read_text(errors='replace')
            assert p.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-3000:]
            if check:r.update(check(d))
            r['passed']=True
        except Exception as e:r['error']=repr(e)
        r['seconds']=round(time.monotonic()-t,3);results.append(r);(out/'results.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(r),flush=True)
    def dump():
        body='matrix meta=(r(N),r(n_failed),r(N_eff),r(n_native_retry),r(n_native_recovered),r(n_simplex),r(n_native),r(n_simplex_retry),r(n_simplex_recovered),r(n_native_attempted))\nmata: savebin("meta.bin",st_matrix("meta"))\n'
        for key in MATS:body+=f'mata: savebin("{key}.bin",st_matrix("r({key})"))\n'
        return body+'mata: savebin("scores.bin",st_data(.,tokens("res res_status")))\n'
    fixtures=[('analytic',np.array([[1,1,1],[1,4,2],[2,1,2]],float),1,1),('singleton',np.array([[3,4,2]],float),1,1)]
    rng=np.random.default_rng(202611032)
    for m,s,q in [(1,1,1),(2,2,2),(3,2,1)]:
        z=rng.uniform(.4,6,(10,m+s+q));z[-1,:m]=z[0,:m]*2;z[-1,m:m+s]=z[0,m:m+s]*.5;z[-1,m+s:]=z[0,m+s:]*2
        for variant in ['base','units','duplicate','near','zeros','permuted']:
            zz=z.copy()
            if variant=='units':zz*=10.**np.linspace(-9,9,m+s+q)
            if variant in ['duplicate','near']:zz[1]=zz[0]*(1 if variant=='duplicate' else 1+1e-8)
            if variant=='zeros':
                if m>1:zz[::2,0]=0
                if s>1:zz[1::2,m]=0
                zz[::2,m+s]=0
                if q>1:zz[:,-1]=0
            if variant=='permuted':zz=zz[::-1]
            fixtures.append((f'random{m}{s}{q}_{variant}',zz,m,s))
    fixtures.append(('zero_bad',np.array([[1,1,0],[2,4,0],[3,1,0]],float),1,1))
    for key,z,m,s in fixtures:
        xs=' '.join(f'x{i+1}' for i in range(m));ys=' '.join(f'y{i+1}' for i in range(s));bs=' '.join(f'b{i+1}' for i in range(z.shape[1]-m-s))
        f=out/(key+'.csv');np.savetxt(f,np.c_[np.arange(1,len(z)+1),z],delimiter=',',header=','.join(['id']+xs.split()+ys.split()+bs.split()),comments='',fmt='%.17g')
        for rts,direction,engine in itertools.product(['crs','vrs','drs','irs'],['yb','xyb','x'],['simplex','ipm','auto']):
            body=f'import delimited using "{f}",clear asdouble\ndea {xs}={ys},model(undesirable) bad({bs}) disposability(weak) weakmatrix dmu(id) gen(res) rts({rts}) direction({direction}) engine({engine})\n'+dump()
            backend=dict(simplex='package_simplex',ipm='stata_linearprogram',auto='simplex_then_native')[engine]
            body+=f'assert "`r(lp_backend)\'"=="{backend}"\nassert "`r(engine)\'"=="{engine}"\nassert "`r(formulation)\'"=="convex_joint_weak_disposal"\n'
            def check(d,z=z,m=m,s=s,direction=direction,rts=rts,key=key):
                r=audit(d,z,m,s,direction,rts)
                if key in ['analytic','singleton','zero_bad']:assert r['missing']==0
                if key=='analytic' and direction=='yb' and rts in ['vrs','irs']:
                    aa=readbin(d/'scores.bin');assert abs(aa[0,0]-1/3)<1e-8;assert r['changed_from_equality']>0
                return r
            run(f'{key}_{rts}_{direction}_{engine}',body,check)
        if a.rscript and key in ['analytic','random111_base','random222_base','zero_bad']:
            (out/'compare.R').write_text(R_CODE)
            for rts,direction in itertools.product(['crs','vrs','drs','irs'],['yb','xyb','x']):
                def compare(d,f=f,m=m,s=s,rts=rts,direction=direction,key=key):
                    p=subprocess.run([a.rscript,str(out/'compare.R'),str(f),str(m),str(s),direction,rts,str(d/'R.csv')],capture_output=True,text=True,timeout=120);(d/'console.txt').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr
                    rr=np.loadtxt(d/'R.csv',ndmin=1);ss=readbin(out/f'{key}_{rts}_{direction}_ipm'/'scores.bin')[:,0];ok=np.isfinite(ss);gap=np.max(abs(rr[ok]-ss[ok]),initial=0);assert gap<2e-7
                    return dict(R_values=int(sum(ok)),R_native_missing=int(sum(~ok)),R_gap=float(gap))
                run(f'R_{key}_{rts}_{direction}','',compare)
    for syntax in ['dea x=y,','dea y,inputs(x)']:
        for rts,short,stage in itertools.product(RTS,[False,True],[1,2]):
            opts='mod(undesirable) bad(b) disp(weak) weakmat' if short else BASE+' weakmatrix'
            run('syntax_'+str(len(results)),SMALL+syntax+f' {opts} rts({rts}) stage({stage}) dmu(id) gen(res)\nassert res_status==0\nassert abs(res-1/3)<1e-8 in 1\n')
    for opt in ['disposability(no)','weakmatrix','disposability(equality) weakmatrix','disposability(weak) fdh','disposability(weak) bootstrap','disposability(weak) reference','disposability(weak) saveref','disposability(weak) saveweights','disposability(weak) gx(x)','disposability(weak) gy(y)','disposability(weak) w(x)','disposability(weak) p(y)','disposability(weak) trace','disposability(weak) ort(input)','disposability(weak) ort(output)','disposability(weak) stage(3)','disposability(weak) direction(no)','disposability(weak) radialtolerance(.1)','disposability(weak) xlo(x)','disposability(weak) contractcheck','disposability(weak) ddfbeta(foo)']:
        run('reject_'+str(len(results)),SMALL+unchanged('dea x=y,model(undesirable) bad(b) '+opt,198))
    for model in ['te','ddf','congestion','sbmundesirable']:
        for opt in ['disp(weak)','disp(equality)','weakmat']:run('wrong_model_'+model+opt,SMALL+unchanged('dea x=y,model('+model+') bad(b) '+opt,198))
    for v,val in [('x','-1'),('y','-1'),('b','-1'),('x','0'),('y','0')]:
        run('domain_'+v+val,SMALL+f'replace {v}={val} in 2\ngen str8 res="old"\n'+unchanged('dea x=y,'+BASE+' gen(res) replacevars',459))
    for opt,rc in [('gen(x) replacevars',110),('gen(b) replacevars',110),('gen(id) replacevars dmu(id)',110),('gen(_dmu)',110),('gen(abcdefghijklmnopqrstuvwxyzabc)',7)]:
        run('collision_'+str(len(results)),SMALL+unchanged('dea x=y,'+BASE+' '+opt,rc))
    for suffix in ['', '_status']:run('protect'+suffix,SMALL+f'gen str8 res{suffix}="old"\n'+unchanged('dea x=y,'+BASE+' gen(res) dmu(id)',110))
    for tol in ['0','-1','.','.a']:run('tolerance_'+tol,SMALL+unchanged('dea x=y,'+BASE+' tol('+tol+')',198))
    for tol in ['1e-4','1e-9','1e-12']:
        run('accepted_tolerance_'+tol,SMALL+f'dea x=y,{BASE} dmu(id) gen(res) tol({tol})\nassert res_status==0\nassert abs(res-1/3)<1e-8 in 1\n')
    run('empty',SMALL+unchanged('dea x=y if id<0,'+BASE,2000))
    run('missing_id',SMALL+'gen _dmu=1\n'+unchanged('dea x=y,'+BASE,110))
    run('underflow',SMALL+'replace x=1e-200 in 1\nreplace x=1e200 in 2\n'+unchanged('dea x=y,'+BASE,459))
    run('optional_matrix',SMALL+f'dea x=y,{BASE} dmu(id) gen(a)\ndea x=y,{BASE} weakmatrix dmu(id) gen(c)\nassert a==c & a_status==c_status\n')
    run('sample_missing',SMALL+'expand 3\ngen str8 label=string(_n)\nreplace label="" in 2\nreplace x=. in 3\n'+f'dea x=y if id>=1 in 2/8,{BASE} weakmatrix dmu(label) gen(res)\nassert r(N)==5\nmatrix obs=r(obs_index)\nassert obs[1,1]==4 & obs[5,1]==8\nassert missing(res) if _n<4 | _n>8\nassert res_status==0 in 4/8\n')
    run('saving',SMALL+'gen str8 res="old"\n'+f'dea x=y if id>=2,{BASE} gen(res) replacevars saving("space = a, b")\nassert _dmu==id if id>=2\nassert missing(_dmu) if id<2\nassert res_status==0 if id>=2\nlocal avg=r(mean_beta)\npreserve\nuse "space = a, b.dta",clear\nassert _N==2\nconfirm variable b res_status\nsummarize res,meanonly\nassert abs(r(mean)-`avg\')<1e-8\nrestore\n'+unchanged('dea x=y,'+BASE+' dmu(id) gen(res) replacevars saving("space = a, b")',602)+f'dea x=y,{BASE} dmu(id) gen(res) replacevars saving("space = a, b") replace\nassert res_status==0\n')
    run('using',SMALL+'dea x=y using "form A",'+BASE+' gen(res)\nconfirm file "form A.dta"\n')
    run('saving_badpath',SMALL+'gen str8 res="old"\n'+unchanged('dea x=y,'+BASE+' gen(res) replacevars saving("absent/path/result")',603))
    run('saving_write_failure',SMALL+'mkdir target.dta\ngen str8 res="old"\n'+unchanged('dea x=y,'+BASE+' gen(res) replacevars saving(target.dta) replace',608))
    run('matrix_limit',SMALL+'expand 5\nset matsize 10\n'+unchanged('dea x=y,'+BASE+' weakmat',908)+f'dea x=y,{BASE} dmu(id) gen(res)\nassert res_status==0\n')
    # The explicit equality label must remain exactly equivalent to omission.
    for rts,direction,engine in itertools.product(['crs','vrs','drs','irs'],['yb','xyb','x'],['simplex','ipm','auto']):
        run('equality_parity_'+rts+direction+engine,SMALL+f'dea x=y,model(undesirable) bad(b) dmu(id) rts({rts}) direction({direction}) engine({engine}) gen(a)\ndea x=y,model(undesirable) bad(b) dmu(id) rts({rts}) direction({direction}) engine({engine}) disposability(equality) gen(c)\nassert a==c\n')
    for opt in ['disp','dispo','dispos','disposa','disposab','disposabi','disposabil','disposabili','disposabilit','disposability']:
        run('abbrev_'+opt,SMALL+f'dea y,inputs(x) mod(undesirable) bad(b) {opt}(WEAK) rts(BCC) dir(YB) weakmatrix dmu(id) gen(res)\nassert abs(res-1/3)<1e-8 in 1\n')
    for opt in ['weakmat','weakmatr','weakmatri','weakmatrix']:
        run('abbrev_'+opt,SMALL+f'dea x=y,{BASE} {opt} dmu(id) gen(res)\nmatrix Z=r(active)\nassert rowsof(Z)==3\n')
    for fault in ['invalid_primal','all_unavailable']:
        private=out/('fault_'+fault);shutil.copytree(pkg,private);p=private/'_dea_weak.ado';t=p.read_text()
        if fault=='invalid_primal':t=t.replace('checks=_dea_lpcheck(sol,E,e,U,u,lo,hi)','if(k==1) sol=J(1,v,1e8)\n            checks=_dea_lpcheck(sol,E,e,U,u,lo,hi)')
        else:t=t.replace('D[k,2..4]=(backend,sretry,retry)','rc=1\n        D[k,2..4]=(backend,sretry,retry)')
        p.write_text(t)
        checks='assert r(n_failed)==1\nassert missing(res) in 1\n' if fault=='invalid_primal' else 'assert r(n_failed)==3\nassert r(N_eff)==0\nassert missing(r(mean_beta))\nassert missing(res)\n'
        run('injected_'+fault,SMALL+f'dea x=y,{BASE} dmu(id) gen(res) weakmat\n'+checks,source=private)
    # Private faults establish unavailable/ retry behavior without modifying production.
    for fault in ['missing','retry']:
        private=out/('fault_'+fault);shutil.copytree(pkg,private);p=private/'_dea_weak.ado';t=p.read_text()
        if fault=='missing':t=t.replace('D[k,2..4]=(backend,sretry,retry)','if(k==1) rc=1\n        D[k,2..4]=(backend,sretry,retry)')
        else:t=t.replace('rc=_dea_nonte_optimize(S,eps,1e-12)','rc=1')
        p.write_text(t)
        checks='assert r(n_failed)==1\nassert missing(res) in 1\nassert res_status==1 in 1\nassert res_status==0 in 2/3\n' if fault=='missing' else 'assert r(n_native_retry)==3\nassert r(n_native_recovered)==3\nassert res_status==0\nassert abs(res-1/3)<1e-8 in 1\n'
        run('injected_'+fault,SMALL+f'dea x=y,{BASE} engine(ipm) rts(vrs) weakmat dmu(id) gen(res)\n'+checks,source=private)
    for fault in ['simplex_retry','auto_fallback']:
        private=out/('fault_'+fault);shutil.copytree(pkg,private);p=private/'_dea_weak.ado';t=p.read_text()
        if fault=='simplex_retry':t=t.replace('rc=res.status','rc=4',1)
        else:t=t.replace('if(engine=="ipm" | (engine=="auto" & rc!=0))','rc=4\n        if(engine=="ipm" | (engine=="auto" & rc!=0))',1)
        p.write_text(t)
        checks='assert r(n_simplex_retry)==3\nassert r(n_simplex_recovered)==3\n' if fault=='simplex_retry' else 'assert r(n_native_attempted)==3\nassert r(n_native)==3\nassert r(n_simplex)==0\n'
        run('injected_'+fault,SMALL+f'dea x=y,{BASE} engine(auto) rts(vrs) weakmat dmu(id) gen(res)\n'+checks+'assert res_status==0\nassert abs(res-1/3)<1e-8 in 1\n',source=private)
    summary=dict(total=len(results),passed=sum(r['passed'] for r in results),failed=sum(not r['passed'] for r in results))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary);return bool(summary['failed'])

if __name__=='__main__':raise SystemExit(main())
