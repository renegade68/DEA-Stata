#!/usr/bin/env python3
"""Specialized-model failures, component statistics, saved samples and independent LPs.

Requires NumPy/SciPy, licensed Stata 16+, and a complete pre-R06 ado snapshot.
Faults exist only in isolated test copies. Infeasible copies append 0=1;
original positive self-referencing data have an explicit feasible witness.
"""
from _support import package_files
from pathlib import Path
import argparse, hashlib, itertools, json, re, shutil, subprocess, time
import numpy as np
from scipy.optimize import linprog
from test_core_retries import WRITER, readbin
from verify_bootstrap_lp import solve as radial_solve
REPO=Path(__file__).resolve().parents[1]
DATA='''input double id x1 x2 y1 y2 b1 b2
1 1 2 1 2 1 2
2 2 1 3 3 2 3
3 3 4 2 5 3 4
4 4 3 5 4 1 3
5 6 5 7 6 2 4
6 7 4 5 8 3 2
end
'''
INIT='scalar failure_mode=1\nscalar fault_cong=0\nscalar fault_sbmu=0\nscalar fault_crs=0\nscalar fault_vrs=0\nscalar fault_value=.\n'
def instrument(src,dst):
 shutil.copytree(src,dst);p=dst/'dea.ado';s=p.read_text()
 for func,flag,rhs in [('_dea_congestion_lp','cong','b_eq'),('_dea_sbmu_solve','sbmu','b_eq_v')]:
  start=s.index('void '+func+'(');end=s.index('\nend\n',start);a=s[start:end]
  cond=f'(st_numscalar("fault_{flag}")==-1 | st_numscalar("fault_{flag}")==k)'
  anchor='        rc = _dea_nonte_optimize(S,tol)';assert a.count(anchor)==1
  a=a.replace(anchor,f'        if (st_numscalar("failure_mode")==2 & {cond}) S.setEquality((A_eq\\J(1,cols(A_eq),0)),({rhs}\\1))\n'+anchor)
  anchor='if (rc == 0)';assert a.count(anchor)==1
  a=a.replace(anchor,f'if (rc == 0 & !(st_numscalar("failure_mode")==1 & {cond}))')
  anchor='    // Store in Stata variable';assert a.count(anchor)==1
  a=a.replace(anchor,f'    for (k=1;k<=n;k++) if (st_numscalar("failure_mode")==3 & {cond}) scores[k]=st_numscalar("fault_value")\n\n'+anchor)
  s=s[:start]+a+s[end:]
 p.write_text(s)
 p=dst/'_dea_ipm.ado';full=p.read_text();start=full.index('void _dea_ipm_envelopment(');end=full.index('real scalar _dea_ipm_slacks(',start);s=full[start:end];anchor='        st_store(obs_k, st_varindex(eff_var), eff_score)';assert s.count(anchor)==1
 s=s.replace(anchor,'        if ((rts=="crs" & (st_numscalar("fault_crs")==-1 | st_numscalar("fault_crs")==k)) |\n            (rts=="vrs" & (st_numscalar("fault_vrs")==-1 | st_numscalar("fault_vrs")==k))) eff_score=st_numscalar("fault_value")\n'+anchor);p.write_text(full[:start]+s+full[end:])

FIXTURES={'regular':DATA,'duplicate':DATA+'expand 2 if id==2\nreplace id=_n\n',
          'collinear':DATA+'replace x2=2*x1\nreplace y2=3*y1\nreplace b2=2*b1\n',
          'units':DATA+'replace x1=x1*100\nreplace y2=y2/10\nreplace b1=b1*10\n'}
FIELDS=['N','N_eff','mean_eff','min_eff','max_eff','n_failed']
CFIELDS=['n_congested','mean_vrs','n_failed_primary','n_failed_crs','n_failed_vrs','n_failed_cong','n_failed_se']

def command(model,route,rts='crs',engine='ipm',sample='',save=False,extra=''):
    opts=f'dmu(id) gen(e) tolerance(1e-7) {extra}'
    if route.startswith('core'):
        io=f'x1 x2=y1 y2 {sample},' if route=='coreA' else f'y1 y2 {sample},inputs(x1 x2)'
        opts+=f' model({model}) rts({rts}) engine({engine})'
        if save:opts+=' saving("saved.dta")'
        name='dea'
    else:
        io=f'x1 x2=y1 y2 {sample},' if route=='wrapA' else f'{sample},inputs(x1 x2) outputs(y1 y2)'
        name='dea_'+model
        if model=='congestion':opts+=f' engine({engine})'
        else:opts+=f' rts({rts})'
    if model=='sbmundesirable':opts+=' bad(b1 b2)'
    return f'{name} {io} {opts}\n'

def checks(model,route,save=False):
    fields=FIELDS+(CFIELDS if model=='congestion' else [])
    s='matrix Stats=('+','.join(f'r({v})' for v in fields)+')\n'
    s+='matrix E=r(efficiency)\nlocal actual_model "`r(model)\'"\n'
    if route.startswith('wrap'):
        s+='assert r(n_efficient)==r(N_eff)\nassert r(mean_efficiency)==r(mean_eff)\n'
    s+=f'assert "`actual_model\'"=="{model}"\n'
    s+='capture matrix list r(stage_status)\nassert _rc==111\ncapture matrix list r(solution_source)\nassert _rc==111\n'
    s+='''quietly count if selected
assert Stats[1,1]==r(N) & rowsof(E)==r(N)
mata: assert(all(st_matrix("E"):==st_data(.,"e","selected")))
quietly count if selected & !missing(e) & abs(e-1)<=1e-7
assert Stats[1,2]==r(N)
quietly summarize e if selected
assert (Stats[1,3]==r(mean) | reldif(Stats[1,3],r(mean))<1e-12) & (Stats[1,4]==r(min) | reldif(Stats[1,4],r(min))<1e-12) & (Stats[1,5]==r(max) | reldif(Stats[1,5],r(max))<1e-12)
assert missing(e) if !selected
'''
    vars='e'
    if model=='congestion':
        vars+=' e_crs e_vrs e_cong_eff e_congestion e_se'
        s+='''quietly count if selected & missing(e,e_crs,e_vrs,e_cong_eff,e_se)
assert Stats[1,6]==r(N)
quietly count if selected & e<1-1e-7
assert Stats[1,7]==r(N)
quietly summarize e_vrs if selected
assert Stats[1,8]==r(mean) | reldif(Stats[1,8],r(mean))<1e-12
'''
        for i,v in enumerate(['e','e_crs','e_vrs','e_cong_eff','e_se'],9):
            s+=f'quietly count if selected & missing({v})\nassert Stats[1,{i}]==r(N)\n'
        s+='''assert e==e_congestion
assert e==e_vrs/e_cong_eff if selected & !missing(e_vrs,e_cong_eff) & e_cong_eff>1e-7
assert missing(e) if selected & (missing(e_vrs,e_cong_eff) | e_cong_eff<=1e-7)
assert e_se==e_crs/e_vrs if selected & !missing(e_crs,e_vrs) & e_vrs>1e-7
assert missing(e_se) if selected & (missing(e_crs,e_vrs) | e_vrs<=1e-7)
'''
    else:s+='quietly count if selected & missing(e)\nassert Stats[1,6]==r(N)\n'
    s+=f'mata: savebin("results.bin",st_data(.,tokens("{vars}"),"selected"))\n'
    s+='mata: savebin("data.bin",st_data(.,tokens("x1 x2 y1 y2 b1 b2"),"selected"))\n'
    if save:
        saved='id x1 x2 y1 y2 '+vars+(' b1 b2' if model=='sbmundesirable' else '')
        s+=f'''preserve
keep if selected
keep {saved}
tempfile expected
save `expected'
use "saved.dta",clear
cf _all using `expected',all
assert _N==Stats[1,1]
restore
'''
    return s

def lp(data,model,rts,k,infeasible=False):
    x,y,b=data[:,:2],data[:,2:4],data[:,4:];n=len(x);m=s=nb=2
    if model=='congestion':
        c=np.r_[np.zeros(n),1.];ae=np.vstack([np.c_[x.T,-x[k]],np.c_[y.T,np.zeros(s)],np.r_[np.ones(n),0.]])
        be=np.r_[np.zeros(m),y[k],1.];bounds=[(0,None)]*n+[(0,1)]
        witness=np.zeros(n+1);witness[k]=witness[-1]=1
    else:
        nv=n+m+s+nb+1;c=np.zeros(nv);c[-1]=1;c[n:n+m]=-1/(m*x[k]);ae=[];be=[]
        for vals,offset,sign in [(x,n,1),(y,n+m,-1),(b,n+m+s,1)]:
            for j in range(vals.shape[1]):
                row=np.zeros(nv);row[:n]=vals[:,j];row[offset+j]=sign;row[-1]=-vals[k,j];ae.append(row);be.append(0.)
        row=np.zeros(nv);row[-1]=1;row[n+m:n+m+s]=1/((s+nb)*y[k]);row[n+m+s:n+m+s+nb]=1/((s+nb)*b[k]);ae.append(row);be.append(1.)
        if rts=='vrs':row=np.zeros(nv);row[:n]=1;row[-1]=-1;ae.append(row);be.append(0.)
        ae=np.array(ae);be=np.array(be);bounds=(0,None);witness=np.zeros(nv);witness[k]=witness[-1]=1
    assert np.all(data>0) and np.max(abs(ae@witness-be))<1e-12
    if infeasible:ae=np.vstack([ae,np.zeros(len(c))]);be=np.r_[be,1.]
    fit=linprog(c,A_eq=ae,b_eq=be,bounds=bounds,method='highs',options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    if infeasible:assert fit.status==2,fit.message
    else:assert fit.success,fit.message
    return fit

def oracle(d,model,rts):
    data=readbin(d/'data.bin');values=readbin(d/'results.bin');maxgap=0.;unavailable=0
    for k in range(len(data)):
        fit=lp(data,model,rts,k);actual=values[k,3] if model=='congestion' else values[k,0]
        if np.isnan(actual):unavailable+=1
        else:
            gap=abs(actual-fit.fun);maxgap=max(maxgap,gap);assert gap<1e-6,(k,actual,fit.fun,gap)
        if model=='congestion':
            x=data[:,:2];y=data[:,2:4]
            crs=radial_solve(x[k],y[k],x,y,'input','crs')
            vrs=radial_solve(x[k],y[k],x,y,'input','vrs')
            for column,expected in [(0,vrs/fit.fun),(1,crs),(2,vrs),(5,crs/vrs)]:
                actual=values[k,column]
                if np.isfinite(actual):
                    gap=abs(actual-expected);maxgap=max(maxgap,gap);assert gap<1e-6,(k,column,actual,expected,gap)
    return dict(independent_lp_checks=len(data)*(3 if model=='congestion' else 1),max_oracle_gap=maxgap,native_unavailable=unavailable)

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--baseline-ado',type=Path,required=True);ap.add_argument('--only');args=ap.parse_args()
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir();hashes={}
    for n in package_files(REPO):
        raw=(REPO/n).read_bytes();(pkg/Path(n).name).write_bytes(raw);hashes[n]=hashlib.sha256(raw).hexdigest()
    audit=out/'instrumented_ado';instrument(pkg,audit);baseline=args.baseline_ado.resolve()
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,baseline_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in baseline.iterdir() if p.is_file()}),indent=2)+'\n')
    rows=[]
    def run(name,body,source=pkg,check=None):
        if args.only and not re.search(args.only,name):return None
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 160\n'
        for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/key).mkdir();head+=f'sysdir set {key} "{d/key}"\n'
        (d/'test.do').write_text(head+f'adopath ++ "{source}"\n'+INIT+WRITER+body+'\ndi "TEST_PASS"\n');start=time.monotonic();row=dict(test=name,passed=False)
        try:
            p=subprocess.run([args.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=240)
            log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else p.stderr.decode(errors='replace')
            assert p.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2600:]
            if check:row.update(check(d))
            row['passed']=True
        except Exception as e:row['error']=repr(e)
        row['seconds']=round(time.monotonic()-start,3);rows.append(row);(out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(row),flush=True)
        return d
    # Original numerical results, independent LP scores and truthful natural failures.
    for fixture,model,rts,engine,selected in itertools.product(FIXTURES,['congestion','sbmundesirable'],['crs','vrs'],['simplex','ipm','auto'],[False,True]):
        name=f'history_{fixture}_{model}_{rts}_{engine}_{selected}'
        if args.only and not re.search(args.only,name):continue
        sample='if id!=2 in 2/6' if selected else '';condition='id!=2 & inrange(_n,2,6)' if selected else '1'
        body=FIXTURES[fixture]+f'gen byte selected={condition}\n'
        cmd=command(model,'coreA',rts,engine,sample,True)
        variables='e e_crs e_vrs e_cong_eff e_congestion e_se' if model=='congestion' else 'e'
        old=run(name+'_before',body+cmd+f'mata: savebin("results.bin",st_data(.,tokens("{variables}"),"selected"))\n',baseline)
        def compare(d,old=old,model=model,rts=rts):
            metrics=oracle(d,model,rts);a=readbin(d/'results.bin');b=readbin(old/'results.bin')
            finite=np.isfinite(a)
            if model=='sbmundesirable':
                # R08 normalizes columns and applies the requested residual test.
                # Keep sample/failure parity and independently validate each model;
                # changed final bits are bounded by the existing oracle threshold.
                np.testing.assert_array_equal(np.isfinite(a),np.isfinite(b))
                np.testing.assert_allclose(a,b,rtol=0,atol=1e-6,equal_nan=True)
            else:np.testing.assert_array_equal(a[finite],b[finite])
            metrics['max_historical_gap']=float(np.max(abs(a[finite]-b[finite]))) if finite.any() else 0.
            metrics.update(historical_pair=True,exact_numeric_file=(d/'results.bin').read_bytes()==(old/'results.bin').read_bytes(),newly_unavailable_cells=int((np.isnan(a)&np.isfinite(b)).sum()))
            return metrics
        run(name+'_after',body+command(model,'coreB',rts,engine,sample,True)+checks(model,'coreB',True),check=compare)
    # Numerical failure versus a genuinely contradictory test LP (0=1).
    for model,rts in [('congestion','crs'),('sbmundesirable','crs'),('sbmundesirable','vrs')]:
        flag='cong' if model=='congestion' else 'sbmu'
        for route,mode,target,selected in itertools.product(['coreA','coreB','wrapA','wrapB'],[1,2],[2,-1],[False,True]):
            sample='if id!=2 in 2/6' if selected else '';condition='id!=2 & inrange(_n,2,6)' if selected else '1';n=4 if selected else 6;nf=n if target==-1 else 1
            s=DATA+f'gen byte selected={condition}\nscalar failure_mode={mode}\nscalar fault_{flag}={target}\n'
            s+=command(model,route,rts,sample=sample,save=route.startswith('core'))+checks(model,route,route.startswith('core'))
            s+=f'assert Stats[1,6]=={nf}\n'
            if target==-1:s+='assert Stats[1,2]==0 & missing(Stats[1,3],Stats[1,4],Stats[1,5])\n'
            def check(d,model=model,rts=rts,mode=mode,target=target):
                data=readbin(d/'data.bin');targets=range(len(data)) if target==-1 else [target-1]
                if mode==2:
                    for k in targets:lp(data,model,rts,k,True)
                    return dict(independent_infeasible_checks=len(targets))
                return {}
            run(f'fault_{model}_{rts}_{route}_{mode}_{target}_{selected}',s,audit,check)
    # Component failures: preserve remaining scores; count unions rather than sums.
    for route,engine,component,target in itertools.product(['coreA','coreB','wrapA','wrapB'],['simplex','ipm','auto'],['crs','vrs','cong','mixed'],[2,-1]):
        s=DATA+'gen byte selected=1\n'
        if component=='mixed':s+='scalar fault_crs=1\nscalar fault_vrs=2\nscalar fault_cong=3\n'
        else:s+=f'scalar fault_{component}={target}\n'
        s+=command('congestion',route,engine=engine,save=route.startswith('core'))+checks('congestion',route,route.startswith('core'))
        n=6 if target==-1 else 1
        if component=='mixed':s+='assert Stats[1,6]==3 & Stats[1,9]==2 & Stats[1,13]==2\n'
        else:
            s+=f'assert Stats[1,6]=={n}\nassert Stats[1,9]=={0 if component=="crs" else n}\n'
            if component=='crs':s+='assert !missing(e) if selected\n'
            if component=='cong':s+='assert !missing(e_se) if selected\n'
        run(f'component_{route}_{engine}_{component}_{target}',s,audit)
    # Ordinary/extended missing, undefined denominators, and efficient-count boundaries.
    for model,flag,values in [('congestion','cong',['.','.a','.z','0','-1','1e-7','1e-8']),('congestion','vrs',['.','.a','.z','0','1e-7']),('sbmundesirable','sbmu',['.','.a','.z','1','1-5e-8','1+5e-8','1+2e-7'])]:
        for route,value,target in itertools.product(['coreA','wrapB'],values,[2,-1]):
            s=DATA+f'gen byte selected=1\nscalar failure_mode=3\nscalar fault_{flag}={target}\nscalar fault_value={value}\n'
            s+=command(model,route,save=route.startswith('core'))+checks(model,route,route.startswith('core'))
            run(f'value_{model}_{flag}_{route}_{value}_{target}',s,audit)
    # A subsequent TE call must not inherit specialized-model failure counters.
    for model in ['congestion','sbmundesirable']:
        s=DATA+'gen byte selected=1\n'+command(model,'coreA')+'dea x1 x2=y1 y2,dmu(id) engine(simplex) gen(te)\n'
        for v in CFIELDS:s+=f'capture confirm scalar r({v})\nassert _rc!=0\n' if False else f'assert missing(r({v}))\n'
        s+='matrix S=r(stage_status)\nassert rowsof(S)==6\n'
        run('no_stale_'+model,s)
    summary=dict(total=len(rows),passed=sum(x['passed'] for x in rows),historical_pairs=sum(x.get('historical_pair',False) for x in rows),exact_numeric_files=sum(x.get('exact_numeric_file',False) for x in rows),independent_lp_checks=sum(x.get('independent_lp_checks',0) for x in rows),independent_infeasible_checks=sum(x.get('independent_infeasible_checks',0) for x in rows),natural_unavailable=sum(x.get('native_unavailable',0) for x in rows),max_oracle_gap=max((x.get('max_oracle_gap',0) for x in rows),default=0))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
    raise SystemExit(summary['total']!=summary['passed'])

if __name__=='__main__':main()
