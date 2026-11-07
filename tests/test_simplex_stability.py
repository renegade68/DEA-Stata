#!/usr/bin/env python3
"""Basis/pivot regressions and independently certified LPs. Stata 16+, NumPy.

Certificates use exact Fraction arithmetic on the actual supplied doubles.
All package copies, do-files, outputs and failed runs remain in --output.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import time
import numpy as np

REPO=Path(__file__).resolve().parents[1]

def mata(a):
    a=np.asarray(a,dtype=float)
    if not a.size:return f'J({a.shape[0]},{a.shape[1]},0)'
    return '('+'\\'.join(','.join(format(v,'.17g') for v in r) for r in a)+')'

HELPERS=r'''
void stability_helpers() {
 struct RSM_Basis scalar B
 struct RSM_Prob scalar P
 struct RSM_Result scalar R
 real matrix M,MI,old
 real scalar i,t,p,step
 real colvector x,d
 B.basic=(2\1);B.Binv=I(2);B.n_upd=0
 for (i=7;i<=12;i++) {
  t=10^(-i)
  assert(_rsm_leave((0\5e-8),(1\1),B,1,t)==1)
  assert(_rsm_leave((1\1),(1\1),B,1,t)==2)
  x=(1e-9\2e-6);d=(1\1000)
  p=_rsm_leave(x,d,B,0,t);step=x[p]/d[p]
  assert(min(x-d*step)>=-min((1e-8,t))-1e-20)
  x=(-.75*min((1e-8,t))\.5*min((1e-8,t))*1000);d=(1\1000)
  p=_rsm_leave(x,d,B,0,t);step=max((0,x[p]))/d[p]
  assert(min(x-d*step)>=-min((1e-8,t))-1e-20)
 }
 M=(1,1\2,2);assert(hasmissing(_rsm_inv(M)))
 M=(1,1\1,1+2^(-30));MI=_rsm_inv(M)
 assert(!hasmissing(MI));assert(max(abs(M*MI-I(2)))<=1e-8)
 M=(2^(-40),0\0,2^40);MI=_rsm_inv(M)
 assert(!hasmissing(MI));assert(max(abs(M*MI-I(2)))==0)
 M=(1,1e-9\0,1);MI=_rsm_inv(M)
 assert(max(abs(M*MI-I(2)))<1e-14)
 P.A=(1,1\2,2);P.m_con=2;B.basic=(1\2);B.Binv=I(2)
 assert(!_rsm_refac(B,P));assert(hasmissing(B.Binv))
 B.Binv=I(2);B.n_upd=0;old=B.Binv
 assert(!_rsm_pivot(B,(1e-14\1),1));assert(all(B.Binv:==old));assert(B.n_upd==0)
 // Tiny nonzero eta coefficients must not be discarded.
 assert(_rsm_pivot(B,(1\5e-13),1));assert(B.Binv[2,1]==-5e-13)
 // An invalid initial basis does not establish infeasibility. Exact duplicate
 // recovery may instead certify this feasible LP in its original coordinates.
 _rsm_std_form(P,J(0,2,0),J(0,1,0),(1,1\2,2),(1\2),(1,1))
 P.init_basic=(1\2);R=_rsm_two_phase(P,1e-10,"")
 assert(R.status==0 & R.objval==1 & all(R.x:==(1\0)))
}
void stability_case(string scalar path,real matrix AE,real colvector BE,
 real matrix AU,real colvector BU,real rowvector C,real scalar tol) {
 struct RSM_Result scalar R
 real scalar f,i
 real rowvector x,diag
 R=_rsm_solve_lp(AU,BU,AE,BE,C,tol)
 x=(rows(R.x)>=cols(C) ? R.x[1..cols(C)]' : J(1,cols(C),.))
 diag=_dea_lpcheck(x,AE,BE,AU,BU,J(1,cols(C),0),J(1,cols(C),.))
 f=fopen(path,"w")
 fput(f,sprintf("%g,%26.17e,%g",R.status,R.objval,R.iters))
 for(i=1;i<=cols(x);i++) fput(f,strofreal(x[i],"%26.17e"))
 fput(f,R.msg);fclose(f)
}
'''

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata',required=True);ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--source',type=Path,help='Optional pre-change package for comparison')
    args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    pkg=out/'ado';shutil.copytree(args.source or REPO/'ado',pkg)
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in pkg.iterdir() if p.is_file()}
    (out/'manifest.json').write_text(json.dumps(hashes,indent=2)+'\n')
    cases=[]
    def add(name,ae,be,c,au=None,bu=None,**kw):
        ae=np.asarray(ae,dtype=float);c=np.asarray(c,dtype=float)
        cases.append(dict(name=name,ae=ae,be=np.asarray(be,dtype=float),c=c,
            au=np.zeros((0,len(c))) if au is None else np.asarray(au,dtype=float),
            bu=np.zeros(0) if bu is None else np.asarray(bu,dtype=float),**kw))
    for tol in [1e-7,1e-10,1e-12]:
        add('duplicate_'+str(tol),[[1,1],[2,2]],[1,2],[1,2],objective=1,tol=tol)
        add('zero_rows_'+str(tol),[[0,0],[1,1],[0,0]],[0,1,0],[1,2],objective=1,tol=tol)
        add('redundant_unbounded_'+str(tol),[[0,0],[0,0]],[0,0],[-1,0],status=2,tol=tol)
        add('inconsistent_'+str(tol),[[1,1],[2,2]],[1,3],[1,2],status=1,tol=tol)
        # Beale cycling example: exact optimum -1 at (1,0,1,0).
        add('beale_'+str(tol),np.zeros((0,4)),[],[-10,57,9,24],
            [[.5,-5.5,-2.5,9],[.5,-1.5,-.5,1],[1,0,0,0]],[0,0,1],objective=-1,tol=tol)
        for exp in [-40,0,40]:
            k=2.**exp
            add('tiny_pivot_%s_%s'%(tol,exp),np.zeros((0,1)),[],[-1],[[k]],[k],objective=-1,tol=tol)
        add('objective_units_'+str(tol),np.zeros((0,1)),[],[-2.**-40],[[2.**-40]],[1],objective=-1,tol=tol)
    rng=np.random.default_rng(19178);certificates=[]
    # Known primal/dual optima; duplicate rows, degenerate bases, row and
    # variable unit changes. Dyadic transformations keep certificates exact.
    for seed in range(40):
        m=int(rng.integers(2,6));n=m+int(rng.integers(2,6))
        A=np.c_[np.eye(m),rng.integers(-4,5,size=(m,n-m))]
        x=np.r_[rng.integers(0,5,size=m),np.zeros(n-m)]
        y=rng.integers(-3,4,size=m);s=np.r_[np.zeros(m),rng.integers(0,5,size=n-m)]
        b=A@x;c=y@A+s
        if seed%2==0:A=np.r_[A,A[:1],2*A[:1]];b=np.r_[b,b[:1],2*b[:1]];y=np.r_[y,0,0]
        for scale in ['base','row','column','combined']:
            row=np.ones(len(b));col=np.ones(n)
            if scale in ['row','combined']:row=2.**rng.integers(-35,36,size=len(b))
            if scale in ['column','combined']:col=2.**rng.integers(-20,21,size=n)
            aa=A*row[:,None]*col;bb=b*row;cc=c*col;xx=x/col;yy=y/row
            af=[[F(v) for v in r] for r in aa];bf=list(map(F,bb));cf=list(map(F,cc));xf=list(map(F,xx));yf=list(map(F,yy))
            assert all(sum(v*z for v,z in zip(r,xf))==rhs for r,rhs in zip(af,bf))
            assert min(xf)>=0
            reduced=[cf[j]-sum(yf[i]*af[i][j] for i in range(len(af))) for j in range(n)]
            assert min(reduced)>=0 and sum(a*z for a,z in zip(cf,xf))==sum(a*z for a,z in zip(bf,yf))
            obj=sum(a*z for a,z in zip(cf,xf));name=f'cert_{seed}_{scale}'
            certificates.append(dict(name=name,objective=str(obj),x=list(map(str,xf)),dual=list(map(str,yf)),reduced=list(map(str,reduced))))
            add(name,aa,bb,cc,objective=float(obj),tol=1e-10)
    (out/'certificates.json').write_text(json.dumps(certificates,indent=2)+'\n')
    head='version 16\nclear all\nset more off\nset linesize 160\n'
    for d in ['PLUS','PERSONAL','SITE','OLDPLACE']:(out/d).mkdir();head+=f'sysdir set {d} "{out/d}"\n'
    helpers=HELPERS[HELPERS.index('void stability_case'): ] if args.source else HELPERS
    head+=f'adopath ++ "{pkg}"\nrun "{pkg/"_dea_simplex.ado"}"\nmata:\n'+helpers
    # Baseline comparison keeps the new helper assertions separate.
    if not args.source:head+='stability_helpers()\n'
    for k in cases:
        head+=f'stability_case("{out/(k["name"]+".txt")}",{mata(k["ae"])},{mata(k["be"].reshape(-1,1))},{mata(k["au"])},{mata(k["bu"].reshape(-1,1))},{mata(k["c"].reshape(1,-1))},{k["tol"]})\n'
    head+='end\ndi "STABILITY_COMPLETE"\n';(out/'test.do').write_text(head)
    start=time.monotonic();p=subprocess.run([args.stata,'-b','do','test.do'],cwd=out,capture_output=True,timeout=600)
    (out/'process_stderr.txt').write_bytes(p.stderr)
    log=(out/'test.log').read_text(errors='replace') if (out/'test.log').exists() else ''
    completed=p.returncode==0 and 'STABILITY_COMPLETE' in log.splitlines()
    results=[]
    for k in cases:
        r=dict(name=k['name'],passed=False);f=out/(k['name']+'.txt')
        if f.exists():
            lines=f.read_text().splitlines();st,obj,it=lines[0].split(',')
            st=int(st);x=np.array([float(v) if v.strip()!='.' else np.nan for v in lines[1:1+len(k['c'])]])
            r.update(status=st,iterations=int(it),message=lines[-1])
            if 'objective' in k:
                if st==0:
                    obj=float(obj);gap=abs(obj-k['objective'])/(1+abs(k['objective']))
                    e=k['ae']@x-k['be'];ineq=np.maximum(k['au']@x-k['bu'],0)
                    den=np.maximum(abs(k['be']),abs(k['ae'])@np.maximum(1,abs(x)))
                    u=np.maximum(abs(k['bu']),abs(k['au'])@np.maximum(1,abs(x)))
                    er=np.divide(abs(e),den,out=np.zeros_like(e),where=den!=0).max(initial=0)
                    ir=np.divide(ineq,u,out=np.zeros_like(ineq),where=u!=0).max(initial=0)
                    br=np.maximum(-x,0).max(initial=0)
                    r.update(objective=obj,relative_objective_gap=gap,max_scaled_residual=float(max(er,ir,br)))
                    r['passed']=bool(gap<=1e-8 and max(er,ir,br)<=k['tol']+1e-15)
            else:r['passed']=st==k['status']
        results.append(r)
    summary=dict(completed=completed,cases=len(results),passed=sum(r['passed'] for r in results),seconds=time.monotonic()-start,failed=[r for r in results if not r['passed']])
    (out/'results.json').write_text(json.dumps(results,indent=2)+'\n');(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2));
    if not completed:print(log[-3500:])
    assert completed and summary['passed']==len(results)

if __name__=='__main__':main()
