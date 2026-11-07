#!/usr/bin/env python3
"""Check certified prices against exact rationals and exercise cycle-gated recovery.

Requires NumPy/SciPy, licensed Stata16+, and the packaged hexadecimal fixture.
Observational instrumentation is confined to a private output copy.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import argparse,hashlib,json,shutil,subprocess
import numpy as np
from test_core_retries import WRITER,readbin
from test_simplex_stability import mata
from test_precision_recovery import decode

REPO=Path(__file__).resolve().parents[1]
def fm(A):return tuple(tuple(F(float(v)) for v in row) for row in np.atleast_2d(A))
def fv(a):return tuple(F(float(v)) for v in np.ravel(a))
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def trans(A):return tuple(zip(*A))
@lru_cache(None)
def solve(A,b):
    n=len(A)
    if any(len(row)!=n for row in A) or len(b)!=n:return None
    T=[list(row)+[rhs] for row,rhs in zip(A,b)]
    for col in range(n):
        p=next((i for i in range(col,n) if T[i][col]),None)
        if p is None:return None
        T[col],T[p]=T[p],T[col];v=T[col][col];T[col]=[x/v for x in T[col]]
        for i in range(n):
            if i==col:continue
            v=T[i][col]
            if v:T[i]=[a-v*b for a,b in zip(T[i],T[col])]
    return tuple(row[-1] for row in T)

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    root=args.output.resolve();root.mkdir(parents=True,exist_ok=False);pkg=root/'ado';shutil.copytree(REPO/'ado',pkg)
    path=REPO/'tests/fixtures/precision_pricing.json';fixture=json.loads(path.read_text());cases=fixture['pricing'];lps=fixture['lp'];assert len(cases)==1588 and len(lps)==18
    source=(pkg/'_dea_simplex.ado').read_text();h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    (root/'manifest.json').write_text(json.dumps(dict(source_sha256={p.name:h(p) for p in pkg.iterdir() if p.is_file()},fixture_sha256=h(path),numpy=np.__version__,stata=args.stata),indent=2)+'\n')
    batches=[]
    def run(name,body,ado=pkg):
        out=root/name;out.mkdir();head='version 16\nclear all\nset more off\nset processors 1\n'
        for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:(out/key).mkdir();head+=f'sysdir set {key} "{out/key}"\n'
        head+=WRITER+f'adopath ++ "{ado}"\nrun "{ado}/_dea_simplex.ado"\n'
        (out/'test.do').write_text(head+body+'\ndi "PRICING_PASS"\nexit,clear\n')
        p=subprocess.run([args.stata,'-b','do','test.do'],cwd=out,capture_output=True,timeout=900)
        (out/'process.json').write_text(json.dumps(dict(returncode=p.returncode,stderr=p.stderr.decode(errors='replace')),indent=2)+'\n')
        log=(out/'test.log').read_text(errors='replace') if (out/'test.log').exists() else ''
        assert p.returncode==0 and 'PRICING_PASS' in log.splitlines(),(name,log[-4000:])
        batches.append(dict(name=name,passed=True));return out
    body='mata:\n'+r'''
void price_probe(string scalar name,real matrix A,real rowvector c,real colvector bb,real scalar tol) {
 struct RSM_Prob scalar P
 struct RSM_Basis scalar B
 real matrix Q,Y
 real rowvector detail
 P.A=A;P.b=J(rows(A),1,0);P.c=c;P.m_con=rows(A);P.n_orig=P.n_var=cols(A)
 B.basic=bb;B.Binv=J(rows(A),rows(A),.);B.n_upd=0
 Q=_rsm_cert_price_bounds(P,B,Y,detail)
 savebin(name+"_bounds.bin",Q);savebin(name+"_pair.bin",Y);savebin(name+"_detail.bin",detail)
 savebin(name+"_enter.bin",(_rsm_cert_enter(P,B,1,tol),_rsm_cert_enter(P,B,0,tol)))
 savebin(name+"_input.bin",A\c)
}
'''
    inputs=[]
    for i,f in enumerate(cases):
        A=decode(f['A'],len(f['c']));c=np.array([float.fromhex(x) for x in f['c']])[None,:];bb=np.array(f['bb'])[:,None];inputs.append((A,c,bb))
        body+=f'price_probe("b{i:04}",{mata(A)},{mata(c)},{mata(bb)},{f["tol"]:.17g})\n'
    body+='assert(_rsm_cert_divup(0,1)==0)\nassert(missing(_rsm_cert_divup(1,0)))\nassert(missing(_rsm_cert_divup(-1,2)))\nassert(missing(_rsm_cert_divup(.,1)))\nassert(missing(_rsm_cert_divup(1,.)))\nsavebin("division.bin",(_rsm_cert_divup(1,3),_rsm_cert_divup(1,7)))\nend\n'
    out=run('bounds',body);rows=[]
    for i,(f,(A,c,bb)) in enumerate(zip(cases,inputs)):
        af=fm(A);cf=fv(c);basic=tuple(int(v)-1 for v in bb.ravel());M=trans(tuple(tuple(row[k] for k in basic) for row in af));y=solve(M,tuple(cf[k] for k in basic)) if len(set(basic))==len(basic) else None
        name=f'b{i:04}';Q=readbin(out/(name+'_bounds.bin'));Y=readbin(out/(name+'_pair.bin'));detail=readbin(out/(name+'_detail.bin')).ravel();enter=readbin(out/(name+'_enter.bin')).ravel().astype(int).tolist()
        assert np.array_equal(readbin(out/(name+'_input.bin')),np.vstack((A,c)))
        available=bool(np.isfinite(Q).all());row=dict(label=f['label'],passed=True,bounds_available=available,bland=enter[0],dantzig=enter[1],previous_pair_counterexample=f['previous_pair_cutoff_counterexample'])
        if y is None:assert not available and enter==[-1,-1]
        else:
            prices=tuple(ci-dot(y,col) for ci,col in zip(cf,trans(af)));cut=F(float(min(1e-8,f['tol'])*min(1,max(abs(c.ravel())))));nb=[k for k in range(len(cf)) if k not in basic];improve=[k for k in nb if prices[k]<-cut]
            if available:
                py=tuple(F(float(hi))+F(float(lo)) for hi,lo in Y);error=max(abs(a-b) for a,b in zip(y,py));assert error<=F(float(detail[1])) and 0<=detail[0]<1
                for k,price in enumerate(prices):assert F(float(Q[k,2]))<=price<=F(float(Q[k,3])),(f['label'],k)
            else:assert enter==[-1,-1]
            for choice in enter:
                if choice==0:assert not improve,(f['label'],'false optimal')
                if choice>0:assert choice-1 in nb and prices[choice-1]<-cut,(f['label'],'false entering')
            if enter[0]>0:assert enter[0]==improve[0]+1
        rows.append(row)
    division=readbin(out/'division.bin').ravel();assert F(float(division[0]))>=F(1,3) and F(float(division[1]))>=F(1,7)
    assert sum(x['bounds_available'] for x in rows)>=1388 and sum(x['bland']>=0 for x in rows)>=1150
    assert sum(x['previous_pair_counterexample'] for x in rows)==12
    (out/'results.json').write_text(json.dumps(rows,indent=2)+'\n')

    # Observe actual calls in a private copy; numerical outputs must match.
    observed=root/'observed_ado';shutil.copytree(pkg,observed)
    a=source.index('real scalar _rsm_cert_enter(');b=source.index('\n}',a)+2;part=source[a:b]
    part=part.replace(' real matrix bounds,Y',' external real scalar _test_price_calls\n real matrix bounds,Y',1).replace(' if(missing(tol)',' _test_price_calls++\n if(missing(tol)',1)
    instrumented=source[:a]+part+source[b:]
    a=instrumented.index('struct RSM_Result scalar _rsm_precision_two_phase(');b=instrumented.index('\n}',a)+2;part=instrumented[a:b]
    part=part.replace('    real colvector second_basis','    external real rowvector _test_retry_modes\n    real colvector second_basis',1).replace('    cycle_seen=0','    _test_retry_modes[certified+1]=_test_retry_modes[certified+1]+1\n    cycle_seen=0',1)
    instrumented=instrumented[:a]+part+instrumented[b:];(observed/'_dea_simplex.ado').write_text(instrumented)
    body='mata:\n_test_price_calls=0\n_test_retry_modes=(0,0)\n'+r'''
void lp_probe(string scalar name,real matrix E,real colvector b,real matrix U,real colvector u,real rowvector c,real scalar tol) {
 struct RSM_Prob scalar P,T
 struct RSM_Result scalar R,S
 external real scalar _test_price_calls
 external real rowvector _test_retry_modes
 real scalar enabled
 real matrix log
 _rsm_std_form(P,U,u,E,b,c)
 for(enabled=0;enabled<=1;enabled++) {
  _test_price_calls=0;_test_retry_modes=(0,0);_rsm_audit_reset(enabled);T=P
  R=_rsm_two_phase(T,tol,"")
  assert(R.recovery[5]+R.recovery[6]==R.iters & R.recovery[7]==R.iters & R.recovery[8]==R.status)
  log=_rsm_audit_read()
  if(enabled) {
   assert(rows(log)==1 & cols(log)==12 & all(log[1,1..10]:==R.recovery))
   assert(S.status==R.status & S.objval==R.objval & S.iters==R.iters & all(S.x:==R.x))
   savebin(name+"_result.bin",(R.status,R.objval,R.iters));savebin(name+"_x.bin",R.x)
   savebin(name+"_record.bin",log);savebin(name+"_counts.bin",(_test_price_calls,_test_retry_modes))
  }
  else {assert(rows(log)==0);S=R;}
 }
 _rsm_audit_reset(0)
 savebin(name+"_eqmat.bin",E);savebin(name+"_eqrhs.bin",b)
 savebin(name+"_ineqmat.bin",U);savebin(name+"_ineqrhs.bin",u);savebin(name+"_cost.bin",c)
}
'''
    lp_inputs=[]
    for i,f in enumerate(lps):
        n=len(f['xstar']);values=[decode(f[key],nc) for key,nc in [('E',n),('b',1),('U',n),('u',1),('c',n)]];E,b,U,u,c=values;lp_inputs.append(values)
        x=tuple(map(F,f['xstar']));y=tuple(map(F,f['ystar']));assert min(x)>=0 and all(np.all(A>0) for A in values)
        assert all(dot(ar,x)==rhs for ar,rhs in zip(fm(E),fv(b))) and all(dot(ar,x)<=rhs for ar,rhs in zip(fm(U),fv(u)))
        assert min(ci-dot(y,col) for ci,col in zip(fv(c),trans(fm(E))))>=0 and dot(fv(c),x)==dot(y,fv(b))==6
        for key,A in zip(['E','b','U','u','c'],values):body+=f'{key}={mata(A)}\n'
        body+=f'lp_probe("p{i:02}",E,b,U,u,c,{f["tol"]:.17g})\n'
    body+='end\n';plain=run('lp_plain',body);obs=run('lp_observed',body,observed);parity=0;lp_rows=[]
    for p in plain.glob('*.bin'):
        if p.name.endswith('_counts.bin'):continue
        assert p.read_bytes()==(obs/p.name).read_bytes(),p.name;parity+=1
    for i,(f,values) in enumerate(zip(lps,lp_inputs)):
        prefix=f'p{i:02}';E,b,U,u,c=values
        for suffix,A in zip(['eqmat','eqrhs','ineqmat','ineqrhs','cost'],values):assert np.array_equal(A,readbin(plain/(prefix+'_'+suffix+'.bin')))
        result=readbin(plain/(prefix+'_result.bin')).ravel();status=int(result[0]);counts=readbin(obs/(prefix+'_counts.bin')).ravel();assert status==f['expected_status'] and result[2]==f['expected_iterations']
        assert counts[2]==int(f['expected_certified_path']) and bool(counts[0]>0)==f['expected_certified_path']
        if status==0:
            xf=readbin(plain/(prefix+'_x.bin')).ravel();assert np.array_equal(xf,np.array([float.fromhex(v) for v in f['expected_x']]))
            x=fv(xf[:E.shape[1]]);assert min(x)>=0 and dot(fv(c),x)==6 and result[1]==6
            assert all(dot(ar,x)==rhs for ar,rhs in zip(fm(E),fv(b))) and all(dot(ar,x)<=rhs for ar,rhs in zip(fm(U),fv(u)))
        lp_rows.append(dict(label=f['label'],passed=True,status=status,iterations=int(result[2]),price_calls=int(counts[0]),precision_paths=counts[1:].astype(int).tolist()))
    (root/'lp_results.json').write_text(json.dumps(lp_rows,indent=2)+'\n')
    summary=dict(all_passed=True,pricing_conditions=len(cases),price_bounds_available=sum(x['bounds_available'] for x in rows),certified_bland_decisions=sum(x['bland']>=0 for x in rows),deferred_bland_decisions=sum(x['bland']<0 for x in rows),prior_pair_cutoff_counterexamples=12,no_false_certified_decisions=True,division_controls=7,positive_LP_conditions=len(lps),LP_executions=4*len(lps),LP_exact_successes=sum(x['status']==0 for x in lp_rows),LP_unresolved_controls=sum(x['status']!=0 for x in lp_rows),instrumented_numeric_files_identical=parity,certified_path_positive_controls=sum(x['expected_certified_path'] for x in lps),certified_path_skipped_controls=sum(not x['expected_certified_path'] for x in lps),single_aggregate_record=True,batches=batches)
    (root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
