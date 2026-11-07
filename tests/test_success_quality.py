#!/usr/bin/env python3
"""Successful-basis quality regression: exact original coordinates and conservative bounds.

Self-contained packaged fixtures. Requires NumPy/SciPy and licensed Stata16+.
Private off/strict policies and deliberate perturbations are test controls only.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import argparse,hashlib,itertools,json,shutil,subprocess,time
import numpy as np
from test_core_retries import WRITER,readbin
from test_simplex_stability import mata
from test_precision_recovery import decode
REPO=Path(__file__).resolve().parents[1]
j=lambda p:json.loads(p.read_text())
h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
fm=lambda A:tuple(tuple(F(float(v)) for v in row) for row in np.atleast_2d(A))
fv=lambda a:tuple(F(float(v)) for v in np.ravel(a))
dot=lambda a,b:sum((x*y for x,y in zip(a,b)),F(0))
trans=lambda A:tuple(zip(*A))
r=None
args=None
runs=[]
def run(name,body):
    out=r/'checks'/name;out.mkdir()
    head='version 16\nclear all\nset more off\nset processors 1\nset linesize 180\n'
    for k in ['PLUS','PERSONAL','SITE','OLDPLACE']:
        (out/k).mkdir();head+=f'sysdir set {k} "{out/k}"\n'
    head+=WRITER+f'adopath ++ "{r}/ado"\n'
    (out/'test.do').write_text(head+body+'\ndi "SUCCESS_QUALITY_PASS"\nexit,clear\n')
    start=time.monotonic()
    p=subprocess.run([args.stata,'-b','do','test.do'],cwd=out,capture_output=True,timeout=900)
    log=(out/'test.log').read_text(errors='replace') if (out/'test.log').exists() else ''
    meta=dict(name=name,returncode=p.returncode,seconds=time.monotonic()-start,stderr=p.stderr.decode(errors='replace'))
    (out/'process.json').write_text(json.dumps(meta,indent=2)+'\n')
    assert p.returncode==0 and 'SUCCESS_QUALITY_PASS' in log.splitlines(),(name,log[-5000:])
    runs.append(meta);print(name,'passed',flush=True);return out

def check_accepted():
    fixture=j(REPO/'tests/fixtures/success_quality.json')
    script=(REPO/'tests/fixtures/success_quality.do').read_text().replace('@ADO@',str(r/'ado'))
    rows=[]
    for policy in [0,1,2]:
        out=run('accepted_'+str(policy),script.replace('@POLICY@',str(policy)))
        corrected=preserved=0
        for z in fixture:
            id=z['id'];tag='off' if policy==0 else 'retain'
            result=readbin(out/f'{id}_result.bin').ravel();x=readbin(out/f'{id}_x.bin').ravel()
            assert [float(v).hex() for v in result]==z[tag+'_result'],(id,policy,'result')
            assert [float(v).hex() for v in x]==z[tag+'_x'],(id,policy,'x')
            assert result[0]==0 and result[2]==float.fromhex(z['off_result'][2])
            if policy:
                meta=readbin(out/f'{id}_success.bin');assert meta.shape==(1,14) and meta[0,0]==z['expected_outcome']
                if z['expected_outcome']==2:
                    if z['kind']=='public':
                        data=readbin(out/f'{id}_input.bin');E=fm(data[1:,:-1]);b=fv(data[1:,-1]);c=fv(data[0,:-1])
                    else:
                        E=fm(readbin(out/f'{id}_eqmat.bin'));b=fv(readbin(out/f'{id}_eqrhs.bin'));c=fv(readbin(out/f'{id}_cost.bin'))
                    xx=fv(x)
                    assert min(xx)>=0 and all(dot(ar,xx)==bb for ar,bb in zip(E,b))
                    yy=tuple(map(F,z['problem']['ystar']))
                    assert all(cc-dot(yy,col)>=0 for cc,col in zip(c,trans(E)))
                    assert dot(c,xx)==dot(yy,b)==6;corrected+=1
                else:
                    assert z['retain_x']==z['off_x'] and z['retain_result']==z['off_result'];preserved+=1
        assert len(fixture)==139
        if policy:
            assert (corrected,preserved)==(137,2)
            assert readbin(out/'boundary_result.bin')[0,1]==6 and readbin(out/'boundary_success.bin')[0,0]==2
        rows.append(dict(policy=policy,conditions=139,exact_corrected=corrected,certified_preserved=preserved))
    (r/'ACCEPTED.json').write_text(json.dumps(dict(passed=True,rows=rows,boundary=1),indent=2)+'\n')

def check_controls():
    body=f'run "{r}/ado/_dea_simplex.ado"\n'
    body+='mata:\n'+r'''
    void control(string scalar name,real matrix E,real colvector b,real rowvector c) {
     struct RSM_Prob scalar P
     struct RSM_Result scalar R
     real scalar policy,record,f
     string scalar id
     for(policy=0;policy<=2;policy++) for(record=0;record<=1;record++) {
      _rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
      _rsm_success_reset(policy,record);_rsm_audit_reset(record)
      R=_rsm_two_phase(P,1e-10,"")
      id=name+"_p"+strofreal(policy)+"_r"+strofreal(record)
      savebin(id+"_result.bin",(R.status,R.objval,R.iters));savebin(id+"_x.bin",R.x)
      savebin(id+"_success.bin",_rsm_success_read());savebin(id+"_counts.bin",_rsm_success_counts_read())
      savebin(id+"_audit.bin",_rsm_audit_read())
      f=fopen(id+"_message.txt","w");fput(f,R.msg);fclose(f)
     }
    }
    control("zero_basic",I(2),(0\1),(1,1))
    control("alternative",(1,1),(1),(1,1))
    control("duplicate",(1,1\2,2),(1\2),(1,2))
    control("negative_duplicate",(1,1\-2,-2),(1\-2),(1,2))
    control("zero_rows",(0,0\1,1\0,0),(0\1\0),(1,2))
    control("sum_row",(1,0\0,1\1,1),(1\2\3),(1,1))
    control("limit4096",J(1,4094,1),1,J(1,4094,1))
    control("limit4097",J(1,4095,1),1,J(1,4095,1))
    // Fault/guard controls on the screen, separate from solver-selected bases.
    P=RSM_Prob();_rsm_std_form(P,J(0,2,0),J(0,1,0),I(2),(1\1),(1,1))
    R=RSM_Result();R.status=0;R.x=(1\1);R.objval=2;R.iters=2
    savebin("invalid_basis.bin",_rsm_success_check(P,R,(1\1),1e-10))
    P.A[.,2]=P.A[.,1]
    savebin("singular_basis.bin",_rsm_success_check(P,R,(1\2),1e-10))
    i=(1::800);E=(J(800,1,1),i,i:^2)
    _rsm_std_form(P,J(0,3,0),J(0,1,0),E,J(800,1,1),(1,1,1))
    R.x=(1\0\0);R.objval=1
    savebin("map_budget.bin",_rsm_success_check(P,R,(1\2),1e-10))
    end
    '''
    out=run('controls',body)
    rows=[]
    for name in ['zero_basic','alternative','duplicate','negative_duplicate','zero_rows','sum_row','limit4096','limit4097']:
     for policy in [0,1,2]:
      p=f'{name}_p{policy}';a=readbin(out/f'{p}_r1_result.bin').ravel();meta=readbin(out/f'{p}_r1_success.bin');off=readbin(out/f'{p}_r0_success.bin')
      assert off.shape[0]==0
      for suffix in ['result.bin','x.bin','counts.bin','message.txt']:assert h(out/f'{p}_r0_{suffix}')==h(out/f'{p}_r1_{suffix}'),(p,suffix)
      if policy!=2 or name not in ['limit4097']:assert a[0]==0,(p,a)
      else:assert a[0]==3,(p,a)
      if a[0]==0:
       for suffix in ['result.bin','x.bin','message.txt']:assert h(out/f'{p}_r1_{suffix}')==h(out/f'{name}_p0_r1_{suffix}'),(p,suffix)
      rows.append(dict(name=name,policy=policy,status=int(a[0]),iterations=int(a[2]),outcomes=meta[:,0].tolist(),reasons=meta[:,1].tolist(),dual_solves=int(readbin(out/f'{p}_r1_counts.bin')[0,5]),primal_solves=int(readbin(out/f'{p}_r1_counts.bin')[0,6])))
    assert readbin(out/'zero_basic_p1_r1_success.bin')[0,3]==2 and sum(readbin(out/'zero_basic_p1_r1_x.bin').ravel()>0)==1
    assert readbin(out/'invalid_basis.bin')[0,1]==2
    assert readbin(out/'singular_basis.bin')[0,1]==4
    assert readbin(out/'map_budget.bin')[0,1]==3 and readbin(out/'map_budget.bin')[0,5]<=1000000
    (r/'CONTROLS.json').write_text(json.dumps(dict(passed=True,solver_conditions=8,executions=48,guard_controls=3,recording_invariant=True,zero_basic_proven=True,rows=rows),indent=2)+'\n')
    print('48 solver controls + 3 direct guard controls passed',flush=True)

def check_arithmetic():
    rng=np.random.default_rng(61104031);cases=[]
    for i in range(1200):
     m=int(rng.integers(1,10));n=int(rng.integers(m,17))
     A=rng.integers(-9,10,(m,n)).astype(float)*np.exp2(rng.integers(-20,21,(m,n)))
     L=rng.integers(-5,6,n).astype(float)*np.exp2(rng.integers(-12,13,n))
     x=L+rng.integers(0,5,n).astype(float)*np.exp2(rng.integers(-12,13,n))
     y=rng.integers(-6,7,m).astype(float)*np.exp2(rng.integers(-10,11,m))
     b=A@x;c=y@A
     c+=np.where(x==L,rng.integers(0,5,n),0)
     if i%3==0:x[0]=np.nextafter(x[0],np.inf)
     if i%3==1:b[0]+=max(1,abs(b[0]))*float(rng.choice([1e-7,1e-9,1e-11]))
     cases.append(dict(A=A,b=b.reshape(-1,1),c=c.reshape(1,-1),x=x.reshape(-1,1),y=y.reshape(1,-1),L=L.reshape(-1,1),tol=float(rng.choice([1e-7,1e-10,1e-12])),kind='mixed_affine'))
    for tol in [1e-7,1e-8,1e-10,1e-12]:
     cutoff=min(1e-8,tol)
     for factor in [0,.5,.99,.999999,1,1.000001,1.01,2]:
      for sign in [1,-1]:
       x=1+factor*2*cutoff/(1-cutoff)
       cases.append(dict(A=np.array([[1.]]),b=np.array([[1.]]),c=np.array([[float(sign)]]),x=np.array([[x]]),y=np.array([[float(sign)]]),L=np.array([[0.]]),tol=tol,kind='cutoff_boundary'))
    for power in [-401,-400,-399,-150,-101,-100,-99,0,99,100,101,150,399,400,401]:
     scale=2.**power
     for sign in [1,-1]:
      cases.append(dict(A=np.array([[scale]]),b=np.array([[scale]]),c=np.array([[float(sign)]]),x=np.array([[1.]]),y=np.array([[sign/scale]]),L=np.array([[-2.]]),tol=1e-10,kind='domain_boundary'))
    cases.append(dict(A=np.eye(2),b=np.array([[2.**30],[2.**30]]),c=np.array([[1.,-1.]]),x=np.array([[2.**30],[2.**30]]),y=np.array([[1.,-1.]]),L=np.array([[-2.**30],[-2.**30]]),tol=1e-10,kind='objective_cancellation'))
    for power in [-400,400]:
     scale=2.**power
     cases.append(dict(A=np.array([[scale]]),b=np.array([[0.]]),c=np.array([[0.]]),x=np.array([[scale]]),y=np.array([[scale]]),L=np.array([[0.]]),tol=1e-10,kind='triple_product_guard'))
    body=f'run "{r}/ado/_dea_simplex.ado"\nmata:\n'
    body+='''void measure(string scalar id,real matrix A,real colvector b,real rowvector c,
    real colvector x,real rowvector y,real colvector L,real scalar tol) {
     struct RSM_Prob scalar P
     real rowvector fast,slow,standard
     _rsm_std_form(P,J(0,cols(A),0),J(0,1,0),A,b,c)
     // Direct arithmetic tests use submitted signed rows, not std_form sign flips.
     P.A=A;P.b=b;P.c=c;P.n_var=P.n_orig=cols(A);P.m_con=rows(A)
     fast=_rsm_quality_fast(P,x,y,L,tol)
     slow=_rsm_affine_quality(P,x,y,J(1,rows(A),0),L,tol)
     savebin(id+"_fast.bin",fast);savebin(id+"_slow.bin",slow)
    }
    '''
    for i,z in enumerate(cases):
     body+=''.join(k+'='+mata(z[k])+'\n' for k in ['A','b','c','x','y','L'])+f'measure("q{i:04}",A,b,c,x,y,L,{z["tol"]:.17g})\n'
    body+='end\n';out=run('arithmetic',body);results=[]
    for i,z in enumerate(cases):
     A=fm(z['A']);b=fv(z['b']);c=fv(z['c']);x=fv(z['x']);y=fv(z['y']);L=fv(z['L'])
     p=tuple(bb-dot(ar,x) for ar,bb in zip(A,b));rc=tuple(cc-dot(y,col) for cc,col in zip(c,trans(A)))
     delta=tuple(xx-ll for xx,ll in zip(x,L));O=dot(c,x);D=dot(b,y)+dot(L,rc);gap=O-D
     assert gap==dot(rc,delta)-dot(y,p)
     pu=sum(abs(yy)*abs(pp) for yy,pp in zip(y,p));cu=sum(abs(rr)*abs(dd) for rr,dd in zip(rc,delta));du=max(F(0),-min(rc))
     cutoff=F(float(min(1e-8,z['tol'])));price=F(float(min(1e-8,z['tol'])));cs=min(F(1),max(abs(cc) for cc in c))
     exact_accept=abs(gap)<=cutoff*(1+abs(O)) and pu<=cutoff*(1+abs(O)) and cu<=cutoff*(1+abs(O))
     record=dict(id=f'q{i:04}',kind=z['kind'],objective=str(O),gap=str(gap),pu=str(pu),cu=str(cu),du=str(du))
     for method in ['fast','slow']:
      q=readbin(out/f'q{i:04}_{method}.bin').ravel();available=bool(np.isfinite(q).all());passed=available and q[0]==1
      if available:
       assert F(float(q[9]))>=abs(gap) and F(float(q[10]))>=pu and F(float(q[11]))>=cu and F(float(q[6]))>=du,(i,method)
      if passed:assert exact_accept and du<=price*cs,(i,method,'false acceptance')
      record[method]=dict(available=available,accepted=bool(passed))
     if record['fast']['accepted']:assert record['slow']['accepted'],(i,'fast passed but EFT did not')
     results.append(record)
    fixture=[{k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in z.items()} for z in cases]
    (r/'ARITHMETIC_FIXTURE.json').write_text(json.dumps(fixture,indent=2)+'\n')
    summary=dict(passed=True,conditions=len(cases),kinds=dict(Counter(z['kind'] for z in cases)),fast_available=sum(z['fast']['available'] for z in results),fast_accepted=sum(z['fast']['accepted'] for z in results),slow_accepted=sum(z['slow']['accepted'] for z in results),all_available_bounds_enclose_exact=True,no_false_acceptance=True,rows=results)
    (r/'ARITHMETIC.json').write_text(json.dumps(summary,indent=2)+'\n');print({k:v for k,v in summary.items() if k!='rows'},flush=True)

def check_public_units():
    cases=[];body=f'run "{r}/ado/_dea_simplex.ado"\n'
    for k,shift,scale,sense,policy in itertools.product([5,25],[(0,0,0),(-6,0,0),(2,-4,1),(-2,-1,-3)],['none','row'],['min','max'],[0,1,2]):
     d=2.**(-k);E=np.array([[1,2,4],[1,2+d,4-d],[1+d,2,4-d]],dtype=float);c=np.array([[1.,2.,3.]])
     L=np.array(shift,dtype=float).reshape(1,-1);xstar=np.ones(3)+L.ravel()
     b=E@xstar;ub=L+4
     assert all(dot(ar,fv(xstar))==bb for ar,bb in zip(fm(E),fv(b)))
     if sense=='max':c=-c
     objective=dot(fv(c),fv(xstar));id=f'q{len(cases):03}'
     case=dict(id=id,k=k,lb=L.tolist(),ub=ub.tolist(),E=E.tolist(),b=b.tolist(),c=c.tolist(),xstar=xstar.tolist(),objective=str(objective),scale=scale,sense=sense,policy=policy);cases.append(case)
     body+='mata: _rsm_success_reset('+str(policy)+',1)\nmata: st_matrix("E",'+mata(E)+');st_matrix("B",'+mata(b.reshape(-1,1))+');st_matrix("C",'+mata(c)+');st_matrix("LB",'+mata(L)+');st_matrix("UB",'+mata(ub)+')\n'
     body+=f'_dea_lpcontract, c(C) aeq(E) beq(B) lb(LB) ub(UB) sense({sense}) engine(simplex) scaling({scale}) tolerance(1e-10)\n'
     body+=f'mata: savebin("{id}_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")));savebin("{id}_x.bin",st_matrix("r(x)"));savebin("{id}_success.bin",_rsm_success_read());savebin("{id}_fast.bin",_rsm_success_fastread());savebin("{id}_attempts.bin",st_matrix("r(attempts)"))\n'
    out=run('public_units',body);rows=[]
    for z in cases:
     id=z['id'];res=readbin(out/f'{id}_result.bin').ravel();x=readbin(out/f'{id}_x.bin').ravel();meta=readbin(out/f'{id}_success.bin');exact=all(dot(ar,fv(x))==bb for ar,bb in zip(fm(z['E']),fv(z['b'])))
     if z['policy']>0:
      assert res[0]==0 and meta.shape==(1,14) and meta[0,0] in [1,2],(id,res,meta)
      if meta[0,0]==2:assert exact and dot(fv(z['c']),fv(x))==F(z['objective']),id
     else:assert res[0]==0,id
     assert min(x-np.ravel(z['lb']))>=0 and min(np.ravel(z['ub'])-x)>=0
     rows.append(dict(id=id,policy=z['policy'],exact=exact,objective_error=float(abs(dot(fv(z['c']),fv(x))-F(z['objective']))),outcome=int(meta[0,0]) if len(meta) else None,iterations=int(res[2])))
    for i in range(0,len(rows),3):
     assert rows[i+1]['outcome']==rows[i+2]['outcome'] and rows[i]['iterations']==rows[i+1]['iterations']==rows[i+2]['iterations']
     for suffix in ['result.bin','x.bin']:assert h(out/f'{rows[i+1]["id"]}_{suffix}')==h(out/f'{rows[i+2]["id"]}_{suffix}')
     if rows[i+1]['outcome']==1:
      for suffix in ['result.bin','x.bin']:assert h(out/f'{rows[i]["id"]}_{suffix}')==h(out/f'{rows[i+1]["id"]}_{suffix}')
    (r/'PUBLIC_UNITS_FIXTURE.json').write_text(json.dumps(cases,indent=2)+'\n')
    summary=dict(passed=True,commands=len(rows),base_conditions=len(rows)//3,zero_objective_conditions=sum(F(z['objective'])==0 for z in cases)//3,refined_each_policy=sum(z['policy']==1 and z['outcome']==2 for z in rows),preserved_each_policy=sum(z['policy']==1 and z['outcome']==1 for z in rows),rows=rows)
    (r/'PUBLIC_UNITS.json').write_text(json.dumps(summary,indent=2)+'\n');print({k:v for k,v in summary.items() if k!='rows'},flush=True)

def check_affine_refinement():
    cases=[x for x in j(r/'PUBLIC_UNITS_FIXTURE.json') if x['policy']==1]
    body=f'run "{r}/ado/_dea_simplex.ado"\nrun "{r}/ado/_dea_lpcontract.ado"\nmata:\n'
    body+='''void exercise(string scalar id,real matrix E,real colvector b,real rowvector c,
     real rowvector L,real rowvector ub,string scalar sense) {
     struct RSM_Result scalar R
     real rowvector actual,cost
     real matrix U
     real colvector u
     cost=(sense=="max" ? -c : c)
     U=I(cols(c));u=(ub-L)'
     R=_rsm_solve_lp(U,u,E,b-E*L',cost,1e-10,0)
     assert(R.status==0)
     savebin(id+"_before.bin",(R.status,R.objval,R.iters))
     R.x[1]=R.x[1]+1e-6
     _rsm_success_reset(2,1)
     _rsm_public_success(R,c,sense,E,b,J(0,cols(c),0),J(0,1,0),L,ub,1e-10,actual)
     savebin(id+"_result.bin",(R.status,c*actual',R.iters));savebin(id+"_x.bin",actual)
     savebin(id+"_success.bin",_rsm_success_read())
    }
    '''
    for z in cases:
     E=np.array(z['E']);b=np.array(z['b']).reshape(-1,1);c=np.array(z['c']);L=np.array(z['lb']);ub=np.array(z['ub'])
     body+=''.join(k+'='+mata(a)+'\n' for k,a in [('E',E),('b',b),('c',c),('L',L),('ub',ub)])+f'exercise("{z["id"]}",E,b,c,L,ub,"{z["sense"]}")\n'
    body+='end\n';out=run('affine_refinement',body);rows=[]
    for z in cases:
     id=z['id'];res=readbin(out/f'{id}_result.bin').ravel();x=fv(readbin(out/f'{id}_x.bin'));meta=readbin(out/f'{id}_success.bin')
     assert res[0]==0 and meta[0,0]==2,(id,res,meta)
     assert all(dot(ar,x)==bb for ar,bb in zip(fm(z['E']),fv(z['b']))) and dot(fv(z['c']),x)==F(z['objective']),id
     assert res[2]==readbin(out/f'{id}_before.bin').ravel()[2]
     rows.append(dict(id=id,lower=z['lb'],objective=z['objective'],exact=True,iterations_preserved=True))
    (r/'AFFINE_REFINEMENT.json').write_text(json.dumps(dict(passed=True,conditions=len(rows),actual_bases=True,controlled_primal_perturbation=1e-6,all_exact=True,rows=rows),indent=2)+'\n');print('32 actual-basis affine refinement controls passed',flush=True)

def check_bounds():
    body=f'run "{r}/ado/_dea_simplex.ado"\n';ids=[]
    for policy,scale,sense,lower in itertools.product([0,1,2],['none','row'],['min','max'],[-2,1]):
     id=f'p{policy}_{scale}_{sense}_{lower}';ids.append(id)
     body+=f'mata: _rsm_success_reset({policy},1)\nmatrix C=(1,2)\nmatrix E=(1,1)\nmatrix B=(3)\nmatrix LB=({lower},{lower})\nmatrix UB=(4,4)\n'
     body+=f'_dea_lpcontract, c(C) aeq(E) beq(B) lb(LB) ub(UB) sense({sense}) engine(simplex) scaling({scale}) tolerance(1e-10)\n'
     body+=f'mata: savebin("{id}_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")));savebin("{id}_x.bin",st_matrix("r(x)"));savebin("{id}_success.bin",_rsm_success_read())\n'
    for policy,kind in itertools.product([0,1,2],['mixed','integer','limited']):
     id=f'p{policy}_{kind}';ids.append(id);integer='1' if kind!='integer' else '1 2';limit='maxnodes(1)' if kind=='limited' else ''
     body+=f'mata: _rsm_success_reset({policy},1)\nmatrix C=(1,1)\nmatrix E=(1,1)\nmatrix B=(1.5)\nmatrix LB=(-1,-1)\nmatrix UB=(3,3)\n'
     body+=f'milp_bb, c(C) aeq(E) beq(B) lb(LB) ub(UB) sense(min) integer({integer}) lpengine(simplex) lptol(1e-10) {limit}\n'
     body+=f'mata: savebin("{id}_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(nodes)")));savebin("{id}_x.bin",st_matrix("r(x)"));savebin("{id}_success.bin",_rsm_success_read())\n'
    out=run('bounds',body);rows=[]
    for id in ids:
     result=readbin(out/f'{id}_result.bin').ravel();x=readbin(out/f'{id}_x.bin').ravel();meta=readbin(out/f'{id}_success.bin')
     if '_min_' in id or '_max_' in id:
      lower=int(id.rsplit('_',1)[1]);expect=(2 if '_min_' in id else 7) if lower==-2 else (4 if '_min_' in id else 5)
      assert result[0]==0 and result[1]==expect and sum(x)==3 and min(x)>=lower and max(x)<=4,(id,result,x)
     else:
      kind=id.split('_')[1];expect={'mixed':0,'integer':1,'limited':4}[kind];assert result[0]==expect,(id,result)
      if kind=='mixed':assert result[1]==1.5 and sum(x)==1.5 and x[0]==np.floor(x[0])
     if id[1]!='0':
      original='p0'+id[2:]
      for suffix in ['result.bin','x.bin']:assert h(out/f'{id}_{suffix}')==h(out/f'{original}_{suffix}'),id
     rows.append(dict(id=id,status=int(result[0]),counts=meta[:,0].tolist()))
    (r/'BOUNDS.json').write_text(json.dumps(dict(passed=True,LP_commands=24,MILP_commands=9,all_existing_results_preserved=True,rows=rows),indent=2)+'\n');print('24 bounded LP + 9 MILP commands passed',flush=True)

def check_adapter_cache():
    body='clear\ninput double x str2 rel double rhs\n1 "" .\n1 ">=" 2\nend\nmatrix C=(1)\nmatrix E=(1)\nmatrix B=(2)\n'
    helpers=['_rsm_public_success','_dea_lpcontract_simplex','_dea_lpcontract_dispatch','_dea_lpcontract_ready']
    for helper in helpers:
        body+='lp x,min engine(simplex)\nassert r(status)==0 & r(objective)==2\n'
        body+=f'mata: mata drop {helper}()\n'
        body+='capture noisily _dea_lpcontract,c(C) aeq(E) beq(B) engine(simplex)\nassert _rc==498\n'
        body+='lp x,min engine(simplex)\nassert r(status)==0 & r(objective)==2\n'
        body+=f'mata: assert(_dea_lpcontract_ready()==9 & findexternal("{helper}()")!=NULL)\n'
        body+=f'mata: savebin("{helper}_x.bin",st_matrix("r(x)"))\n'
    out=run('adapter_cache',body)
    for helper in helpers:assert readbin(out/(helper+'_x.bin'))[0,0]==2
    (r/'ADAPTER_CACHE.json').write_text(json.dumps(dict(passed=True,helpers=helpers,public_reloads=4,internal_rejections=4,rejection_code=498),indent=2)+'\n')

def main():
    global r,args
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stata',required=True);parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args();r=args.output.resolve();r.mkdir(parents=True,exist_ok=False)
    (r/'checks').mkdir();shutil.copytree(REPO/'ado',r/'ado')
    (r/'manifest.json').write_text(json.dumps(dict(source_sha256={p.name:h(p) for p in (r/'ado').iterdir() if p.is_file()},fixture_sha256={p.name:h(p) for p in (REPO/'tests/fixtures').glob('success_quality.*')}),indent=2)+'\n')
    for fn in [check_accepted,check_controls,check_arithmetic,check_public_units,check_affine_refinement,check_bounds,check_adapter_cache]:fn()
    (r/'summary.json').write_text(json.dumps(dict(passed=True,batches=len(runs),accepted_conditions=139,boundary=1,arithmetic_conditions=1297,public_unit_commands=96,affine_controls=32,bounded_LP=24,MILP=9,solver_controls=48,direct_guards=3,public_adapter_reloads=4,internal_adapter_rejections=4,runs=runs),indent=2)+'\n')
if __name__=='__main__':main()
