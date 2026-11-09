"""Independent binary-rational witnesses for fuzzy center LPs.

Synthetic fixtures only; no Git history, research paths or installed R needed.
"""
from fractions import Fraction as F
import itertools
import numpy as np
from test_core_retries import readbin
from test_fuzzy_models import fixture
def mat(a):
 a=np.asarray(a,float)
 if a.ndim==1:a=a[None,:]
 if not a.size:return f'J({a.shape[0]},{a.shape[1]},0)'
 return '('+'\\'.join(','.join(format(float(v),'.17g') if np.isfinite(v) else '.' for v in row) for row in a)+')'
def fr(A):return [[F(float(v)) for v in row] for row in A]
def vec(p):return [F(float(v)) for v in readbin(p).ravel()]
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def integer(a):
 assert all(float(v).is_integer() for v in a) and a[0] in [-1,0,1]
 assert all(0<=v<2**26 for v in a[1:])
 return int(a[0])*sum(int(v)*2**(26*k) for k,v in enumerate(a[1:]))
def power(e):return F(2)**int(e)
def cert(d,tag):
 def b(k):return readbin(d/(tag+'_'+k+'.bin'))
 G=fr(b('G'));h=list(map(F,b('h').ravel()));T=fr(b('T'));info=b('info').ravel();ne=int(info[0]);assert info[1]==1
 n=len(T);assert len(T[0])==n and all(T[i][k]==0 for i in range(n) for k in range(n) if i!=k)
 assert all(q[i]>0 and (q[i].numerator&(q[i].numerator-1))==0 and (q[i].denominator&(q[i].denominator-1))==0 for i,q in enumerate(T))
 x=[F(integer(a),integer(b('den').ravel())) for a in b('num')]
 ids=[int(v)-1 for v in b('basis').ravel()];ex=b('ex').ravel()
 dual=[F(integer(a),integer(b('bden').ravel()))*power(ex[-1]-ex[i]) for i,a in enumerate(b('dual'))]
 if (d/(tag+'_cost.bin')).exists():c=[F(integer(a))*power(b('ce').item()) for a in b('cost')]
 else:c=[-v for v in vec(d/(tag+'_c.bin'))]
 assert len(ids)==n and len(set(ids))==n
 residual=[dot(row,x)-v for row,v in zip(G,h)]
 assert all(v==0 for v in residual[:ne]) and all(v<=0 for v in residual[ne:])
 assert all(v>=0 for k,v in zip(ids,dual) if k>=ne)
 assert all(dot([G[k][j] for k in ids],dual)==c[j] for j in range(n))
 assert all(residual[k]*v==0 for k,v in zip(ids,dual))
 assert dot(c,x)==dot([h[k] for k in ids],dual)
 face=[int(v)-1 for v in b('face').ravel()]
 assert face==[k for k,v in zip(ids,dual) if k>=ne and v>0]
 all_dual=[F(0)]*len(G)
 for i,v in zip(ids,dual):all_dual[i]=v
 return dict(dual=all_dual,G=G,h=h,T=T,x=x,c=c,ne=ne,face=face,raw=[T[i][i]*x[i] for i in range(n)],objective=-dot(c,x),ops=int(info[4]),limbs=int(info[5]))

DUMP='void dump135(string scalar tag,struct _F135_Center scalar R) {\n struct _C58_Cert scalar f\n savebin(tag+"_diag.bin",R.diag);savebin(tag+"_scores.bin",(R.score,R.radial))\n savebin(tag+"_weights.bin",R.weights);savebin(tag+"_lambda.bin",R.lambda)\n if(!R.ok)return\n f=R.primary\n savebin(tag+"_G.bin",R.G);savebin(tag+"_h.bin",R.gh);savebin(tag+"_T.bin",R.T);savebin(tag+"_c.bin",R.c*R.T)\n savebin(tag+"_num.bin",f.numerator);savebin(tag+"_den.bin",f.denominator);savebin(tag+"_dual.bin",f.dual_numerator)\n savebin(tag+"_bden.bin",f.basis_denominator);savebin(tag+"_ex.bin",f.row_exponent);savebin(tag+"_basis.bin",f.basis);savebin(tag+"_face.bin",f.face)\n savebin(tag+"_info.bin",(rows(R.E),f.status,f.reason,f.solver_status,f.ops,f.maxlimbs))\n}\nvoid probe135(string scalar tag,real matrix X,real matrix Y,string scalar rts,string scalar ort) {\n struct _F135_Center scalar R\n real scalar k\n savebin(tag+"_X.bin",X);savebin(tag+"_Y.bin",Y)\n for(k=1;k<=rows(X);k++) {R=_f135_stage(X,Y,k,rts,ort,1e-10);dump135(tag+"_"+strofreal(k),R);}\n}\n'

def grid_body():
    jobs=[];body='quietly _dea_fuzzy_load\nmata:\n'+DUMP
    body+='ctx135=_c57_context_save()\n_dea_rec_abort();_rsm_audit_reset(0);_rsm_success_reset(1,0);_dea_link_reset(2,0)\n'
    for m,s in [(1,1),(2,1),(1,2),(2,2)]:
     a,_=fixture(m,s,'near');X0=a[:,1:1+m];Y0=a[:,1+m:1+m+s]
     for variant in ['raw','reverse','dyadic','duplicate','zero_columns','columns_reverse']:
      X=X0.copy();Y=Y0.copy();ids=list(range(1,6))
      if variant=='reverse':X=X[::-1];Y=Y[::-1];ids=ids[::-1]
      if variant=='dyadic':X=X*np.array([2.**(16-32*j) for j in range(m)]);Y=Y*np.array([2.**(-24+48*j) for j in range(s)])
      if variant=='duplicate':X=np.r_[X,X[:1]];Y=np.r_[Y,Y[:1]];ids+=[1]
      if variant=='zero_columns':X=np.c_[X,np.zeros(len(X))];Y=np.c_[Y,np.zeros(len(Y))]
      if variant=='columns_reverse':X=X[:,::-1];Y=Y[:,::-1]
      for rts,ori in itertools.product(['crs','vrs','drs','irs'],['input','output']):
       name=f'{m}{s}_{variant}_{rts}_{ori}';jobs.append(dict(name=name,dataset=f'{m}{s}',variant=variant,ids=ids,rts=rts,orientation=ori,X=X.tolist(),Y=Y.tolist()))
       body+=f'probe135("{name}",{mat(X)},{mat(Y)},"{rts}","{ori}")\n'
    body+='''_c57_context_restore(ctx135)
    after135=_c57_context_save();names135=_c57_context_names()
    for(j135=1;j135<=cols(names135);j135++) {
     assert(asarray(ctx135,names135[j135]+":existed")==asarray(after135,names135[j135]+":existed"))
     if(asarray(ctx135,names135[j135]+":existed")) asserteq(asarray(ctx135,names135[j135]),asarray(after135,names135[j135]))
    }
    end
    '''
    return body,jobs

def check_grid(d,jobs):
    rows=[];invariants={}
    for job in jobs:
     name=job['name'];X,Y=fr(np.array(job['X'])),fr(np.array(job['Y']));n=len(X);m=len(X[0]);s=len(Y[0]);rts=job['rts'];ori=job['orientation']
     w=[F(1),F(-1)] if rts=='vrs' else [F(-1)] if rts=='drs' else [F(1)] if rts=='irs' else [];p=m+s+len(w)
     assert np.array_equal(readbin(d/(name+'_X.bin')),np.array(job['X'])) and np.array_equal(readbin(d/(name+'_Y.bin')),np.array(job['Y']))
     for k in range(n):
      tag=name+'_'+str(k+1);diag=readbin(d/(tag+'_diag.bin')).ravel();assert diag[0]==1,(tag,diag.tolist())
      q=cert(d,tag);T=q['T'];ex=X[k]+[F(0)]*(s+len(w));ey=[F(0)]*m+Y[k]+[F(0)]*len(w);intercept=[F(0)]*(m+s)+w
      E=ex if ori=='input' else ey;U=[[-v for v in x]+y+w for x,y in zip(X,Y)];c=[-v-z for v,z in zip(ey,intercept)] if ori=='input' else [v-z for v,z in zip(ex,intercept)]
      mapped=lambda a:[v*T[j][j] for j,v in enumerate(a)]
      assert q['ne']==1 and q['G']==[mapped(E)]+[mapped(v) for v in U]+[[-F(i==j) for j in range(p)] for i in range(p)]
      assert q['h']==[F(1)]+[F(0)]*(n+p) and q['c']==mapped([-v for v in c])
      radial=-q['objective'] if ori=='input' else q['objective'];score=radial if ori=='input' else 1/radial;lam=q['dual'][1:n+1]
      assert 0<score<=1 and all(v>=0 for v in lam)
      if rts=='vrs':assert sum(lam)==1
      if rts=='drs':assert sum(lam)<=1
      if rts=='irs':assert sum(lam)>=1
      assert all(dot([r[j] for r in X],lam)<=X[k][j]*(radial if ori=='input' else 1) for j in range(m))
      assert all(dot([r[j] for r in Y],lam)>=Y[k][j]*(radial if ori=='output' else 1) for j in range(s))
      out=readbin(d/(tag+'_scores.bin')).ravel();gap=abs(F(out[0])-score)/score;assert gap<=F(2e-16) and out[1]==float(radial)
      key=(job['dataset'],rts,ori,job['ids'][k]);invariants.setdefault(key,score);assert invariants[key]==score
      rows.append(dict(case=name,row=k+1,exact_score=str(score),score=float(out[0]),relative_export_error=float(gap),ops=q['ops'],limbs=q['limbs']))
    return dict(profiles=len(jobs),certificates=len(rows),raw_unique_scores=len(invariants),max_relative_export_error=max(r['relative_export_error'] for r in rows),max_ops=max(r['ops'] for r in rows),max_limbs=max(r['limbs'] for r in rows))
