"""Independent Fraction check of original LP optimality, faces and exports.
No numerical optimizer or Mata arithmetic is used to authorize a result.
"""
from pathlib import Path
from fractions import Fraction as F
import json,struct
import numpy as np
from test_core_retries import readbin

def integer(row):
 row=list(row);assert row[0] in [-1,0,1]
 assert all(v==int(v) and 0<=v<2**26 for v in row[1:])
 return int(row[0])*sum(int(v)*2**(26*i) for i,v in enumerate(row[1:]))
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def vec(a):return [F(float(x)) for x in np.ravel(a)]
def matrix(a):return [vec(row) for row in a]
def load(d,tag):return {k:readbin(d/(tag+'_'+k+'.bin')) for k in ['c','E','b','U','ub','lo','hi']}
def model(p):
 n=p['c'].size;upper=np.flatnonzero(np.isfinite(p['hi'].ravel()))
 return matrix(np.vstack([p['E'],p['U'],-np.eye(n),np.eye(n)[upper]])),vec(np.r_[p['b'].ravel(),p['ub'].ravel(),-p['lo'].ravel(),p['hi'].ravel()[upper]]),len(p['E'])
def certificate(d,tag):
 p=load(d,tag);read=lambda k:readbin(d/(tag+'_'+k+'.bin'));status=int(read('meta').ravel()[0]);assert status in [-1,0,1]
 if status==-1:return None
 D=integer(read('den').ravel());Db=integer(read('basis_den').ravel());assert D>0 and Db>0
 x=[F(integer(row),D) for row in read('num')];G,h,ne=model(p);n=len(x);ids=read('basis').ravel().astype(int)-1
 assert len(ids)==n and len(set(ids))==n and min(ids)>=0 and max(ids)<len(G)
 Q=[integer(row) for row in read('dual')];es=read('exponent').ravel();assert len(es)==n+1 and all(v==int(v) for v in es)
 y=[F(0)]*len(G)
 for k,i in enumerate(ids):y[i]=F(Q[k],Db)*F(2)**int(es[-1]-es[k])
 residual=[dot(row,x)-b for row,b in zip(G,h)];sign=lambda v:int(v>0)-int(v<0)
 assert [sign(v) for v in residual]==read('primal_sign').ravel().tolist()
 assert [sign(y[i]) for i in ids]==read('dual_sign').ravel().tolist()
 primal=all(v==0 if i<ne else v<=0 for i,v in enumerate(residual));dual=all(v>=0 for v in y[ne:])
 c=vec(p['c']);assert all(sum(row[j]*yy for row,yy in zip(G,y))==c[j] for j in range(n))
 optimal=primal and dual;assert (status==1)==optimal
 if not optimal:return None
 objective=dot(c,x);assert objective==dot(h,y)
 face=[i for i in ids if i>=ne and y[i]>0]
 assert sorted(face)==sorted((read('face').ravel().astype(int)-1).tolist())
 return dict(p=p,x=x,objective=objective,G=G,h=h,ne=ne,face=face)
def face(p,proof):
 E=p['E'];b=p['b']
 if proof['face']:
  E=np.vstack([E,np.array([[float(v) for v in proof['G'][i]] for i in proof['face']])]);b=np.vstack([b,np.array([float(proof['h'][i]) for i in proof['face']]).reshape(-1,1)])
 return E,b
def same(A,B):assert A.shape==B.shape and np.array_equal(A,B,equal_nan=True)
def compare(got,want,context):
 assert struct.pack('<d',float(got))==struct.pack('<d',float(want)),(context,got,float(want))

def audit_pipeline(d,X,Y,weight):
 n,m=X.shape;s=Y.shape[1];sx=X.max(0);sy=Y.max(0);xs=X/sx;ys=Y/sy;v=n+m+s+1
 payload={k:readbin(d/(k+'.bin')) for k in ['scores','input_total','input_technical','input_congestion','output_slack','target_inputs','target_outputs','lambda_stage1','lambda_stage2','projection_min','projection_max','congestion_min','congestion_max','range_status']}
 count=fields=0
 for i in range(n):
  a=d/'audit'/f'dmu_{i+1}';tags=sorted(f.name.removesuffix('_meta.bin') for f in a.glob('*_meta.bin') if '_attempt' not in f.name and '_repair' not in f.name);proofs={}
  for tag in tags:proofs[tag]=certificate(a,tag);count+=1
  radial=proofs.get('radial');slack=proofs.get('slack');base=proofs.get('decomposition')
  if radial:
   p=radial['p'];same(p['E'],np.vstack([np.c_[xs.T,np.eye(m),np.zeros((m,s+1))],np.c_[ys.T,np.zeros((s,m)),-np.eye(s),-ys[i]],np.r_[np.ones(n),np.zeros(m+s+1)]]));same(p['b'],np.r_[xs[i],np.zeros(s),1].reshape(-1,1))
  if slack:
   E,b=face(radial['p'],radial);same(slack['p']['E'],E);same(slack['p']['b'],b)
   w=np.r_[sx,sy] if weight=='raw' else 1/np.r_[xs[i],ys[i]];w=w/max(w);same(slack['p']['c'],np.r_[np.zeros(n),w,0].reshape(1,-1))
   assert slack['x'][v-1]==radial['objective']
  if base:
   ids=readbin(a/'slack_basis.bin').ravel().astype(int)-1
   B=np.array([[float(x) for x in slack['G'][j]] for j in ids]);hb=np.array([float(slack['h'][j]) for j in ids]).reshape(-1,1)
   same(base['p']['E'],np.vstack([np.c_[B,np.zeros((v,n))],np.c_[-ys.T,np.zeros((s,v-n)),ys.T],np.r_[np.zeros(v),np.ones(n)]]))
   same(base['p']['b'],np.vstack([hb,np.zeros((s,1)),[[1]]]))
   same(base['p']['U'],np.vstack([np.c_[np.zeros((m,v)),xs.T],np.c_[xs.T,np.zeros((m,v-n)),-xs.T]]))
   same(base['p']['ub'],np.r_[xs[i],np.zeros(m)].reshape(-1,1))
   w=sx if weight=='raw' else 1/xs[i];w=w/max(w)
   for j,value in enumerate(base['p']['c'].ravel()[-n:]):
    exact=dot(vec(xs[j]),vec(w));assert abs(F(float(value))-exact)<=(2*m+1)*F(2)**-52*abs(exact)
   point=base['x'];lam=point[:n];mu=point[-n:];sl=point[n:n+m];so=point[n+m:n+m+s];xx=matrix(xs.T);yy=matrix(ys.T);xi=[row[i] for row in xx]
   assert point[:v]==slack['x'];phi=point[n+m+s];tech=[dot(row,[u-l for u,l in zip(mu,lam)]) for row in xx];co=[x-y for x,y in zip(sl,tech)]
   assert all(0<=t<=a for t,a in zip(tech,sl))
   total=sum(x/y for x,y in zip(sl,xi))/m;technical=sum(x/y for x,y in zip(tech,xi))/m
   expected=dict(scores=[total-technical,phi,1/phi,total,technical,F(0)],lambda_stage1=lam,lambda_stage2=mu,input_total=[x*F(float(u)) for x,u in zip(sl,sx)],input_technical=[x*F(float(u)) for x,u in zip(tech,sx)],input_congestion=[x*F(float(u)) for x,u in zip(co,sx)],output_slack=[x*F(float(u)) for x,u in zip(so,sy)],target_inputs=[dot(row,lam)*F(float(u)) for row,u in zip(xx,sx)],target_outputs=[dot(row,lam)*F(float(u)) for row,u in zip(yy,sy)])
   for key,nums in expected.items():
    for x,y in zip(payload[key][i],nums):compare(x,y,(i,key));fields+=1
   Fd,bd=face(base['p'],base);Fs,bs=face(slack['p'],slack)
   for kind,dim,E,b in [('projection',m+s,Fs,bs),('congestion',m,Fd,bd)]:
    for j in range(dim):
     for sign in [-1,1]:
      tag=f'{kind}_{j+1}_{sign}';proof=proofs.get(tag)
      if not proof:continue
      same(proof['p']['E'],E);same(proof['p']['b'],b)
      q=np.r_[(xs.T if j<m else ys.T)[j if j<m else j-m],np.zeros(v-n)] if kind=='projection' else np.r_[np.zeros(len(point)-n),xs[:,j]]
      same(proof['p']['c'],(sign*q).reshape(1,-1))
      value=sign*proof['objective'];key=kind+('_min' if sign<0 else '_max')
      scale=(np.r_[sx,sy][j] if kind=='projection' else sx[j])
      if kind=='congestion':value=xi[j]-value;key=kind+('_min' if sign>0 else '_max')
      if payload['range_status'][i,0 if kind=='projection' else 1]==0:compare(payload[key][i,j],value*F(float(scale)),(i,key,j));fields+=1
 return dict(exact_final_basis_checks=count,exact_numeric_fields_bitwise_verified=fields,range_reference='Fraction primal/dual certificates on original exact faces;floating range oracle not used for exact acceptance')
