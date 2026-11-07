#!/usr/bin/env python3
"""Additive Shapley / historical allocation conventions, independent DEA duals,
actual LP capture, sample contracts and forced-failure propagation. Stata 16+,
NumPy/SciPy; optional Benchmarking reference. No general Nash optimizer claim.
"""
from _support import package_files
from pathlib import Path
import argparse,hashlib,itertools,json,math,re,shutil,subprocess,time
import numpy as np
from test_core_retries import WRITER,readbin
from test_nonte_tolerance import instrument,CAPTURE
from test_lp_residuals import residual
from scipy.optimize import linprog
from test_output_protection import unchanged
from test_io_syntax import fixture as io_fixture,QUALIFIERS,CAPTURE as RET_CAPTURE,VERIFY
REPO=Path(__file__).resolve().parents[1]
DATA='input double(id x y)\n1 1 1\n2 2 1\n3 4 1\nend\n'
SUFFIX={'cooperative':'indiv shapley contrib rank','nash':'disagree coop bargain gain'}
COMMON=['n','n_lp','n_failed','n_boundary','complete','tolerance']
SPEC={'cooperative':['v_grand','detail_shown'],'nash':['total_gain','allocated_gain','unallocated_gain','mean_gain','nash_product','log_product','n_zero_gain','n_capped','product_underflow','allocation_complete']}

def fixture(m,s,variant):
 rng=np.random.default_rng(9113+10*m+s);a=rng.uniform(.2,7,(6,m+s));a[5,:m]=2*a[0,:m];a[5,m:]=.5*a[0,m:]
 if variant=='zero' and m>1:a[1,0]=0
 if variant=='zero' and s>1:a[2,m]=0
 if variant=='zerocol' and m>1:a[:,m-1]=0
 if variant=='zerocol' and s>1:a[:,-1]=0
 if variant=='units':a*=10.**np.linspace(-12,12,m+s)
 if variant=='tiny':a*=1e-12
 if variant=='duplicate':a[1]=a[0]
 if variant=='near':a[1]=a[0]*(1+1e-9)
 selected=np.ones(len(a));z=np.c_[np.arange(1,len(a)+1),a,selected]
 if variant=='sample':z[1,-1]=0;z=np.vstack([z,z[0]]);z[-1,0]=7;z[-1,1:1+m]*=.001;z[-1,-1]=0
 if variant=='permuted':z=z[[5,2,0,4,1,3]]
 return z,['id']+[f'x{i+1}' for i in range(m)]+[f'y{i+1}' for i in range(s)]+['selected']

def dual_score(x,y,k,rts):
 # Independent multiplier model (production uses bounded envelopment).
 n,m=x.shape;ss=y.shape[1]
 c=np.r_[np.zeros(m),-y[k]];E=np.r_[x[k],np.zeros(ss)][None,:]
 U=np.c_[-x,y];bounds=[(0,None)]*(m+ss)
 if rts=='vrs':
  U=np.c_[U,-np.ones(n)];E=np.c_[E,0];c=np.r_[c,1];bounds.append((None,None))
 f=linprog(c,A_ub=U,b_ub=np.zeros(n),A_eq=E,b_eq=[1],bounds=bounds,method='highs-ipm',options=dict(presolve=False,primal_feasibility_tolerance=1e-10,dual_feasibility_tolerance=1e-10,ipm_optimality_tolerance=1e-12))
 assert f.success,f.message
 assert max(U@f.x,default=0)<1e-8 and abs((E@f.x)[0]-1)<1e-8
 assert -1e-8<=-f.fun<=1+1e-8
 return -f.fun


def verify(d,m,s,mode,tol):
 a=readbin(d/'data.bin');x=a[:,1:1+m];y=a[:,1+m:];x=x/np.where(x.max(0)==0,1,x.max(0));y=y/np.where(y.max(0)==0,1,y.max(0));n=len(a)
 R=readbin(d/'results.bin');D=readbin(d/'lp_diagnostics.bin');M=readbin(d/'meta.bin').ravel();rts=(d/'rts.txt').read_text()
 want=np.array([dual_score(x,y,i,rts) for i in range(n)])
 good=np.isfinite(D[:,1]);gap=np.max(abs(D[good,1]-want[good]),initial=0);assert gap<3e-6,(gap,R,want)
 assert np.array_equal(good,D[:,2]==1);assert np.array_equal(R[:,0],D[:,1],equal_nan=True)
 assert M[0]==M[1]==n and M[2]==sum(~good) and M[3]==sum(D[:,3]) and M[4]==good.all() and M[5]==tol
 assert np.all((D[good,1]>=0)&(D[good,1]<=1))
 assert np.array_equal(D[good,3],D[good,0]!=D[good,1]);assert np.max(abs(D[good,0]-D[good,1]),initial=0)<=min(tol,1e-8)
 e=R[:,0];shapley_checks=0
 if mode=='cooperative':
  assert np.array_equal(R[:,0],R[:,1],equal_nan=True);assert np.all(R[good,2]==0) and np.isnan(R[~good,2]).all()
  if good.all():
   assert abs(M[6]-sum(e))<1e-13;assert np.array_equal(R[:,3],1+(e[:,None]<e[None,:]).sum(1))
   # Independently enumerate the original characteristic function, not e_i directly.
   v=np.array([sum(e[j] for j in range(n) if k&(1<<j)) for k in range(1<<n)])
   val=np.zeros(n)
   for i in range(n):
    for k in range(1<<n):
     if not k&(1<<i):
      t=k.bit_count();val[i]+=math.factorial(t)*math.factorial(n-t-1)/math.factorial(n)*(v[k|(1<<i)]-v[k])
   assert np.max(abs(val-R[:,1]))<2e-14;shapley_checks=n
  else:assert np.isnan(M[6]) and np.isnan(R[:,3]).all()
 else:
  c=np.where(e>=1-tol,1,e/(.8+.2*e));assert np.allclose(c,R[:,1],atol=2e-15,rtol=0,equal_nan=True)
  if good.all():
   G=sum(c-e);w=np.maximum(.01,1-e);raw=G*w/sum(w);b=np.minimum(1,e+raw);g=b-e
   assert np.allclose(R[:,2],b,atol=3e-15,rtol=0);assert np.array_equal(R[:,3],R[:,2]-e)
   assert np.all(R[:,3]>=0);assert np.allclose(M[6:10],[G,sum(g),G-sum(g),np.mean(g)],atol=4e-15,rtol=0)
   assert M[12]==sum(g==0) and M[13]==sum(e+raw>1) and M[15]==1
   if np.any(g==0):assert M[10]==0 and np.isnan(M[11]) and M[14]==0
   else:assert abs(M[11]-sum(np.log(g)))<1e-11 and abs(M[10]-np.exp(M[11]))<1e-15
  else:assert np.isnan(R[:,2:]).all() and np.isnan(M[6:15]).all() and M[15]==0
 return dict(independent_dual_LPs=n,max_score_gap=float(gap),unavailable=int(sum(~good)),additive_shapley_values=shapley_checks)

def verify_lp(d,tol):
 # Self reference provides an explicit witness even when a presolver declares
 # a nearly duplicated VRS model infeasible. Retain that presolve status, and
 # solve the *same* LP using HiGHS IPM without presolve at tighter tolerance.
 counts=dict(lp_checks=0,accepted=0,rejected=0,rejected_feasible=0,restarts=0,
             oracle_presolve_infeasible=0,max_scaled=0.,max_relative_objective_gap=0.)
 for p in sorted(d.glob('lp*_meta.bin')):
  base=p.name[:-9];status,obj,requested,actual,sign=readbin(p).ravel()
  assert requested==min(tol,1e-8) and actual<=requested and sign==1
  c=readbin(d/(base+'_c.bin')).ravel();E=readbin(d/(base+'_E.bin'));U=readbin(d/(base+'_U.bin'));B=readbin(d/(base+'_B.bin'));x=readbin(d/(base+'_x.bin')).ravel()
  n=len(c)-1;k=int(base[2:])-1;assert k<n
  witness=np.zeros(n+1);witness[k]=1;witness[-1]=1
  lb=np.nan_to_num(B[0],nan=-np.inf);ub=np.nan_to_num(B[1],nan=np.inf)
  assert np.max(U[:,:-1]@witness-U[:,-1],initial=0)==0
  assert np.max(abs(E[:,:-1]@witness-E[:,-1]),initial=0)==0
  assert np.all(witness>=lb) and np.all(witness<=ub)
  def norm(M):
   if not len(M):return None,None
   z=np.maximum(abs(M).max(1),1e-300);return M[:,:-1]/z[:,None],M[:,-1]/z
  ae,be=norm(E);au,bu=norm(U);kwargs=dict(A_eq=ae,b_eq=be,A_ub=au,b_ub=bu,bounds=list(zip(lb,ub)))
  first=linprog(c,**kwargs,method='highs',options=dict(primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9))
  fit=linprog(c,**kwargs,method='highs-ipm',options=dict(presolve=False,primal_feasibility_tolerance=1e-10,dual_feasibility_tolerance=1e-10,ipm_optimality_tolerance=1e-12))
  assert fit.success,(p.name,fit.message)
  rr=residual(fit.x,E[:,:-1],E[:,-1],U[:,:-1],U[:,-1],lb,ub)
  assert max(rr[3:])<3e-9 and -1e-9<=fit.fun<=1+1e-9,(p.name,rr,fit.fun)
  (d/(base+'_oracle.json')).write_text(json.dumps(dict(presolve_status=int(first.status),verification_status=int(fit.status),objective=float(fit.fun),scaled_residual=float(max(rr[3:])),self_witness_exact_float_residual=0))+'\n')
  counts['oracle_presolve_infeasible']+=int(first.status==2);counts['lp_checks']+=1
  attempts=sorted(d.glob('attempt'+base[2:]+'_*.bin'));assert 1<=len(attempts)<=2
  aa=readbin(attempts[0]).ravel();assert aa[0]==min(tol,1e-8) and aa[1]==min(tol,1e-8)
  if len(attempts)==2:
   counts['restarts']+=1;bb=readbin(attempts[1]).ravel();assert aa[2]==1 and aa[3]==0 and bb[1]==min(tol,1e-12)
  if status==0:
   counts['accepted']+=1;assert np.isfinite(x).all() and np.isfinite(obj)
   rr=residual(x,E[:,:-1],E[:,-1],U[:,:-1],U[:,-1],lb,ub);scaled=float(max(rr[3:]));assert scaled<=min(tol,1e-8)+5e-16
   gap=abs(obj-fit.fun);assert gap<1e-6 and abs(obj-c@x)<1e-8
   counts['max_scaled']=max(counts['max_scaled'],scaled);counts['max_relative_objective_gap']=max(counts['max_relative_objective_gap'],float(gap))
  else:counts['rejected']+=1;counts['rejected_feasible']+=1
 return counts


def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--only');ap.add_argument('--rscript');ap.add_argument('--baseline-ado',type=Path);a=ap.parse_args()
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir();h={}
 for n in package_files(REPO):
  raw=(REPO/n).read_bytes();(pkg/Path(n).name).write_bytes(raw);h[n]=hashlib.sha256(raw).hexdigest()
 audit=out/'audit_ado';instrument(pkg,audit)
 (out/'manifest.json').write_text(json.dumps(dict(source_sha256=h,stata=a.stata),indent=2)+'\n');rows=[]
 def run(name,body,source=pkg,check=None):
  if a.only and not re.search(a.only,name):return
  d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset type double\nset processors 1\nset linesize 180\n'
  for f in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/f).mkdir();head+=f'sysdir set {f} "{d/f}"\n'
  head+=f'adopath ++ "{source}"\n'+WRITER+CAPTURE+'scalar r08calls=0\n';(d/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n');r=dict(test=name,passed=False);start=time.monotonic()
  try:
   p=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=180)
   log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else ''
   if p.returncode or not log:(d/'process.json').write_text(json.dumps(dict(returncode=p.returncode,stdout=p.stdout.decode(errors='replace'),stderr=p.stderr.decode(errors='replace'))))
   assert p.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2500:]
   if check:r.update(check(d))
   r['passed']=True
  except Exception as e:r['error']=repr(e)
  r['seconds']=round(time.monotonic()-start,3);rows.append(r);(out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(r),flush=True)
  return d
 specs=[]
 for m,s in [(1,1),(2,1),(1,2),(2,2),(3,2)]:
  for var in ['base','units','tiny','zero','zerocol','duplicate','near','sample','permuted']:
   z,cols=fixture(m,s,var);key=f'{m}{s}_{var}';f=out/('data_'+key+'.csv');np.savetxt(f,z,delimiter=',',header=','.join(cols),comments='',fmt='%.17g');specs.append((key,m,s,f))
   xs=' '.join(c for c in cols if c.startswith('x'));ys=' '.join(c for c in cols if c.startswith('y'))
   for mode,rts,tol in itertools.product(SUFFIX,['crs','vrs'],['1e-6','1e-8','1e-10']):
    body=f'import delimited using "{f}",clear asdouble\ndea_{mode} {xs}={ys} if selected, dmu(id) rts({rts}) tol({tol}) gen(res) matrix detail debug\n'
    body+='matrix meta=('+','.join('r('+k+')' for k in COMMON+SPEC[mode])+')\nmata: savebin("meta.bin",st_matrix("meta"))\n'
    for mat in ['results','lp_diagnostics']:body+=f'mata: savebin("{mat}.bin",st_matrix("r({mat})"))\n'
    body+='foreach v of varlist res_* {\nassert missing(`v\') if !selected\n}\nkeep if selected\n'
    body+=f'mata: savebin("data.bin",st_data(.,tokens("id {xs} {ys}")))\n'
    body+=f'file open rr using rts.txt,write\nfile write rr "{rts}"\nfile close rr\n'
    run(key+'_'+mode+'_'+rts+'_'+tol,body,check=lambda p,m=m,s=s,mode=mode,tol=tol:verify(p,m,s,mode,float(tol)))
   if var in ['base','units','near']:
    for mode,rts,tol in itertools.product(SUFFIX,['crs','vrs'],['1e-4','1e-7','1e-10']):
     body=f'import delimited using "{f}",clear asdouble\ndea_{mode} {xs}={ys},rts({rts}) tol({tol}) matrix\nassert r08calls==r(n)\n'
     run('lp_'+key+'_'+mode+'_'+rts+'_'+tol,body,source=audit,check=lambda d,tol=tol:verify_lp(d,float(tol)))
 for mode in SUFFIX:
  cmd=f'dea_{mode},inputs(x) outputs(y) gen(res)'
  for opt in ['tol(0)','tol(-1)','tol(.)','tol(.a)','rts(irs)']:
   run(mode+'_invalid_'+re.sub(r'\W','_',opt),DATA+'gen str12 res_'+SUFFIX[mode].split()[0]+'="old"\n'+unchanged(cmd+' replacevars '+opt,198))
  for v,val in itertools.product(['x','y'],[-1,0]):run(f'{mode}_domain_{v}_{val}',DATA+f'replace {v}={val} in 1\n'+unchanged(cmd,459))
  run(mode+'_empty',DATA+unchanged(cmd.replace(',', ' if 0,',1),2000))
  run(mode+'_string_quantity',DATA+'tostring x,replace\n'+unchanged(cmd,109))
  for kind in ['numeric','string']:
   body=DATA+('replace id=. in 1\n' if kind=='numeric' else 'tostring id,replace\nreplace id="" in 1\n')+cmd+' dmu(id) matrix\nassert r(n)==2\nassert missing(res_'+SUFFIX[mode].split()[0]+') in 1\n'
   run(mode+'_missing_id_'+kind,body)
  run(mode+'_abbreviated_options',DATA+cmd+' rts(vrs) tol(1e-8) mat det debug\n'+RET_CAPTURE+'\n'+cmd+' rts(vrs) tolerance(1e-8) matrix detail debug replacevars\n'+VERIFY)
  run(mode+'_singleton',DATA+'keep in 1\n'+cmd+' matrix detail\nassert r(complete)==1\nassert abs(res_'+SUFFIX[mode].split()[0]+'-1)<1e-8\n')
  run(mode+'_reload_workspace',DATA+'matrix user=(9,8)\nscalar __dea_coop_v_grand=7654321\nlocal before : all matrices\n'+cmd+' matrix\nlocal after : all matrices\nassert "`before\'"=="`after\'"\nassert __dea_coop_v_grand==7654321\ndiscard\n'+cmd+' replacevars matrix\nassert __dea_coop_v_grand==7654321\nlocal after : all matrices\nassert "`before\'"=="`after\'"\n')
  run(mode+'_large24', 'set obs 24\ngen x=_n\ngen y=1\n'+cmd+' matrix detail\nassert r(n)==24 & r(n_lp)==24 & r(complete)==1\nassert abs(res_'+SUFFIX[mode].split()[0]+'-1/x)<1e-7\n'+('assert r(detail_shown)==0\n' if mode=='cooperative' else ''))
  run(mode+'_numeric_id',DATA+'replace id=1234567890123+_n\n'+cmd+' dmu(id)\n',check=lambda d:dict(id_labels_checked=3) if all(str(1234567890123+i) in (d/'test.log').read_text() for i in range(1,4)) else (_ for _ in ()).throw(AssertionError('ID precision')))
  for rts,label,form in itertools.product(['crs','vrs'],QUALIFIERS,['A','B']):
   q=QUALIFIERS[label];outs=' '.join('res_'+s for s in SUFFIX[mode].split())
   c=f'dea_{mode} '+('x1 x2=y1 y2,' if form=='A' else ',inputs(x1 x2) outputs(y1 y2)')+f' dmu(sid) gen(res) rts({rts}) matrix detail debug'
   body=io_fixture('nash')+f'gen byte selected=0\nreplace selected=1 {q}\ntempfile expected\npreserve\nkeep if selected\n'+c+'\n'+RET_CAPTURE+'keep obsid '+outs+'\nsave `expected\'\nrestore\n'+c.replace(',',f' {q},',1)+'\n'+VERIFY+'foreach v of varlist '+outs+' {\nassert missing(`v\') if !selected\n}\nkeep if selected\ncf obsid '+outs+' using `expected\',all\n'
   run(mode+'_sample_'+rts+'_'+label+'_'+form,body)
 for fault in ['first','all','bound','nonfinite','above','below','outside']:
  target=out/('fault_ado_'+fault);shutil.copytree(pkg,target)
  if fault in ['bound','nonfinite']:
   p=target/'_dea_nonte_lp.ado';s=p.read_text().replace('        x=S.parameters()', '        x=S.parameters()\n        '+('x[1]=-1000' if fault=='bound' else 'x[1]=.'));p.write_text(s)
  else:
   p=target/'_dea_games.ado';s=p.read_text()
   if fault in ['first','all']:s=s.replace('        S=LinearProgram();',('        if(i==1) continue\n' if fault=='first' else '        continue\n')+'        S=LinearProgram();')
   else:s=s.replace('v=S.value();D[i,1]=v', 'v='+{'above':'1+5e-10','below':'-5e-10','outside':'1+1e-5'}[fault]+';D[i,1]=v')
   p.write_text(s)
  for mode in SUFFIX:
   b=DATA+f'dea_{mode} x=y,gen(res) matrix detail tol(1e-7)\n'
   if fault in ['above','below']:
    b+='assert r(n_boundary)==3 & r(n_failed)==0\nassert res_'+SUFFIX[mode].split()[0]+'=='+('1' if fault=='above' else '0')+'\n'
   else:
    b+=f'assert r(n_failed)=={1 if fault=="first" else 3} & r(complete)==0\n'
    if mode=='cooperative':b+='assert missing(r(v_grand))\nassert missing(res_rank)\nassert missing(res_shapley,res_contrib) in 1\n'
    else:b+='assert missing(res_coop) in 1\nassert missing(res_bargain,res_gain)\nassert missing(r(mean_gain),r(nash_product),r(total_gain))\nassert r(allocation_complete)==0\n'
   run('fault_'+fault+'_'+mode,b,source=target)
 pure='quietly run "'+str(pkg/'_dea_games.ado')+'"\nmata:\n'
 pure+='m=J(1,10,.);R=_dea_games_allocate((1\\.5\\.25),1e-7,m)\nassert(m[3]>0 & m[8]==1 & m[5]==0 & m[7]==1 & m[6]>=.)\nassert(R[1,4]==0 & all(R[,4]:>=0))\n'
 pure+='R=_dea_games_allocate(J(1000,1,.5),1e-7,m)\nassert(m[5]==0 & m[9]==1 & m[7]==0 & m[6]<-1000)\n'
 pure+='R=_dea_games_allocate((.a\\.5),1e-7,m)\nassert(R[1,2]>=. & all(R[,3..4]:>=.) & m[10]==0 & m[5]>=.)\n'
 pure+='R=_dea_games_allocate(J(3,1,1),1e-7,m)\nassert(m[1]==0 & m[2]==0 & m[3]==0 & m[5]==0 & m[7]==3 & m[8]==0 & m[9]==0 & m[10]==1)\n'
 pure+='R=_dea_games_allocate(J(3,1,0),1e-7,m)\nassert(m[1]==0 & m[5]==0 & m[7]==3 & m[10]==1)\n'
 pure+='assert(_dea_games_coalition((.\\.5),0)==0 & _dea_games_coalition((.\\.5),1)>=. & _dea_games_coalition((.\\.5),2)==.5)\n'
 pure+='assert(all(_dea_games_rank((.5\\.5\\.25)):==(1\\1\\3)))\nend\n'
 run('pure_boundary_and_underflow',pure)
 if a.baseline_ado:
  baseline=out/'before_ado';shutil.copytree(a.baseline_ado,baseline)
  countpkg=out/'old_count_ado';shutil.copytree(baseline,countpkg)
  pp=countpkg/'dea_cooperative.ado';ss=pp.read_text().replace('(void) S.optimize()', 'st_numscalar("r08calls",st_numscalar("r08calls")+1)\n    (void) S.optimize()');pp.write_text(ss)
  run('old_repeated_lp_count',DATA+'dea_cooperative x=y\nassert r08calls==15\n',source=countpkg)
  run('old_scalar_overwrite',DATA+'scalar __dea_coop_v_grand=7654321\ndea_cooperative x=y\ncapture confirm scalar __dea_coop_v_grand\nassert _rc==111\n',source=baseline)
  run('old_limit24','set obs 24\ngen x=_n\ngen y=1\n'+unchanged('dea_cooperative x=y',198),source=baseline)
  for mode,fn in [('cooperative','_dea_coop_efficiency'),('nash','_dea_nash_eff')]:
   target=out/('old_fault_'+mode);shutil.copytree(baseline,target);p=target/f'dea_{mode}.ado';s=p.read_text();pos=s.index('real scalar '+fn+'(');pos=s.index('{',pos)+1;s=s[:pos]+'\n    if(x0[1]==1) return(.)\n'+s[pos:];p.write_text(s)
   body=DATA+f'dea_{mode} x=y,gen(res)\n'
   body+=('assert missing(res_indiv) in 1\nassert res_shapley==0 in 1\nassert r(v_grand)<.\n' if mode=='cooperative' else 'assert missing(res_disagree) in 1\nassert res_coop==1 in 1\nassert !missing(res_bargain,res_gain) in 2/3\n')
   run('old_missing_'+mode,body,source=target)
 if a.rscript:
  script='''library(Benchmarking)
a<-commandArgs(TRUE);d<-read.csv(a[1]);m<-as.integer(a[2]);s<-as.integer(a[3]);d<-d[d$selected==1,];x<-as.matrix(d[,2:(m+1),drop=FALSE]);y<-as.matrix(d[,(m+2):(m+s+1),drop=FALSE]);x<-sweep(x,2,apply(x,2,function(v)max(c(v,1e-300))),"/");y<-sweep(y,2,apply(y,2,function(v)max(c(v,1e-300))),"/");write.csv(data.frame(crs=dea(x,y,RTS="crs",ORIENTATION="in")$eff,vrs=dea(x,y,RTS="vrs",ORIENTATION="in")$eff),a[4],row.names=FALSE);cat(as.character(packageVersion("Benchmarking")))
'''
  (out/'compare.R').write_text(script)
  def rc(d):
   count=0;gap=0;missing=0
   for key,m,s,f in specs:
    if not key.endswith('_base'):continue
    dest=d/(key+'.csv');p=subprocess.run([a.rscript,str(out/'compare.R'),str(f),str(m),str(s),str(dest)],capture_output=True,text=True,timeout=180);(d/(key+'.console')).write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr;ref=np.genfromtxt(dest,delimiter=',',names=True)
    for mode,rts in itertools.product(SUFFIX,['crs','vrs']):
     q=out/(key+'_'+mode+'_'+rts+'_1e-10')
     if not q.exists():continue
     got=readbin(q/'results.bin')[:,0];ok=np.isfinite(got);missing+=sum(~ok);g=max(abs(got[ok]-ref[rts][ok]),default=0);assert g<3e-6;gap=max(gap,g);count+=sum(ok)
   assert count>0,'R comparison needs base grid cases (remove --only restriction)'
   return dict(R_score_comparisons=int(count),R_max_gap=float(gap),R_comparison_unavailable=int(missing))
  run('R_Benchmarking','',check=rc)
 summary=dict(total=len(rows),passed=sum(r['passed'] for r in rows),failed=sum(not r['passed'] for r in rows));(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary);return bool(summary['failed'])
if __name__=='__main__':raise SystemExit(main())
