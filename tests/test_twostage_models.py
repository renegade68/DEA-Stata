#!/usr/bin/env python3
"""Independent two-stage LP formulas, shared-weight decomposition and public contracts.
Stata 16+, NumPy/SciPy; optional original R DJL/Benchmarking comparisons.
"""
from _support import package_files
from pathlib import Path
import argparse,csv,hashlib,io,itertools,json,re,shutil,subprocess,time
import numpy as np
from scipy.optimize import linprog
from test_core_retries import WRITER,readbin
from test_output_protection import unchanged
from test_io_syntax import CAPTURE,VERIFY,QUALIFIERS,fixture as io_fixture
REPO=Path(__file__).resolve().parents[1]
PUBLISHED_DATA='"dmu","x1","x2","z1","z2","y1","y2"\n1,1178744,673512,7451757,856735,984143,681687\n2,1381822,1352755,10020274,1812894,1228502,834754\n3,1177494,592790,4776548,560244,293613,658428\n4,601320,594259,3174851,371863,248709,177331\n5,6699063,3531614,37392862,1753794,7851229,3925272\n6,2627707,668363,9747908,952326,1713598,415058\n7,1942833,1443100,10685457,643412,2239593,439039\n8,3789001,1873530,17267266,1134600,3899530,622868\n9,1567746,950432,11473162,546337,1043778,264098\n10,1303249,1298470,8210389,504528,1697941,554806\n11,1962448,672414,7222378,643178,1486014,18259\n12,2592790,650952,9434406,1118489,1574191,909295\n13,2609941,1368802,13921464,811343,3609236,223047\n14,1396002,988888,7396396,465509,1401200,332283\n15,2184944,651063,10422297,749893,3355197,555482\n16,1211716,415071,5606013,402881,854054,197947\n17,1453797,1085019,7695461,342489,3144484,371984\n18,757515,547997,3631484,995620,692731,163927\n19,159422,182338,1141950,483291,519121,46857\n20,145442,53518,316829,131920,355624,26537\n21,84171,26224,225888,40542,51950,6491\n22,15993,10502,52063,14574,82141,4181\n23,54693,28408,245910,49864,0.1,18980\n24,163297,235094,476419,644816,142370,16976\n25,1544215,828963,7832893,667964,1602873,477733\n'
MATS=['efficiencies','status','weights_input','weights_intermediate','weights_output','decomposition_status','follower_status','follower_route','t121_quality','t122_quality1','t122_quality2','lambda_stage1','lambda_stage2']

def fit(c,A,b,E=None,h=None,bounds=(0,None)):
 return linprog(c,A_ub=A,b_ub=b,A_eq=E,b_eq=h,bounds=bounds,method='highs',options=dict(primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9))

def stage(x,y,k,rts,ort):
 n,m=x.shape;s=y.shape[1]
 # Independent multiplier dual, unlike production's envelopment formulation.
 if ort=='input':
  c=np.r_[np.zeros(m),-y[k],0 if rts=='vrs' else []];eq=np.r_[x[k],np.zeros(s),0 if rts=='vrs' else []]
 else:
  c=np.r_[x[k],np.zeros(s),1 if rts=='vrs' else []];eq=np.r_[np.zeros(m),y[k],0 if rts=='vrs' else []]
 A=np.c_[-x,y];bd=[(0,None)]*(m+s)
 if rts=='vrs':A=np.c_[A,-np.ones(n)];bd.append((None,None));c[-1]=1
 f=fit(c,A,np.zeros(n),eq[None,:],np.ones(1),bd)
 assert f.success,f.message
 return -f.fun if ort=='input' else 1/f.fun

def joint(x,z,y,k,mode,lock=None):
 n,m=x.shape;d=z.shape[1];s=y.shape[1]
 ex=np.r_[x[k],np.zeros(d+s)];ez=np.r_[np.zeros(m),z[k],np.zeros(s)];ey=np.r_[np.zeros(m+d),y[k]]
 A=np.vstack([np.c_[-x,z,np.zeros((n,s))],np.c_[np.zeros((n,m)),-z,y]])
 if mode=='primary':c=-ey;E=ex[None,:];h=[1]
 elif mode=='decomp':c=-ez;E=np.vstack([ex,ey]);h=[1,lock]
 elif mode=='stage1':c=-ey;E=np.vstack([ex,ez]);h=[1,lock] # paper normalization, differs from production
 else:c=-ez;E=np.vstack([ex,ey-lock*ez]);h=[1,0]
 return fit(c,A,np.zeros(2*n),E,np.array(h))

def fixture(shape,variant):
 m,d,s=shape;rng=np.random.default_rng(82919+100*m+10*d+s)
 a=rng.uniform(.3,8,(6,m+d+s));a[5,:m]=a[0,:m]*2;a[5,m:m+d]=a[0,m:m+d];a[5,m+d:]=a[0,m+d:]*.5
 if variant=='near':a[1]=a[0]*(1+1e-8)
 if variant=='zero':
  for start,size in [(0,m),(m,d),(m+d,s)]:
   if size>1:a[2,start]=0
 if variant=='zerocol':
  for start,size in [(0,m),(m,d),(m+d,s)]:
   if size>1:a[:,start+size-1]=0
 if variant=='units':a*=10.**np.linspace(-10,10,m+d+s)
 if variant=='tiny':a*=1e-10
 z=np.c_[np.arange(1,7),a,np.ones(6)]
 if variant=='duplicate':z=np.vstack([z,z[2]]);z[-1,0]=7
 if variant=='sample':z=np.vstack([z,z[0]]);z[-1,0]=7;z[-1,1:1+m]*=.01;z[-1,-1]=0;z[1,-1]=0
 if variant=='permuted':z=z[[5,2,0,4,1,3]]
 return z,['id']+[f'x{i+1}' for i in range(m)]+[f'z{i+1}' for i in range(d)]+[f'y{i+1}' for i in range(s)]+['selected']

def verify(path,shape,kind,method,rts,ort,leader,tol):
 m,d,s=shape;a=readbin(path/'data.bin');n=len(a);x=a[:,1:1+m];z=a[:,1+m:1+m+d];y=a[:,1+m+d:]
 scales=[np.where(v.max(0)==0,1,v.max(0)) for v in [x,z,y]];x,z,y=[v/scale for v,scale in zip([x,z,y],scales)]
 M={p.stem:readbin(p) for p in path.glob('*.bin')};R=M['efficiencies'];meta=M['meta'].ravel()
 assert np.array_equal(M['status'],np.isfinite(R));assert np.array_equal(meta[:3],np.sum(~np.isfinite(R),axis=0))
 rel=kind=='network' and method=='multiplicative' and rts=='crs' and ort=='input';nc=method=='noncooperative'
 counts=dict(independent_LPs=0,joint_LPs=0,max_score_gap=0.,max_decomp_gap=0.,unavailable_scores=int(np.isnan(R).sum()),decomp_fallback=0,follower_failed=0)
 gap=0
 if nc:
  from twostage_witness import check
  truths=check(path.name,a)
  counts['floating_follower_oracle_max_gap']=0.
  counts['exact_follower_rows']=0
 for k in range(n):
  b=np.array([stage(x,z,k,rts,ort),stage(z,y,k,rts,ort)]);counts['independent_LPs']+=2
  if rel:
   f=joint(x,z,y,k,'primary');assert f.success;counts['joint_LPs']+=1
   if np.isfinite(R[k,0]):
    gap=max(gap,abs(R[k,0]+f.fun));st=M['decomposition_status'][k,0];assert st in [1,2]
    if st==2:
     f=joint(x,z,y,k,'decomp',R[k,0]);counts['joint_LPs']+=1;assert f.success
     counts['max_decomp_gap']=max(counts['max_decomp_gap'],abs(R[k,1]+f.fun));assert abs(R[k,1]+f.fun)<3e-6
    else:counts['decomp_fallback']+=1
   else:assert M['decomposition_status'][k,0]==0 and np.isnan(R[k]).all()
  else:
   if kind=='stackelberg':
    bs=M['bench'][k];assert np.allclose(bs[np.isfinite(bs)],b[np.isfinite(bs)],atol=3e-6,rtol=1e-8)
   for j,(xx,yy) in enumerate([(x,z),(z,y)]):
    lam=M[f'lambda_stage{j+1}'][k]
    val=M['bench'][k,j] if kind=='stackelberg' else R[k,j+1]
    if not np.isfinite(val):assert np.isnan(lam).all();continue
    assert np.isfinite(lam).all() and lam.min()>-3e-6
    assert np.max(xx.T@lam-(val*xx[k] if ort=='input' else xx[k]))<3e-6
    assert np.max((yy[k] if ort=='input' else yy[k]/val)-yy.T@lam)<3e-6
    if rts=='vrs':assert abs(sum(lam)-1)<3e-6
   if nc:
    active=np.count_nonzero(z.max(0)>0)
    route=M['follower_route'][k,0]
    assert route in [0,1,2,3] and (route in [0,3] or route==(2 if active==1 else 1))
    assert (route==3)==(M['t121_quality'][k,9]==1)
    j=0 if leader=='stage1' else 1;lead=M['bench'][k,j];follow=R[k,2-j];st=M['follower_status'][k,0]
    if route==3:
     # The certified leader need not equal the old native benchmark bit-for-bit.
     # Independent benchmarks and linked outputs have separate verification.
     assert abs(R[k,j+1]-b[j])<3e-6
     lead=R[k,j+1]
    else:assert R[k,j+1]==lead or np.isnan(R[k,j+1]) and np.isnan(lead)
    if st==2:
     f=joint(x,z,y,k,leader,lead);assert f.success;counts['joint_LPs']+=1
     expected=-f.fun/lead if leader=='stage1' else -f.fun
     if route==3:
      # Retain the floating comparison as a diagnostic. The near fixture shows
      # both non-dyadic normalization loss and floating feasibility effects.
      counts['floating_follower_oracle_max_gap']=max(counts['floating_follower_oracle_max_gap'],abs(follow-expected))
      expected=float(truths[k][2-j])
      assert list(R[k])==[float(v) for v in truths[k]]
      counts['exact_follower_rows']+=1
     gap=max(gap,abs(follow-expected));assert follow<=b[1-j]+3e-6
    else:assert st in [0,1] and np.isnan(follow) and np.isnan(R[k,0]);counts['follower_failed']+=int(st==1)
   else:
    used=np.isfinite(R[k,1:]);gap=max(gap,np.max(abs(R[k,1:][used]-b[used]),initial=0))
  if np.isfinite(R[k]).all():
   expected=np.mean(R[k,1:]) if method=='additive' else np.prod(R[k,1:]);assert abs(expected-R[k,0])<3e-6
  if 'weights_input' in M:
   w=np.r_[M['weights_input'][k]*scales[0],M['weights_intermediate'][k]*scales[1],M['weights_output'][k]*scales[2]]
   if np.isnan(w).all():continue
   assert np.isfinite(w).all() and min(w)>-3e-6
   vx=x@w[:m];wz=z@w[m:m+d];uy=y@w[m+d:]
   assert max((wz-vx)/np.maximum(1,abs(vx)+abs(wz)))<3e-6 and max((uy-wz)/np.maximum(1,abs(uy)+abs(wz)))<3e-6
   if np.isfinite(R[k,1]):assert abs(wz[k]/vx[k]-R[k,1])<3e-6
   if np.isfinite(R[k,2]):assert abs(uy[k]/wz[k]-R[k,2])<3e-6
   if rel:assert abs(vx[k]-1)<3e-6
   elif leader=='stage1':assert abs(wz[k]-1)<3e-6
   else:assert abs(vx[k]-1)<3e-6
 assert gap<3e-6,gap;counts['max_score_gap']=gap
 if rel:assert meta[3]==counts['decomp_fallback']
 if nc:assert meta[4]==counts['follower_failed']
 return counts

R_SCRIPT=r'''
library(DJL);library(Benchmarking)
a=commandArgs(TRUE);dat=read.csv(a[1]);m=as.integer(a[2]);d=as.integer(a[3]);s=as.integer(a[4]);z=as.matrix(dat[dat$selected==1,]);x=z[,2:(1+m),drop=FALSE];w=z[,(2+m):(1+m+d),drop=FALSE];y=z[,(2+m+d):(1+m+d+s),drop=FALSE]
x=sweep(x,2,pmax(apply(x,2,max),1e-300),'/');w=sweep(w,2,pmax(apply(w,2,max),1e-300),'/');y=sweep(y,2,pmax(apply(y,2,max),1e-300),'/')
res=data.frame(id=z[,1]);for(r in c('crs','vrs'))for(o in c('in','out')) {e1=dea(x,w,RTS=r,ORIENTATION=o)$eff;e2=dea(w,y,RTS=r,ORIENTATION=o)$eff;if(o=='out'){e1=1/e1;e2=1/e2};res[[paste(r,o,'s1',sep='_')]]=e1;res[[paste(r,o,'s2',sep='_')]]=e2}
f=dm.network.dea(xdata.s1=x,zdata=w,ydata.s2=y,type='co');res$co_s1=f$eff.s1;res$co_s2=f$eff.s2
for(l in c('1st','2nd')) {f=dm.network.dea(xdata.s1=x,zdata=w,ydata.s2=y,type='nc',leader=l);res[[paste0(l,'_s1')]]=f$eff.s1;res[[paste0(l,'_s2')]]=f$eff.s2}
write.csv(res,a[5],row.names=FALSE)
'''

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--only');ap.add_argument('--rscript');a=ap.parse_args()
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir();h={}
 paths=[REPO/n for n in package_files(REPO)]
 for p in paths:
  if p.suffix in ['.ado','.do','.sthlp']:raw=p.read_bytes();(pkg/p.name).write_bytes(raw);h[p.name]=hashlib.sha256(raw).hexdigest()
 (out/'manifest.json').write_text(json.dumps(dict(source_sha256=h,stata=a.stata),indent=2)+'\n');rows=[]
 def run(name,body,check=None):
  if a.only and not re.search(a.only,name):return
  path=out/name;path.mkdir();head='version 16\nclear all\nset more off\nset processors 1\nset type double\nset linesize 160\n'
  for f in ['PLUS','PERSONAL','SITE','OLDPLACE']:(path/f).mkdir();head+=f'sysdir set {f} "{path/f}"\n'
  head+=f'adopath ++ "{pkg}"\n'+WRITER;(path/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n');r=dict(test=name,passed=False);t=time.monotonic()
  try:
   p=subprocess.run([a.stata,'-b','do','test.do'],cwd=path,capture_output=True,timeout=180)
   if p.returncode or not (path/'test.log').exists():(path/'process.json').write_text(json.dumps(dict(returncode=p.returncode,stdout=p.stdout.decode(errors='replace'),stderr=p.stderr.decode(errors='replace')),indent=2))
   log=(path/'test.log').read_text(errors='replace');assert p.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2000:]
   if check:r.update(check(path))
   r['passed']=True
  except Exception as e:r['error']=repr(e)
  r['seconds']=round(time.monotonic()-t,3);rows.append(r);(out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(r),flush=True)
 profiles=[('network',method,rts,ort,'stage1') for method,rts,ort in itertools.product(['multiplicative','independent','additive'],['crs','vrs'],['input','output'])]
 profiles += [('stackelberg','independent',rts,'input',l) for rts,l in itertools.product(['crs','vrs'],['stage1','stage2'])]+[('stackelberg','noncooperative','crs','input',l) for l in ['stage1','stage2']]
 specs=[]
 for sh,var in itertools.product([(1,1,1),(2,1,2),(1,2,1),(2,2,2)],['base','units','tiny','zero','zerocol','duplicate','sample','permuted','near']):
  dat,cols=fixture(sh,var);specs.append((f'{sh[0]}x{sh[1]}x{sh[2]}_{var}',sh,dat,cols))
 paper=np.genfromtxt(io.StringIO(PUBLISHED_DATA),delimiter=',',skip_header=1)
 specs.append(('published_base',(2,2,2),np.c_[paper,np.ones(len(paper))],['id','x1','x2','z1','z2','y1','y2','selected']))
 for key,sh,dat,cols in specs:
  f=out/('data_'+key+'.csv');np.savetxt(f,dat,delimiter=',',header=','.join(cols),comments='',fmt='%.17g');m,d,s=sh;xs=' '.join(f'x{i+1}' for i in range(m));zs=' '.join(f'z{i+1}' for i in range(d));ys=' '.join(f'y{i+1}' for i in range(s))
  for profile,tol in itertools.product(profiles,['1e-7','1e-10']):
   kind,method,rts,ort,leader=profile;rel=kind=='network' and method=='multiplicative' and rts=='crs' and ort=='input';nc=method=='noncooperative'
   cmd=f'dea_{kind} if selected, inputs({xs}) intermediates({zs}) outputs({ys}) method({method}) rts({rts}) '+(f'ort({ort})' if kind=='network' else f'leader({leader})')+f' dmu(id) gen(res) tol({tol}) matrix detail'
   body=f'import delimited using "{f}",clear varnames(1) asdouble\n'+cmd+'\n'
   body+='matrix meta=(r(n_failed_overall),r(n_failed_stage1),r(n_failed_stage2),'+('r(n_decomp_fallback),.' if kind=='network' else '.,r(n_failed_follower_lp)')+')\nmata: savebin("meta.bin",st_matrix("meta"))\n'
   for name in MATS:
    if name.startswith('weights_') and not(rel or nc):continue
    if name=='decomposition_status' and not rel:continue
    if name in ['follower_status','follower_route','t121_quality'] and not nc:continue
    if name.startswith('t122_') and (rel or nc):continue
    if name.startswith('lambda_') and rel:continue
    body+=f'mata: savebin("{name}.bin",st_matrix("r({name})"))\n'
   expected='relational_crs' if rel else 'noncooperative_crs' if nc else 'independent_ordering' if kind=='stackelberg' else 'independent_mean' if method=='additive' else 'independent_product'
   body+=f'assert "`r(formulation)\'"=="{expected}"\n'
   body+='matrix expected=r(efficiencies)\nforeach v of varlist res_* {\n assert missing(`v\') if !selected\n}\nkeep if selected\n'
   if kind=='stackelberg':body+='mata: savebin("bench.bin",st_data(.,tokens("res_indep1 res_indep2")))\n'+f'mata: assert(all(st_data(.,"res_leader"):==st_matrix("expected")[,{2 if leader=="stage1" else 3}]))\n'
   body+=f'mata: savebin("data.bin",st_data(.,tokens("id {xs} {zs} {ys}")))\n'
   name=key+'_'+'_'.join(profile)+'_'+tol
   run(name,body,lambda p,sh=sh,profile=profile,tol=tol:verify(p,sh,*profile,float(tol)))
 run('smallmiddle_available', 'input double(id x z y)\n1 1e8 1 1\n2 1 1 1\nend\ndea_network x=y,intermediates(z) gen(res) matrix detail\nassert !missing(res_overall,res_stage1,res_stage2)\nassert abs(res_overall-1e-8)<1e-12 in 1\nassert abs(res_stage1-1e-8)<1e-12 in 1\nassert abs(res_stage2-1)<1e-4 in 1\nassert abs(res_overall-res_stage1*res_stage2)<1e-12\n')
 data='input double(id x z y)\n1 1 1 1\n2 2 3 2\n3 3 2 3\nend\n'
 for kind in ['network','stackelberg']:
  cmd=f'dea_{kind},inputs(x) intermediates(z) outputs(y) gen(res)'
  for opt in ['tol(0)','tol(-1)','tol(.)','tol(.a)']:
   run(kind+'_invalid_'+re.sub(r'\W','_',opt),data+'gen str12 res_overall="keep"\n'+unchanged(cmd+' replacevars '+opt,198))
  for var,change in itertools.product(['x','z','y'],[-1,0]):run(f'{kind}_domain_{var}_{change}',data+f'replace {var}={change} in 1\ngen str12 res_overall="keep"\n'+unchanged(cmd+' replacevars',459))
  run(kind+'_empty',data+'gen str12 res_overall="keep"\n'+unchanged(cmd.replace(',', ' if 0,',1)+' replacevars',2000))
  run(kind+'_missing_id',data+'replace id=. in 1\n'+cmd+' dmu(id)\nassert r(n)==2\nassert missing(res_overall) in 1\n')
  run(kind+'_singleton',data+'keep in 1\n'+cmd+' matrix detail\nassert abs(res_overall-1)<3e-7\n')
  run(kind+'_workspace_reload',data+'matrix user=(9,8)\nlocal before : all matrices\n'+cmd+' matrix\nlocal after : all matrices\nassert "`before\'"=="`after\'"\n'+f'discard\n'+cmd+' replacevars matrix\nlocal after : all matrices\nassert "`before\'"=="`after\'"\n')
 run('unsupported_nc_vrs',data+'gen str12 res_overall="keep"\n'+unchanged('dea_stackelberg, inputs(x) intermediates(z) outputs(y) method(noncooperative) rts(vrs) gen(res) replacevars',198))
 for method,rts,leader in [('independent','crs','stage1'),('independent','vrs','stage2'),('noncooperative','crs','stage1'),('noncooperative','crs','stage2')]:
  for label,q in QUALIFIERS.items():
   cmd=f'dea_stackelberg,inputs(x1 x2) intermediates(z1 z2) outputs(y1 y2) dmu(sid) gen(res) method({method}) rts({rts}) leader({leader}) matrix detail debug'
   outs='res_leader res_follower res_overall res_indep1 res_indep2'
   body=io_fixture('network')+f'gen byte selected=0\nreplace selected=1 {q}\ntempfile expected\npreserve\nkeep if selected\n'+cmd+'\n'+CAPTURE+'keep obsid '+outs+'\nsave `expected\'\nrestore\n'+cmd.replace(',',f' {q},',1)+'\n'+VERIFY+'foreach v of varlist '+outs+' {\nassert missing(`v\') if !selected\n}\nkeep if selected\ncf obsid '+outs+' using `expected\',all\n'
   run('stack_sample_'+method+'_'+rts+'_'+leader+'_'+label,body)
 for leader in ['stage1','stage2']:
  cmd='dea_stackelberg, inputs(x) intermediates(z) outputs(y) method(noncooperative) leader('+leader+') gen(res)'
  run('nc_protect_'+leader,data+'gen str12 res_overall="keep"\n'+unchanged(cmd,110))
  run('nc_model_collision_'+leader,data+'rename x res_leader\n'+unchanged(cmd.replace('inputs(x)','inputs(res_leader)')+' replacevars',110))
  run('nc_replace_'+leader,data+cmd+'\nrename res_overall expected\ndrop res_*\ngen str12 res_overall="old"\n'+cmd+' replacevars matrix detail\nassert res_overall==expected\n')
 if a.rscript:
  (out/'compare.R').write_text(R_SCRIPT)
  def rcheck(path):
   gaps=[];checked=0;unavailable=0
   for key,sh,dat,cols in specs:
    if not key.endswith('_base'):continue
    dest=path/(key+'.csv');p=subprocess.run([a.rscript,str(out/'compare.R'),str(out/('data_'+key+'.csv')),*map(str,sh),str(dest)],capture_output=True,text=True,timeout=180);(path/(key+'.console.txt')).write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr
    ref=np.genfromtxt(dest,delimiter=',',names=True)
    for profile in profiles:
     kind,method,rts,ort,leader=profile;d0=out/(key+'_'+'_'.join(profile)+'_1e-10')
     if not d0.exists():continue
     got=readbin(d0/'efficiencies.bin');unavailable+=int(np.isnan(got).sum());tag=rts+'_'+('in' if ort=='input' else 'out')
     if method=='noncooperative':tag='1st' if leader=='stage1' else '2nd';want=np.c_[ref[tag+'_s1'],ref[tag+'_s2']]
     elif kind=='network' and method=='multiplicative' and rts=='crs' and ort=='input':
      # DJL chooses an arbitrary centralized optimum; only the system score is unique.
      want=ref['co_s1']*ref['co_s2'];good=np.isfinite(got[:,0]);gap=max(abs(got[good,0]-want[good]),default=0);assert gap<3e-6;gaps.append(gap);checked+=sum(good);continue
     else:want=np.c_[ref[tag+'_s1'],ref[tag+'_s2']]
     good=np.isfinite(got[:,1:]);gap=max(abs(got[:,1:][good]-want[good]),default=0);assert gap<3e-6;gaps.append(gap);checked+=int(good.sum())
   return dict(R_score_comparisons=int(checked),R_max_gap=max(gaps,default=0),Stata_unavailable_component_scores=unavailable)
  run('R_original_packages','',rcheck)
 summary=dict(total=len(rows),passed=sum(r['passed'] for r in rows),failed=sum(not r['passed'] for r in rows));(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary);return bool(summary['failed'])
if __name__=='__main__':raise SystemExit(main())
