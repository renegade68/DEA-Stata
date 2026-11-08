#!/usr/bin/env python3
"""Portable public recovery/certificate checks using private observation copies."""
from pathlib import Path
from fractions import Fraction as F
import argparse,csv,hashlib,itertools,json,shutil,subprocess,time
import numpy as np
from test_core_retries import WRITER,readbin
REPO=Path(__file__).resolve().parents[1]
SAVE='''
void _c63_record(string scalar tag,real rowvector c,real matrix G,
 real colvector h,real scalar ne,struct _C58_Cert scalar fit) {
 savebin(tag+"_c.bin",c);savebin(tag+"_G.bin",G);savebin(tag+"_h.bin",h)
 savebin(tag+"_meta.bin",(fit.status,fit.reason,ne,fit.ops,fit.maxlimbs))
 savebin(tag+"_num.bin",fit.numerator);savebin(tag+"_den.bin",fit.denominator)
 savebin(tag+"_basis.bin",fit.basis);savebin(tag+"_dual.bin",fit.dual_numerator)
 savebin(tag+"_basis_den.bin",fit.basis_denominator);savebin(tag+"_exponent.bin",fit.row_exponent)
 savebin(tag+"_face.bin",fit.face);savebin(tag+"_trace.bin",fit.trace62)
 savebin(tag+"_ray.bin",fit.ray);savebin(tag+"_ray_num.bin",fit.ray_num);savebin(tag+"_ray_den.bin",fit.ray_den)
}
void _c70_cost_record(string scalar tag,transmorphic scalar C,real scalar ce,real scalar n) {
 real matrix packed
 real rowvector a
 real scalar j
 packed=J(n,129,0)
 for(j=1;j<=n;j++) {a=asarray(C,(1,j));packed[j,1..cols(a)]=a;}
 savebin(tag+"_exact_cost.bin",packed);savebin(tag+"_exact_exponent.bin",ce)
}
'''
def integer(row):
 row=list(row);assert row[0] in [-1,0,1]
 assert all(v==int(v) and 0<=v<2**26 for v in row[1:])
 return int(row[0])*sum(int(v)*2**(26*i) for i,v in enumerate(row[1:]))
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def check(prefix):
 read=lambda k:readbin(Path(str(prefix)+'_'+k+'.bin'))
 meta=read('meta').ravel();ne=int(meta[2]);G=[[F(float(v)) for v in row] for row in read('G')]
 h=[F(float(v)) for v in read('h').ravel()];c=[F(float(v)) for v in read('c').ravel()]
 if Path(str(prefix)+'_exact_cost.bin').exists():c=[F(integer(row))*F(2)**int(read('exact_exponent')[0,0]) for row in read('exact_cost')]
 proof=ray=0
 if meta[0]==1:
  D=integer(read('den').ravel());assert D>0
  x=[F(integer(row),D) for row in read('num')]
  assert all(dot(a,x)==b for a,b in zip(G[:ne],h[:ne])) and all(dot(a,x)<=b for a,b in zip(G[ne:],h[ne:]))
  ids=[int(v)-1 for v in read('basis').ravel()];ex=list(map(int,read('exponent').ravel()));Db=integer(read('basis_den').ravel());dual=[F(0)]*len(G)
  factor_path=Path(str(prefix)+'_original_objective_factor.bin');factor=integer(read('original_objective_factor').ravel()) if factor_path.exists() else 1
  assert factor>0
  for i,row in enumerate(read('dual')):dual[ids[i]]=factor*F(integer(row),Db)*F(2)**(ex[-1]-ex[i])
  assert all(v>=0 for v in dual[ne:])
  assert [sum((a[j]*v for a,v in zip(G,dual)),F(0)) for j in range(len(c))]==c
  assert dot(c,x)==dot(h,dual)
  assert [int(v)-1 for v in read('face').ravel()]==[i for i in ids if i>=ne and dual[i]>0]
  proof=1
 z=read('trace').ravel()
 assert z[21]==sum(z[16:21]) and z[12]==sum(z[[6,8,10,11,24]]) and z[12]<=16
 if z[15]==1:
  D=integer(read('ray_den').ravel());assert D>0
  x=[F(integer(row),D) for row in read('ray_num')];d=[F(integer(row)) for row in read('ray')]
  assert all(dot(a,x)==b for a,b in zip(G[:ne],h[:ne])) and all(dot(a,x)<=b for a,b in zip(G[ne:],h[ne:]))
  assert all(dot(a,d)==0 for a in G[:ne]) and all(dot(a,d)<=0 for a in G[ne:]) and dot(c,d)>0
  assert meta[0]!=1;ray=1
 return proof,ray
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir()
 source={}
 for f in sorted((REPO/'ado').iterdir()):
  if f.is_file() and f.suffix in ['.ado','.do','.sthlp']:
   shutil.copy2(f,pkg/f.name);source['ado/'+f.name]=hashlib.sha256(f.read_bytes()).hexdigest()
 audit=out/'audit_ado';shutil.copytree(pkg,audit)
 f=audit/'_dea_crosseff_simplex.ado';s=f.read_text()
 for call,stage,obj,G,h,ne,fit in [
  ('_c58_fit(c*T,E*T,b,U*T,ub,tol,primary,G,h)','p','-c*T','G','h','rows(E)','primary'),
  ('_c70_fit(c2,F,fb,U*T,ub,tol,exactC,exact_exponent,objective_diag,secondary,G2,h2)','s','-c2','G2','h2','rows(F)','secondary')]:
  assert s.count(' '+call)==1
  tag='st_global("CROSSEFF_AUDIT_DIR")+"/dmu"+strofreal(index,"%02.0f")+"_'+stage+'"'
  s=s.replace(' '+call,' '+call+f'\n _c63_record({tag},{obj},{G},{h},{ne},{fit})',1)
  if stage=='s':s=s.replace(' '+call,' '+call+f'\n _c70_cost_record({tag},exactC,exact_exponent,cols(c2))',1)
 s=s.rsplit('\nend',1)[0]+SAVE+'\nend\n';f.write_text(s)
 for name,needle in [('_dea_crosseff_simplex.ado',' n=cols(c);T=I(n);x=J(1,n,.)'),('_dea_crosseff_cert.ado',' fit=_c58b_solve(-c,E,b,U,ub'),('_dea_crosseff_exactobjective.ado',' fit=_c70b_solve(-c,E,b,U,ub')]:
  f=audit/name;s=f.read_text();assert needle in s;s=s.replace(needle,' assert(tol==st_numscalar("CROSSEFF_REQUESTED_TOL"))\n'+needle,1);f.write_text(s)
 code='version 16\nclear all\nset more off\nset processors 1\nset type double\nset linesize 180\n'
 for k in ['PLUS','PERSONAL','SITE','OLDPLACE']:(out/k).mkdir();code+=f'sysdir set {k} "{out/k}"\n'
 code+=f'adopath ++ "{audit}"\n'+WRITER
 basic=[[1,1,1],[2,2,4],[3,4,2]];raydata=[[1,1,1],[2,2,2]]
 near=json.loads((REPO/'tests/fixtures/crosseff_recovery.json').read_text())['near'];cases=[]
 specs=[('native_cold',basic,1,1,'ipm','crs','input','arbitrary','1e-7','legacy')]
 for tol,rts,ori,method in itertools.product(['1e-4','1e-7','1e-12'],['crs','vrs'],['input','output'],['arbitrary','benevolent','aggressive']):
  specs.append((f'near_{tol}_{rts}_{ori}_{method}',near,2,2,'simplex',rts,ori,method,tol,'legacy'))
 specs += [('ray',raydata,1,1,'simplex','vrs','input','benevolent','1e-7','legacy'),('limzhu',near,2,2,'simplex','vrs','input','benevolent','1e-10','limzhu')]
 for name,data,m,n,engine,rts,ori,method,tol,rule in specs:
  d=out/name;d.mkdir();names=['id']+[f'x{i+1}' for i in range(m)]+[f'y{i+1}' for i in range(n)]
  with (d/'fixture.csv').open('w') as fh:w=csv.writer(fh);w.writerow(names);w.writerows(data)
  xs=' '.join(names[1:1+m]);ys=' '.join(names[1+m:])
  code+=f'import delimited using "{d/"fixture.csv"}",clear varnames(1) asdouble\nglobal CROSSEFF_AUDIT_DIR "{d}"\nscalar CROSSEFF_REQUESTED_TOL={tol}\n'
  code+=f'dea_crosseff {xs}={ys},dmu(id) rts({rts}) ort({ori}) method({method}) vrsc({rule}) eng({engine}) tol({tol}) gen(res) matrix\n'
  code+=f'assert r(tolerance)=={tol}\nassert r(diagnostic_version)==6201\nassert "`r(actual_engine)\'"=="{engine}"\n'
  if name=='native_cold':code+='mata: assert(findexternal("_rsm_solve_lp()")==NULL)\n'
  if name=='ray':code+='assert r(n_exact_unbounded)>0\nassert r(n_failed_secondary)>0\n'
  for key in ['primary_recovery_diagnostics','secondary_recovery_diagnostics']:
   code+=f'mata: savebin("{d/key}.bin",st_matrix("r({key})"))\n'
  code+='matrix cc=(r(n_recovery_requested),r(n_recovery_finite),r(n_recovery_unknown),r(n_recovery_skipped),r(n_recovery_equality),r(n_recovery_auxiliary),r(n_exact_unbounded),r(recovery_integer_ops),r(recovery_steps),r(certificate_integer_ops))\n'
  code+=f'mata: savebin("{d}/counts.bin",st_matrix("cc"))\n'
  code+=f'dea_crosseff {xs}={ys},dmu(id) rts({rts}) ort({ori}) eng({engine}) tol({tol}) gen(arbitrary)\nassert arbitrary_self==res_self\n'
  cases.append(dict(test=name,m=m,s=n,engine=engine))
 code+='di "CROSSEFF_RECOVERY_PASS"\nexit,clear\n';(out/'test.do').write_text(code);start=time.monotonic()
 q=subprocess.run([args.stata,'-b','do','test.do'],cwd=out,capture_output=True,timeout=1800)
 log=(out/'test.log').read_text(errors='replace');assert q.returncode==0 and 'CROSSEFF_RECOVERY_PASS' in log.splitlines(),'\n'.join(log.splitlines()[-25:])
 rows=[];proofs=rays=faces=0
 for case in cases:
  d=out/case['test'];A=np.vstack([readbin(d/'primary_recovery_diagnostics.bin'),readbin(d/'secondary_recovery_diagnostics.bin')])
  assert A.shape[1]==31 and np.all(A[:,19]==np.sum(A[:,14:19],axis=1)) and np.all(A[:,12]==np.sum(A[:,[6,8,10,11]],axis=1)) and np.all(A[:,12]<=16)
  expected=[sum(A[:,0]),sum(A[:,1]==1),sum(A[:,1]==3),sum(A[:,1]==-1),sum(A[:,7]),sum(A[:,9]),sum(A[:,13]==1),sum(A[:,19]-A[:,14]),sum(A[:,12]-A[:,6]),sum(A[:,19])]
  np.testing.assert_array_equal(expected,readbin(d/'counts.bin').ravel())
  for f in d.glob('dmu*_meta.bin'):
   pref=Path(str(f).removesuffix('_meta.bin'));a,b=check(pref);proofs+=a;rays+=b
   if pref.name.endswith('_s'):
    pp=Path(str(pref).removesuffix('_s')+'_p');read=lambda pref,k:readbin(Path(str(pref)+'_'+k+'.bin'))
    ne=int(read(pp,'meta')[0,2]);face=read(pp,'face').ravel().astype(int)-1;N=int(read(pref,'meta')[0,2])
    np.testing.assert_array_equal(read(pref,'G')[:N],np.vstack([read(pp,'G')[:ne],read(pp,'G')[face]]));np.testing.assert_array_equal(read(pref,'h')[:N],np.vstack([read(pp,'h')[:ne],read(pp,'h')[face]]));faces+=1
  rows.append(dict(**case,passed=True))
 for name,sha in source.items():assert hashlib.sha256((REPO/name).read_bytes()).hexdigest()==sha
 (out/'manifest.json').write_text(json.dumps(dict(source_sha256=source),indent=2)+'\n');(out/'results.json').write_text(json.dumps(rows,indent=2)+'\n')
 summary=dict(passed=True,cases=len(rows),failed=0,exact_finite_certificates=proofs,exact_rays=rays,original_face_links=faces,source_files_preserved=True,actual_tolerance_return_and_delivery_asserted=True,seconds=time.monotonic()-start)
 (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
