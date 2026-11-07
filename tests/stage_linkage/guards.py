from runlib import *
from test_recoverylog import SNAP
import shutil,difflib
body=f'run "{r}/candidate/ado/_dea_simplex.ado"\nmata:\n'+r'''
void direct_guards() {
 struct RSM_Prob scalar P,Q
 struct RSM_Result scalar R,S
 real matrix XP,B
 real rowvector detail,res
 real scalar ok,i
 P.m_con=2;P.n_orig=2;P.n_var=2;P.A=I(2);P.b=(0\1)
 R.basiscols=(1\2);R.status=0;R.x=P.b;R.objval=0;R.iters=0
 ok=_dea_link_basis(P,R,XP,B,detail);assert(ok==1 & B[1,1]==0 & B[1,2]==0)
 savebin("zero_bounds.bin",B)
 P.b=(-2^-30\1);ok=_dea_link_basis(P,R,XP,B,detail);assert(ok==1 & B[1,2]<0)
 savebin("negative_bounds.bin",B)
 P.b=(2^-30\1);ok=_dea_link_basis(P,R,XP,B,detail);assert(ok==1 & B[1,1]>0)
 savebin("positive_bounds.bin",B)
 for(i=1;i<=10;i++) {
  Q=P;S=R
  if(i==1) Q.m_con=33
  if(i==2) Q.n_orig=4093
  if(i==3) Q.m_con=1.5
  if(i==4) Q.A=J(1,2,0)
  if(i==5) Q.b=J(1,1,0)
  if(i==6) S.basiscols=(1\1)
  if(i==7) S.basiscols=(1\3)
  if(i==8) S.basiscols=(1\.)
  if(i==9) Q.A=J(2,2,1)
  if(i==10) Q.A=I(2)*2^500
  assert(_dea_link_basis(Q,S,XP,B,detail)==0)
 }
 Q=P;Q.m_con=33;assert(hasmissing(_dea_link_residual(Q,R.x)))
 Q=P;Q.n_orig=4093;assert(hasmissing(_dea_link_residual(Q,R.x)))
 Q=P;Q.b=J(1,1,0);assert(hasmissing(_dea_link_residual(Q,R.x)))
 Q=P;Q.A=J(1,2,0);assert(hasmissing(_dea_link_residual(Q,R.x)))
 Q=P;assert(hasmissing(_dea_link_residual(Q,J(3,1,1))))
 Q=P;assert(hasmissing(_dea_link_residual(Q,(.\1))))
 savebin("direct_count.bin",19)
}
direct_guards()
end
'''
run('guards',body,r/'candidate/ado')
faults={
 'retry_failure':('T=_rsm_two_phase(trial,linktol,"");meta[11]=T.status;','T=_rsm_two_phase(trial,linktol,"");T.status=3;meta[11]=T.status;',6),
 'new_basis_unavailable':('if(!_dea_link_basis(P,T,XP,bounds,detail))','if(1)',7),
 'excessive_change':('meta[22]=safe','safe=oldrad*.9\n meta[22]=safe',8),
 'bad_primal':('x=T.x;x[I.idx_theta]=safe','x=T.x;x[I.idx_theta]=safe;x[2]=x[2]+1',9),
 'stage2_failure':('res2 = _rsm_two_phase(prob2, linktol, trace)','res2 = _rsm_two_phase(prob2, linktol, trace)\n    res2.status=3',None),
 'stage2_bad_certificate':('res2 = _rsm_two_phase(prob2, linktol, trace)','res2 = _rsm_two_phase(prob2, linktol, trace)\n    if(linkmeta[15]==1) {res2.x[1]=res2.x[1]+1;res2.basiscols=J(0,1,.);}',-1)}
z=next(z for z in j(r/'GRID_DATA.json') if z['kind']=='near_parallel');rows=[]
for name,(anchor,replacement,reason) in faults.items():
 pkg=r/'faults'/name/'ado';shutil.copytree(r/'candidate/ado',pkg);p=pkg/'_dea_simplex.ado';old=p.read_text();assert old.count(anchor)==1;s=old.replace(anchor,replacement,1);p.write_text(s);(pkg.parent/'fault.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),s.splitlines(True))))
 body=f'run "{pkg}/_dea_simplex.ado"\nmata: _dea_link_reset(2,1)\ninput double id double x1 double x2 double y1 double y2\n'
 for k,row in enumerate(z['data'],1):body+=str(k)+' '+' '.join(format(x,'.17g') for x in row)+'\n'
 body+='end\ndea x1 x2=y1 y2,dmu(id) engine(simplex) rts(crs) ort(output) stage(2) reference gen(eff)\n'+SNAP.replace('PREFIX',name)+f'mata: savebin("{name}_link.bin",_dea_link_read())\n'
 d=run('fault_'+name,body,pkg);M=readbin(d/f'{name}_link.bin');score=readbin(d/f'{name}_scores.bin')
 if reason is not None and reason>=0:
  assert (M[:,14]==0).all() and (M[[0,2,3,4,5],23]==reason).all()
  for end in ['scores','m_lambda','m_slack_in','m_slack_out','m_stage_status','m_numerical_diagnostics']:
   assert h(d/f'{name}_{end}.bin')==j(r/'REFERENCE.json')['grid_baseline'][f'near_parallel_crs_output_s2_{end}.bin']
 elif reason is None:
  assert h(d/f'{name}_scores.bin')==h(r/'checks/grid_repair/near_parallel_crs_output_s2_scores.bin')
  assert np.isnan(readbin(d/f'{name}_m_lambda.bin')).all() and np.isnan(readbin(d/f'{name}_m_slack_in.bin')).all()
  assert (readbin(d/f'{name}_m_stage_status.bin')[:,1]==3).all()
 else:
  assert h(d/f'{name}_scores.bin')==h(r/'checks/grid_repair/near_parallel_crs_output_s2_scores.bin')
  L=readbin(d/f'{name}_m_lambda.bin');S=readbin(d/f'{name}_m_stage_status.bin')
  assert np.isnan(L[[0,2,3,4,5]]).all() and np.isfinite(L[1]).all()
  assert (S[[0,2,3,4,5],1]==3).all() and S[1,1]==0 and (M[[0,2,3,4,5],28]<=0).all()
 rows.append(dict(case=name,passed=True,reason=reason))
dump('GUARD_REVIEW.json',dict(passed=True,direct_checks=19,fault_commands=6,fault_DMUs=36,rows=rows,production_changed=False))
print('Guard/fault checks passed.',flush=True)
