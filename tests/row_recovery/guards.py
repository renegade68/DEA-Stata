from runlib import *
import shutil,difflib
body=f'run "{r}/candidate/ado/_dea_simplex.ado"\nmata:\n'+r'''
void guard() {
 struct RSM_Prob scalar P,Q
 struct RSM_Result scalar R,S
 real colvector kr,kc,L
 real rowvector d,m
 real scalar ok
 _rsm_std_form(P,J(0,2,0),J(0,1,0),(1,0\0,1\1,1),(1\2\3),(1,1))
 assert(_rsm_rowmap(P,kr,kc,d))
 R.status=0;R.x=(1\2);R.objval=3;R.iters=0;R.basiscols=(1\2)
 Q=P;S=R;assert(_rsm_row_accept(Q,S,kr,1e-10));savebin("good.bin",S.x)
 Q=P;S=R;assert(!_rsm_row_accept(Q,S,(2\3),1e-10))
 Q=P;S=R;Q.b[3]=3.0000000000000004;assert(!_rsm_row_accept(Q,S,kr,1e-10))
 Q=P;S=R;Q.A[3,1]=1.0000000000000002;assert(!_rsm_row_accept(Q,S,kr,1e-10))
 Q=P;S=R;S.basiscols=(1\1);assert(!_rsm_row_accept(Q,S,kr,1e-10))
 Q=P;S=R;S.basiscols=1;assert(!_rsm_row_accept(Q,S,kr,1e-10))
 Q=P;S=R;S.status=3;assert(!_rsm_row_accept(Q,S,kr,1e-10))
 Q=P;S=R;S.x[1]=.;assert(!_rsm_row_accept(Q,S,kr,1e-10))
 Q=P;S=R;S.x=(0\0);assert(_rsm_row_accept(Q,S,kr,1e-10));assert(all(S.x:==R.x))
 // Empty bases: preserve alternatives; correct bound/objective; reject false triviality.
 _rsm_std_form(P,J(0,2,0),J(0,1,0),J(0,2,0),J(0,1,0),(1,0))
 L=(-2\3);R.status=0;R.x=(-2\7);R.objval=-2;R.basiscols=J(0,1,.)
 S=R;m=_rsm_success_check(P,S,R.basiscols,1e-10,L);assert(m[1]==1 & all(S.x:==R.x))
 S=R;S.x=(3\1);S.objval=100;m=_rsm_success_check(P,S,R.basiscols,1e-10,L);assert(m[1]==2 & all(S.x:==L) & S.objval==-2)
 Q=P;Q.c[1]=-1;S=R;m=_rsm_success_check(Q,S,R.basiscols,1e-10,L);assert(m[1]==0)
 _rsm_std_form(Q,J(0,2,0),J(0,1,0),(0,0),1e-200,(1,0));S=R;m=_rsm_success_check(Q,S,R.basiscols,1e-10,L);assert(m[1]==-1)
 _rsm_std_form(Q,J(0,2,0),J(0,1,0),(1e-200,0),0,(1,0));S=R;m=_rsm_success_check(Q,S,R.basiscols,1e-10,L);assert(m[1]==-1)
 S=R;S.x=J(2,1,.);m=_rsm_success_check(P,S,R.basiscols,1e-10,L);assert(m[1]==-1)
 // Full original nonbasic column must satisfy the dual inequality.
 _rsm_std_form(P,J(0,3,0),J(0,1,0),(1,0,1\0,1,1\1,1,2),(1\2\3),(1,1,0))
 R.x=(1\2\0);R.objval=3;R.basiscols=(1\2);assert(_rsm_rowmap(P,kr,kc,d));assert(!_rsm_row_accept(P,R,kr,1e-10))
 savebin("guards.bin",16)
}
guard()
end
'''
run('row_guards',body)
z=next(z for z in j(r/'GRID_FIXTURE.json') if z['id']=='g0162');E=np.array(z['E']);b=np.array(z['b']).reshape(-1,1);c=np.array(z['c']).reshape(1,-1)
source=(r/'candidate/ado/_dea_simplex.ado').read_text();records=[]
for fault in ['map_budget','wrong_rows','missing_x','wrong_basis','full_check','iteration_limit']:
 pkg=r/('fault_'+fault)/'ado';shutil.copytree(r/'candidate/ado',pkg);s=source
 a=s.index('struct RSM_Result scalar _rsm_row_retry(');b1=s.index('\n}',a)+2;part=s[a:b1]
 if fault=='map_budget':part=part.replace('_rsm_rowmap(original,kr,kc,detail)','_rsm_rowmap(original,kr,kc,detail,0)')
 if fault=='wrong_rows':part=part.replace(' if(_rsm_row_accept',' kr=J(0,1,.)\n if(_rsm_row_accept')
 if fault=='missing_x':part=part.replace(' if(_rsm_row_accept',' result.x=J(original.n_orig,1,.)\n if(_rsm_row_accept')
 if fault=='wrong_basis':part=part.replace(' if(_rsm_row_accept',' result.basiscols=J(0,1,.)\n if(_rsm_row_accept')
 if fault=='full_check':part=part.replace('if(_rsm_row_accept(original,result,kr,tol))','if(0)')
 if fault=='iteration_limit':
  part=part.replace(' real colvector kr,kc,allcols',' external real scalar _t35_faultrow\n real colvector kr,kc,allcols').replace(' result=_rsm_two_phase_core',' _t35_faultrow=1\n result=_rsm_two_phase_core').replace(' if(_rsm_row_accept',' _t35_faultrow=0\n if(_rsm_row_accept')
 s=s[:a]+part+s[b1:]
 if fault=='iteration_limit':
  for name in ['_rsm_two_phase_once','_rsm_fp_once','_rsm_precision_once']:
   a=s.index('struct RSM_Result scalar '+name+'(');b1=s.index('\n}',a)+2;part=s[a:b1];k=part.index('{')+1;part=part[:k]+'\n external real scalar _t35_faultrow\n'+part[k:]
   anchor='    max_iter = max((10000, 20 * prob.m_con))';assert part.count(anchor)==1;part=part.replace(anchor,anchor+'\n    if(_t35_faultrow) max_iter=1');s=s[:a]+part+s[b1:]
 assert s!=source;(pkg/'_dea_simplex.ado').write_text(s)
 script=f'run "{pkg}/_dea_simplex.ado"\nmata: _t35_faultrow=0\n'+DIRECT+f'mata: trial("fault",{mata(E)},{mata(b)},{mata(c)},1e-10,1)\n'
 out=run('row_fault_'+fault,script,pkg);res=readbin(out/'fault_result.bin').ravel();audit=readbin(out/'fault_audit.bin');old=(3,0,j(r/'REFERENCE.json')['g0162_previous_iterations'])
 assert res[0]==3 and res[2]>=old[2] and audit.shape==(1,12) and audit[0,6]==res[2]
 if fault=='map_budget':assert res[2]==old[2]
 records.append(dict(fault=fault,status=int(res[0]),iterations=int(res[2]),previous_iterations=int(old[2]),source_sha256={p.name:h(p) for p in pkg.iterdir() if p.is_file()}))
dump('GUARDS_REVIEW.json',dict(passed=True,direct_controls=16,fault_paths=6,false_successes=0,one_top_level_audit=True,rows=records))
