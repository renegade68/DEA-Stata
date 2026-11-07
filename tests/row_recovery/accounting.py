from runlib import *
import shutil,difflib
pkg=r/'accounted/ado';shutil.copytree(r/'candidate/ado',pkg);p=pkg/'_dea_simplex.ado';old=p.read_text();s=old
for name in ['_rsm_two_phase_once','_rsm_fp_once','_rsm_precision_once']:
 a=s.index('struct RSM_Result scalar '+name+'(');b=s.index('\n}',a)+2;part=s[a:b];start=part.index('{')+1
 part=part[:start]+'\n external real scalar _t35_work\n'+part[start:]
 part=part.replace('return(result)','_t35_work=_t35_work+result.iters;return(result)');s=s[:a]+part+s[b:]
s=s.replace('struct RSM_Result scalar _rsm_two_phase_core(', 'struct RSM_Result scalar _t35_core_impl(',1)
s=s.replace('struct RSM_Result scalar _rsm_row_retry(', 'struct RSM_Result scalar _t35_row_impl(',1)
helper=r'''
struct RSM_Result scalar _rsm_two_phase_core(struct RSM_Prob scalar P,
 real scalar tol,string scalar trace,real scalar allow_duplicates,| real scalar check_success,real scalar allow_rows) {
 external real scalar _t35_depth,_t35_maxdepth
 external real matrix _t35_calls
 struct RSM_Result scalar R
 real scalar m
 if(args()<5) check_success=1
 if(args()<6) allow_rows=1
 m=P.m_con;_t35_depth++;_t35_maxdepth=max((_t35_maxdepth,_t35_depth))
 R=_t35_core_impl(P,tol,trace,allow_duplicates,check_success,allow_rows)
 _t35_calls=_t35_calls\(_t35_depth,allow_duplicates,check_success,allow_rows,m,R.status,R.iters)
 _t35_depth--
 return(R)
}
struct RSM_Result scalar _rsm_row_retry(struct RSM_Prob scalar P,
 real scalar tol,string scalar trace,real scalar check_success) {
 external real scalar _t35_rows,_t35_rowdepth,_t35_maxrowdepth
 struct RSM_Result scalar R
 _t35_rows++;_t35_rowdepth++;_t35_maxrowdepth=max((_t35_maxrowdepth,_t35_rowdepth))
 R=_t35_row_impl(P,tol,trace,check_success)
 _t35_rowdepth--
 return(R)
}
'''
a=s.index('struct RSM_Result {');b=s.index('\n}',a)+2;s=s[:b]+'\n'+helper+s[b:]
p.write_text(s);(r/'accounting.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),s.splitlines(True))))
cases=j(r/'GRID_FIXTURE.json');body=f'run "{pkg}/_dea_simplex.ado"\n'+DIRECT
for z in cases:
 id=z['id']+'_p1';n=len(z['c']);E=np.array(z['E']).reshape(-1,n);b=np.array(z['b']).reshape(-1,1);c=np.array(z['c']).reshape(1,-1)
 body+='mata: _t35_work=_t35_depth=_t35_maxdepth=_t35_rows=_t35_rowdepth=_t35_maxrowdepth=0;_t35_calls=J(0,7,.)\n'
 body+=f'mata: trial("{id}",{mata(E)},{mata(b)},{mata(c)},1e-10,1)\n'
 body+=f'mata: savebin("{id}_actual.bin",(_t35_work,_t35_rows,_t35_depth,_t35_maxdepth,_t35_rowdepth,_t35_maxrowdepth));savebin("{id}_calls.bin",_t35_calls)\n'
out=run('accounting',body,pkg);records=[];parity=0
for z in cases:
 id=z['id']+'_p1';a=readbin(out/(id+'_actual.bin')).ravel();calls=readbin(out/(id+'_calls.bin'));res=readbin(out/(id+'_result.bin')).ravel();audit=readbin(out/(id+'_audit.bin'))
 assert a[0]==res[2] and a[1]<=1 and a[2]==a[4]==0 and a[3]<=2 and a[5]<=1,(id,a,res)
 assert audit.shape==(1,12) and audit[0,6]==res[2] and audit[0,4]+audit[0,5]==res[2]
 assert (calls[:,0]==1).sum()==1 and np.all(calls[calls[:,0]==2,1]==0) and np.all(calls[calls[:,0]==2,3]==0)
 for suffix in ['result.bin','x.bin','basis.bin','success.bin','counts.bin','fast.bin','audit.bin','recovery.bin','message.txt']:
  assert h(out/(id+'_'+suffix))==h(r/'checks/grid'/(id+'_'+suffix)),(id,suffix);parity+=1
 records.append(dict(id=z['id'],actual_iterations=int(a[0]),row_retry_calls=int(a[1]),core_calls=len(calls),max_depth=int(a[3])))
dump('ACCOUNTING_REVIEW.json',dict(passed=True,conditions=len(cases),row_retry_calls=sum(z['row_retry_calls'] for z in records),max_core_depth=max(z['max_depth'] for z in records),max_row_retry_depth=1,one_top_level_audit=True,actual_iteration_sums_equal=True,unchanged_files=parity,rows=records,source_sha256={p.name:h(p) for p in pkg.iterdir() if p.is_file()}))
print('Actual iteration accounting verified',len(records),parity,flush=True)
