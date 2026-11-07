from runlib import *
import difflib,shutil
pkg=r/'observed/ado';shutil.copytree(r/'candidate/ado',pkg)
p=pkg/'_dea_simplex.ado';old=p.read_text();s=old
helper=WRITER.split('mata:\n',1)[1].rsplit('end',1)[0].replace('savebin','_t34_save')
helper+=r'''
void _t34_observe(struct RSM_Prob scalar P,struct RSM_Result scalar R,
 real colvector basic,real rowvector meta,real colvector L,real scalar tol) {
 real scalar seq
 string scalar id
 seq=st_numscalar("ROW_SEQ")+1;st_numscalar("ROW_SEQ",seq)
 id=st_global("ROW_CASE")+"_obs"+strofreal(seq,"%04.0f")
 _t34_save(id+"_A.bin",P.A[,1..P.n_orig]);_t34_save(id+"_b.bin",P.b)
 _t34_save(id+"_c.bin",P.c[1..P.n_orig]);_t34_save(id+"_basis.bin",basic)
 _t34_save(id+"_x.bin",R.x);_t34_save(id+"_meta.bin",meta)
 _t34_save(id+"_L.bin",L);_t34_save(id+"_tol.bin",tol)
}
'''
a=s.index('struct RSM_Result {');b=s.index('\n}',a)+2;s=s[:b]+'\n'+helper+'\n'+s[b:]
anchor=' _rsm_success_counts[1]=_rsm_success_counts[1]+1';assert s.count(anchor)==1
s=s.replace(anchor,' if(args()==5) _t34_observe(original,candidate,basic,meta,L,tol)\n else _t34_observe(original,candidate,basic,meta,J(original.n_orig,1,0),tol)\n'+anchor,1)
anchor='capture mata: mata drop _rsm_success_reset()'
s=s.replace(anchor,'capture mata: mata drop _t34_save()\ncapture mata: mata drop _t34_observe()\n'+anchor,1)
p.write_text(s)
(r/'observed.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),s.splitlines(True),fromfile='baseline/_dea_simplex.ado',tofile='observed/_dea_simplex.ado')))
dump('OBSERVATION_SOURCE.json',dict(passed=True,numerical_operations_replaced=False,modified=['_dea_simplex.ado'],source_sha256={p.name:h(p) for p in pkg.iterdir() if p.is_file()}))
print('Read-only numerical observer built in isolated source copy.',flush=True)
