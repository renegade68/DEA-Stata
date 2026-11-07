from runlib import *
from test_recoverylog import SNAP
import shutil,difflib
pkg=r/'observed/ado';shutil.copytree(r/'candidate/ado',pkg);p=pkg/'_dea_simplex.ado';old=p.read_text();s=old
helper=WRITER.split('mata:\n',1)[1].rsplit('end',1)[0].replace('savebin','_t39_save')
helper+=r'''
string scalar _t39_tag(string scalar phase) {
 return(st_global("T39_CASE")+"_d"+strofreal(st_numscalar("T39_SEQ"),"%02.0f")+"_"+phase)
}
void _t39_snap(string scalar phase,struct RSM_Prob scalar P,struct RSM_Result scalar R) {
 string scalar id
 if(st_numscalar("T39_SEQ")>st_numscalar("T39_N")) return
 id=_t39_tag(phase)
 _t39_save(id+"_A.bin",P.A[,1..P.n_orig]);_t39_save(id+"_b.bin",P.b);_t39_save(id+"_c.bin",P.c[1..P.n_orig])
 _t39_save(id+"_x.bin",R.x);_t39_save(id+"_basis.bin",R.basiscols);_t39_save(id+"_result.bin",(R.status,R.objval,R.iters))
}
void _t39_bounds(string scalar phase,real matrix XP,real matrix bounds,real rowvector detail) {
 string scalar id
 if(st_numscalar("T39_SEQ")>st_numscalar("T39_N")) return
 id=_t39_tag(phase)
 _t39_save(id+"_xp.bin",XP);_t39_save(id+"_bounds.bin",bounds);_t39_save(id+"_detail.bin",detail)
}
'''
a=s.index('struct RSM_Result {');b=s.index('\n}',a)+2;s=s[:b]+'\n'+helper+'\n'+s[b:]
s=s.replace('capture mata: mata drop _dea_link_reset()', '\n'.join('capture mata: mata drop '+n+'()' for n in ['_t39_save','_t39_tag','_t39_snap','_t39_bounds'])+'\ncapture mata: mata drop _dea_link_reset()',1)
s=s.replace('mode=_dea_link_mode;linktol=tol;oldrad=_dea_extract_theta(R,I)','mode=_dea_link_mode;linktol=tol;oldrad=_dea_extract_theta(R,I)\n st_numscalar("T39_SEQ",st_numscalar("T39_SEQ")+1)\n _t39_snap("old",P,R)',1)
s=s.replace('meta[6]=1;meta[7]=sum(bounds[,2]:<0);','_t39_bounds("old",XP,bounds,detail)\n meta[6]=1;meta[7]=sum(bounds[,2]:<0);',1)
s=s.replace('meta[12]=1;meta[13]=min(bounds[,1]);','_t39_snap("retry",P,T)\n _t39_bounds("retry",XP,bounds,detail)\n meta[12]=1;meta[13]=min(bounds[,1]);',1)
s=s.replace('_dea_link_apply(link_original,info1,res1,k,use_stage,tol,linktol,linkmeta)','_dea_link_apply(link_original,info1,res1,k,use_stage,tol,linktol,linkmeta)\n    _t39_snap("new",link_original,res1)',1)
s=s.replace('struct RSM_Prob scalar prob1, prob2, link_original','struct RSM_Prob scalar prob1, prob2, link_original, t39_stage2',1)
anchor='_dea_lp_stage2(prob2, info2, X, Y, x0, y0, ort, rts, theta_star)';assert s.count(anchor)==2;s=s.replace(anchor,anchor+'\n    t39_stage2=prob2')
anchor='res2 = _rsm_two_phase(prob2, linktol, trace)';assert s.count(anchor)==1;s=s.replace(anchor,anchor+'\n    _t39_snap("s2",t39_stage2,res2)',1)
p.write_text(s);(r/'observed.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),s.splitlines(True))))
data={z['kind']:z for z in j(r/'GRID_DATA.json')};commands=[z for z in j(r/'GRID_COMMANDS.json') if z['stage']==2]
body=f'run "{pkg}/_dea_simplex.ado"\n'
for z in commands:
 v=data[z['dataset']];id=z['id'];body+=f'global T39_CASE "{id}"\nscalar T39_SEQ=0\nscalar T39_N={z["DMUs"]}\nmata: _dea_link_reset(2,1)\nclear\ninput double id '+' '.join('double '+n for n in v['names'])+'\n'
 for k,row in enumerate(v['data'],1):body+=str(k)+' '+' '.join(format(x,'.17g') for x in row)+'\n'
 body+=f'end\ndea {" ".join(v["inputs"])}={" ".join(v["outputs"])},dmu(id) engine(simplex) rts({z["rts"]}) ort({z["ort"]}) stage(2) reference gen(eff)\n'+SNAP.replace('PREFIX',id)+f'mata: savebin("{id}_public_link.bin",st_matrix("r(stage_linkage)"))\n'+f'mata: savebin("{id}_link.bin",_dea_link_read())\n'
out=run('observed',body,pkg);parity=0
for z in commands:
 for p in (r/'checks/grid_repair').glob(z['id']+'_*.bin'):
  if p.name.endswith('_rc.bin'):continue
  assert h(p)==h(out/p.name),p.name;parity+=1
dump('OBSERVED_REVIEW.json',dict(passed=True,commands=len(commands),identical_files=parity,numerical_statements_replaced=False,source_sha256={p.name:h(p) for p in pkg.iterdir() if p.is_file()}))
print('Observer parity:',parity,flush=True)
