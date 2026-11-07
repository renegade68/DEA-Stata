from runlib import *
from test_recoverylog import SNAP
import itertools
data=next(z for z in j(r/'GRID_DATA.json') if z['kind']=='near_parallel');commands=[]
for engine,rts,ort,stage,ti in itertools.product(['simplex','auto','ipm'],['crs','vrs'],['input','output'],[1,2],range(5)):
 commands.append(dict(id=f'{engine}_{rts}_{ort}_s{stage}_t{ti}',engine=engine,rts=rts,ort=ort,stage=stage,tol=[1e-6,1e-7,1e-9,1e-10,1e-12][ti],policy=1))
for policy,stage in itertools.product([0,2],[1,2]):commands.append(dict(id=f'policy{policy}_s{stage}',engine='simplex',rts='crs',ort='output',stage=stage,tol=1e-7,policy=policy))
pkg=r/'candidate/ado';body=f'run "{pkg}/_dea_simplex.ado"\n'
for z in commands:
 id=z['id'];body+=f'mata: _rsm_success_reset({z["policy"]},1);_dea_link_reset(2,1)\nclear\ninput double id double x1 double x2 double y1 double y2\n'
 for k,row in enumerate(data['data'],1):body+=str(k)+' '+' '.join(format(x,'.17g') for x in row)+'\n'
 body+=f'end\ndea x1 x2=y1 y2,dmu(id) engine({z["engine"]}) rts({z["rts"]}) ort({z["ort"]}) stage({z["stage"]}) tolerance({z["tol"]:.17g}) reference gen(eff)\nscalar case_rc=_rc\nmata: savebin("{id}_rc.bin",st_numscalar("case_rc"))\n'+SNAP.replace('PREFIX',id)+f'mata: savebin("{id}_link.bin",_dea_link_read())\n'
 if z['engine']!='ipm':body+='mata: assert(all(st_matrix("r(stage_linkage)")[.,30]:=='+format(z['tol'],'.17g')+'))\n'
d=run('sensitivity_candidate',body,pkg);ref=j(r/'REFERENCE.json');parity=ipm=auto=0
for name,sha in ref['sensitivity_candidate'].items():assert h(d/name)==sha,name;parity+=1
for z in commands:
 id=z['id']
 if z['engine']=='ipm':
  for name,sha in ref['sensitivity_baseline'].items():
   if name.startswith(id+'_'):assert h(d/name)==sha;ipm+=1
 if z['engine']=='auto':
  for name,sha in ref['sensitivity_baseline'].items():
   if not name.startswith(id+'_'):continue
   if name.endswith('_s_auto_scale_ratio.bin'):assert h(d/name)==sha;continue
   assert h(d/name)==h(d/name.replace('auto_','simplex_',1));auto+=1
 if z['stage']==1:assert h(d/f'{id}_scores.bin')==h(d/f'{id.replace("_s1","_s2")}_scores.bin')
dump('SENSITIVITY_REVIEW.json',dict(passed=True,commands=len(commands),reference_files_identical=parity,IPM_identical_files=ipm,auto_simplex_identical_files=auto,stage_score_parity=True))
