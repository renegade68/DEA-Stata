from runlib import *
from test_recoverylog import SNAP
import itertools
data=j(r/"GRID_DATA.json")
commands=[]
for z,rts,ort,stage in itertools.product(data,['crs','vrs','drs','irs'],['input','output'],[1,2]):
 commands.append(dict(id=f'{z["kind"]}_{rts}_{ort}_s{stage}',dataset=z['kind'],rts=rts,ort=ort,stage=stage,DMUs=len(z['data'])))
dump('GRID_COMMANDS.json',commands)
for mode in ['off','audit','repair','slow']:
 pkg=r/'candidate/ado';body=f'run "{pkg}/_dea_simplex.ado"\n'
 for z in commands:
  v=next(d for d in data if d['kind']==z['dataset']);id=z['id']
  body+=f'mata: _dea_link_reset({dict(off=0,audit=1,repair=2,slow=2)[mode]},1,{int(mode!="slow")})\n'
  body+='clear\ninput double id '+' '.join('double '+n for n in v['names'])+'\n'
  for k,row in enumerate(v['data'],1):body+=str(k)+' '+' '.join(format(x,'.17g') for x in row)+'\n'
  body+=f'end\ncapture noisily dea {" ".join(v["inputs"])}={" ".join(v["outputs"])},dmu(id) engine(simplex) rts({z["rts"]}) ort({z["ort"]}) stage({z["stage"]}) reference gen(eff)\nscalar case_rc=_rc\nmata: savebin("{id}_rc.bin",st_numscalar("case_rc"))\n'
  body+='if case_rc==0 {\n'+SNAP.replace('PREFIX',id)+f'mata: savebin("{id}_public_link.bin",st_matrix("r(stage_linkage)"))\n'+'}\n'
  body+=f'mata: savebin("{id}_link.bin",_dea_link_read())\n'
 run('grid_'+mode,body,pkg)
print('Grid complete:',len(data),'datasets;',len(commands),'commands per mode; 4 modes.',flush=True)
