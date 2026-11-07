from runlib import *
from test_recoverylog import SNAP
import itertools
observed='--observed' in sys.argv;pkg=r/('observed' if observed else 'candidate')/'ado'
base=np.array([[1,3,1,2],[2,2,3,1],[3,4,2,3],[4,6,4,5],[5,5,5,4],[6,7,6,7]],dtype=float)
data=[]
for kind in ['base','input_sum','output_sum','both_sum','duplicate_DMU','proportional','near_parallel','units']:
 a=base.copy();inputs=['x1','x2'];outputs=['y1','y2'];names=inputs+outputs
 if kind in ['input_sum','both_sum']:a=np.column_stack([a,a[:,0]+a[:,1]]);names.append('x3');inputs=inputs+['x3']
 if kind in ['output_sum','both_sum']:a=np.column_stack([a,a[:,2]+a[:,3]]);names.append('y3');outputs=outputs+['y3']
 if kind=='duplicate_DMU':a=np.vstack([a,a[:2]])
 if kind=='proportional':a[:,1]=2*a[:,0];a[:,3]=2*a[:,2]
 if kind=='near_parallel':a[:,1]=a[:,0]+2.**(-25)*np.array([0,1,2,0,1,2]);a[:,3]=a[:,2]+2.**(-25)*np.array([0,2,1,2,1,0])
 if kind=='units':a=a*np.array([2.**-10,2.**10,2.**-6,2.**6])
 assert (a>0).all();data.append(dict(kind=kind,names=names,inputs=inputs,outputs=outputs,data=a.tolist()))
body=f'run "{pkg}/_dea_simplex.ado"\nmata: _dea_link_reset(0,0)\n';commands=[]
for z,rts,ort,stage in itertools.product(data,['crs','vrs'],['input','output'],[1,2]):
 for policy in ([1] if observed else [0,1,2]):
  id=f'{z["kind"]}_{rts}_{ort}_s{stage}_p{policy}'
  body+=f'global ROW_CASE "{id}"\nscalar ROW_SEQ=0\nmata: _rsm_success_reset({policy},1)\nclear\ninput double id '+' '.join('double '+n for n in z['names'])+'\n'
  for i,row in enumerate(z['data'],1):body+=str(i)+' '+' '.join(format(v,'.17g') for v in row)+'\n'
  body+='end\n'+f'capture noisily dea {" ".join(z["inputs"])}={" ".join(z["outputs"])},dmu(id) engine(simplex) rts({rts}) ort({ort}) stage({stage}) gen(eff)\nscalar case_rc=_rc\nmata: savebin("{id}_rc.bin",st_numscalar("case_rc"))\n'
  body+='if case_rc==0 {\n'+SNAP.replace('PREFIX',id)+f'mata: savebin("{id}_success.bin",_rsm_success_read());savebin("{id}_counts.bin",_rsm_success_counts_read())\n'+'}\n'
  commands.append(dict(id=id,dataset=z['kind'],rts=rts,ort=ort,stage=stage,policy=policy,DMUs=len(z['data'])))
out=run('dea_observed' if observed else 'dea',body,pkg)
dump('DEA_DATA.json',data);dump('DEA_OBSERVED_COMMANDS.json' if observed else 'DEA_COMMANDS.json',commands)
print('DEA commands',len(commands),flush=True)
