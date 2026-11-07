from runlib import *
observed='--observed' in sys.argv;reduced='--reduced' in sys.argv;assert not (observed and reduced)
pkg=r/('observed' if observed else 'candidate')/'ado'
cases=j(r/('REDUCED_FIXTURE.json' if reduced else 'GRID_FIXTURE.json'));body=f'run "{pkg}/_dea_simplex.ado"\n';commands=[]
for z in cases:
 for scale in ['none','row']:
  for policy in ([1] if observed or reduced else [0,1,2]):
   id=z['id']+'_'+scale+'_p'+str(policy);n=len(z['c']);names=' '.join('x'+str(i+1) for i in range(n))
   body+=f'global ROW_CASE "{id}"\nscalar ROW_SEQ=0\nmata: _rsm_success_reset({policy},1)\nclear\ninput '+' '.join('double '+name for name in names.split())+' str2 rel double rhs\n'
   body+=' '.join(format(v,'.17g') for v in z['c'])+' "" .\n'
   for ar,b in zip(z['E'],z['b']):body+=' '.join(format(v,'.17g') for v in ar)+f' "=" {b:.17g}\n'
   body+=f'end\ncapture noisily lp {names},min engine(simplex) scaling({scale}) tolerance(1e-10)\nscalar case_rc=_rc\nmata: savebin("{id}_rc.bin",st_numscalar("case_rc"))\n'
   body+='if case_rc==0 {\n'+f'mata: savebin("{id}_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")));savebin("{id}_x.bin",st_matrix("r(x)"));savebin("{id}_success.bin",_rsm_success_read());savebin("{id}_attempts.bin",st_matrix("r(attempts)"))\n'+'}\n'
   commands.append(dict(id=id,fixture=z['id'],scale=scale,policy=policy))
out=run('public_observed' if observed else 'public_reduced' if reduced else 'public',body,pkg)
dump('PUBLIC_OBSERVED_COMMANDS.json' if observed else 'PUBLIC_REDUCED_COMMANDS.json' if reduced else 'PUBLIC_COMMANDS.json',commands)
print('Public commands',len(commands),flush=True)
