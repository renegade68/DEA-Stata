from runlib import *
from test_recoverylog import SNAP
import itertools
label=sys.argv[1] if len(sys.argv)>1 else 'accounting';z=next(z for z in j(r/'GRID_DATA.json') if z['kind']=='near_parallel');body=f'run "{r}/candidate/ado/_dea_simplex.ado"\n';commands=[]
for rts,stage,on in itertools.product(['crs','vrs'],[1,2],[False,True]):
 id=f'{rts}_s{stage}_{int(on)}';body+='mata: _dea_link_reset(2,1)\nclear\ninput double id double x1 double x2 double y1 double y2\n'
 for k,row in enumerate(z['data'],1):body+=str(k)+' '+' '.join(format(x,'.17g') for x in row)+'\n'
 body+=f'end\ndea x1 x2=y1 y2,dmu(id) engine(simplex) rts({rts}) ort(output) stage({stage}) reference gen(eff) '+('recoverylog' if on else '')+'\n'+SNAP.replace('PREFIX',id)+f'mata: savebin("{id}_link.bin",_dea_link_read())\n'
 commands.append(dict(id=id,rts=rts,stage=stage,on=on))
d=run(label,body,r/'candidate/ado');missing=0;rows=[];parity=0
for z in commands:
 if not z['on']:continue
 id=z['id'];off=id[:-1]+'0';M=readbin(d/f'{id}_link.bin');A=readbin(d/f'{id}_m_recovery_log.bin');C=readbin(d/f'{id}_m_recovery_stats.bin').ravel()
 missing+=int(np.isnan(A[:,19]).sum())
 expected=6*z['stage']+(18 if z['rts']=='vrs' else 0)+int(np.isfinite(M[:,10]).sum())
 assert len(A)==C[0]==C[7]==expected and C[6]==0
 assert np.array_equal(A[:,6],np.arange(1,len(A)+1))
 assert np.array_equal(A[:,11]+A[:,12],A[:,13])
 assert np.array_equal(A[:,[7,8,11,12,13]].sum(axis=0),C[1:6])
 for p in d.glob(off+'_*.bin'):
  assert h(p)==h(d/p.name.replace(off+'_',id+'_',1)),p.name;parity+=1
 rows.append(dict(**z,LP_calls=expected,retries=int(np.isfinite(M[:,10]).sum()),missing_caller_status=int(np.isnan(A[:,19]).sum())))
out=dict(passed=missing==0,commands=8,missing_caller_status=missing,logging_numeric_identical_files=parity,rows=rows)
dump(label.upper()+'_REVIEW.json',out);print(out,flush=True)
if label=='accounting':assert missing==0
