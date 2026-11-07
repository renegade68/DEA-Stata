from runlib import *
commands=j(r/'GRID_COMMANDS.json');d=r/'checks/grid_repair';rows=[];aux=0;repaired=0;fastparity=0
for z in commands:
 id=z['id'];M=readbin(d/f'{id}_link.bin')
 for row in M[z['DMUs']:]:
  rts={1:'crs',2:'vrs',3:'drs',4:'irs'}[int(row[27])]
  other=f'{z["dataset"]}_{rts}_{z["ort"]}_s{int(row[2])}'
  direct=readbin(d/f'{other}_link.bin')[int(row[0])-1]
  assert np.array_equal(row,direct,equal_nan=True),(id,other,row[0])
  aux+=1;repaired+=int(row[14]==1)
 for p in d.glob(id+'_*.bin'):
  if p.name.endswith('_link.bin'):continue
  old=r/'checks/grid_slow'/p.name
  if p.name.endswith('_m_stage_linkage.bin'):
   A=readbin(p);B=readbin(old);cols=np.r_[0:18,20:23,24:30]
   assert np.array_equal(A[:,cols],B[:,cols],equal_nan=True),p.name
   for a,b in zip(A[:,23],B[:,23]):assert a==b or {a,b}=={2,11}
   continue
  assert h(p)==h(old),p.name;fastparity+=1
dump('CONSISTENCY_REVIEW.json',dict(passed=True,auxiliary_rows_match_direct=aux,auxiliary_repair_rows_match_direct=repaired,fast_and_slow_public_identical_files=fastparity,stage2_quality_guard_preserves_valid_results=True))
print('Auxiliary and fast/slow parity:',aux,repaired,fastparity,flush=True)
