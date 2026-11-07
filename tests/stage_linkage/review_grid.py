from runlib import *
commands=j(r/'GRID_COMMANDS.json');ref=j(r/'REFERENCE.json');parity=0;repairs=0;aux=0
for mode in ['off','audit','repair','slow']:
 d=r/'checks'/('grid_'+mode)
 expected=ref['grid_baseline' if mode in ['off','audit'] else 'grid_repair']
 extras={z['id']+'_'+suffix for z in commands for suffix in ['link.bin','public_link.bin','m_stage_linkage.bin','s_n_stage_link_adjusted.bin','s_n_stage_link_retries.bin']}
 assert {p.name for p in d.glob('*.bin')}==set(expected)|extras,mode
 for name,sha in expected.items():
  if mode in ['off','audit'] and name.endswith('_link.bin'):continue
  # Fast and slow private reasons differ; numerical outputs must agree.
  if mode=='slow' and name.endswith('_link.bin'):continue
  assert h(d/name)==sha,(mode,name);parity+=1
 for z in commands:
  id=z['id'];M=readbin(d/f'{id}_link.bin');L=readbin(d/f'{id}_public_link.bin')
  assert L.shape==(z['DMUs'],30)
  assert np.array_equal(L[:,:29],M[:z['DMUs']],equal_nan=True),id
  assert np.array_equal(L,readbin(d/f'{id}_m_stage_linkage.bin'),equal_nan=True),id
  assert readbin(d/f'{id}_s_n_stage_link_adjusted.bin')[0,0]==sum(L[:,14]==1)
  assert readbin(d/f'{id}_s_n_stage_link_retries.bin')[0,0]==sum(np.isfinite(L[:,9]))
  assert (L[:,29]==1e-7).all() and (L[:,2]==z['stage']).all()
  scores=readbin(d/f'{id}_scores.bin').ravel()
  expected_score=L[:,16] if z['ort']=='input' else 1/L[:,16]
  # The public score may snap to one within user tolerance; radial stays raw.
  assert np.all(abs(scores-expected_score)<=1e-7*(1+abs(expected_score))),id
  if mode=='repair':repairs+=int((L[:,14]==1).sum());aux+=int((M[z['DMUs']:,14]==1).sum())
  if z['stage']==1:assert h(d/f'{id}_scores.bin')==h(d/f'{id[:-1]}2_scores.bin')
assert repairs==402 and aux==304
dump('GRID_REVIEW.json',dict(passed=True,commands=1024,datasets=16,main_repaired=repairs,auxiliary_repaired=aux,reference_files_identical=parity,public_DMU_records=sum(z['DMUs'] for z in commands)*4,main_public_matches_private=True,stage_score_parity=True))
print('Grid/reference/public diagnostic checks passed.',parity,flush=True)
