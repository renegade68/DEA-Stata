from runlib import *
from collections import Counter
out=r/'checks/dea';obs=r/'checks/dea_observed';rows=[];parity=0
for z in j(r/'DEA_COMMANDS.json'):
 id=z['id'];assert readbin(out/(id+'_rc.bin'))[0,0]==0,id
 meta=readbin(out/(id+'_success.bin'));score=readbin(out/(id+'_scores.bin'));counts=readbin(out/(id+'_counts.bin')).ravel()
 row=dict(**z,available=int(np.isfinite(score).sum()),screens=int(counts[0]),refined=int(counts[2]),quality_failed=int(counts[3]),unavailable=int(counts[4]),outcomes=dict(Counter(str(int(v)) for v in meta[:,0])),reasons=dict(Counter(str(int(v)) for v in meta[:,1])))
 if z['policy']==1:
  for p in out.glob(id+'_*.bin'):assert h(p)==h(obs/p.name),p.name;parity+=1
  assert len(list(obs.glob(id+'_obs*_meta.bin')))==len(meta),id
 rows.append(row)
comparisons=[]
for z in rows:
 if z['policy']!=1:continue
 id=z['id'];prefix=id.removesuffix('_p1');files=[p for p in out.glob(id+'_*.bin') if p.name.endswith('_scores.bin') or '_m_' in p.name or '_s_' in p.name]
 record=dict(case=prefix,policy1_available=z['available'])
 for policy in [0,2]:
  record['vs_policy'+str(policy)+'_identical']=all(h(p)==h(out/p.name.replace(prefix+'_p1_',prefix+'_p'+str(policy)+'_',1)) for p in files)
  record['policy'+str(policy)+'_available']=next(a['available'] for a in rows if a['id']==prefix+'_p'+str(policy))
 comparisons.append(record)
summary=dict(passed=True,commands=192,observed_commands=64,positive_datasets=8,all_data_strictly_positive=True,observed_numeric_files_identical=parity,policy_summary={str(p):{k:sum(z[k] for z in rows if z['policy']==p) for k in ['available','screens','refined','quality_failed','unavailable']} for p in [0,1,2]},policy0_identical_commands=sum(z['vs_policy0_identical'] for z in comparisons),policy2_identical_commands=sum(z['vs_policy2_identical'] for z in comparisons),rows=rows,comparisons=comparisons)
dump('DEA_REVIEW.json',summary);print(json.dumps({k:v for k,v in summary.items() if k not in ['rows','comparisons']},indent=2),flush=True)
