from runlib import *
import itertools,shutil,difflib
z=next(z for z in j(r/'GRID_DATA.json') if z['kind']=='near_parallel')
def data():
 s='clear\ninput double id double x1 double x2 double y1 double y2 byte selected\n'
 for k,row in enumerate(z['data'],1):s+=str(k*2-1)+' -1 -1 -1 -1 0\n'+str(k*2)+' '+' '.join(format(x,'.17g') for x in row)+' 1\n'
 return s+'end\ngenerate str20 sentinel="original data"\n'
body=f'run "{r}/candidate/ado/_dea_simplex.ado"\n';commands=[]
for rts,ort,stage,form in itertools.product(['crs','vrs','drs','irs'],['input','output'],[1,2],['A','B']):
 id=f'{rts}_{ort}_s{stage}_{form}';commands.append(id)
 body+=data()+'mata: _dea_link_reset(2,0)\n'+f'mata: savebin("{id}_before.bin",st_data(.,tokens("id x1 x2 y1 y2 selected")))\n'
 syntax='x1 x2=y1 y2 if selected,' if form=='A' else 'y1 y2 if selected,inputs(x1 x2)'
 body+=f'dea {syntax} dmu(id) engine(simplex) rts({rts}) ort({ort}) stage({stage}) reference gen(eff) recoverylog\n'
 body+='matrix L=r(stage_linkage)\nmata: assert(rows(_dea_link_read())==0)\n'
 body+='assert r(N)==6 & rowsof(L)==6 & colsof(L)==30\n'
 body+='mata: assert(all(st_matrix("L")[.,1]:==2*(1::6)));assert(sum(st_matrix("L")[.,15])==st_numscalar("r(n_stage_link_adjusted)"));assert(sum(st_matrix("L")[.,10]:<.)==st_numscalar("r(n_stage_link_retries)"))\n'
 body+=f'mata: savebin("{id}_public.bin",st_matrix("L"));savebin("{id}_after.bin",st_data(.,tokens("id x1 x2 y1 y2 selected")));savebin("{id}_score.bin",st_data(.,"eff","selected"))\nassert sentinel=="original data"\nassert missing(eff) if !selected\n'
 # The result is independent of the private recorder capacity and VRS auxiliaries.
 body+='mata: _dea_link_reset(2,1);_dea_link_records=J(10000,29,.)\n'
 body+=f'dea {syntax} dmu(id) engine(simplex) rts({rts}) ort({ort}) stage({stage}) reference gen(eff) replacevars\n'
 body+='mata: assert(all(st_matrix("r(stage_linkage)"):==st_matrix("L")));assert(rows(_dea_link_read())==10000)\n'
 body+='matrix drop L\n'
# Absent metadata must not leak from a preceding simplex command.
body+=data()+'dea x1 x2=y1 y2 if selected,dmu(id) engine(ipm) gen(ipm)\ncapture matrix list r(stage_linkage)\nassert _rc==111\nassert missing(r(n_stage_link_adjusted)) & missing(r(n_stage_link_retries))\n'
body+='dea x1 x2=y1 y2 if selected,dmu(id) engine(simplex) model(sbm) gen(sbm)\ncapture matrix list r(stage_linkage)\nassert _rc==111\nassert missing(r(n_stage_link_adjusted)) & missing(r(n_stage_link_retries))\n'
out=run('public',body);parity=0
for id in commands:
 assert h(out/f'{id}_before.bin')==h(out/f'{id}_after.bin')
 refid='near_parallel_'+id[:-2];L=readbin(out/f'{id}_public.bin');M=readbin(r/'checks/grid_repair'/f'{refid}_public_link.bin')
 assert np.array_equal(L[:,1:],M[:,1:],equal_nan=True),id
 assert h(out/f'{id}_score.bin')==h(r/'checks/grid_repair'/f'{refid}_scores.bin');parity+=1
# Force primal/dual unavailability without changing the distributed package.
faults=[]
for fault in ['dual_only','both_failed']:
 pkg=r/'public_faults'/fault/'ado';shutil.copytree(r/'candidate/ado',pkg);p=pkg/'_dea_simplex.ado';old=p.read_text();s=old
 anchor='res1 = _rsm_two_phase(prob1, tol, trace)';assert s.count(anchor)==1;s=s.replace(anchor,anchor+'\n    res1.status=3',1)
 if fault=='both_failed':
  anchor='res_m = _rsm_two_phase(prob_m, tol, trace)';assert s.count(anchor)==1;s=s.replace(anchor,anchor+'\n    res_m.status=3',1)
 p.write_text(s);(pkg.parent/'fault.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),s.splitlines(True))))
 body=data()
 for stage in [1,2]:
  body+=f'dea x1 x2=y1 y2 if selected,dmu(id) engine(simplex) stage({stage}) gen(e{stage})\n'
  body+='assert r(n_stage_link_adjusted)==0 & r(n_stage_link_retries)==0\nmata: L=st_matrix("r(stage_linkage)");assert(rows(L)==6 & all(L[.,24]:==12) & all(L[.,15]:==0));assert(all(L[.,6]:>=.) & all(L[.,12]:>=.) & all(L[.,29]:>=.))\n'
  if fault=='dual_only':body+='mata: assert(all(L[.,16]:<.) & all(L[.,16]:==L[.,17]));assert(all(st_matrix("r(solution_source)")[.,1]:==2))\n'
  else:body+='mata: assert(all(L[.,16]:>=.) & all(L[.,17]:>=.));assert(st_numscalar("r(n_radial_failed)")==6)\n'
 if fault=='both_failed':body+='dea x1 x2=y1 y2 if selected,dmu(id) engine(auto) gen(auto_eff)\ncapture matrix list r(stage_linkage)\nassert _rc==111\nassert missing(r(n_stage_link_adjusted)) & missing(r(n_stage_link_retries))\n'
 run('public_'+fault,body,pkg);faults.append(fault)
dump('PUBLIC_REVIEW.json',dict(passed=True,main_commands=len(commands),capacity_repeats=len(commands),selected_DMUs_per_command=6,source_observation_indices_preserved=True,private_record_cap_independent=True,original_data_unchanged=True,main_auxiliary_separated=True,missing_metadata_leak_commands=2,dual_failure_commands=4,auto_fallback_metadata_clear_commands=1,faults=faults,grid_score_parity=parity))
