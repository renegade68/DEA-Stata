from runlib import *
names=['GRID','EXACT','FAST','GUARD','ACCOUNTING','CONSISTENCY','PUBLIC','SENSITIVITY','OBSERVED']
reports={n:j(r/(n+'_REVIEW.json')) for n in names}
for n,z in reports.items():assert z['passed'],n
assert reports['EXACT']['repaired_main_stage2_records']==201
assert reports['EXACT']['exact_witness_checks']==reports['EXACT']['basis_error_enclosures_verified']==402
assert reports['EXACT']['enumeration_fallbacks']==0
assert reports['FAST']['cases']==535 and reports['FAST']['accepted_outside_trigger_bound']==0
assert reports['GUARD']['direct_checks']==19 and reports['GUARD']['fault_commands']==6
assert reports['ACCOUNTING']['missing_caller_status']==0
summary=dict(passed=True,grid=reports['GRID'],exact={k:v for k,v in reports['EXACT'].items() if k!='records'},fast={k:v for k,v in reports['FAST'].items() if k!='rows'},guards=reports['GUARD'],accounting=reports['ACCOUNTING'],public=reports['PUBLIC'],sensitivity=reports['SENSITIVITY'],consistency=reports['CONSISTENCY'],observer_identical_files=reports['OBSERVED']['identical_files'])
dump('FINAL_REVIEW.json',summary)
print('Stage-link numerical,public and exact-certificate regression passed.',flush=True)
