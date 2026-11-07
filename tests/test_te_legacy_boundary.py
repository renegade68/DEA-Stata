#!/usr/bin/env python3
"""Keep public ordinary-TE domain checks separate from legacy IDEA's inner LPs."""
from pathlib import Path
import argparse,itertools,json,subprocess,shutil
from test_tolerance_validation import BASE
from test_recoverylog import SNAP
from test_core_retries import WRITER
REPO=Path(__file__).resolve().parents[1]
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--baseline-ado',type=Path,required=True);a=ap.parse_args()
 out=a.output.resolve();out.mkdir(parents=True);pkg=out/'ado';shutil.copytree(REPO/'ado',pkg);runs=[];cases=[];nfiles=0
 profiles={'ordinal_input':('x','ordinal(x)'), 'ordinal_output':('x','ordinal(y)'), 'ordinal_both':('x','ordinal(x y)'), 'bounded':('xs','bounded(xs)'), 'bounded_ordinal_output':('xs','bounded(xs) ordinal(y)'), 'numeric_ordinal_output':('x','xlo(lo) xhi(hi) ordinal(y)')}
 for profile,(x,extra) in profiles.items():
  commands=[]
  for engine,rts,ori,stage,form in itertools.product(['simplex','ipm','auto'],['crs','vrs'],['input','output'],[1,2],['A','B']):
   head=f'{x}=y,' if form=='A' else f'y, inputs({x})';prefix=f'c{len(cases):03}'
   cmd=f'dea {head} {extra} engine({engine}) rts({rts}) ort({ori}) stage({stage}) dmu(id) gen(eff) tolerance(1e-8)'
   cases.append(dict(prefix=prefix,profile=profile,engine=engine,rts=rts,orientation=ori,stage=stage,form=form,command=cmd));commands.append((prefix,cmd))
  for mode,source in [('baseline',a.baseline_ado.resolve()),('candidate',pkg)]:
   d=out/(profile+'_'+mode);d.mkdir();body=''
   for prefix,cmd in commands:
    body+='clear\n'+BASE+cmd+'\n'+SNAP.replace('PREFIX',prefix)
    # An internal IDEA exception must not leak into the next public command.
    if mode=='candidate':
     body+='replace y=0\ncapture noisily dea x=y, dmu(id) gen(eff) replacevars\nassert _rc==459\n'
     body+='replace y=1\nreplace x=-1 in 1\ncapture noisily dea x=y, dmu(id) gen(eff) replacevars recoverylog\nassert _rc==459\nmata: assert(_dea_rec_on==0)\n'
   (d/'test.do').write_text(f'clear all\nset more off\nadopath ++ "{source}"\n'+WRITER+body+'di "TEST_PASS"\nexit,clear\n')
   p=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=300);log=(d/'test.log').read_text() if (d/'test.log').exists() else '';ok=p.returncode==0 and 'TEST_PASS' in log.splitlines();runs.append(dict(test=d.name,passed=ok));(out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n');assert ok,(d.name,log[-3000:]);print(runs[-1],flush=True)
  before=out/(profile+'_baseline');after=out/(profile+'_candidate')
  assert sorted(p.name for p in before.glob('*.bin'))==sorted(p.name for p in after.glob('*.bin'))
  for p in before.glob('*.bin'):assert p.read_bytes()==(after/p.name).read_bytes(),p.name;nfiles+=1
 summary=dict(legacy_combinations=len(cases),legacy_executions=2*len(cases),subsequent_invalid_public_commands=2*len(cases),numeric_files_identical=nfiles,batches=len(runs),all_passed=True)
 (out/'cases.json').write_text(json.dumps(cases,indent=2)+'\n');(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)
if __name__=='__main__':main()
