#!/usr/bin/env python3
"""Self-contained bounded exact-row/empty-basis/recovery regression.
Requires NumPy/SciPy,Python3 and licensed Stata16+. See ROW_RECOVERY.txt.
"""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,sys,time
REPO=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stata',required=True);p.add_argument('--output',required=True,type=Path);args=p.parse_args()
 r=args.output.resolve();r.mkdir(parents=True,exist_ok=False)
 for n in ['checks','logs']:(r/n).mkdir()
 shutil.copytree(REPO/'ado',r/'candidate/ado')
 for f in (REPO/'tests/row_recovery').glob('*.py'):shutil.copy2(f,r/f.name)
 fixture=REPO/'tests/fixtures/row_recovery.json';data=json.loads(fixture.read_text())
 for name,key in [('GRID_FIXTURE.json','grid'),('EXACT_ORACLES.json','oracles'),('REFERENCE.json','reference')]: (r/name).write_text(json.dumps(data[key],indent=2)+'\n')
 (r/'config.json').write_text(json.dumps(dict(stata=args.stata,tests=str(REPO/'tests')),indent=2)+'\n')
 h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();source={p.name:h(p) for p in (r/'candidate/ado').iterdir() if p.is_file()}
 (r/'BUILD.json').write_text(json.dumps(dict(source_sha256=source),indent=2)+'\n')
 (r/'manifest.json').write_text(json.dumps(dict(source_sha256=source,fixtures_sha256={fixture.name:h(fixture)},component_sha256={p.name:h(p) for p in (REPO/'tests/row_recovery').glob('*.py')}),indent=2)+'\n')
 jobs=[('grid.py',[]),('review_grid.py',[]),('math_checks.py',[]),('observe.py',[]),('public.py',[]),('public.py',['--observed']),('review_public.py',[]),('dea.py',[]),('dea.py',['--observed']),('review_dea.py',[]),('accounting.py',[]),('guards.py',[]),('affine_bounds.py',[]),('review.py',[])]
 rows=[]
 for name,extra in jobs:
  start=time.monotonic();label=name.removesuffix('.py')+('_observed' if extra else '')
  with (r/'logs'/(label+'.console')).open('w') as f:result=subprocess.run([sys.executable,'-B',str(r/name),*extra],stdout=f,stderr=subprocess.STDOUT,timeout=1800)
  row=dict(name=label,returncode=result.returncode,seconds=time.monotonic()-start);rows.append(row);(r/'RUNS.json').write_text(json.dumps(rows,indent=2)+'\n');print(row,flush=True)
  assert result.returncode==0,(name,'See the saved console and Stata log.')
 review=json.loads((r/'FINAL_REVIEW.json').read_text());assert review['passed']
 summary=dict(all_passed=True,**review,components=rows)
 (r/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(review,indent=2))
if __name__=='__main__':main()
