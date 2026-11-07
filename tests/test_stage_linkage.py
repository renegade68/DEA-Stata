#!/usr/bin/env python3
"""Bounded DEA stage-link regression; see STAGE_LINKAGE.txt."""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,sys,time
REPO=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stata',required=True);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
 r=a.output.resolve();r.mkdir(parents=True,exist_ok=False)
 for n in ['checks','logs']:(r/n).mkdir()
 shutil.copytree(REPO/'ado',r/'candidate/ado')
 for f in (REPO/'tests/stage_linkage').glob('*.py'):shutil.copy2(f,r/f.name)
 fixture=REPO/'tests/fixtures/stage_linkage.json';z=json.loads(fixture.read_text())
 for name,key in [('GRID_DATA.json','data'),('GRID_COMMANDS.json','commands'),('FAST_FIXTURE.json','fast'),('REFERENCE.json','reference')]: (r/name).write_text(json.dumps(z[key],indent=2)+'\n')
 (r/'config.json').write_text(json.dumps(dict(stata=a.stata,tests=str(REPO/'tests')),indent=2)+'\n')
 h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
 (r/'manifest.json').write_text(json.dumps(dict(source_sha256={p.name:h(p) for p in (r/'candidate/ado').iterdir() if p.is_file()},fixture_sha256=h(fixture),components={p.name:h(p) for p in (REPO/'tests/stage_linkage').glob('*.py')}),indent=2)+'\n')
 rows=[]
 for name in ['grid','review_grid','observe','exact_review','fast_math','guards','accounting','consistency','public','sensitivity','review']:
  start=time.monotonic()
  with (r/'logs'/(name+'.console')).open('w') as f:result=subprocess.run([sys.executable,'-B',str(r/(name+'.py'))],stdout=f,stderr=subprocess.STDOUT,timeout=1800)
  row=dict(name=name,returncode=result.returncode,seconds=time.monotonic()-start);rows.append(row);(r/'RUNS.json').write_text(json.dumps(rows,indent=2)+'\n');print(row,flush=True)
  assert result.returncode==0,(name,'Inspect saved console and Stata log.')
 summary=json.loads((r/'FINAL_REVIEW.json').read_text());assert summary['passed']
 (r/'summary.json').write_text(json.dumps(dict(all_passed=True,**summary,components=rows),indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
