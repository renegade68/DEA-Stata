from pathlib import Path
from fractions import Fraction as F
import hashlib,json,re,subprocess,sys,time
import numpy as np
r=Path(__file__).resolve().parent;sys.dont_write_bytecode=True;config=json.loads((r/'config.json').read_text());sys.path.insert(0,config['tests'])
from test_core_retries import WRITER,readbin
from test_simplex_stability import mata
j=lambda p:json.loads(p.read_text());h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
dump=lambda n,x:(r/n).write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n')
def run(name,body,pkg=None):
 pkg=pkg or r/'candidate/ado';out=r/'checks'/name;out.mkdir()
 head='version 16\nclear all\nset more off\nset processors 1\nset linesize 180\n'
 for k in ['PLUS','PERSONAL','SITE','OLDPLACE']:(out/k).mkdir();head+=f'sysdir set {k} "{out/k}"\n'
 script=head+WRITER+f'adopath ++ "{pkg}"\n'+body+'\ndi "STAGE_LINK_PASS"\nexit,clear\n'
 (out/'test.do').write_text(script);start=time.monotonic()
 p=subprocess.run([config['stata'],'-b','do','test.do'],cwd=out,capture_output=True,timeout=1800)
 meta=dict(returncode=p.returncode,seconds=time.monotonic()-start,stdout=p.stdout.decode(errors='replace'),stderr=p.stderr.decode(errors='replace'))
 (out/'process.json').write_text(json.dumps(meta,indent=2)+'\n')
 log=(out/'test.log').read_text(errors='replace') if (out/'test.log').exists() else ''
 assert p.returncode==0 and 'STAGE_LINK_PASS' in log.splitlines(),(name,log[-4500:])
 (out/'summary.json').write_text(json.dumps(dict(passed=True,marker='STAGE_LINK_PASS',source_sha256={p.name:h(p) for p in pkg.iterdir() if p.is_file()},**meta),indent=2)+'\n')
 print(name,'completed',round(meta['seconds'],3),flush=True);return out
