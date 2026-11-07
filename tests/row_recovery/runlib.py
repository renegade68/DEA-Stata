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
 script=head+WRITER+f'adopath ++ "{pkg}"\n'+body+'\ndi "ROW_RECOVERY_PASS"\nexit,clear\n'
 (out/'test.do').write_text(script);start=time.monotonic()
 p=subprocess.run([config['stata'],'-b','do','test.do'],cwd=out,capture_output=True,timeout=1800)
 meta=dict(returncode=p.returncode,seconds=time.monotonic()-start,stdout=p.stdout.decode(errors='replace'),stderr=p.stderr.decode(errors='replace'))
 (out/'process.json').write_text(json.dumps(meta,indent=2)+'\n')
 log=(out/'test.log').read_text(errors='replace') if (out/'test.log').exists() else ''
 assert p.returncode==0 and 'ROW_RECOVERY_PASS' in log.splitlines(),(name,log[-4500:])
 (out/'summary.json').write_text(json.dumps(dict(passed=True,marker='ROW_RECOVERY_PASS',source_sha256={p.name:h(p) for p in pkg.iterdir() if p.is_file()},**meta),indent=2)+'\n')
 print(name,'completed',round(meta['seconds'],3),flush=True);return out

DIRECT=r'''
mata:
void trial(string scalar id,real matrix E,real colvector b,real rowvector c,real scalar tol,real scalar policy) {
 struct RSM_Prob scalar P,Q
 struct RSM_Result scalar R
 real scalar f
 _rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c);Q=P
 _rsm_success_reset(policy,1);_rsm_audit_reset(1)
 R=_rsm_two_phase(P,tol,"")
 savebin(id+"_result.bin",(R.status,R.objval,R.iters));savebin(id+"_x.bin",R.x)
 savebin(id+"_basis.bin",R.basiscols);savebin(id+"_success.bin",_rsm_success_read())
 savebin(id+"_counts.bin",_rsm_success_counts_read());savebin(id+"_fast.bin",_rsm_success_fastread())
 savebin(id+"_A.bin",Q.A[,1..Q.n_orig]);savebin(id+"_b.bin",Q.b);savebin(id+"_c.bin",Q.c[1..Q.n_orig])
 savebin(id+"_audit.bin",_rsm_audit_read());savebin(id+"_recovery.bin",R.recovery)
 f=fopen(id+"_message.txt","w");fput(f,R.msg);fclose(f)
}
end
'''
def direct(cases,name):
 body=f'run "{r}/candidate/ado/_dea_simplex.ado"\n'+DIRECT
 for z in cases:
  for policy in [0,1,2]:
   id=z['id']+'_p'+str(policy);E=np.array(z['E']).reshape(-1,len(z['c']));b=np.array(z['b']).reshape(-1,1);c=np.array(z['c']).reshape(1,-1)
   body+=f'capture noisily mata: trial("{id}",{mata(E)},{mata(b)},{mata(c)},{z.get("tol",1e-10):.17g},{policy})\nscalar case_rc=_rc\nmata: savebin("{id}_rc.bin",st_numscalar("case_rc"))\n'
 return run(name,body)
