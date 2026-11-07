#!/usr/bin/env python3
"""Panel native LP tolerance/residual audit and component failure preservation."""
from _support import package_files
from pathlib import Path
import argparse,hashlib,itertools,json,re,shutil,subprocess,time
from test_nonte_tolerance import instrument,verify,CAPTURE
from test_core_retries import WRITER
REPO=Path(__file__).resolve().parents[1]
DATA='''input double(id period x1 x2 y1 y2)
1 1 2 4 2 1
2 1 4 2 1 2
3 1 3 3 1 1
1 3 2 4 2.4 1.2
2 3 4 2 1.2 2.4
3 3 3 3 1.2 1.2
1 7 2 4 2.6 1.3
2 7 4 2 1.3 2.6
3 7 3 3 1.3 1.3
end
'''
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--only');a=ap.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir();h={}
    for n in package_files(REPO):
        raw=(REPO/n).read_bytes();(pkg/Path(n).name).write_bytes(raw);h[n]=hashlib.sha256(raw).hexdigest()
    audit=out/'audit_ado';instrument(pkg,audit);(out/'manifest.json').write_text(json.dumps(dict(source_sha256=h),indent=2)+'\n');results=[]
    def run(name,body,source=audit,check=None):
        if a.only and not re.search(a.only,name):return
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset processors 1\nset type double\n'
        for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/key).mkdir();head+=f'sysdir set {key} "{d/key}"\n'
        head+=f'adopath ++ "{source}"\n'+WRITER+CAPTURE+'scalar r08calls=0\n';(d/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n');r=dict(test=name,passed=False);t=time.monotonic()
        try:
            p=subprocess.run([a.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=180);log=(d/'test.log').read_text(errors='replace')
            assert p.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2000:]
            if check:r.update(check(d))
            r['passed']=True
        except Exception as e:r['error']=repr(e)
        r['seconds']=round(time.monotonic()-t,3);results.append(r);(out/'results.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(r),flush=True)
    def cmd(kind,rts,ori,tol):
        return f'dea_{"window" if kind=="window" else "malmquist"} x1 x2=y1 y2,dmu(id) time(period) gen(res) rts({rts}) ort({ori}) tol({tol}) '+('width(2)' if kind=='window' else ('global' if kind=='global' else ''))+'\n'
    for kind,ori,tol,variant in itertools.product(['adjacent','global','window'],['input','output'],['1e-4','1e-7','1e-10'],['base','units','zero']):
        for rts in (['crs','vrs','drs','irs'] if kind=='window' else ['crs','vrs']):
            body=DATA+f'scalar expected_tol={tol}\n'
            if variant=='units':body+='replace x1=x1*1e-12\nreplace x2=x2*1e9\nreplace y1=y1*1e8\nreplace y2=y2*1e-9\n'
            if variant=='zero':body+='replace x1=0 if id==1\nreplace y1=0 if id==2\n'
            body+=cmd(kind,rts,ori,tol)+f'assert r(tolerance)=={tol}\nassert r08calls>0\n'
            run(f'lp_{kind}_{rts}_{ori}_{tol}_{variant}',body,check=lambda d,tol=tol:verify(d,float(tol)))
    for fault in ['all','third','bound','equality','nonfinite']:
        target=out/f'fault_{fault}';shutil.copytree(pkg,target);p=target/'_dea_nonte_lp.ado';s=p.read_text()
        if fault in ['all','third']:
            anchor='    c=S.getCoefficients(); n=cols(c)'
            change='    return(1)\n' if fault=='all' else '    st_numscalar("fault_calls",st_numscalar("fault_calls")+1)\n    if (st_numscalar("fault_calls")==3) return(1)\n'
            assert anchor in s;s=s.replace(anchor,change+anchor)
        else:
            change={'bound':'x[1]=-1000','equality':'x[1]=x[1]+1000','nonfinite':'x[1]=.'}[fault]
            s=s.replace('        x=S.parameters()','        x=S.parameters()\n        '+change)
        p.write_text(s)
        for kind,ori in itertools.product(['adjacent','global','window'],['input','output']):
            body=DATA+'keep if period<=3\nscalar fault_calls=0\n'+cmd(kind,'vrs',ori,'1e-7')
            if fault!='third':
                if kind=='window':body+='assert r(n_obs)==0 & r(n_failed)==6 & r(n_attempted)==6 & r(n_failed_obs)==6\nassert missing(r(mean_eff),r(sd_eff))\nassert missing(res_mean,res_sd)\n'
                else:body+='assert r(n_dmus)==3 & r(n_failed)==3 & r(n_failed_ec)==3 & r(n_failed_tc)==3 & r(n_failed_mpi)==3\nassert missing(r(mean_ec),r(mean_tc),r(mean_mpi))\nassert missing(res_ec,res_tc,res_mpi)\n'
            elif kind=='adjacent':
                # Period growth can independently make VRS cross distances infeasible.
                body+='assert !missing(res_ec) if period==3\nassert missing(res_tc,res_mpi) if period==3 & id==1\n'
            elif kind=='global':body+='assert !missing(res_mpi) if period==3\nassert missing(res_ec,res_tc) if period==3 & id==1\nassert r(n_failed)==1 & r(n_failed_ec)==1 & r(n_failed_mpi)==0\n'
            else:body+='assert r(n_obs)==5 & r(n_failed)==1 & r(n_attempted)==6 & r(n_failed_obs)==1\nassert missing(res_mean) if period==1 & id==3\n'
            run(f'fault_{fault}_{kind}_{ori}',body,source=target)
    # Analytic CRS growth, both directions, adjacent and fixed global pool.
    for kind,ori in itertools.product(['adjacent','global'],['input','output']):
        body='input double(id period x1 x2 y1 y2)\n1 1 1 2 1 2\n2 1 2 4 1 2\n1 3 1 2 2 4\n2 3 2 4 2 4\n1 7 1 2 4 8\n2 7 2 4 4 8\nend\n'+cmd(kind,'crs',ori,'1e-9')
        body+='assert abs(res_ec-1)<1e-7 & abs(res_tc-2)<1e-7 & abs(res_mpi-2)<1e-7 if period>1\nassert abs(res_mpi[3]*res_mpi[5]-4)<1e-7\n'
        run(f'analytic_growth_{kind}_{ori}',body,source=pkg)
    for kind in ['malmquist','window']:
        body=f'run "{pkg}/dea_{kind}.ado"\nrun "{pkg}/dea_{kind}.ado"\n'+DATA
        body+=cmd('adjacent' if kind=='malmquist' else 'window','crs','input','1e-9')
        body+='assert '+('r(n_failed_total)' if kind=='malmquist' else 'r(n_failed)')+'==0\n'
        run('reload_'+kind,body,source=pkg)
    summary=dict(total=len(results),passed=sum(r['passed'] for r in results),failed=sum(not r['passed'] for r in results),LPs=sum(r.get('lp_checks',0) for r in results),rejected_feasible=sum(r.get('rejected_feasible',0) for r in results));(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary);return bool(summary['failed'])
if __name__=='__main__':raise SystemExit(main())
