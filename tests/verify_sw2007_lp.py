#!/usr/bin/env python3
"""Audit sequential SW2007 pseudo-frontiers, inverse-score correction and inference.

Uses the production draws recorded by a test-only source copy and independent
SciPy HiGHS LPs. NumPy checks covariance and both interval types from all B2 fits.
"""
from _support import package_files
import argparse,difflib,hashlib,json,subprocess,time
from pathlib import Path
import numpy as np
import scipy
from verify_bootstrap_lp import solve,read_matrix
from test_sw2007 import DATA,REPO


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir()
    names=package_files(REPO);hashes={}
    for name in names:
        data=(REPO/name).read_bytes();(pkg/Path(name).name).write_bytes(data);hashes[name]=hashlib.sha256(data).hexdigest()
    f=pkg/'_dea_sw2007_seq.ado';orig=f.read_text()
    record='''    st_matrix("AUDIT_X",X)
    st_matrix("AUDIT_Y",Y)
    st_matrix("AUDIT_ORIG",theta)
    if (cols(st_matrix("AUDIT_DRAWS"))==0) {
        st_matrix("AUDIT_DRAWS",d)
        st_matrix("AUDIT_BOOT",result)
    }
    else {
        st_matrix("AUDIT_DRAWS",(st_matrix("AUDIT_DRAWS"),d))
        st_matrix("AUDIT_BOOT",(st_matrix("AUDIT_BOOT"),result))
    }
'''
    needle='    st_store(ix,dest,result)';assert orig.count(needle)==1
    mod=orig.replace(needle,record+needle);f.write_text(mod)
    (out/'instrumentation.patch').write_text(''.join(difflib.unified_diff(orig.splitlines(True),mod.splitlines(True))))
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,numpy=np.__version__,scipy=scipy.__version__),indent=2))
    results=[]
    for dim in [1,2]:
        for rts in ['crs','vrs','drs','irs']:
            for ori in ['input','output']:
                name=f'{dim}d_{rts}_{ori}';work=out/name;work.mkdir();io='x=y' if dim==1 else 'x x2=y y2'
                s=f'clear all\nset more off\nset linesize 120\nadopath ++ "{pkg}"\n'+DATA+f'''dea_sw2007 {io}, envars(z) rts({rts}) ort({ori}) tolerance(1e-8) reps(30) bcreps(20) seed(127) gen(a)
matrix AUDIT_B=e(b_reps)
matrix AUDIT_POINT=e(b)
matrix AUDIT_V=e(V)
matrix AUDIT_CI=e(ci_percentile)
matrix AUDIT_PIVOT=e(ci_pivot)
mkmat a_eff_bc,matrix(AUDIT_BC)
'''
                for m in ['X','Y','ORIG','DRAWS','BOOT','B','POINT','V','CI','PIVOT','BC']:
                    s+=f'preserve\nclear\nsvmat double AUDIT_{m},names(v)\nexport delimited using "{m}.csv", replace\nrestore\n'
                (work/'test.do').write_text(s+'di "TEST_PASS"\n');row=dict(test=name,passed=False);start=time.monotonic()
                try:
                    subprocess.run([args.stata,'-b','do','test.do'],cwd=work,timeout=240,capture_output=True)
                    log=(work/'test.log').read_text(errors='replace')
                    if 'TEST_PASS' not in log.splitlines():raise RuntimeError(log[-1800:])
                    x,y,theta,draws,boot,B,point,V,CI,pivot,bc=[read_matrix(work/f'{m}.csv') for m in ['X','Y','ORIG','DRAWS','BOOT','B','POINT','V','CI','PIVOT','BC']]
                    theta=theta[:,0];n,bcount=draws.shape;expected=np.empty_like(boot)
                    initial=np.array([solve(x[i],y[i],x,y,ori,rts) for i in range(n)])
                    for b in range(bcount):
                        ratio=1/(theta*draws[:,b]);xx=x/ratio[:,None] if ori=='input' else x;yy=y*ratio[:,None] if ori=='output' else y
                        for i in range(n):expected[i,b]=1/solve(x[i],y[i],xx,yy,ori,rts)
                    expected_bc=1/(2/theta-expected.mean(axis=1))
                    expected_ci=np.quantile(B,[.025,.975],axis=0,method='inverted_cdf')
                    expected_pivot=2*point-expected_ci[::-1]
                    row.update(lp_count=n*(1+bcount),max_initial_gap=float(abs(initial-theta).max()),max_lp_gap=float(abs(expected-boot).max()),max_bc_gap=float(abs(expected_bc-bc[:,0]).max()),max_cov_gap=float(abs(np.cov(B,rowvar=False)-V).max()),max_ci_gap=float(abs(expected_ci-CI).max()),max_pivot_gap=float(abs(expected_pivot-pivot).max()))
                    row['passed']=max(row[k] for k in row if k.startswith('max_'))<1e-7
                    np.savetxt(work/'independent_boot.csv',expected,delimiter=',',fmt='%.17g')
                except Exception as exc:row['error']=str(exc)
                row['seconds']=round(time.monotonic()-start,2);results.append(row);(out/'results.json').write_text(json.dumps(results,indent=2));print(json.dumps(row),flush=True)
    summary=dict(total=len(results),passed=sum(r['passed'] for r in results),lp_count=sum(r.get('lp_count',0) for r in results));(out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
    raise SystemExit(not all(r['passed'] for r in results))

if __name__=='__main__':main()
