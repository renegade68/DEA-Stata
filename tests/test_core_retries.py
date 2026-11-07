#!/usr/bin/env python3
"""Native retry integration: per-DMU evidence, initial-solution parity and scope.
Requires licensed Stata 16+, NumPy/SciPy. Faults exist only in private copies.
"""
from _support import package_files
import argparse, hashlib, itertools, json, re, shutil, struct, subprocess, time
from pathlib import Path
import numpy as np
from test_ipm_stages import FIXTURES, script as projection_script, oracle
from test_lp_residuals import residual
from test_saveref import instrument as auto_instrument

REPO=Path(__file__).resolve().parents[1]
SIMPLE='input double id x y\n1 1 1\n2 2 4\n3 4 2\nend\n'
WRITER='''mata:
void savebin(string scalar file,real matrix A) {
 real matrix B
 real scalar f
 B=bufio();bufbyteorder(B,2);f=fopen(file,"w")
 fbufput(B,f,"%4b",rows(A));fbufput(B,f,"%4b",cols(A));fbufput(B,f,"%8z",A);fclose(f)
}
end
'''

def readbin(p):
    b=p.read_bytes();r,c=struct.unpack('<ii',b[:8]);a=np.frombuffer(b[8:],dtype='<f8').reshape(r,c).copy()
    a[a>8e307]=np.nan;return a


def same_bits(a,b):
    # Accepted solution fields must be finite. Byte comparison also detects
    # signed-zero differences that ordinary numeric equality would miss.
    assert a.shape==b.shape and np.isfinite(a).all() and np.isfinite(b).all()
    assert a.tobytes()==b.tobytes()


def check_trace(t,tol):
    assert t.shape[1] in (75,150);counts=np.zeros(3,dtype=int)
    for k in range(len(t)):
        for start in range(0,t.shape[1],25):
            final=t[k,start];h=t[k,start+1:start+25].reshape(4,6)
            if np.isnan(final): assert np.isnan(h).all();continue
            assert np.isfinite(h[0,0]) and h[0,2]==tol
            used=np.isfinite(h[:,0]);assert 1<=used.sum()<=4
            for v in range(4):
                if not used[v]:assert np.isnan(h[v]).all();continue
                assert h[v,2]==(min(tol,1e-12) if v==3 else tol)
                assert h[v,3]>=0
                if h[v,1]==0:
                    assert h[v,0]==0 and h[v,4]<=tol and np.isfinite(h[v,5])
            if final==0:
                assert h[np.where(used)[0][-1],1]==0
                assert not np.any(h[:,0]==2)
                if h[0,1]==0: assert used.sum()==1
            counts[0]+=used[1:].sum()
            counts[1]+=int(final==0 and h[0,1]!=0)
            counts[2]+=int(final==91)
    return counts


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--baseline-ado',type=Path,required=True);ap.add_argument('--only')
    args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    pkg=out/'ado';pkg.mkdir();hashes={}
    for n in package_files(REPO):
        raw=(REPO/n).read_bytes();(pkg/Path(n).name).write_bytes(raw);hashes[n]=hashlib.sha256(raw).hexdigest()
    (out/'manifest.json').write_text(json.dumps(hashes,indent=2)+'\n');results=[]
    baseline=args.baseline_ado.resolve()
    (out/'baseline_manifest.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in baseline.iterdir() if p.is_file()},indent=2)+'\n')

    def run(name,body,source=pkg,check=None):
        if args.only and not re.search(args.only,name):return None
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 120\n'
        for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:
            (d/key).mkdir();head+=f'sysdir set {key} "{d/key}"\n'
        head+=f'adopath ++ "{source}"\n'+WRITER
        (d/'test.do').write_text(head+body+'\ndi "TEST_PASS"\n');row=dict(test=name,passed=False);start=time.monotonic()
        try:
            p=subprocess.run([args.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=240)
            log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else ''
            assert p.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-2400:] or p.stderr.decode(errors='replace')
            if check:row.update(check(d))
            row['passed']=True
        except Exception as e:row['error']=repr(e)
        row['seconds']=round(time.monotonic()-start,3);results.append(row)
        (out/'results.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(row),flush=True);return d

    def clone(name):
        dst=out/('source_'+name);shutil.copytree(pkg,dst);return dst

    # Pair both stages against the immediately preceding reviewed snapshot.
    # Only initially accepted LPs are equality controls; recovered rows may differ.
    for fixture,rts,ori,tol in itertools.product(FIXTURES,['crs','vrs','drs','irs'],['input','output'],[1e-7,1e-10]):
        name=f'parity_{fixture}_{rts}_{ori}_{tol}'
        if args.only and not re.search(args.only,name):continue
        def make(new):
            s=FIXTURES[fixture]
            qual='if id>=2 in 2/8' if fixture=='duplicate' else ''
            for stage in [1,2]:
                spec=f'x1 x2=y1 y2 {qual},' if stage==1 else f'y1 y2 {qual}, inputs(x1 x2)'
                s+=f'dea {spec} dmu(id) engine(ipm) rts({rts}) ort({ori}) stage({stage}) tolerance({tol}) gen(e{stage}) reference replacevars\n'
                mats=['efficiency','lambda']+(['slack_in','slack_out'] if stage==2 else [])
                if rts=='vrs':mats+=['vrsfrontier']
                if new:
                    mats+=['stage_status','ipm_trace']
                    s+=f'matrix C=(r(n_ipm_retries),r(n_ipm_recovered),r(n_ipm_conflicts))\nmata: savebin("counts{stage}.bin",st_matrix("C"))\n'
                    if rts=='vrs':
                        mats+=['vrs_ipm_trace'];s+=f'matrix C=(r(n_vrs_ipm_retries),r(n_vrs_ipm_recovered),r(n_vrs_ipm_conflicts))\nmata: savebin("vc{stage}.bin",st_matrix("C"))\n'
                for mat in mats:s+=f'mata: savebin("{mat}{stage}.bin",st_matrix("r({mat})"))\n'
            return s
        b=run(name+'_before',make(False),baseline)
        def parity(d,b=b,rts=rts,tol=tol):
            first=0;projected=0;recovered=0;unavailable=0
            for stage in [1,2]:
                t=readbin(d/f'ipm_trace{stage}.bin');st=readbin(d/f'stage_status{stage}.bin')
                np.testing.assert_array_equal(check_trace(t,tol),readbin(d/f'counts{stage}.bin')[0])
                eff=readbin(d/f'efficiency{stage}.bin');old=readbin(b/f'efficiency{stage}.bin')
                mask=t[:,2]==0;first+=int(mask.sum());same_bits(eff[mask],old[mask])
                assert np.array_equal(st[:,0],np.isnan(eff[:,0]).astype(float))
                recovered+=int(((t[:,0]==0)&(t[:,2]!=0)).sum())
                if stage==1:
                    assert np.isnan(t[:,25:50]).all()
                    same_bits(readbin(d/'lambda1.bin')[mask],readbin(b/'lambda1.bin')[mask])
                else:
                    np.testing.assert_array_equal(t[:,25],st[:,1])
                    mask=mask&(t[:,27]==0);projected+=int(mask.sum())
                    for mat in ['lambda','slack_in','slack_out']:
                        same_bits(readbin(d/f'{mat}2.bin')[mask],readbin(b/f'{mat}2.bin')[mask])
                    recovered+=int(((t[:,25]==0)&(t[:,27]!=0)).sum())
                    unavailable+=int((st[:,1]!=0).sum())
                if rts=='vrs':
                    v=readbin(d/f'vrs_ipm_trace{stage}.bin')
                    np.testing.assert_array_equal(check_trace(v,tol),readbin(d/f'vc{stage}.bin')[0])
                    np.testing.assert_array_equal(t[:,:25],v[:,50:75])
                    vm=np.all(v[:,[2,52,102]]==0,axis=1)
                    same_bits(readbin(d/f'vrsfrontier{stage}.bin')[vm],readbin(b/f'vrsfrontier{stage}.bin')[vm])
            return dict(initial_radial_bitwise=first,initial_projection_bitwise=projected,recovered_phases=recovered,unavailable_projections=unavailable)
        run(name+'_after',make(True),check=parity)

    # Initial-only failures exercise real recovery, with unchanged native calls.
    for phase in ['primal','stage2','auxiliary']:
        target=clone(phase);p=target/'_dea_ipm.ado';s=p.read_text()
        begin={'primal':'void _dea_ipm_envelopment(','stage2':'real scalar _dea_ipm_slacks(','auxiliary':'void _dea_ipm_vrsfrontier('}[phase]
        end={'primal':'real scalar _dea_ipm_slacks(','stage2':'real scalar _dea_ipm_dual_input(','auxiliary':'// 5. UTILITY:'}[phase]
        lo=s.index(begin);hi=s.index(end,lo);part=s[lo:hi]
        anchor={'primal':'if (S.converged() && S.errorcode() == 0)','stage2':'if (S2.converged() && S2.errorcode() == 0)','auxiliary':'if (S.converged() & S.errorcode()==0)'}[phase]
        assert part.count(anchor)==(2 if phase=='primal' else 1)
        part=part.replace(anchor,anchor.replace('if (','if (0 & ',1));p.write_text(s[:lo]+part+s[hi:])
        for rts,ori,stage in itertools.product(['crs','vrs','drs','irs'],['input','output'],[1,2]):
            if phase=='auxiliary' and rts!='vrs':continue
            answer=(.25 if rts in ['crs','irs'] else .5) if ori=='output' else (.25 if rts in ['crs','drs'] else 1/3)
            body=SIMPLE+f'dea x=y, dmu(id) engine(ipm) rts({rts}) ort({ori}) stage({stage}) tolerance(1e-7) gen(eff) reference\n'
            body+=f'assert r(n_radial_failed)==0 & abs(eff[3]-({answer}))<1e-6\n'
            if phase=='auxiliary':
                body+='assert r(n_vrs_ipm_retries)>0 & r(n_vrs_ipm_recovered)>0\nmata: savebin("trace.bin",st_matrix("r(vrs_ipm_trace)"))\n'
            else:
                if phase=='primal' or stage==2:body+='assert r(n_ipm_retries)>0 & r(n_ipm_recovered)>0\n'
                body+='mata: savebin("trace.bin",st_matrix("r(ipm_trace)"))\n'
            run(f'initial_fault_{phase}_{rts}_{ori}_{stage}',body,target,check=lambda d:dict(trace_counts=check_trace(readbin(d/'trace.bin'),1e-7).tolist()))

    # Dual solver entry points preserve optional-argument compatibility and
    # reciprocal output scores after a real native retry.
    target=clone('dual');p=target/'_dea_ipm.ado';s=p.read_text()
    lo=s.index('real scalar _dea_ipm_dual_input(');hi=s.index('void _dea_ipm_vrsfrontier(',lo);part=s[lo:hi]
    assert part.count('if (lp.converged() == 1 & lp.errorcode() == 0)')==2
    part=part.replace('if (lp.converged() == 1 & lp.errorcode() == 0)','if (0 & lp.converged() == 1 & lp.errorcode() == 0)');p.write_text(s[:lo]+part+s[hi:])
    for rts,ori in itertools.product(['crs','vrs','drs','irs'],['input','output']):
        answer=(.25 if rts in ['crs','irs'] else .5) if ori=='output' else (.25 if rts in ['crs','drs'] else 1/3)
        body=f'''run "{target/'_dea_ipm.ado'}"
mata:
audit=J(1,25,.)
e=_dea_ipm_dual_{ori}((1\\2\\4),(1\\4\\2),3,"{rts}",1e-7,audit)
assert(abs(e-({answer}))<1e-6 & audit[1]==0 & audit[3]!=0)
savebin("trace.bin",(audit,J(1,50,.)))
end
'''
        run(f'dual_retry_{rts}_{ori}',body,target,check=lambda d:dict(trace_counts=check_trace(readbin(d/'trace.bin'),1e-7).tolist()))

    # Contradictory primal reports must not disappear through dual fallback.
    target=clone('conflict');p=target/'_dea_ipm.ado';s=p.read_text();lo=s.index('void _dea_ipm_envelopment(');hi=s.index('real scalar _dea_ipm_slacks(',lo);part=s[lo:hi]
    part=part.replace('primal_status=_dea_ipm_finish(primal_status,S.errorcode(),','primal_ok=0;eff_score=.;lam_k=J(1,n,.)\n        primal_status=_dea_ipm_finish(2,2,');p.write_text(s[:lo]+part+s[hi:])
    for stage in [1,2]:
        body=SIMPLE+f'''dea x=y, dmu(id) engine(ipm) rts(vrs) stage({stage}) tolerance(1e-7) gen(eff)
assert r(n_radial_failed)==3 & r(n_ipm_conflicts)==3 & r(n_failed)==3
assert missing(eff)
matrix T=r(ipm_trace)
mata: assert(all(st_matrix("T")[.,1]:==91) & all(st_matrix("T")[.,51..75]:>=.))
''';run(f'conflict_{stage}',body,target)

    # Even if an unbounded report is followed only by rejected candidates,
    # a finite dual result must not hide that unresolved primal evidence.
    blocked=out/'source_unbounded_rejected';shutil.copytree(target,blocked)
    p=blocked/'_dea_native_retry.ado';text=p.read_text().replace('x=S.parameters();d=', 'x=S.parameters();x[1]=x[1]+100;d=');p.write_text(text)
    run('unbounded_then_rejected',SIMPLE+'''dea x=y, dmu(id) engine(ipm) rts(crs) stage(1) tolerance(1e-7) gen(eff)
assert r(n_radial_failed)==3 & r(n_failed)==3
assert missing(eff)
matrix T=r(ipm_trace)
mata: assert(any(st_matrix("T")[.,1]:==90) & all(st_matrix("T")[.,51..75]:>=.))
''',blocked)

    # r() scope and row alignment after nested routines, I/O and repeated calls.
    auto=out/'source_auto';auto_instrument(pkg,auto)
    for mode in ['saving','bootstrap','contractcheck','auto','small_matsize']:
        for stage in [1,2]:
            source=auto if mode=='auto' else pkg
            extras={'saving':'saveref saving("result.dta")','bootstrap':'bootstrap reps1(50)','contractcheck':'contractcheck','auto':'','small_matsize':''}[mode]
            eng='auto' if mode=='auto' else 'ipm'
            body=SIMPLE+'gen byte useit=id!=2\n'+('set matsize 10\n' if mode=='small_matsize' else '')
            body+=f'dea x=y if useit, dmu(id) engine({eng}) rts(vrs) stage({stage}) tolerance(1e-7) gen(a) reference {extras}\n'
            body+='matrix T=r(ipm_trace)\nmatrix V=r(vrs_ipm_trace)\nmatrix C=(r(n_ipm_retries),r(n_ipm_recovered),r(n_ipm_conflicts),r(n_vrs_ipm_retries),r(n_vrs_ipm_recovered),r(n_vrs_ipm_conflicts))\n'
            body+='assert rowsof(T)==2 & colsof(T)==75 & rowsof(V)==2 & colsof(V)==150\nassert missing(a)==!useit\n'
            if mode=='small_matsize':body+='assert c(matsize)==10\n'
            body+=f'dea y if useit, inputs(x) dmu(id) engine(ipm) rts(vrs) stage({stage}) tolerance(1e-7) gen(b) reference replacevars\n'
            body+='matrix U=r(ipm_trace)\nmatrix W=r(vrs_ipm_trace)\nmatrix D=(r(n_ipm_retries),r(n_ipm_recovered),r(n_ipm_conflicts),r(n_vrs_ipm_retries),r(n_vrs_ipm_recovered),r(n_vrs_ipm_conflicts))\nmata: assert(all(st_matrix("T"):==st_matrix("U")) & all(st_matrix("V"):==st_matrix("W")) & all(st_matrix("C"):==st_matrix("D")))\nassert a==b\n'
            body+='dea x=y if useit, dmu(id) engine(simplex) rts(crs) gen(c) replacevars\ncapture matrix list r(ipm_trace)\nassert _rc==111\ncapture matrix list r(vrs_ipm_trace)\nassert _rc==111\n'
            run(f'preservation_{mode}_{stage}',body,source)

    # A previous engine may already be in Mata memory when the updated ado
    # entry point is loaded. Its old optional-argument API must not leak through.
    for rts,stage in itertools.product(['crs','vrs'],[1,2]):
        body=f'run "{baseline/"_dea_ipm.ado"}"\n'+SIMPLE
        body+=f'dea x=y, dmu(id) engine(ipm) rts({rts}) stage({stage}) tolerance(1e-7) gen(eff)\n'
        body+='assert r(n_radial_failed)==0\nmatrix T=r(ipm_trace)\nassert rowsof(T)==3 & colsof(T)==75\n'
        if stage==1:body+='mata: assert(all(st_matrix("T")[.,26..50]:>=.))\n'
        if rts=='vrs':body+='matrix V=r(vrs_ipm_trace)\nassert rowsof(V)==3 & colsof(V)==150\n'
        run(f'cached_old_engine_{rts}_{stage}',body)

    summary=dict(cases=len(results),failed=[r['test'] for r in results if not r['passed']],
                 initial_radial_bitwise=sum(r.get('initial_radial_bitwise',0) for r in results),
                 initial_projection_bitwise=sum(r.get('initial_projection_bitwise',0) for r in results),
                 recovered_phases=sum(r.get('recovered_phases',0) for r in results),
                 unavailable_projections=sum(r.get('unavailable_projections',0) for r in results))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');assert not summary['failed'],summary

if __name__=='__main__':main()
