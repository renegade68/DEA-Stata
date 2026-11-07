#!/usr/bin/env python3
"""TE returned-projection diagnostics and small positive score preservation.
Requires Stata 16+, NumPy, and the immediately preceding R06 ado snapshot.
This suite checks numerical residuals, not exact LP feasibility/optimality.
"""
from _support import package_files
from pathlib import Path
import argparse, hashlib, itertools, json, re, shutil, subprocess, time
import numpy as np
from test_core_retries import WRITER, readbin
from test_simplex_stages import FIXTURES as SIMPLE, INIT, instrument as simplex_fault
from test_ipm_stages import FIXTURES as IPM, instrument as ipm_fault
from test_saveref import instrument as auto_fault
REPO=Path(__file__).resolve().parents[1]
FIXTURES={k:SIMPLE[k] for k in ['slack','random','scaled','duplicate']}
FIXTURES['collinear']=SIMPLE['slack']+'replace x2=2*x1\nreplace y2=3*y1\n'
FIXTURES['boundary']=IPM['p03_vrs']
FIELDS='radial_initial radial_stage2 radial_adjustment stage2_status checked within_tolerance input_abs output_abs rts_abs lambda_sum lambda_gap lambda_bound slack_in_bound slack_out_bound input_scaled output_scaled rts_scaled lambda_scaled slack_in_scaled slack_out_scaled max_scaled input_normalized output_normalized'

def export(stage,diagnostics=True):
    names=['efficiency','lambda','stage_status','solution_source']
    if stage==2:names+=['slack_in','slack_out']
    s=''.join(f'mata: savebin("{n}.bin",st_matrix("r({n})"))\n' for n in names)
    s+='capture matrix list r(ipm_trace)\nif !_rc {\nmata: savebin("ipm_trace.bin",st_matrix("r(ipm_trace)"))\n}\n'
    if diagnostics:
        s+='matrix D=r(numerical_diagnostics)\nlocal dn : colnames D\n'
        s+=f'assert "`dn\'"=="{FIELDS}"\n'
        s+='mata: savebin("diagnostics.bin",st_matrix("D"))\n'
        s+='matrix C=(r(n_projection_checked),r(n_projection_outside_tol),r(n_projection_unchecked))\nmata: savebin("counts.bin",st_matrix("C"))\n'
    s+='mata: savebin("data.bin",st_data(.,tokens("x1 x2 y1 y2"),"selected"))\n'
    s+='assert sentinel=="preserve me"\nassert missing(e) if !selected\n'
    return s

def audit(d,ori,rts,tol,stage):
    q=readbin(d/'diagnostics.bin');e=readbin(d/'efficiency.bin')[:,0];st=readbin(d/'stage_status.bin');data=readbin(d/'data.bin');l=readbin(d/'lambda.bin');n=len(e)
    assert q.shape==(n,23)
    np.testing.assert_array_equal(q[:,3],st[:,1]);np.testing.assert_array_equal(np.isnan(q[:,0]),np.isnan(e))
    np.testing.assert_allclose(q[:,0],e if ori=='input' else 1/e,rtol=2e-15,atol=0,equal_nan=True)
    use=(st[:,1]==0);np.testing.assert_array_equal(np.isfinite(q[:,1]),use)
    np.testing.assert_array_equal(q[use,2],q[use,1]-q[use,0])
    assert np.isnan(q[~use,5:]).all() and np.all(q[~use,4]==0)
    np.testing.assert_array_equal(readbin(d/'counts.bin')[0],[sum(q[:,4]),sum(q[:,5]==0),sum(use & (q[:,4]==0))])
    if stage==1:assert not use.any();return dict(checked_rows=0,outside_rows=0,max_scaled=0.)
    xs=np.max(abs(data[:,:2]),axis=0);ys=np.max(abs(data[:,2:]),axis=0);xs[xs==0]=1;ys[ys==0]=1
    x=data[:,:2]/xs;y=data[:,2:]/ys;sin=readbin(d/'slack_in.bin');sout=readbin(d/'slack_out.bin')
    for k in np.where(use)[0]:
        assert q[k,4]==1
        si=sin[k]/xs;so=sout[k]/ys;lam=l[k];rad=q[k,1]
        bx=rad*x[k] if ori=='input' else x[k];by=y[k] if ori=='input' else rad*y[k]
        rx=lam@x+si-bx;ry=lam@y-so-by
        dx=np.maximum(abs(bx),abs(x).T@np.maximum(1,abs(lam))+np.maximum(1,abs(si)))
        dy=np.maximum(abs(by),abs(y).T@np.maximum(1,abs(lam))+np.maximum(1,abs(so)))
        rawx=lam@data[:,:2]+sin[k]-(rad*data[k,:2] if ori=='input' else data[k,:2])
        rawy=lam@data[:,2:]-sout[k]-(data[k,2:] if ori=='input' else rad*data[k,2:])
        total=sum(lam);gap=total-1
        rt=abs(gap) if rts=='vrs' else max(0,gap) if rts=='drs' else max(0,-gap) if rts=='irs' else 0
        bounds=[max(0,-min(v)) for v in [lam,sin[k],sout[k]]]
        sb=[max(np.maximum(0,-v)/np.maximum(1,abs(v))) for v in [lam,si,so]]
        scaled=[max(abs(rx)/dx),max(abs(ry)/dy),rt/max(1,sum(np.maximum(1,abs(lam))))]+sb
        # Direct original-coordinate arithmetic differs from scaled-row evaluation
        # by rounding; use a forward rounding allowance for these absolute residuals.
        for col,value,scale in [(6,max(abs(rawx)),max(abs(lam)@abs(data[:,:2])+abs(sin[k])+abs(rad*data[k,:2]))),(7,max(abs(rawy)),max(abs(lam)@abs(data[:,2:])+abs(sout[k])+abs(rad*data[k,2:])) )]:
            assert abs(q[k,col]-value)<=2e-14*max(1,scale),(k,col,q[k,col],value,scale)
        np.testing.assert_allclose(q[k,8:14],[rt,total,gap]+bounds,rtol=1e-12,atol=3e-14)
        np.testing.assert_allclose(q[k,14:21],scaled+[max(scaled)],rtol=2e-8,atol=3e-15)
        np.testing.assert_allclose(q[k,21:23],[max(abs(rx)),max(abs(ry))],rtol=2e-8,atol=3e-14)
        assert q[k,5]==(q[k,20]<=tol)
    return dict(checked_rows=int(use.sum()),outside_rows=int(sum(q[:,5]==0)),max_scaled=float(max(q[use,20],default=0.)))

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--baseline-ado',type=Path,required=True);ap.add_argument('--only');args=ap.parse_args()
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);pkg=out/'ado';pkg.mkdir();hashes={}
    for n in package_files(REPO):
        p=REPO/n;(pkg/p.name).write_bytes(p.read_bytes());hashes[n]=hashlib.sha256(p.read_bytes()).hexdigest()
    baseline=args.baseline_ado.resolve();(out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes,baseline_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in baseline.iterdir() if p.is_file()}),indent=2)+'\n')
    sf=out/'simplex_fault';simplex_fault(pkg,sf);ipf=out/'ipm_fault';ipm_fault(pkg,ipf);af=out/'auto_fault';auto_fault(pkg,af);rows=[]
    def run(name,body,source=pkg,check=None):
        if args.only and not re.search(args.only,name):return None
        d=out/name;d.mkdir();head='version 16\nclear all\nset more off\nset linesize 200\n'
        for key in ['PLUS','PERSONAL','SITE','OLDPLACE']:(d/key).mkdir();head+=f'sysdir set {key} "{d/key}"\n'
        (d/'test.do').write_text(head+f'adopath ++ "{source}"\n'+INIT+WRITER+body+'\ndi "TEST_PASS"\n');t=time.monotonic();row=dict(test=name,passed=False)
        try:
            p=subprocess.run([args.stata,'-b','do','test.do'],cwd=d,capture_output=True,timeout=240)
            log=(d/'test.log').read_text(errors='replace') if (d/'test.log').exists() else p.stderr.decode(errors='replace')
            assert p.returncode==0 and 'TEST_PASS' in log.splitlines(),log[-3200:]
            if check:row.update(check(d))
            row['passed']=True
        except Exception as e:row['error']=repr(e)
        row['seconds']=round(time.monotonic()-t,3);rows.append(row);(out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(row),flush=True);return d
    # Unchanged accepted scores, references, slacks, statuses, and native traces.
    for fixture,engine,rts,ori,tol,stage in itertools.product(FIXTURES,['simplex','ipm','auto'],['crs','vrs','drs','irs'],['input','output'],[1e-7,1e-10],[1,2]):
        name=f'grid_{fixture}_{engine}_{rts}_{ori}_{tol}_{stage}'
        if args.only and not re.search(args.only,name):continue
        qual='if id>=2 in 2/8' if fixture=='duplicate' else '';condition='id>=2 & inrange(_n,2,8)' if qual else '1'
        body=FIXTURES[fixture]+f'gen byte selected={condition}\n'
        opts=f'dmu(id) engine({engine}) rts({rts}) ort({ori}) tolerance({tol}) stage({stage}) reference gen(e)'
        old=run(name+'_before',body+f'dea x1 x2=y1 y2 {qual}, {opts}\n'+export(stage,False),baseline)
        def check(d,old=old,ori=ori,rts=rts,tol=tol,stage=stage):
            files=[p.name for p in old.glob('*.bin')];assert all((old/f).read_bytes()==(d/f).read_bytes() for f in files),[f for f in files if (old/f).read_bytes()!=(d/f).read_bytes()]
            z=audit(d,ori,rts,tol,stage);z.update(historical_pair=True,exact_numeric_files=len(files));return z
        run(name+'_after',body+f'dea y1 y2 {qual}, inputs(x1 x2) {opts}\n'+export(stage),check=check)
    # Optional recovery must report the actual RHS, without replacing initial scores.
    for rts,ori,tol,band in itertools.product(['crs','vrs','drs','irs'],['input','output'],[1e-7,1e-10],[1e-7,1e-6]):
        body=FIXTURES['boundary']+'gen byte selected=1\n'
        body+=f'dea x1 x2=y1 y2,dmu(id) engine(ipm) rts({rts}) ort({ori}) tolerance({tol}) radialtolerance({band}) reference gen(e)\n'
        body+='matrix P=r(projection)\nmata: assert(all(st_matrix("P")[.,1..3]:==st_matrix("r(numerical_diagnostics)")[.,1..3]))\n'+export(2)
        run(f'recovery_{rts}_{ori}_{tol}_{band}',body,check=lambda d,ori=ori,rts=rts,tol=tol:audit(d,ori,rts,tol,2))
    # Diagnostics exist even without reference; they agree with an explicit request.
    for engine,ori,rts in itertools.product(['simplex','ipm','auto'],['input','output'],['ccr','bcc','nirs','ndrs']):
        s=FIXTURES['slack']+f'dea x1 x2=y1 y2,dmu(id) engine({engine}) ort({ori}) rts({rts}) gen(e)\n'
        s+='matrix D=r(numerical_diagnostics)\nscalar n=r(n_projection_checked)\n'
        s+=f'dea x1 x2=y1 y2,dmu(id) engine({engine}) ort({ori}) rts({rts}) gen(e1) reference replacevars\nmata: assert(all(st_matrix("D"):==st_matrix("r(numerical_diagnostics)")))\nassert e==e1\nassert n==r(n_projection_checked)\n'
        run(f'no_reference_{engine}_{ori}_{rts}',s)
    # Failure, dual-only score, and skipped stage 2 remain distinguishable.
    for engine,ori,stage,mode in itertools.product(['simplex','ipm'],['input','output'],[1,2],['primal','both','stage2']):
        s=FIXTURES['slack']+'gen byte selected=1\n'
        s+=f'scalar fault_primal={int(mode in ["primal","both"])}\nscalar fault_dual={int(mode=="both")}\nscalar fault_stage2={int(mode=="stage2")}\n'
        s+=f'dea x1 x2=y1 y2,dmu(id) engine({engine}) ort({ori}) stage({stage}) reference gen(e)\n'+export(stage)
        run(f'fault_{engine}_{ori}_{stage}_{mode}',s,sf if engine=='simplex' else ipf,check=lambda d,ori=ori,stage=stage:audit(d,ori,'crs',1e-7,stage))
    # Deliberate small-score correction: quantify the change, including failed stage 2.
    for ori,power,stage in itertools.product(['input','output'],[8,9,10,11],[1,2]):
        values=f'1 1 1\n2 1e{power} 1' if ori=='input' else f'1 1 1e{power}\n2 1 1'
        s='input double id x y\n'+values+'\nend\n'
        cmd=f'dea x=y,dmu(id) engine(simplex) ort({ori}) stage({stage}) tolerance(1e-12) gen(e) reference\n'
        run(f'small_before_{ori}_{power}_{stage}',s+cmd+('assert e[2]==0\n' if power==9 or (ori=='output' and power>8) else f'assert abs(e[2]-1e-{power})<1e-18\n'),baseline)
        s+=cmd+f'assert e[2]>0 & abs(e[2]-1e-{power})<1e-18\n'
        s+='matrix D=r(numerical_diagnostics)\nassert D[2,1]>0\n'
        if stage==2:
            s+='scalar fault_stage2=1\n'+cmd.replace('gen(e)','gen(f)')+'assert e==f\nassert r(n_stage2_failed)==2\n'
        run(f'small_after_{ori}_{power}_{stage}',s,sf)
    # A presentation cleanup may not destroy a solver-accepted primal vector.
    tiny={
        'lambda':'input double id x1 x2 y1 y2\n1 1 1 1 1\n2 1 1 1e-9 1e-9\nend\n',
        'slack_in':'input double id x1 x2 y1 y2\n1 1 1 1 1\n2 2 2.000000002 1 1\nend\n',
        'slack_out':'input double id x1 x2 y1 y2\n1 1 1 1 1.000000001\n2 2 2 1 1\nend\n',
    }
    for fixture,rts,tol,stage in itertools.product(tiny,['crs','vrs','drs','irs'],[1e-7,1e-12],[1,2]):
        name=f'cleanup_{fixture}_{rts}_{tol}_{stage}'
        if args.only and not re.search(args.only,name):continue
        body=tiny[fixture]+'gen byte selected=1\ngen str20 sentinel="preserve me"\n'
        cmd=f'dea x1 x2=y1 y2,dmu(id) engine(simplex) rts({rts}) ort(input) stage({stage}) tolerance({tol}) reference gen(e)\n'
        old=run(name+'_before',body+cmd+export(stage,False),baseline)
        def check(d,old=old,rts=rts,tol=tol,stage=stage,fixture=fixture):
            z=audit(d,'input',rts,tol,stage);assert z['outside_rows']==0
            cols=['lambda']+(['slack_in','slack_out'] if stage==2 else [])
            a=np.column_stack([readbin(d/(n+'.bin')) for n in cols]);b=np.column_stack([readbin(old/(n+'.bin')) for n in cols]);changed=np.any(~((a==b)|(np.isnan(a)&np.isnan(b))),axis=1)
            np.testing.assert_array_equal(readbin(d/'stage_status.bin'),readbin(old/'stage_status.bin'))
            if fixture=='lambda' and rts in ['crs','drs'] and tol==1e-12:
                assert readbin(d/'lambda.bin')[1,0]>0 and readbin(old/'lambda.bin')[1,0]==0
                assert readbin(d/'efficiency.bin')[1,0]>0
            if tol==1e-7:assert not changed.any()
            z.update(cleanup_case=True,restored_reference_rows=int(changed.sum()));return z
        run(name+'_after',body+cmd+export(stage),check=check)
    # Loading an old engine first must not retain an obsolete Mata signature.
    for engine in ['simplex','ipm']:
        s=f'run "{baseline/("_dea_"+engine+".ado")}"\n'+FIXTURES['slack']+f'dea x1 x2=y1 y2,dmu(id) engine({engine}) gen(e)\n'
        s+='matrix D=r(numerical_diagnostics)\nassert colsof(D)==23\n'
        run('cached_'+engine,s)
    # A non-TE result must not inherit diagnostics from a previous TE command.
    for model,opt in [('fdh',''),('sbm',''),('congestion',''),('sbmundesirable','bad(x2)')]:
        s=FIXTURES['slack']+'dea x1 x2=y1 y2,dmu(id) gen(e)\n'
        s+=f'dea x1=y1 y2,dmu(id) model({model}) {opt} gen(f)\ncapture matrix list r(numerical_diagnostics)\nassert _rc==111\nassert missing(r(n_projection_checked))\n'
        run('no_stale_'+model,s)
    # Saving/using, VRS auxiliaries and bootstrap must preserve the main audit.
    for engine,ori,rts,kind in itertools.product(['simplex','ipm','auto'],['input','output'],['crs','vrs'],['saving','bootstrap']):
        cmd=f'dea x1 x2=y1 y2,dmu(id) engine({engine}) ort({ori}) rts({rts}) reference'
        body=FIXTURES['slack']+cmd+' gen(control)\nmatrix D=r(numerical_diagnostics)\n'
        extra='bootstrap reps1(60)' if kind=='bootstrap' else ''
        body+='set seed 5297\n'+cmd+f' gen(e) replacevars {extra} saving("saved.dta")\n'
        body+='assert e==control\nmata: assert(all(st_matrix("D"):==st_matrix("r(numerical_diagnostics)")))\n'
        body+='preserve\nuse "saved.dta",clear\nassert _N==rowsof(D)\n'
        body+='capture confirm variable numerical_diagnostics\nassert _rc==111\nrestore\n'
        run(f'integration_{engine}_{ori}_{rts}_{kind}',body)
    # Force auto fallback in a private dispatcher copy; compare with explicit IPM.
    for ori,rts,stage in itertools.product(['input','output'],['crs','vrs'],[1,2]):
        cmd=f'dea x1 x2=y1 y2,dmu(id) ort({ori}) rts({rts}) stage({stage}) reference'
        body=FIXTURES['slack']+cmd+' engine(ipm) gen(control)\nmatrix D=r(numerical_diagnostics)\n'
        body+=cmd+' engine(auto) gen(e) replacevars\nassert "`r(engine)\'"=="ipm"\n'
        body+='assert e==control\nmata: assert(all(st_matrix("D"):==st_matrix("r(numerical_diagnostics)")))\n'
        run(f'integration_fallback_{ori}_{rts}_{stage}',body,af)
    # Direct malformed/perturbed returned candidates exercise the audit, without
    # claiming that a production solver actually returned these vectors.
    for mode,setup,check in [
        ('valid','','assert(D[5]==1 & D[6]==1 & D[21]==0)'),
        ('input','si=.01','assert(D[5]==1 & D[6]==0 & abs(D[7]-10)<1e-10 & D[15]>0)'),
        ('output','so=.01','assert(D[5]==1 & D[6]==0 & abs(D[8]-.001)<1e-14 & D[16]>0)'),
        ('lambda','lam=(-.1,1.1)','assert(D[5]==1 & D[6]==0 & D[12]==.1 & D[18]>0)'),
        ('slack_in','si=-.01','assert(D[5]==1 & D[6]==0 & D[13]==10 & D[19]>0)'),
        ('slack_out','so=-.01','assert(D[5]==1 & D[6]==0 & abs(D[14]-.001)<1e-14 & D[20]>0)'),
        ('missing','lam=(.,0)','assert(D[5]==0 & D[6]>=.)'),
        ('overflow','XR=J(2,1,1e308);radial=10','assert(D[5]==0 & D[6]>=.)'),
        ('failed','status=3','assert(D[5]==0 & D[6]>=.)'),
        ('skipped','status=.;radial=.','assert(D[5]==0 & D[6]>=.)'),
    ]:
        for ori in ['input','output']:
            # Overflow the expanded side for either orientation.
            setup=setup.replace('XR=J(2,1,1e308)','XR=YR=J(2,1,1e308)')
            body=f'run "{pkg/"_dea_ipm.ado"}"\nmata:\n'
            body+='X=(1\\2);Y=(1\\1);XR=X:*1000;YR=Y:*.1;lam=(1,0);si=so=0;radial=1;status=0\n'
            body+=setup+'\n'+f'D=_dea_teaudit(X,Y,XR,YR,1,"{ori}","vrs",1000,.1,lam,si,so,1,radial,status,1e-7)\n'+check+'\nend\n'
            run('helper_'+mode+'_'+ori,body)
    for rts,total in itertools.product(['crs','vrs','drs','irs'],[.5,1.,1.5]):
        violation=abs(total-1) if rts=='vrs' else max(0,total-1) if rts=='drs' else max(0,1-total) if rts=='irs' else 0
        body=f'run "{pkg/"_dea_ipm.ado"}"\nmata:\n'
        body+=f'X=Y=(1\\1);D=_dea_teaudit(X,Y,X,Y,1,"input","{rts}",1,1,({total},0),0,0,1,1,0,1e-7)\n'
        body+=f'assert(D[9]=={violation} & D[10]=={total} & D[11]=={total-1})\nassert(abs(D[17]-{violation}/max((1,max((1,{total}))+1)))<1e-15)\nend\n'
        run(f'helper_rts_{rts}_{total}',body)
    summary=dict(cleanup_cases=sum(r.get('cleanup_case',False) for r in rows),restored_reference_rows=sum(r.get('restored_reference_rows',0) for r in rows),total=len(rows),passed=sum(r['passed'] for r in rows),historical_pairs=sum(r.get('historical_pair',False) for r in rows),exact_numeric_files=sum(r.get('exact_numeric_files',0) for r in rows),checked_rows=sum(r.get('checked_rows',0) for r in rows),outside_rows=sum(r.get('outside_rows',0) for r in rows),max_scaled=max((r.get('max_scaled',0) for r in rows),default=0))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
    raise SystemExit(not rows or summary['total']!=summary['passed'])
if __name__=='__main__':main()
