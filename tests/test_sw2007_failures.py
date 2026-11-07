#!/usr/bin/env python3
"""Fault-injection checks of whole-replication rejection and bounded retry.

All changes are confined to private test copies. The numerical kernel is
unaltered except for explicitly invalidating one result or a regression sample.
"""
from _support import package_files
import argparse,difflib,hashlib,json,subprocess
from pathlib import Path
from test_sw2007 import DATA,REPO
from test_output_protection import unchanged


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stata',required=True);ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    names=package_files(REPO)
    hashes={n:hashlib.sha256((REPO/n).read_bytes()).hexdigest() for n in names}
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256=hashes),indent=2));results=[]
    for stage in ['dea','regression']:
        for permanent in [False,True]:
            name=f'{stage}_{"permanent" if permanent else "recover"}';work=out/name;work.mkdir();pkg=work/'ado';pkg.mkdir()
            for n in names:(pkg/Path(n).name).write_bytes((REPO/n).read_bytes())
            f=pkg/'_dea_sw2007_seq.ado';original=f.read_text()
            if stage=='dea':
                needle='    st_store(ix,dest,result)'
                cond='1' if permanent else 'st_numscalar("TEST_ATTEMPTS")<3'
                injection='    st_numscalar("TEST_ATTEMPTS",st_numscalar("TEST_ATTEMPTS")+1)\n'+f'    if ({cond}) result[1]=.\n'
            else:
                needle="        capture qui truncreg `draw' `envars' if `touse', ll(`tau') iterate(`maxiter')"
                cond='1' if permanent else "`reg_try'<=2"
                injection=f"        if {cond} qui replace `draw' = . if `touse'\n"
            assert original.count(needle)==1;modified=original.replace(needle,injection+needle);f.write_text(modified)
            (work/'instrumentation.patch').write_text(''.join(difflib.unified_diff(original.splitlines(True),modified.splitlines(True))))
            cmd='dea_sw2007 x=y, envars(z) rts(crs) reps(10) bcreps(10) seed(127) gen(a) replacevars'
            body=DATA+'scalar TEST_ATTEMPTS=0\ngen a_eff=991\ngen a_eff_bc=992\n'
            if permanent:body+=unchanged(cmd,498)+'\nassert "`e(cmd)\'"==""\n'
            else:
                body+=cmd+'\nassert e(bc_success)==10 & e(boot_success)==10\n'
                body+=f'assert e({"bc_failed" if stage=="dea" else "boot_failed"})==2\nassert a_eff_bc>0 & a_eff_bc<1\n'
            (work/'test.do').write_text(f'clear all\nset more off\nadopath ++ "{pkg}"\n'+body+'\ndi "TEST_PASS"\n')
            subprocess.run([args.stata,'-b','do','test.do'],cwd=work,timeout=120,capture_output=True)
            log=(work/'test.log').read_text(errors='replace');row=dict(test=name,passed='TEST_PASS' in log.splitlines())
            if not row['passed']:row['tail']=log[-1800:]
            results.append(row);(out/'results.json').write_text(json.dumps(results,indent=2));print(json.dumps(row),flush=True)
    raise SystemExit(not all(r['passed'] for r in results))

if __name__=='__main__':main()
