#!/usr/bin/env python3
"""Raw-point failure reasons, withheld exports and restored test-only work limits.

Synthetic integer packets are generated locally. No private data, Git history,
R, or numerical tolerance changes are required. Work-limit injection is confined
to this test's Stata process; the complete helper is reloaded afterwards.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import struct
import subprocess

REPO = Path(__file__).resolve().parents[1]
MISSING = 2. ** 1023


def packet(value):
    sign = (value > 0) - (value < 0)
    value = abs(value)
    digits = []
    while value:
        digits.append(value % 2**26)
        value //= 2**26
    return [sign, *(digits or [0])]


def writebin(path, matrix):
    n, m = len(matrix), len(matrix[0])
    path.write_bytes(struct.pack('<ii', n, m) +
                     struct.pack('<' + 'd' * (n*m), *(v for row in matrix for v in row)))


def readbin(path):
    data = path.read_bytes()
    n, m = struct.unpack('<ii', data[:8])
    values = struct.unpack('<' + 'd' * (n*m), data[8:])
    return [list(values[i*m:(i+1)*m]) for i in range(n)]


def cases():
    rows = []
    def add(name, nums, denominator=1, exponents=None, free=0, reason=0,
            ops=2, limbs=2, limit=10000000, cap=512, injected=False):
        rows.append(dict(name=name, nums=nums, denominator=denominator,
                         exponents=exponents or [0] * (len(nums)-free), free=free,
                         reason=reason, ops=ops, limbs=limbs, limit=limit,
                         cap=cap, injected=injected))
    add('zero', [0], ops=1)
    add('one', [1])
    add('third', [1], 3)
    add('tie_even', [2**53+1], 2**53, limbs=4)
    add('tie_odd', [2**53+3], 2**53, limbs=4)
    add('minsub', [1], 2**1074, limbs=43)
    add('underflow', [1], 2**1075, reason=4, limbs=43)
    add('overflow', [1], exponents=[1023], reason=3, ops=81, limbs=41)
    add('late_overflow', [1, 1], exponents=[0, 1023], reason=3, ops=83, limbs=41)
    add('late_underflow', [1, 1], exponents=[0, -1075], reason=4, ops=87, limbs=43)
    add('free_positive', [1, 3, 1], free=1, ops=6)
    add('free_negative', [1, 1, 3], free=1, ops=6)
    add('free_zero', [1, 1, 1], free=1, ops=5)
    add('free_cancel', [1, 2**500+1, 2**500], free=1, ops=25, limbs=20)
    add('zero_denominator', [1], 0, reason=1, ops=1)
    add('bad_positive_shift', [1], exponents=[.5], reason=14, ops=0, limbs=0)
    add('bad_negative_shift', [1], exponents=[-.5], reason=14, ops=1)
    add('wide_positive_shift', [1], exponents=[13312], reason=12, ops=513, limbs=514)
    add('wide_negative_shift', [1], exponents=[-13312], reason=12, ops=514, limbs=514)
    add('zero_positive_shift', [0], exponents=[13312], ops=1)
    add('work_zero', [1], reason=11, ops=1, limit=0, injected=True)
    add('work_one', [1], reason=11, limit=1, injected=True)
    add('work_two', [1], limit=2, injected=True)
    add('limb_one', [1], reason=12, ops=1, cap=1, injected=True)
    add('combined_limit', [1], reason=12, ops=1, limit=0, cap=1, injected=True)
    add('late_budget', [0, 1], reason=11, ops=3, limit=1, injected=True)
    add('free_budget', [1, 3, 1], free=1, reason=11, ops=5, limit=2, injected=True)
    add('scale_cap', [1], exponents=[26], reason=12, limbs=3, cap=2, injected=True)
    return rows


IO = r'''
mata:
void point_save(string scalar path,real matrix A) {
 real matrix B
 real scalar f
 B=bufio();bufbyteorder(B,2);f=fopen(path,"w")
 fbufput(B,f,"%4b",rows(A));fbufput(B,f,"%4b",cols(A));fbufput(B,f,"%8z",A);fclose(f)
}
real matrix point_load(string scalar path) {
 real matrix B,A
 real scalar f,n,m
 B=bufio();bufbyteorder(B,2);f=fopen(path,"r")
 n=fbufget(B,f,"%4b",1,1);m=fbufget(B,f,"%4b",1,1)
 A=fbufget(B,f,"%8z",n,m);fclose(f);return(A)
}
void point_assert_defaults() {
 struct _C58_Work scalar W
 W=_c84_work();assert(W.limit==10000000 & W.limbs==512 & W.ops==0 & W.error==0)
}
void point_context(string scalar prefix) {
 string rowvector names
 pointer scalar p
 real scalar i
 names=_c57_context_names()
 for(i=1;i<=cols(names);i++) {
  p=findexternal(names[i]);assert(p!=NULL)
  point_save(prefix+strofreal(i)+".bin",*p)
 }
}
void point_case(string scalar name,real scalar injected) {
 struct _C58_Cert scalar fit
 struct _C58_Work scalar W
 real rowvector v,p,x,d
 pointer scalar budget,cap
 fit.numerator=point_load(name+"/N.bin");fit.denominator=point_load(name+"/D.bin")
 v=point_load(name+"/ve.bin");p=point_load(name+"/parameters.bin")
 if(injected) {
  budget=findexternal("_point_test_limit");cap=findexternal("_point_test_limbs")
  *budget=p[2];*cap=p[3]
 }
 W=_c84_work();point_save(name+"/work.bin",(W.limit,W.limbs,W.ops,W.error))
 point_context(name+"/before_");x=_c84_point(v,p[1],fit,d);point_context(name+"/after_")
 point_save(name+"/point.bin",x);point_save(name+"/diag.bin",d)
 point_save(name+"/echo_N.bin",fit.numerator);point_save(name+"/echo_D.bin",fit.denominator)
}
end
'''
OVERRIDE = r'''
capture mata: mata drop _c84_work()
mata:
struct _C58_Work scalar _c84_work() {
 external real scalar _point_test_limit,_point_test_limbs
 struct _C58_Work scalar W
 W.ops=0;W.limit=_point_test_limit;W.limbs=_point_test_limbs;W.maxlimbs=0;W.error=0
 return(W)
}
end
'''


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stata', required=True)
    ap.add_argument('--output', required=True, type=Path)
    args = ap.parse_args()
    out = args.output.resolve()
    if out == REPO or REPO in out.parents:
        ap.error('Choose an output directory outside the package.')
    out.mkdir(parents=True, exist_ok=False)
    spec = cases()
    for row in spec:
        d = out / row['name']
        d.mkdir()
        nums = [packet(n) for n in row['nums']]
        width = max(map(len, nums))
        writebin(d/'N.bin', [n+[0]*(width-len(n)) for n in nums])
        writebin(d/'D.bin', [packet(row['denominator'])])
        writebin(d/'ve.bin', [row['exponents']])
        writebin(d/'parameters.bin', [[row['free'], row['limit'], row['cap']]])
    code = 'version 16\nclear all\nset more off\nset processors 1\nset type double\n'
    for name in ['PLUS', 'PERSONAL', 'SITE', 'OLDPLACE']:
        (out/name).mkdir()
        code += f'sysdir set {name} "{out/name}"\n'
    code += f'adopath ++ "{REPO/"ado"}"\n'
    warm = 'clear\ninput id x y\n1 1 1\n2 2 1\nend\nquietly dea_crosseff x=y,dmu(id) eng(simplex) matrix\n'
    code += warm + 'mata: assert(_c84_ready()==8405)\n' + IO
    code += 'mata: point_save("warm_before.bin",st_matrix("r(ce_matrix)"))\n'
    for row in spec:
        if not row['injected']:
            code += f'mata: point_case("{row["name"]}",0)\n'
    code += 'mata: _point_test_limit=10000000;_point_test_limbs=512\n' + OVERRIDE
    for row in spec:
        if row['injected']:
            code += f'mata: point_case("{row["name"]}",1)\n'
    code += f'run "{REPO/"ado/_dea_crosseff_raw.ado"}"\n'
    code += 'mata: point_assert_defaults()\n'
    code += warm + 'mata: point_save("warm_after.bin",st_matrix("r(ce_matrix)"))\n'
    code += 'display "CROSSEFF_POINT_PASS"\nexit,clear\n'
    (out/'test.do').write_text(code)
    proc = subprocess.run([args.stata, '-b', 'do', 'test.do'], cwd=out,
                          capture_output=True, timeout=300)
    log = (out/'test.log').read_text(errors='replace')
    assert proc.returncode == 0 and 'CROSSEFF_POINT_PASS' in log.splitlines(), log[-3000:]
    verified = []
    for row in spec:
        d = out/row['name']
        diag = readbin(d/'diag.bin')[0]
        assert diag == [int(row['reason']==0), row['reason'], row['ops'], row['limbs']], row['name']
        got = readbin(d/'point.bin')[0]
        dim = len(row['exponents'])
        if row['reason']:
            assert got == [MISSING]*dim, row['name']
        else:
            want = []
            for i, e in enumerate(row['exponents']):
                n = row['nums'][i] - (row['nums'][dim] if row['free'] and i==dim-1 else 0)
                q = Fraction(n, row['denominator']) * Fraction(2)**int(e)
                want.append(float(q))
            assert struct.pack('<'+'d'*dim, *got)==struct.pack('<'+'d'*dim, *want), row['name']
        assert readbin(d/'work.bin') == [[row['limit'], row['cap'], 0, 0]], row['name']
        for n in ['N', 'D']:
            assert (d/(n+'.bin')).read_bytes() == (d/('echo_'+n+'.bin')).read_bytes()
        for i in range(1, 20):
            assert (d/f'before_{i}.bin').read_bytes() == (d/f'after_{i}.bin').read_bytes()
        verified.append(dict(case=row['name'], passed=True, reason=row['reason'],
                             test_only_limit_injected=row['injected']))
    assert (out/'warm_before.bin').read_bytes() == (out/'warm_after.bin').read_bytes()
    summary = dict(passed=True, cases=len(spec), default_work_cases=sum(not r['injected'] for r in spec),
                   injected_work_cases=sum(r['injected'] for r in spec),
                   withheld_cases=sum(bool(r['reason']) for r in spec),
                   state_matrix_byte_comparisons=19*len(spec), exact_fraction_points=True,
                   defaults_restored=True, command_before_after_identical=True,
                   scope='synthetic conversion/helper regression; invalid and extended inputs included; not LP infeasibility proofs')
    (out/'results.json').write_text(json.dumps(verified, indent=2)+'\n')
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
