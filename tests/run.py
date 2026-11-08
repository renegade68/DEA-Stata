#!/usr/bin/env python3
"""Run public DEA regression suites sequentially, without editing their source."""
import argparse
import ast
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from _support import executable, package_files, public_files

REPO = Path(__file__).resolve().parents[1]


def suites(repo=REPO):
    result = {}
    for name in public_files(repo):
        p = Path(name)
        if len(p.parts) != 2 or p.parts[0] != 'tests' or p.suffix != '.py':
            continue
        if not p.stem.startswith(('test_', 'verify_')) or p.stem == 'test_portability':
            continue
        tree = ast.parse((Path(repo) / p).read_text(encoding='utf-8'))
        flags, required, dependencies = set(), set(), set()
        for n in ast.walk(tree):
            if isinstance(n, ast.Import):
                dependencies.update(x.name.split('.')[0] for x in n.names)
            elif isinstance(n, ast.ImportFrom) and n.module:
                dependencies.add(n.module.split('.')[0])
            if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) and n.func.attr == 'add_argument':
                options = {v.value for v in n.args if isinstance(v, ast.Constant) and isinstance(v.value, str)}
                flags.update(options)
                if any(k.arg == 'required' and isinstance(k.value, ast.Constant) and k.value.value is True for k in n.keywords):
                    required.update(options)
        if '--stata' not in flags or '--output' not in flags:
            continue
        if p.stem in ('test_output_protection', 'test_sample_selection'):
            required.add('--baseline-ado')
        if p.stem == 'test_sample_selection':
            # Its historical crossgame control uses the shared independent LP verifier.
            dependencies.update({'numpy', 'scipy'})
        key = p.stem.removeprefix('test_')
        result[key] = dict(file=name, flags=flags, required=required - {'--stata', '--output'},
                           dependencies=dependencies & {'numpy', 'scipy'})
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--list', action='store_true', help='List suites and historical-input requirements')
    ap.add_argument('--suite', action='append', help='Suite name; repeat to run sequentially. Default: smoke')
    ap.add_argument('--stata', help='Overrides DEA_STATA and automatic discovery')
    ap.add_argument('--rscript', help='Opt in to R checks; overrides DEA_RSCRIPT')
    ap.add_argument('--with-r', action='store_true', help='Opt in to R checks using discovery or DEA_RSCRIPT')
    ap.add_argument('--output', type=Path, help='New output directory; default: a unique OS temporary directory')
    ap.add_argument('--baseline-ado', type=Path, help='Reviewed historical ado directory, never the current package')
    ap.add_argument('--only', help='Existing case-name regex, for suites that support --only')
    a = ap.parse_args()
    registry = suites()
    if a.list:
        print('smoke: selected help examples plus fixed raw-recovery regression')
        for name, info in sorted(registry.items()):
            print(name + (': requires ' + ', '.join(sorted(info['required'])) if info['required'] else ''))
        return 0
    jobs = []
    for name in a.suite or ['smoke']:
        if name == 'smoke':
            jobs.extend([('help_syntax', ['--only', '^help_(dea|dea_lp|dea_crosseff)$']), ('crosseff_raw', [])])
        elif name in registry:
            jobs.append((name, []))
        else:
            ap.error('Unknown suite. Use --list.')
    if len({name for name, _ in jobs}) != len(jobs):
        ap.error('Each suite may be selected only once.')
    try:
        package_files(REPO)
        stata = executable('stata', a.stata)
        rscript = executable('rscript', a.rscript) if a.with_r or a.rscript else None
    except ValueError as e:
        ap.error(str(e))
    baseline = a.baseline_ado.expanduser().resolve() if a.baseline_ado else None
    if baseline and (baseline == (REPO / 'ado').resolve() or not (baseline / 'dea.ado').is_file()):
        ap.error('--baseline-ado must name a separate reviewed historical ado directory.')
    for name, extra in jobs:
        info = registry[name]
        missing = info['required'] - ({'--baseline-ado'} if baseline else set())
        if missing:
            ap.error(name + ' needs ' + ', '.join(sorted(missing)) + '. It is not a ZIP-only current-package test.')
        for dep in info['dependencies']:
            if importlib.util.find_spec(dep) is None:
                ap.error('Missing ' + dep + '. Install tests/requirements.txt into your Python environment.')
        if a.only:
            if '--only' not in info['flags'] or '--only' in extra:
                ap.error('--only is not supported for this selection.')
            extra += ['--only', a.only]
        if baseline and '--baseline-ado' in info['flags']:
            extra += ['--baseline-ado', str(baseline)]
        if rscript and '--rscript' in info['flags']:
            extra += ['--rscript', rscript]
    if rscript and not any('--rscript' in registry[name]['flags'] for name, _ in jobs):
        ap.error('The selected suites do not have R comparisons.')
    if a.output:
        out = a.output.expanduser().resolve()
        if out == REPO or REPO in out.parents:
            ap.error('Choose an output directory outside the package so results cannot enter the public inventory.')
        out.mkdir(parents=True, exist_ok=False)
    else:
        out = Path(tempfile.mkdtemp(prefix='dea-tests-')).resolve()
    print('Local results: ' + str(out), flush=True)
    summary = dict(schema=1, scope='selected regression suites', results=[],
                   raw_outputs_private=True, r_requested=bool(rscript),
                   requested_suites=[name for name, _ in jobs])
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    for name, extra in jobs:
        cmd = [sys.executable, '-B', str(REPO / registry[name]['file']), '--stata', stata,
               '--output', str(out / name), *extra]
        print('Running ' + name, flush=True)
        with (out / (name + '.console.txt')).open('w', encoding='utf-8') as f:
            p = subprocess.run(cmd, cwd=REPO, env=env, stdout=f, stderr=subprocess.STDOUT)
        summary['results'].append(dict(suite=name, passed=p.returncode == 0, returncode=p.returncode,
                                       case_filter_used='--only' in extra))
        summary['passed'] = all(x['passed'] for x in summary['results']) and len(summary['results']) == len(jobs)
        # This report deliberately excludes commands, absolute paths and log excerpts.
        (out / 'shareable-summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
        print(name + (': PASS' if p.returncode == 0 else ': FAIL (inspect the local console file)'), flush=True)
        if p.returncode:
            return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
