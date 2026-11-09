#!/usr/bin/env python3
"""Standard-library checks of package discovery, privacy gates and runner setup."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stdout, redirect_stderr
import io
from unittest.mock import patch
import _support as support
import check_publication as privacy
import run as runner


class Portability(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='dea-portability-')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name) / 'package space 한글'
        self.root.mkdir()
        for name, value in [('ado/dea.ado', '* test source\n'), ('tests/sample.txt', 'safe\n')]:
            p = self.root / name
            p.parent.mkdir(exist_ok=True)
            p.write_text(value, encoding='utf-8')
        self.inventory = self.root / 'tests/public_files.json'
        self.write_inventory(['ado/dea.ado', 'tests/public_files.json', 'tests/sample.txt'])

    def write_inventory(self, names):
        self.inventory.write_text(json.dumps(sorted(names)), encoding='utf-8')

    def test_zip_inventory_and_unlisted_scratch(self):
        (self.root / 'ado/scratch.ado').write_text('* not public')
        with patch('subprocess.check_output', side_effect=AssertionError('Git must not be called')):
            self.assertEqual(support.package_files(self.root), ['ado/dea.ado'])

    def test_missing_inventory_member(self):
        (self.root / 'ado/dea.ado').unlink()
        with self.assertRaises(ValueError):
            support.package_files(self.root)

    def test_inventory_traversal(self):
        self.write_inventory(['../external'])
        with self.assertRaises(ValueError):
            support.public_files(self.root)

    def test_inventory_symlink(self):
        p = self.root / 'ado/dea.ado'
        p.unlink()
        p.symlink_to(self.root / 'tests/sample.txt')
        with self.assertRaises(ValueError):
            support.public_files(self.root)

    def test_inventory_duplicate(self):
        self.write_inventory(['ado/dea.ado', 'ado/dea.ado'])
        with self.assertRaises(ValueError):
            support.public_files(self.root)

    def test_historical_zip_rejected(self):
        with self.assertRaisesRegex(ValueError, 'historical comparison'):
            support.require_git_history(self.root, 'HEAD')

    def make_executable(self, name):
        p = Path(self.tmp.name) / name
        p.write_text('#!/bin/sh\nexit 0\n')
        p.chmod(0o755)
        return p

    def test_executable_precedence_and_spaces(self):
        explicit = self.make_executable('explicit 한글')
        env_path = self.make_executable('environment')
        env = {'DEA_STATA': str(env_path), 'PATH': ''}
        self.assertEqual(support.executable('stata', str(explicit), env, 'linux'), str(explicit.resolve()))
        self.assertEqual(support.executable('stata', environ=env, platform='linux'), str(env_path.resolve()))

    def test_executable_path_discovery(self):
        p = self.make_executable('stata-mp')
        self.assertEqual(support.executable('stata', environ={'PATH': str(p.parent)}, platform='linux'), str(p.resolve()))

    def test_invalid_explicit_never_falls_back(self):
        p = self.make_executable('stata-mp')
        with self.assertRaises(ValueError):
            support.executable('stata', 'absent-command', {'PATH': str(p.parent)}, 'linux')

    def test_non_executable_rejected(self):
        p = self.make_executable('not-executable')
        p.chmod(0o644)
        with self.assertRaises(ValueError):
            support.executable('stata', str(p), {}, 'linux')

    def test_privacy_home_paths_and_tokens(self):
        values = ['/' + 'Users' + '/example/private', '/' + 'home' + '/example/file',
                  'C:' + '\\' + 'Users' + '\\' + 'example' + '\\' + 'file',
                  'gh' + 'p_' + 'A' * 36,
                  '-----BEGIN ' + 'PRIVATE KEY-----',
                  'https:' + '//user:password' + chr(64) + 'example.test/file']
        for value in values:
            with self.subTest(value_type=values.index(value)):
                report, _ = privacy.scan_data('tests/sample.txt', value.encode())
                self.assertTrue(report)
                self.assertNotIn(value, json.dumps(report))

    def test_privacy_contact_and_matrix_product(self):
        address = 'sample' + chr(64) + 'example.test'
        self.assertTrue(privacy.scan_data('tests/sample.txt', address.encode())[0])
        self.assertFalse(privacy.scan_data('tests/sample.py', b'value = q@a - g.fun\n')[0])

    def test_private_terms(self):
        self.assertTrue(privacy.scan_data('tests/sample.txt', b'PrivateProject', ['privateproject'])[0])
        report, _ = privacy.scan_data('tests/PrivateProject.txt', b'safe', ['privateproject'])
        self.assertNotIn('PrivateProject', json.dumps(report))

    def test_temporary_home_metadata(self):
        path = '/' + 'private/var/folders/' + 'ab/' + 'x' * 24 + '/T/example'
        self.assertTrue(privacy.scan_data('tests/sample.txt', path.encode())[0])

    def run_runner(self, args, fake_run):
        with patch.object(sys, 'argv', ['run.py', *args]), patch.object(runner, 'executable', return_value='/local/fake-stata'), \
                patch.object(runner.subprocess, 'run', side_effect=fake_run), redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            return runner.main()

    def test_runner_sequential_failure_and_private_summary(self):
        out = Path(self.tmp.name) / 'local results'
        seen = []
        def fake_run(cmd, **kwargs):
            seen.append(cmd)
            return subprocess.CompletedProcess(cmd, 1)
        code = self.run_runner(['--suite', 'saveref', '--suite', 'help_syntax', '--output', str(out)], fake_run)
        self.assertEqual(code, 1)
        self.assertEqual(len(seen), 1)
        self.assertIn('--stata', seen[0])
        report = (out / 'shareable-summary.json').read_text()
        self.assertNotIn(str(out), report)
        self.assertNotIn('/local/', report)
        self.assertFalse(json.loads(report)['passed'])

    def test_runner_success_and_environment_independence(self):
        out = Path(self.tmp.name) / 'success results'
        def fake_run(cmd, **kwargs):
            self.assertEqual(kwargs['env']['PYTHONDONTWRITEBYTECODE'], '1')
            return subprocess.CompletedProcess(cmd, 0)
        self.assertEqual(self.run_runner(['--suite', 'saveref', '--output', str(out)], fake_run), 0)
        self.assertTrue(json.loads((out / 'shareable-summary.json').read_text())['passed'])

    def test_runner_missing_historical_input_before_execution(self):
        out = Path(self.tmp.name) / 'unused'
        def fake_run(*args, **kwargs):
            self.fail('Should not launch any suite')
        with self.assertRaises(SystemExit) as error:
            self.run_runner(['--suite', 'projection_recovery', '--output', str(out)], fake_run)
        self.assertEqual(error.exception.code, 2)
        self.assertFalse(out.exists())

    def test_runner_output_inside_package_rejected(self):
        def fake_run(*args, **kwargs):
            self.fail('Should not launch any suite')
        with self.assertRaises(SystemExit):
            self.run_runner(['--suite', 'saveref', '--output', str(runner.REPO / 'private-results')], fake_run)

    def test_binary_rejected(self):
        self.assertEqual(privacy.scan_data('tests/sample.dat', b'\xff')[0][0]['kind'], 'unreviewed_binary')

    def test_registry_marks_historical_inputs(self):
        registry = runner.suites()
        self.assertEqual(len(registry), 90)
        self.assertFalse(registry["twostage_quality"]["required"])
        self.assertFalse(registry["games_quality"]["required"])
        for name in ("crosseff_costunits", "crosseff_costunits_flow", "crosseff_costunits_link", "crossgame_recovery"):
            self.assertEqual(registry[name]["required"], set())
        self.assertIn('crosseff_decoder', registry)
        self.assertIn('crosseff_dyadic', registry)
        self.assertFalse(registry['crosseff_dyadic']['required'])
        for name in ('sample_selection', 'output_protection', 'projection_recovery'):
            self.assertIn('--baseline-ado', registry[name]['required'])
        self.assertFalse(registry['crosseff_raw']['required'])

    @unittest.skipUnless(shutil.which('git'), 'Git is optional for ZIP use; needed only for index-audit tests')
    def test_index_bytes_and_email_metadata(self):
        def git(*args):
            return subprocess.check_output(['git', '--no-optional-locks', '-C', str(self.root), *args], stderr=subprocess.PIPE)
        git('-c', 'init.defaultBranch=main', 'init', '-q')
        address = 'sample' + chr(64) + 'example.test'
        (self.root / 'tests/sample.txt').write_text(address)
        git('add', '--', 'ado/dea.ado', 'tests/sample.txt', 'tests/public_files.json')
        (self.root / 'tests/sample.txt').write_text('safe\n')
        self.assertTrue(privacy.audit(self.root)['passed'])
        self.assertFalse(privacy.audit(self.root, staged=True)['passed'])
        git('add', '--', 'tests/sample.txt')
        env = {k: 'Test Author' for k in ('GIT_AUTHOR_NAME', 'GIT_COMMITTER_NAME')}
        env.update({k: 'test' + chr(64) + 'users.noreply.github.com' for k in ('GIT_AUTHOR_EMAIL', 'GIT_COMMITTER_EMAIL')})
        with patch.dict(os.environ, env):
            self.assertTrue(privacy.audit(self.root, staged=True, check_identity=True)['passed'])
            with patch.dict(os.environ, {'GIT_COMMITTER_EMAIL': address}):
                self.assertFalse(privacy.audit(self.root, staged=True, check_identity=True)['passed'])
        (self.root / 'private.txt').write_text('not approved')
        git('add', '--', 'private.txt')
        self.assertFalse(privacy.audit(self.root, staged=True)['passed'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
