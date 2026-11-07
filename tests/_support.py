"""Package discovery and local executable configuration for public tests."""
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys


def public_files(repo):
    repo = Path(repo).resolve()
    names = json.loads((repo / 'tests/public_files.json').read_text(encoding='utf-8'))
    if not isinstance(names, list) or not names or names != sorted(set(names)):
        raise ValueError('Invalid public file inventory')
    for name in names:
        if not isinstance(name, str):
            raise ValueError('Invalid public file name')
        p = PurePosixPath(name)
        if (p.is_absolute() or '..' in p.parts or '\\' in name or ':' in name
                or p.as_posix() != name or '.git' in p.parts):
            raise ValueError('Unsafe public file name')
        full = repo / name
        if full.is_symlink() or any(q.is_symlink() for q in full.parents if q != repo and repo in q.parents):
            raise ValueError('Public inventory must not contain symlinks')
        if not full.is_file() or repo not in full.resolve().parents:
            raise ValueError('Missing or external public file: ' + name)
    return names


def package_files(repo):
    """The same reviewed runtime files in a clone and a source ZIP, without Git."""
    names = [n for n in public_files(repo) if n.startswith('ado/')]
    if not names or any(len(PurePosixPath(n).parts) != 2 or
                        PurePosixPath(n).suffix not in ('.ado', '.sthlp', '.do') for n in names):
        raise ValueError('Invalid runtime file inventory')
    return names


def require_git_history(repo, ref):
    """Historical comparisons never silently become current-versus-current tests."""
    repo = Path(repo).resolve()
    if not (repo / '.git').exists() or not shutil.which('git'):
        raise ValueError('This historical comparison needs a Git clone, or an explicit '
                         '--baseline-ado where supported. Current-package tests work from ZIP.')
    proc = subprocess.run(['git', '--no-optional-locks', '-C', str(repo),
                           'rev-parse', '--verify', '--end-of-options', str(ref) + '^{commit}'],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode:
        raise ValueError('Historical commit is unavailable; obtain the required Git history '
                         'or supply --baseline-ado where supported.')


def executable(kind, explicit=None, environ=None, platform=None):
    """Explicit option > environment > PATH > standard macOS application location."""
    env = os.environ if environ is None else environ
    platform = sys.platform if platform is None else platform
    variable = {'stata': 'DEA_STATA', 'rscript': 'DEA_RSCRIPT'}[kind]
    configured = explicit or env.get(variable)
    if configured:
        candidates = [configured]
    else:
        candidates = ['stata-mp', 'stata-se', 'stata'] if kind == 'stata' else ['Rscript']
        if platform == 'darwin':
            if kind == 'stata':
                candidates += ['/Applications/Stata/' + app + '.app/Contents/MacOS/' + binary
                               for app, binary in [('StataMP', 'stata-mp'), ('StataSE', 'stata-se'),
                                                   ('Stata', 'stata')]]
            else:
                candidates += ['/usr/local/bin/Rscript', '/opt/homebrew/bin/Rscript']
    for candidate in candidates:
        p = Path(candidate).expanduser()
        located = str(p.resolve()) if p.is_file() else shutil.which(str(candidate), path=env.get('PATH', ''))
        if located and Path(located).is_file() and os.access(located, os.X_OK):
            return str(Path(located).resolve())
    raise ValueError('Executable not found or not executable. Set ' + variable +
                     ' or pass --' + kind + '. No program was installed or downloaded.')
