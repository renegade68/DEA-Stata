#!/usr/bin/env python3
"""Audit the public inventory or staged Git blobs without printing matched secrets."""
import argparse
import getpass
import hashlib
import io
import json
import os
from pathlib import Path
import re
import subprocess
import tokenize
from _support import public_files

REPO = Path(__file__).resolve().parents[1]
# Existing public author contact, restricted to attributed ado/help files.
# This permits attribution only, not private contact details in test data.
PUBLIC_CONTACT_SHA256 = {'392fee16f489779ba61d300984fedd3a21385ddbb455a66df71513da468c2098'}
PATTERNS = {
    'home_directory': re.compile(r'''(?:/(?:Users|home)/|[A-Za-z]:[\\/]+Users[\\/]+)[^/\\\s"'<>]+'''),
    'user_temporary_directory': re.compile(r'''/private/var/folders/[A-Za-z0-9]{2}/[A-Za-z0-9_-]{8,}'''),
    'private_key': re.compile(r'-----BEGIN (?:RSA |EC |OPENSSH |DSA )?PRIVATE KEY-----'),
    'access_token': re.compile(r'\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,}|AKIA[A-Z0-9]{16})\b'),
    'url_credentials': re.compile(r'''https?://[^/\s:@]+:[^/\s@]+@'''),
}
EMAIL = re.compile(r'(?<![\w.])[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}(?![\w.])')


def scan_data(name, data, private_terms=()):
    findings, contacts = [], 0
    try:
        text = data.decode('utf-8')
    except UnicodeDecodeError:
        return [dict(file=name, line=0, kind='unreviewed_binary')], 0
    # Preserve line numbers; matrix multiplication in Python is not an email.
    email_lines = text.splitlines()
    if name.endswith('.py'):
        email_lines = [''] * len(email_lines)
        try:
            for t in tokenize.generate_tokens(io.StringIO(text).readline):
                if t.type in (tokenize.STRING, tokenize.COMMENT):
                    for i, part in enumerate(t.string.splitlines()):
                        email_lines[t.start[0] - 1 + i] += part
        except (tokenize.TokenError, IndentationError):
            findings.append(dict(file=name, line=0, kind='invalid_python_tokenization'))
    for i, line in enumerate([name] + text.splitlines()):
        for kind, pattern in PATTERNS.items():
            if pattern.search(line):
                findings.append(dict(file=name, line=i, kind=kind))
        if any(term.casefold() in line.casefold() for term in private_terms if term):
            findings.append(dict(file=name, line=i, kind='private_term'))
    for i, line in enumerate(email_lines, 1):
        for match in EMAIL.finditer(line):
            digest = hashlib.sha256(match.group().lower().encode()).hexdigest()
            if name.startswith('ado/') and digest in PUBLIC_CONTACT_SHA256:
                contacts += 1
            else:
                findings.append(dict(file=name, line=i, kind='contact_address'))
    if any(pattern.search(name) for pattern in PATTERNS.values()) or any(
            term.casefold() in name.casefold() for term in private_terms if term):
        for row in findings:
            row['file'] = '(redacted relative filename)'
    return findings, contacts


def git(repo, *args):
    return subprocess.check_output(['git', '--no-optional-locks', '-C', str(repo), *args])


def audit(repo, staged=False, private_terms=(), check_identity=False):
    repo = Path(repo).resolve()
    names = public_files(repo)
    allowed = set(names)
    findings, contacts = [], 0
    if staged:
        if not (repo / '.git').exists():
            raise ValueError('--staged requires this package to be a Git working tree.')
        # Read the index, including its inventory, rather than trusting working-tree bytes.
        names = [p.decode('utf-8') for p in git(repo, 'ls-files', '-z').split(b'\0') if p]
        try:
            index_inventory = json.loads(git(repo, 'show', ':tests/public_files.json'))
        except subprocess.CalledProcessError as e:
            raise ValueError('Stage the reviewed public inventory before the final audit.') from e
        if set(names) != set(index_inventory) or set(index_inventory) != allowed:
            findings.append(dict(file='tests/public_files.json', line=0, kind='index_inventory_mismatch'))
        allowed = set(index_inventory)
        for entry in git(repo, 'ls-files', '--stage', '-z').split(b'\0'):
            if entry and entry.split(b' ', 1)[0] not in (b'100644', b'100755'):
                findings.append(dict(file='(index entry)', line=0, kind='nonregular_index_file'))
    for name in names:
        if name not in allowed:
            findings.append(dict(file=name, line=0, kind='unreviewed_file'))
            continue
        data = git(repo, 'show', ':' + name) if staged else (repo / name).read_bytes()
        rows, count = scan_data(name, data, private_terms)
        findings.extend(rows)
        contacts += count
    identity = 'not_checked'
    if check_identity:
        for which in ('GIT_AUTHOR_IDENT', 'GIT_COMMITTER_IDENT'):
            value = git(repo, 'var', which).decode('utf-8')
            match = re.search(r'<([^<>]+)>', value)
            if not match or not match[1].lower().endswith('@users.noreply.github.com'):
                findings.append(dict(file='(commit metadata)', line=0, kind=which.lower() + '_not_noreply'))
        identity = 'checked_without_disclosing_identity'
    return dict(schema=1, passed=not findings, scope='index' if staged else 'public_inventory',
                files_checked=len(names), findings=findings, existing_public_contact_occurrences=contacts,
                commit_identity=identity,
                limitations='Pattern-based scan; does not prove absence of all personal data or inspect Git history.')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--staged', action='store_true', help='Audit every index blob and the exact index inventory')
    ap.add_argument('--check-identity', action='store_true', help='Require noreply author and committer emails')
    ap.add_argument('--private-term', action='append', default=[], help='Additional local-only text to exclude; repeat as needed')
    a = ap.parse_args()
    terms = list(a.private_term)
    for term in (getpass.getuser(), Path.home().name):
        if len(term) >= 4 and term not in ('root', 'runner', 'user', 'codex'):
            terms.append(term)
    try:
        report = audit(REPO, a.staged, terms, a.check_identity)
    except (ValueError, OSError, subprocess.CalledProcessError) as e:
        # The error may contain local paths; print only its class on a sharing surface.
        print(json.dumps(dict(passed=False, kind=type(e).__name__, message='Audit could not complete; check the local inventory and Git setup.')))
        return 2
    print(json.dumps(report, indent=2))
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
