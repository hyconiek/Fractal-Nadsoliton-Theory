#!/usr/bin/env python3
"""Check the new handoff manifest and deterministic bounded science fixtures."""
from pathlib import Path
import hashlib, subprocess, sys
base=Path(__file__).resolve().parent
manifest=base/'MANIFEST.sha256'
rows=manifest.read_text().strip().splitlines()
expected={}
for row in rows:
    digest,path=row.split('  ',1)
    assert path not in expected, path
    expected[path]=digest
files={p.relative_to(base).as_posix() for p in base.rglob('*') if p.is_file() and p.name!='MANIFEST.sha256'}
assert files==set(expected),f'manifest coverage missing={files-set(expected)} stale={set(expected)-files}'
for path,digest in expected.items():
    assert hashlib.sha256((base/path).read_bytes()).hexdigest()==digest,path
result=subprocess.run([sys.executable,str(base/'verify_cut_shell.py')],cwd=base,capture_output=True,text=True,check=True)
print('PACKAGE PASS',len(rows),'hash entries and deterministic checks')
print(result.stdout.strip())
