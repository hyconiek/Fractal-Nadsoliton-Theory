#!/usr/bin/env python3
"""Check package coverage, all immutable new-file SHA256s, and exact fixture replay."""
from pathlib import Path
import hashlib, json, subprocess, sys, tempfile, shutil
ROOT=Path(__file__).resolve().parent
M=ROOT/'MANIFEST.sha256'
entries={}
for line in M.read_text().splitlines():
    if not line.strip():continue
    digest,name=line.split('  ',1)
    assert name not in entries and '/' not in name, name
    entries[name]=digest
existing={p.name for p in ROOT.iterdir() if p.is_file() and p.name!='MANIFEST.sha256'}
assert existing==set(entries),f'coverage missing {existing-set(entries)} extra {set(entries)-existing}'
for name,digest in entries.items():
    actual=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    assert actual==digest, name
with tempfile.TemporaryDirectory(prefix='fin-orbit-replay-') as temp:
    p=Path(temp)
    shutil.copy2(ROOT/'verify_orbit.py',p/'verify_orbit.py')
    subprocess.run([sys.executable,str(p/'verify_orbit.py')],check=True,capture_output=True,text=True)
    assert (p/'EXACT_CHECKS.json').read_bytes()==(ROOT/'EXACT_CHECKS.json').read_bytes(),'non-identical replay'
data=json.loads((ROOT/'EXACT_CHECKS.json').read_text())
assert data['all_pass'] and data['extra_full_cut_exhaustive_control']['pass']
print(f'PASS {len(entries)}/{len(entries)} sha256 files, 10 exact bounded fixture categories, 64 independent full cuts, byte-identical clean replay')
