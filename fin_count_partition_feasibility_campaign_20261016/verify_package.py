#!/usr/bin/env python3
"""Independent on-disk package SHA256 and file coverage validator."""
from pathlib import Path
import hashlib
P=Path(__file__).resolve().parent
items={}
for ln in (P/'MANIFEST.sha256').read_text().splitlines():
 if not ln.strip():continue
 digest,name=ln.split('  ',1)
 assert name not in items and (P/name).is_file(),name
 assert hashlib.sha256((P/name).read_bytes()).hexdigest()==digest,name
 items[name]=digest
actual={str(f.relative_to(P)) for f in P.rglob('*') if f.is_file()}-{'MANIFEST.sha256','PACKAGE_VERIFICATION.json'}
assert actual==set(items),f'missing {actual-set(items)} extra {set(items)-actual}'
print(f'PASS {len(items)}/{len(items)} SHA256 manifest members, coverage complete')
