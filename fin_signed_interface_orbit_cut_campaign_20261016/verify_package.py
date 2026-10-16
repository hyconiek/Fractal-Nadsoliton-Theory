#!/usr/bin/env python3
"""Verify this incremental package's manifest. No scientific status inferred."""
from pathlib import Path
import hashlib
root=Path(__file__).resolve().parent
manifest=root/'MANIFEST.sha256'
found=[]
for row in manifest.read_text().splitlines():
 if not row.strip():continue
 digest,rel=row.split('  ',1)
 p=root/rel
 assert p.is_file(),f'file missing: {rel}'
 actual=hashlib.sha256(p.read_bytes()).hexdigest()
 assert actual==digest,f'hash mismatch: {rel}'
 found.append(rel)
actuals=sorted(str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p.name!='MANIFEST.sha256')
assert sorted(found)==actuals,f'coverage mismatch manifest={len(found)}, actual={len(actuals)}'
print(f'PASS {len(found)}/{len(actuals)} file hashes and complete coverage')
