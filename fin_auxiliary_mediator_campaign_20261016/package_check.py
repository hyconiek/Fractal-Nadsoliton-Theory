#!/usr/bin/env python3
from pathlib import Path
import hashlib
root=Path(__file__).parent
rows=(root/'MANIFEST.sha256').read_text().splitlines()
for row in rows:
    digest,name=row.split('  ',1)
    p=root/name
    assert p.is_file(), f'missing: {name}'
    assert hashlib.sha256(p.read_bytes()).hexdigest()==digest, f'hash mismatch {name}'
actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and p.name not in ('MANIFEST.sha256','PACKAGE_VERIFICATION.json','__pycache__') and '__pycache__' not in p.parts}
expected={row.split('  ',1)[1] for row in rows}
assert actual==expected, ('uncovered files',actual-expected,expected-actual)
print(f'HM manifest: {len(rows)}/{len(rows)} valid with coverage')
