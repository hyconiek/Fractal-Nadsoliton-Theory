#!/usr/bin/env python3
"""Check local manifest; this does not validate the mathematical results."""
from pathlib import Path
import hashlib
p=Path(__file__).resolve().parent
lines=(p/'MANIFEST.sha256').read_text().splitlines()
seen=set()
for line in lines:
 h,f=line.split('  ',1)
 assert f not in seen, f
 seen.add(f)
 data=(p/f).read_bytes()
 assert hashlib.sha256(data).hexdigest()==h, f
actual={q.relative_to(p).as_posix() for q in p.rglob('*') if q.is_file() and q.name!='MANIFEST.sha256'}
assert seen==actual,(sorted(actual-seen),sorted(seen-actual))
print(f'MANIFEST PASS {len(lines)}/{len(lines)} files with complete coverage')
