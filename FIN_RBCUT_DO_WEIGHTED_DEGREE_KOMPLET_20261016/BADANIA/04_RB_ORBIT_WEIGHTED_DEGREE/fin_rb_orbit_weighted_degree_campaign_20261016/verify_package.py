#!/usr/bin/env python3
"""Verify this output manifest. Integrity != analytical proof."""
from pathlib import Path
import hashlib
D=Path(__file__).resolve().parent
lines=(D/'MANIFEST.sha256').read_text().splitlines()
assert lines and len(set(lines))==len(lines)
covered=set()
for line in lines:
    h,rel=line.split('  ',1)
    p=D/rel
    assert p.is_file() and '..' not in Path(rel).parts,(rel,'missing')
    assert hashlib.sha256(p.read_bytes()).hexdigest()==h,rel
    covered.add(rel)
actual={p.name for p in D.iterdir() if p.is_file()}-{'MANIFEST.sha256'}
assert covered==actual,(actual-covered,covered-actual)
print(f'PASS manifest {len(lines)}/{len(actual)} with complete coverage')
