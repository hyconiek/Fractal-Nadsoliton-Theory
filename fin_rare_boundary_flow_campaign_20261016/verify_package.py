#!/usr/bin/env python3
"""Read-only RB file manifest checker. Does not independently certify mathematical proofs."""
import hashlib,pathlib,json,sys
P=pathlib.Path(__file__).resolve().parent
lines=(P/'MANIFEST.sha256').read_text().splitlines()
listed={}
for line in lines:
 sha,rel=line.split('  ',1)
 p=P/rel
 assert rel not in listed,rel
 assert p.is_file(),rel
 assert hashlib.sha256(p.read_bytes()).hexdigest()==sha,rel
 listed[rel]=sha
actual={p.relative_to(P).as_posix() for p in P.rglob('*') if p.is_file() and p.name not in {'MANIFEST.sha256','PACKAGE_VERIFICATION.json'}}
assert set(listed)==actual,(sorted(actual-set(listed)),sorted(set(listed)-actual))
print(f'PASS {len(listed)}/{len(actual)} manifest entries and full coverage')
