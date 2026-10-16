#!/usr/bin/env python3
from pathlib import Path
from hashlib import sha256
root=Path(__file__).resolve().parent
manifest=root/'MANIFEST.sha256'; entries={}
for line in manifest.read_text().splitlines():
    digest,rel=line.split('  ',1)
    assert rel not in entries and '..' not in Path(rel).parts
    entries[rel]=digest
actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p.name!='MANIFEST.sha256'}
assert actual==set(entries),{'missing':sorted(actual-set(entries)),'extra':sorted(set(entries)-actual)}
for rel,digest in entries.items():
    assert sha256((root/rel).read_bytes()).hexdigest()==digest,rel
for line in (root/'INPUTS.sha256').read_text().splitlines():
    digest,path=line.split('  ',1)
    p=Path(path)
    assert p.is_file() and sha256(p.read_bytes()).hexdigest()==digest,path
print(f'PASS package {len(entries)}/{len(actual)}; inputs checked')
