"""Verify current RI handoff contents against its separately authored manifest."""
from pathlib import Path
from hashlib import sha256
root=Path(__file__).resolve().parent
manifest=(root/'MANIFEST.sha256').read_text().splitlines()
seen=set()
for line in manifest:
    h,rel=line.split('  ',1)
    assert rel not in seen
    seen.add(rel)
    f=root/rel
    assert f.is_file(), f'missing {rel}'
    assert sha256(f.read_bytes()).hexdigest()==h, f'hash error {rel}'
existing={str(f.relative_to(root)) for f in root.rglob('*') if f.is_file() and '__pycache__' not in f.parts and f.name!='MANIFEST.sha256'}
assert existing==seen, f'extra {existing-seen}, missing manifest {seen-existing}'
print(f'MANIFEST PASS {len(seen)}/{len(seen)}')
