#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parent
entries=(root/'MANIFEST.sha256').read_text().splitlines()
for entry in entries:
 hash_,rel=entry.split('  ',1)
 path=root/rel
 assert path.is_file(), f'missing {rel}'
 assert hashlib.sha256(path.read_bytes()).hexdigest()==hash_, f'bad sha {rel}'
actual={str(x.relative_to(root)) for x in root.rglob('*') if x.is_file() and x.name not in ('MANIFEST.sha256','PACKAGE_VERIFICATION.json')}
listed={line.split('  ',1)[1] for line in entries}
assert actual==listed,(sorted(actual-listed),sorted(listed-actual))
for fn in ['CLAIM_REGISTER.json','SOURCE_PROVENANCE.json','EXACT_DIAGNOSTICS.json','BA-001/CLASS_AND_FLUX_CONTRACT.json','BA-003/MODEL_FAMILY.json']:
 json.loads((root/fn).read_text())
assert json.loads((root/'EXACT_DIAGNOSTICS.json').read_text())['check_count']==24
print('PACKAGE PASS',len(entries),'SHA entries, complete coverage; exact checks 24')
