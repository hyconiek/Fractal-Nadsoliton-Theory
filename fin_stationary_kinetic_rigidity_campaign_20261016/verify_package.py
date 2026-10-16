#!/usr/bin/env python3
from pathlib import Path
from hashlib import sha256
root=Path(__file__).resolve().parent
rows=[line.split('  ',1) for line in (root/'MANIFEST.sha256').read_text().splitlines() if line.strip()]
expected={filename:digest for digest,filename in rows}
actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p.name!='MANIFEST.sha256' and '__pycache__' not in p.parts}
assert actual==set(expected),(sorted(actual-set(expected)),sorted(set(expected)-actual))
for filename,digest in expected.items():
 assert sha256((root/filename).read_bytes()).hexdigest()==digest,filename
print(f'MANIFEST PASS {len(rows)}/{len(rows)}; complete file coverage')
