#!/usr/bin/env python3
"""Validate the scoped NOC manifest; package verification JSON is excluded."""
from pathlib import Path
import hashlib,json,zipfile
R=Path(__file__).resolve().parent
rows=[]
for line in (R/'MANIFEST.sha256').read_text().splitlines():
 digest,rel=line.split('  ',1)
 p=R/rel
 assert p.is_file(),f'Missing {rel}'
 actual=hashlib.sha256(p.read_bytes()).hexdigest()
 assert actual==digest,f'Hash mismatch {rel}'
 rows.append(rel)
assert len(set(rows))==len(rows)
actual_files={x.relative_to(R).as_posix() for x in R.rglob('*') if x.is_file() and not x.name.endswith('.pyc') and '__pycache__' not in x.parts}-{'MANIFEST.sha256','PACKAGE_VERIFICATION.json'}
assert actual_files==set(rows),(sorted(actual_files-set(rows)),sorted(set(rows)-actual_files))
archive=R.parent/'FIN_NOC_001_NOC_005_ONLY_HANDOFF_20261015.zip'
if archive.is_file():
 assert zipfile.ZipFile(archive).testzip() is None
 zip_test='PASS'
else:
 zip_test='NOT_PRESENT_IN_EXTRACTED_DIRECTORY'
out={'status':'PASS','manifest_entries':len(rows),'zip_test':zip_test,'zip_bytes':archive.stat().st_size if archive.is_file() else None,'scope':'NOC-only; manifest covers files excluding itself and package checker output; historical sources not repackaged'}
(R/'PACKAGE_VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n')
print('PASS',len(rows),'manifest entries, ZIP integrity:',zip_test)
