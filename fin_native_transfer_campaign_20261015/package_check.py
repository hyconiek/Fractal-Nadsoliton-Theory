#!/usr/bin/env python3
import hashlib,pathlib,sys
root=pathlib.Path(__file__).resolve().parent
pairs=[]
for line in (root/'MANIFEST.sha256').read_text().splitlines():
    if not line.strip(): continue
    expected,relative=line.split('  ',1)
    path=root/relative
    assert path.is_file(),relative
    actual=hashlib.sha256(path.read_bytes()).hexdigest()
    assert actual==expected,(relative,expected,actual)
    pairs.append(relative)
expected_files=set(p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and p.name not in ('MANIFEST.sha256','PACKAGE_VERIFICATION.json','package_check.py'))
assert set(pairs)==expected_files,{'missing':sorted(expected_files-set(pairs)),'surplus':sorted(set(pairs)-expected_files)}
print(f'PASS {len(pairs)}/{len(pairs)} source-manifest hashes and complete coverage; package metadata/checker excluded explicitly')
