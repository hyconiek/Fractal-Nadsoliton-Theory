#!/usr/bin/env python3
"""Validate scope/manifest, with optional external-source hashes."""
from pathlib import Path
import hashlib,json,sys
D=Path(__file__).resolve().parent
claims=json.loads((D/'OBJECT_CONTRACT.json').read_text())
reg=json.loads((D/'CLAIM_REGISTER.json').read_text())
assert claims['task']=='BR-001' and claims['decision']=='OPEN_INPUT'
assert claims['object_supplied'] is False and claims['admitted_object_id'] is None
assert claims['do_not_execute']==['BR-002','BR-003']
assert reg['scientific_claims_added']==0 and reg['status']=='OPEN_INPUT'
assert len(reg['claims'])==4
for line in (D/'MANIFEST.sha256').read_text().splitlines():
    digest,name=line.split('  ',1); p=D/name
    assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==digest,name
for textfile in ['HANDOFF.md','ADMISSION_OR_OPEN_INPUT.md','REPORT.md']:
    t=(D/textfile).read_text()
    assert 'OPEN_INPUT' in t and 'BR-002' in t and 'BR-003' in t
print('PASS: BR-001 decision, deferred gates, package hashes, basic report consistency')
if '--sources' in sys.argv:
    root=D.parent
    for line in (D/'INPUTS.sha256').read_text().splitlines():
        digest,rem=line.split('  ',1)
        rel=rem.split('    # ',1)[0]
        p=root/rel
        assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==digest,rel
    print('PASS: available external source hashes')
