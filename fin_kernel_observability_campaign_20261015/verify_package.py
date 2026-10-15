#!/usr/bin/env python3
"""Verify deterministic KO source archive file coverage and local task layout."""
from pathlib import Path
import hashlib,json,sys
R=Path(__file__).resolve().parent
MAN='MANIFEST.sha256'
SKIP={MAN,'PACKAGE_VERIFICATION.json'}
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
files=sorted(x for x in R.rglob('*') if x.is_file() and '__pycache__' not in x.parts)
paths=[str(x.relative_to(R)) for x in files]
assert MAN in paths
entries={}
for line in (R/MAN).read_text().splitlines():
 if not line.strip():continue
 h,name=line.split('  ',1)
 assert name not in entries,(name,'duplicate')
 entries[name]=h
actual={str(p.relative_to(R)):digest(p) for p in files if str(p.relative_to(R)) not in SKIP}
assert entries==actual,{'only_manifest':sorted(set(entries)-set(actual)),'only_actual':sorted(set(actual)-set(entries)),'mismatched':[x for x in entries if x in actual and entries[x]!=actual[x]]}
for n in ['HANDOFF.md','REPORT.md','NONCONCLUSIONS.md','NEXT_OPEN_ATOM.md','REPLAY.md','SOURCE_INPUTS.sha256','PROVENANCE.json','CLAIM_REGISTER.json','INDEPENDENT_VERIFICATION.json','PROPOSED_AGENTS_PATCH.md']:
 assert n in actual,n
for i in range(1,6):
 t=f'KO-{i:03}'
 for n in ['REPORT.md','RESULTS.json','PROOF_OR_OPEN_LEMMA.md','NONCONCLUSIONS.md','REPLAY.md','INPUTS.sha256']:
  assert f'{t}/{n}' in actual, f'{t}/{n}'
assert json.loads((R/'INDEPENDENT_VERIFICATION.json').read_text())['assertions_passed']==38
assert json.loads((R/'KO-004/RESULTS.json').read_text())['scientific_verdict']=='DEFERRED_NOT_REQUIRED'
result={'status':'PASS','manifest_file_count':len(entries),'task_cards':5,'self_consistency_assertions_passed':38,'proof_scope':'synthetic finite/exact rational/interval calculations; not a mechanized limit theorem'}
(R/'PACKAGE_VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
