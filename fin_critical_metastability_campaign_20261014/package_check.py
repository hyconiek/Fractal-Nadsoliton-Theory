#!/usr/bin/env python3
"""Check source and stage contracts. Can optionally verify manifest if present."""
from pathlib import Path
import hashlib,json
R=Path(__file__).resolve().parent
checks={}
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for i in range(1,6):
 d=R/f'MT-{i:03d}'
 for n in ['REPORT.md','RESULTS.json','PROOF_OR_OPEN_LEMMA.md','NONCONCLUSIONS.md','REPLAY.md','INPUTS.sha256']:
  checks[f'{d.name}/{n}']= (d/n).is_file()
 for line in (d/'INPUTS.sha256').read_text().splitlines():
  h,p=line.split('  ',1)
  f=(d/p).resolve()
  checks[f'{d.name} source hash {p}']= f.is_file() and sha(f)==h
for n in ['HANDOFF.md','CLAIM_REGISTER.json','NEXT_OPEN_ATOM.md','NONCONCLUSIONS.md','REPLAY.md','RESEARCH_SUMMARY_PL.md','PROVENANCE.json','CLEAN_REPLAY_VERIFICATION.json','INDEPENDENT_VERIFICATION.json','NUMERICAL_DIAGNOSTICS.json']:
 checks[f'root {n}']=(R/n).is_file()
checks['nominal 57/57']=json.loads((R/'NUMERICAL_DIAGNOSTICS.json').read_text())['checks_passed']==57
checks['independent 58/58']=json.loads((R/'INDEPENDENT_VERIFICATION.json').read_text())['passed']==58
checks['clean replay 2/2']=json.loads((R/'CLEAN_REPLAY_VERIFICATION.json').read_text())['byte_identical_outputs']==2
if (R/'MANIFEST.sha256').exists():
 for line in (R/'MANIFEST.sha256').read_text().splitlines():
  h,p=line.split('  ',1)
  checks[f'manifest {p}']=(R/p).is_file() and sha(R/p)==h
print('PACKAGE CHECK',sum(checks.values()),'/',len(checks))
if not all(checks.values()):
 print('FAILS', [k for k,v in checks.items() if not v]); raise SystemExit(1)
