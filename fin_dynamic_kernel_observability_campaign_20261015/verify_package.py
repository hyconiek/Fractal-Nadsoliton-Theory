#!/usr/bin/env python3
"""Reproducible package completeness, manifest, JSON and contract checks.
This does not certify the analytical DKO/MT/KO theorem hypotheses.
"""
from pathlib import Path
import hashlib, json, sys
R=Path(__file__).resolve().parent
M=R/'MANIFEST.sha256'
req_roots=['HANDOFF.md','REPORT.md','CLAIM_REGISTER.json','NONCONCLUSIONS.md','NEXT_OPEN_ATOM.md','REPLAY.md','SOURCE_INPUTS.sha256','PROPOSED_AGENTS_PATCH.md','RESEARCH_SUMMARY_PL.md','PROVENANCE.json','INDEPENDENT_VERIFICATION.json']
req_stage=['REPORT.md','RESULTS.json','PROOF_OR_OPEN_LEMMA.md','NONCONCLUSIONS.md','REPLAY.md','INPUTS.sha256']
errors=[]
for p in req_roots:
 if not (R/p).is_file():errors.append('missing root '+p)
for k in range(1,6):
 for name in req_stage:
  file=f'DKO-{k:03d}/{name}'
  if not (R/file).is_file():errors.append('missing stage '+file)
rows=[]
for row in M.read_text().splitlines():
 if not row.strip():continue
 try:h,name=row.split('  ',1)
 except ValueError:errors.append('bad manifest line '+row);continue
 p=R/name
 if not p.is_file():errors.append('missing manifest target '+name)
 elif hashlib.sha256(p.read_bytes()).hexdigest()!=h:errors.append('hash mismatch '+name)
 rows.append(name)
actual={str(p.relative_to(R)) for p in R.rglob('*') if p.is_file() and p.name not in {'MANIFEST.sha256','PACKAGE_VERIFICATION.json'}}
if set(rows)!=actual:errors.append('manifest coverage error: '+str(sorted(actual-set(rows)))[:1000]+str(sorted(set(rows)-actual))[:1000])
json_parsed=0
for p in R.rglob('*.json'):
 if p.name=='PACKAGE_VERIFICATION.json':continue
 try:json.loads(p.read_text());json_parsed+=1
 except Exception as e: errors.append('bad json '+str(p.relative_to(R))+': '+str(e))
checks=json.loads((R/'INDEPENDENT_VERIFICATION.json').read_text()).get('n_checks')
if checks!=32:errors.append('independent verification n_checks !=32; got '+str(checks))
res={'status':'PASS' if not errors else 'FAIL','manifest_entries':len(rows),'files_without_manifests':len(actual),'parsed_json':json_parsed,'independent_bounded_assertions':checks,'stage_count':5,'errors':errors,'claim_ceiling':'NOT a mechanized proof of stochastic theorem or real-world apparatus'}
(R/'PACKAGE_VERIFICATION.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
if errors:sys.exit(1)
