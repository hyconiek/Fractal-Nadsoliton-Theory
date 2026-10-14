#!/usr/bin/env python3
"""Check provenance, hash every handoff source and output, build self-contained ZIP."""
from pathlib import Path
import json,hashlib,zipfile
R=Path(__file__).resolve().parent
DEST=Path('/mnt/data/FIN_RG_F01_F05_FULL_HANDOFF_20261014.zip')
def sha(data):return hashlib.sha256(data).hexdigest()
prior=R/'source_context/FIN_RG_E01_E05_FULL_HANDOFF_20261014.zip'
with zipfile.ZipFile(prior) as archive:
 base='fin_generative_operational_campaign_20261014/'
 manifest=archive.read(base+'MANIFEST.sha256').decode().splitlines()
 matched=0;errors=[]
 for line in manifest:
  hash_expected,rel=line.split('  ',1)
  fname=base+rel
  try:match=sha(archive.read(fname))==hash_expected
  except KeyError:match=False
  matched+=int(match)
  if not match:errors.append(rel)
 assert matched==len(manifest),(matched,len(manifest),errors)
 assert len(archive.namelist())==len(manifest)+1,(len(archive.namelist()),len(manifest))
 source={'source_package':'user-mounted FIN_RG_E01_E05_FULL_HANDOFF_20261014.zip',
         'source_sha256':sha(prior.read_bytes()),'verified_embedded_source_manifest':f'{matched}/{len(manifest)}',
         'input_plan':'source_context/FIN_N2_CALIBRATION_NEXT_CHAIN_USER.md','input_plan_sha256':sha((R/'source_context/FIN_N2_CALIBRATION_NEXT_CHAIN_USER.md').read_bytes()),
         'source_model':'source_model.py copied byte-for-byte from previously supplied RG-E ZIP, which itself describes it as a formula reimplementation rather than original safe_seed source bytes',
         'repository_baseline_claim':'08b6d235 (user intake; not live git verified)',
         'user_pasted_E_intake':'supplied as chat message; there is no original separately downloaded review/certificate in this F source bundle',
         'gain_cover':'0.00119114 after output TV allowances, inherited conditional user-review statement; proof file GAIN_BOX_CERTIFICATE unavailable locally for independent redo',
         'empirical_records':'NONE; all synthetic test vectors are explicitly labelled hypothetical'}
 (R/'SOURCE_PROVENANCE.json').write_text(json.dumps(source,indent=2,sort_keys=True)+'\n')

# For each task no missing standard artifact.
for i in range(1,6):
 d=R/f'RG-F{i:02}'
 for f in ['REPORT.md','RESULTS.json','NONCONCLUSIONS.md','REPLAY.md','INPUTS.sha256']:
  assert (d/f).is_file(),d/f
 for line in (d/'INPUTS.sha256').read_text().splitlines():
  h,path=line.split('  ',1)
  assert sha((R/path).read_bytes())==h,(d,path)

needed=['HANDOFF.md','CLAIM_REGISTER.json','RESEARCH_SUMMARY_PL.md','NONCONCLUSIONS.md','REPLAY.md','SOURCE_PROVENANCE.json',
        'INDEPENDENT_VERIFICATION.json','CLEAN_REPLAY_VERIFICATION.json','RG-F05/READINESS.json','RG-F04/STATISTICAL_PROOF.md',
        'RG-F04/DECISION_RULE.json','RG-F04/DATA_SCHEMA.json','decision_rule.py','positive_process.py',
        'run_calibration.py','source_context/FIN_RG_E01_E05_FULL_HANDOFF_20261014.zip']
for f in needed:assert (R/f).is_file(),f
v=json.loads((R/'INDEPENDENT_VERIFICATION.json').read_text())
assert v['checks_passed']==v['checks_total'] and v['checks_total']>=90
assert json.loads((R/'CLEAN_REPLAY_VERIFICATION.json').read_text())['all_sha256_equal']

# Manifest covers all actual files except itself, excludes machine-generated caches.
files=sorted([p for p in R.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='MANIFEST.sha256'])
manifest_out=''.join(f'{sha(p.read_bytes())}  {p.relative_to(R).as_posix()}\n' for p in files)
(R/'MANIFEST.sha256').write_text(manifest_out)
files.append(R/'MANIFEST.sha256')
with zipfile.ZipFile(DEST,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
 for p in sorted(files):z.write(p,arcname=R.name+'/'+p.relative_to(R).as_posix())
with zipfile.ZipFile(DEST) as z:
 assert z.testzip() is None
 assert len(z.namelist())==len(files)
 for line in manifest_out.splitlines():
  h,rel=line.split('  ',1)
  assert sha(z.read(R.name+'/'+rel))==h
summary={'zip':str(DEST),'bytes':DEST.stat().st_size,'archive_sha256':sha(DEST.read_bytes()),
         'file_members':len(files),'manifest_entries':len(files)-1,'manifest_all_match':True,
         'previous_source_manifest':f'{matched}/{len(manifest)}','independent_checks':v['checks_total'],
         'clean_replay_files':json.loads((R/'CLEAN_REPLAY_VERIFICATION.json').read_text())['number_of_compared_outputs'],
         'scientific_verdict':'NO_NEW_PHYSICAL_INPUT; DESIGN_ONLY','architect_RG_F06':'NOT_EXECUTED'}
print(json.dumps(summary,indent=2))
