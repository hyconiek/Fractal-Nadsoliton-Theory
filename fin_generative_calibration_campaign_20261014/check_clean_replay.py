#!/usr/bin/env python3
"""Clean workspace replay, numerical/byte reproduction in this environment."""
import tempfile,shutil,subprocess,hashlib,json,os
from pathlib import Path
R=Path(__file__).resolve().parent
SELECT=['RG-F01/PROCESS_CONTRACT.json','RG-F01/COMPOSITION_TESTS.json',
 'RG-F02/SYMMETRY_AND_ERASURE_CONTROLS.json','RG-F03/TOTAL_LAW_ERROR_BUDGET.json',
 'RG-F03/PREDICTION_ENCLOSURES.json','RG-F04/DECISION_RULE.json',
 'RG-F04/SYNTHETIC_COMPOSITE_TESTS.json','INDEPENDENT_VERIFICATION.json']
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
with tempfile.TemporaryDirectory(prefix='fin-f-clean-') as temp:
 copy=Path(temp)/R.name
 shutil.copytree(R,copy,ignore=shutil.ignore_patterns('__pycache__','.git'))
 env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
 for script in ['run_calibration.py','decision_rule.py','independent_verify.py']:
  p=subprocess.run(['python',script],cwd=copy,env=env,capture_output=True,text=True,check=True)
  print(script,':',p.stdout.strip()[-140:])
 hashes={s:{'original_sha256':digest(R/s),'clean_sha256':digest(copy/s),
             'byte_match':digest(R/s)==digest(copy/s)} for s in SELECT}
 status={'number_of_compared_outputs':len(SELECT),'all_sha256_equal':all(v['byte_match'] for v in hashes.values()),
   'within_environment_replay_only':True,'cross_environment_identical_bytes_not_claimed':True,'file_hashes':hashes}
 (R/'CLEAN_REPLAY_VERIFICATION.json').write_text(json.dumps(status,indent=2,sort_keys=True)+'\n')
 print('clean replay:',status['number_of_compared_outputs'],'/',sum(v['byte_match'] for v in hashes.values()),'byte-identical')
 if not status['all_sha256_equal']:raise SystemExit('clean replay FAILED')
