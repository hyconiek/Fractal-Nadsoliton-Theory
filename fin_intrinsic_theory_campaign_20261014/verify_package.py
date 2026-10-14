from pathlib import Path
import json, hashlib, sys
R=Path(__file__).resolve().parent
required=['HANDOFF.md','CLAIM_REGISTER.json','NONCONCLUSIONS.md','REPLAY.md','REMOTE_SOURCES.json','NUMERICAL_REPLAY.json','INDEPENDENT_VERIFICATION.json','theory_core.py','fold_diagnostics.py','independent_verify.py']
for t in range(1,6):
    d=R/f'TH-{t:03d}'
    for f in ['REPORT.md','RESULTS.json','PROOF_OR_OPEN_LEMMA.md','NONCONCLUSIONS.md','REPLAY.md','INPUTS.sha256']:
        required.append(str((d/f).relative_to(R)))
missing=[x for x in required if not (R/x).exists()]
json_bad=[]
for p in R.rglob('*.json'):
    try: json.load(open(p))
    except Exception as e: json_bad.append([str(p.relative_to(R)),repr(e)])
statuses={}
for t in range(1,6): statuses[f'TH-{t:03d}']=json.load(open(R/f'TH-{t:03d}'/'RESULTS.json'))['execution_status']
ind=json.load(open(R/'INDEPENDENT_VERIFICATION.json'))
out={'status':'PASS' if not missing and not json_bad and ind['status']=='PASS' else 'FAIL','missing':missing,'invalid_json':json_bad,'task_execution_statuses':statuses,'independent_verification':f"{ind['passed']}/{ind['total']}"}
(R/'PACKAGE_VERIFICATION.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2))
if out['status']!='PASS': sys.exit(1)
