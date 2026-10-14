#!/usr/bin/env python3
"""Manifest and output structural checks; not a theorem prover."""
import hashlib,json,subprocess
from pathlib import Path
r=Path(__file__).resolve().parent
checks=[]
def ck(name,predicate):checks.append({'name':name,'passed':bool(predicate)})
for task in ['KI-001','KI-002','KI-003','KI-004','KI-005']:
    p=r/task
    for fn in ['REPORT.md','RESULTS.json','PROOF_OR_OPEN_LEMMA.md','NONCONCLUSIONS.md','REPLAY.md','INPUTS.sha256']:
        ck(task+'/'+fn,(p/fn).exists() and (p/fn).stat().st_size>30)
    try:
        ob=json.loads((p/'RESULTS.json').read_text());ck(task+' status',ob.get('execution_status')=='EXECUTED' and not ob.get('fitted_objects'))
    except:ck(task+' status',False)
res=json.loads((r/'INDEPENDENT_VERIFICATION.json').read_text())
ck('independent check flags',res['all_pass'] and res['passed']==res['total']==41)
k=json.loads((r/'KI-002'/'DERIVED_SOFT_CLOCKS.json').read_text())
ck('mobilities positive',all(m['mobility_m_nominal']>0 for m in k['models'].values()))
ck('numerical detailed balance',k['finite_detailed_balance_max_log_defect_nominal']<5e-12)
ck('nontrivial kinetic clocks',k['models']['HB']['mobility_m_nominal']>k['models']['BAR']['mobility_m_nominal']>k['models']['MET']['mobility_m_nominal'])
ck('local diagnostic present',(r/'KI-003'/'LOCAL_DRIFT_DIAGNOSTICS.json').exists())
manifest=r/'MANIFEST.sha256'
if manifest.exists():
    ok=True; cnt=0
    for line in manifest.read_text().splitlines():
        if not line.strip():continue
        sha,name=line.split('  ',1);p=r/name
        if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=sha:ok=False
        cnt+=1
    ck('manifest integrity',ok)
else:cnt=0;ck('manifest integrity',False)
summary={'checks':checks,'passed':sum(x['passed'] for x in checks),'total':len(checks),'all_pass':all(x['passed'] for x in checks),'manifest_entries':cnt}
print(json.dumps({k:summary[k] for k in ['passed','total','all_pass','manifest_entries']},indent=2))
if not summary['all_pass']:print([x['name'] for x in checks if not x['passed']]);raise SystemExit(1)
