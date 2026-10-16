#!/usr/bin/env python3
"""Mechanical BA review/FP plan packaging; --check is read-only."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
PLAN=ROOT/'FIN_COUNT_PARTITION_FEASIBILITY_NEXT_CHAIN_20261016.md'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def reject(x):raise ValueError(x)

def main():
    p=argparse.ArgumentParser(); p.add_argument('--check',action='store_true'); a=p.parse_args()
    required=['INTAKE_20261016.md','PROOF_AUDIT_AND_SCOPE.md','COUNT_PARTITION_FEASIBILITY_SEED.md',
        'PARTITION_DECISION.md','RESEARCH_SUMMARY_PL.md','CLAIM_REGISTER.json','NONCONCLUSIONS.md',
        'HANDOFF.md','REPLAY.md','verify_intake.py','verification.json','finalize_review.py']
    assert all((HERE/n).is_file() for n in required)
    for path in HERE.glob('*.json'):json.loads(path.read_text(),parse_constant=reject)
    for path in HERE.glob('*.py'):ast.parse(path.read_text(),filename=str(path))
    for path in [*HERE.glob('*.md'),PLAN]:
        assert all(line==line.rstrip() for line in path.read_text().splitlines()),path
        for target in re.findall(r'\[[^\]]+\]\(([^)]+)\)',path.read_text()):
            if not target.startswith(('http://','https://','#')):assert (path.parent/target.split('#',1)[0]).exists(),(path,target)
    saved=json.loads((HERE/'verification.json').read_text())
    assert saved['status']=='PASS_SCOPED_BA_INTAKE_CHECKS' and len(saved['checks'])==22 and all(c['passed'] for c in saved['checks'])
    assert saved['source_manifest']==saved['source_manifest_after'] and saved['source_manifest']['matched']==52
    assert saved['replay']['byte_equal']
    for rel,h in saved['input_hashes'].items():assert sha(ROOT/rel)==h,rel
    for folder in ['fin_reversible_bracket_averaging_campaign_20261016',*saved['actual_dependencies']]:
        for row in (ROOT/folder/'MANIFEST.sha256').read_text().splitlines():
            h,n=row.split('  ',1); assert sha(ROOT/folder/n)==h,(folder,n)
    for n,h in saved['source_manifest']['excluded_hashes'].items():assert sha(ROOT/'fin_reversible_bracket_averaging_campaign_20261016'/n)==h
    cards=PLAN.read_text().split('### FP-')[1:]; assert len(cards)==4
    fields=['Question','Inputs/dependencies','Method','Deliverables','Acceptance','Failure/stop','Must not conclude','Priority/class']
    for i,card in enumerate(cards,1):assert card.startswith(f'{i:03d} ') and all(field in card for field in fields),i
    result={'date':'2026-10-16','status':'PASS_MECHANICAL_BA_REVIEW_AND_FP_PLAN_CHECKS',
        'required_files':required,'plan_sha256':sha(PLAN),'FP_cards':4,
        'source_universal_IBA':'OPEN','closed_partition_rates':'NOT_CONSTRUCTED',
        'scope':'mechanical syntax/link/hash/scope checks, not formal proof or FP execution'}
    if a.check:
        assert result==json.loads((HERE/'PACKAGE_CHECKS.json').read_text())
        names=[]
        for row in (HERE/'MANIFEST.sha256').read_text().splitlines():
            h,n=row.split('  ',1); assert n not in names and sha(HERE/n)==h,n; names.append(n)
        assert set(names)=={path.name for path in HERE.iterdir() if path.is_file() and path.name!='MANIFEST.sha256'}
    else:
        (HERE/'PACKAGE_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
        files=sorted(path for path in HERE.iterdir() if path.is_file() and path.name!='MANIFEST.sha256')
        (HERE/'MANIFEST.sha256').write_text(''.join(f'{sha(path)}  {path.name}\n' for path in files))
    print(json.dumps({'status':result['status'],'read_only':a.check,'manifest_entries':len(required)+1,'FP_cards':4},indent=2))

if __name__=='__main__':main()
