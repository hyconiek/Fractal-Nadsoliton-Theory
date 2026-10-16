#!/usr/bin/env python3
"""Mechanical RC review packaging; --check does not write anything."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re
import subprocess

HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
PLAN=ROOT/'FIN_RATE_ACTIVITY_SOURCE_NEXT_CHAIN_20261016.md'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def reject(x):raise ValueError(x)

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--check',action='store_true'); args=parser.parse_args()
    required=['INTAKE_20261016.md','PROOF_AUDIT_AND_ADDITIONAL_CONTROLS.md','ACTIVITY_DECISION.md',
        'REPORT_PL.md','NONCONCLUSIONS.md','CLAIM_REGISTER.json','HANDOFF.md','REPLAY.md','verify_intake.py','verification.json','finalize_review.py']
    assert all((HERE/n).is_file() for n in required)
    for path in HERE.glob('*.json'):json.loads(path.read_text(),parse_constant=reject)
    for path in HERE.glob('*.py'):ast.parse(path.read_text(),filename=str(path))
    for path in [*HERE.glob('*.md'),PLAN]:
        assert all(line==line.rstrip() for line in path.read_text().splitlines()),path
        for target in re.findall(r'\[[^\]]+\]\(([^)]+)\)',path.read_text()):
            if not target.startswith(('http://','https://','#')):assert (path.parent/target.split('#',1)[0]).exists(),(path,target)
    saved=json.loads((HERE/'verification.json').read_text())
    assert saved['status']=='PASS_SCOPED_RATE_COUPLING_INTAKE_CHECKS' and len(saved['checks'])==17 and all(c['passed'] for c in saved['checks'])
    assert saved['source_manifest']==saved['source_manifest_after'] and saved['source_manifest']['matched']==10
    assert saved['exact_fixtures']==1116 and all(saved['replay']['output_byte_equal'].values())
    for path,h in saved['input_hashes'].items():assert sha(ROOT/path)==h,path
    old=subprocess.run(['git','show',saved['source_inputs']['baseline_AGENTS_verified_at_git']+':AGENTS.md'],cwd=ROOT,capture_output=True,check=True).stdout
    assert hashlib.sha256(old).hexdigest()==saved['source_inputs']['declared']['AGENTS.md']
    assert (ROOT/'AGENTS.md').read_bytes().startswith(old)
    source=ROOT/'fin_rate_coupling_derivation_20261016'
    for line in (source/'MANIFEST.sha256').read_text().splitlines():
        h,n=line.split('  ',1); assert sha(source/n)==h,n
    for rel in saved['previous_reviews']:
        for line in (ROOT/rel/'MANIFEST.sha256').read_text().splitlines():
            h,n=line.split('  ',1); assert sha(ROOT/rel/n)==h,(rel,n)
    cards=PLAN.read_text().split('### AS-')[1:]; assert len(cards)==3
    fields=['Question','Inputs/dependencies','Method','Deliverables','Acceptance','Failure/stop','Must not conclude','Priority/class']
    for i,card in enumerate(cards,1):assert card.startswith(f'{i:03d} ') and all(f in card for f in fields),i
    result={'date':'2026-10-16','status':'PASS_MECHANICAL_REVIEW_AND_CONDITIONAL_PLAN_CHECKS','required_files':required,
        'plan_sha256':sha(PLAN),'conditional_AS_cards':3,'scope':'mechanical syntax/link/hash checks, not formal proof or source admission'}
    if args.check:
        assert result==json.loads((HERE/'PACKAGE_CHECKS.json').read_text())
        names=[]
        for row in (HERE/'MANIFEST.sha256').read_text().splitlines():
            h,n=row.split('  ',1); assert n not in names and sha(HERE/n)==h,n; names.append(n)
        assert set(names)=={p.name for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256'}
    else:
        (HERE/'PACKAGE_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
        files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256')
        (HERE/'MANIFEST.sha256').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in files))
    print(json.dumps({'status':result['status'],'read_only':args.check,'manifest_entries':len(required)+1},indent=2))

if __name__=='__main__':main()
