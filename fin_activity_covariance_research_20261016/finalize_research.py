#!/usr/bin/env python3
"""Mechanical research/handoff packaging; never edits proofs or prior inputs."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
ORDER=ROOT/'FIN_STATIONARY_KINETIC_RIGIDITY_RESEARCH_ORDER_20261016.md'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def reject(x):raise ValueError(x)

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--check',action='store_true'); args=parser.parse_args()
    required=['THEOREMS.md','REPORT_PL.md','ACTIVITY_PROVIDER_CONTRACT.json','CLAIM_REGISTER.json',
        'ADMISSION_OR_OPEN_INPUT.json','PREMISE_LEDGER.md','DEPENDENCY_DELTA.md',
        'HANDOFF.md','NONCONCLUSIONS.md','PROPOSED_AGENTS_PATCH.md','REPLAY.md','verify_research.py','RESULTS.json','finalize_research.py']
    assert all((HERE/n).is_file() for n in required)
    for f in HERE.glob('*.json'):json.loads(f.read_text(),parse_constant=reject)
    for f in HERE.glob('*.py'):ast.parse(f.read_text(),filename=str(f))
    for f in [*HERE.glob('*.md'),ORDER]:
        assert all(line==line.rstrip() for line in f.read_text().splitlines()),f
        for target in re.findall(r'\[[^\]]+\]\(([^)]+)\)',f.read_text()):
            if not target.startswith(('http://','https://','#')):assert (f.parent/target.split('#',1)[0]).exists(),(f,target)
    saved=json.loads((HERE/'RESULTS.json').read_text())
    assert saved['status']=='PASS_BOUNDED_ACTIVITY_COVARIANCE_RESEARCH_CHECKS'
    assert saved['count']==len(saved['checks'])==22 and all(c['passed'] for c in saved['checks'])
    assert saved['implementation_sha256']==sha(HERE/'verify_research.py')
    for rel,h in saved['inputs'].items():assert sha(ROOT/rel)==h,rel
    # Preserve the complete current rate review, not just a selected source file.
    prior=ROOT/'fin_rate_coupling_review_20261016'
    for row in (prior/'MANIFEST.sha256').read_text().splitlines():
        h,n=row.split('  ',1); assert sha(prior/n)==h,n
    cards=ORDER.read_text().split('### KR-')[1:]; assert len(cards)==4
    for i,card in enumerate(cards,1):assert card.startswith(f'{i:03d} '),i
    contract=json.loads((HERE/'ACTIVITY_PROVIDER_CONTRACT.json').read_text())
    claims=json.loads((HERE/'CLAIM_REGISTER.json').read_text())
    assert contract['AS_002']==contract['AS_003']=='NOT_ACTIVATED'
    assert claims['independent_activity_source']=='NOT_DERIVED'
    result={'date':'2026-10-16','status':'PASS_MECHANICAL_RESEARCH_AND_HANDOFF_CHECKS',
        'required_files':required,'bounded_checks':22,'order_sha256':sha(ORDER),
        'KR_cards':4,'source_proof_status':'analytical, not formal machine proof',
        'next_class_process_theorem':'NOT_EXECUTED','source_bridge':'OPEN; AS-002/003 NOT_ACTIVATED'}
    input_rows=''.join(f'{h}  {rel}\n' for rel,h in sorted(saved['inputs'].items()))
    if args.check:
        assert (HERE/'INPUTS.sha256').read_text()==input_rows
        assert result==json.loads((HERE/'PACKAGE_CHECKS.json').read_text())
        names=[]
        for row in (HERE/'MANIFEST.sha256').read_text().splitlines():
            h,n=row.split('  ',1); assert n not in names and sha(HERE/n)==h,n; names.append(n)
        assert set(names)=={f.name for f in HERE.iterdir() if f.is_file() and f.name!='MANIFEST.sha256'}
    else:
        (HERE/'INPUTS.sha256').write_text(input_rows)
        (HERE/'PACKAGE_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
        files=sorted(f for f in HERE.iterdir() if f.is_file() and f.name!='MANIFEST.sha256')
        (HERE/'MANIFEST.sha256').write_text(''.join(f'{sha(f)}  {f.name}\n' for f in files))
    print(json.dumps({'status':result['status'],'read_only':args.check,'manifest_entries':len(required)+2,'KR_cards':4},indent=2))

if __name__=='__main__':main()
