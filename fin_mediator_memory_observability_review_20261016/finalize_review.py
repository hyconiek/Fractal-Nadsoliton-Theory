#!/usr/bin/env python3
"""Mechanical review packaging, with read-only --check; not a proof checker."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
PLAN=ROOT/'FIN_POST_MM_KINETIC_BRIDGE_ACTIVATION_CHAIN_20261016.md'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def reject(x):raise ValueError(x)

def main():
    p=argparse.ArgumentParser(); p.add_argument('--check',action='store_true'); a=p.parse_args()
    required=['INTAKE_20261016.md','PROOF_COMPLETIONS_AND_CORRECTIONS.md','BRIDGE_DECISION.md',
        'RESEARCH_SUMMARY_PL.md','NONCONCLUSIONS.md','HANDOFF.md','REPLAY.md','CLAIM_REGISTER.json',
        'certify_margin.py','MARGIN_CERTIFICATE.json','verify_intake.py','verification.json','finalize_review.py']
    assert all((HERE/n).is_file() for n in required)
    for f in HERE.glob('*.json'):json.loads(f.read_text(),parse_constant=reject)
    for f in HERE.glob('*.py'):ast.parse(f.read_text(),filename=str(f))
    for f in [*HERE.glob('*.md'),PLAN]:
        assert all(line==line.rstrip() for line in f.read_text().splitlines()),f
        for target in re.findall(r'\[[^\]]+\]\(([^)]+)\)',f.read_text()):
            if not target.startswith(('http://','https://','#')):assert (f.parent/target.split('#',1)[0]).exists(),(f,target)
    saved=json.loads((HERE/'verification.json').read_text())
    assert saved['status']=='PASS_SCOPED_MM_INTAKE_CHECKS' and len(saved['checks'])==28 and all(c['passed'] for c in saved['checks'])
    assert saved['source_manifest']['matched']==58 and saved['replay']['source_flags']==35
    assert saved['source_manifest']==saved['source_manifest_after']
    for rel,h in saved['input_hashes'].items():assert sha(ROOT/rel)==h,rel
    source=ROOT/'fin_mediator_memory_observability_campaign_20261016'
    for line in (source/'MANIFEST.sha256').read_text().splitlines():
        h,name=line.split('  ',1); assert sha(source/name)==h,name
    for name,h in saved['source_manifest']['excluded_hashes'].items():assert sha(source/name)==h,name
    from certify_margin import build
    assert build()==json.loads((HERE/'MARGIN_CERTIFICATE.json').read_text())
    cards=PLAN.read_text().split('### KB-')[1:]
    fields=['Title','Scientific question','Why it matters','Exact inputs','Dependencies','Method','Deliverables',
        'Acceptance criterion','Failure/refutation criterion','What becomes possible','Must NOT conclude','Estimated computational class','Priority']
    assert len(cards)==3
    for i,card in enumerate(cards,1):assert card.startswith(f'{i:03d} ') and all(f in card for f in fields),i
    result={'date':'2026-10-16','status':'PASS_MECHANICAL_REVIEW_AND_CONDITIONAL_PLAN_CHECKS',
        'required_files':required,'KB_cards':3,'plan_sha256':sha(PLAN),'json_ast_links_whitespace':'PASS',
        'scope':'Mechanical packaging/dependency checks, not formal stochastic proof or KB execution'}
    if a.check:
        assert result==json.loads((HERE/'PACKAGE_CHECKS.json').read_text())
        names=[]
        for row in (HERE/'MANIFEST.sha256').read_text().splitlines():
            h,name=row.split('  ',1); assert name not in names and sha(HERE/name)==h,name; names.append(name)
        assert set(names)=={f.name for f in HERE.iterdir() if f.is_file() and f.name!='MANIFEST.sha256'}
    else:
        (HERE/'PACKAGE_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
        files=sorted(f for f in HERE.iterdir() if f.is_file() and f.name!='MANIFEST.sha256')
        (HERE/'MANIFEST.sha256').write_text(''.join(f'{sha(f)}  {f.name}\n' for f in files))
    print(json.dumps({'status':result['status'],'read_only':a.check,'manifest_entries':len(required)+1,'conditional_KB_cards':3},indent=2))

if __name__=='__main__':main()
