#!/usr/bin/env python3
"""Mechanical combined-review checks, not an all-N cut or stochastic proof."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
PLAN=ROOT/'FIN_SIGNED_INTERFACE_ORBIT_CUT_NEXT_CHAIN_20261016.md'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def reject(value):raise ValueError(value)
def agents_prefix(expected):
    digest=hashlib.sha256()
    for line in (ROOT/'AGENTS.md').read_bytes().splitlines(keepends=True):
        digest.update(line)
        if digest.hexdigest()==expected:return True
    return False

def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args()
    required=['INTAKE_20261016.md','PROOF_AUDIT_AND_CORRECTIONS.md',
        'SIGNED_INTERFACE_CUT_DERIVATIONS.md','CRITICAL_DELETION_DRIFT_OBSTRUCTION.md',
        'CUT_DECISION.md','RESEARCH_SUMMARY_PL.md','CLAIM_REGISTER.json','INPUT_REGISTRY.json',
        'NONCONCLUSIONS.md','HANDOFF.md','REPLAY.md','verify_intake.py','verification.json','finalize_review.py']
    assert all((HERE/name).is_file() for name in required)
    for path in HERE.glob('*.json'):json.loads(path.read_text(),parse_constant=reject)
    for path in HERE.glob('*.py'):ast.parse(path.read_text(),filename=str(path))
    for path in [*HERE.glob('*.md'),PLAN]:
        assert all(line==line.rstrip() for line in path.read_text().splitlines()),path
        for target in re.findall(r'\[[^\]]+\]\(([^)]+)\)',path.read_text()):
            if not target.startswith(('http://','https://','#')):assert (path.parent/target.split('#',1)[0]).exists(),(path,target)
    registry=json.loads((HERE/'INPUT_REGISTRY.json').read_text())
    for key,rel in registry.items():
        if key!='note':assert (ROOT/rel).is_dir(),rel
    saved=json.loads((HERE/'verification.json').read_text())
    assert saved['status']=='PASS_SCOPED_RB_AND_WEIGHTED_DEGREE_INTAKE'
    assert len(saved['checks'])==22 and all(c['passed'] for c in saved['checks'])
    assert saved['RB_manifest']==saved['RB_manifest_after'] and saved['bundle_manifest']==saved['bundle_manifest_after']
    assert saved['RB_manifest']['matched']==47 and saved['bundle_manifest']['matched']==75
    assert sum(v['matched_members'] for v in saved['ZIP_member_binding'].values())==71
    replay=saved['replay']['campaigns'];assert all(v['scientific_output_byte_equal'] for v in replay.values())
    assert replay['repair']['logs']['verify_package.py']['exit_code']==1
    assert replay['repair']['logs']['verify_package.py']['independent_manifest_fallback']['matched']==17
    assert all(v['logs']['verify_package.py']['exit_code']==0 for k,v in replay.items() if k!='repair')
    for rel,h in saved['input_hashes'].items():
        assert agents_prefix(h) if rel=='AGENTS.md' else sha(ROOT/rel)==h,rel
    for key in ('RB','shell','repair','orbit','degree'):
        folder=ROOT/registry[key]
        for row in (folder/'MANIFEST.sha256').read_text().splitlines():
            h,rel=row.split('  ',1);assert sha(folder/rel)==h,(key,rel)
    for row in (ROOT/registry['bundle']/'MANIFEST_SHA256.txt').read_text().splitlines():
        h,rel=row.split('  ',1);assert sha(ROOT/registry['bundle']/rel)==h,rel
    cards=PLAN.read_text().split('### SC-')[1:]
    assert len(cards)==4 and all(card.startswith(f'{i:03d} ') for i,card in enumerate(cards,1))
    assert 'SC-001–SC-004 NOT EXECUTED' in PLAN.read_text()
    assert 'FIN combined RB / RBCUT intake and signed-interface cut frontier' in (ROOT/'AGENTS.md').read_text()
    result={'date':'2026-10-16','status':'PASS_MECHANICAL_COMBINED_RB_REVIEW_AND_SC_PLAN',
        'required_files':required,'plan_sha256':sha(PLAN),'SC_cards':4,'bounded_intake_checks':22,
        'source_preservation':'verified original source and outer bundle manifests',
        'uniform_actual_cut_bound':'OPEN','actual_rates':'NOT_CONSTRUCTED','universal_IBA':'OPEN',
        'scope':'mechanical JSON/AST/link/hash/replay-disposition checks; not a mathematical proof'}
    if args.check:
        assert result==json.loads((HERE/'PACKAGE_CHECKS.json').read_text())
        names=[]
        for row in (HERE/'MANIFEST.sha256').read_text().splitlines():
            h,name=row.split('  ',1);assert name not in names and sha(HERE/name)==h,name;names.append(name)
        assert set(names)=={p.name for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256'}
    else:
        (HERE/'PACKAGE_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
        files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256')
        (HERE/'MANIFEST.sha256').write_text(''.join(f'{sha(path)}  {path.name}\n' for path in files))
        names=[path.name for path in files]
    print(json.dumps({'status':result['status'],'read_only':args.check,'manifest_entries':len(names),'SC_cards':4},indent=2))

if __name__=='__main__':main()
