#!/usr/bin/env python3
"""Mechanical SC integration/research-order checks; not a proof checker."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
RESEARCH=ROOT/'fin_balanced_phase_cut_research_20261016'
PLAN=ROOT/'FIN_BALANCED_PHASE_CURRENT_NEXT_CHAIN_20261016.md'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def reject(value):raise ValueError(value)
def prefix(h):
    d=hashlib.sha256()
    for line in (ROOT/'AGENTS.md').read_bytes().splitlines(keepends=True):
        d.update(line)
        if d.hexdigest()==h:return True
    return False
def manifest(folder):
    names=[]
    for row in (folder/'MANIFEST.sha256').read_text().splitlines():
        h,n=row.split('  ',1);assert n not in names and sha(folder/n)==h,(folder,n);names.append(n)
    assert set(names)=={p.name for p in folder.iterdir() if p.is_file() and p.name!='MANIFEST.sha256'}
    return len(names)

def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args()
    required=['INTAKE_20261016.md','PROOF_AUDIT.md','CURRENT_DECISION.md','CLAIM_REGISTER.json',
        'NONCONCLUSIONS.md','REPLAY.md','HANDOFF.md','COMPLETION_AUDIT.md',
        'verify_intake.py','verification.json','finalize_review.py']
    assert all((HERE/n).is_file() for n in required)
    for path in HERE.glob('*.json'):json.loads(path.read_text(),parse_constant=reject)
    for path in HERE.glob('*.py'):ast.parse(path.read_text(),filename=str(path))
    for path in [*HERE.glob('*.md'),PLAN]:
        assert all(x==x.rstrip() for x in path.read_text().splitlines()),path
        for target in re.findall(r'\[[^\]]+\]\(([^)]+)\)',path.read_text()):
            if not target.startswith(('http:','https:','#')):assert (path.parent/target).exists(),(path,target)
    saved=json.loads((HERE/'verification.json').read_text())
    assert saved['status']=='PASS_SCOPED_SC_INTAKE_WITH_NEW_BALANCED_PHASE_RESEARCH'
    assert len(saved['checks'])==11 and all(c['passed'] for c in saved['checks'])
    assert saved['source_manifest']==saved['source_manifest_after'] and saved['source_manifest']['matched']==48
    assert saved['replay']['byte_equal'] and saved['replay']['producer_assertions']==148
    assert len(saved['source_inputs']['matched'])==10 and len(saved['source_inputs']['unavailable'])==2
    for rel,h in saved['input_hashes'].items():assert prefix(h) if rel=='AGENTS.md' else sha(ROOT/rel)==h,rel
    source=ROOT/'fin_signed_interface_orbit_cut_campaign_20261016';source_names=[]
    for row in (source/'MANIFEST.sha256').read_text().splitlines():
        h,name=row.split('  ',1);assert name not in source_names and sha(source/name)==h,name;source_names.append(name)
    assert len(source_names)==48 and set(source_names)=={str(p.relative_to(source)) for p in source.rglob('*') if p.is_file() and p.name!='MANIFEST.sha256'}
    research_count=manifest(RESEARCH)
    data=json.loads((RESEARCH/'RESULTS.json').read_text());assert len(data['checks'])==20 and all(c['passed'] for c in data['checks'])
    for rel,h in data['input_hashes'].items():assert sha(ROOT/rel)==h,rel
    cards=PLAN.read_text().split('### PC-')[1:]
    assert len(cards)==4 and all(card.startswith(f'{i:03d} ') for i,card in enumerate(cards,1))
    assert 'PC-001–PC-004 NOT EXECUTED' in PLAN.read_text()
    assert '## FIN SC intake and explicit balanced-phase current frontier' in (ROOT/'AGENTS.md').read_text()
    result={'started':'2026-10-16','final_validation':'2026-10-17',
        'status':'PASS_MECHANICAL_SC_REVIEW_RESEARCH_AND_PC_ORDER','required_files':required,
        'plan_sha256':sha(PLAN),'PC_cards':4,'research_manifest_entries':research_count,
        'source_manifest_entries':48,'source_assertions':148,'intake_controls':11,'new_research_controls':20,
        'source_preservation':'actual before/after manifest verified',
        'coupling':'OPEN','full_generator':'NOT_CONSTRUCTED','unconditional_IBA':'OPEN',
        'scope':'mechanical evidence packaging, not proof by tests or physical validation'}
    if args.check:
        assert result==json.loads((HERE/'PACKAGE_CHECKS.json').read_text())
        count=manifest(HERE)
    else:
        (HERE/'PACKAGE_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
        files=sorted(path for path in HERE.iterdir() if path.is_file() and path.name!='MANIFEST.sha256')
        (HERE/'MANIFEST.sha256').write_text(''.join(f'{sha(path)}  {path.name}\n' for path in files));count=len(files)
    print(json.dumps({'status':result['status'],'read_only':args.check,'manifest_entries':count,'PC_cards':4},indent=2))

if __name__=='__main__':main()
