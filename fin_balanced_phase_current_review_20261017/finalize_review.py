#!/usr/bin/env python3
"""Mechanical corrected-PC intake and physical-search packaging checks."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
PHYSICS=ROOT/'fin_q12_physical_bridge_research_20261017'
PLAN=ROOT/'FIN_PC_PHYSICAL_COMPATIBILITY_NEXT_CHAIN_20261017.md'
SOURCE=ROOT/'fin_balanced_phase_current_campaign_20261016'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def reject(value):raise ValueError(value)
def prefix(h):
    digest=hashlib.sha256()
    for line in (ROOT/'AGENTS.md').read_bytes().splitlines(keepends=True):
        digest.update(line)
        if digest.hexdigest()==h:return True
    return False
def check_manifest(folder,excluded=()):
    names=[]
    for line in (folder/'MANIFEST.sha256').read_text().splitlines():
        h,n=line.split('  ',1);assert n not in names and sha(folder/n)==h,(folder,n);names.append(n)
    actual={str(p.relative_to(folder)) for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    assert actual==set(names)|{'MANIFEST.sha256'}|set(excluded),(folder,actual-set(names))
    return len(names)
def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args()
    required=['INTAKE_20261017.md','PROOF_AUDIT.md','FINITE_ORBIT_HALL_REDUCTION.md',
        'NEXT_DECISION.md','NONCONCLUSIONS.md','REPLAY.md','HANDOFF.md','CLAIM_REGISTER.json',
        'COMPLETION_AUDIT.md','verify_intake.py','verification.json','finalize_review.py']
    assert all((HERE/n).is_file() for n in required)
    for path in HERE.glob('*.json'):json.loads(path.read_text(),parse_constant=reject)
    for path in HERE.glob('*.py'):ast.parse(path.read_text(),filename=str(path))
    for path in [*HERE.glob('*.md'),PLAN]:
        assert all(line==line.rstrip() for line in path.read_text().splitlines()),path
        for target in re.findall(r'\[[^\]]+\]\(([^)]+)\)',path.read_text()):
            if not target.startswith(('http:','https:','#')):assert (path.parent/target).exists(),(path,target)
    saved=json.loads((HERE/'verification.json').read_text())
    assert saved['status']=='PASS_SCOPED_PC_INTAKE_AND_PHYSICAL_BOUNDARY'
    assert len(saved['checks'])==17 and all(c['passed'] for c in saved['checks'])
    assert saved['source_manifest']==saved['source_manifest_after'] and saved['source_manifest']['matched']==26
    assert saved['replay']['byte_equal'] and 'derived' in saved['replay']['wrapper_scope']
    for rel,h in saved['input_hashes'].items():assert prefix(h) if rel=='AGENTS.md' else sha(ROOT/rel)==h,rel
    source_count=check_manifest(SOURCE,('check_package.py','PACKAGE_VERIFICATION.json'))
    for name,h in saved['source_manifest']['excluded_hashes'].items():assert sha(SOURCE/name)==h,name
    physical_count=check_manifest(PHYSICS)
    phy=json.loads((PHYSICS/'RESULTS.json').read_text());assert len(phy['checks'])==16 and all(c['passed'] for c in phy['checks'])
    for rel,h in phy['input_hashes'].items():assert sha(ROOT/rel)==h,rel
    text=PLAN.read_text()
    for name in ('DH-001','DH-002','DH-003','QP-001','QP-002','QP-003'):assert '### '+name+' ' in text,name
    assert 'no task below is executed' in text
    assert '## FIN actual balanced-phase-current PC intake and physical compatibility boundary' in (ROOT/'AGENTS.md').read_text()
    result={'date':'2026-10-17','status':'PASS_MECHANICAL_CORRECTED_PC_INTAKE_AND_BRIDGE_ASSESSMENT',
        'required_files':required,'source_manifest_entries':source_count,'physical_manifest_entries':physical_count,
        'intake_controls':17,'physical_controls':16,'actual_browser_captures':3,'plan_sha256':sha(PLAN),
        'full_Hall':'OPEN','PC003':'NOT_ACTIVATED','physical_source':'NOT_FOUND',
        'scope':'mechanical JSON/AST/link/hash checks and exact bounded replays, not physical validation'}
    if args.check:
        assert result==json.loads((HERE/'PACKAGE_CHECKS.json').read_text());count=check_manifest(HERE)
    else:
        (HERE/'PACKAGE_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
        files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256')
        (HERE/'MANIFEST.sha256').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in files));count=len(files)
    print(json.dumps({'status':result['status'],'read_only':args.check,'manifest_entries':count},indent=2))

if __name__=='__main__':main()
