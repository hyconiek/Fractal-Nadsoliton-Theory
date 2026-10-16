#!/usr/bin/env python3
"""Mechanical physical-search package checks, not a source-of-physics proof."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def reject(value):raise ValueError(value)
def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args()
    required=['Q12_SPECTATOR_ACTIVITY_OBSTRUCTION.md','THERMAL_RESET_RESOURCE_BOUNDARY.md',
        'LITERATURE_DATABASE.md','LITERATURE_REGISTER.json','BROWSER_READ_LOG.jsonl','browser_capture.py',
        'PHYSICS_GAPS_AND_SEARCH_PL.md','TOOL_ACCESS_STATUS.json','NONCONCLUSIONS.md','REPLAY.md',
        'HANDOFF.md','verify_bridges.py','RESULTS.json','finalize_research.py']
    assert all((HERE/n).is_file() for n in required)
    for path in HERE.glob('*.py'):ast.parse(path.read_text(),filename=str(path))
    for path in HERE.glob('*.json'):json.loads(path.read_text(),parse_constant=reject)
    for path in HERE.glob('*.md'):
        assert all(line==line.rstrip() for line in path.read_text().splitlines()),path
        for target in re.findall(r'\[[^\]]+\]\(([^)]+)\)',path.read_text()):
            if not target.startswith(('https:','http:','#')):assert (path.parent/target).exists(),(path,target)
    saved=json.loads((HERE/'RESULTS.json').read_text());assert len(saved['checks'])==16 and all(c['passed'] for c in saved['checks'])
    assert saved['status']=='PASS_BOUNDED_Q12_PHYSICAL_COMPATIBILITY_CHECKS' and saved['empirical_evidence'] is False
    for rel,h in saved['input_hashes'].items():assert sha(ROOT/rel)==h,rel
    log=[json.loads(line) for line in (HERE/'BROWSER_READ_LOG.jsonl').read_text().splitlines()]
    assert len(log)==3 and len({row['label'] for row in log})==3
    papers=json.loads((HERE/'LITERATURE_REGISTER.json').read_text())['works'];assert len(papers)==5
    result={'date':'2026-10-17','status':'PASS_MECHANICAL_Q12_PHYSICAL_COMPATIBILITY_PACKAGE',
        'required_files':required,'bounded_controls':16,'actual_browser_captures':3,'annotated_primary_works':5,
        'PDF_hash_scope':'no original PDF byte hashes claimed','physical_source':'NOT_FOUND',
        'conditional_implementation':'ANALYTICAL_TARGET_DESIGNED_NOT_AUTONOMOUS_SOURCE',
        'scope':'JSON/AST/link/hash/capture provenance, not physical validation'}
    if args.check:
        assert result==json.loads((HERE/'PACKAGE_CHECKS.json').read_text());names=[]
        for row in (HERE/'MANIFEST.sha256').read_text().splitlines():
            h,n=row.split('  ',1);assert n not in names and sha(HERE/n)==h,n;names.append(n)
        assert set(names)=={p.name for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256'}
    else:
        (HERE/'PACKAGE_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
        files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256')
        (HERE/'MANIFEST.sha256').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in files));names=[p.name for p in files]
    print(json.dumps({'status':result['status'],'read_only':args.check,'manifest_entries':len(names)},indent=2))

if __name__=='__main__':main()
