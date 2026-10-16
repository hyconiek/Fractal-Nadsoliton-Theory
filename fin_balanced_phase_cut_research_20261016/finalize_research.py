#!/usr/bin/env python3
"""Mechanical packaging of new research, not a mathematical proof checker."""
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
    required=['ACTUAL_COUNTERCUT_TO_STRONG_SHIELDING.md','BALANCED_PHASE_DEFECT_THEOREM.md',
        'COMPACT_PHASE_PREFIX_LOCAL_LIMIT.md','EXACT_COMPACT_CURRENT_MARGINALS_AND_GLUING.md',
        'CURRENT_CARRY_CONTRACT.json','REPORT_PL.md','NONCONCLUSIONS.md','REPLAY.md','HANDOFF.md',
        'verify_research.py','RESULTS.json','finalize_research.py']
    assert all((HERE/n).is_file() for n in required)
    for path in HERE.glob('*.py'):ast.parse(path.read_text(),filename=str(path))
    for path in HERE.glob('*.json'):json.loads(path.read_text(),parse_constant=reject)
    for path in HERE.glob('*.md'):
        assert all(x==x.rstrip() for x in path.read_text().splitlines()),path
        for target in re.findall(r'\[[^\]]+\]\(([^)]+)\)',path.read_text()):
            if not target.startswith(('http:','https:','#')):assert (path.parent/target).exists(),(path,target)
    data=json.loads((HERE/'RESULTS.json').read_text())
    assert data['status']=='PASS_BOUNDED_BALANCED_PHASE_RESEARCH_CONTROLS'
    assert len(data['checks'])==20 and all(c['passed'] for c in data['checks'])
    for rel,h in data['input_hashes'].items():assert sha(ROOT/rel)==h,rel
    result={'date':'2026-10-16','status':'PASS_MECHANICAL_BALANCED_PHASE_RESEARCH_PACKAGE',
        'required_files':required,'bounded_controls':20,'all_N_proof_by_tests':False,
        'strong_shielding':'REFUTED_ON_ACTUAL_CUT_SEQUENCE','original_SC':'OPEN',
        'current_marginals':'DERIVED_ON_FIXED_COMPACT_BANDS','edge_coupling':'OPEN','full_generator':'NOT_CONSTRUCTED'}
    if args.check:
        assert result==json.loads((HERE/'PACKAGE_CHECKS.json').read_text())
        names=[]
        for line in (HERE/'MANIFEST.sha256').read_text().splitlines():
            h,name=line.split('  ',1);assert name not in names and sha(HERE/name)==h;names.append(name)
        assert set(names)=={p.name for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256'}
    else:
        (HERE/'PACKAGE_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
        files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256')
        (HERE/'MANIFEST.sha256').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in files));names=[p.name for p in files]
    print(json.dumps({'status':result['status'],'read_only':args.check,'manifest_entries':len(names)},indent=2))

if __name__=='__main__':main()
