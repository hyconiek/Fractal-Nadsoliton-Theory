#!/usr/bin/env python3
"""Complete new DH review without editing incoming science or AGENTS."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');a=ap.parse_args()
    required=['INTAKE_20261017.md','PROOF_COMPLETIONS_AND_SCOPE.md','CLUSTERED_COMPLEMENT_SEED.md',
              'CLAIM_REGISTER.json','REPORT.md','HANDOFF.md','NONCONCLUSIONS.md','REPLAY.md',
              'verify_intake.py','verification.json','finalize_review.py']
    assert all((HERE/n).is_file() for n in required)
    d=json.loads((HERE/'verification.json').read_text())
    assert d['status']=='PASS_SCOPED_DH_DENSE_COMPLEMENT_INTAKE'
    assert all(c['passed'] for c in d['checks'])
    assert d['source_manifest']==d['source_manifest_after'] and d['source_manifest']['matched']==29
    assert d['replay']['logs']['verify_dh.py']['exit_code']==0
    for rel,h in d['input_hashes'].items():assert sha(ROOT/rel)==h,rel
    source=ROOT/'fin_dense_orbit_fibre_transport_campaign_20261017'
    for l in (source/'MANIFEST.sha256').read_text().splitlines():
        h,n=l.split('  ',1);assert sha(source/n)==h,n
    for p in HERE.glob('*.py'):ast.parse(p.read_text())
    for p in HERE.glob('*.json'):json.loads(p.read_text())
    for p in HERE.glob('*.md'):
        for t in re.findall(r'\[[^\]]+\]\(([^)]+)\)',p.read_text()):
            if not t.startswith(('http:','https:','#')):assert (p.parent/t).exists(),(p,t)
    files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256')
    lines=''.join(f'{sha(p)}  {p.name}\n' for p in files)
    if a.check:assert (HERE/'MANIFEST.sha256').read_text()==lines
    else:(HERE/'MANIFEST.sha256').write_text(lines)
    print(json.dumps({'status':'PASS_SCOPED_DH_INTAKE_PACKAGE','entries':len(files),'bounded_checks':len(d['checks']),'read_only':a.check}))
if __name__=='__main__':main()
