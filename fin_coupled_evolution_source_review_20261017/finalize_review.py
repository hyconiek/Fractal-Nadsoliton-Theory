#!/usr/bin/env python3
"""Source-preserving packaging of the new UE intake; --check is read-only."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SOURCE=ROOT/'fin_coupled_evolution_source_campaign_20261017'
REQUIRED=['INTAKE_20261017.md','PROOF_AUDIT_AND_COMPLETIONS.md','LOCAL_SOURCE_BINDING.json',
          'REPORT_PL.md','CLAIM_REGISTER.json','NONCONCLUSIONS.md','NEXT_OPEN_ATOM.md',
          'HANDOFF.md','REPLAY.md','verify_intake.py','verification.json','finalize_review.py']
BASE='ffb62af2671f62aa582c688a90edd998e4a9621b29e35fd240e1863038c9e0d9'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def strict_load(p):
    def reject(x):raise ValueError(x)
    return json.loads(p.read_text(),parse_constant=reject)
def check_source():
    names=[]
    for line in (SOURCE/'MANIFEST.sha256').read_text().splitlines():
        h,n=line.split('  ',1);p=Path(n)
        assert not p.is_absolute() and '..' not in p.parts
        assert p.as_posix() not in names and sha(SOURCE/p)==h,n;names.append(p.as_posix())
    actual={p.relative_to(SOURCE).as_posix() for p in SOURCE.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    assert actual==set(names)|{'MANIFEST.sha256'}
    assert len(names)==22
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');a=ap.parse_args()
    assert all((HERE/n).is_file() for n in REQUIRED)
    result=strict_load(HERE/'verification.json')
    assert result['status']=='PASS_SCOPED_UE_INTAKE'
    assert len(result['checks'])==30 and all(c['passed'] for c in result['checks'])
    assert result['source_manifest']==result['source_manifest_after']
    assert result['source_manifest']['matched']==22 and len(result['original_local_references'])==3
    assert result['replay']['returncode']==0
    assert set(result['replay']['changed_fields'])<= {'runtime'}
    for p,h in result['input_hashes'].items():assert sha(ROOT/p)==h,p
    for p,h in result['original_local_references'].items():assert sha(SOURCE/p)==h,p
    assert sha(ROOT/'AGENTS.md')==BASE
    check_source()
    for p in HERE.glob('*.json'):strict_load(p)
    for p in HERE.glob('*.py'):ast.parse(p.read_text(),filename=str(p))
    for c in strict_load(HERE/'CLAIM_REGISTER.json')['claims']:assert (HERE/c['evidence']).is_file()
    root=ROOT/'FIN_COUPLED_EVOLUTION_SOURCE_INTEGRATION_20261017.md'
    for p in [*HERE.glob('*.md'),root]:
        assert all(l==l.rstrip() for l in p.read_text().splitlines()),p
        for t in re.findall(r'\[[^\]]+\]\(([^)]+)\)',p.read_text()):
            if not t.startswith(('http:','https:','#')):assert (p.parent/t).exists(),(p,t)
    index={'date':'2026-10-17','source_manifest_sha256':sha(SOURCE/'MANIFEST.sha256'),
           'root_integration_sha256':sha(root),'AGENTS_sha256':BASE,
           'actual_input_hashes':result['input_hashes'],
           'structural_binding_sha256':sha(HERE/'LOCAL_SOURCE_BINDING.json'),
           'scope':'local new review only; original artifacts and historical statements unchanged'}
    indexpath=HERE/'INTEGRATION_BINDING.json'
    if a.check:assert strict_load(indexpath)==index
    else:indexpath.write_text(json.dumps(index,indent=2)+'\n')
    files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256')
    lines=''.join(f'{sha(p)}  {p.name}\n' for p in files)
    if a.check:assert (HERE/'MANIFEST.sha256').read_text()==lines
    else:(HERE/'MANIFEST.sha256').write_text(lines)
    print(json.dumps({'status':'PASS_SOURCE_PRESERVING_UE_INTEGRATION','read_only':a.check,
                      'review_entries':len(files),'source_entries':22,'intake_controls':30,
                      'AGENTS':'UNCHANGED','strict_physical_source':'NOT_ESTABLISHED'},indent=2))
if __name__=='__main__':main()
