#!/usr/bin/env python3
"""Complete the NS intake package without modifying source/predecessor packages."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SOURCE=ROOT/'fin_nonlinear_stability_transfer_campaign_20261017'
PLAN=ROOT/'FIN_COHERENT_MODE_RADIATION_NEXT_CHAIN_20261017.md'
INTEGRATION=ROOT/'FIN_NONLINEAR_STABILITY_TRANSFER_INTEGRATION_20261017.md'
BASE='ffb62af2671f62aa582c688a90edd998e4a9621b29e35fd240e1863038c9e0d9'
REQUIRED=['INTAKE_20261017.md','REPORT_PL.md','PROOF_COMPLETIONS_AND_SCOPE.md',
          'PHASE_AND_RADIATION_SEED.md','CLAIM_REGISTER.json','NEXT_DECISION.md',
          'NONCONCLUSIONS.md','HANDOFF.md','REPLAY.md','verify_intake.py',
          'verification.json','COMPLETION_AUDIT_20261019.md','finalize_review.py']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def strict_json(p):
    def reject(v):raise ValueError(v)
    return json.loads(p.read_text(),parse_constant=reject)
def prose_links(text):
    # Equations in code blocks can contain E_[-R,R](0), which is not a link.
    prose=[];fence=None
    for line in text.splitlines():
        stripped=line.lstrip()
        marker='```' if stripped.startswith('```') else '~~~' if stripped.startswith('~~~') else None
        if marker:
            if fence is None:fence=marker
            elif fence==marker:fence=None
            continue
        if fence is None and not line.startswith(('    ','\t')):prose.append(line)
    clean=re.sub(r'`[^`]*`','', '\n'.join(prose))
    return re.findall(r'\[[^\]]+\]\(([^)]+)\)',clean)
def manifest(folder):
    names=[]
    for line in (folder/'MANIFEST.sha256').read_text().splitlines():
        h,n=line.split('  ',1);p=Path(n)
        assert not p.is_absolute() and '..' not in p.parts
        assert p.as_posix() not in names and sha(folder/p)==h,(folder,n)
        names.append(p.as_posix())
    actual={p.relative_to(folder).as_posix() for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    assert actual==set(names)|{'MANIFEST.sha256'},(folder,actual-set(names))
    return len(names)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');args=ap.parse_args()
    assert all((HERE/n).is_file() for n in REQUIRED)
    saved=strict_json(HERE/'verification.json')
    assert saved['status']=='PASS_SCOPED_NS_INTAKE'
    assert len(saved['checks'])==33 and all(c['passed'] for c in saved['checks'])
    assert saved['source_manifest']==saved['source_manifest_after']
    assert saved['source_manifest']['matched']==16 and manifest(SOURCE)==16
    assert sum(r['status']=='MATCHED_BYTES' for r in saved['original_references'])==4
    assert sum(r['status']=='UNAVAILABLE_ORIGINAL' for r in saved['original_references'])==1
    for r in saved['original_references']:
        if r['status']=='MATCHED_BYTES':assert sha(ROOT/r['local'])==r['sha256']
    assert saved['replay']['exit_code']==0 and saved['replay']['json_equal']
    assert saved['replay']['byte_equal'] and not saved['replay']['changed_controls']
    assert saved['verification_date']=='2026-10-19'
    assert re.fullmatch('[0-9a-f]{40}',saved['verification_HEAD'])
    for rel,h in saved['input_hashes'].items():assert sha(ROOT/rel)==h,rel
    for rel,record in saved['actual_dependencies'].items():
        assert manifest(ROOT/rel)==record['matched']
        assert sha(ROOT/rel/'MANIFEST.sha256')==record['manifest_sha256']
    assert sha(ROOT/'AGENTS.md')==BASE
    for p in HERE.glob('*.json'):strict_json(p)
    for p in HERE.glob('*.py'):ast.parse(p.read_text(),filename=str(p))
    for c in strict_json(HERE/'CLAIM_REGISTER.json')['claims']:assert (HERE/c['evidence']).is_file()
    assert prose_links('    E_[-R,R](0)-E_[-R,R](T)')==[]
    assert prose_links('[proof](PROOF_COMPLETIONS_AND_SCOPE.md)')==['PROOF_COMPLETIONS_AND_SCOPE.md']
    for p in [*HERE.glob('*.md'),PLAN,INTEGRATION]:
        assert all(l==l.rstrip() for l in p.read_text().splitlines()),p
        for target in prose_links(p.read_text()):
            if not target.startswith(('http:','https:','#')):assert (p.parent/target).exists(),(p,target)
    text=PLAN.read_text()
    for task in ('RP-001','RP-002','RP-003','RP-004','RP-005'):assert '## '+task+' ' in text
    assert 'RP-001–RP-005 NOT EXECUTED' in text
    index={'review_started':'2026-10-17','completion_date':'2026-10-19',
           'verification_HEAD':saved['verification_HEAD'],'source_manifest_sha256':sha(SOURCE/'MANIFEST.sha256'),
           'verification_sha256':sha(HERE/'verification.json'),
           'next_plan_sha256':sha(PLAN),'root_integration_sha256':sha(INTEGRATION),
           'actual_input_hashes':saved['input_hashes'],'AGENTS_sha256':BASE,
           'source_mapping':'bound quartic Phi; nonzero Phi6 comparison-only',
           'RP':'UNEXECUTED','physics_or_ToE':'NOT_ESTABLISHED'}
    p=HERE/'INTEGRATION_BINDING.json'
    if args.check:assert strict_json(p)==index
    else:p.write_text(json.dumps(index,indent=2,allow_nan=False)+'\n')
    files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256')
    lines=''.join(f'{sha(p)}  {p.name}\n' for p in files)
    if args.check:assert (HERE/'MANIFEST.sha256').read_text()==lines
    else:(HERE/'MANIFEST.sha256').write_text(lines)
    assert manifest(HERE)==len(files)
    print(json.dumps({'status':'PASS_SOURCE_PRESERVING_NS_INTEGRATION','read_only':args.check,
                      'review_entries':len(files),'source_entries':16,'intake_controls':33,
                      'verification_date':saved['verification_date'],'RP':'UNEXECUTED',
                      'AGENTS':'UNCHANGED'},indent=2))
if __name__=='__main__':main()
