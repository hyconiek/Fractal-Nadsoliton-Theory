#!/usr/bin/env python3
"""Package the scoped RP intake, preserving all prior evidence and contexts."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SOURCE=ROOT/'fin_coherent_mode_radiation_campaign_20261017'
PLAN=ROOT/'FIN_RP004_RETARDED_FLUX_CONTINUATION_20261019.md'
STATUS=ROOT/'FIN_COHERENT_MODE_RADIATION_INTEGRATION_20261019.md'
REQUIRED=['INTAKE_20261019.md','REPORT_PL.md','PROOF_COMPLETIONS_AND_SCOPE.md',
          'RETARDED_FLUX_PROOF_ROUTE.md','CLAIM_REGISTER.json','NEXT_DECISION.md',
          'NONCONCLUSIONS.md','HANDOFF.md','REPLAY.md','COMPLETION_AUDIT.md',
          'verify_intake.py','verification.json','finalize_review.py']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def strict_json(p):
    def reject(x):raise ValueError(x)
    return json.loads(p.read_text(),parse_constant=reject)
def manifest(folder):
    names=[]
    for row in (folder/'MANIFEST.sha256').read_text().splitlines():
        h,n=row.split('  ',1);p=Path(n)
        assert not p.is_absolute() and '..' not in p.parts
        assert p.as_posix() not in names and sha(folder/p)==h,(folder,n)
        names.append(p.as_posix())
    actual={p.relative_to(folder).as_posix() for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    assert actual==set(names)|{'MANIFEST.sha256'},(folder,actual-set(names))
    return len(names)
def prose_links(text):
    lines=[];fence=None
    for line in text.splitlines():
        stripped=line.lstrip()
        marker='```' if stripped.startswith('```') else '~~~' if stripped.startswith('~~~') else None
        if marker:
            if fence is None:fence=marker
            elif fence==marker:fence=None
            continue
        if fence is None and not line.startswith(('    ','\t')):lines.append(line)
    text=re.sub(r'`[^`]*`','', '\n'.join(lines))
    return re.findall(r'\[[^\]]+\]\(([^)]+)\)',text)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');a=ap.parse_args()
    assert all((HERE/n).is_file() for n in REQUIRED)
    saved=strict_json(HERE/'verification.json')
    assert saved['status']=='PASS_SCOPED_RP_ON_SHELL_INTAKE'
    assert len(saved['checks'])==34 and all(c['passed'] for c in saved['checks'])
    assert saved['source_manifest']==saved['source_manifest_after']
    assert manifest(SOURCE)==saved['source_manifest']['entries']==19
    assert len(saved['matched_Library_payloads'])==4
    for r in saved['matched_Library_payloads']:assert sha(ROOT/r['local'])==r['sha256']
    for rel,h in saved['input_hashes'].items():assert sha(ROOT/rel)==h,rel
    for rel,d in saved['dependencies'].items():
        assert manifest(ROOT/rel)==d['entries']
        assert sha(ROOT/rel/'MANIFEST.sha256')==d['manifest_sha256']
    assert saved['replay']['returncode']==0 and 'PASS 16/16 bounded controls' in saved['replay']['stdout']
    assert saved['uniform_Gamma_lower_rational']=='125/73728'
    assert saved['RP004']=='OPEN_NONLINEAR_OUTWARD_FLUX'
    agents=(ROOT/'AGENTS.md').read_bytes();context=saved['AGENTS_context']
    assert len(agents)>=context['bytes']
    assert hashlib.sha256(agents[:context['bytes']]).hexdigest()==context['sha256']
    prior=context['prior_verified_prefix_bytes']
    assert hashlib.sha256(agents[:prior]).hexdigest()==context['prior_verified_prefix_sha256']
    assert '## FIN RP intake: on-shell coupling admitted, nonlinear flux still open' in agents[:context['bytes']].decode()
    later=len(agents)-context['bytes']
    for p in HERE.glob('*.json'):strict_json(p)
    for p in HERE.glob('*.py'):ast.parse(p.read_text(),filename=str(p))
    for c in strict_json(HERE/'CLAIM_REGISTER.json')['claims']:assert (HERE/c['evidence']).is_file()
    for p in [*HERE.glob('*.md'),PLAN,STATUS]:
        assert all(l==l.rstrip() for l in p.read_text().splitlines()),p
        for target in prose_links(p.read_text()):
            if not target.startswith(('http:','https:','#')):assert (p.parent/target).exists(),(p,target)
    for task in ('RP004-A','RP004-B','RP004-C','RP004-D'):assert '## '+task+' ' in PLAN.read_text()
    assert 'UNEXECUTED CONTINUATION ORDER' in PLAN.read_text()
    index={'date':'2026-10-19','source_manifest_sha256':sha(SOURCE/'MANIFEST.sha256'),
           'verification_sha256':sha(HERE/'verification.json'),'actual_instruction_context':context,
           'next_continuation_sha256':sha(PLAN),'root_integration_sha256':sha(STATUS),
           'scientific_inputs':saved['input_hashes'],'verdict':'ON_SHELL_ONLY; NONLINEAR_FLUX_OPEN',
           'next_tasks':'NOT_EXECUTED','natural_kernel_or_physics_source':'NOT_ESTABLISHED'}
    p=HERE/'INTEGRATION_BINDING.json'
    if a.check:assert strict_json(p)==index
    else:
        assert later==0,'Do not seal an unreviewed later instruction context'
        p.write_text(json.dumps(index,indent=2,allow_nan=False)+'\n')
    files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256')
    lines=''.join(f'{sha(p)}  {p.name}\n' for p in files)
    if a.check:assert (HERE/'MANIFEST.sha256').read_text()==lines
    else:(HERE/'MANIFEST.sha256').write_text(lines)
    assert manifest(HERE)==len(files)
    print(json.dumps({'status':'PASS_SCOPED_RP_INTEGRATION','read_only':a.check,
                      'review_entries':len(files),'source_entries':19,'bounded_checks':34,
                      'instruction_later_tail_bytes':later,'RP004':'OPEN_NONLINEAR_OUTWARD_FLUX',
                      'continuation':'NOT_EXECUTED'},indent=2))
if __name__=='__main__':main()
