#!/usr/bin/env python3
"""Seal/check this new ID intake only; no old source/strategy evidence is edited."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
SOURCE=ROOT/'fin_intrinsic_entropy_source_bound_20261019'
REQUIRED=['INTAKE_20261019.md','PROOF_AUDIT_AND_SCOPE.md','SUBSAMPLING_ERROR_FLOOR.md',
          'REPORT_PL.md','CLAIM_REGISTER.json','NONCONCLUSIONS.md','HANDOFF.md','REPLAY.md',
          'STRATEGY_RECONCILIATION.md','LITERATURE_SCOPE.md','COMPLETION_AUDIT.md',
          'verify_intake.py','verification.json','finalize_review.py']
ROOTDOCS=['FIN_INTRINSIC_ENTROPY_SOURCE_INTEGRATION_20261019.md',
          'FIN_POST_ID001_SOURCE_DECISION_20261019.md','FIN_RESEARCH_EXECUTION_STATUS_20261019.json']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def strict_json(p):
    def reject(v):raise ValueError(v)
    return json.loads(p.read_text(),parse_constant=reject)
def manifest(folder):
    names=[]
    for line in (folder/'MANIFEST.sha256').read_text().splitlines():
        h,n=line.split('  ',1);p=Path(n)
        assert not p.is_absolute() and '..' not in p.parts
        assert p.as_posix() not in names and sha(folder/p)==h,(folder,n);names.append(p.as_posix())
    actual={str(p.relative_to(folder)) for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    assert actual==set(names)|{'MANIFEST.sha256'},(folder,actual-set(names))
    return len(names)
def links(text):
    lines=[];fence=None
    for line in text.splitlines():
        t=line.lstrip();mark='```' if t.startswith('```') else '~~~' if t.startswith('~~~') else None
        if mark:
            if fence is None:fence=mark
            elif fence==mark:fence=None
            continue
        if fence is None and not line.startswith(('    ','\t')):lines.append(line)
    return re.findall(r'\[[^\]]+\]\(([^)]+)\)',re.sub(r'`[^`]*`','','\n'.join(lines)))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');a=ap.parse_args()
    assert all((HERE/n).is_file() for n in REQUIRED)
    d=strict_json(HERE/'verification.json')
    assert d['status']=='PASS_SCOPED_ID001_INTAKE_WITH_SUBSAMPLING_FLOOR'
    assert len(d['checks'])==23 and all(v['passed'] for v in d['checks'])
    assert d['source_manifest']==d['source_manifest_after'] and manifest(SOURCE)==14
    assert d['freshness']['returncode']==0 and d['replay']['returncode']==0
    assert d['replay']['max_numeric_drift']<1e-9
    for rel,h in d['input_hashes'].items():assert sha(ROOT/rel)==h,rel
    for line in (ROOT/'fin_strategy_audit_20261019/MANIFEST.sha256').read_text().splitlines():
        h,n=line.split('  ',1);assert sha(ROOT/n)==h,n
    assert d['architecture_manifest']['entries']==47
    ctx=(ROOT/'AGENTS.md').read_bytes();record=d['AGENTS_context']
    assert len(ctx)>=record['bytes'] and hashlib.sha256(ctx[:record['bytes']]).hexdigest()==record['sha256']
    assert hashlib.sha256(ctx[:record['older_verified_prefix_bytes']]).hexdigest()==record['older_verified_prefix_sha256']
    later=len(ctx)-record['bytes']
    for p in HERE.glob('*.json'):strict_json(p)
    for p in HERE.glob('*.py'):ast.parse(p.read_text())
    for claim in strict_json(HERE/'CLAIM_REGISTER.json')['claims']:assert (HERE/claim['evidence']).is_file()
    for p in [*HERE.glob('*.md'),*[ROOT/n for n in ROOTDOCS if n.endswith('.md')]]:
        assert all(line==line.rstrip() for line in p.read_text().splitlines()),p
        for t in links(p.read_text()):
            if not t.startswith(('http:','https:','#')):assert (p.parent/t).exists(),(p,t)
    status=strict_json(ROOT/ROOTDOCS[2]);tasks=status['tasks']
    assert tasks['ID-001']['execution']=='EXECUTED_AND_SCOPED_REVIEW_ACCEPTED'
    assert status['blocker_delta']['F01']=='OPEN_SOURCE' and status['default_automatic_start']==[]
    for name in ('RP004-A','RP004-B','RP004-C','MH-001'):assert tasks[name]['execution']=='NOT_EXECUTED'
    original=strict_json(ROOT/'research_tasks/TASKS.json')
    assert all(t['execution_status']=='NOT_EXECUTED' for t in original['tasks'])
    readpaths=[*ROOTDOCS,'FIN_RESEARCH_STRATEGY.md','FIN_RESEARCH_EXECUTION_PLAN.md',
       'FIN_BLOCKERS_REGISTER.md','FIN_PRIOR_RESEARCH_MAP.md',
       'fin_strategy_audit_20261019/AUDIT_SCOPE.md','fin_strategy_audit_20261019/INTRINSIC_EVOLUTION_ASSESSMENT.md',
       'fin_strategy_audit_20261019/ONTOLOGY_STEERING.md','fin_strategy_audit_20261019/REQUEST_ADDITIONAL.txt',
       'research_tasks/TASKS.json','research_tasks/RP004_SHARED_CONTRACT.json']
    readpaths += [str(p.relative_to(ROOT)) for p in sorted((ROOT/'research_tasks').glob('*.md'))]
    index={'date':'2026-10-19','actual_AGENTS_context':record,'verification_sha256':sha(HERE/'verification.json'),
           'source_manifest_sha256':sha(SOURCE/'MANIFEST.sha256'),
           'historical_architecture_manifest_sha256':sha(ROOT/'fin_strategy_audit_20261019/MANIFEST.sha256'),
           'strategy_and_status_binding':{n:sha(ROOT/n) for n in readpaths},
           'reading_scope':'named latest strategy/audit/task md and exact ID inputs; not another whole-repository proof audit',
           'F02':'SCOPED_CLASS_OBSTRUCTION_PROVED','F01':'OPEN_SOURCE','new_campaigns_activated':[]}
    p=HERE/'INTEGRATION_BINDING.json'
    if a.check:assert strict_json(p)==index
    else:
        assert later==0,'Do not seal an unreviewed later instruction revision'
        p.write_text(json.dumps(index,indent=2,allow_nan=False)+'\n')
    files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256')
    data=''.join(f'{sha(p)}  {p.name}\n' for p in files)
    if a.check:assert (HERE/'MANIFEST.sha256').read_text()==data
    else:(HERE/'MANIFEST.sha256').write_text(data)
    assert manifest(HERE)==len(files)
    print(json.dumps({'status':'PASS_SCOPED_ID001_INTEGRATION','read_only':a.check,
          'review_entries':len(files),'source_entries':14,'architecture_entries':47,'intake_controls':23,
          'F01':'OPEN_SOURCE','new_campaigns_activated':[],'later_instruction_bytes':later},indent=2))
if __name__=='__main__':main()
