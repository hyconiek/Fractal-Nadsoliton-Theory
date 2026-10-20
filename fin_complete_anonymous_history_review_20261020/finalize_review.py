#!/usr/bin/env python3
"""Seal/check only this new AH intake, its rational certificate and two root indices."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re
import verify_intake as audit

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = ROOT/'fin_complete_anonymous_history_campaign_20261020'
ROOTDOCS = ['FIN_COMPLETE_ANONYMOUS_HISTORY_INTEGRATION_20261020.md',
            'FIN_POST_AH_SOURCE_AND_VALIDATION_NEXT_CHAIN_20261020.md']
REQUIRED = ['INTAKE_20261020.md','PROOF_COMPLETIONS_AND_SCOPE.md','REPORT_PL.md',
            'CLAIM_REGISTER.json','NONCONCLUSIONS.md','PROOF_OR_OPEN_LEMMA.md',
            'HANDOFF.md','REPLAY.md','COMPLETION_AUDIT.md','verify_intake.py',
            'verification.json','finalize_review.py']


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def load(p):
    def reject(v): raise ValueError(v)
    def pairs(items):
        d={}
        for k,v in items:
            assert k not in d,('duplicate key',k)
            d[k]=v
        return d
    return json.loads(p.read_text(),parse_constant=reject,object_pairs_hook=pairs)


def manifest(folder,root_relative=False):
    out=[]
    for line in (folder/'MANIFEST.sha256').read_text().splitlines():
        h,n=line.split(None,1)
        assert not Path(n).is_absolute() and '..' not in Path(n).parts and n not in out
        assert sha((ROOT if root_relative else folder)/n)==h,(folder,n)
        out.append(n)
    return out


def links(text):
    lines=[]; fence=None
    for line in text.splitlines():
        s=line.lstrip()
        marker='```' if s.startswith('```') else '~~~' if s.startswith('~~~') else None
        if marker:
            fence=marker if fence is None else None if fence==marker else fence
        elif fence is None and not line.startswith(('    ','\t')):
            lines.append(line)
    return re.findall(r'\[[^\]]+\]\(([^)]+)\)',re.sub(r'`[^`]*`','','\n'.join(lines)))


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--check',action='store_true'); args=parser.parse_args()
    assert all((HERE/n).is_file() for n in REQUIRED)
    d=load(HERE/'verification.json')
    assert d['status']=='PASS_SCOPED_AH_INTAKE_WITH_EXACT_ARITHMETIC_SUPPORT'
    assert d['bounded_check_count']==len(d['bounded_checks'])==len(set(d['bounded_checks']))==49
    assert d['source_manifest']=={'entries':15,'complete_coverage':True}
    assert len(d['source_preserved'])==d['source_file_count']==16
    assert audit.snapshot(SOURCE)==d['source_preserved']
    assert sha(ROOT/'AGENTS.md')==d['live_AGENTS_sha256']
    assert d['copy_only_replay']['verify_package.py']['returncode']==0
    assert d['copy_only_replay']['verify.py']['returncode']==0
    assert d['copy_only_replay']['result_comparison']['JSON_equal']
    assert sum(r['state']=='MATCH' for r in d['original_inputs'])==6 and len(d['original_inputs'])==7
    assert sum(r['state']=='MISSING_ORIGINAL' for r in d['original_inputs'])==1
    assert not d['RP004_audit_executed'] and d['F01']==d['F03']==d['F06']=='OPEN_SOURCE'
    cert=d['example_budget_certificate']; audit.check_certificate(cert)
    cp=HERE/'FINITE_PREFIX_BUDGET_CERTIFICATE.json'
    if args.check: assert load(cp)==cert
    else: cp.write_text(json.dumps(cert,indent=2)+'\n')
    assert cert['m']==569231067 and cert['L']==1138462136
    architecture=manifest(ROOT/'fin_strategy_audit_20261019',root_relative=True)
    assert len(architecture)==47
    for name in ['fin_post_rp004_history_review_20261020', 'fin_rp004_flux_review_20261019',
                 'fin_intrinsic_entropy_source_review_20261019']:
        manifest(ROOT/name)
    assert len(manifest(SOURCE))==15
    claims=load(HERE/'CLAIM_REGISTER.json')
    assert claims['next_automatic_research']==[] and not claims['new_physical_model']
    assert all((HERE/r['evidence']).is_file() for r in claims['claims'])
    for p in HERE.glob('*.json'): load(p)
    for p in HERE.glob('*.py'): ast.parse(p.read_text())
    for p in [*HERE.glob('*.md'),*[ROOT/n for n in ROOTDOCS]]:
        assert all(line==line.rstrip() for line in p.read_text().splitlines()),p
        for target in links(p.read_text()):
            if not target.startswith(('http:','https:','#')): assert (p.parent/target).exists(),(p,target)
    extra=['AGENTS.md','GEN001-gen956.md','FIN_RELATIONAL_HISTORY_NEXT_CHAIN_20261020.md',
           'FIN_POST_RP004_RESEARCH_DECISION_20261019.md','FIN_POST_RP004_EXECUTION_STATUS_20261019.json',
           'FIN_POST_RP004_HISTORY_INTEGRATION_20261020.md',
           'FIN_POST_ID001_SOURCE_DECISION_20261019.md',
           'fin_post_rp004_history_review_20261020/INTAKE_20261020.md',
           'fin_post_rp004_history_review_20261020/PROOF_AUDIT_AND_SCOPE.md',
           'fin_post_rp004_history_review_20261020/MEAN_BRACKET_CORRECTION.md',
           'fin_post_rp004_history_review_20261020/MANIFEST.sha256',
           'fin_rp004_flux_review_20261019/MANIFEST.sha256',
           'fin_intrinsic_entropy_source_review_20261019/MANIFEST.sha256',
           'research_tasks/00_EXECUTION_RULES.md']
    extra += [r['current_counterpart'] for r in d['original_inputs'] if r['state']=='MATCH']
    names=sorted(set(architecture+extra+[str((SOURCE/n).relative_to(ROOT)) for n in d['source_preserved']]))
    input_text=''.join(f'{sha(ROOT/n)}  {n}\n' for n in names)
    binding={'date':'2026-10-20','entry_HEAD':'0c464f01','source_entries':15,'source_present_files':16,
             'original_input_matches':6,'original_input_missing':1,'bounded_controls':49,
             'bound_input_count':len(names),'live_AGENTS_sha256':d['live_AGENTS_sha256'],
             'verification_sha256':sha(HERE/'verification.json'),
             'certificate_sha256':sha(cp),'root_documents':{n:sha(ROOT/n) for n in ROOTDOCS},
             'AH_verdict':claims['verdict'],'F01_F03_F06':'OPEN_SOURCE',
             'RP004_audit_executed':False,'source_or_AGENTS_rewritten':False}
    if args.check:
        assert (HERE/'INPUTS.sha256').read_text()==input_text
        assert load(HERE/'INTEGRATION_BINDING.json')==binding
    else:
        (HERE/'INPUTS.sha256').write_text(input_text)
        (HERE/'INTEGRATION_BINDING.json').write_text(json.dumps(binding,indent=2)+'\n')
    files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256')
    mtext=''.join(f'{sha(p)}  {p.name}\n' for p in files)
    if args.check: assert (HERE/'MANIFEST.sha256').read_text()==mtext
    else: (HERE/'MANIFEST.sha256').write_text(mtext)
    assert len(manifest(HERE))==len(files)
    print(json.dumps({'status':'PASS_SCOPED_AH_INTEGRATION_NOT_A_FORMAL_PROOF',
                      'read_only':args.check,'review_entries':len(files),'source_entries':15,
                      'bounded_controls':49,'input_bindings':len(names),
                      'example_m':cert['m'],'example_L':cert['L'],
                      'source_law':'OPEN_SOURCE','RP004_audit_executed':False},indent=2))


if __name__=='__main__': main()
