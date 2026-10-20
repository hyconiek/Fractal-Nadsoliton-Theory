#!/usr/bin/env python3
"""Seal/check this new history intake; never rewrite old packages or AGENTS."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = ROOT/'FIN_POST_RP004_RESEARCH_ONLY_20261019'
ROOTDOCS = ['FIN_POST_RP004_HISTORY_INTEGRATION_20261020.md', 'FIN_RELATIONAL_HISTORY_NEXT_CHAIN_20261020.md']
REQUIRED = ['INTAKE_20261020.md', 'REPORT_PL.md', 'PROOF_AUDIT_AND_SCOPE.md',
            'MEAN_BRACKET_CORRECTION.md', 'CLAIM_REGISTER.json', 'NONCONCLUSIONS.md',
            'HANDOFF.md', 'PROOF_OR_OPEN_LEMMA.md', 'REPLAY.md', 'COMPLETION_AUDIT.md',
            'verify_intake.py', 'verification.json', 'finalize_review.py']


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def strict_json(p):
    def reject(x): raise ValueError(x)
    return json.loads(p.read_text(), parse_constant=reject)


def verify_manifest(folder, root_relative=False):
    result = []
    for line in (folder/'MANIFEST.sha256').read_text().splitlines():
        h, name = line.split(None, 1)
        name = name.lstrip('*')
        assert not Path(name).is_absolute() and '..' not in Path(name).parts and name not in result
        assert sha((ROOT if root_relative else folder)/name) == h, (folder, name)
        result.append(name)
    return result


def links(text):
    out, fence = [], None
    for line in text.splitlines():
        s = line.lstrip()
        marker = '```' if s.startswith('```') else '~~~' if s.startswith('~~~') else None
        if marker:
            fence = marker if fence is None else None if fence == marker else fence
        elif fence is None and not line.startswith(('    ', '\t')):
            out.append(line)
    return re.findall(r'\[[^\]]+\]\(([^)]+)\)', re.sub(r'`[^`]*`', '', '\n'.join(out)))


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--check', action='store_true'); args = parser.parse_args()
    assert all((HERE/n).is_file() for n in REQUIRED)
    d = strict_json(HERE/'verification.json')
    assert d['status'] == 'PASS_SCOPED_FIVE_STUDY_HISTORY_INTAKE'
    assert d['bounded_check_count'] == len(d['bounded_checks']) == 38
    assert len(set(d['bounded_checks'])) == 38
    assert d['consolidated_manifest']['entries'] == 75 and d['source_file_count'] == 77
    assert sum(r['entries'] for r in d['nested_manifests'].values()) == 55
    actual = {str(p.relative_to(SOURCE)):sha(p) for p in SOURCE.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    assert actual == d['source_preserved']
    assert sha(ROOT/'AGENTS.md') == d['live_AGENTS_sha256']
    assert len(d['copy_only_replay']) == 6 and all(r['returncode']==0 for r in d['copy_only_replay'].values())
    assert all(r.get('output_comparison', {}).get('max_numeric_drift', 0)<1e-10 for r in d['copy_only_replay'].values())
    rows = d['provenance']['input_reference_rows']
    assert len(rows)==26
    assert sum(r['state']=='MATCH' for r in rows)==23
    assert sum(r['state']=='MATCH_ARCHIVED_SNAPSHOT_NOT_LIVE' for r in rows)==1
    assert sum(r['state']=='MISSING_ORIGINAL' for r in rows)==2
    assert not d['RP004_full_proof_audited_here']
    assert d['F01']==d['F03']==d['F06']=='OPEN_SOURCE'
    architecture = verify_manifest(ROOT/'fin_strategy_audit_20261019', root_relative=True)
    assert len(architecture)==47
    for name in ['fin_rp004_flux_review_20261019', 'fin_intrinsic_entropy_source_review_20261019']:
        verify_manifest(ROOT/name)
    for p in HERE.glob('*.py'): ast.parse(p.read_text())
    for p in HERE.glob('*.json'): strict_json(p)
    claims = strict_json(HERE/'CLAIM_REGISTER.json')
    assert claims['full_anonymous_history_identification']=='PROPOSED_TARGET_NOT_EXECUTED'
    assert not claims['new_physical_model']
    assert all((HERE/c['evidence']).is_file() for c in claims['claims'])
    for p in [*HERE.glob('*.md'), *[ROOT/n for n in ROOTDOCS]]:
        assert all(line==line.rstrip() for line in p.read_text().splitlines()), p
        for target in links(p.read_text()):
            if not target.startswith(('http:','https:','#')):
                assert (p.parent/target).exists(), (p,target)
    extra = ['AGENTS.md','GEN001-gen956.md','FIN_POST_RP004_RESEARCH_DECISION_20261019.md',
             'FIN_POST_RP004_EXECUTION_STATUS_20261019.json',
             'FIN_POST_ID001_SOURCE_DECISION_20261019.md','FIN_RESEARCH_EXECUTION_STATUS_20261019.json',
             'fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py',
             'fin_rp004_flux_review_20261019/MANIFEST.sha256',
             'fin_intrinsic_entropy_source_review_20261019/MANIFEST.sha256',
             'fin_toe_bridge_research_20261015/verify_bridges.py',
             'fin_strategy_audit_20261019/INTRINSIC_EVOLUTION_ASSESSMENT.md']
    extra += [r['local_counterpart'] for r in rows if r['state']!='MISSING_ORIGINAL']
    names = sorted(set(architecture+extra+[str((SOURCE/p).relative_to(ROOT)) for p in actual]))
    text = ''.join(f'{sha(ROOT/n)}  {n}\n' for n in names)
    binding = {'date':'2026-10-20','entry_HEAD':'0c464f01','source_preserved_files':77,
               'consolidated_entries':75,'nested_direct_entries':55,'bounded_controls':38,
               'input_rows':26,'matching_input_rows_including_archived_strategy':24,'missing_original_input_rows':2,
               'bound_input_count':len(names),'live_AGENTS_sha256':d['live_AGENTS_sha256'],
               'verification_sha256':sha(HERE/'verification.json'),
               'root_integration_documents':{n:sha(ROOT/n) for n in ROOTDOCS},
               'source_open':True,'RP004_independent_audit_pending':True,'AH_tasks_executed':False,
               'AGENTS_or_old_plans_rewritten':False}
    if args.check:
        assert (HERE/'INPUTS.sha256').read_text()==text
        assert strict_json(HERE/'INTEGRATION_BINDING.json')==binding
    else:
        (HERE/'INPUTS.sha256').write_text(text)
        (HERE/'INTEGRATION_BINDING.json').write_text(json.dumps(binding,indent=2)+'\n')
    files = sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256')
    mtext = ''.join(f'{sha(p)}  {p.name}\n' for p in files)
    if args.check: assert (HERE/'MANIFEST.sha256').read_text()==mtext
    else: (HERE/'MANIFEST.sha256').write_text(mtext)
    assert len(verify_manifest(HERE))==len(files)
    print(json.dumps({'status':'PASS_SCOPED_HISTORY_INTEGRATION','read_only':args.check,
                      'review_entries':len(files),'source_files':77,'bound_inputs':len(names),
                      'bounded_controls':38,'source_law':'OPEN_SOURCE','AH_tasks_executed':False},indent=2))


if __name__=='__main__': main()
