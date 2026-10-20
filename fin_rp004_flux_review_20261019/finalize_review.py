#!/usr/bin/env python3
"""Seal/check only this new proof submission and its two root status documents."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = ROOT/'fin_rp004_flux_continuation_20261019'
ROOTDOCS = ['FIN_POST_RP004_RESEARCH_DECISION_20261019.md',
            'FIN_POST_RP004_EXECUTION_STATUS_20261019.json']
REQUIRED = ['INTAKE_20261019.md', 'MOVING_CONE_THEOREM.md', 'BOOTSTRAP_COMPLETION.md',
            'OUTWARD_FLUX_THEOREM.md', 'STRATEGY_RECONCILIATION.md', 'REPORT_PL.md',
            'PROOF_OR_OPEN_LEMMA.md', 'CLAIM_REGISTER.json', 'COMPLETE_FLUX_BUDGET.json',
            'PREPARATION_AND_ERROR_CONTRACT.json', 'COMPLETION_AUDIT.md',
            'NONCONCLUSIONS.md', 'HANDOFF.md', 'REPLAY.md', 'verify_intake.py',
            'verification.json', 'finalize_review.py']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pairs(items):
    d = {}
    for k, v in items:
        assert k not in d, ('duplicate JSON key', k)
        d[k] = v
    return d


def load(path):
    def reject(value):
        raise ValueError(value)
    return json.loads(path.read_text(), parse_constant=reject, object_pairs_hook=pairs)


def verify_manifest(folder, root_relative=False):
    declared = []
    for line in (folder/'MANIFEST.sha256').read_text().splitlines():
        h, name = line.split(None, 1)
        name = name.lstrip('*')
        p = Path(name)
        assert not p.is_absolute() and '..' not in p.parts and name not in declared
        assert sha((ROOT if root_relative else folder)/p) == h, (folder, name)
        declared.append(name)
    return declared


def local_links(text):
    out = []
    fence = None
    for line in text.splitlines():
        stripped = line.lstrip()
        mark = '```' if stripped.startswith('```') else '~~~' if stripped.startswith('~~~') else None
        if mark:
            fence = mark if fence is None else None if fence == mark else fence
        elif fence is None and not line.startswith(('    ', '\t')):
            out.append(line)
    return re.findall(r'\[[^\]]+\]\(([^)]+)\)', re.sub(r'`[^`]*`', '', '\n'.join(out)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    assert all((HERE/p).is_file() for p in REQUIRED)
    data = load(HERE/'verification.json')
    assert data['status'] == 'PASS_BOUNDED_CONTROLS_WITH_DISCLOSED_SOURCE_OMISSION'
    assert data['intake_bounded_check_count'] == len(data['intake_bounded_checks']) == 26
    assert len(data['copy_only_source_replay']) == 2
    assert all(r['returncode'] == 0 for r in data['copy_only_source_replay'])
    assert data['producer_symbolic_assert_statements'] == 12
    assert data['strategy_freshness']['returncode'] == 0
    assert data['analytic_proofs_not_established_by_flags']
    assert not data['independent_peer_review'] and not data['physical_validation']
    for folder, key in [(SOURCE, 'source_preserved'),
                        (ROOT/'fin_strategy_audit_20261019', 'strategy_preserved'),
                        (ROOT/'research_tasks', 'tasks_preserved')]:
        actual = {str(p.relative_to(folder)): sha(p) for p in folder.rglob('*')
                  if p.is_file() and '__pycache__' not in p.parts}
        assert actual == data[key], folder
    assert sha(ROOT/'AGENTS.md') == data['live_AGENTS_sha256']
    assert not (SOURCE/'finalize.py').exists(), 'New original helper needs a fresh provenance audit'
    source_rows = data['manifest_checks']
    assert sum(r['matching'] for r in source_rows.values()) == 30
    assert sum(r['entries'] for r in source_rows.values()) == 31
    assert source_rows['.']['missing'] == ['finalize.py']
    inputs = data['original_inputs']['.']
    assert len(inputs) == 10 and sum(r['state'] == 'MATCH' for r in inputs) == 8
    architecture = verify_manifest(ROOT/'fin_strategy_audit_20261019', root_relative=True)
    assert len(architecture) == 47
    for name in ['fin_coherent_mode_radiation_review_20261019',
                 'fin_coherent_mode_radiation_campaign_20261017',
                 'fin_nonlinear_stability_transfer_review_20261017',
                 'fin_intrinsic_entropy_source_review_20261019']:
        verify_manifest(ROOT/name)
    for p in HERE.glob('*.json'):
        load(p)
    for p in HERE.glob('*.py'):
        ast.parse(p.read_text())
    claims = load(HERE/'CLAIM_REGISTER.json')
    assert not claims['independent_reviewer_acceptance']
    assert all((HERE/row['evidence']).is_file() for row in claims['claims'])
    for p in [*HERE.glob('*.md'), ROOT/ROOTDOCS[0]]:
        assert all(line == line.rstrip() for line in p.read_text().splitlines()), p
        for link in local_links(p.read_text()):
            if not link.startswith(('http:', 'https:', '#')):
                assert (p.parent/link).exists(), (p, link)
    status = load(ROOT/ROOTDOCS[1])
    assert not status['AGENTS_rewritten'] and not status['historical_plans_and_statuses_rewritten']
    assert status['default_automatic_new_research'] == []
    assert status['blockers']['F01'] == status['blockers']['F03'] == status['blockers']['F06'] == 'OPEN_SOURCE'
    old = load(ROOT/'FIN_RESEARCH_EXECUTION_STATUS_20261019.json')
    assert old['tasks']['RP004-A']['execution'] == 'NOT_EXECUTED'
    contract = load(HERE/'PREPARATION_AND_ERROR_CONTRACT.json')
    assert sha(ROOT/contract['parent_contract']) == contract['parent_contract_sha256']
    extra = ['AGENTS.md', 'FIN_RESEARCH_EXECUTION_STATUS_20261019.json',
             'FIN_POST_ID001_SOURCE_DECISION_20261019.md',
             'fin_intrinsic_entropy_source_review_20261019/STRATEGY_RECONCILIATION.md',
             'fin_intrinsic_entropy_source_review_20261019/SUBSAMPLING_ERROR_FLOOR.md',
             'fin_nonlinear_stability_transfer_review_20261017/PROOF_COMPLETIONS_AND_SCOPE.md',
             'fin_nonlinear_stability_transfer_review_20261017/INTAKE_20261017.md',
             'fin_coherent_mode_radiation_review_20261019/PROOF_COMPLETIONS_AND_SCOPE.md',
             'fin_coherent_mode_radiation_review_20261019/RETARDED_FLUX_PROOF_ROUTE.md',
             'fin_coherent_mode_radiation_review_20261019/INTAKE_20261019.md',
             'FIN_RP004_RETARDED_FLUX_CONTINUATION_20261019.md']
    sources = sorted(set(architecture + extra +
                         [str((SOURCE/p).relative_to(ROOT)) for p in data['source_preserved']]))
    intext = ''.join(f'{sha(ROOT/p)}  {p}\n' for p in sources)
    binding = {'date': '2026-10-19', 'entry_HEAD': '0c464f01',
               'review_status': 'CONDITIONAL_ANALYTICAL_PROOFS_SUPPLIED_SELF_AUDITED',
               'independent_review_done': False, 'AGENTS_sha256': data['live_AGENTS_sha256'],
               'bound_input_count': len(sources),
               'source_present_file_count': len(data['source_preserved']),
               'architecture_manifest_entries': len(architecture),
               'intake_controls': 26, 'verification_sha256': sha(HERE/'verification.json'),
               'root_status_documents': {p: sha(ROOT/p) for p in ROOTDOCS},
               'original_inputs_eight_match_two_different_current_counterparts': True,
               'source_root_finalizer_absence_preserved': True,
               'source_model_and_physical_gates_unchanged': True}
    if args.check:
        assert (HERE/'INPUTS.sha256').read_text() == intext
        assert load(HERE/'INTEGRATION_BINDING.json') == binding
    else:
        (HERE/'INPUTS.sha256').write_text(intext)
        (HERE/'INTEGRATION_BINDING.json').write_text(json.dumps(binding, indent=2)+'\n')
    files = sorted(p for p in HERE.iterdir() if p.is_file() and p.name != 'MANIFEST.sha256')
    manifest_text = ''.join(f'{sha(p)}  {p.name}\n' for p in files)
    if args.check:
        assert (HERE/'MANIFEST.sha256').read_text() == manifest_text
    else:
        (HERE/'MANIFEST.sha256').write_text(manifest_text)
    assert len(verify_manifest(HERE)) == len(files)
    print(json.dumps({'status': 'PASS_NEW_REVIEW_PACKAGE_NOT_A_PDE_PROOF_CERTIFICATE',
                      'read_only': args.check, 'review_entries': len(files),
                      'inputs': len(sources), 'source_matches': 30,
                      'source_absent': ['finalize.py'], 'bounded_controls': 26,
                      'AGENTS_and_old_strategy_preserved': True,
                      'independent_review_done': False}, indent=2))


if __name__ == '__main__':
    main()
