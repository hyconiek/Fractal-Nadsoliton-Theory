#!/usr/bin/env python3
"""Seal/check only this new ES intake/continuation and its two root documents."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SOURCE=ROOT/'fin_evolution_source_research_20261020'
ROOTDOCS=['FIN_EVOLUTION_SOURCE_INTEGRATION_20261020.md','FIN_CAUSAL_SOURCE_TRANSFER_NEXT_CHAIN_20261020.md']
REQUIRED=['INTAKE_20261020.md','EXACT_CAUSAL_FILTER.md','FULL_PATH_TRANSFER.md',
          'FILTERED_BLOCK_DEFECT.md','REPORT_PL.md','CLAIM_REGISTER.json',
          'PROOF_OR_OPEN_LEMMA.md','NONCONCLUSIONS.md','HANDOFF.md','REPLAY.md',
          'COMPLETION_AUDIT.md','verify_intake.py','verification.json','finalize_review.py']


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def load(p):
    def reject(x): raise ValueError(x)
    return json.loads(p.read_text(),parse_constant=reject)


def manifest(folder,root_relative=False):
    names=[]
    for line in (folder/'MANIFEST.sha256').read_text().splitlines():
        h,n=line.split(None,1)
        assert not Path(n).is_absolute() and '..' not in Path(n).parts and n not in names
        assert sha((ROOT if root_relative else folder)/n)==h,(folder,n)
        names.append(n)
    return names


def links(text):
    out=[]; fence=None
    for line in text.splitlines():
        s=line.lstrip(); mark='```' if s.startswith('```') else '~~~' if s.startswith('~~~') else None
        if mark: fence=mark if fence is None else None if fence==mark else fence
        elif fence is None and not line.startswith(('    ','\t')): out.append(line)
    return re.findall(r'\[[^\]]+\]\(([^)]+)\)',re.sub(r'`[^`]*`','','\n'.join(out)))


def main():
    p=argparse.ArgumentParser(); p.add_argument('--check',action='store_true'); a=p.parse_args()
    assert all((HERE/n).is_file() for n in REQUIRED)
    d=load(HERE/'verification.json')
    assert d['status']=='PASS_SCOPED_ES_INTAKE_AND_NEW_CAUSAL_TRANSFER_CONTROLS'
    assert d['bounded_control_count']==len(d['bounded_controls'])==len(set(d['bounded_controls']))==83
    assert len(d['source_preserved'])==d['source_file_count']==16
    assert d['source_manifest_entries']==15 and d['source_input_matches']==14
    actual={str(f.relative_to(SOURCE)):sha(f) for f in SOURCE.rglob('*') if f.is_file() and '__pycache__' not in f.parts}
    assert actual==d['source_preserved']
    assert sha(ROOT/'AGENTS.md')==d['live_AGENTS_sha256']
    assert d['copy_only_replay']['check_package.py --check']['returncode']==0
    assert d['copy_only_replay']['verify_research.py']['returncode']==0
    assert d['copy_only_replay']['result_comparison']['JSON_equal']
    for line in (SOURCE/'INPUTS.sha256').read_text().splitlines():
        h,n=line.split(None,1); assert sha(ROOT/n)==h
    assert len(manifest(SOURCE))==15
    architecture=manifest(ROOT/'fin_strategy_audit_20261019',root_relative=True)
    assert len(architecture)==47
    for n in ['fin_complete_anonymous_history_review_20261020','fin_post_rp004_history_review_20261020',
              'fin_intrinsic_entropy_source_review_20261019','fin_rp004_flux_review_20261019']:
        manifest(ROOT/n)
    claims=load(HERE/'CLAIM_REGISTER.json')
    assert not claims['actual_independent_F01_provider_exported'] and not claims['new_physical_law_adopted']
    assert claims['F01']==claims['F03']==claims['F06']=='OPEN_SOURCE'
    assert all((HERE/c['evidence']).is_file() for c in claims['claims'])
    for f in HERE.glob('*.json'): load(f)
    for f in HERE.glob('*.py'): ast.parse(f.read_text())
    for f in [*HERE.glob('*.md'),*[ROOT/n for n in ROOTDOCS]]:
        assert all(line==line.rstrip() for line in f.read_text().splitlines()),f
        for target in links(f.read_text()):
            if not target.startswith(('http:','https:','#')): assert (f.parent/target).exists(),(f,target)
    extra=['AGENTS.md','FIN_RESEARCH_STRATEGY.md','research_tasks/00_EXECUTION_RULES.md',
           'FIN_POST_AH_SOURCE_AND_VALIDATION_NEXT_CHAIN_20261020.md',
           'FIN_EVOLUTION_SOURCE_RESEARCH_STATUS_20261020.md',
           'fin_record_information_review_20261015/TAG_FEEDBACK_REMOVAL_LEMMA.md',
           'fin_generative_operational_review_20261014/EXACT_N2_PROCESS_AND_SHARP_SEPARATION.md',
           'fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py',
           'fin_rp004_independent_audit_20261019/HANDOFF.md',
           'fin_rp004_independent_audit_20261019/REPORT.md']
    names=sorted(set(architecture+extra+[str((SOURCE/n).relative_to(ROOT)) for n in actual]+d['source_input_paths']))
    inputs=''.join(f'{sha(ROOT/n)}  {n}\n' for n in names)
    binding={'date':'2026-10-20','entry_HEAD':'0c464f01','source_entries':15,'preserved_source_files':16,
             'source_input_matches':14,'new_bounded_controls':83,'bound_input_count':len(names),
             'verification_sha256':sha(HERE/'verification.json'),'AGENTS_sha256':d['live_AGENTS_sha256'],
             'root_documents':{n:sha(ROOT/n) for n in ROOTDOCS},'CF001_002_003':'CONDITIONAL_DERIVATIONS_SUPPLIED',
             'F01':'OPEN_SOURCE','actual_independent_source_supplied':False,'source_or_AGENTS_rewritten':False}
    if a.check:
        assert (HERE/'INPUTS.sha256').read_text()==inputs
        assert load(HERE/'INTEGRATION_BINDING.json')==binding
    else:
        (HERE/'INPUTS.sha256').write_text(inputs)
        (HERE/'INTEGRATION_BINDING.json').write_text(json.dumps(binding,indent=2)+'\n')
    files=sorted(f for f in HERE.iterdir() if f.is_file() and f.name!='MANIFEST.sha256')
    mtext=''.join(f'{sha(f)}  {f.name}\n' for f in files)
    if a.check: assert (HERE/'MANIFEST.sha256').read_text()==mtext
    else: (HERE/'MANIFEST.sha256').write_text(mtext)
    assert len(manifest(HERE))==len(files)
    print(json.dumps({'status':'PASS_ES_INTAKE_AND_CONDITIONAL_CAUSAL_TRANSFER_PACKAGE',
                      'read_only':a.check,'review_entries':len(files),'source_entries':15,
                      'new_controls':83,'F01':'OPEN_SOURCE'},indent=2))


if __name__=='__main__': main()
