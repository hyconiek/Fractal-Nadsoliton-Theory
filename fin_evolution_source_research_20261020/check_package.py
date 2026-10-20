#!/usr/bin/env python3
"""Seal/check this new research packet without promoting its unmet source objective."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
ROOTDOC='FIN_EVOLUTION_SOURCE_RESEARCH_STATUS_20261020.md'
REQUIRED=['PROOF.md','PREMISE_LEDGER.md','REPORT_PL.md','SOURCE_ATOM.md',
          'LITERATURE_SCOPE.md','READING_AND_ACCESS.json','NONCONCLUSIONS.md',
          'CLAIM_REGISTER.json','HANDOFF.md','REPLAY.md','verify_research.py','RESULTS.json','check_package.py']


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def load(p):
    def reject(x): raise ValueError(x)
    return json.loads(p.read_text(),parse_constant=reject)


def manifest(folder,root_relative=False):
    count=0; seen=set()
    for line in (folder/'MANIFEST.sha256').read_text().splitlines():
        h,n=line.split(None,1)
        assert not Path(n).is_absolute() and '..' not in Path(n).parts and n not in seen
        assert sha((ROOT if root_relative else folder)/n)==h,(folder,n)
        seen.add(n); count+=1
    return count


def local_links(text):
    out=[]; fence=None
    for line in text.splitlines():
        s=line.lstrip(); marker='```' if s.startswith('```') else '~~~' if s.startswith('~~~') else None
        if marker: fence=marker if fence is None else None if fence==marker else fence
        elif fence is None and not line.startswith(('    ','\t')): out.append(line)
    return re.findall(r'\[[^\]]+\]\(([^)]+)\)',re.sub(r'`[^`]*`','','\n'.join(out)))


def main():
    p=argparse.ArgumentParser(); p.add_argument('--check',action='store_true'); args=p.parse_args()
    assert all((HERE/n).is_file() for n in REQUIRED)
    results=load(HERE/'RESULTS.json'); claims=load(HERE/'CLAIM_REGISTER.json')
    assert results['status']=='PASS_BOUNDED_RESET_HISTORY_SOURCE_PROBE'
    assert results['bounded_control_count']==len(results['bounded_controls'])==len(set(results['bounded_controls']))==122
    assert results['F01']==claims['F01']=='OPEN_SOURCE_NOT_DERIVED'
    assert not claims['new_physical_postulate_adopted'] and not claims['Markov_inverse_and_AH_theorems_refuted']
    for rel,h in results['input_hashes'].items(): assert sha(ROOT/rel)==h,rel
    assert all((HERE/c['evidence']).is_file() for c in claims['claims'])
    for n in ['fin_complete_anonymous_history_review_20261020','fin_post_rp004_history_review_20261020',
              'fin_intrinsic_entropy_source_review_20261019','fin_rp004_flux_review_20261019']:
        manifest(ROOT/n)
    assert manifest(ROOT/'fin_strategy_audit_20261019',root_relative=True)==47
    for f in HERE.glob('*.json'): load(f)
    for f in HERE.glob('*.py'): ast.parse(f.read_text())
    for f in [*HERE.glob('*.md'),ROOT/ROOTDOC]:
        assert all(line==line.rstrip() for line in f.read_text().splitlines()),f
        for target in local_links(f.read_text()):
            if not target.startswith(('http:','https:','#')): assert (f.parent/target).exists(),(f,target)
    extra=['FIN_RESEARCH_STRATEGY.md','FIN_POST_AH_SOURCE_AND_VALIDATION_NEXT_CHAIN_20261020.md',
           'FIN_COMPLETE_ANONYMOUS_HISTORY_INTEGRATION_20261020.md',
           'fin_unified_strict_evolution_audit_20261017/STRICT_EVOLUTION_SOURCE_AUDIT.md',
           'fin_unified_strict_evolution_audit_20261017/MULTISCALE_FRACTAL_DYNAMICS.md',
           'fin_toe_bridge_research_20261015/REPOSITORY_FRONTIER_AUDIT.md']
    names=sorted(set(list(results['input_hashes'])+extra))
    inputs=''.join(f'{sha(ROOT/n)}  {n}\n' for n in names)
    binding={'date':'2026-10-20','bounded_controls':122,'results_sha256':sha(HERE/'RESULTS.json'),
             'root_status_sha256':sha(ROOT/ROOTDOC),'input_count':len(names),
             'F01':'OPEN_SOURCE_NOT_DERIVED','positive_source_object_exported':False,
             'new_counterfamily':'diagnostic reset-history family, not FIN adoption',
             'independent_provider_condition':'SOURCE_ATOM.md, unpaid',
             'external_pdf_used_as_proof_dependency':False,
             'old_sources_and_AGENTS_rewritten':False}
    if args.check:
        assert (HERE/'INPUTS.sha256').read_text()==inputs
        assert load(HERE/'INTEGRATION_BINDING.json')==binding
    else:
        (HERE/'INPUTS.sha256').write_text(inputs)
        (HERE/'INTEGRATION_BINDING.json').write_text(json.dumps(binding,indent=2)+'\n')
    files=sorted(f for f in HERE.iterdir() if f.is_file() and f.name!='MANIFEST.sha256')
    mtext=''.join(f'{sha(f)}  {f.name}\n' for f in files)
    if args.check: assert (HERE/'MANIFEST.sha256').read_text()==mtext
    else: (HERE/'MANIFEST.sha256').write_text(mtext)
    assert manifest(HERE)==len(files)
    print(json.dumps({'package':'PASS_BOUNDED_RESEARCH_PACKAGE_ONLY','entries':len(files),
                      'controls':122,'F01':'OPEN_SOURCE_NOT_DERIVED','read_only':args.check},indent=2))


if __name__=='__main__': main()
