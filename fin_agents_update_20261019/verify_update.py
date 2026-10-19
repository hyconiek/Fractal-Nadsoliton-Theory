#!/usr/bin/env python3
"""Verify a new instruction prefix without modifying historical FIN packages."""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
PACKAGES={
    'fin_nonlinear_stability_transfer_review_20261017':'finalize_review.py',
    'fin_coupled_evolution_source_review_20261017':'finalize_review.py',
    'fin_unified_strict_evolution_audit_20261017':'finalize_audit.py',
    'fin_dense_orbit_fibre_transport_review_20261017':'finalize_review.py',
    'fin_intermediate_clock_bridge_research_20261017':'finalize_research.py',
    'fin_balanced_phase_current_review_20261017':'finalize_review.py',
    'fin_q12_physical_bridge_research_20261017':'finalize_research.py',
}
def digest(data):return hashlib.sha256(data).hexdigest()
def sha(path):return digest(path.read_bytes())
def read_json(path):
    def reject(v):raise ValueError(v)
    return json.loads(path.read_text(),parse_constant=reject)
def check_prefix(data,c):
    n=c['base_prefix_bytes'];m=c['updated_prefix_bytes']
    if not (0<n<m<=len(data)):raise ValueError('truncated instruction revision')
    if digest(data[:n])!=c['base_prefix_sha256']:raise ValueError('historical instruction prefix changed')
    if digest(data[:m])!=c['updated_prefix_sha256']:raise ValueError('reviewed instruction revision changed')
    if digest(data[n:m])!=c['append_sha256']:raise ValueError('append changed')
    return len(data)-m
def manifest(folder):
    names=[]
    for line in (folder/'MANIFEST.sha256').read_text().splitlines():
        h,n=line.split('  ',1);p=Path(n)
        assert not p.is_absolute() and '..' not in p.parts
        assert p.as_posix() not in names and sha(folder/p)==h,(folder,n)
        names.append(p.as_posix())
    actual={p.relative_to(folder).as_posix() for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    assert actual==set(names)|{'MANIFEST.sha256'},folder
    return {'entries':len(names),'manifest_sha256':sha(folder/'MANIFEST.sha256')}
def archived_checks(prefix):
    # Only --check finalizers execute; no source writer follows these symlinks.
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
    logs={}
    with tempfile.TemporaryDirectory(prefix='fin-agents-archive-context-') as temp:
        tree=Path(temp)
        for p in ROOT.iterdir():
            if p.name.startswith('.') or p.name==HERE.name:continue
            destination=tree/p.name
            if p.name in PACKAGES:
                shutil.copytree(p,destination,ignore=shutil.ignore_patterns('__pycache__'))
            elif p.name=='AGENTS.md':
                destination.write_bytes(prefix)
            else:destination.symlink_to(p,target_is_directory=p.is_dir())
        assert (tree/'AGENTS.md').read_bytes()==prefix
        for name,script in PACKAGES.items():
            assert sha(tree/name/script)==sha(ROOT/name/script)
            run=subprocess.run(['python3',str(tree/name/script),'--check'],cwd=tree,
                               env=env,capture_output=True,text=True,timeout=30)
            logs[name]={'returncode':run.returncode,'stdout':run.stdout,'stderr':run.stderr,
                        'scope':'ARCHIVED_CONTEXT with exact historical AGENTS, not current instruction validation'}
            assert run.returncode==0,(name,run.stderr)
    return {'scope':'ARCHIVED_CONTEXT_ONLY','temporary_tree_removed':True,
            'baseline_AGENTS_sha256':digest(prefix),'finalizers':logs}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--replay-archives',action='store_true')
    ap.add_argument('--output',type=Path);ap.add_argument('--check',action='store_true')
    ap.add_argument('--seal',action='store_true');args=ap.parse_args()
    c=read_json(HERE/'UPDATE_CONTRACT.json');data=(ROOT/'AGENTS.md').read_bytes()
    extra=check_prefix(data,c);n=c['base_prefix_bytes'];m=c['updated_prefix_bytes']
    neg={}
    for name,mutated in (
        ('changed_historical_byte',bytes([data[0]^1])+data[1:]),
        ('changed_new_instruction_byte',data[:n]+bytes([data[n]^1])+data[n+1:]),
        ('truncated_revision',data[:m-1])):
        try:check_prefix(mutated,c)
        except ValueError:neg[name]='REJECTED'
        else:raise AssertionError(name)
    assert check_prefix(data[:m]+b'\nUnreviewed later append.\n',c)>0
    evidence={name:manifest(ROOT/name) for name in PACKAGES}
    refs=[ROOT/'FIN_COHERENT_MODE_RADIATION_NEXT_CHAIN_20261017.md',
          ROOT/'fin_nonlinear_stability_transfer_review_20261017/CLAIM_REGISTER.json',
          ROOT/'fin_nonlinear_stability_transfer_review_20261017/PROOF_COMPLETIONS_AND_SCOPE.md']
    ref_hashes={str(p.relative_to(ROOT)):sha(p) for p in refs}
    assert 'RP-001–RP-005 NOT EXECUTED' in refs[0].read_text()
    ns=read_json(refs[1]);assert ns['actual_FIN_parameter_membership'] is False
    assert ns['new_PM_PB_PP_EV_adopted'] is False
    archives={};live_check={}
    if args.replay_archives:
        direct=subprocess.run(['python3',str(ROOT/'fin_nonlinear_stability_transfer_review_20261017/finalize_review.py'),'--check'],
                              cwd=ROOT,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),capture_output=True,text=True,timeout=30)
        assert direct.returncode!=0 and 'AGENTS.md' in direct.stderr
        live_check={'expected_result':'HISTORICAL_FULL_FILE_HASH_MISMATCH',
                    'returncode':direct.returncode,'diagnostic':'AGENTS.md assertion; not falsely advertised as a live pass'}
        archives=archived_checks(data[:n])
    assert (ROOT/'AGENTS.md').read_bytes()==data
    assert {name:manifest(ROOT/name) for name in PACKAGES}==evidence
    result={'date':'2026-10-19',
            'status':'PASS_APPEND_ONLY_FIN_AGENTS_UPDATE' if extra==0 else 'REVIEWED_PREFIX_WITH_UNREVIEWED_LATER_APPEND',
            'actual_live_AGENTS_sha256':digest(data),'base_prefix_sha256':digest(data[:n]),
            'reviewed_revision_sha256':digest(data[:m]),'unreviewed_tail_bytes':extra,
            'negative_controls':neg,'scientific_packages_unchanged':evidence,
            'reference_hashes':ref_hashes,'direct_live_archive_check':live_check,'archived_replay':archives,
            'scope':'current instruction-prefix verification plus separately labelled historical-context package checks; no new science',
            'RP':'NOT_ACTIVATED','new_physical_claims':False}
    if args.output:args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if args.seal:
        assert args.replay_archives and extra==0 and args.output and args.output.resolve()==(HERE/'RESULTS.json').resolve()
        files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256')
        (HERE/'MANIFEST.sha256').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in files))
    if args.check:
        manifest(HERE)
        saved=read_json(HERE/'RESULTS.json')
        assert saved['reviewed_revision_sha256']==c['updated_prefix_sha256']
        assert saved['scientific_packages_unchanged']==evidence and saved['reference_hashes']==ref_hashes
        assert saved['archived_replay']['scope']=='ARCHIVED_CONTEXT_ONLY'
        assert len(saved['archived_replay']['finalizers'])==7
        assert all(v['returncode']==0 for v in saved['archived_replay']['finalizers'].values())
        for p in HERE.glob('*.json'):read_json(p)
        ast.parse(Path(__file__).read_text())
    print(json.dumps({'status':result['status'],'historical_prefix_unchanged':True,
                      'unreviewed_tail_bytes':extra,'scientific_package_manifests':len(evidence),
                      'historical_context_finalizers':len(archives.get('finalizers',{})),
                      'saved_package_checked':args.check,'RP':'NOT_ACTIVATED'},indent=2))
if __name__=='__main__':main()
