#!/usr/bin/env python3
"""Check an append-only live instruction update; preserve old FIN inputs."""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
N = 1252727
BASE = '300f3864fde7843bd46f53840baebfc8190e3501c5d19de394a7266d3c1367e4'
PACKAGES = {
    'fin_evolution_source_review_20261020': ('finalize_review.py','--check'),
    'fin_evolution_source_research_20261020': ('check_package.py','--check'),
    'fin_complete_anonymous_history_review_20261020': ('finalize_review.py','--check'),
    'fin_post_rp004_history_review_20261020': ('finalize_review.py','--check'),
    'fin_rp004_flux_review_20261019': ('finalize_review.py','--check'),
    'fin_intrinsic_entropy_source_review_20261019': ('finalize_review.py','--check'),
    'fin_strategy_audit_20261019': ('verify_strategy.py','--inputs-only'),
}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def sha(p):
    return digest(p.read_bytes())


def load(p):
    return json.loads(p.read_text(),parse_constant=lambda s: (_ for _ in ()).throw(ValueError(s)))


def check_prefix(data, c):
    n,m = c['historical_bytes'],c['updated_bytes']
    if not 0 < n < m <= len(data):
        raise ValueError('truncated instruction revision')
    if digest(data[:n]) != c['historical_sha256']:
        raise ValueError('historical instruction prefix changed')
    if digest(data[n:m]) != c['append_sha256'] or digest(data[:m]) != c['updated_sha256']:
        raise ValueError('new reviewed instruction prefix changed')
    return len(data)-m


def snapshot(folder):
    files = sorted(p for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    listing = ''.join(f'{sha(p)}  {p.relative_to(folder).as_posix()}\n' for p in files)
    return {'files':len(files),'tree_sha256':digest(listing.encode()),'manifest_sha256':sha(folder/'MANIFEST.sha256')}


def manifest(folder, root_relative=False):
    names = []
    for row in (folder/'MANIFEST.sha256').read_text().splitlines():
        h,name = row.split(None,1)
        p = Path(name)
        assert not p.is_absolute() and '..' not in p.parts and name not in names
        assert sha((ROOT if root_relative else folder)/name) == h,(folder,name)
        names.append(name)
    return len(names)


def limits():
    resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))


def run(script, args, cwd):
    r = subprocess.run(['python3',str(script),*args],cwd=cwd,
        env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'),
        capture_output=True,text=True,timeout=60,preexec_fn=limits)
    return {'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr}


def replay(prefix):
    logs = {}
    with tempfile.TemporaryDirectory(prefix='fin-causal-historical-context-') as d:
        tree = Path(d)
        for p in ROOT.iterdir():
            if p.name.startswith('.') or p.name == HERE.name:
                continue
            if p.name in PACKAGES:
                shutil.copytree(p,tree/p.name,ignore=shutil.ignore_patterns('__pycache__'))
            elif p.name == 'AGENTS.md':
                (tree/p.name).write_bytes(prefix)
            else:
                (tree/p.name).symlink_to(p,target_is_directory=p.is_dir())
        assert sha(tree/'AGENTS.md') == BASE
        for name,(script,*args) in PACKAGES.items():
            assert sha(tree/name/script) == sha(ROOT/name/script)
            logs[name] = run(tree/name/script,args,tree)
            logs[name]['scope'] = 'ARCHIVED_CONTEXT_ONLY; historical AGENTS input, not live instructions or new science'
            assert logs[name]['returncode'] == 0,(name,logs[name])
    return logs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--record',action='store_true')
    ap.add_argument('--check',action='store_true')
    ap.add_argument('--replay-archives',action='store_true')
    ap.add_argument('--output',type=Path)
    args = ap.parse_args()
    data = (ROOT/'AGENTS.md').read_bytes()
    append = (HERE/'AGENTS_APPEND.md').read_bytes()
    expected = {'date':'2026-10-20','historical_bytes':N,'historical_sha256':BASE,
                'append_bytes':len(append),'append_sha256':digest(append),
                'updated_bytes':N+len(append),'updated_sha256':digest(data[:N]+append)}
    if args.record:
        assert digest(data[:N]) == BASE and data[N:] == append
        (HERE/'UPDATE_CONTRACT.json').write_text(json.dumps(expected,indent=2)+'\n')
    contract = load(HERE/'UPDATE_CONTRACT.json')
    assert contract == expected
    extra = check_prefix(data,contract)
    negative = {}
    for label,mutant in [('old_byte',bytes([data[0]^1])+data[1:]),
                         ('new_byte',data[:N]+bytes([data[N]^1])+data[N+1:]),
                         ('truncation',data[:contract['updated_bytes']-1])]:
        try:
            check_prefix(mutant,contract)
        except ValueError:
            negative[label] = 'REJECTED'
        else:
            raise AssertionError(label)
    assert check_prefix(data[:contract['updated_bytes']]+b'\nUnreviewed later append\n',contract) > 0
    protected = {name:snapshot(ROOT/name) for name in PACKAGES}
    counts = {name:manifest(ROOT/name,name == 'fin_strategy_audit_20261019') for name in PACKAGES}
    research = ROOT/'fin_reversible_source_transfer_research_20261020'
    science = {'manifest_entries':manifest(research),'manifest_sha256':sha(research/'MANIFEST.sha256')}
    rp_paths = ['fin_rp004_flux_review_20261019/'+n for n in
                ['MOVING_CONE_THEOREM.md','BOOTSTRAP_COMPLETION.md',
                 'OUTWARD_FLUX_THEOREM.md','COMPLETE_FLUX_BUDGET.json']]
    names = sorted([name+'/MANIFEST.sha256' for name in PACKAGES] + rp_paths +
        ['fin_reversible_source_transfer_research_20261020/MANIFEST.sha256',
         'FIN_CAUSAL_SOURCE_TRANSFER_NEXT_CHAIN_20261020.md',
         'FIN_POST_RP004_RESEARCH_DECISION_20261019.md',
         'research_tasks/00_EXECUTION_RULES.md','research_tasks/RP004_SHARED_CONTRACT.json',
         'fin_coherent_mode_radiation_review_20261019/PROOF_COMPLETIONS_AND_SCOPE.md',
         'fin_nonlinear_stability_transfer_review_20261017/PROOF_COMPLETIONS_AND_SCOPE.md'])
    input_text = ''.join(f'{sha(ROOT/name)}  {name}\n' for name in names)
    live,archives = {},{}
    if args.replay_archives:
        live = run(ROOT/'fin_evolution_source_review_20261020/finalize_review.py',['--check'],ROOT)
        assert live['returncode'] != 0 and 'AGENTS.md' in live['stderr'],live
        live['expected_status'] = 'HISTORICAL_FULL_FILE_HASH_REJECTED; not a scientific failure or live PASS'
        archives = replay(data[:N])
    assert data == (ROOT/'AGENTS.md').read_bytes()
    assert protected == {name:snapshot(ROOT/name) for name in PACKAGES}
    result = {'date':'2026-10-20','status':'PASS_SCOPED_APPEND_ONLY_INSTRUCTION_UPDATE' if extra == 0 else 'REVIEWED_PREFIX_WITH_LATER_UNREVIEWED_TAIL',
              'unreviewed_tail_bytes':extra,'historical_prefix_sha256':digest(data[:N]),
              'reviewed_update_sha256':contract['updated_sha256'],'negative_controls':negative,
              'protected_packages':protected,'protected_manifest_entries':counts,
              'new_RS_package':science,'direct_live_old_checker':live,
              'RP_proofs_availability_binding_only':{name:sha(ROOT/name) for name in rp_paths},
              'archived_context_only_replays':archives,'F01_F03_F06':'OPEN_SOURCE',
              'CSP001_003':'NOT_ACTIVATED','new_independent_RP_audit_executed':False,
              'new_physical_law':False,'commit_or_push':False}
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    if args.record:
        assert args.replay_archives and args.output and args.output.resolve() == HERE/'RESULTS.json'
        (HERE/'INPUTS.sha256').write_text(input_text)
        files = sorted(p for p in HERE.iterdir() if p.is_file() and p.name != 'MANIFEST.sha256')
        (HERE/'MANIFEST.sha256').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in files))
    if args.check:
        manifest(HERE)
        assert (HERE/'INPUTS.sha256').read_text() == input_text
        saved = load(HERE/'RESULTS.json')
        assert saved['protected_packages'] == protected and saved['new_RS_package'] == science
        assert len(saved['archived_context_only_replays']) == 7
        assert all(x['returncode'] == 0 for x in saved['archived_context_only_replays'].values())
        if args.replay_archives:
            assert result == saved
        names = {line.split(None,1)[1] for line in (HERE/'MANIFEST.sha256').read_text().splitlines()}
        assert names|{'MANIFEST.sha256'} == {p.name for p in HERE.iterdir() if p.is_file()}
    for f in HERE.glob('*.json'):
        load(f)
    ast.parse(Path(__file__).read_text())
    print(json.dumps({'status':result['status'],'historical_prefix_unchanged':True,
                      'archived_context_replays':len(archives),'protected_packages':len(protected),
                      'F01':'OPEN_SOURCE','CSP':'NOT_ACTIVATED'}))


if __name__ == '__main__':
    main()
