#!/usr/bin/env python3
"""Bounded explicit-path rg search, not an exhaustive source inventory."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
QUERIES = {
    'pair_reversible': r'(two.time|pair.?wise|two.point|dwuczas|pair semigroup).{0,140}(reversib|Markov|semigroup)|reversib.{0,140}(lumpab|two.time|pair.?wise)',
    'compression': r'PL.?2P|P L.?2 P|second.derivative|compression defect|variance.{0,30}compression',
    'finite_lag_transfer': r'trace.{0,90}(semigroup|defect)|semigroup.{0,90}trace|two.lag|doubling|sampled.record.{0,90}defect',
    'source_clock': r'ACTIVITY_SOURCE|SOURCE_GATE|relative.activation.barrier|spectral.inclusion|11.fold|eleven.fold',
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def selected_paths():
    index = json.loads((ROOT/'fin_strategy_audit_20261019/SEARCH_FILE_INDEX.json').read_text())
    paths = {p for p in index['selected_paths'] if (ROOT/p).is_file()}
    for pattern in ('fin_*20261020/*.md', 'fin_auxiliary_mediator*/*.md',
                    'fin_intermediate_clock*/*.md', 'FIN_*SOURCE*.md',
                    'fin_transition_hierarchy/*.md', 'fin_transition_hierarchy/*.py'):
        paths.update(p.relative_to(ROOT).as_posix() for p in ROOT.glob(pattern) if p.is_file())
    result = []
    excluded = {'inputs', 'source_context', 'REPLAY', '__pycache__', 'embedded', 'archives'}
    for name in sorted(paths):
        p = Path(name)
        if p.suffix not in ('.md', '.tex', '.py') or p.name == 'AGENTS.md':
            continue
        if HERE.name in p.parts or excluded.intersection(p.parts):
            continue
        if p.name in ('AGENTS_INPUT.md', 'AGENTS_latest.md'):
            continue
        if p.suffix == '.tex' and (ROOT/p.with_suffix('.md')).is_file():
            continue
        result.append(name)
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--refresh', action='store_true', help='explicitly reselect the indexed search scope')
    args = ap.parse_args()
    saved = HERE/'SEARCH_RESULTS.json'
    names = sorted(json.loads(saved.read_text())['searched_files']) if saved.is_file() and not args.refresh else selected_paths()
    before = {p:sha(ROOT/p) for p in names}
    matches = {}
    for label, pattern in QUERIES.items():
        lines = []
        for start in range(0, len(names), 100):
            batch = names[start:start+100]
            run = subprocess.run(['rg', '-n', '-i', '--max-columns', '400', '--', pattern, *batch],
                                 cwd=ROOT, text=True, capture_output=True, timeout=60)
            if run.returncode not in (0, 1):
                raise RuntimeError(run.stderr)
            lines.extend(run.stdout.splitlines())
        matches[label] = {'pattern':pattern, 'matching_lines':len(lines), 'lines':lines}
    assert before == {p:sha(ROOT/p) for p in names}, 'search inputs changed during scan'
    out = {'date':'2026-10-20', 'selected_existing_files':len(names),
           'scope':'bounded indexed text search; not all repo, not all proofs read',
           'selection_index_sha256':sha(ROOT/'fin_strategy_audit_20261019/SEARCH_FILE_INDEX.json'),
           'searched_files':before, 'queries':matches,
           'scientific_novelty_claim':'local interface delta only; no literature-priority or global absence claim'}
    args.output.write_text(json.dumps(out, indent=2, ensure_ascii=False)+'\n')
    print(json.dumps({'selected_existing_files':len(names),
                      'matching_lines':{k:v['matching_lines'] for k,v in matches.items()}}))


if __name__ == '__main__':
    main()
