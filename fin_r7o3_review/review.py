"""Read-only R7O3 intake; new audit outputs live exclusively in this directory.

Uses the inspected fixed-witness jet formulas with an independently implemented
outward binary64 backend. No witness generation or stored PASS is a proof gate.
Run with normal Python (not -O, since internal assertions are proof gates).
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import argparse
import hashlib
import json
import multiprocessing as mp
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = ROOT / 'FIN_R7O3_TARGETP_HANDOFF_20260920'
PRE = ROOT / 'FIN_R7N_HANDOFF_20260920'
TAU = F(67, 250)
sys.path.insert(0, str(ROOT))


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def load(p):
    return json.loads(p.read_text())


def lines(p):
    return [json.loads(s) for s in p.read_text().splitlines() if s.strip()]


def save(name, value):
    target = HERE / name
    temp = target.with_suffix(target.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2) + '\n')
    temp.replace(target)


def cell(raw):
    c = tuple(tuple(map(F, pair)) for pair in raw)
    assert len(c) == 4 and all(len(p) == 2 and p[0] < p[1] for p in c)
    return c


def tree(node, expected, leafmap, parent_id, used):
    assert cell(node['cell']) == expected
    if node['kind'] == 'SAFE':
        lid = node['leaf_id']
        assert lid not in used
        cert = leafmap[lid]
        assert cert['original_index'] == parent_id
        assert cell(cert['cell']) == expected
        used.add(lid)
        return
    assert node['kind'] == 'SPLIT'
    axis = node['axis']
    assert type(axis) is int and 0 <= axis < 4
    split = F(node['split'])
    lo, hi = expected[axis]
    assert lo < split < hi
    left, right = list(expected), list(expected)
    left[axis], right[axis] = (lo, split), (split, hi)
    tree(node['left'], tuple(left), leafmap, parent_id, used)
    tree(node['right'], tuple(right), leafmap, parent_id, used)


def compact_partition(rows):
    trie = {}
    for row in rows:
        path = row['path']
        assert len(path) % 2 == 0
        node = trie
        for i in range(0, len(path), 2):
            tag = path[i:i+2]
            assert tag[0] in '0123' and tag[1] in 'LR' and 'leaf' not in node
            node = node.setdefault(tag, {})
        assert not node
        node['leaf'] = cell(row['cell'])
    def walk(node):
        if 'leaf' in node:
            assert len(node) == 1
            return node['leaf']
        assert len(node) == 2
        axes = {int(k[0]) for k in node}
        assert len(axes) == 1
        axis = next(iter(axes))
        left, right = walk(node[f'{axis}L']), walk(node[f'{axis}R'])
        assert left[axis][1] == right[axis][0]
        assert all(left[k] == right[k] for k in range(4) if k != axis)
        out = list(left)
        out[axis] = (left[axis][0], right[axis][1])
        return tuple(out)
    return walk(trie)


def accepted_dependency(path, register):
    data = load(register)
    want = data['source_sha256'][str(path.relative_to(ROOT))]
    assert sha(path) == want, str(path)


def registry():
    manifest = {}
    for line in (SOURCE / 'MANIFEST.sha256').read_text().splitlines():
        if not line.strip():
            continue
        want, rel = line.split(maxsplit=1)
        rel = rel.lstrip('*')
        p = (SOURCE / rel).resolve()
        assert p.is_relative_to(SOURCE) and rel not in manifest
        assert sha(p) == want, rel
        manifest[rel] = want
    certs = lines(SOURCE / 'certificates/active_leaf_certificates.jsonl')
    leafmap = {r['leaf_id']: r for r in certs}
    assert len(certs) == len(leafmap) == 12425
    sys.path.insert(0, str(SOURCE / 'src'))
    from verification_core import validate_cert, tree_record_hash
    for cert in certs:
        assert not validate_cert(cert), cert['leaf_id']
        assert cert['threshold'] == '67/250'
        cell(cert['cell'])
    frozen = load(SOURCE / 'inherited/frozen_r7n_residual_5432.json')['parents']
    prior = load(PRE / 'checkpoints/R7N-021_t_refine_second_cheap_v1.json')
    original = prior['refined_failed']
    parents = load(SOURCE / 'parent_registry.json')['parents']
    trees = lines(SOURCE / 'certificates/reconstructed_parent_trees.jsonl')
    assert len(original) == len(frozen) == len(parents) == len(trees) == 5432
    assert {p['original_index'] for p in parents} == set(range(5432))
    assert {p['original_index'] for p in trees} == set(range(5432))
    assert {p['original_index'] for p in frozen} == set(range(5432))
    pmap = {p['original_index']: p for p in parents}
    tmap = {p['original_index']: p for p in trees}
    fmap = {p['original_index']: p for p in frozen}
    used = set()
    for i, old in enumerate(original):
        p, t, f = pmap[i], tmap[i], fmap[i]
        assert p['original_path'] == f['original_path'] == t['original_path'] == old['path']
        assert cell(p['original_cell']) == cell(f['original_cell']) == cell(old['cell'])
        assert tree_record_hash(t) == t['tree_sha256']
        local = set()
        tree(t['tree'], cell(old['cell']), leafmap, i, local)
        assert not used.intersection(local)
        assert len(local) == p['safe_terminal_count'] and p['unresolved_terminal_count'] == 0
        assert {leafmap[k]['certificate_sha256'] for k in local} == set(p['leaf_certificate_sha256s'])
        used.update(local)
    assert used == set(leafmap)
    # Reconstruct the exact accepted predecessor cells, not their claimed SAFE labels.
    base = load(PRE / 'checkpoints/R7N-020_trace_e2_compression_v1.json')
    safe = base['trace_terminals'] + base['e2_terminals'] + base['compression_terminals']
    predecessor_files = ['checkpoints/R7N-020_trace_e2_compression_v1.json']
    for name in ['R7N-021_refine_once_v1.json', 'R7N-021_t_refine_once_v2.json', 'R7N-021_t_refine_second_cheap_v1.json']:
        rel = 'checkpoints/' + name
        predecessor_files.append(rel)
        safe += load(PRE / rel)['refined_safe']
    inv = load(ROOT / 'fin_r7n_review/inventory.json')['archive_sha256']
    for rel in predecessor_files:
        assert sha(PRE / rel) == inv[rel]
    accepted_dependency(ROOT / 'fin_r7n_review/partial_leaf_replay.json', ROOT / 'fin_r7n_review/verification.json')
    paid = load(ROOT / 'fin_r7n_review/partial_leaf_replay.json')
    assert paid['complete'] and not paid['failed'] and paid['processed'] == paid['total'] == 13231
    assert len(safe) == 13231
    expected = {r['path']: (cell(r['cell']), True) for r in safe}
    assert len(expected) == len(safe)
    for row in original:
        assert row['path'] not in expected
        expected[row['path']] = (cell(row['cell']), False)
    compact = lines(SOURCE / 'inherited/r7n_final_compact_partition.jsonl')
    actual = {r['path']: (cell(r['cell']), r['safe']) for r in compact}
    assert len(compact) == len(actual) == 18663 and actual == expected
    assert compact_partition(compact) == cell(base['root_hull'])
    # Join to the already independently accepted tail premises, not new tail proofs.
    from global_join_checker import check
    join = check()
    assert join['global_join_pass']
    for row in join['tail_hashes']:
        assert sha(PRE / row['path']) == row['got'] == inv[row['path']]
    tail_register = ROOT / 'fin_rank7_intake_review/consolidated_verification.json'
    tail_evidence = load(tail_register)['evidence_sha256']
    tail_paths = [ROOT/'fin_rank7_intake_review'/name for name in
                  ['FR1_independent_geometry.json','FR1_replay.json','FR42_replay.json']]
    for p in tail_paths:
        assert sha(p) == tail_evidence[str(p.relative_to(ROOT))]
    # Check the mathematical feature enclosures against the accepted exact provider inputs.
    import mpmath as mpm
    iv = mpm.iv
    iv.dps = 65
    def point(q):
        q = F(q)
        return iv.mpf(q.numerator) / q.denominator
    def endpoint(t):
        sign, man, exponent, _ = t
        return (-1 if sign else 1) * F(man) * F(2) ** exponent
    spectrum = load(ROOT / 'fin_handoff_audit/results.json')['exact']['laplacian_intervals']
    assert spectrum == load(PRE / 'inputs/FR223_20260916/inputs/fin_handoff_audit/results.json')['exact']['laplacian_intervals']
    L = [point(F(a)) + (point(F(b))-point(F(a))) * iv.mpf([0,1]) for a,b in spectrum]
    obs = load(SOURCE / 'inherited/spectral_obs_rounded9.json')['obs']
    assert len(obs) == 7 and all(len(r) == 4 for r in obs)
    scales = [iv.sqrt(L[k]/n) for k,n in [(3,6),(4,6),(5,6),(6,12)]]
    for j, row in zip([0,4,6,2,3,5,1], obs):
        for a, k in enumerate([3,4,5,6]):
            v = scales[a] * (iv.cos(2*iv.pi*k*j/12) if k != 6 else (-1)**j)
            # Exact cosine zeros must not be tested against tiny interval pi roundoff.
            if k != 6 and (k*j) % 6 == 3:
                v = point(0)
            lo, hi = map(endpoint, v._mpi_)
            p, q = map(F, row[a])
            assert p <= lo <= hi <= q, (j, k, row[a])
    out = dict(manifest_entries=len(manifest), manifest_pass=True,
               source_manifest_sha256=sha(SOURCE/'MANIFEST.sha256'),
               parents=5432, active_leaves=12425, prior_safe=13231,
               original_compact_cells=18663, expanded_compact_leaves=25656,
               unresolved=0, exact_geometry=True, predecessor_identity=True,
               spectral_enclosures_checked=True, tail_join=join,
               source_sha256=manifest,
               dependency_sha256={str(p.relative_to(ROOT)):sha(p) for p in
                  [ROOT/'fin_r7n_review/verification.json',ROOT/'fin_r7n_review/partial_leaf_replay.json',
                   ROOT/'fin_r7n_review/inventory.json',tail_register,
                   ROOT/'fin_handoff_audit/results.json',
                   ROOT/'fin_r7o2_review/intervals_fast.py'] + tail_paths +
                   [PRE/rel for rel in predecessor_files]})
    save('registry.json', out)
    print(json.dumps({k:v for k,v in out.items() if k not in ['source_sha256','tail_join']},indent=2),flush=True)


def backend(fast=True):
    sys.path.insert(0, str(SOURCE / 'src'))
    import intervals
    if fast:
        from fin_r7o2_review.intervals_fast import Interval
        intervals.QI = Interval
    import fixed_witness_checker as ck
    if fast:
        ck.enable_outward_rounding = lambda digits: None
    return ck


def worker_init():
    global CHECKER
    CHECKER = backend(True)


def one(item):
    i, cert = item
    result = CHECKER.certify_fixed(cert)
    result.update(index=i, leaf_id=cert['leaf_id'], certificate_sha256=cert['certificate_sha256'])
    return result


def replay(limit=None, workers=2, rational=False):
    certs = lines(SOURCE / 'certificates/active_leaf_certificates.jsonl')
    if rational:
        chosen = sorted(set([0,1,2,3,4,10878,len(certs)-1] +
                            [i for i,c in enumerate(certs) if c['pd_method']=='GERSHGORIN']))
        ck = backend(False)
        output = []
        for i in chosen:
            r = ck.certify_fixed(certs[i],9)
            r.update(index=i, leaf_id=certs[i]['leaf_id'], exact_match=all(ck.exact_match(certs[i],r).values()))
            output.append(r)
            print('rational',i,r['ok'],r['exact_match'],flush=True)
        save('rational_sample.json',dict(count=len(output),passed=sum(r['ok'] and r['exact_match'] for r in output),certificates=output))
        return
    assert 1 <= workers <= 2
    if limit is not None:
        assert 0 < limit <= len(certs)
    total = len(certs) if limit is None else limit
    out = []
    start = time.monotonic()
    with mp.Pool(workers, initializer=worker_init) as pool:
        for r in pool.imap(one, enumerate(certs[:total]), chunksize=5):
            out.append(r)
            if len(out)%100==0 or len(out)==total:
                record=dict(total=len(certs),processed=len(out),complete=len(out)==len(certs),
                            failed=[r['index'] for r in out if not r['ok']],
                            seconds=time.monotonic()-start,workers=workers,
                            input_sha256=sha(SOURCE/'certificates/active_leaf_certificates.jsonl'),
                            arithmetic='Binary64 nextafter-outward operations with exact rational public endpoints',
                            certificates=out)
                save('leaf_replay.json' if limit is None else 'pilot_replay.json',record)
                print(len(out),'/',total,'failures',len(record['failed']),'seconds',round(record['seconds'],1),flush=True)


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Proof checker must not run with assertions disabled')
    ap = argparse.ArgumentParser()
    ap.add_argument('stage',choices=['registry','replay'])
    ap.add_argument('--limit',type=int)
    ap.add_argument('--workers',type=int,default=2)
    ap.add_argument('--rational',action='store_true')
    args=ap.parse_args()
    registry() if args.stage=='registry' else replay(args.limit,args.workers,args.rational)
