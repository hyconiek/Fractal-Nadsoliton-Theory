"""Adjudicate freshly recomputed evidence. Does not replace a full jet replay."""
from fractions import Fraction as F
from itertools import permutations
from collections import Counter
from pathlib import Path
import json
import platform
import subprocess
import sys
import os

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fin_r7o3_review.review import HERE, ROOT, SOURCE, TAU, load, lines, sha, save


def add(a,b):
    return a[0]+b[0],a[1]+b[1]


def neg(a):
    return -a[1],-a[0]


def mul(a,b):
    v=[x*y for x in a for y in b]
    return min(v),max(v)


def det(A):
    out=(F(0),F(0))
    for p in permutations(range(len(A))):
        term=(F(1),F(1))
        for i,j in enumerate(p):
            term=mul(term,A[i][j])
        odd=sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))%2
        out=add(out,neg(term) if odd else term)
    return out


def proof_gate(cert, record):
    assert record['ok'] and cert['threshold']=='67/250'
    assert cert['certificate_sha256']==record['certificate_sha256']
    assert record['leaf_id']==cert['leaf_id']
    den=cert['basis_den']
    assert type(den) is int and den>0
    assert len(cert['basis_num'])==4 and all(len(r)==3 for r in cert['basis_num'])
    assert all(type(x) is int for r in cert['basis_num'] for x in r)
    B=[[F(x,den) for x in r] for r in cert['basis_num']]
    ids=cert['rank_minor_rows']
    assert len(set(ids))==3 and all(type(i) is int and 0<=i<4 for i in ids)
    minor=det([[(x,x) for x in B[i]] for i in ids])
    assert minor[0]==minor[1]!=0
    assert minor[0]==F(record['rank_minor_exact'])==F(cert['rank_minor_exact'])
    G=[[sum(B[k][i]*B[k][j] for k in range(4)) for j in range(3)] for i in range(3)]
    assert G==[list(map(F,r)) for r in record['gram_matrix_exact']]
    assert G==[list(map(F,r)) for r in cert['gram_matrix_exact']]
    E=[[tuple(map(F,p)) for p in row] for row in record['moment_entry_enclosures']]
    assert len(E)==3 and all(len(r)==3 for r in E)
    assert all(len(E[i][j])==2 and E[i][j][0]<=E[i][j][1] and E[i][j]==E[j][i]
               for i in range(3) for j in range(3))
    K=[[add((TAU*G[i][j],TAU*G[i][j]),neg(E[i][j])) for j in range(3)] for i in range(3)]
    # New exact-rational matrix check, independent of saved PD diagnostics.
    minors=[det([row[:n] for row in K[:n]])[0] for n in [1,2,3]]
    gersh=[K[i][i][0]-sum(max(map(abs,K[i][j])) for j in range(3) if i!=j) for i in range(3)]
    assert min(minors)>0 or min(gersh)>0
    return 'SYLVESTER' if min(minors)>0 else 'GERSHGORIN'


def completion_gate(record, certs):
    assert record['complete'] and not record['failed']
    assert record['processed']==record['total']==len(certs)==12425
    assert len(record['certificates'])==len(certs)
    assert [r['index'] for r in record['certificates']]==list(range(len(certs)))
    assert len({c['leaf_id'] for c in certs})==len(certs)


def main():
    registry=load(HERE/'registry.json')
    provenance=load(HERE/'provenance.json')
    assert provenance['archive_manifest_match'] and provenance['files_compared']==registry['manifest_entries']
    archive=ROOT/provenance['archive_path']
    archive_present=archive.exists()
    if archive_present:
        assert sha(archive)==provenance['archive_sha256']
    # Extracted proof inputs below remain mandatory and individually hashed.
    # A subsequently removed archive is not reconstructed or treated as present.
    assert registry['exact_geometry'] and registry['unresolved']==0
    assert registry['parents']==5432 and registry['active_leaves']==12425
    assert registry['prior_safe']==13231 and registry['expanded_compact_leaves']==25656
    assert sha(SOURCE/'MANIFEST.sha256')==registry['source_manifest_sha256']
    for rel,want in registry['source_sha256'].items():
        assert sha(SOURCE/rel)==want,rel
    for rel,want in registry['dependency_sha256'].items():
        assert sha(ROOT/rel)==want,rel
    certs=lines(SOURCE/'certificates/active_leaf_certificates.jsonl')
    parent_map={p['original_index']:p for p in load(SOURCE/'inherited/frozen_r7n_residual_5432.json')['parents']}
    checker_actual=sha(SOURCE/'src/fixed_witness_checker.py')
    checker_declared=Counter(c['checker_source_sha256'] for c in certs)
    for c in certs:
        assert c['original_parent_sha256']==parent_map[c['original_index']]['original_cell_sha256']
    replay=load(HERE/'leaf_replay.json')
    assert replay['input_sha256']==sha(SOURCE/'certificates/active_leaf_certificates.jsonl')
    completion_gate(replay,certs)
    counts=Counter()
    for c,r in zip(certs,replay['certificates']):
        counts[proof_gate(c,r)]+=1
    sample=load(HERE/'rational_sample.json')
    assert sample['passed']==sample['count']==len(sample['certificates'])
    assert all(r['ok'] and r['exact_match'] for r in sample['certificates'])
    selected={r['index'] for r in sample['certificates']}
    assert 10878 in selected
    assert {i for i,c in enumerate(certs) if c['pd_method']=='GERSHGORIN'}<=selected
    mathcheck=load(HERE/'analytic_checks.json')
    assert mathcheck['pass'] and F(mathcheck['sigma_upper'])<TAU
    assert mathcheck['current_spectral_provider_recomputed']
    for rel,want in mathcheck['provider_sha256'].items():
        assert sha(ROOT/rel)==want,rel
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
    test=subprocess.run([sys.executable,'-m','unittest','fin_r7o3_review.test_review',
                         'fin_r7o2_review.test_review','-v'],cwd=ROOT,env=env,
                         capture_output=True,text=True,timeout=90)
    assert test.returncode==0,test.stdout+test.stderr
    evidence=[p for p in HERE.iterdir() if p.suffix in ['.py','.md','.json'] and p.name!='verification.json']
    out=dict(status='PASS',date='2026-09-22',target_P_global=True,target_S_global=False,
             threshold='67/250',eigenvalue_order='descending; second largest',
             domain='Supplied C4 shared nonnegative fields J3,J4,J5,J6',
             manifest_entries=registry['manifest_entries'],source_unchanged=True,
             archive_sha256=provenance['archive_sha256'],
             archive_present_at_finalization=archive_present,
             archive_scope='Historical ZIP/extracted identity check; present archive rechecked only if available. Extracted inputs always fully checked.',
             parents_closed=5432,formula_leaves_recomputed=12425,unresolved=0,
             prior_accepted_safe_cells=13231,expanded_compact_leaves=25656,
             exact_pd_rechecks=dict(counts),rational_crosscheck_count=sample['count'],
             certificate_checker_hashes=dict(checker_declared),replayed_checker_sha256=checker_actual,
             checker_provenance_scope='Declared historical checker hashes are not assumed to identify the shipped checker; fresh replay uses and binds the actual shipped source',
             replay_seconds=replay['seconds'],test_log=test.stderr,
             tail_proof_status='Previously accepted FR1/FR42 reused with identity and domain checks; not newly re-proved',
             nonconclusions=['Target S','full-X7 ceiling or exhaustion','global minimizing-orbit uniqueness',
                             'gain/clock source','selector/QW-2191','laboratory evidence',
                             'legacy completion/role transfer','SM/GR','L_total','ToE'],
             check_scope='Integrity, complete replay ledger, and exact PD recheck; run review.py replay for all jet formulas',
             environment=dict(python=platform.python_version()),
             evidence_sha256={p.name:sha(p) for p in sorted(evidence)})
    if '--record' in sys.argv:
        save('verification.json',out)
    print(json.dumps({k:v for k,v in out.items() if k not in ['evidence_sha256','test_log']},indent=2))


if __name__=='__main__':
    if not __debug__:
        raise RuntimeError('Assertions must be enabled')
    main()
