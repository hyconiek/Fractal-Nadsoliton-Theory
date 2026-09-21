"""Integrity and acceptance gate after the full R7O2 mathematical replay."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import hashlib,json,os,subprocess,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
SOURCE=ROOT/'FIN_R7O2_CONTINUATION_HANDOFF_20260920'
def load(name):return json.loads((HERE/name).read_text())
def det(A):return A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])-A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0])
def main():
    n=0
    for line in (SOURCE/'MANIFEST.sha256').read_text().splitlines():
        if not line.strip():continue
        want,rel=line.split(maxsplit=1);p=SOURCE/rel.lstrip('*')
        assert hashlib.sha256(p.read_bytes()).hexdigest()==want,rel;n+=1
    reg=load('registry.json');data=load('leaf_replay.json');rows=load('safe_leaves.json')
    assert data['complete'] and data['processed']==data['total']==len(rows)==7340 and not data['failed']
    assert len(data['certificates'])==len(rows)
    for raw,c in zip(rows,data['certificates']):
        assert c['original_index']==raw['original_index'] and c['path']==raw['path'] and c['ok']
        assert c['basis_num']==raw['basis_num'] and c['basis_den']==raw['basis_den']
        assert len(c['center_c'])==3;list(map(F,c['center_c']))
        B=[[F(x,c['basis_den']) for x in r] for r in c['basis_num']]
        assert F(c['rank_minor']) and F(c['rank_minor']) in [det([B[i] for i in inds]) for inds in combinations(range(4),3)]
        bounds=[tuple(map(F,p)) for p in c['pd_bounds']]
        assert all(a<=b for a,b in bounds)
        assert all(a>0 for a,b in bounds) or F(c['gersh_lower'])>0
    missing=reg['missing'];repair=reg['repair_parents']
    assert len(missing)==1972 and len(repair)==54 and reg['claimed_closed']==3406
    parent_res=F(reg['unprocessed_fraction'])+F(reg['repair_parent_fraction'])
    leaf_res=F(reg['unprocessed_fraction'])+F(reg['unresolved_terminal_fraction'])
    previous=json.loads((ROOT/'fin_r7n_review/partial_geometry.json').read_text())
    assert F(reg['closed_parent_fraction'])+parent_res==F(previous['residual_fraction'])
    for name in ['sample_binary_replay.json','sample_rational_replay.json']:
        sample=load(name);assert sample['processed']==20 and not sample['failed']
        assert [c['center_c'] for c in sample['certificates']]==[c['center_c'] for c in data['certificates'][:20]]
    env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1'
    tests=subprocess.run([sys.executable,'-m','unittest','fin_r7o2_review.test_review','-v'],cwd=ROOT,env=env,capture_output=True,text=True,timeout=60)
    assert tests.returncode==0,tests.stderr;assert 'Ran 6 tests' in tests.stderr
    depnames=['src/target_p_trace_cover.py','src/compression_interval_probe.py','src/generic_threshold_shifted.py',
              'inputs/FR223_20260916/src/intervals.py','inputs/FR223_20260916/src/off_face.py',
              'inputs/FR223_20260916/src/off_face_local.py','inputs/FR223_20260916/src/boundary_ising.py',
              'inputs/FR223_20260916/inputs/fin_handoff_audit/results.json']
    prior=json.loads((ROOT/'fin_r7n_review/inventory.json').read_text())['archive_sha256'];deps={}
    for name in depnames:
        p=ROOT/'FIN_R7N_HANDOFF_20260920'/name;h=hashlib.sha256(p.read_bytes()).hexdigest();assert h==prior[name];deps[name]=h
    out=dict(status='PASS',date='2026-09-20',manifest_entries=n,source_archive_unchanged=True,
             formula_leaves_replayed=7340,processed_parents=3460,certified_closed_parents=3406,
             partial_parents=54,unresolved_terminals=171,unprocessed_parents=1972,
             conservative_parent_residual=str(parent_res),terminal_residual=str(leaf_res),
             conservative_parent_residual_decimal=float(parent_res),terminal_residual_decimal=float(leaf_res),
             compact_hull_certified_fraction_decimal=float(1-leaf_res),global_target_P=False,
             new_tests=6,test_log=tests.stderr,dependency_sha256=deps,
             fixed_witness_cross_backend_sample=20,
             mathematical_scope='Partial shared-field C4 Target P only. No new parent or unresolved-leaf research was run.',
             checker_scope='Completed replay records plus original source/dependency integrity; full formulas are rerun by review.py replay.',
             evidence_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.iterdir()) if p.suffix in ['.py','.md','.json'] and p.name!='verification.json'})
    if '--record' in sys.argv:(HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['test_log','dependency_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
