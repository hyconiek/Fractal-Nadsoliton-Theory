"""Check integrity/structure of the newly replayed evidence and record the handoff.

This does not replace the scientific replay commands documented in README.md.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import re
import sys
import os
import subprocess

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
def load(name):return json.loads((HERE/name).read_text())
def interval(x):
    a,b=map(F,x);assert a<=b;return a,b
def shifted_ok(p):
    assert interval(p['c2'])[0]>0
    assert interval(p['P'])[1]<=0 or interval(p['P1'])[0]>=0

def validate_repair(repair):
    assert repair['status']=='PASS'
    bypath={x['path']:x for x in repair['leaves']};assert len(bypath)==len(repair['leaves'])
    seen=set()
    def visit(path,box):
        if path in bypath:
            row=bypath[path];assert row['box']==[[str(x) for x in v] for v in box]
            assert row['status']=='INTERVAL_CERTIFIED';shifted_ok(row['bounds']);seen.add(path);return
        children=[p for p in bypath if p.startswith(path)];assert children
        axes={int(p[len(path)]) for p in children};assert len(axes)==1
        axis=next(iter(axes));assert 0<=axis<4
        lo,hi=box[axis];mid=(lo+hi)/2
        left=list(box);right=list(box);left[axis]=(lo,mid);right[axis]=(mid,hi)
        visit(path+str(axis)+'L',left);visit(path+str(axis)+'R',right)
    visit('',[tuple(map(F,p)) for p in repair['original_box']]);assert seen==set(bypath)
    return len(seen)

def main():
    env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1'
    controls=subprocess.run([sys.executable,'-m','unittest','fin_rank7_intake_review.test_integration','-v'],
                            cwd=ROOT,env=env,capture_output=True,text=True,timeout=60)
    assert controls.returncode==0,controls.stdout+controls.stderr
    assert 'Ran 7 tests' in controls.stderr
    inv=load('continuation_inventory.json')
    for rel,want in inv['source_sha256'].items():
        f=ROOT/'FIN_rank7_CONTINUATION_HANDOFF_20260916_FR223'/rel
        assert hashlib.sha256(f.read_bytes()).hexdigest()==want,rel
    tests=load('replay_latest_fresh.json');assert len(tests)==29
    count=0
    for row in tests:
        assert row['returncode']==0
        assert 'skipped=' not in row['stderr']
        count+=int(re.search(r'Ran (\d+) tests?',row['stderr'])[1])
    assert count==119
    union=load('FR223_union_replay.json');assert union['finished'] and union['count']==106
    failures=[]
    for row in union['rows']:
        if row['status']!='INTERVAL_CERTIFIED':failures.append(row['name']);continue
        p=row['proof']
        if row['kind']=='shifted':shifted_ok(p)
        else:
            assert all(p['signs'].values())
            for k in ['boundary_schur','endpoint_schur']:
                assert interval(p[k][0])[1]<0
                assert all(interval(v)[0]>0 for v in p[k][1:])
    assert len(failures)==7 and union['passed']==99
    repairs=load('FR223_subdivision_repairs.json');assert {x['name'] for x in repairs}==set(failures)
    total_leaves=0
    for repair in repairs:
        total_leaves+=validate_repair(repair)
    assert total_leaves==18
    phases=load('phase_recertification.json')
    for kind in ['quartic','full']:
        rows=phases[kind]['roots'];assert len(rows)==60
        for row in rows:
            radius=F(row['radius']);assert F(row['contraction_upper'])<1
            assert all(-radius<interval(v)[0]<=interval(v)[1]<radius for v in row['krawczyk_offsets'])
            signs=[]
            for v in row['ldl_pivots']:
                lo,hi=interval(v);assert lo>0 or hi<0;signs.append(hi<0)
            assert sum(signs)==row['negative_index']
    tail=load('FR42_replay.json');assert tail['status']=='PASS' and tail['leaves']==637
    assert F(tail['all_66_pair_diameter_squared_upper'])<5
    assert F(tail['perturbation_upper'])<F(tail['boundary_reserve'])
    geo=load('FR1_independent_geometry.json');assert geo['status']=='PASS'
    assert all(interval(row['strict_gap'])[0]>0 for row in geo['rows'].values())
    global_checks=load('global_rechecks.json');assert interval(global_checks['energy_at_exact_g_3_71835'])[1]<0
    paths=[p for p in HERE.iterdir() if p.suffix in ['.py','.md','.json'] and p.name!='consolidated_verification.json']
    summary=dict(status='PASS',date='2026-09-19',tests=119,test_files=29,
                 additional_integration_controls=7,integration_control_log=controls.stderr,
                 source_archives_unchanged=True,FR223_direct_masks=99,FR223_repaired_masks=7,
                 FR223_repair_leaves=18,FR223_total_domains=106,FR42_leaves=637,
                 quartic_local_roots=60,full_local_roots=60,global_4D_theorem=False,
                 claim_register='ACCEPTED_RESULTS.md',
                 check_type='Integrity/structure check of freshly recomputed evidence; full scientific replay commands are separate.',
                 evidence_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)})
    if '--record' in sys.argv:
        (HERE/'consolidated_verification.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['evidence_sha256','integration_control_log']},indent=2))

if __name__=='__main__':main()
