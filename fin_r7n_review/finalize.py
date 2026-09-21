"""Integrity and acceptance gates for the completed independent R7N replay.

This is a checkpoint validator, not a substitute for the full replay commands.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import os
import platform
import subprocess
import sys

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
def load(name):return json.loads((HERE/name).read_text())

def phase_gate(records,geometry,collars):
    expected={'quartic':27272,'16':54341,'20':25292}
    for k,n in expected.items():
        r=records[k];assert r['complete'] and r['processed']==r['total']==n
        assert not r['failed'] and F(r['minimum_margin'])>0
    for kind in ['quartic','full']:
        assert geometry[kind]['complete_partition']
        assert geometry[kind]['root_ids']==list(range(60))
        assert collars[kind]['count']==60 and len(collars[kind]['roots'])==60
        assert all(F(r['q_upper'])<1 for r in collars[kind]['roots'])
        assert collars[kind]['indices']=={'0':12,'1':24,'2':18,'3':6}
    return {'quartic':60,'full':60}

def main():
    records={k:load('leaf_replay_'+k+'.json') for k in ['quartic','16','20']}
    exact=phase_gate(records,load('phase_geometry.json'),load('collar_replay.json'))
    partial=load('partial_leaf_replay.json');geo=load('partial_geometry.json')
    assert partial['complete'] and not partial['failed'] and partial['processed']==partial['total']==13231
    assert geo['safe_leaves']==13231 and geo['residual_leaves']==5432
    assert 0<F(geo['residual_fraction'])<1
    inv=load('inventory.json');assert not inv['MANIFEST.sha256']['missing'] and not inv['MANIFEST.sha256']['changed']
    assert all(x['matches'] for x in inv['resolved_original_inputs'])
    for name,want in inv['archive_sha256'].items():
        assert hashlib.sha256((ROOT/'FIN_R7N_HANDOFF_20260920'/name).read_bytes()).hexdigest()==want,name
    env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1'
    tests=subprocess.run([sys.executable,'-m','unittest','fin_r7n_review.test_review','-v'],
                         cwd=ROOT,env=env,text=True,capture_output=True,timeout=60)
    assert tests.returncode==0,tests.stdout+tests.stderr
    assert 'Ran 8 tests' in tests.stderr
    files=[p for p in HERE.iterdir() if p.suffix in ['.py','.json','.md'] and p.name!='verification.json']
    data=dict(status='PASS',date='2026-09-20',exact_phase_counts=exact,
              phase_gradient_cells_recomputed=106905,local_roots_and_collars=120,
              partial_target_P_cells_recomputed=13231,unresolved_target_P_cells=5432,
              target_P_global=False,target_S_global=False,source_archive_unchanged=True,
              original_manifest_entries=740,new_unit_tests=8,test_output=tests.stderr,
              environment={'python':platform.python_version()},
              source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
              verification_scope='Checks completed replay records and archive integrity; full mathematical replay commands remain in README.',
              old_baseline_tests='119+7 retained as baseline evidence, not counted as new tests in this intake.',
              raw_portable_checker='Local 120-second timeout; independent replay completed separately.')
    if '--record' in sys.argv:(HERE/'verification.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({k:v for k,v in data.items() if k not in ['source_sha256','test_output']},indent=2))

if __name__=='__main__':main()
