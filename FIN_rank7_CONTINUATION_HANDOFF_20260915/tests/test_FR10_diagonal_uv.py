from fractions import Fraction as F
import json
from pathlib import Path
import frontier_local_boxes as flb
ROOT=Path(__file__).resolve().parents[1]

def test_fr10_acceptance_box():
    R=flb.raw_box(F(1,8192),F(1,4800),F(1,4800),F(1,3072))
    assert R['status']=='INTERVAL_CERTIFIED'
    assert all(R['signs'].values())
    assert R['boundary_ok'] and R['endpoint_ok']

def test_fr10_larger_e_negative_control():
    R=flb.raw_box(F(1,8192),F(1,4800),F(1,4800),F(1,3052))
    assert R['status']=='FAILED'
    assert not R['boundary_ok']

def test_fr10_larger_uv_negative_control():
    R=flb.raw_box(F(1,8192),F(1,4700),F(1,4700),F(1,3072))
    assert R['status']=='FAILED'
    assert not R['boundary_ok']

def test_fr10_record_scope_and_persistence():
    rec=flb.write_fr10()
    assert rec['status']=='INTERVAL_CERTIFIED'
    assert rec['domain']['e=1-q_even']=='<=1/3072'
    saved=json.loads((ROOT/'results/FR10_diagonal_uv_local_box.json').read_text())
    assert 'No full positive-orthant' in saved['scope']
