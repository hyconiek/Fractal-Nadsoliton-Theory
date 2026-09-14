from fractions import Fraction as F
import json
from pathlib import Path

import off_face_local as loc

ROOT=Path(__file__).resolve().parents[1]


def test_acceptance_radius_is_interval_certified():
    R=loc.raw_local_cone(F(1,8192),F(1,10000))
    assert R['status']=='INTERVAL_CERTIFIED'
    assert all(R['signs'].values())
    assert R['boundary_ok'] and R['endpoint_ok']
    assert R['c2'].v.lo>0


def test_larger_dyadic_radius_is_rejected_by_same_checker():
    R=loc.raw_local_cone(F(1,4096),F(1,10000))
    assert R['status']=='FAILED'
    assert not R['boundary_ok']
    # This is a checker negative-control, not a counterexample to the theorem.


def test_characteristic_anchor_numerically_contains_zero():
    R=loc.raw_local_cone(F(1,8192),F(1,10000))
    # Interval dependency widens the exact identities, but must contain them.
    assert R['P'].v.lo<=0<=R['P'].v.hi
    assert R['P1'].v.lo<=0<=R['P1'].v.hi


def test_write_certificate_and_scope():
    rec=loc.write_certificate()
    assert rec['status']=='INTERVAL_CERTIFIED'
    assert rec['domain']['rho']=='1/8192'
    assert 'does not prove the full off-face orthant ceiling' in rec['nontransfer']
    saved=json.loads((ROOT/'certificates/R7P-068_local_cone.json').read_text())
    assert saved['status']=='INTERVAL_CERTIFIED'
