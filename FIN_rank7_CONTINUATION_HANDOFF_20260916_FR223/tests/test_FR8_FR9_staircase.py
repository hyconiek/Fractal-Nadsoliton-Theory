from fractions import Fraction as F
import frontier_local_boxes as flb


def test_fr8_diagonal_rv_arm():
    R=flb.raw_box(F(1,6500),F(1,8192),F(1,4600),F(1,100000))
    assert R['status']=='INTERVAL_CERTIFIED'
    assert R['boundary_ok'] and R['endpoint_ok']


def test_fr9_diagonal_ru_arm():
    R=flb.raw_box(F(1,6500),F(1,2432),F(1,8192),F(1,100000))
    assert R['status']=='INTERVAL_CERTIFIED'
    assert R['boundary_ok'] and R['endpoint_ok']


def test_fr9_exploratory_negative_control():
    R=flb.raw_box(F(1,6500),F(1,2404),F(1,8192),F(1,100000))
    assert R['status']=='FAILED'
