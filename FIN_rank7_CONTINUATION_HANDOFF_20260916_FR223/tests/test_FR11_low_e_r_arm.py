from fractions import Fraction as F
import frontier_local_boxes as flb

def test_fr11_low_e_r_arm_passes():
    R=flb.raw_box(F(1,7600),F(1,4800),F(1,4800),F(1,81920))
    assert R['status']=='INTERVAL_CERTIFIED'
    assert R['boundary_ok'] and R['endpoint_ok']

def test_fr11_larger_e_is_negative_control():
    R=flb.raw_box(F(1,7600),F(1,4800),F(1,4800),F(1,75000))
    assert R['status']=='FAILED'
    assert not R['boundary_ok']

def test_fr11_larger_r_is_negative_control():
    R=flb.raw_box(F(1,7500),F(1,4800),F(1,4800),F(1,81920))
    assert R['status']=='FAILED'
    assert not R['boundary_ok']
