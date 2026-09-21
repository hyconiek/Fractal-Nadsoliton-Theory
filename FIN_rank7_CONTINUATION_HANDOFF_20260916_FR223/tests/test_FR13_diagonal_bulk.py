from fractions import Fraction as F
import frontier_local_boxes as flb

def test_fr13_passes():
    R=flb.raw_box(F(1,6200),F(1,6800),F(1,6800),F(1,400))
    assert R['status']=='INTERVAL_CERTIFIED'
    assert R['boundary_ok'] and R['endpoint_ok']

def test_fr13_negative_control():
    R=flb.raw_box(F(1,6200),F(1,6800),F(1,6800),F(1,399))
    assert R['status']=='FAILED'
    assert not R['boundary_ok']
