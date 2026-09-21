from fractions import Fraction as F
import frontier_local_boxes as flb

def test_fr12_passes():
    R=flb.raw_box(F(1,4096),F(1,8192),F(1,8192),F(1,1024))
    assert R['status']=='INTERVAL_CERTIFIED'
    assert R['boundary_ok'] and R['endpoint_ok']

def test_fr12_negative_control():
    R=flb.raw_box(F(1,4096),F(1,8192),F(1,8192),F(1,850))
    assert R['status']=='FAILED'
    assert not R['boundary_ok']
