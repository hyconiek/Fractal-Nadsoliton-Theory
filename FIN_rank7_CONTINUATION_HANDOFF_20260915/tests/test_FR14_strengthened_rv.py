from fractions import Fraction as F
import frontier_local_boxes as flb

def test_fr14_passes():
    R=flb.raw_box(F(1,6400),F(1,8192),F(1,4600),F(1,100000))
    assert R['status']=='INTERVAL_CERTIFIED'

def test_fr14_negative_control():
    R=flb.raw_box(F(1,6300),F(1,8192),F(1,4600),F(1,100000))
    assert R['status']=='FAILED'
