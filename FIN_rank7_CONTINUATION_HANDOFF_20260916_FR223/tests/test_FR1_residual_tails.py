from fractions import Fraction as F
from src.frontier_residual_tails import build, a_tail_record, s_tail_record, t_tail_record, q_lower_bound_record


def test_q_lower_bound_exact_record():
    r=q_lower_bound_record()
    assert r['status']=='EXACT_PROVED'
    assert 'q*(1+a)^2 >= 1+a^2' in r['compact_form']


def test_large_J3_tail_strict():
    r=a_tail_record()
    assert r['status']=='INTERVAL_CERTIFIED'
    assert r['a_threshold']=='1/30'
    assert F(r['strict_gap_interval'][0])>0


def test_large_J4_tail_strict():
    r=s_tail_record()
    assert r['status']=='INTERVAL_CERTIFIED'
    assert r['s_threshold']=='1/128'
    assert F(r['strict_gap_interval'][0])>0


def test_large_J5_tail_strict_and_large_improvement():
    r=t_tail_record()
    assert r['status']=='INTERVAL_CERTIFIED'
    assert r['t_threshold']=='1/9'
    assert F(r['strict_gap_interval'][0])>0
    assert F(1,9) > F(1,2**11)


def test_new_residual_hull():
    r=build()['new_residual_outer_box']['necessary_conditions_for_unresolved_point']
    assert len(r)==3
