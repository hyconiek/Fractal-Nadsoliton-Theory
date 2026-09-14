from fractions import Fraction as F
from quantum_bridge_audit import run, DELTA_SMALL, CURRENT_SMALL, D0

def test_provenance_gap_exact():
    assert DELTA_SMALL-CURRENT_SMALL == F(1,2000000000000000)

def test_loading_factor_exact():
    r=run()['R7P-117']
    assert F(r['conditional_floor_factor']) == F(49,50)**2
    assert F(r['conditional_floor_lower']) == F(49,50)**2*D0

def test_no_causal_promotion():
    r=run()
    assert 'No causal' in r['R7P-115']['verdict']
    assert r['R7P-120']['classification'].endswith('NO_CAUSAL_SPECTRAL_BRIDGE')
