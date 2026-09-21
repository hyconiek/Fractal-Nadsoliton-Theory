from fractions import Fraction as F
from src.frontier_j6_tail import build_record, shifted_boundary_cover, diameter_bound, local_box_checks

def test_FR5_global_tail_certificate():
    r=build_record(False)
    assert r['status']=='INTERVAL_CERTIFIED_REPLAYED'
    assert r['boundary_shift']['cover']['unresolved_count']==0
    assert r['boundary_shift']['cover']['terminal_leaves']==508
    assert r['boundary_shift']['cover']['reason_counts']=={'SAFE_A_NONPOS':356,'SAFE_B_NONNEG':150,'LOCAL_FR13':2}
    assert F(r['covariance_perturbation']['perturbation_upper']) < F(1,10**6)

def test_FR5_feature_diameter_and_local_quarantine():
    d=diameter_bound(); assert F(d['diameter2_upper'])<5
    q=local_box_checks(); assert q['local_box_inside_FR13_projection']; assert q['odd_mass_inside_FR13_e_radius']
