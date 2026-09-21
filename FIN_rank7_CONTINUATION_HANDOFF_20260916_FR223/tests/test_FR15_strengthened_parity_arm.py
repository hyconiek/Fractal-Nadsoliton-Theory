from src.frontier_fr15 import fr15_record

def test_FR15_strengthened_parity_arm():
    r=fr15_record(); assert r['status']=='INTERVAL_CERTIFIED_REPLAYED'
    assert r['strict_checks']['boundary_schur_pass']; assert r['strict_checks']['P1_zero_endpoint_schur_pass']
    assert r['negative_control']['e=1/13000']['status']=='FAILED'
