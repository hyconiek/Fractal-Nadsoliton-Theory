from src.frontier_fr16 import fr16_record

def test_FR16_wall_adapted_rv():
    r=fr16_record(); assert r['status']=='INTERVAL_CERTIFIED_REPLAYED'
    assert r['strict_checks']['boundary_schur_pass']; assert r['strict_checks']['P1_zero_endpoint_schur_pass']
    assert r['negative_control']['rx=1/4800']['status']=='FAILED'
