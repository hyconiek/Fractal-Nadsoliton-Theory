from off_face_tail import run

def test_tail_is_strictly_certified():
    r=run(); assert r['strict']; assert r['global_4D_ceiling']=='NOT_PROVED'

def test_residual_is_not_hidden():
    assert 'remains unresolved' in run()['residual']
