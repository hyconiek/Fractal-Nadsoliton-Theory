import importlib.util, pathlib, numpy as np
P=pathlib.Path(__file__).with_name('safe_seed.py')
s=importlib.util.spec_from_file_location('safe_seed',P); m=importlib.util.module_from_spec(s); s.loader.exec_module(m)
_,_,_,_,A=m.build_rank7()
def test_A():
    assert np.max(abs(A-A.T))<1e-12
    assert np.max(abs(A.sum(1)))<1e-12
    assert np.linalg.matrix_rank(A,tol=1e-10)==7
    assert np.ptp(np.diag(A))<1e-12
    assert np.linalg.eigvalsh(A).min()>-1e-12
def test_generators_N2():
    st=m.compositions(2); pi=m.count_pi(st,m.G_FROZEN,A)
    for kin in ['heat_bath','metropolis','barker']:
        q=m.count_generator(st,m.G_FROZEN,A,kin)
        assert m.generator_checks(st,q,pi)['passed']
def test_all_C12_sectors_at_degenerate_g0():
    st=m.compositions(2); pi=m.count_pi(st,0,A); q=m.count_generator(st,0,A)
    _,secs=m.sector_spectra(st,q,pi)
    assert set(map(int,secs.keys()))==set(range(12))
    assert all(v['dimension']>0 for v in secs.values())
    assert max(v['max_eigenpair_residual'] for v in secs.values())<1e-10
def test_rho3_fixture():
    st=m.compositions(3); pi=m.count_pi(st,m.G_FROZEN,A); q=m.count_generator(st,m.G_FROZEN,A)
    _,secs=m.sector_spectra(st,q,pi)
    rho=-max(secs['4']['eigenvalues'])
    assert abs(rho-0.13143978619564534)<1e-12
