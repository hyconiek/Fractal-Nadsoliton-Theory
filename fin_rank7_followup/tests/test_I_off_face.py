import math, sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'src'))
import off_face as o


def test_R7P065_eta_and_schur_equivalence():
    d=o.strict_target_data(); assert float(eval(d['eta_global_lower_interval'][0]))>0
    rng=np.random.default_rng(65067)
    _,sigma,_,_,C=o.constants()
    for _ in range(40):
        J=rng.uniform(0,5,4); p=o.p_from_fields(J)
        q=p[::2].sum()
        def stats(F,pp):
            m=pp@F; X=F-m; return m,X.T@(pp[:,None]*X)
        mp,Cp=stats(C[::2],p[::2]/q); mm,Cm=stats(C[1::2],p[1::2]/(1-q))
        W=q*Cp+(1-q)*Cm; b=np.sqrt(q*(1-q))*(mp-mm)
        M=W+np.outer(b,b)
        eta=1-b[3]**2/sigma; assert eta>0
        Mt=W[:3,:3]+np.outer(b[:3],b[:3])/eta
        assert sum(np.linalg.eigvalsh(M)>sigma+1e-11)==sum(np.linalg.eigvalsh(Mt)>sigma+1e-11)


def test_R7P067_compactification_matches_finite_fields():
    rng=np.random.default_rng(67067)
    for _ in range(100):
        J=rng.uniform(0,12,4); x=np.exp(-J)
        p1=o.p_from_fields(J); p2=o.p_from_compact(x)
        np.testing.assert_allclose(p1,p2,rtol=0,atol=3e-14)
    t=o.compactification_theorem(); assert t['exponent_rows'][0]==['0','0','0','0']
    # Irrational exponents must remain explicit, not split into independent variables.
    assert any('sqrt(3)' in a for row in t['exponent_rows'] for a in row)


def test_compact_boundary_anchor_and_known_equality():
    _,sigma,_,_,_=o.constants(); L=o.L_global()
    Jstar=math.atanh(math.sqrt(1-sigma/(L[3]/6))); xstar=math.exp(-Jstar)
    p=o.p_from_compact([0,0,0,0]); assert np.isclose(p.sum(),1) and p[0]==1
    M,_=o.covariance_direct_compact([xstar,1,1,0]); eig=np.linalg.eigvalsh(M)
    assert abs(eig[-2]-sigma)<2e-13 and abs(eig[-1]-sigma)<2e-13


def test_independent_field_and_compact_covariances():
    for J in ([.3,.7,1.1,2.0],[2.2,.1,3.0,5.0],[.01,.02,.03,.04]):
        M1,p1=o.covariance_direct_fields(J); M2,p2=o.covariance_direct_compact(np.exp(-np.array(J)))
        np.testing.assert_allclose(M1,M2,rtol=0,atol=4e-14); np.testing.assert_allclose(p1,p2,rtol=0,atol=4e-14)

def test_reoptimized_parity_slope_is_interval_certified_negative():
    c=o.reoptimized_parity_slope_certificate()
    lo,hi=c['display']
    assert lo < hi < -0.13128285839
    assert lo > -0.13128285841


def test_aligned_compact_chart_matches_boundary_variables():
    c=o.aligned_compact_chart()
    assert c['boundary_match'].startswith('At y=0')
    assert 'sqrt(3)' in ''.join(c['odd_aggregate_weights'])

def test_aligned_compact_probabilities_match_finite_fields():
    rng=np.random.default_rng(67068)
    for _ in range(100):
        J=rng.uniform(0,8,4)
        r,s,t,y=np.exp([-2*J[0],-1.5*J[1],-.5*J[2],-2*J[3]])
        np.testing.assert_allclose(o.p_from_aligned_compact(r,s,t,y),o.p_from_fields(J),rtol=0,atol=4e-14)
