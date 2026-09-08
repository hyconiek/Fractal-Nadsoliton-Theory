"""New 30-round goal: source-level audit of actual projected learning laws.

Historical files are read, not rewritten. No PDF or plotting.
"""
from fractions import Fraction as F
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm


def strict():
    return np.array([[0. if i==j else math.cos(.18575*min(abs(i-j),12-abs(i-j))+.1625)/
        (1+min(abs(i-j),12-abs(i-j))**1.8) for j in range(12)] for i in range(12)])


def off(M):
    value=np.array(M,copy=True);np.fill_diagonal(value,0)
    return value


def projected_rhs(K,rho,eta,gamma):
    return eta*(off(rho.real)-gamma*K),-1j*(K@rho-rho@K)


def lyapunov(K,rho,gamma):
    return gamma*np.sum(K*K)/2-np.trace(K@rho).real


def spectral_projectors_circulant(n=12):
    basis=np.exp(2j*np.pi*np.outer(np.arange(n),np.arange(n))/n)/math.sqrt(n)
    groups=[[0]]+[[k,n-k] for k in range(1,n//2)]+[[n//2]]
    return [basis[:,idx]@basis[:,idx].conj().T for idx in groups]


def pure_time_average(K,psi):
    """Strict-only pinching; its exact sector multiplicities are certified.

    This must not silently group the spectrum of another changing kernel.
    """
    if K.shape!=(12,12) or not np.array_equal(K,strict()):
        raise ValueError('This averaging routine is certified only for frozen strict K')
    projections=spectral_projectors_circulant(len(K))
    return sum((P@np.outer(psi,psi.conj())@P).real for P in projections)


def pure_average_witness(K,gamma):
    n=len(K);x=np.arange(n);psi=np.zeros(n,dtype=complex)
    for k in [0,n//2]:
        v=np.ones(n)/math.sqrt(n) if k==0 else (-1.)**x/math.sqrt(n)
        weight=1/n+gamma*float(v@K@v)
        if weight<0:raise ValueError('Negative modal weight')
        psi+=math.sqrt(weight)*v
    for k in range(1,n//2):
        co=math.sqrt(2/n)*np.cos(2*np.pi*k*x/n)
        si=math.sqrt(2/n)*np.sin(2*np.pi*k*x/n)
        weight=1/n+gamma*float(co@K@co)
        if weight<0:raise ValueError('Negative modal weight')
        psi+=math.sqrt(weight)*(co+1j*si)
    return psi


def expected_learning(T,eta=F(1,100),gamma=F(1,100)):
    factor=1-eta*gamma
    # Signal phase uniform, independent Gaussian noise variance 1/4.
    variance_outer=F(1,8)+F(1,4)+F(1,16)
    stationary_variance=eta**2*variance_outer/(1-factor**2)
    return factor,variance_outer,stationary_variance


def off_double_commutator_map(K):
    n=len(K);pairs=list(itertools.combinations(range(n),2));columns=[]
    for k in range(n):
        D=np.zeros((n,n));D[k,k]=1
        C=K@(K@D-D@K)-(K@D-D@K)@K
        columns.append([C[i,j] for i,j in pairs])
    return np.array(columns).T


def certify_strict_spectrum():
    """Reuse the exact transcendental enclosures from the earlier audit."""
    import importlib.util
    path=Path(__file__).resolve().parent.parent/'fin_replication_consistency'/'certify.py'
    spec=importlib.util.spec_from_file_location('fin_weight_enclosures',path)
    cert=importlib.util.module_from_spec(spec);spec.loader.exec_module(cert)
    weights=cert.strict_weights();root3=cert.root_interval(3,2)
    halfroot=cert.scale(F(1,2),root3)
    one=(F(1),F(1));half=(F(1,2),F(1,2));zero=(F(0),F(0))
    table=[one,halfroot,half,zero,cert.scale(-1,half),cert.scale(-1,halfroot),
           cert.scale(-1,one),cert.scale(-1,halfroot),cert.scale(-1,half),zero,half,halfroot]
    eigen=[]
    for k in range(7):
        value=cert.scale((-1)**k,weights[5])
        for d in range(1,6):
            value=cert.add(value,cert.scale(2,cert.multiply(weights[d-1],table[(k*d)%12])))
        q=10**14;lo,hi=value
        lower=lo*q;upper=hi*q
        eigen.append((F(lower.numerator//lower.denominator,q),
                      F(-((-upper.numerator)//upper.denominator),q)))
    gaps=[]
    for i,j in itertools.combinations(range(7),2):
        lo1,hi1=eigen[i];lo2,hi2=eigen[j]
        gap=max(lo1-hi2,lo2-hi1)
        assert gap>0
        gaps.append(gap)
    density_lower=min(F(1,12)+F(1,20)*lo for lo,hi in eigen)
    assert density_lower>0
    return dict(eigenvalue_intervals=[[str(lo),str(hi)] for lo,hi in eigen],
                minimum_intersector_gap_lower=str(min(gaps)),
                density_minimum_eigenvalue_lower=str(density_lower),
                assertion='Seven distinct real Fourier eigenvalues; multiplicities 1,2,2,2,2,2,1.')


def replay_actual_source(seed=8621):
    import contextlib
    import importlib.util
    import io
    path=Path(__file__).resolve().parent.parent/'nadsoliton_neural_analysis.py'
    spec=importlib.util.spec_from_file_location('archived_fin_neural_source',path)
    source=importlib.util.module_from_spec(spec);spec.loader.exec_module(source)
    oldstate=np.random.get_state();log=io.StringIO()
    try:
        np.random.seed(seed)
        with contextlib.redirect_stdout(log):
            learned,target=source.run_hebbian_emergence_simulation()
    finally:np.random.set_state(oldstate)
    wrong=source.analyze_network_properties(-target,target)
    n=len(target);teacher=off(np.cos((np.arange(n)[:,None]-np.arange(n)[None,:])*math.pi/4))
    return dict(seed=seed,iterations=30000,stdout=log.getvalue(),
        learned_teacher_correlation=float(np.corrcoef(learned.ravel(),teacher.ravel())[0,1]),
        learned_legacy_target_correlation=float(np.corrcoef(learned.ravel(),target.ravel())[0,1]),
        learned_strict_correlation=float(np.corrcoef(learned.ravel(),strict().ravel())[0,1]),
        negative_target_test=dict(correlation=float(wrong['hebbian_correlation']),
                                 emitted_claim=wrong['conclusion']))


def run():
    K=strict();n=len(K);gamma=.05;eta=.2
    factor,variance_outer,stationary_variance=expected_learning(30000)
    signal=np.cos((np.arange(n)[:,None]-np.arange(n)[None,:])*math.pi/4)
    target=.5*off(signal)
    rho=np.eye(n)/n+gamma*K
    dK,drho=projected_rhs(K,rho,eta,gamma)
    assert np.linalg.eigvalsh(rho)[0]>0
    assert np.linalg.norm(dK)<1e-14 and np.linalg.norm(drho)<1e-14
    psi=pure_average_witness(K,gamma)
    averaged=pure_time_average(K,psi)
    average_error=float(np.linalg.norm(averaged-rho))
    instant_error=float(np.linalg.norm(off(np.outer(psi,psi.conj()).real)-gamma*K))
    assert abs(np.vdot(psi,psi).real-1)<1e-13 and average_error<1e-13 and instant_error>.01
    # Degenerate projection cannot be replaced by a scalar occupation*identity.
    P=np.eye(2);z=np.array([1.,0.]);correct=P@np.outer(z,z)@P
    scalar=np.vdot(z,z).real*P
    degeneracy_error=float(np.linalg.norm(correct-scalar))
    assert degeneracy_error==1
    rng=np.random.default_rng(8621)
    Z=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));state=Z@Z.conj().T;state/=np.trace(state)
    W=off(rng.normal(size=(n,n)));W=(W+W.T)/2
    KW,rW=projected_rhs(W,state,eta,gamma)
    derivative=gamma*np.sum(W*KW)-np.trace(KW@state+W@rW).real
    prediction=-eta*np.linalg.norm(off(state.real)-gamma*W)**2
    assert abs(derivative-prediction)<1e-13
    # Zero-learning motion need not be a stationary state: exact two-level orbit.
    k=.7;g=.2;smallK=np.array([[0.,k],[k,0.]])
    smallrho=np.array([[.7,g*k],[g*k,.3]],complex)
    assert np.linalg.eigvalsh(smallrho)[0]>0
    period=math.pi/k
    orbit=[]
    for time in [0.,period/8,period/4]:
        U=expm(-1j*time*smallK);r=U@smallrho@U.conj().T
        update,motion=projected_rhs(smallK,r,eta,g)
        assert np.linalg.norm(update)<1e-14
        orbit.append(dict(t=time,kernel_derivative=float(np.linalg.norm(update)),
                          state_derivative=float(np.linalg.norm(motion)),population=float(r[0,0].real)))
    rank=np.linalg.matrix_rank(off_double_commutator_map(K),tol=1e-11)
    assert rank==n-1
    spectral_certificate=certify_strict_spectrum()
    actual_source=replay_actual_source()
    return dict(programs=[f'ST{x}' for x in range(8621,8631)],completed_rounds=list(range(1,11)),
        actual_rule=dict(name='projected leaky Hebb, not PCA/Oja',decay_factor=str(factor),
            initial_memory_after_30000=float(factor)**30000,
            outer_product_variance=str(variance_outer),stationary_entry_variance=str(stationary_variance),
            conditional_Lyapunov_exponent=math.log(float(factor)),
            teacher_covariance_offdiagonal=target.tolist()),
        entropy_obstruction=dict(maximum_12_state_entropy=math.log(12),requested_4_bit_target=math.log(16),
            unavoidable_relative_error=1-math.log(12)/math.log(16)),
        degeneracy_scalarization_counterexample_norm=degeneracy_error,
        projected_mixed_fixed_point=dict(gamma=gamma,min_density_eigenvalue=float(np.linalg.eigvalsh(rho)[0]),
            K_derivative=float(np.linalg.norm(dK)),rho_derivative=float(np.linalg.norm(drho)),
            gamma_max=float(-1/(n*np.linalg.eigvalsh(K)[0]))),
        pure_time_average=dict(norm=float(np.vdot(psi,psi).real),average_covariance_error=average_error,
                               instantaneous_kernel_mismatch=instant_error),
        lyapunov_identity=dict(actual_derivative=float(derivative),predicted_derivative=float(prediction)),
        zero_learning_nonstationary_orbit=orbit,
        exact_spectral_certificate=spectral_certificate,
        actual_archived_source_replay=actual_source,
        strict_invariant_set_test=dict(off_double_commutator_diagonal_rank=int(rank),expected_rank=n-1,
            min_positive_edge=float(K[K>0].min()),
            statement='Maximum-principle proof, not sampled rank alone, licenses the dense-positive rigidity result.'))


if __name__=='__main__':
    result=run();Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='actual_rule'},indent=2))
