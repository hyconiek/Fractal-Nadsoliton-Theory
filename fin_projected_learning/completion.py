"""Final batch of the second thirty-round goal: inverse and compatible laws.

No PDF/plotting. The proofs and source limitations are in REPORT.md.
"""
from fractions import Fraction as F
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

import research as r


def lap(K):return np.diag(K.sum(axis=1))-K


def state_distance(rho):
    p=np.diag(rho).real
    return (p[:,None]+p[None,:])/2-rho.real


def cyclic_average(M):
    n=len(M)
    return sum(np.roll(np.roll(M,j,axis=0),j,axis=1) for j in range(n))/n


def regular_projection(G):
    n=len(G);G=r.off((G+G.T)/2)
    degree=G.sum(axis=1)
    multipliers=(degree-np.sum(degree)/(2*n-2))/(n-2)
    return r.off(G-multipliers[:,None]-multipliers[None,:])


def precision_loss(L,C):
    for matrix in [L,C]:
        if not np.allclose(matrix,matrix.T,rtol=0,atol=1e-13):
            raise ValueError('Real symmetric positive-definite matrices are required')
        try:np.linalg.cholesky(matrix)
        except np.linalg.LinAlgError as exc:raise ValueError('Positive definiteness is required') from exc
    sign,logdet=np.linalg.slogdet(L)
    if sign<=0:raise ValueError('Positive definiteness is required')
    return float(np.trace(L@C)-logdet)


def fast_reduced_rhs(rho,gamma):
    K=r.off(rho.real)/gamma
    return -1j*(K@rho-rho@K)


def run():
    K=r.strict();n=12;gamma=.05;eta=.2;s=float(K[0].sum())
    rho=np.eye(n)/n+gamma*K
    seeded=[]
    for t in [0.,1.,10.,100.]:
        learned=(1-math.exp(-eta*gamma*t))*K
        dK,dr=r.projected_rhs(learned,rho,eta,gamma)
        expected=eta*gamma*math.exp(-eta*gamma*t)*K
        assert np.linalg.norm(dK-expected)<1e-13 and np.linalg.norm(dr)<1e-13
        seeded.append(dict(t=t,kernel_scale=1-math.exp(-eta*gamma*t),state_derivative=float(np.linalg.norm(dr))))
    psi=r.pure_average_witness(K,gamma);pure=np.outer(psi,psi.conj())
    cyclic_seed=[]
    for t in [0.,1.,10.]:
        scale=1-math.exp(-eta*gamma*t)
        integrated=t-scale/(eta*gamma)
        U=expm(-1j*integrated*K);state=U@pure@U.conj().T
        mismatch=float(np.linalg.norm(r.off(cyclic_average(state.real))-gamma*K))
        assert mismatch<1e-13
        cyclic_seed.append(dict(t=t,cyclic_source_error=mismatch,
            raw_source_error=float(np.linalg.norm(r.off(state.real)-gamma*K))))
    q=state_distance(rho)
    centered=np.eye(n)-np.ones((n,n))/n
    assert np.linalg.eigvalsh(centered@q@centered)[-1]<1e-12
    v=np.sqrt(2/n)*np.cos(2*np.pi*np.arange(n)/n)
    incompatible=float(v@K@v)
    assert incompatible>0
    # A viable interior step of fixed-degree positive projected learning.
    psi=np.zeros(n);psi[0]=1/math.sqrt(2);psi[1]=-1/math.sqrt(2)
    state=np.outer(psi,psi).astype(complex);A=lap(K)
    U=expm(-.03j*A);rotated=U@state@U.conj().T
    before=r.lyapunov(K,state,gamma);after_rotation=r.lyapunov(K,rotated,gamma)
    h=.001;direction=regular_projection(r.off(rotated.real)-gamma*K)
    repaired=K+h*direction
    assert repaired[~np.eye(n,dtype=bool)].min()>0
    assert np.max(abs(repaired.sum(axis=1)-s))<1e-13
    after=r.lyapunov(repaired,rotated,gamma)
    bound=-(1/h-gamma/2)*np.linalg.norm(repaired-K)**2
    assert abs(before-after_rotation)<1e-13 and after-after_rotation<=bound+1e-13
    # Restore a full Green matrix before taking a positive precision inverse.
    lifted=100*np.eye(n)+K;precision=np.linalg.inv(lifted)
    mass2=1/(100+s);parent=precision-mass2*np.eye(n)
    offprecision=r.off(precision)
    assert offprecision[~np.eye(n,dtype=bool)].max()<0
    assert np.max(abs(parent.sum(axis=1)))<1e-14
    green_error=float(np.linalg.norm(np.linalg.inv(parent+mass2*np.eye(n))-lifted))
    assert green_error<1e-10
    # Fixed-covariance precision descent; the quoted positive sign is ascent.
    C=np.array([[2.,.3],[.3,1.]])
    L0=np.array([[1.4,.2],[.2,.9]])
    gradient=C-np.linalg.inv(L0)
    ascent_derivative=float(np.sum(gradient*gradient))
    sol=solve_ivp(lambda t,y:(np.linalg.inv(y.reshape(2,2))-C).ravel(),
                  (0,30),L0.ravel(),rtol=1e-11,atol=1e-13)
    if not sol.success:raise RuntimeError(sol.message)
    final=sol.y[:,-1].reshape(2,2)
    precision_error=float(np.linalg.norm(final-np.linalg.inv(C)))
    assert precision_error<1e-6
    # A genuine joint covariance/precision functional is flat on reciprocal pairs.
    joint=[]
    for scale in [.4,1.,3.]:
        L=scale*L0;cov=np.linalg.inv(L)
        value=np.trace(L@cov)-np.linalg.slogdet(L)[1]-np.linalg.slogdet(cov)[1]-2
        joint.append(float(value));assert abs(value)<1e-13
    # Fast-learning limit, initialized on the algebraic response manifold.
    rng=np.random.default_rng(8648)
    Z=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3));rho0=Z@Z.conj().T;rho0/=np.trace(rho0)
    g=.5;K0=r.off(rho0.real)/g;T=.2
    reduced=solve_ivp(lambda t,y:fast_reduced_rhs(y.reshape(3,3),g).ravel(),
                      (0,T),rho0.ravel(),rtol=1e-11,atol=1e-13)
    if not reduced.success:raise RuntimeError(reduced.message)
    reduced_state=reduced.y[:,-1].reshape(3,3)
    reduced_energy=lambda R:np.linalg.norm(r.off(R.real))**2/(2*g)
    energy_error=abs(reduced_energy(reduced_state)-reduced_energy(rho0))
    fast=[]
    M=max(np.linalg.norm(K0),1/g)
    for rate in [20.,100.,500.]:
        initial=np.concatenate((K0.ravel().astype(complex),rho0.ravel()))
        def rhs(t,y):
            Kt=y[:9].reshape(3,3).real;Rt=y[9:].reshape(3,3)
            dK,dR=r.projected_rhs(Kt,Rt,rate,g)
            return np.concatenate((dK.ravel().astype(complex),dR.ravel()))
        exact=solve_ivp(rhs,(0,T),initial,rtol=1e-10,atol=1e-12)
        if not exact.success:raise RuntimeError(exact.message)
        Kt=exact.y[:9,-1].reshape(3,3).real;Rt=exact.y[9:,-1].reshape(3,3)
        lag=float(np.linalg.norm(Kt-r.off(Rt.real)/g));lag_bound=2*M/(rate*g*g)
        loss=float(r.lyapunov(K0,rho0,g)-r.lyapunov(Kt,Rt,g))
        loss_bound=4*M*M*T/(rate*g*g)
        assert lag<=lag_bound and -1e-12<=loss<=loss_bound
        fast.append(dict(eta=rate,lag=lag,lag_bound=lag_bound,dissipation=loss,
            dissipation_bound=loss_bound,reduced_state_error=float(np.linalg.norm(Rt-reduced_state))))
    # Normal ordering is not covariant under an arbitrary spectral basis rotation.
    H=np.array([[1.,1.],[1.,-1.]])/math.sqrt(2);E=np.diag([1.,0.])
    covariance_defect=float(np.linalg.norm(r.off(H@E@H.T)-H@r.off(E)@H.T))
    assert covariance_defect>.5
    gibbs=[]
    for gsmall in [.01,.003,.001]:
        beta=n*gsmall;state=expm(beta*K);state/=np.trace(state)
        residual=float(np.linalg.norm(r.off(state)-gsmall*K))
        assert residual>0
        gibbs.append(dict(gamma=gsmall,beta=beta,stationary_residual=residual))
    spectrum=r.certify_strict_spectrum()
    return dict(completed_computational_rounds=list(range(21,30)),
        encoded_state_learning=seeded,
        pure_cyclically_averaged_escape=dict(checks=cyclic_seed,
            scope='Exact trajectory of a DIFFERENT spatially averaged learning rule; spectral occupations remain supplied.'),
        Dirichlet_source=dict(strict_CND_violation=incompatible,
            supplied_distance_max_centered_eigenvalue=float(np.linalg.eigvalsh(centered@q@centered)[-1])),
        regular_positive_repair=dict(state_step_energy_error=abs(before-after_rotation),
            projected_energy_drop=after-after_rotation,upper_bound=bound,
            minimum_weight=float(repaired[~np.eye(n,dtype=bool)].min()),
            row_sum_error=float(np.max(abs(repaired.sum(axis=1)-s)))),
        positive_precision_parent=dict(diagonal_lift=100.,minimum_precision_eigenvalue=float(np.linalg.eigvalsh(precision)[0]),
            maximum_offdiagonal_precision=float(offprecision[~np.eye(n,dtype=bool)].max()),
            screened_mass_squared=mass2,parent_row_sum_error=float(np.max(abs(parent.sum(axis=1)))),
            exact_green_reconstruction_error=green_error),
        precision_flow=dict(ascent_F_derivative=ascent_derivative,scalar_ascent_boundary_time=math.log(2)-.5,
            corrected_descent_error=precision_error,initial_loss=precision_loss(L0,C),final_loss=precision_loss(final,C)),
        joint_bootstrap=dict(reciprocal_pair_functional_values=joint,
            chain_rule_counterexample='C(L)=L^-1 makes the frozen-covariance gradient zero, but total gradient of n-logdetL is -L^-1.'),
        Gibbs_obstruction=dict(distinct_strict_eigenvalues=len(spectrum['eigenvalue_intervals']),
            maximum_scalar_intersections=2,near_zero_gamma_examples=gibbs,
            scope='Exact no-go at positive gamma; no uniform residual lower bound as gamma tends to zero.'),
        fast_learning=dict(reduced_energy_error=energy_error,checks=fast),
        basis_covariance=dict(Hadamard_normal_ordering_defect=covariance_defect))


if __name__=='__main__':
    result=run();Path(__file__).with_name('completion_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
