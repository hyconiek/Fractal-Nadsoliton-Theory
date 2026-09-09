"""A state-programmed CPTP simulation, not the original self-feedback ODE.

Fresh program copies C, interaction schedule, and overall rate are supplied.
The target W is encoded in C=I/n+gamma W, not generated from no input.
"""
from fractions import Fraction as F
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

try:
    from . import research as r
except ImportError:
    import research as r


def unitary(n,tau):
    V,S,D,P=r.structures(n)
    plus=(np.eye(n*n)+S)/2-D
    base=np.exp(-.5j*tau)
    return base*np.eye(n*n)+(np.exp(.5j*tau)-base)*plus+(np.exp(.5j*(n-1)*tau)-base)*P


def channel(C,tau):
    n=len(C);U=unitary(n,tau).reshape(n,n,n,n)
    return np.einsum('iapb,bc,jaqc->ijpq',U,C,U.conj(),optimize=True).reshape(n*n,n*n)


def choi(superoperator,n):
    return superoperator.reshape(n,n,n,n).transpose(0,2,1,3).reshape(n*n,n*n)


def variance_operator(C):
    n=len(C);h=r.source(C)
    return (np.eye(n)+(n-2)*C.T)/4-h@h


def noise(C,sigma):
    n=len(C);V,*_=r.structures(n);h=r.source(C)
    centered=V-np.kron(h,np.eye(n))
    B=variance_operator(C)
    return r.marginal(centered@np.kron(sigma,C)@centered,n)-(B@sigma+sigma@B)/2


def perron_one_step_loss(n,c0,tau):
    """Exact, for any program with <uniform|program|uniform>=c0."""
    return ((1-c0)*(1-4/n**2)*np.sin(tau/2)**2
            +c0*4*(n-1)/n**2*np.sin((n-2)*tau/4)**2)


def perron_leakage_coefficient(n,c0):
    return (n*n-4)/(4*n*n)+(n-2)*(n-4)/(4*n)*c0


def loading_data(W):
    """Optimum of the declared uniform-variance bound, not all algorithms."""
    n=len(W)
    if n!=12:raise ValueError('The simplified Perron formula is certified here for n=12')
    lam=np.linalg.eigvalsh(W);ell=-lam[0];s=lam[-1]
    gamma=1/(n*ell);C=np.eye(n)/n+gamma*W
    variance_cost=66*ell**2+30*ell*s-s*s
    leakage_cost=55*ell**2+20*ell*s
    return dict(maximum_feasible_loading=float(gamma),minimum_reference_eigenvalue=float(lam[0]),
        perron_eigenvalue=float(s),optimal_variance_cost=float(variance_cost),
        optimal_perron_leakage_cost=float(leakage_cost),
        check_variance_cost=float(np.linalg.eigvalsh(variance_operator(C))[-1]/gamma**2),
        scope='Optimizes program loading/variance bounds for the specified canonical processor; not universal sample complexity.')


def loading_certificate():
    import importlib.util
    root=Path(__file__).resolve().parents[1]
    spec=importlib.util.spec_from_file_location('prior_learning',root/'fin_projected_learning/research.py')
    previous=importlib.util.module_from_spec(spec);spec.loader.exec_module(previous)
    data=previous.certify_strict_spectrum()
    intervals=[tuple(map(F,pair)) for pair in data['eigenvalue_intervals']]
    minimum=min(range(7),key=lambda k:intervals[k][0])
    assert all(intervals[minimum][1]<intervals[j][0] for j in range(7) if j!=minimum)
    assert all(intervals[0][0]>intervals[j][1] for j in range(1,7))
    ell=(-intervals[minimum][1],-intervals[minimum][0]);s=intervals[0]
    gamma=(1/(12*ell[1]),1/(12*ell[0]))
    variance=(66*ell[0]**2+30*ell[0]*s[0]-s[1]**2,
              66*ell[1]**2+30*ell[1]*s[1]-s[0]**2)
    leakage=(55*ell[0]**2+20*ell[0]*s[0],55*ell[1]**2+20*ell[1]*s[1])
    assert variance[0]>0 and leakage[0]>0
    return dict(minimum_mode=minimum,perron_mode=0,
        loading_interval=[str(x) for x in gamma],
        optimal_variance_cost_interval=[str(x) for x in variance],
        optimal_leakage_cost_interval=[str(x) for x in leakage],
        source='Recomputed outward rational strict Fourier eigenvalue enclosures; analytic optimization in PROOF.md.')


def certified_strict_bound(copies,time=F(1)):
    if not isinstance(copies,int) or copies<1:raise ValueError('Positive integer copy budget required')
    time=F(time)
    if time<0:raise ValueError('Positive time required')
    # n=12, gamma=1/20, ||W||<=||W||F<3, ||T(C)||<3/20,
    # ||V||=11/2, 0<=B_C<=5I/6. All scalar arithmetic is rational.
    bound=F(2000,3)*time**2/copies+F(5324108,3)*time**3/copies**2
    return min(F(2),bound)


def run():
    W=r.strict();n=12;gamma=.05;g=1/gamma;C=np.eye(n)/n+gamma*W
    target=expm(1j*W);S_target=np.kron(target,target.conj())
    J_target=choi(S_target,n)/n
    B=variance_operator(C)
    V,*_=r.structures(n)
    direct_B=r.marginal(V@V@np.kron(np.eye(n),C),n)-r.source(C)@r.source(C)
    assert np.linalg.norm(B-direct_B)<1e-13 and np.linalg.eigvalsh(B)[0]>0
    assert np.linalg.norm(noise(C,np.eye(n)))<1e-12
    rho0=np.zeros((n,n));rho0[0,0]=1
    target_rho=target@rho0@target.conj().T
    perron=np.ones(n)/np.sqrt(n);P_perron=np.outer(perron,perron)
    c0=float(perron@C@perron)
    asymptotic_perron=g*g*perron_leakage_coefficient(n,c0)
    records=[]
    for N in [100,1000,10000,100000,1000000]:
        one=channel(C,g/N);total=np.linalg.matrix_power(one,N)
        J=choi(total,n)/n
        J=(J+J.conj().T)/2
        delta=J-J_target
        rho=(total@rho0.ravel()).reshape(n,n)
        rho=(rho+rho.conj().T)/2
        perron_output=(total@P_perron.ravel()).reshape(n,n)
        perron_loss=float(1-np.vdot(perron,perron_output@perron).real)
        records.append(dict(copies=N,
            normalized_choi_trace_distance=float(sum(abs(np.linalg.eigvalsh(delta)))/2),
            selected_input_trace_distance=float(sum(abs(np.linalg.eigvalsh(rho-target_rho)))/2),
            trace_preservation_error=float(np.linalg.norm(np.trace(choi(total,n).reshape(n,n,n,n),axis1=0,axis2=2)-np.eye(n))),
            minimum_normalized_choi_eigenvalue=float(np.linalg.eigvalsh(J)[0]),
            perron_survival_loss=perron_loss,
            perron_loss_times_copies=perron_loss*N,
            rigorous_diamond_error_upper=str(certified_strict_bound(N)),
            rigorous_diamond_error_upper_decimal=float(certified_strict_bound(N))))
    # A small fixed phase verifies the entire Stinespring channel construction.
    tau=.003;U=unitary(n,tau)
    direct=r.marginal(U@np.kron(rho0,C)@U.conj().T,n)
    assert np.linalg.norm((channel(C,tau)@rho0.ravel()).reshape(n,n)-direct)<1e-13
    assert np.linalg.norm(U-expm(1j*tau*V))<1e-13
    assert certified_strict_bound(1000000)<F(7,10000)
    one_perron=(channel(C,tau)@P_perron.ravel()).reshape(n,n)
    assert abs((1-np.vdot(perron,one_perron@perron).real)-perron_one_step_loss(n,c0,tau))<1e-13
    return dict(status='Conditional sample-programmed quantum simulation with an explicit finite-copy error bound.',
        encoding='C=I/12+(1/20)W; every program copy already contains the target covariance.',
        supplied_resources=['Fresh independent program copies','Canonical two-body interaction V','Interaction timing and overall g=20'],
        variance_eigenvalue_range=[float(np.linalg.eigvalsh(B)[0]),float(np.linalg.eigvalsh(B)[-1])],
        exact_one_million_copy_diamond_bound=str(certified_strict_bound(1000000)),
        asymptotic_perron_loss_times_copies=float(asymptotic_perron),
        optimal_loading=loading_data(W),outward_loading_certificate=loading_certificate(),
        records=records,
        scope=['This is not the finite-eta deterministic self-feedback law.',
            'Finite-copy noise is retained; exact finite quantum programming is not claimed.',
            'No non-target-encoded state source or physical calibration is supplied.'])


if __name__=='__main__':
    result=run();Path(__file__).with_name('collision_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
