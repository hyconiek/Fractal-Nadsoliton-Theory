"""Rounds 11--19: temporal source tests and positive clock rigidity.

Only numerical replay here. General proofs, scopes and removal examples are
in REPORT.md. No plotting or document compilation.
"""
from fractions import Fraction as F
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from research import strict,lap,noise_curve


def telegraph(Q,C,kappa):
    n=len(Q);I=np.eye(n)
    return np.block([[Q+C-kappa*I,kappa*I],[kappa*I,Q-C-kappa*I]])


def collapse_matrices(n):
    return np.hstack((np.eye(n),np.eye(n))),np.vstack((np.eye(n),np.eye(n)))/2


def observed_semigroup(Q,C,kappa,t):
    collapse,lift=collapse_matrices(len(Q))
    return collapse@expm(t*telegraph(Q,C,kappa))@lift


def bernstein(lam,drift,atoms):
    return drift*lam+sum(rate*(-math.expm1(-lam*time)) for time,rate in atoms)


def shape_gap(lam1,lam2,time):
    # Inputs used here are well separated; exact positivity is analytic.
    return -math.expm1(-lam1*time)/lam1+math.expm1(-lam2*time)/lam2


def run():
    W=strict();Q=-lap(W);s,R,D=noise_curve(W);C=s*D;n=len(W)
    kappa=.7;L=telegraph(Q,C,kappa);collapse,lift=collapse_matrices(n)
    off=L.copy();np.fill_diagonal(off,0)
    assert off.min()>=0 and np.max(abs(L.sum(axis=1)))<1e-13
    first=collapse@L@lift
    second=collapse@L@L@lift
    first_error=float(np.linalg.norm(first-Q))
    second_error=float(np.linalg.norm(second-Q@Q-C@C))
    assert first_error<1e-13 and second_error<1e-13
    semigroup_defect=float(np.linalg.norm(observed_semigroup(Q,C,kappa,.5)-
        observed_semigroup(Q,C,kappa,.2)@observed_semigroup(Q,C,kappa,.3)))
    assert semigroup_defect>1e-5
    z=.8
    direct=collapse@np.linalg.inv(z*np.eye(2*n)-L)@lift
    reduced=np.linalg.inv(z*np.eye(n)-Q-C@np.linalg.inv((z+2*kappa)*np.eye(n)-Q)@C)
    schur_error=float(np.linalg.norm(direct-reduced))
    assert schur_error<1e-13
    normC=float(np.linalg.norm(C,2));t=.7
    fast=[]
    for rate in [1.,5.,20.,100.]:
        error=float(np.linalg.norm(observed_semigroup(Q,C,rate,t)-expm(t*Q),2))
        upper=math.expm1(normC**2*t/(2*rate))
        assert error<=upper+1e-13
        fast.append(dict(kappa=rate,operator_error=error,analytic_upper_bound=upper))
    eigen=np.linalg.eigvalsh(-Q);positive=eigen[eigen>1e-10]
    lo,hi=float(positive[0]),float(positive[-1])
    assert hi-lo>.1
    atom_time=.4;intensity=.3;drift=.8
    gap=(bernstein(lo,drift,[(atom_time,intensity)])/lo-
         bernstein(hi,drift,[(atom_time,intensity)])/hi)
    predicted=intensity*shape_gap(lo,hi,atom_time)
    assert gap>0 and abs(gap-predicted)<1e-14
    cutoff=.2
    tail_upper=gap/shape_gap(lo,hi,cutoff)
    assert tail_upper>=intensity
    # One nonzero spectral value does not force deterministic time.
    single_lambda=1.;atom_rate=.5;atom_tau=1.
    adjusted_drift=1-atom_rate*(-math.expm1(-1.))
    single=bernstein(1.,adjusted_drift,[(1.,atom_rate)])
    pair=bernstein(2.,adjusted_drift,[(1.,atom_rate)])
    assert abs(single-1)<1e-15 and pair<2
    tiny=[]
    for tau in [.1,.01,.001]:
        rho=1/tau
        f1=bernstein(lo,0.,[(tau,rho)])
        f2=bernstein(hi,0.,[(tau,rho)])
        tiny.append(dict(clock_step=tau,jump_rate=rho,shape_defect=f1/lo-f2/hi,
                         first_mode_error=abs(f1-lo)))
    minR=float(np.linalg.eigvalsh(R)[0])
    commutator=float(np.linalg.norm((-Q)@D-D@(-Q)))
    assert minR<0 and commutator<1e-13
    return dict(completed_rounds=list(range(11,20)),programs=[f'ST{i}' for i in range(8601,8610)],
        stationary_law_does_not_fix_time=dict(variance='1/4',contraction='3/5',
            discrete_AR_covariances=[str(F(1,4)*F(3,5)**k) for k in range(6)],
            discrete_iid_covariances=['1/4']+['0']*5,
            Poisson_AR_centered_linear_rate='2/5',Poisson_AR_centered_square_rate='16/25',
            matched_reset_rates=['2/5','2/5']),
        telegraph=dict(kappa=kappa,first_derivative_error=first_error,
            second_derivative_identity_error=second_error,
            extra_second_derivative_norm=float(np.linalg.norm(C@C)),
            semigroup_defect=semigroup_defect,schur_resolvent_error=schur_error,
            fast_switching=fast),
        moment_extremes=dict(zero_variance='delta_0',unit_variance='(delta_-1+delta_1)/2',
            flat_polynomial='u^2-1/4',flat_moment_identity=str(F(1,16)-F(1,2)*F(1,4)+F(1,16))),
        subordination=dict(low_nonzero_eigenvalue=lo,high_eigenvalue=hi,
            positive_atom_shape_gap=gap,integral_prediction=predicted,
            cutoff=cutoff,atom_tail_upper_bound=tail_upper,
            one_mode_counterexample=dict(drift=adjusted_drift,f1=single,f2=pair),
            small_jump_ambiguity=tiny),
        commuting_is_not_heat_mixture=dict(min_eigenvalue_R0=minR,trace_R0=float(np.trace(R)),
            A_D_commutator=commutator))


if __name__=='__main__':
    result=run()
    Path(__file__).with_name('temporal_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
