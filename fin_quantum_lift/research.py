"""New goal after ST8652: conditional microscopic lift of the supplied source.

Tensor replication and a mean-field pair law are explicit additional premises.
This does not identify them as already derived FIN physics.
"""
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm,logm


def strict(n=12):
    if n!=12:raise ValueError('The strict reference here has twelve labels')
    return np.array([[0. if i==j else math.cos(.18575*min(abs(i-j),12-abs(i-j))+.1625)/
        (1+min(abs(i-j),12-abs(i-j))**1.8) for j in range(n)] for i in range(n)])


def legacy_cycle():
    """Same intermediate legacy reference as fin_projected_learning/geometry.py.

    Signed weights are Hamiltonian data here, not positive transition rates.
    """
    return np.array([[0. if i==j else
        4*math.log(2)*math.cos(math.pi*min(abs(i-j),12-abs(i-j))/4+math.pi/6)/
        (1+.01*min(abs(i-j),12-abs(i-j))) for j in range(12)] for i in range(12)])


def source(rho):
    A=rho.real.copy();np.fill_diagonal(A,0);return A


def structures(n):
    swap=np.zeros((n*n,n*n));D=np.zeros_like(swap);omega=np.zeros(n*n)
    for i in range(n):
        D[i*n+i,i*n+i]=1;omega[i*n+i]=1
        for j in range(n):swap[i*n+j,j*n+i]=1
    P=np.outer(omega,omega)/n
    V=(n*P+swap-2*D)/2
    return V,swap,D,P


def lift_from_observables(n):
    V=np.zeros((n*n,n*n))
    for i in range(n):
        for j in range(i+1,n):
            X=np.zeros((n,n));X[i,j]=X[j,i]=1/math.sqrt(2)
            V+=np.kron(X,X)
    return V


def marginal(R,n):
    return np.trace(R.reshape(n,n,n,n),axis1=1,axis2=3)


def partial_transpose(R,n):
    return R.reshape(n,n,n,n).transpose(0,3,2,1).reshape(n*n,n*n)


def product_evolution(C,t):
    """t is the pair phase g*time, not an absolute physical clock."""
    n=len(C);V,S,D,P=structures(n)
    U=expm(1j*t*V)
    R=U@np.kron(C,C)@U.conj().T
    return marginal(R,n),R


def product_closed(C,t):
    """Exact for real symmetric C with uniform diagonal and diag(C²)."""
    n=len(C);I=np.eye(n)
    if not (np.allclose(C,C.T) and np.max(abs(C.imag))<1e-14
            and np.allclose(np.diag(C),np.ones(n)/n)
            and np.allclose(np.diag(C@C),np.ones(n)*np.trace(C@C)/n)):
        raise ValueError('The declared constant-diagonal/constant-square-diagonal class is required')
    A=2*(math.cos(t)-1)/n
    B=2*(math.cos((n/2-1)*t)-math.cos(t))/n
    return C+A*(C-I/n)+B*(C@C-np.trace(C@C)*I/n)


def product_acceleration(C):
    n=len(C)
    return -2/n*(C-np.eye(n)/n)+(2-n/2)*(C@C-np.trace(C@C)*np.eye(n)/n)


def antisymmetric_completion(C):
    """Explicit target-dependent stationary two-body completion, not a broadcaster."""
    n=len(C)
    if n<3:raise ValueError('This affine construction requires n>=3')
    V,S,D,P=structures(n);A=(np.eye(n*n)-S)/2
    R=2/(n-2)*A@(np.kron(C,np.eye(n))+np.kron(np.eye(n),C)-np.eye(n*n)/(n-1))@A
    return R


def entropy(C):
    p=np.linalg.eigvalsh(C);return float(-sum(x*math.log(x) for x in p if x>1e-14))


def local_collective_commutant_rank(n):
    """Numerical rank only; the proof uses the isolated Omega eigenspace."""
    V,*_=structures(n);I=np.eye(n);basis=[]
    for i in range(n):
        H=np.zeros((n,n),complex);H[i,i]=1;basis.append(H)
        for j in range(i+1,n):
            H=np.zeros((n,n),complex);H[i,j]=H[j,i]=1;basis.append(H)
            H=np.zeros((n,n),complex);H[i,j]=1j;H[j,i]=-1j;basis.append(H)
    columns=[]
    for H in basis:
        L=np.kron(H,I)+np.kron(I,H);A=V@L-L@V
        columns.append(np.r_[A.real.ravel(),A.imag.ravel()])
    return int(np.linalg.matrix_rank(np.array(columns).T,tol=1e-10))


def run():
    W=strict();n=12;gamma=.05;C=np.eye(n)/n+gamma*W
    V,S,D,P=structures(n);Rprod=np.kron(C,C)
    exact_lift_error=float(np.linalg.norm(V-lift_from_observables(n)))
    contraction=marginal(V@np.kron(np.eye(n),C),n)
    assert np.linalg.norm(contraction-source(C))<1e-13
    hartree_residual=float(np.linalg.norm(source(C)@C-C@source(C)))
    product_stationary_residual=float(np.linalg.norm(V@Rprod-Rprod@V))
    accel=-marginal(V@(V@Rprod-Rprod@V)-(V@Rprod-Rprod@V)@V,n)
    assert np.linalg.norm(accel-product_acceleration(C))<1e-13
    times=[.03,.2,math.pi/5,math.pi,2*math.pi]
    rows=[]
    for t in times:
        rho,R=product_evolution(C,t)
        residual=float(np.linalg.norm(rho-product_closed(C,t)))
        assert residual<2e-13
        rows.append(dict(pair_phase=t,closed_form_error=residual,
            marginal_distance=float(np.linalg.norm(rho-C)),entropy_gain=entropy(rho)-entropy(C),
            uniform_population_error=float(np.max(abs(np.diag(rho)-1/n))),
            exact_antisymmetric_source_stays_stationary=True))
    RA=antisymmetric_completion(C)
    anti_eig=np.linalg.eigvalsh(RA)
    pt_eig=np.linalg.eigvalsh(partial_transpose(RA,n))
    assert anti_eig[0]>-1e-13 and abs(np.trace(RA)-1)<1e-13
    assert np.linalg.norm(marginal(RA,n)-C)<1e-13
    assert np.linalg.norm(V@RA-RA@V)<1e-13
    assert pt_eig[0]<-1/n+1e-12 and np.sum(pt_eig<-1e-10)==1
    H=logm(C);L=np.kron(H,np.eye(n))+np.kron(np.eye(n),H)
    legacy=legacy_cycle();C_legacy=np.eye(n)/n+.001*legacy
    legacy_product=np.kron(C_legacy,C_legacy)
    legacy_stationary=antisymmetric_completion(C_legacy)
    legacy_out,_=product_evolution(C_legacy,.2)
    assert np.linalg.eigvalsh(C_legacy)[0]>1/22
    assert np.linalg.norm(legacy_out-product_closed(C_legacy,.2))<1e-12
    return dict(status='Research checkpoint; conditional quantum lift, not sourced fundamental physics.',
        parameters=dict(n=n,gamma=gamma),
        observable_lift_residual=exact_lift_error,hartree_stationary_residual=hartree_residual,
        exact_product_stationary_residual=product_stationary_residual,
        marginal_acceleration_norm=float(np.linalg.norm(accel)),
        product_log_collective_commutator=float(np.linalg.norm(V@L-L@V)),
        local_collective_commutant_numerical_ranks={str(k):local_collective_commutant_rank(k) for k in [2,3,4,5]},
        predicted_ranks={str(k):k*k-1 for k in [3,4,5]},
        product_trajectory=rows,
        antisymmetric_completion=dict(trace=float(np.trace(RA)),
            minimum_eigenvalue=float(anti_eig[0]),marginal_residual=float(np.linalg.norm(marginal(RA,n)-C)),
            stationary_residual=float(np.linalg.norm(V@RA-RA@V)),
            partial_transpose_negative_eigenvalues=pt_eig[pt_eig<-1e-10].tolist(),
            negativity=float(-sum(pt_eig[pt_eig<0])),
            minimum_input_density_eigenvalue=float(np.linalg.eigvalsh(C)[0]),
            sufficient_completion_density_floor=1/(2*(n-1))),
        separate_legacy_cycle_check=dict(loading=.001,
            minimum_program_eigenvalue=float(np.linalg.eigvalsh(C_legacy)[0]),
            product_commutator_norm=float(np.linalg.norm(V@legacy_product-legacy_product@V)),
            correlated_commutator_norm=float(np.linalg.norm(V@legacy_stationary-legacy_stationary@V)),
            common_source_amplitude='4 ln(2); phase angles are rational multiples of pi',
            arithmetic_scope='Canonical legacy ratios are algebraic; the strict P512 transcendence obstruction is not transferred.'),
        scope=['The tensor-pair microscopic law is a declared completion, not a strict-source theorem.',
            'Stationary correlated completion is target dependent, not a universal quantum broadcasting channel.',
            'Mean-field stationary does not imply finite-N product stationary.'])


if __name__=='__main__':
    result=run();Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
