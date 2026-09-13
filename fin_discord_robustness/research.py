"""New goal: robustness beyond one canonical microscopic realization.

The task concerns two 12-level factors and stated equilibrium/source
premises. Exchange asymmetry is not the internal Z12 orientation selector.
"""
from fractions import Fraction as F
import importlib.util
import json
from pathlib import Path

import numpy as np
import sympy as sp


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('previous_quantum_lift',ROOT/'fin_quantum_lift/research.py')
q=importlib.util.module_from_spec(spec);spec.loader.exec_module(q)


def polynomial_certificate(n=12):
    a,x=sp.symbols('a x',real=True)
    eig=[a+sp.Rational(n-1,2),a+sp.Rational(1,2),a-sp.Rational(1,2),-a-sp.Rational(1,2)]
    def projector(index):
        return sp.prod((x-eig[j])/(eig[index]-eig[j]) for j in range(4) if j!=index)
    poly=sp.factor(-sp.Rational(1,2)+projector(1)+sp.Rational(n,2)*projector(0))
    target=[sp.Rational(n-1,2),sp.Rational(1,2),-sp.Rational(1,2),-sp.Rational(1,2)]
    assert all(sp.simplify(poly.subs(x,e)-t)==0 for e,t in zip(eig,target))
    denominator=sp.factor(sp.denom(poly))
    roots=sp.solve(denominator,a)
    assert set(roots)=={-sp.Rational(1,2),-sp.Rational(n,4)}
    old_exception=-sp.Rational(n-2,2*n)
    assert denominator.subs(a,old_exception)!=0
    return dict(n=n,reconstruction_polynomial=str(poly),denominator=str(denominator),
        genuine_transfer_exceptions=[str(z) for z in roots],
        old_singular_block_value=str(old_exception),old_value_now_paid=True,
        exact_four_band_reconstruction=True)


def exception_bounds(n=12):
    k=F(n-2,2*n);rows=[]
    for a in [F(-1,2),F(-n,4)]:
        energies=[a+F(n-1,2),a+F(1,2),a-F(1,2),-a-F(1,2)]
        inverse=max(1/k,1/abs(a+k),F(n))
        norm=max(abs(v) for v in energies)
        assert inverse<=n and norm<=F(n-1,2)
        rows.append(dict(exchange_shift=str(a),inverse_cross_block_norm=str(inverse),
            Hamiltonian_norm=str(norm),canonical_bound_remains_valid=True))
    return rows


def tradeoff_certificate():
    spec=importlib.util.spec_from_file_location('strict_spectral_provider',ROOT/'fin_projected_learning/research.py')
    previous=importlib.util.module_from_spec(spec);spec.loader.exec_module(previous)
    source=previous.certify_strict_spectrum();eigen=[tuple(map(F,row)) for row in source['eigenvalue_intervals']]
    gap=F(source['minimum_intersector_gap_lower'])/20
    difference=(eigen[0][0]-eigen[6][1])/20
    assert eigen[0][1]<F(5,3)
    program_minimum=min(F(1,12)+F(11,100)*lo for lo,hi in eigen)
    assert program_minimum>F(1,125)
    floor=gap*difference/7392;residual_coefficient=gap/616
    assert floor>0
    return dict(strict_loading='1/20',density_gap_lower=str(gap),mode_weight_difference_lower=str(difference),
        uniform_CQ_distance_lower=str(floor),uniform_CQ_distance_lower_decimal=float(floor),
        stationarity_residual_coefficient=str(residual_coefficient),
        amplified_program_eigenvalue_lower=str(program_minimum),
        required_exchange_asymmetry_if_stationary_and_CQ=str(2*floor),
        inequality='d_CQ_A(R) + D_tr(R,Swap R Swap)/2 + coefficient*||[H,R]||_1 >= floor',
        premise='H=V+c P_sym+B_anti; average of the two labelled marginals is the stated strict C.',
        attainability_claimed=False)


def exceptional_state(n,a):
    V,S,D,P=q.structures(n);d=n*n
    anti=np.zeros(d);anti[1]=1/np.sqrt(2);anti[n]=-1/np.sqrt(2)
    if a==-.5:
        symmetric=np.zeros(d);symmetric[2]=symmetric[2*n]=1/np.sqrt(2)
    else:symmetric=np.eye(n).ravel()/np.sqrt(n)
    R=np.eye(d)/d+(np.outer(symmetric,anti)+np.outer(anti,symmetric))/(2*d)
    H=V+a*S
    assert np.linalg.norm(H@R-R@H)<1e-13
    assert np.linalg.norm(V@R-R@V)>1e-4
    return H,R


def asymmetric_cq_construction(gamma=.05):
    n=12;W=q.strict();C=np.eye(n)/n+gamma*W
    program=np.eye(n)/n+gamma*F(11,5)*W
    program=np.asarray(program,dtype=float)
    assert np.linalg.eigvalsh(program)[0]>0
    R=np.zeros((n*n,n*n));branches=[]
    for i in range(n):
        first=np.zeros((n,n));first[i,i]=1
        mask=np.ones(n);mask[i]=0
        second=n/(n-1)*mask[:,None]*program*mask[None,:]
        branches.append((first,second));R+=np.kron(first,second)/n
    V,S,D,P=q.structures(n);H=V-S/2;Rbar=(R+S@R@S)/2
    assert np.linalg.norm(H@R)<1e-13 and np.linalg.norm(R@H)<1e-13
    assert np.linalg.norm(V@Rbar-Rbar@V)<1e-13
    left=q.marginal(R,n);right=q.marginal(S@R@S,n)
    assert np.linalg.norm(left-np.eye(n)/n)<1e-13
    assert np.linalg.norm(right-(np.eye(n)/n+2*gamma*W))<1e-13
    assert np.linalg.norm((left+right)/2-C)<1e-13
    return R,Rbar,program,branches


def stationary_broadcast_channels(rho):
    n=len(rho);V,S,D,P=q.structures(n)
    plus=(np.eye(n*n)+S)/2-D;minus=(np.eye(n*n)-S)/2
    left=np.kron(rho,np.eye(n));scale=2/(n-1)
    Eplus=scale*plus@left@plus;Eminus=scale*minus@left@minus
    balanced=(Eplus+Eminus)/2
    cq=np.zeros((n*n,n*n),complex)
    for i in range(n):
        basis=np.zeros((n,n));basis[i,i]=1;Q=np.eye(n)-basis
        cq+=np.kron(basis,Q@rho@Q)/(n-1)
    assert np.linalg.norm(balanced-(cq+S@cq@S)/2)<1e-13
    return Eplus,Eminus,balanced,cq


def broadcasting_certificate(n=12):
    p=F(n-2,2*(n-1))
    choi_pt=(1-(n+1)*p)/(n*n)
    return dict(marginal_depolarizing_factor=str(p),marginal_identity_coefficient=str(F(1,2*(n-1))),
        deterministic_CPTP_channels=True,
        family='r E_plus+(1-r)E_minus, 0<=r<=1',
        exact_separability_condition='r=1/2 for every input density',
        partial_transpose_minor_offdiagonal='(2r-1)(rho_ii+rho_jj)/(2(n-1))',
        uniform_diagonal_negativity_lower='abs(2r-1)/n',
        local_Choi_partial_transpose_eigenvalue=str(choi_pt),
        one_fixed_party_LOCC_without_quantum_transfer_ruled_out=bool(choi_pt<0),
        LOCC_scope='Universal channel on an unknown input, potentially entangled with a reference; not preparation of one fixed known separable target.',
        physical_source_claimed=False)


def run():
    polynomial=polynomial_certificate();bounds=exception_bounds();tradeoff=tradeoff_certificate()
    n=12;V,S,D,P=q.structures(n);Ps=(np.eye(n*n)+S)/2;Pa=np.eye(n*n)-Ps
    W=q.strict();C=np.eye(n)/n+.05*W
    u=np.ones(n)/np.sqrt(n);v=(-1.)**np.arange(n)/np.sqrt(n)
    A=np.outer(u,u)-np.outer(v,v);B=np.outer(u,v)+np.outer(v,u)
    # CQ on the first factor, same marginals, but its swap average need not be CQ.
    cq=np.kron(C,C)+.0001*np.kron(A,B)
    average=(cq+S@cq@S)/2
    assert np.linalg.eigvalsh(cq)[0]>0
    assert np.linalg.norm(cq@np.kron(C,np.eye(n))-np.kron(C,np.eye(n))@cq)<1e-13
    created=np.linalg.norm(average@np.kron(C,np.eye(n))-np.kron(C,np.eye(n))@average)
    assert created>1e-6
    rng=np.random.default_rng(8658);raw=rng.normal(size=(n*n,n*n))+1j*rng.normal(size=(n*n,n*n));raw=(raw+raw.conj().T)/2
    H=V+17*Ps+Pa@raw@Pa
    arbitrary=rng.normal(size=(n*n,4))+1j*rng.normal(size=(n*n,4));R=arbitrary@arbitrary.conj().T;R/=np.trace(R)
    Rbar=(R+S@R@S)/2
    identity_error=np.linalg.norm(V@Rbar-Rbar@V-Ps@(H@R-R@H)@Ps)
    assert identity_error<1e-13
    counterexamples=[]
    for a in [-.5,-3.]:
        H0,R0=exceptional_state(n,a)
        counterexamples.append(dict(a=a,stationary_residual=float(np.linalg.norm(H0@R0-R0@H0)),
            canonical_residual=float(np.linalg.norm(V@R0-R0@V)),
            minimum_eigenvalue=float(np.linalg.eigvalsh(R0)[0]),
            strict_marginals_claimed=False))
    asymmetric,symmetric,program,branches=asymmetric_cq_construction()
    asymmetry=float(sum(abs(np.linalg.eigvalsh(asymmetric-S@asymmetric@S)))/2)
    assert asymmetry>2*float(F(tradeoff['uniform_CQ_distance_lower']))
    assert np.linalg.norm(asymmetric-S@asymmetric@S)>1e-3
    # Passive channels agree after twirling; coherent controlled-U is a
    # stronger resource and distinguishes the models at the same time.
    from scipy.linalg import expm
    t=np.pi
    control_canonical=float(np.trace(expm(-1j*t*V)@symmetric).real)
    control_shifted=float(np.trace(expm(-1j*t*(V-S/2))@asymmetric).real)
    assert abs(control_canonical)<1e-12 and abs(control_shifted-1)<1e-12
    positive,negative,balanced,cq_state=stationary_broadcast_channels(program)
    np.testing.assert_allclose(balanced,symmetric,atol=1e-13)
    np.testing.assert_allclose(cq_state,asymmetric,atol=1e-13)
    expected=(np.eye(n)+(n-2)*program)/(2*(n-1))
    channel_rows=[]
    for r in [0.,.25,.5,.75,1.]:
        output=r*positive+(1-r)*negative
        pt=q.partial_transpose(output,n);ev=np.linalg.eigvalsh(pt)
        negativity=float(-sum(ev[ev<0]))
        assert np.linalg.norm(q.marginal(output,n)-expected)<1e-13
        assert np.linalg.norm(V@output-output@V)<1e-13
        assert negativity+1e-12>=abs(2*r-1)/n
        if r==.5:assert negativity<1e-12
        channel_rows.append(dict(mixing_weight=r,negativity=negativity,
            swap_expectation=float(np.trace(S@output).real)))
    return dict(status='Mixed-flow discord robustness, pure-flow tradeoff and stationary-channel identifiability proved; no fundamental source closure.',
        polynomial_transfer=polynomial,exceptional_block_bounds=bounds,
        exact_tradeoff=tradeoff,
        numerical_checks=dict(pure_family_symmetrization_identity_error=float(identity_error),
            classical_input_swap_average_marginal_commutator=float(created),
            CQ_input_stationarity_claimed=False),
        actual_transfer_counterexamples=counterexamples,
        asymmetric_CQ_witness=dict(first_marginal='I/12',second_marginal='I/12+W/10',
            average_marginal='I/12+W/20',same_individual_marginals_claimed=False,
            program_loading='11/100',preparation_success_probability='11/12',
            alternative_deterministic_Luders_instrument=True,
            stationary_exchange_shift='-1/2',one_sided_CQ_on_first=True,
            trace_distance_exchange_asymmetry=asymmetry,
            symmetric_state_nonclassicality_commutator=float(np.linalg.norm(symmetric@np.kron(C,np.eye(n))-np.kron(C,np.eye(n))@symmetric)),
            controlled_U_X_expectations_at_time_pi=[control_canonical,control_shifted],
            passive_equivalence_scope='Same fixed outcome-wise swap-covariant instruments and uncontrolled free-evolution channels; not controlled-Hamiltonian or energy-query access.'),
        exact_broadcasting_certificate=broadcasting_certificate(),
        stationary_channel_mixture_checks=channel_rows,
        scope=['The mixed-flow family has a uniform non-CQ equilibrium bound for all a,b.',
            'The pure-flow family requires a state-symmetry premise or the explicit asymmetry/residual tradeoff.',
            'Exchange of the two factors is not the Z12 orientation problem or a physical selector source.'])


if __name__=='__main__':
    result=run();Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
