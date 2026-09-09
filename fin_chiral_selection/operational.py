"""ST8652: exact ensemble obstruction for the original PURE feedback law.

No-signalling interpretation additionally assumes standard remote preparation,
Born readout and branchwise evolution of conditionally prepared pure states.
The algebraic failure of an affine density-map extension needs no relativity.
"""
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

try:
    from . import research as r
except ImportError:
    import research as r


def pure_rhs(t,state,n,eta,gamma,propagation='K'):
    psi=state[:n];K=state[n:].reshape(n,n)
    if propagation not in ['K','laplacian']:raise ValueError(propagation)
    H=K if propagation=='K' else K-np.diag(K.real.sum(axis=1))
    return np.r_[1j*H@psi,
        (eta*(r.off(np.outer(psi,psi.conj()).real)-gamma*K.real)).ravel()]


def pure_trajectory(psi,K,times,eta=.2,gamma=.05,propagation='K'):
    n=len(psi)
    solution=solve_ivp(lambda t,y:pure_rhs(t,y,n,eta,gamma,propagation),
        (0,max(times)),np.r_[psi,K.ravel()].astype(complex),t_eval=times,
        method='DOP853',rtol=3e-12,atol=3e-14,max_step=.002)
    if not solution.success:raise RuntimeError(solution.message)
    return solution.y[:n].T,solution.y[n:].T.reshape(-1,n,n)


def exact_defect():
    import sympy as s
    I=s.eye(2);X=s.Matrix([[0,1],[1,0]]);Y=s.Matrix([[0,-s.I],[s.I,0]])
    P=s.Matrix([[9,3],[3,1]])/10;Q=I-P
    def off(A):return A-s.diag(*A.diagonal())
    comm=lambda A,B:A*B-B*A
    defect=s.I*(comm(off(P),P)+comm(off(Q),Q))
    assert P*P==P and Q*Q==Q and P*Q==s.zeros(2) and P+Q==I
    assert defect==s.Rational(12,25)*Y
    # Canonical vertex readout, with arbitrary couplings to two extra sites.
    kk=s.symbols('k01 k02 k03 k12 k13 k23',real=True)
    K=s.Matrix([[0,kk[0],kk[1],kk[2]],[kk[0],0,kk[3],kk[4]],
                [kk[1],kk[3],0,kk[5]],[kk[2],kk[4],kk[5],0]])
    D=s.zeros(4);laplacian_defect=s.zeros(4)
    for small in [P,Q]:
        big=s.zeros(4);big[:2,:2]=small;A=off(big)
        D+=-2*comm(A,comm(K,big))-comm(K,comm(A,big))
        effective=A-s.diag(*list(A*s.ones(4,1)))
        laplacian_defect+=s.I*comm(effective,big)
    diagonal=[s.simplify(D[i,i]) for i in range(4)]
    assert diagonal==[-s.Rational(72,25)*kk[0],s.Rational(72,25)*kk[0],0,0]
    expected_lap=s.zeros(4);expected_lap[:2,:2]=defect
    assert laplacian_defect==expected_lap
    Z=s.diag(1,-1);Z0=s.diag(1,0);Z1=s.diag(0,1)
    second_moment=(s.kronecker_product(P,P)+s.kronecker_product(Q,Q)
                   -s.kronecker_product(Z0,Z0)-s.kronecker_product(Z1,Z1))/2
    contrast=s.trace(s.kronecker_product(X,Z)*second_moment)
    assert contrast==s.Rational(12,25)
    # For n=12 and eta=1/5: Delta rho''(0)=Y/125,
    # Delta rho(t)=t²Y/250+O(t³).
    return dict(projectors=[[[str(v) for v in P.row(i)] for i in range(2)],
                            [[str(v) for v in Q.row(i)] for i in range(2)]],
        exact_pure_projectors=True,identical_ensemble_density=True,
        pair_acceleration_defect_over_eta='12/25 * sigma_y',
        general_density_leading_coefficient='6*eta/(25*n) * t^2 * sigma_y',
        strict12_eta_one_fifth_density_leading_coefficient='1/250 * t^2 * sigma_y',
        canonical_vertex_probability_rotated_minus_basis='-12*eta*K01/(25*n) * t^3 + O(t^4)',
        generic_extra_coupling_cancellation_exact=True,
        replacing_K_propagation_by_graph_Laplacian_keeps_second_order_defect=True,
        two_copy_pair_ensemble_X_tensor_Z_contrast='12/25',
        applies_to_every_common_initial_real_symmetric_K=True,
        independent_of_gamma=True)


def finite_strict_certificate():
    """An explicit finite positive probability gap, not only Taylor existence.

    The norm bounds are conservative exact rational inequalities. The selected
    strict weights are independently re-enclosed on every invocation.
    """
    from fractions import Fraction as F
    import importlib.util
    root=Path(__file__).resolve().parents[1]
    spec=importlib.util.spec_from_file_location('strict_weights',root/'fin_replication_consistency/certify.py')
    cert=importlib.util.module_from_spec(spec);spec.loader.exec_module(cert)
    weights=cert.strict_weights()
    initial_norm_squared_upper=12*(2*sum(hi*hi for lo,hi in weights[:5])+weights[5][1]**2)
    eta=F(1,5);gamma=F(1,20);T=F(1,10000);B=F(3);n=12
    assert initial_norm_squared_upper<(B-eta*T)**2
    L=eta*(1+gamma*B);J=eta*(2*B+gamma*L)
    third_density_derivative_bound=2*J+12*B*L+8*B**3
    remainder_bound=F(2,3*n)*third_density_derivative_bound*T**3
    leading=F(1,250)*T*T;lower=leading-remainder_bound
    positive_edges_lower=min(lo for lo,hi in weights)-L*T
    assert lower>0 and positive_edges_lower>0
    second_density_bound=2*L+4*B*B
    third_kernel_bound=eta*(second_density_bound+gamma*J)
    fourth_density_bound=2*third_kernel_bound+12*B*J+6*L*second_density_bound+2*B*third_density_derivative_bound
    assert weights[0][0]>F(2,5)
    vertex_leading_lower=eta*F(2,5)/25*T**3
    vertex_remainder=F(1,6*n)*fourth_density_bound*T**4
    vertex_lower=vertex_leading_lower-vertex_remainder
    assert vertex_lower>0
    return dict(time=str(T),kernel_frobenius_upper=str(B),
        kernel_derivative_frobenius_upper=str(L),
        kernel_second_derivative_frobenius_upper=str(J),
        pure_density_third_derivative_frobenius_upper=str(third_density_derivative_bound),
        readout_leading_term=str(leading),readout_absolute_remainder_bound=str(remainder_bound),
        readout_gap_lower=str(lower),readout_gap_lower_decimal=float(lower),
        any_quantum_channel_four_input_worst_trace_error_lower=str(F(n,4)*lower),
        pure_density_fourth_derivative_frobenius_upper=str(fourth_density_bound),
        vertex_probability_gap_lower=str(vertex_lower),
        vertex_probability_gap_lower_decimal=float(vertex_lower),
        vertex_remainder_bound=str(vertex_remainder),
        strict_positive_edges_proved=True,
        proof_scope='Four differing pure branches; same strict K0, eta=1/5, gamma=1/20; y+ AND canonical vertex readouts.')


def finite_two_site_certificate():
    from fractions import Fraction as F
    # Alternating Taylor inequalities, x=1/100 and t=1.
    x=F(1,100)
    exp_lower=1-x+x*x/2-x**3/6
    exp_upper=exp_lower+x**4/24
    delta_lower=12*(100*exp_lower-99)
    delta_upper=12*(100*exp_upper-99)
    assert 0<delta_lower<delta_upper<1
    trace_lower=F(2,5)*(delta_lower-delta_upper**3/6)
    edge_lower=7*exp_lower-6
    assert trace_lower>F(239,10000) and edge_lower>F(93,100)
    return dict(time='1',eta='1/5',gamma='1/20',initial_edge='1',
        delta_lower=str(delta_lower),delta_upper=str(delta_upper),
        trace_distance_lower=str(trace_lower),trace_distance_lower_decimal=float(trace_lower),
        positive_edge_lower=str(edge_lower),
        any_quantum_channel_worst_trace_error_greater_than='239/20000')


def two_site_exact(t,eta=.2,gamma=.05,k0=1.,a=.6,b=.8):
    decay=eta*gamma
    one_minus=-math.expm1(-decay*t)
    common=2*k0*one_minus/decay
    extra=a/gamma*(t-one_minus/decay)
    bloch=np.array([0,b*math.cos(common)*math.sin(extra),-b*math.sin(common)*math.sin(extra)])
    kernel_minus=k0*math.exp(-decay*t)-a/(2*gamma)*one_minus
    return dict(bloch=bloch.tolist(),trace_distance=float(np.linalg.norm(bloch)/2),
        smallest_branch_edge=kernel_minus,common_angle=common,extra_angle=extra)


def run():
    n=12;K=r.strict();times=[.0125,.025,.05,.1]
    basis=np.eye(n);states=[basis[0],basis[1],(3*basis[0]+basis[1])/math.sqrt(10),
                           (-basis[0]+3*basis[1])/math.sqrt(10)]
    paths=[pure_trajectory(psi,K,times) for psi in states]
    Y=np.zeros((n,n),complex);Y[0,1]=-1j;Y[1,0]=1j
    y_plus=(basis[0]+1j*basis[1])/math.sqrt(2)
    readout_projector=np.outer(y_plus,y_plus.conj())
    records=[]
    for index,t in enumerate(times):
        densities=[np.outer(path[0][index],path[0][index].conj()) for path in paths]
        difference=(densities[2]+densities[3]-densities[0]-densities[1])/n
        # Remaining ten basis branches cancel exactly and need not be integrated.
        readout=np.trace(readout_projector@difference).real
        records.append(dict(t=t,density_defect_frobenius=float(np.linalg.norm(difference)),
            trace_distance=float(np.sum(abs(np.linalg.eigvalsh(difference)))/2),
            y_plus_probability_difference=float(readout),
            predicted_leading_probability=t*t/250,
            probability_over_t_squared=float(readout/t**2),
            vertex_zero_probability_basis_minus_rotated=float(-difference[0,0].real),
            vertex_leading_probability=float(.2*K[0,1]/25*t**3),
            taylor_matrix_scaled_error=float(np.linalg.norm(difference/t**2-Y/250))))
    norm_error=max(abs(np.linalg.norm(psi)**2-1) for path in paths for psi in path[0])
    edge_margin=min(M.real[~np.eye(n,dtype=bool)].min() for path in paths for M in path[1])
    assert norm_error<1e-10 and edge_margin>0
    assert abs(records[0]['probability_over_t_squared']-.004)<2e-6
    # Independent exact finite-time two-site formula versus the original pure ODE.
    X=np.array([[0.,1.],[1.,0.]])
    states2=[np.array([1.,0.]),np.array([0.,1.]),np.array([3.,1.])/math.sqrt(10),
             np.array([-1.,3.])/math.sqrt(10)]
    paths2=[pure_trajectory(psi,X,[1.]) for psi in states2]
    rho=[np.outer(path[0][0],path[0][0].conj()) for path in paths2]
    diff2=(rho[2]+rho[3]-rho[0]-rho[1])/2
    exact2=two_site_exact(1.)
    bvec=exact2['bloch'];Y2=Y[:2,:2];Z=np.diag([1.,-1.])
    expected2=(bvec[1]*Y2+bvec[2]*Z)/2
    err=float(np.linalg.norm(diff2-expected2));assert err<1e-10
    return dict(status='Exact non-affinity proved; standard no-signalling interpretation is premise-explicit.',
        exact_certificate=exact_defect(),finite_strict_certificate=finite_strict_certificate(),strict12_numerical_witness=records,
        pure_state_norm_error=float(norm_error),positive_edge_sampled_margin=float(edge_margin),
        two_site_closed_form=exact2,two_site_full_ODE_agreement=err,
        finite_two_site_certificate=finite_two_site_certificate(),
        physical_scope='Not an experiment; not a universal prohibition of mean-field, stochastic, or nonstandard composite laws.')


if __name__=='__main__':
    data=run();Path(__file__).with_name('operational_results.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2))
