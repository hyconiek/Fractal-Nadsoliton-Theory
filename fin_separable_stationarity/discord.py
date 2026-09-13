"""Exact block obstruction and conservative discord-distance certificates."""
from fractions import Fraction as F
import importlib.util
import json
from pathlib import Path

import numpy as np
import sympy as sp

try:
    from . import research as r
except ImportError:
    import research as r


def exact_block_certificate(n=12):
    if n<4 or n%2:raise ValueError('Even alternating-mode certificate required')
    u=sp.ones(n,1)/sp.sqrt(n);v=sp.Matrix([(-1)**i for i in range(n)])/sp.sqrt(n)
    P=u*u.T+v*v.T;J=sp.diag(*[(-1)**i for i in range(n)])
    A=(u*v.T+v*u.T)/2-J/n;k=sp.Rational(n-2,2*n)
    assert sp.simplify(A*A-((sp.eye(n)-P)/n**2+k*k*P))==sp.zeros(n)
    inverse=-n*J+(n+1/k)*(u*v.T+v*u.T)
    assert sp.simplify(A*inverse)==sp.eye(n)
    shifted=A-k*v*u.T
    assert shifted.rank()==n-1
    return dict(n=n,inverse_norm=str(n),smallest_singular_value=str(sp.Rational(1,n)),
        exceptional_exchange_shift=str(-k),exceptional_block_rank=n-1,
        exact_inverse_and_square_identity=True)


def distance_certificate():
    spec=importlib.util.spec_from_file_location('strict_spectrum',r.ROOT/'fin_projected_learning/research.py')
    previous=importlib.util.module_from_spec(spec);spec.loader.exec_module(previous)
    source=previous.certify_strict_spectrum();eig=[tuple(map(F,p)) for p in source['eigenvalue_intervals']]
    spec=importlib.util.spec_from_file_location('strict_weights',r.ROOT/'fin_replication_consistency/certify.py')
    weights=importlib.util.module_from_spec(spec);spec.loader.exec_module(weights)
    weight_boxes=weights.strict_weights()
    assert weight_boxes[0][0]>F(2,5) and all(lo>0 for lo,hi in weight_boxes)
    assert eig[0][1]<F(5,3)
    gap=F(source['minimum_intersector_gap_lower'])/20
    difference=(eig[0][0]-eig[6][1])/20
    universal=gap*difference/F(7392)
    assert universal>0
    return dict(explicit_separable_state_CQ_trace_distance_strictly_greater_than='1/15400',
        strict_first_edge_lower='2/5',all_strict_edge_lower_bounds_positive=True,
        strict_density_gap_lower=str(gap),uniform_alternating_weight_difference_lower=str(difference),
        every_canonical_stationary_state_CQ_trace_distance_lower=str(universal),
        universal_lower_decimal=float(universal),
        scope='Distance to the one-sided classical-quantum set; not entropic discord, negativity or laboratory evidence.')


def legacy_certificate():
    field=sp.QQ.algebraic_field(sp.sqrt(2),sp.sqrt(3))
    weights=[sp.cos(sp.pi*d/4+sp.pi/6)/(1+sp.Rational(d,100)) for d in range(1,7)]
    eigen=[sp.simplify(2*sum(weights[d-1]*sp.cos(sp.pi*k*d/6) for d in range(1,6))+(-1)**k*weights[5]) for k in range(7)]
    for mode in [0,6]:
        for other in range(7):
            if other!=mode:assert field.from_sympy(eigen[mode]-eigen[other])!=field.zero
    program_floor=F(1,12)-F(11,5000)*33
    assert program_floor>0
    return dict(output_loading='1/1000',amplified_loading='11/5000',
        coarse_program_eigenvalue_lower=str(program_floor),
        common_factor_removed='4 ln(2)',simple_modes_0_and_6_exactly_verified=True,
        field='Q(sqrt(2),sqrt(3))',
        scope='Separate legacy separable construction and qualitative no-zero-discord result; strict numeric bounds are not transferred.')


def run():
    n=12;W=r.q.strict();C=np.eye(n)/n+.05*W
    R,_,_,_=r.stationary_separable(C)
    comm=R@np.kron(C,np.eye(n))-np.kron(C,np.eye(n))@R
    comm_norm=float(sum(np.linalg.svd(comm,compute_uv=False)))
    c0=float(np.linalg.eigvalsh(C)[-1])
    numerical_distance_bound=comm_norm/(4*(1+c0))
    assert abs(comm[0,n])>1/6600
    assert numerical_distance_bound>1/15400
    return dict(status='Canonical stationary nonclassicality proved; separability separately constructed.',
        exact_block=exact_block_certificate(),distance_certificate=distance_certificate(),
        separate_legacy_certificate=legacy_certificate(),
        numerical_witness=dict(commutator_trace_norm=comm_norm,
            commutator_entry=complex(comm[0,n]).real,
            CQ_distance_lower_from_computed_norm=numerical_distance_bound),
        scope=['The invertible-block theorem uses simple unequal marginal eigenvalues.',
            'The full Hartree-invisible interaction class is not covered automatically.',
            'Classical shared randomness in preparation does not imply zero quantum discord of its output.'])


if __name__=='__main__':
    data=run();Path(__file__).with_name('discord_results.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2))
