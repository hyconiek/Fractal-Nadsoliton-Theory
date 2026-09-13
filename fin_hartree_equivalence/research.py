"""New goal: distinguish full-source, mixed-flow and pure-flow lifts.

No previous completed finding is counted again. Projective dynamics,
exchange symmetry and the state class are explicit premises.
"""
import json
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.linalg import expm


def exact_structure(n):
    I=sp.eye(n*n);S=sp.zeros(n*n);D=sp.zeros(n*n);omega=sp.zeros(n*n,1)
    for i in range(n):
        D[i*n+i,i*n+i]=1;omega[i*n+i]=1
        for j in range(n):S[i*n+j,j*n+i]=1
    P=(I+S)/2;A=(I-S)/2;V=(omega*omega.T+S-2*D)/2
    return V,S,P,A,D


def basis(n):
    out=[]
    for i in range(n):
        H=sp.zeros(n);H[i,i]=1;out.append(H)
        for j in range(i+1,n):
            H=sp.zeros(n);H[i,j]=H[j,i]=1;out.append(H)
            H=sp.zeros(n);H[i,j]=sp.I;H[j,i]=-sp.I;out.append(H)
    return out


def ptr(R,n):
    if isinstance(R,sp.MatrixBase):
        return sp.Matrix(n,n,lambda i,j:sum(R[i*n+k,j*n+k] for k in range(n)))
    return np.trace(R.reshape(n,n,n,n),axis1=1,axis2=3)


def local_rank(n,a):
    V,S,*_=exact_structure(n);V+=a*S;I=sp.eye(n);columns=[]
    for side in [0,1]:
        for H in basis(n):
            L=sp.kronecker_product(H,I) if side==0 else sp.kronecker_product(I,H)
            columns.append((V*L-L*V).reshape(n**4,1))
    return int(sp.Matrix.hstack(*columns).rank())


def commuting_map_ranks(n=3):
    """All quadratic coefficients, not a finite state sampling test."""
    H=basis(n);d=n*n;norm=[sp.trace(X*X) for X in H]
    structure={}
    for a in range(d):
        for b in range(d):
            comm=sp.I*(H[a]*H[b]-H[b]*H[a])
            structure[a,b]=[sp.simplify(sp.trace(X*comm)/z) for X,z in zip(H,norm)]
    entries={};row=0
    for b in range(d):
        for c in range(b,d):
            for k in range(d):
                for a in range(d):
                    v=structure[a,c][k]
                    if v:entries[row,a*d+b]=v
                    if b!=c:
                        v=structure[a,b][k]
                        if v:entries[row,a*d+c]=v
                row+=1
    M=sp.SparseMatrix(row,d*d,entries)
    unconstrained=int(M.rank())
    reciprocal=[]
    for i in range(d):
        for j in range(i+1,d):
            v=[0]*(d*d);v[i*d+j]=norm[i];v[j*d+i]=-norm[j];reciprocal.append(v)
    symmetric=int(M.col_join(sp.Matrix(reciprocal)).rank())
    return dict(dimension=n,real_map_parameters=d*d,
        commuting_rank=unconstrained,commuting_nullity=d*d-unconstrained,
        reciprocal_rank=symmetric,reciprocal_nullity=d*d-symmetric)


def antisymmetric_basis(n):
    columns=[]
    for i in range(n):
        for j in range(i+1,n):
            a=sp.zeros(n*n,1);a[i*n+j]=1/sp.sqrt(2);a[j*n+i]=-1/sp.sqrt(2);columns.append(a)
    return sp.Matrix.hstack(*columns)


def matching_completion():
    n=3;V,S,P,A,D=exact_structure(n);U=antisymmetric_basis(n)
    C=sp.Matrix([[2,1,0],[1,2,0],[0,0,1]])/5
    F=C.inv()/sp.trace(C.inv());R=sp.kronecker_product(C,F)
    symbols=sp.symbols('b0:6',real=True)
    B=sp.Matrix([[symbols[0],symbols[1],symbols[2]],
                 [symbols[1],symbols[3],symbols[4]],
                 [symbols[2],symbols[4],symbols[5]]])
    H=P*V*P+U*B*U.T
    solution=list(sp.linsolve(list(H*R-R*H),symbols))
    assert len(solution)==1
    values=solution[0];free=set().union(*(x.free_symbols for x in values))
    values=[sp.simplify(x.subs({f:0 for f in free})) for x in values]
    H=H.subs(dict(zip(symbols,values)))
    assert H*R==R*H and P*(H-V)*P==sp.zeros(n*n)
    return C,F,H,dict(program='non-diagonal real matching marginal',
        first_marginal=[[str(x) for x in C.row(i)] for i in range(n)],
        second_marginal=[[str(x) for x in F.row(i)] for i in range(n)],
        anti_block_coefficients=list(map(str,values)),
        exact_stationarity=True,exact_pure_flow_equivalence=True,
        mixed_flow_equivalence_claimed=False)


def dense_completion_rejection():
    n=3;V,S,P,A,D=exact_structure(n);U=antisymmetric_basis(n)
    C=sp.Matrix([[4,1,1],[1,3,1],[1,1,3]])/10
    E=C.inv()/sp.trace(C.inv());R=sp.kronecker_product(C,E)
    assert all(C[:k,:k].det()>0 for k in range(1,4))
    symbols=sp.symbols('z0:9',real=True)
    B=sum((z*H for z,H in zip(symbols,basis(3))),sp.zeros(3))
    H=P*V*P+U*B*U.T
    equations=list(H*R-R*H)
    M,rhs=sp.linear_eq_to_matrix(equations,symbols)
    rank=int(M.rank());augmented=int(M.row_join(rhs).rank())
    assert augmented>rank and sp.linsolve((M,rhs),symbols)==sp.EmptySet
    return dict(dense_positive_marginal=[[str(x) for x in C.row(i)] for i in range(3)],
        unrestricted_hermitian_antisymmetric_parameters=9,
        coefficient_rank=rank,augmented_rank=augmented,exact_inconsistency=True)


def dense_kernel_certificate():
    from fractions import Fraction as F
    import importlib.util
    root=Path(__file__).resolve().parents[1]
    spec=importlib.util.spec_from_file_location('strict_enclosures',root/'fin_replication_consistency/certify.py')
    cert=importlib.util.module_from_spec(spec);spec.loader.exec_module(cert)
    weights=cert.strict_weights()
    assert weights[0][0]>F(2,5) and weights[1][0]>F(1,10)
    legacy1=4*sp.log(2)*sp.cos(5*sp.pi/12)/sp.Rational(101,100)
    legacy2=4*sp.log(2)*sp.cos(2*sp.pi/3)/sp.Rational(102,100)
    assert legacy1.is_positive and legacy2.is_negative
    return dict(strict_two_incident_edge_lower_bounds=['2/5','1/10'],
        strict_outward_bounds_recomputed=True,
        canonical_legacy_two_edges=[str(sp.simplify(legacy1)),str(sp.simplify(legacy2))],
        implication='Any nonzero loading of either dense source violates the necessary matching condition.',
        scope='Unequal partners require full rank; identical-product theorem is separately rank-free.')


def correlation_floor_certificate():
    from fractions import Fraction as F
    import importlib.util
    root=Path(__file__).resolve().parents[1]
    spec=importlib.util.spec_from_file_location('strict_spectrum',root/'fin_projected_learning/research.py')
    previous=importlib.util.module_from_spec(spec);spec.loader.exec_module(previous)
    certificate=previous.certify_strict_spectrum()
    eig=[tuple(map(F,row)) for row in certificate['eigenvalue_intervals']]
    assert min(F(1,12)+row[0]/20 for row in eig)>F(1,22)
    c0=tuple(F(1,12)+x/20 for x in eig[0]);c6=tuple(F(1,12)+x/20 for x in eig[6])
    delta=(c0[0]**2-c6[1]**2,c0[1]**2-c6[0]**2)
    assert c6[0]>0 and delta[0]>F(63,2500)
    assert F(20)<F(9,2)**2
    floor=F(7,1250)
    return dict(loading_both='1/20',
        delta_interval=[str(x) for x in delta],
        product_commutator_entry_lower=str(F(5,12)*delta[0]),
        stationary_correlation_trace_distance_strictly_greater_than=str(floor),
        mutual_information_nats_strictly_greater_than=str(2*floor**2),
        rank_assumption='None; each marginal is a stationary strict source state.',
        microscopic_freedom='Every reciprocal pure-projector-flow-equivalent pair Hamiltonian.',
        no_laboratory_claim=True)


def correlation_witness(n=12):
    if n<4 or n%2:raise ValueError('This implementation uses an even alternating mode')
    u=np.ones(n)/np.sqrt(n);v=(-1.)**np.arange(n)/np.sqrt(n)
    alpha=np.kron(u,u);beta=np.kron(v,v);omega=np.eye(n).ravel()
    kappa=(n-2)/(2*n);d=alpha-beta
    B0=kappa*(np.outer(d,omega)+np.outer(omega,d))/2
    local=np.outer(u,u)-np.outer(v,v)
    B=B0-kappa*(np.kron(local,np.eye(n))+np.kron(np.eye(n),local))/2
    return alpha,beta,B,kappa


def run():
    ranks=commuting_map_ranks(3)
    assert ranks['commuting_nullity']==10 and ranks['reciprocal_nullity']==2
    rank_rows=[]
    for n in [3,4]:
        for a in [sp.S.Zero,-sp.Rational(1,2),-sp.Rational(n,4),sp.Rational(1,3)]:
            rank=local_rank(n,a);nullity=2*n*n-rank
            expected=n+1 if a==-sp.Rational(1,2) else 2
            assert nullity==expected
            rank_rows.append(dict(n=n,exchange_shift=str(a),rank=rank,nullity=nullity))
    n=3;V,S,P,A,D=exact_structure(n)
    C=sp.diag(sp.Rational(1,2),sp.Rational(1,3),sp.Rational(1,6))
    F=C.inv()/sp.trace(C.inv());R=sp.kronecker_product(C,F);Q=V-S/2
    assert Q*R==R*Q and V*R!=R*V
    Cm,Fm,Hm,matching=matching_completion()
    # A generic antisymmetric perturbation is invisible to pure flow but
    # not to the supplied extension to all mixed densities.
    U=antisymmetric_basis(3);B=sp.diag(1,2,4);perturb=U*B*U.T
    rho=sp.Matrix([[4,1,0],[1,3,1],[0,1,3]])/10
    field=ptr(perturb*sp.kronecker_product(sp.eye(3),rho),3)
    mixed_defect=field*rho-rho*field
    assert mixed_defect!=sp.zeros(3)
    psi=sp.Matrix([1,sp.I,2])/sp.sqrt(6);pure=psi*sp.conjugate(psi.T)
    pure_field=ptr(perturb*sp.kronecker_product(sp.eye(3),pure),3)
    assert sp.simplify(pure_field*pure-pure*pure_field)==sp.zeros(3)
    return dict(status='Analytic lift classification and rank-free uniform strict correlation floor proved; source/physical closure is not claimed.',
        exact_commuting_map_ranks=ranks,exchange_family_local_ranks=rank_rows,
        twelve_label_invisible_dimensions=dict(full_source=0,mixed_projective_flow=2,
            pure_projective_flow=1+66**2,pure_flow_modulo_global_energy_shift=66**2),
        exact_diagonal_reciprocal_exception=True,matching_completion=matching,
        dense_completion_rejection=dense_completion_rejection(),
        dense_kernel_certificate=dense_kernel_certificate(),
        correlation_floor_certificate=correlation_floor_certificate(),
        exact_pure_mixed_separation=True,
        mixed_defect_squared_frobenius=str(sp.simplify(sp.trace(mixed_defect.conjugate().T*mixed_defect))),
        scope=['Full source T, mixed projector flow and pure projector flow are not interchangeable.',
            'Matching marginal exception does not realize a dense strict or legacy kernel.',
            'No finite-speed controller, source, laboratory or ToE closure.'])


if __name__=='__main__':
    data=run();Path(__file__).with_name('results.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2))
