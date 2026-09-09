"""Second checkpoint: covariance rank and degeneracy of learning minima.

Numerical searches select witnesses only. Exact minor enclosures and the
analytic arguments in REPORT.md license the theorem-level conclusions.
"""
from fractions import Fraction as F
import importlib.util
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import null_space
from scipy.integrate import solve_ivp

import research as r


def source_certifier():
    path=Path(__file__).resolve().parent.parent/'fin_replication_consistency'/'certify.py'
    spec=importlib.util.spec_from_file_location('fin_geometry_enclosures',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def legacy(cyclic=True):
    def weight(i,j):
        if i==j:return 0.
        d=min(abs(i-j),12-abs(i-j)) if cyclic else abs(i-j)
        return 4*math.log(2)*math.cos(math.pi*d/4+math.pi/6)/(1+.01*d)
    return np.array([[weight(i,j) for j in range(12)] for i in range(12)])


def determinant_interval(low,high,scale):
    n=len(low);lower=0;upper=0
    for perm in itertools.permutations(range(n)):
        inversions=sum(perm[i]>perm[j] for i in range(n) for j in range(i+1,n))
        lo=hi=1
        for i,j in enumerate(perm):
            products=[lo*int(low[i,j]),lo*int(high[i,j]),hi*int(low[i,j]),hi*int(high[i,j])]
            lo,hi=min(products),max(products)
        if inversions%2:lower-=hi;upper-=lo
        else:lower+=lo;upper+=hi
    return F(lower,scale**n),F(upper,scale**n)


def rank_witness(K,low,high,scale):
    best=None
    for tail in itertools.combinations(range(1,12),5):
        I=(0,)+tail;J=tuple(x for x in range(12) if x not in I)
        singular=np.linalg.svd(K[np.ix_(I,J)],compute_uv=False)[-1]
        if best is None or singular>best[0]:best=(singular,I,J)
    sigma,I,J=best
    lo,hi=determinant_interval(low[np.ix_(I,J)],high[np.ix_(I,J)],scale)
    assert lo>0 or hi<0
    return dict(rows=I,columns=J,smallest_singular_value=float(sigma),
        determinant_interval=[str(lo),str(hi)],real_covariance_rank_lower_bound=6,
        complex_density_rank_lower_bound=3)


def skew_basis(n):
    result=[]
    for i,j in itertools.combinations(range(n),2):
        X=np.zeros((n,n));X[i,j]=1;X[j,i]=-1;result.append(X)
    return result


def tangent_ranks(K):
    n=len(K);basis=skew_basis(n)
    comm=[X@K-K@X for X in basis]
    full_orbit=np.array([C.ravel() for C in comm]).T
    diag_map=np.array([np.diag(C) for C in comm]).T
    u=np.ones(n)/math.sqrt(n)
    fix=np.array([X@u for X in basis]).T
    Z=null_space(fix)
    restricted_orbit=full_orbit@Z
    restricted_diag=diag_map@Z
    return dict(full_orbit_rank=int(np.linalg.matrix_rank(full_orbit,tol=1e-10)),
        full_diagonal_rank=int(np.linalg.matrix_rank(diag_map,tol=1e-10)),
        fixed_uniform_group_dimension=Z.shape[1],
        fixed_uniform_orbit_rank=int(np.linalg.matrix_rank(restricted_orbit,tol=1e-10)),
        fixed_uniform_diagonal_rank=int(np.linalg.matrix_rank(restricted_diag,tol=1e-10)))


def zero_diagonal_basis(S):
    """Numerical realization of the trace-zero Rayleigh-basis induction."""
    n=len(S)
    if n==1 or np.linalg.norm(S)<1e-13:return np.eye(n)
    values,vectors=np.linalg.eigh(S)
    lo,hi=values[0],values[-1]
    if not lo<0<hi:raise ValueError('Nonzero traceless compression must straddle zero')
    v=math.sqrt(hi/(hi-lo))*vectors[:,0]+math.sqrt(-lo/(hi-lo))*vectors[:,-1]
    B=np.column_stack((v,null_space(v.reshape(1,n))))
    transformed=B.T@S@B
    tail=zero_diagonal_basis(transformed[1:,1:])
    block=np.eye(n);block[1:,1:]=tail
    return B@block


def covariance_rank_minimizers(K):
    n=len(K);fourier=np.exp(2j*np.pi*np.outer(np.arange(n),np.arange(n))/n)/math.sqrt(n)
    eigen=np.fft.fft(K[0]).real
    gamma=-1/(n*eigen[n//2])
    c=1/n+gamma*eigen
    states=[]
    for signs in itertools.product((-1,1),repeat=5):
        rho=c[0]*np.outer(fourier[:,0],fourier[:,0].conj())
        for k,sign in enumerate(signs,1):
            idx=k if sign==1 else n-k
            rho+=2*c[k]*np.outer(fourier[:,idx],fourier[:,idx].conj())
        states.append(rho)
    residual=max(np.linalg.norm(r.off(rho.real)-gamma*K) for rho in states)
    ranks=[int(np.linalg.matrix_rank(rho,tol=1e-11)) for rho in states]
    assert residual<1e-12 and set(ranks)=={6}
    return gamma,c,states,dict(gamma=gamma,count=len(states),ranks=sorted(set(ranks)),
        fixed_point_residual=float(residual),minimal_entropy=float(-sum(x*math.log(x) for x in
            [c[0]]+[2*c[k] for k in range(1,6)] if x>0)))


def positive_circulant_permutations(K):
    n=len(K);values=np.fft.fft(K[0]).real;accepted=[]
    for order in itertools.permutations(range(1,6)):
        spectrum=values.copy()
        for k,index in enumerate(order,1):spectrum[k]=spectrum[n-k]=values[index]
        weights=np.fft.ifft(spectrum).real
        if weights[1:].min()>1e-12:
            accepted.append(dict(order=order,min_weight=float(weights[1:].min()),weights=weights[1:7].tolist()))
    return accepted


def certify_circulant_assignments():
    cert=source_certifier()
    values=[tuple(F(x) for x in I) for I in r.certify_strict_spectrum()['eigenvalue_intervals']]
    h3=cert.scale(F(1,2),cert.root_interval(3,2))
    one=(F(1),F(1));half=(F(1,2),F(1,2));zero=(F(0),F(0))
    table=[one,h3,half,zero,cert.scale(-1,half),cert.scale(-1,h3),
           cert.scale(-1,one),cert.scale(-1,h3),cert.scale(-1,half),zero,half,h3]
    accepted=[];rejected=[]
    for order in itertools.permutations(range(1,6)):
        weights=[]
        for d in range(1,7):
            value=cert.add(values[0],cert.scale((-1)**d,values[6]))
            for k,index in enumerate(order,1):
                value=cert.add(value,cert.scale(2,cert.multiply(values[index],table[(k*d)%12])))
            weights.append(cert.scale(F(1,12),value))
        if all(lo>0 for lo,hi in weights):
            accepted.append(dict(order=order,min_weight_lower=str(min(lo for lo,hi in weights))))
        else:
            negative=[(d+1,hi) for d,(lo,hi) in enumerate(weights) if hi<0]
            if not negative:raise AssertionError('Unresolved interval sign in assignment census')
            d,upper=min(negative,key=lambda item:item[1])
            rejected.append(dict(order=order,negative_distance=d,negative_upper=str(upper)))
    assert len(accepted)==2 and len(rejected)==118
    return dict(accepted=accepted,rejected=rejected,unresolved=0)


def regular_isospectral_witness(K,time=.002):
    """A curve of equilibrium points, NOT a trajectory of the learning ODE."""
    n=len(K);basis=np.array(skew_basis(n));ones=np.ones(n)
    fixed=np.array([X@ones for X in basis]).T[:n-1]
    seed=np.array([(k%7)-3 for k in range(len(basis))],float)
    def tangent(M,h):
        comm=basis@M-M@basis
        diagonal=np.array([np.diag(C) for C in comm]).T[:n-1]
        constraints=np.vstack((fixed,diagonal))
        x=h-constraints.T@np.linalg.solve(constraints@constraints.T,constraints@h)
        X=np.einsum('k,kij->ij',x,basis)
        return X@M-M@X
    seed/=np.linalg.norm(tangent(K,seed))
    solution=solve_ivp(lambda t,y:tangent(y.reshape(n,n),seed).ravel(),
        (0,time),K.ravel(),rtol=1e-11,atol=1e-13)
    if not solution.success:raise RuntimeError(solution.message)
    other=solution.y[:,-1].reshape(n,n)
    w0=K[0,1:7]
    offdiag=other[~np.eye(n,dtype=bool)]
    new_weight=float(np.max(np.min(abs(offdiag[:,None]-w0[None,:]),axis=1)))
    assert new_weight>1e-5
    return other,dict(parameter=time,kernel_displacement=float(np.linalg.norm(other-K)),
        diagonal_error=float(np.max(abs(np.diag(other)))),
        row_sum_error=float(np.max(abs(other.sum(axis=1)-K.sum(axis=1)))),
        eigenvalue_error=float(np.max(abs(np.linalg.eigvalsh(other)-np.linalg.eigvalsh(K)))),
        minimum_weight=float(offdiag.min()),new_weight_distance_from_original_values=new_weight,
        scope='Numerical representative of an analytically proved 39-dimensional local equilibrium manifold; not learning evolution.')


def run():
    K=r.strict();cert=source_certifier();w=cert.strict_weights();v=cert.legacy_weights()
    rank={}
    for name,M,weights,cycle in [('strict',K,w,True),('legacy_cycle',legacy(),v,True),
                                 ('legacy_line',legacy(False),v,False)]:
        lo,hi=cert.integer_matrix(weights,cycle)
        rank[name]=rank_witness(M,lo,hi,cert.SCALE)
    gamma,c,states,minimum=covariance_rank_minimizers(K)
    real=np.eye(12)/12+gamma*K
    reflection=np.eye(12)[[(-i)%12 for i in range(12)]]
    pair_error=max(np.linalg.norm(reflection@states[i]@reflection.T-states[-1-i]) for i in range(32))
    assert pair_error<1e-12
    ranks=tangent_ranks(K)
    assert ranks==dict(full_orbit_rank=61,full_diagonal_rank=11,fixed_uniform_group_dimension=55,
                      fixed_uniform_orbit_rank=50,fixed_uniform_diagonal_rank=11)
    # A distinct positive circulant representative: relabel vertices by a unit mod12.
    perm=[5*i%12 for i in range(12)];other=K[np.ix_(perm,perm)]
    assert np.linalg.norm(other-K)>1 and other[other>0].min()>0
    g=.05;rho=np.eye(12)/12+g*K;rho2=np.eye(12)/12+g*other
    bound=-(np.trace(rho@rho)-1/12)/(2*g)
    values=[r.lyapunov(M,R,g) for M,R in [(K,rho),(other,rho2)]]
    assert max(abs(value-bound) for value in values)<1e-13
    circles=positive_circulant_permutations(K)
    assignment_certificate=certify_circulant_assignments()
    witness,witness_data=regular_isospectral_witness(K)
    stationary_rho=np.eye(12)/12+.05*witness
    dK,dr=r.projected_rhs(witness,stationary_rho,.2,.05)
    witness_data['learning_derivative_residual']=float(np.linalg.norm(dK)+np.linalg.norm(dr))
    # Fast-learning counterexample to invariance of the Markov-rate cone.
    initial=np.zeros(12,complex);initial[0]=1/math.sqrt(2);initial[1]=-1/math.sqrt(2)
    state=np.outer(initial,initial.conj());rate=1000.;g=.05
    initial_full=np.concatenate((K.ravel().astype(complex),state.ravel()))
    def fast_rhs(theta,y):
        M=y[:144].reshape(12,12).real;R=y[144:].reshape(12,12)
        update,motion=r.projected_rhs(M,R,rate,g)
        return np.concatenate(((update/rate).ravel().astype(complex),(motion/rate).ravel()))
    sol=solve_ivp(fast_rhs,(0,1),initial_full,rtol=1e-11,atol=1e-13)
    if not sol.success:raise RuntimeError(sol.message)
    edge=float(sol.y[1,-1].real)
    analytic_upper=F(-3233,80000)+F(19,3000)
    assert edge<0 and edge<float(analytic_upper)
    cone_exit=dict(eta=rate,gamma=g,dimensionless_model_time=.001,K01=edge,
        analytic_upper_bound=str(analytic_upper),scope='Exact inequality proves negativity; integration is a check.')
    # Learning can lose regularity before any edge crosses zero.
    early=solve_ivp(fast_rhs,(0,.1),initial_full,rtol=1e-11,atol=1e-13)
    if not early.success:raise RuntimeError(early.message)
    earlyK=early.y[:144,-1].reshape(12,12).real
    degree=earlyK.sum(axis=1);Aearly=np.diag(degree)-earlyK
    duality=dict(positive_min_weight=float(earlyK[~np.eye(12,dtype=bool)].min()),
        degree_spread=float(np.ptp(degree)),
        Hamiltonian_Laplacian_commutator=float(np.linalg.norm(earlyK@Aearly-Aearly@earlyK)),
        regular_witness_intertime_commutator=float(np.linalg.norm(K@witness-witness@K)))
    assert duality['positive_min_weight']>0 and duality['Hamiltonian_Laplacian_commutator']>1e-3
    # Sparse pure stationary point admitted by positivity projection.
    signs=np.array([1.]*6+[-1.]*6)/math.sqrt(12)
    pure=np.outer(signs,signs);clipped=np.maximum(r.off(pure)/g,0.)
    residual=r.off(pure)-g*clipped
    projected=np.where(clipped>0,residual,np.maximum(residual,0.))
    assert np.linalg.norm(projected)<1e-14 and np.linalg.norm(clipped@pure-pure@clipped)<1e-13
    return dict(completed_rounds=list(range(11,21)),programs=[f'ST{i}' for i in range(8631,8641)],
        off_diagonal_rank_witnesses=rank,stationary_minimum_rank=minimum,
        reflection_pair_residual=float(pair_error),
        isospectral_minimizer_tangent_ranks=ranks,local_minimizer_dimension=50,
        local_regular_minimizer_dimension=39,
        distinct_equal_minima=dict(kernel_difference=float(np.linalg.norm(other-K)),
            functional_values=[float(x) for x in values],purity_bound=float(bound)),
        regular_isospectral_witness=witness_data,
        positive_circulant_assignment_scan=dict(count=len(circles),assignments=circles,
            scope='Numerical discovery, followed by the separate exact interval census.'),
        exact_circulant_assignment_certificate=assignment_certificate,
        positive_cone_exit=cone_exit,
        adaptive_duality_test=duality,
        projected_cone_sparse_equilibrium=dict(component_sizes=[6,6],density_rank=1,
            projected_learning_residual=float(np.linalg.norm(projected)),
            state_commutator=float(np.linalg.norm(clipped@pure-pure@clipped))))


if __name__=='__main__':
    result=run();Path(__file__).with_name('geometry_results.json').write_text(json.dumps(result,indent=2)+'\n')
    summary={k:v for k,v in result.items() if k!='exact_circulant_assignment_certificate'}
    summary['exact_circulant_assignment_counts']={k:len(result['exact_circulant_assignment_certificate'][k])
                                                for k in ['accepted','rejected']}
    print(json.dumps(summary,indent=2))
