"""ST8579--ST8590: replication consistency and the scope of FIN coupling laws.

Proofs are in report.tex. Floating checks and exhaustive finite scans are
labelled separately from analytic statements. No old research file is changed.
"""
from __future__ import annotations
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import expm_multiply


def strict():
    n=12
    return np.array([[0. if i==j else
        math.cos(.18575*min(abs(i-j),n-abs(i-j))+.16250)/
        (1+min(abs(i-j),n-abs(i-j))**1.8)
        for j in range(n)] for i in range(n)])


def legacy(cyclic=True):
    n=12
    def value(i,j):
        if i==j: return 0.
        d=min(abs(i-j),n-abs(i-j)) if cyclic else abs(i-j)
        return 4*math.log(2)*math.cos(math.pi*d/4+math.pi/6)/(1+.01*d)
    return np.array([[value(i,j) for j in range(n)] for i in range(n)])


def lap(W):
    return np.diag(W.sum(axis=1))-W


def finite_extension(W,N,c):
    """Antipodal common events between opposite parity classes.

    n must be divisible by four. Each singleton antipodal rate is reduced by
    c times the number of opposite-parity partners. Other rates are unchanged.
    Rejects genuinely negative rates rather than clipping them into validity.
    """
    n=len(W)
    if n%4: raise ValueError('Antipode must preserve parity')
    states=list(itertools.product(range(n),repeat=N))
    ix={x:i for i,x in enumerate(states)}
    rows=[];cols=[];vals=[]
    for row,x in enumerate(states):
        exits={}
        for k in range(N):
            opposite=sum((x[k]-x[l])%2 for l in range(N) if l!=k)
            for j in range(n):
                if j==x[k]: continue
                rate=W[x[k],j]-(c*opposite if j==(x[k]+n//2)%n else 0.)
                if rate < 0:
                    raise ValueError(f'Negative singleton rate {rate} at {x}')
                y=list(x);y[k]=j; target=ix[tuple(y)]
                exits[target]=exits.get(target,0.)+rate
        for k,l in itertools.combinations(range(N),2):
            if (x[k]-x[l])%2:
                y=list(x);y[k]=(y[k]+n//2)%n;y[l]=(y[l]+n//2)%n
                target=ix[tuple(y)]
                exits[target]=exits.get(target,0.)+c
        exits[row]=-sum(exits.values())
        for col,value in exits.items():
            if value:
                rows.append(row);cols.append(col);vals.append(value)
    return coo_matrix((vals,(rows,cols)),shape=(len(states),len(states))).tocsr(),states


def projection(states,fewer):
    ix={x:i for i,x in enumerate(fewer)}
    m=len(fewer[0]);rows=np.arange(len(states))
    cols=[ix[x[:m]] for x in states]
    return coo_matrix((np.ones(len(states)),(rows,cols)),
                     shape=(len(states),len(fewer))).tocsr()


def sparse_max(M):
    return float(np.max(abs(M.data))) if M.nnz else 0.


def pair_rates(G,states):
    """N=2 common-change rates b_ij for every pair of initial vertices."""
    n=max(max(x) for x in states)+1
    B=np.zeros((n,n))
    for r,x in enumerate(states):
        row=G.getrow(r)
        B[x[0],x[1]]=sum(v for j,v in zip(row.indices,row.data)
                       if states[j][0]!=x[0] and states[j][1]!=x[1])
    return B


def finite_bound(exit_rate,diagonal_budget,m):
    if m<1 or exit_rate<0 or diagonal_budget<0:
        raise ValueError('Invalid rate, diagonal budget or clone count')
    if exit_rate==0:return 0.
    return exit_rate/(2*m)*(1+math.sqrt(1+4*m*(m-1)*diagonal_budget/exit_rate))


def hopfield_scan(W):
    states=np.array(list(itertools.product((-1.,1.),repeat=len(W))))
    fields=states@W.T
    updated=np.where(fields>0,1.,np.where(fields<0,-1.,states))
    # States follow lexicographic binary order.
    powers=2**np.arange(len(W)-1,-1,-1)
    successor=((updated+1)/2@powers).astype(int)
    energies=-.5*np.sum(states*fields,axis=1)
    fixed=np.flatnonzero(successor==np.arange(len(states)))
    visited=np.zeros(len(states),bool);cycle_lengths=[]
    for start in range(len(states)):
        if visited[start]:continue
        path=[];local={};x=start
        while not visited[x] and x not in local:
            local[x]=len(path);path.append(x);x=successor[x]
        if x in local:cycle_lengths.append(len(path)-local[x])
        visited[path]=True
    return dict(fixed_points=len(fixed),cycle_length_histogram={str(k):cycle_lengths.count(k)
        for k in sorted(set(cycle_lengths))},minimum_abs_field=float(np.min(abs(fields))),
        lowest_energy=float(energies.min()),
        numerical_ground_state_count=int(np.count_nonzero(abs(energies-energies.min())<1e-10)),
        largest_synchronous_energy_increase=float(np.max(energies[successor]-energies)))


def run():
    W=strict();s=float(W.sum(axis=1)[0]);a=float(W[0,6]);c=a/2
    G2,x2=finite_extension(W,2,c)
    G3,x3=finite_extension(W,3,c)
    P=projection(x3,x2)
    B=pair_rates(G2,x2)
    assert np.max(abs(np.diag(B)))==0 and B[0,1]>0
    proj_error=sparse_max(G3@P-P@G2)
    assert proj_error<1e-12
    symmetry=sparse_max(G2-G2.T)
    assert symmetry<1e-12
    singleton=np.array([[int(x[0]==i) for i in range(12)] for x in x2],float)
    marginal_error=float(np.max(abs(G2@singleton-singleton@(-lap(W)))))
    assert marginal_error<1e-12
    ix={x:i for i,x in enumerate(x2)}
    symmetry_checks=[]
    for sign in [-1,1]:
        for shift in range(12):
            perm=[ix[tuple((sign*v+shift)%12 for v in x)] for x in x2]
            symmetry_checks.append(sparse_max(G2[perm,:][:,perm]-G2))
    assert max(symmetry_checks)<1e-12
    tables=[]
    for delta in [0.,1e-8,1e-5,.01]:
        tables.append(dict(delta=delta,bounds=[dict(m=m,bound=finite_bound(s,delta,m))
             for m in [1,2,4,16,64,1024]],infinite_bound=math.sqrt(s*delta)))
    minimum_invalid_rate=a-3*c
    assert minimum_invalid_rate<0
    # Explicitly test rejection at small n, not an enormous invalid strict matrix.
    small=np.ones((4,4))-np.eye(4)
    try:finite_extension(small,4,.5)
    except ValueError: rejected=True
    else:rejected=False
    assert rejected
    Gi,_=finite_extension(W,2,0.)
    t=.4
    initial=np.zeros(len(x2)); initial[ix[(0,0)]]=1
    p=expm_multiply(t*G2.T,initial)
    p0=expm_multiply(t*Gi.T,initial)
    clone_finite_tv=float(np.sum(abs(p-p0))/2)
    assert clone_finite_tv>1e-7
    scans={name:hopfield_scan(V) for name,V in
           [('strict_cycle',W),('legacy_cycle',legacy()),('legacy_line',legacy(False))]}
    V=legacy(False);largest=float(np.linalg.eigvalsh(V)[-1])
    crit=dict(EI_ratio=float(V[V>0].sum()/abs(V[V<0].sum())),
              max_eigenvalue=largest,stable_gain=.5/largest,unstable_gain=2/largest,
              stable_max_real_eigenvalue=-.5,unstable_max_real_eigenvalue=1.)
    # Stronger infinite-replication Gram inequality supersedes the one-group
    # square-root bound above when both-origin replication is available.
    Bcom=s*np.eye(12)+W
    incidence=[]
    for i,j in itertools.combinations(range(12),2):
        v=np.zeros(12);v[i]=v[j]=math.sqrt(W[i,j]);incidence.append(v)
    incidence=np.array(incidence).T
    cp_error=float(np.linalg.norm(incidence@incidence.T-Bcom))
    assert cp_error<1e-12
    # A maximal-activity family, not selected even by the full activity matrix.
    R=W/s
    Gproduct=s*(np.kron(R,R)-np.eye(144))
    Gshift=np.zeros((144,144))
    for row,x in enumerate(x2):
        for r in range(1,12):
            target=ix[tuple((v+r)%12 for v in x)]
            Gshift[row,target]+=W[0,r]
        Gshift[row,row]=-s
    maximum=[]
    for eta in [.25,.75]:
        G=(1-eta)*Gproduct+eta*Gshift
        budget=pair_rates(csr_matrix(G),x2)
        assert np.max(abs(budget-s))<1e-12
        assert np.max(abs(G@singleton-singleton@(-lap(W))))<1e-12
        maximum.append(dict(eta=eta,activity_error=float(np.max(abs(budget-s))),
            detailed_balance_error=float(np.linalg.norm(G-G.T)),
            transition_00_to_12=float(G[ix[(0,0)],ix[(1,2)]])))
    assert maximum[0]['transition_00_to_12']>maximum[1]['transition_00_to_12']>0
    alternating=np.array([(-1)**i for i in range(12)])
    negative_witness=float(alternating@B@alternating)
    assert negative_witness<0
    return dict(campaign='ST8579--ST8590',strict_row_sum=s,strict_antipodal_rate=a,
        common_rate=c,finite_horizon=3,common_pair_matrix=B.tolist(),
        three_to_two_projection_error=proj_error,singleton_error=marginal_error,
        detailed_balance_error=symmetry,dihedral_errors=symmetry_checks,
        uniform_stationary_error=float(np.max(abs(np.asarray(G3.sum(axis=0))))),
        dimensions=[len(x2),len(x3)],nnz=[G2.nnz,G3.nnz],
        four_copy_negative_singleton_rate=minimum_invalid_rate,
        four_copy_target_resolved_bound=a/3,finite_invalid_constructor_rejected=rejected,
        approximate_bounds=tables,clone_initial_time=t,clone_initial_finite_time_tv=clone_finite_tv,
        hopfield_scans=scans,criticality_counterexample=crit,
        infinite_gram_counterexample_quadratic_form=negative_witness,
        unsigned_incidence_cp_error=cp_error,maximal_activity_counterexamples=maximum,
        infinite_uniform_diagonal_budget_bound='b_ij <= delta; B_N <= binom(N,2)*delta',
        scope='Analytic theorems in report.tex; numerical scans are not interval certificates.')


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,
        default=Path(__file__).with_name('results.json'))
    args=parser.parse_args();result=run()
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    for k in ['campaign','strict_antipodal_rate','common_rate','three_to_two_projection_error',
              'singleton_error','detailed_balance_error','four_copy_negative_singleton_rate',
              'clone_initial_finite_time_tv','hopfield_scans','criticality_counterexample']:
        print(k,':',result[k])
