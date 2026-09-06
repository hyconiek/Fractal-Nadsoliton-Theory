"""Active thirty-round FIN goal: exact composition and hierarchy tests.

No PDF, plots or typesetting. Analytic arguments are kept in REPORT.md.
"""
from __future__ import annotations
import itertools
import json
import math
from fractions import Fraction as F
from pathlib import Path

import numpy as np


def strict():
    return np.array([[0. if i==j else math.cos(.18575*min(abs(i-j),12-abs(i-j))+.1625)/
        (1+min(abs(i-j),12-abs(i-j))**1.8) for j in range(12)] for i in range(12)])


def lap(W):return np.diag(W.sum(axis=1))-W


def activity_counterexample():
    v=[0,1,2,3,2,1]*2
    H=[[v[(k-i)%12] for k in range(12)] for i in range(12)]
    norm=sum(x*x for x in v)
    B=[[F(sum(H[i][k]*H[j][k] for k in range(12)),norm)
        for j in range(12)] for i in range(12)]
    assert B[0][1]==F(16,19) and B[0][2]==F(11,19)
    facet=B[0][2]-B[0][1]-B[1][2]+1
    assert facet==F(-2,19)
    values=[]
    for x,y,z in itertools.product((0,1),repeat=3):
        value=x*z-x*y-y*z+y
        assert value>=0
        values.append(value)
    return B,H,dict(gram_denominator=norm,b01=str(B[0][1]),b12=str(B[1][2]),
        b02=str(B[0][2]),violated_facet=str(facet),boolean_facet_values=values)


def pair_generator(W,normalized_activity):
    n=len(W);f=lambda i:(i+n//2)%n
    states=list(itertools.product(range(n),repeat=2));ix={x:i for i,x in enumerate(states)}
    G=np.zeros((n*n,n*n))
    for row,(i,j) in enumerate(states):
        common=W[i,f(i)]*float(normalized_activity[i][j])
        for k in range(n):
            if k!=i:
                rate=W[i,k]-(common if k==f(i) else 0.)
                if rate < 0:raise ValueError('Negative rate')
                G[row,ix[(k,j)]]+=rate
            if k!=j:
                rate=W[j,k]-(common if k==f(j) else 0.)
                if rate < 0:raise ValueError('Negative rate')
                G[row,ix[(i,k)]]+=rate
        G[row,ix[(f(i),f(j))]]+=common
        G[row,row]=-G[row].sum()
    return G,states


def finite_difference_laws(order):
    """Two disjoint rational laws matching moments through the given order."""
    if order<1:raise ValueError('order must be positive')
    laws=[[],[]]
    for k in range(order+2):
        laws[k%2].append((F(-1)+F(2*k,order+1),F(math.comb(order+1,k),2**order)))
    return laws


def moment(law,order):return sum(weight*point**order for point,weight in law)


def wasserstein_one(mu,nu):
    nodes=sorted(set(x for x,_ in mu)|set(x for x,_ in nu))
    difference={x:F(0) for x in nodes}
    for x,w in mu:difference[x]+=w
    for x,w in nu:difference[x]-=w
    cdf=F(0);area=F(0)
    for i in range(len(nodes)-1):
        cdf+=difference[nodes[i]]
        area+=abs(cdf)*(nodes[i+1]-nodes[i])
    return area


def shifted_moment(law,p,d,order):
    return sum(w*(p+d*x)**order for x,w in law)


def recover_moments(values,p,d):
    recovered=[F(1)]
    for n,value in enumerate(values,1):
        lower=sum(F(math.comb(n,k))*p**(n-k)*d**k*recovered[k] for k in range(n))
        recovered.append((value-lower)/d**n)
    return recovered


def selfsimilar_moments(r,innovation,max_order):
    r=F(r);step=1-r;mom=[F(1)]
    for n in range(1,max_order+1):
        rhs=sum(F(math.comb(n,k))*r**k*step**(n-k)*mom[k]*moment(innovation,n-k)
                for k in range(n))
        mom.append(rhs/(1-r**n))
    return mom


def noise_curve(W):
    s=float(W[0].sum());R=W/s
    D=np.zeros_like(W);amplitude=R[0,2]/2
    for i in range(12):
        for j in range(12):
            d=min(abs(i-j),12-abs(i-j))
            if d==1:D[i,j]=amplitude
            elif d==2:D[i,j]=-amplitude
    assert np.max(abs(D.sum(axis=1)))<1e-15
    for u in [-1.,0.,1.]:
        candidate=R+u*D
        assert np.min(candidate)>=0 and np.max(abs(candidate.sum(axis=1)-1))<1e-14
    return s,R,D


def common_tick_generator(R,D,law,N,clock):
    result=np.zeros((len(R)**N,len(R)**N))
    for value,weight in law:
        single=R+float(value)*D
        tensor=single
        for _ in range(1,N):tensor=np.kron(tensor,single)
        result+=float(weight)*tensor
    return clock*(result-np.eye(len(result)))


def run():
    W=strict();B,H,exact=activity_counterexample();G,xs=pair_generator(W,B)
    P=np.array([[x[0]==i for i in range(12)] for x in xs],float)
    marginal=float(np.max(abs(G@P-P@(-lap(W)))))
    reversibility=float(np.max(abs(G-G.T)))
    assert marginal<1e-13 and reversibility<1e-13
    ix={x:i for i,x in enumerate(xs)}
    sym=[]
    for sign,shift in itertools.product((-1,1),range(12)):
        perm=[ix[tuple((sign*k+shift)%12 for k in x)] for x in xs]
        sym.append(float(np.max(abs(G[np.ix_(perm,perm)]-G))))
    assert max(sym)<1e-13
    lawA=[(F(-1,2),F(1,2)),(F(1,2),F(1,2))]
    lawB=[(F(-1,4),F(4,5)),(F(1),F(1,5))]
    assert [moment(lawA,k) for k in range(3)]==[moment(lawB,k) for k in range(3)]
    assert moment(lawB,3)-moment(lawA,3)==F(3,16)
    s,R,D=noise_curve(W)
    G2a=common_tick_generator(R,D,lawA,2,s)
    G2b=common_tick_generator(R,D,lawB,2,s)
    pair_error=float(np.max(abs(G2a-G2b)))
    assert pair_error<1e-13
    def target(law,N):return s*sum(float(w)*(R[0,1]+float(u)*D[0,1])**N for u,w in law)
    triple_difference=target(lawB,3)-target(lawA,3)
    predicted=s*float(F(3,16))*D[0,1]**3
    assert abs(triple_difference-predicted)<1e-15 and predicted>0
    moment_checks=[]
    for order in range(1,13):
        mu,nu=finite_difference_laws(order)
        assert sum(w for _,w in mu)==sum(w for _,w in nu)==1
        equal=[moment(mu,k)==moment(nu,k) for k in range(order+1)]
        gap=moment(mu,order+1)-moment(nu,order+1)
        formula=F((-1)**(order+1)*math.factorial(order+1),2**order)*F(2,order+1)**(order+1)
        assert all(equal) and gap==formula and gap!=0
        distance=wasserstein_one(mu,nu)
        assert distance==F(2,order+1)
        moment_checks.append(dict(order=order,matched_moments=len(equal),next_moment_gap=str(gap),
                                 latent_total_variation=1,wasserstein_one=str(distance)))
    p,d=F(1,3),F(1,20)
    values=[shifted_moment(lawB,p,d,n) for n in range(1,9)]
    recovered=recover_moments(values,p,d)
    assert recovered==[moment(lawB,n) for n in range(9)]
    innovationA=[(F(-1),F(1,2)),(F(1),F(1,2))]
    innovationB=[(F(-1),F(1,4)),(F(0),F(1,2)),(F(1),F(1,4))]
    ifsA=selfsimilar_moments(F(3,5),innovationA,10)
    ifsB=selfsimilar_moments(F(1,3),innovationB,10)
    assert ifsA[2]==ifsB[2]==F(1,4)
    assert ifsA[4]==F(35,272) and ifsB[4]==F(11,80)
    return dict(completed_rounds=list(range(1,11)),programs=[f'ST{i}' for i in range(8591,8601)],
        pair_gram_insufficiency=dict(**exact,strict_antipodal_rate=float(W[0,6]),
            strict_facet_gap=float(-2*W[0,6]/19),singleton_error=marginal,
            reversibility_error=reversibility,dihedral_errors=sym,
            gram_nonnegative_integer_factor=H),
        pair_law_nonuniqueness=dict(lawA=[(str(x),str(w)) for x,w in lawA],
            lawB=[(str(x),str(w)) for x,w in lawB],
            pair_generator_error=pair_error,third_moment_gap='3/16',
            triple_target_rate_A=target(lawA,3),triple_target_rate_B=target(lawB,3),
            triple_difference=triple_difference,exact_coefficient_prediction=predicted,
            clock=s,curve_amplitude=float(D[0,1])),
        arbitrary_finite_order_checks=moment_checks,
        all_order_recovery=dict(test_p=str(p),test_d=str(d),
            recovered_moments=[str(x) for x in recovered],
            raw_moment_noise_amplification=[str(d**(-n)) for n in [2,4,8]],
            strict_curve_noise_amplification=[float(D[0,1]**(-n)) for n in [2,4,8]]),
        selfsimilar_candidates=dict(bernoulli_contraction='3/5',
            bernoulli_moments=[str(x) for x in ifsA],
            ternary_contraction='1/3',ternary_moments=[str(x) for x in ifsB],
            same_second_moment='1/4',different_fourth_moments=['35/272','11/80']))


if __name__=='__main__':
    result=run()
    Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='pair_gram_insufficiency'},indent=2))
