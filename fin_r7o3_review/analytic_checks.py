"""Independent 12-label model reconstruction and cross-formulation controls.

Point tests below are controls, not substitutes for the all-cell Taylor proof.
Exact aggregation and the threshold comparison are separate analytic gates.
"""
from fractions import Fraction as F
from itertools import product
from collections import Counter
from pathlib import Path
import sys
import mpmath as mp

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from fin_r7o3_review.review import HERE,ROOT,SOURCE,TAU,load,lines,save,sha

iv=mp.iv
iv.dps=65


def I(x):
    x=F(x)
    return iv.mpf(x.numerator)/x.denominator


def bounds(x):
    def endpoint(t):
        sign,man,exponent,_=t
        return (-1 if sign else 1)*F(man)*F(2)**exponent
    return tuple(map(endpoint,x._mpi_))


def main():
    # Exact cosines in Q(sqrt(3)), encoded as (rational part, sqrt(3) part).
    c3=[(F(x),F(0)) for x in [1,0,-1,0]*3]
    c4=[(F(x),F(0)) for x in [1,F(-1,2),F(-1,2)]*4]
    c5=[(F(a),F(b)) for a,b in [(1,0),(0,F(-1,2)),(F(1,2),0),(0,0),
         (F(-1,2),0),(0,F(1,2)),(-1,0),(0,F(1,2)),(F(-1,2),0),(0,0),
         (F(1,2),0),(0,F(-1,2))]]
    groups=Counter()
    for j in range(12):
        groups[(1-c3[j][0],F(2,3)*(1-c4[j][0]),
                2*(1-c5[j][0]),-2*c5[j][1],F(1-(-1)**j,2))]+=1
    expected={(F(a),F(b),F(c),F(d),F(e)):m for a,b,c,d,e,m in
              [(0,0,0,0,0,1),(0,1,3,0,0,2),(2,0,4,0,0,1),
               (2,1,1,0,0,2),(1,0,2,0,1,2),(1,1,2,-1,1,2),(1,1,2,1,1,2)]}
    assert dict(groups)==expected
    # Both exponent brackets are proved by exact squared rational comparisons.
    assert F(19,11)**2<3<F(26,15)**2
    raw=load(ROOT/'fin_handoff_audit/results.json')['exact']['laplacian_intervals']
    from fin_projected_learning.research import certify_strict_spectrum
    W=[tuple(map(F,p)) for p in certify_strict_spectrum()['eigenvalue_intervals']]
    rebuilt=[['0','0']]+[[str(W[0][0]-b),str(W[0][1]-a)] for a,b in W[1:]]
    assert rebuilt==raw
    L=[I(a)+(I(b)-I(a))*iv.mpf([0,1]) for a,b in raw]
    sigma=(2*L[3]*(L[4]+L[5])-L[4]*L[5])/(24*L[3])
    assert bounds(sigma)[1]<TAU
    sqrt3=iv.sqrt(3)
    cosines=[[I(a)+I(b)*sqrt3 for a,b in row] for row in [c3,c4,c5]]
    scales=[iv.sqrt(L[k]/n) for k,n in [(3,6),(4,6),(5,6),(6,12)]]
    features=[[scales[k]*cosines[k][j] for k in range(3)]+[scales[3]*(-1)**j]
              for j in range(12)]
    certs=lines(SOURCE/'certificates/active_leaf_certificates.jsonl')
    sample=load(HERE/'rational_sample.json')['certificates']
    checked=0
    for record in sample:
        c=certs[record['index']]
        box=[tuple(map(F,p)) for p in c['cell']]
        B=[[I(F(x,c['basis_den'])) for x in row] for row in c['basis_num']]
        center=list(map(I,c['center_c']))
        z=[[sum((features[j][k]*B[k][a] for k in range(4)),I(0))-center[a]
            for a in range(3)] for j in range(12)]
        pts=list(product(*box))+[tuple((a+b)/2 for a,b in box)]
        for r,s,t,y in pts:
            J=[-iv.log(I(r))/2,-2*iv.log(I(s))/3,-2*iv.log(I(t)),-iv.log(I(y))/2]
            w=[iv.exp(sum((J[k]*(cosines[k][j]-1) for k in range(3)),I(0))
                      +J[3]*((-1)**j-1)) for j in range(12)]
            denominator=sum(w,I(0))
            for a in range(3):
                for b in range(a,3):
                    moment=sum((w[j]*z[j][a]*z[j][b] for j in range(12)),I(0))/denominator
                    lo,hi=bounds(moment)
                    cl,ch=map(F,record['moment_entry_enclosures'][a][b])
                    assert cl<=lo<=hi<=ch,(record['index'],a,b)
            checked+=1
    save('analytic_checks.json',dict(
        **{'pass':True},exact_12_label_aggregation=True,representatives=[0,4,6,2,3,5,1],
        multiplicities=[1,2,1,2,2,2,2],exact_sqrt3_bracket=['19/11','26/15'],
        sigma_interval=list(map(str,bounds(sigma))),sigma_upper=str(bounds(sigma)[1]),
        current_spectral_provider_recomputed=True,
        provider_sha256={str(p.relative_to(ROOT)):sha(p) for p in
                         [ROOT/'fin_projected_learning/research.py',ROOT/'fin_replication_consistency/certify.py']},
        direct_12_label_points=checked,moment_entries_checked=6*checked,
        sample_indices=[r['index'] for r in sample],
        point_check_scope='Independent exponential 12-label formula at all corners and midpoint of selected cells; diagnostic controls, not global coverage',
        proof_scope='Exact aggregation and threshold comparison; global Taylor calculus is inspected separately and every leaf is freshly replayed'))
    print('Exact aggregation and sigma < tau PASS;',checked,'independent 12-label points PASS',flush=True)


if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    main()
