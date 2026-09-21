"""Corrected proof-layer checks, independent of stored PASS booleans.

Exact-decimal phase fixture, true pi/square roots, interval Krawczyk and
interval congruence/LDL inertia. Original input packages are never rewritten.
"""
from fractions import Fraction as F
from pathlib import Path
import itertools
import json
import sys
import numpy as np
import mpmath as mp

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
SOURCE=ROOT/'FIN_rank7_CONTINUATION_HANDOFF_20260916_FR223'
mp.mp.dps=80;mp.iv.dps=70;iv=mp.iv

def I(x):
    f=F(str(x)) if not isinstance(x,F) else x
    return iv.mpf(f.numerator)/f.denominator

def endpoint(t):
    sign,man,exponent,_=t
    return (-1 if sign else 1)*F(man)*F(2)**exponent

def bounds(v): return tuple(endpoint(x) for x in v._mpi_)
def mid(v):
    a,b=bounds(v);return float((a+b)/2)
def serialize(v):
    a,b=bounds(v);q=10**30
    return [str(F((a*q).__floor__(),q)),str(F((b*q).__ceil__(),q))]

def phase_FH(phi,kind):
    amps=list(map(I,['0.1131879146','0.1698528641','0.2269339093']))
    z6=I('-0.3380663037');h=[];d=[];dd=[]
    for j in range(12):
        v=z6*((-1)**j)/iv.sqrt(12);row=[];row2=[]
        for a,k,p in zip(amps,[3,4,5],phi):
            angle=2*iv.pi*k*j/12+p;aa=a/iv.sqrt(3)
            v+=aa*iv.cos(angle);row.append(-aa*iv.sin(angle));row2.append(-aa*iv.cos(angle))
        h.append(v);d.append(row);dd.append(row2)
    if kind=='full':
        w=[iv.exp(v) for v in h];s=sum(w,I(0));p=[v/s for v in w]
        g=[sum((p[j]*d[j][a] for j in range(12)),I(0)) for a in range(3)]
        H=[[sum((p[j]*d[j][a]*d[j][b] for j in range(12)),I(0))-g[a]*g[b]
            +(sum((p[j]*dd[j][a] for j in range(12)),I(0)) if a==b else 0)
            for b in range(3)] for a in range(3)]
    else:
        # Fourier orthogonality makes M2 and M2^2 exactly phase independent.
        # Therefore phase derivatives of K4 equal those of M3/6+M4/24.
        weights=[v*v/2+v**3/6 for v in h]
        weights2=[v+v*v/2 for v in h]
        g=[sum((weights[j]*d[j][a] for j in range(12)),I(0))/12 for a in range(3)]
        H=[[sum((weights2[j]*d[j][a]*d[j][b]
                 +(weights[j]*dd[j][a] if a==b else 0) for j in range(12)),I(0))/12
            for b in range(3)] for a in range(3)]
    return g,H

def det3(A):
    return A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])-A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0])

def phase_root(rec,kind):
    center=[F(str(x)) for x in rec['phase']];radius=F(1,10**7)
    c=[I(x) for x in center]
    box=[I(x)+iv.mpf(['-0.0000001','0.0000001']) for x in center]
    g0,H0=phase_FH(c,kind);_,HB=phase_FH(box,kind)
    Hmid=np.array([[mid(x) for x in row] for row in H0])
    A=np.linalg.inv(Hmid);Ai=[[I(repr(x)) for x in row] for row in A]
    alo,ahi=bounds(det3(Ai));assert ahi<0 or alo>0
    offsets=[];contraction=[]
    delta=iv.mpf(['-0.0000001','0.0000001'])
    for a in range(3):
        v=-sum((Ai[a][j]*g0[j] for j in range(3)),I(0))
        row_norm=F(0)
        for b in range(3):
            e=I(int(a==b))-sum((Ai[a][k]*HB[k][b] for k in range(3)),I(0))
            elo,ehi=bounds(e);row_norm+=max(abs(elo),abs(ehi))
            v+=e*delta
        lo,hi=bounds(v);assert -radius<lo<=hi<radius
        assert row_norm<1;contraction.append(row_norm)
        offsets.append(v)
    # A rational change of basis near numerical eigenvectors is only a proposal;
    # nonsingularity and every LDL pivot sign are subsequently interval checked.
    _,Q=np.linalg.eigh(Hmid);Q=[[I(repr(x)) for x in row] for row in Q]
    qlo,qhi=bounds(det3(Q));assert qhi<0 or qlo>0
    H=[[sum((Q[i][a]*HB[i][j]*Q[j][b] for i in range(3) for j in range(3)),I(0))
        for b in range(3)] for a in range(3)]
    L=[[I(int(i==j)) for j in range(3)] for i in range(3)];D=[];signs=[]
    for k in range(3):
        pivot=H[k][k]-sum((L[k][j]*L[k][j]*D[j] for j in range(k)),I(0))
        lo,hi=bounds(pivot);assert hi<0 or lo>0
        signs.append(-1 if hi<0 else 1);D.append(pivot)
        for i in range(k+1,3):
            L[i][k]=(H[i][k]-sum((L[i][j]*L[k][j]*D[j] for j in range(k)),I(0)))/pivot
    index=signs.count(-1)
    assert index==rec['negative_index']
    return dict(center=[str(x) for x in center],radius=str(radius),
                contraction_upper=str(max(contraction)),preconditioner=[[repr(x) for x in row] for row in A],
                krawczyk_offsets=[serialize(x) for x in offsets],ldl_pivots=[serialize(x) for x in D],negative_index=index)

def phases():
    out={}
    for kind,path in [('quartic','results/R7P-089_quartic_roots.json'),('full','certificates/R7P-092_full_phase_roots.json')]:
        data=json.loads((SOURCE/path).read_text());roots=data['roots'];rows=[]
        for i,rec in enumerate(roots):
            rows.append(phase_root(rec,kind))
            if i%15==0: print(kind,i+1,flush=True)
        assert len(rows)==60
        # Pairwise torus separation: numerical locator distances are enormous
        # relative to the rational boxes; verify using rational 2pi enclosure.
        for a,b in itertools.combinations(rows,2):
            ca=list(map(F,a['center']));cb=list(map(F,b['center']))
            separated=False
            for x,y in zip(ca,cb):
                distances=[bounds(I(x-y)+n*2*iv.pi) for n in [-1,0,1]]
                if all(lo>2*F(a['radius']) or hi<-2*F(a['radius']) for lo,hi in distances):separated=True
            assert separated
        out[kind]=dict(count=60,roots=rows,scope='At least 60 distinct local roots of the exact decimal-amplitude fixture; no global exhaustion.')
        (HERE/'phase_recertification.json').write_text(json.dumps(out,indent=2)+'\n')

def globals_recheck():
    sys.path.insert(0,str(ROOT))
    from fin_projected_learning.research import certify_strict_spectrum
    from fin_handoff_audit.research import spectrum_intervals
    L=[iv.mpf([str(x.lo),str(x.hi)]) for x in spectrum_intervals()]
    original=json.loads((SOURCE/'inputs/fin_handoff_audit/results.json').read_text())['exact']['laplacian_intervals']
    assert original==[[str(x.lo),str(x.hi)] for x in spectrum_intervals()]
    nums=[836365227,1700465,10147081,26989852,24859253,13242432,9756607,13242432,24859253,26989852,10147081,1700465]
    assert sum(nums)==10**9 and min(nums)>0
    p=[I(F(x,10**9)) for x in nums]
    S={k:sum((p[j]*iv.cos(2*iv.pi*k*j/12) for j in range(12)),I(0)) for k in [3,4,5]}
    Q=sum((L[k]/6*S[k]**2 for k in [3,4,5]),I(0))+L[6]/12*sum((p[j]*(-1)**j for j in range(12)),I(0))**2
    D=sum((x*iv.log(12*x) for x in p),I(0));E=D-I('3.71835')*Q/2
    assert bounds(E)[1]<0
    eig=certify_strict_spectrum()
    delta=F(eig['minimum_intersector_gap_lower'])/20
    assert delta==F(27234855667,12500000000000)
    def angular_f(r):
        x=iv.sqrt(L[6]/12)*I(r);t=(iv.exp(2*x)-1)/(iv.exp(2*x)+1)
        return L[3]*(1+t)-L[6]*t/x
    fl=angular_f('0.41421132290');fh=angular_f('0.41421132293')
    assert bounds(fl)[1]<0 and bounds(fh)[0]>0
    x=iv.sqrt(L[6]/12)*I('0.41421132293');t=(iv.exp(2*x)-1)/(iv.exp(2*x)+1)
    lag=L[6]*t/(12*x)
    others=[L[k]/12-lag for k in [4,5]]
    assert all(bounds(v)[1]<0 for v in others)
    out=dict(energy_at_exact_g_3_71835=serialize(E),ratio=serialize(2*D/Q),
             historical_delta_replayed_from_original_W_provider=str(delta),
             angular_radius=['0.41421132290','0.41421132293'],
             angular_endpoint_function=[serialize(fl),serialize(fh)],
             angular_other_sector_upper_checks=[serialize(v) for v in others],
             lower_bracket_dependency='Existing ST448 plus A7<=A_full; not newly replayed full ST448 cover.')
    (HERE/'global_rechecks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':globals()[sys.argv[1]]()
