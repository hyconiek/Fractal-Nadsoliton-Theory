"""R7P-036: interval Krawczyk certificate for a stationary full-X7 index-2 witness at g=5.

All proof decisions use exact rational interval endpoints. The only floating operation is
creation of a rationalized preconditioner from a numerical midpoint Jacobian; after
rationalization, the Krawczyk inclusion is checked with QI arithmetic.
"""
from __future__ import annotations
from fractions import Fraction as F
import json
from pathlib import Path
import numpy as np

from .intervals import QI, exp_interval

ROOT=Path(__file__).resolve().parents[1]
AUDIT=ROOT/'inputs/fin_handoff_audit/results.json'
G=F(5)
CENTER_J=F('1.244948195164464')
CENTER_K=F('0.7795556209881978')
RADIUS=F(1,10_000_000)  # 1e-7


def _laplacian_intervals():
    data=json.loads(AUDIT.read_text())['exact']['laplacian_intervals']
    return [QI(F(lo),F(hi)) for lo,hi in data]


def hyp(J:QI):
    ep=exp_interval(J); em=exp_interval(-J)
    return (ep-em)/2,(ep+em)/2


def xy(J:QI,K:QI):
    sh,ch=hyp(J); z=exp_interval(-2*K); D=ch+z
    return sh/D,(ch-z)/D,(sh,ch,z,D)


def residual_and_jacobian(J:QI,K:QI,a:QI,b:QI):
    x,y,(sh,ch,z,D)=xy(J,K)
    xJ=(QI(1)+z*ch)/(D**2)
    xK=2*z*sh/(D**2)
    yJ=xK
    yK=4*z*ch/(D**2)
    fv=[J-G*a*x,K-G*b*y]
    jac=[[QI(1)-G*a*xJ,-G*a*xK],[-G*b*yJ,QI(1)-G*b*yK]]
    return fv,jac


def midpoint(I:QI): return (I.lo+I.hi)/2

def compact(I:QI,digits=30):
    """Outward-round an exact rational interval to compact denominator 10^digits."""
    Q=10**digits
    lo_num=(I.lo.numerator*Q)//I.lo.denominator
    hi_num=-((-I.hi.numerator*Q)//I.hi.denominator)
    return QI(F(lo_num,Q),F(hi_num,Q))

def qstr(I:QI):
    I=compact(I)
    return [str(I.lo),str(I.hi)]

def fstr(x:F,digits=30):
    Q=10**digits
    # For positive margins we store a conservative downward rational decimal.
    n=(x.numerator*Q)//x.denominator
    return str(F(n,Q))


def certify():
    L=_laplacian_intervals(); l3,l4,l5,l6=L[3],L[4],L[5],L[6]
    a=l3/6; b=l6/12
    XJ=QI(CENTER_J-RADIUS,CENTER_J+RADIUS)
    XK=QI(CENTER_K-RADIUS,CENTER_K+RADIUS)
    mJ,mK=QI(CENTER_J),QI(CENTER_K)
    Fm,Jm=residual_and_jacobian(mJ,mK,a,b)
    _,JX=residual_and_jacobian(XJ,XK,a,b)

    # Rational preconditioner from midpoint interval-Jacobian.
    M=np.array([[float(midpoint(q)) for q in row] for row in Jm],dtype=float)
    Cfloat=np.linalg.inv(M)
    C=[[F(str(Cfloat[i,j])) for j in range(2)] for i in range(2)]

    base=[]
    centers=[CENTER_J,CENTER_K]
    for i in range(2):
        s=QI(centers[i])
        for j in range(2): s=s-C[i][j]*Fm[j]
        base.append(s)
    B=[]
    for i in range(2):
        row=[]
        for k in range(2):
            s=QI(1 if i==k else 0)
            for j in range(2): s=s-C[i][j]*JX[j][k]
            row.append(s)
        B.append(row)
    delta=[QI(-RADIUS,RADIUS),QI(-RADIUS,RADIUS)]
    Kraw=[]
    for i in range(2):
        s=base[i]
        for k in range(2): s=s+B[i][k]*delta[k]
        Kraw.append(s)
    inclusion=[Kraw[i].lo>centers[i]-RADIUS and Kraw[i].hi<centers[i]+RADIUS for i in range(2)]
    assert all(inclusion)
    inclusion_margin=[min(Kraw[i].lo-(centers[i]-RADIUS),(centers[i]+RADIUS)-Kraw[i].hi) for i in range(2)]

    # At the true stationary root, x=J/(5a), y=K/(5b), so no additional
    # transcendental enclosure is needed for the Hessian-inertia check.
    x=XJ/(G*a); y=XK/(G*b)
    one=QI(1); invg=QI(F(1,5))

    # Exact symmetry block structure on the period-4 two-harmonic family:
    # singleton 3-sine; (3-cos,6-alt); and two independent (4,5) blocks.
    H3s=invg-a*(one-y)/2
    H33=invg-a*((one+y)/2-x*x)
    H66=invg-b*(one-y*y)
    det36=H33*H66-a*b*x*x*(one-y)**2
    H44=invg-l4/12
    H55=invg-l5/12
    det45=H44*H55-l4*l5*x*x/144

    # H3s>0 and H36 positive definite; each of the two disjoint H45 blocks
    # has determinant <0, hence exactly one negative eigenvalue. Therefore H7
    # has exactly two negative eigenvalues and no zero eigenvalue.
    assert H3s.lo>0
    assert H33.lo>0 and det36.lo>0
    assert H44.lo>0 and H55.lo>0 and det45.hi<0

    out={
      'id':'R7P-036-stationary-index2',
      'claim_id':'CLM-STATIONARY-INDEX2',
      'domain':'full X7 dual stationary system restricted to the exact invariant two-harmonic family at g=5',
      'quantifiers':'for the fixed strict spectral tuple enclosed by the accepted outward spectral intervals; Krawczyk inclusion is uniform over those intervals',
      'assumptions':['conditional active-gain dual with exact g=5','accepted strict spectral intervals','two-harmonic invariant family h_j=J cos(pi j/2)+K(-1)^j'],
      'proof_type':'parametric rational interval Krawczyk + exact symmetry block inertia',
      'inputs':['inputs/fin_handoff_audit/results.json','proofs/R7P-034_two_harmonic_reduction.md'],
      'conclusion':'there exists a unique stationary root in the stated box and its full H7 inertia is exactly (2 negative, 0 zero, 5 positive); hence unrestricted stationary index<=1 is false',
      'global_pass':False,
      'task':'R7P-036',
      'scientific_status':'INTERVAL_CERTIFIED',
      'gain':['5','5'],
      'spectral_inputs':{f'lambda{k}':qstr(L[k]) for k in [3,4,5,6]},
      'root_box':{'J':[str(CENTER_J-RADIUS),str(CENTER_J+RADIUS)],
                  'K':[str(CENTER_K-RADIUS),str(CENTER_K+RADIUS)]},
      'certified_intervals':[
          {'name':'J_root_box','lo':str(CENTER_J-RADIUS),'hi':str(CENTER_J+RADIUS)},
          {'name':'K_root_box','lo':str(CENTER_K-RADIUS),'hi':str(CENTER_K+RADIUS)}],
      'preconditioner_rational':[[fstr(v) for v in row] for row in C],
      'F_mid_interval':[qstr(q) for q in Fm],
      'Jacobian_box':[[qstr(q) for q in row] for row in JX],
      'Krawczyk_image':[qstr(q) for q in Kraw],
      'strict_inclusion_margin':[fstr(v) for v in inclusion_margin],
      'root_scope':'For every spectral tuple in the accepted strict intervals, one unique root lies in this J,K box; in particular the strict tuple has one unique root.',
      'stationary_lift':'The period-4 reflection-even two-harmonic family has zero residual in all omitted X7 coordinates identically, so the certified 2D root is a full seven-coordinate stationary point.',
      'hessian_blocks':{
          'x_interval':qstr(x),'y_interval':qstr(y),
          'H_3sin':qstr(H3s),
          'H_36_diag_3cos':qstr(H33),'H_36_diag_6alt':qstr(H66),'H_36_det':qstr(det36),
          'H_45_diag_4':qstr(H44),'H_45_diag_5':qstr(H55),'H_45_det_each':qstr(det45)},
      'inertia_conclusion':'H7 has exactly 2 negative, 0 zero, and 5 positive eigenvalues throughout the certified stationary root box for the strict spectral tuple.',
      'claim_consequence':'The unrestricted universal stationary-point conjecture index(H7)<=1 is false (already at exact supplied gain g=5). This does not decide a gain-restricted conjecture near g_eq.',
      'restricted_full_mismatch':{'H4_index':1,'H7_index':2,'extra_negative_sector':'(4sin,5sin)','statement':'The certified stationary point is index 1 in C4 but index 2 in full X7.'}
    }
    path=ROOT/'certificates/R7P-036_stationary_index2_witness.json'
    path.write_text(json.dumps(out,indent=2)+'\n')
    return out

if __name__=='__main__':
    d=certify()
    print(json.dumps(d,indent=2))
