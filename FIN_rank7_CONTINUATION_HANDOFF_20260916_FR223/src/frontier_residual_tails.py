"""Post-handoff frontier: new global tail certificates for the unresolved 4D core.

These lemmas are *new continuation results* after R7P-128.  They do not alter
historical campaign claims in-place.  All proof decisions use the accepted
strict spectral intervals.

Compact variables:
    a = exp(-J3),
    s = exp(-3 J4/2),
    t = exp(-J5/2),
    y = exp(-2 J6).

The key tool is Courant--Fischer plus a low-rank limiting support.  Instead of
bounding the whole covariance trace, remove the single dominant line that
survives on the relevant compact face.  The second eigenvalue is then bounded
by the trace of the compressed covariance, and that trace only sees weights
that vanish with the chosen compact coordinate.
"""
from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
import json

from intervals import QI, sqrt_interval
from boundary_ising import strict_intervals

ROOT = Path(__file__).resolve().parents[1]


def _pair(I: QI):
    return [str(I.lo), str(I.hi)]


def spectral_data():
    L = strict_intervals()
    l3,l4,l5,l6 = L[3],L[4],L[5],L[6]
    sigma = (2*l3*(l4+l5)-l4*l5)/(24*l3)
    return L,l3,l4,l5,l6,sigma


def q_lower_bound_record():
    """Exact strengthened parity lower bound.

    With the exact parity partition sums
      Z+ = 2 e^J4 cosh(J3+J5) + 4 e^(-J4/2) cosh(J3-J5/2),
      Z- = 2 e^J4 + 4 e^(-J4/2) cosh(sqrt(3) J5/2),
    one has Z+ >= cosh(J3) Z- for J3,J4,J5>=0.

    Proof: after factoring 2e^(-J4/2), it suffices at J4=0 because
    A=cosh(J3+J5)-cosh(J3)>=0.  The J4=0 remainder is
      cosh(J3) A0(J5) + sinh(J3) B0(J5),
    B0=2 sinh(x/2)(cosh(x/2)-1)>=0, and
      A0=sum_{n>=3} [4^n+2-2*3^n]/4^n * x^(2n)/(2n)! >=0.
    Here h3=12 and h_{n+1}=4h_n+2*3^n-6>0.
    """
    return {
      'id':'FR1-q-lower-bound',
      'status':'EXACT_PROVED',
      'domain':'J3,J4,J5,J6>=0',
      'statement':'q/(1-q) >= cosh(J3), hence q >= cosh(J3)/(1+cosh(J3)).',
      'compact_form':'with a=exp(-J3), q >= (1+a^2)/(1+a)^2; equivalently q*(1+a)^2 >= 1+a^2.',
      'partition_sums':{
        'Z_even':'2 exp(J4) cosh(J3+J5)+4 exp(-J4/2) cosh(J3-J5/2)',
        'Z_odd':'2 exp(J4)+4 exp(-J4/2) cosh(sqrt(3) J5/2)',
        'odds':'q/(1-q)=exp(2J6) Z_even/Z_odd'
      },
      'proof_atoms':[
        'A=cosh(J3+J5)-cosh(J3)>=0, so increasing J4 only helps after factoring 2 exp(-J4/2).',
        'At J4=0 the difference equals cosh(J3) A0(J5)+sinh(J3) B0(J5).',
        'B0(x)=sinh(x)-2sinh(x/2)=2sinh(x/2)(cosh(x/2)-1)>=0.',
        'A0(x)=cosh(x)+2cosh(x/2)-1-2cosh(sqrt(3)x/2). Its x^2 and x^4 coefficients vanish; for n>=3 the x^(2n) coefficient is proportional to h_n=4^n+2-2*3^n>0.',
        'h_3=12 and h_(n+1)=4 h_n + 2*3^n-6 >0 for n>=3.'
      ],
      'equality_note':'At J5=0 the strengthened ratio is exactly cosh(J3); J6=0 then saturates the q lower bound.'
    }


def a_tail_record():
    """Global large-J3 tail via a rank-one limiting support."""
    _,l3,l4,l5,l6,sigma = spectral_data()
    # Sum of squared distances after projection orthogonal to line F0--F4=F8.
    C1 = l3 + 2*l6 + (l4*l5)/(l4+l5)  # six odd states, each weight <= a
    C2 = 2*l3 + (l4*l5)/(l4+l5)       # j=2,6,10, each weight <= a^2
    a0 = QI(F(1,30))
    bound = C1*a0 + C2*a0*a0
    gap = sigma-bound
    assert gap.lo>0
    return {
      'id':'FR1-large-J3-tail',
      'status':'INTERVAL_CERTIFIED',
      'domain':'a=exp(-J3)<=1/30; arbitrary J4,J5,J6>=0',
      'criterion':'lambda2(M4) <= C1*a + C2*a^2 < sigma_*',
      'C1_interval':_pair(C1),'C2_interval':_pair(C2),
      'a_threshold':'1/30',
      'upper_bound_at_threshold':_pair(bound),
      'sigma_interval':_pair(sigma),'strict_gap_interval':_pair(gap),
      'geometry':'At a=0 only feature values F0 and F4=F8 survive. Project to the orthogonal complement of F4-F0; the limiting covariance is zero on that 3D compression.',
      'weight_bound':'The six odd-label ratios are <=a and labels 2,6,10 are <=a^2; normalization is >=1.',
      'courant_fischer':'lambda2 is at most the largest eigenvalue of this 3D compression, which is at most its trace; variance is bounded by second moment about F0.',
      'exact_coefficients':{
        'C1':'lambda3 + 2 lambda6 + lambda4 lambda5/(lambda4+lambda5)',
        'C2':'2 lambda3 + lambda4 lambda5/(lambda4+lambda5)'
      }
    }


def s_tail_record():
    """Global large-J4 tail via a rank-two limiting support."""
    _,l3,l4,l5,l6,sigma = spectral_data()
    # u is the c3/c5 line occupied by labels 0,3,6,9 when s=0.
    # On u^perp the limiting variation is only c6, whose variance <= lambda6/12.
    Cs = (3*l3*l4 + 2*l3*l5 + 3*l4*l5)/(l3+l5)
    s0 = QI(F(1,128))
    base = l6/12
    bound = base + Cs*s0
    gap = sigma-bound
    assert gap.lo>0
    return {
      'id':'FR1-large-J4-tail',
      'status':'INTERVAL_CERTIFIED',
      'domain':'s=exp(-3J4/2)<=1/128; arbitrary J3,J5,J6>=0',
      'criterion':'lambda2(M4) <= lambda6/12 + Cs*s < sigma_*',
      'Cs_interval':_pair(Cs),'s_threshold':'1/128',
      'lambda6_over_12_interval':_pair(base),
      'upper_bound_at_threshold':_pair(bound),
      'sigma_interval':_pair(sigma),'strict_gap_interval':_pair(gap),
      'geometry':'At s=0 labels 0,3,6,9 lie in span{u,e6}, where u is the c3/c5 line. Compress to u^perp. The c6 variance is globally <=lambda6/12, and all remaining two compressed coordinates are constant on the s=0 support.',
      'weight_bound':'Every state that changes those remaining compressed coordinates carries a factor s; its normalized probability is <= its unnormalized anchor ratio <=s.',
      'exact_coefficient':'Cs=(3 lambda3 lambda4 + 2 lambda3 lambda5 + 3 lambda4 lambda5)/(lambda3+lambda5)'
    }


def t_tail_record():
    """Strong global large-J5 tail using the F0--F5=F7 dominant line."""
    _,l3,l4,l5,l6,sigma = spectral_data()
    rt3 = sqrt_interval(QI(3),50)
    D = 4*l3+9*l4-4*rt3*l5+7*l5+8*l6
    assert D.lo>0
    C1=(18*l3*l4-24*rt3*l3*l5+42*l3*l5+64*l3*l6-9*rt3*l4*l5+18*l4*l5+36*l4*l6+4*l5*l6)/(6*D)
    C2=(3*l3*l4+l3*l5+3*l4*l5+6*l4*l6+2*l5*l6)/D
    C3=(6*l3*l4+14*l3*l5+3*rt3*l4*l5+24*l4*l5+12*l4*l6+28*l5*l6)/(2*D)
    C4=2*(9*l3*l4+3*l3*l5+8*l3*l6+9*l4*l5+8*l5*l6)/(3*D)
    t0=QI(F(1,9))
    bound=C1*t0+C2*t0**2+C3*t0**3+C4*t0**4
    gap=sigma-bound
    assert gap.lo>0
    return {
      'id':'FR1-large-J5-projected-tail',
      'status':'INTERVAL_CERTIFIED',
      'domain':'t=exp(-J5/2)<=1/9; arbitrary J3,J4,J6>=0',
      'criterion':'lambda2(M4) <= C1*t + C2*t^2 + C3*t^3 + C4*t^4 < sigma_*',
      't_threshold':'1/9','sigma_interval':_pair(sigma),
      'coefficient_intervals':{'C1':_pair(C1),'C2':_pair(C2),'C3':_pair(C3),'C4':_pair(C4)},
      'upper_bound_at_threshold':_pair(bound),'strict_gap_interval':_pair(gap),
      'geometry':'As t->0 the slow labels 5 and 7 have the same C4 feature vector. Together with label 0 they occupy one affine line. Project orthogonally to F5-F0; those slow t^(2-sqrt(3)) states then contribute exactly zero to the compression.',
      'weight_groups':{
        't':'labels 2,10','t^2':'labels 3,9','t^3':'labels 1,11,4,8 (using 2+sqrt(3)>=3)','t^4':'label 6',
        'removed_line':'labels 0,5,7'
      },
      'courant_fischer':'lambda2 is bounded by the trace of the 3D compressed covariance; normalized probabilities are <= anchor-relative unnormalized weights because w0=1.',
      'improvement_over_R7P069':'Replaces the old certified t<=2^-11 tail by the much larger t<=1/9 tail.'
    }


def build():
    return {
      'continuation':'post-R7P-128 frontier 1',
      'q_lower_bound':q_lower_bound_record(),
      'large_J3_tail':a_tail_record(),
      'large_J4_tail':s_tail_record(),
      'large_J5_tail':t_tail_record(),
      'new_residual_outer_box':{
        'necessary_conditions_for_unresolved_point':['a=exp(-J3)>1/30','s=exp(-3J4/2)>1/128','t=exp(-J5/2)>1/9'],
        'note':'Certified boundary/face/local-cone regions from the completed campaign are additionally removed; this box is only a coarse outer hull of the new residual.'
      }
    }


def main():
    out=build()
    p=ROOT/'results/FR1_residual_tail_upgrade.json'
    p.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:(v.get('status') if isinstance(v,dict) else None) for k,v in out.items()},indent=2))

if __name__=='__main__':
    main()
