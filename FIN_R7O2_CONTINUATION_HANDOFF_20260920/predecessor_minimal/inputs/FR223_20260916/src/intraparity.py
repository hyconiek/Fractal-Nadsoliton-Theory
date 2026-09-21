"""R7P-057--059 exact intraparity reconstruction.

The two conditional parity laws are always derived from the SAME nonnegative
fields J3,J4,J5,J6.  J6 changes only the parity mixture weight; it cancels
from each conditional distribution.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import json, math
import sympy as sp

from intervals import QI, sqrt_interval
import boundary_ising as bi

ROOT=Path(__file__).resolve().parents[1]


def qij(I): return [str(I.lo),str(I.hi)]


def partition_formulas():
    """Exact unscaled conditional partition sums and shared-field maps."""
    return {
      'Z_even':'2*exp(J4)*cosh(J3+J5)+4*exp(-J4/2)*cosh(J3-J5/2)',
      'Z_odd':'2*exp(J4)+4*exp(-J4/2)*cosh(sqrt(3)*J5/2)',
      'q_even':'exp(J6)*Z_even/(exp(J6)*Z_even+exp(-J6)*Z_odd)',
      'q0':'Z_even/(Z_even+Z_odd)',
      'even_parameters':{'X':'exp(2*J3)','Y':'exp(3*J4/2)','Z':'exp(J5/2)'},
      'odd_parameters':{'Y':'exp(3*J4/2)','S':'exp(sqrt(3)*J5/2)','shared_relation':'S=Z**sqrt(3)'},
    }


def parity_weight_proof():
    """Exact proof skeleton for q0>=1/2 without a cooperativity assumption."""
    # Put A=J3, B=J5/2.  Since exp(J4)>=exp(-J4/2), it suffices to prove
    # F=cosh(A+2B)-1+2cosh(A-B)-2cosh(sqrt(3)B)>=0.
    # dF/dA at A=0 is sinh(2B)-2sinh(B)=2sinh(B)(cosh(B)-1)>=0,
    # and d2F/dA2>0, so F(A,B)>=F(0,B).  The latter has even-power
    # coefficients 4^n+2-2*3^n, zero for n=1,2 and positive for n>=3.
    coeff=[4**n+2-2*3**n for n in range(1,9)]
    return {
      'claim':'Z_even>=Z_odd for J3,J4,J5>=0; equality iff J3=J5=0',
      'reduction':'F(A,B)=cosh(A+2B)-1+2*cosh(A-B)-2*cosh(sqrt(3)*B)',
      'dF_dA_at_A0':'2*sinh(B)*(cosh(B)-1)>=0',
      'd2F_dA2':'cosh(A+2B)+2*cosh(A-B)>0',
      'F0_series_coefficients_formula':'4^n+2-2*3^n',
      'first_coefficients_n1_to_n8':coeff,
      'coefficient_induction':'a_2=0; a_(n+1)=4*a_n+2*3^n-6>0 for n>=2',
      'consequence':'q0>=1/2; for J6>=0, q_even>=q0>=1/2',
      'equality':'J3=J5=0 (arbitrary J4) at J6=0 for q0=1/2',
    }


def conditional_features_symbolic():
    """Actual C4 conditional feature points, with the constant k6 coordinate."""
    l3,l4,l5,l6=sp.symbols('l3 l4 l5 l6',positive=True)
    a3=sp.sqrt(l3/6);a4=sp.sqrt(l4/6);a5=sp.sqrt(l5/6);a6=sp.sqrt(l6/12)
    plus=sp.Matrix([
      [ a3, a4, a5, a6],
      [ a3,-a4/2,-a5/2,a6],
      [-a3, a4,-a5,a6],
      [-a3,-a4/2,a5/2,a6],
    ])
    minus=sp.Matrix([
      [0,a4,0,-a6],
      [0,-a4/2,sp.sqrt(3)*a5/2,-a6],
      [0,-a4/2,-sp.sqrt(3)*a5/2,-a6],
    ])
    return (l3,l4,l5,l6),plus,minus


def odd_domain_symbolics():
    """Exact (u,d) law and covariance for the odd conditional sector."""
    u,d,l4,l5=sp.symbols('u d l4 l5',real=True)
    p0=1-u;pp=(u+d)/2;pm=(u-d)/2
    C=sp.Matrix([
      [3*l4*u*(1-u)/8,-sp.sqrt(3*l4*l5)*(1-u)*d/8],
      [-sp.sqrt(3*l4*l5)*(1-u)*d/8,l5*(u-d**2)/8]
    ])
    return {
      'symbols':(u,d,l4,l5),'probabilities':(p0,pp,pm),'covariance':C,
      'domain':['0<=d<=u<=1','4*(1-u)^2 >= u^2-d^2'],
      'reconstruction':{
        'S2':'(u+d)/(u-d)',
        'Y2':'4*(1-u)^2/(u^2-d^2)',
        'J5':'log(S2)/sqrt(3)',
        'J4':'log(Y2)/3',
      },
      'boundary_interpretation':'d=u is p_minus=0; 4(1-u)^2=u^2-d^2 is J4=0.'
    }


def odd_supremum_proof():
    """Exact/strict-spectral proof of sup lambda1(C_minus)."""
    L=bi.strict_intervals();l4,l5=L[4],L[5]
    signs={
      'lambda5_minus_lambda4':qij(l5-l4),
      '9lambda4_minus_5lambda5':qij(9*l4-5*l5),
      '7lambda4_minus_3lambda5':qij(7*l4-3*l5),
    }
    assert (l5-l4).lo>0 and (9*l4-5*l5).lo>0 and (7*l4-3*l5).lo>0
    # A=3*l4, B=l5, N=8*Cminus, c=(A+B)/4=8*m.
    # mI-C PSD is checked by positive diagonals and determinant.  The determinant
    # is affine in y=d^2.  Its y derivative switches at uc=(3A-B)/(4A)<2/3.
    # On y=u^2 it is (A+B)^2(2u-1)^2/16.  On y=0 it is
    # (A+B-4Bu)(A(2u-1)^2+B)/16.  On the J4=0 lower boundary for u>=2/3,
    # mapping u=2/3+z/3 gives the four Bernstein coefficients below.
    l4s,l5s=sp.symbols('l4 l5',positive=True)
    bern=[
      (l4s+3*l5s)*(9*l4s-5*l5s)/144,
      (3*l4s+l5s)*(7*l4s-3*l5s)/144,
      5*(3*l4s+l5s)**2/144,
      (3*l4s+l5s)**2/16,
    ]
    return {
      'status':'INTERVAL_CERTIFIED_EXACT_FORMULA',
      'supremum':'(3*lambda4+lambda5)/32',
      'attainment':'closure point (p0,p_plus,p_minus)=(1/2,1/2,0)',
      'field_limit':'J4,J5->+infinity with J5=sqrt(3)*J4 asymptotically',
      'finite_attainment':False,
      'spectral_signs':signs,
      'determinant_piecewise':{
        'dy_sign_switch':'u_c=(9*lambda4-lambda5)/(12*lambda4) < 2/3',
        'u_le_uc_at_y_eq_u2':'(3*lambda4+lambda5)^2*(2*u-1)^2/16',
        'uc_to_2over3_at_y0':'(3*lambda4+lambda5-4*lambda5*u)*(3*lambda4*(2*u-1)^2+lambda5)/16',
        '2over3_to_1_lower_boundary_Bernstein':[str(sp.factor(x)) for x in bern],
      },
      'diagonal_bounds':['8*m-3*lambda4*u*(1-u) >= lambda5/4 >0',
                         '8*m-lambda5*(u-d^2) >= (9*lambda4-5*lambda5)/12 >0'],
      'uniqueness':'Only u=d=1/2 makes the determinant zero; elsewhere the PSD inequality is strict.'
    }


def dangerous_symbolics():
    u,y,l4,l5,sigma=sp.symbols('u y l4 l5 sigma',positive=True)
    c11=3*l4*u*(1-u)/8;c22=l5*(u-y)/8;off2=3*l4*l5*(1-u)**2*y/64
    det=sp.factor((sigma-c11)*(sigma-c22)-off2)
    dmin=sp.factor((8*sigma-l5*u)*(8*sigma-3*l4*u*(1-u))/(l5*(3*l4*(1-u)-8*sigma)))
    return u,y,l4,l5,sigma,det,dmin


def dangerous_interval_data():
    L=bi.strict_intervals();l3,l4,l5=L[3],L[4],L[5]
    sigma=(2*l3*(l4+l5)-l4*l5)/(24*l3)
    disc=1-32*sigma/(3*l4+l5);root=sqrt_interval(disc,36)
    uminus=(1-root)/2;uplus=(1+root)/2
    uc=1-8*sigma/(3*l4)
    # Uniform signs required for the cleared determinant equivalence on the full dangerous u interval.
    s1=8*sigma-l5*uplus
    s2=3*l4*(1-uplus)-8*sigma
    # Simple global smaller-eigenvalue bound: lambda_min <= trace/2 <= 3 l4/64+l5/24.
    small_gap=sigma-(3*l4/F(64)+l5/F(24))
    assert disc.lo>0 and disc.hi<1
    assert uplus.hi<uc.lo
    assert s1.lo>0 and s2.lo>0 and small_gap.lo>0
    return {
      'sigma_interval':qij(sigma),'discriminant_interval':qij(disc),
      'u_minus_interval':qij(uminus),'u_plus_interval':qij(uplus),'u_coefficient_switch_interval':qij(uc),
      'uniform_signs':{
        '8sigma_minus_lambda5_uplus':qij(s1),
        '3lambda4_1minusuplus_minus_8sigma':qij(s2),
        'sigma_minus_tracehalf_upper':qij(small_gap),
      },
      'dangerous_u_formula':'u_+- = (1 +- sqrt(1-32*sigma/(3*lambda4+lambda5)))/2',
      'dmin2_formula':'((8*sigma-lambda5*u)*(8*sigma-3*lambda4*u*(1-u)))/(lambda5*(3*lambda4*(1-u)-8*sigma))',
      'dangerous_set':'u in [u_-,u_+] and d^2 >= dmin2(u), intersected with 0<=d<=u and 4(1-u)^2>=u^2-d^2',
      'logic':'lambda_min(C_minus)<sigma globally; on [u_-,u_+] the determinant is strictly decreasing in d^2, so lambda_max>=sigma iff det(sigma I-C_minus)<=0 iff d^2>=dmin2(u).'
    }


def build_results():
    sy=odd_domain_symbolics();u,d,l4,l5=sy['symbols']
    return {
      'R7P-057':{
        'partition':partition_formulas(),'parity_weight':parity_weight_proof(),
        'plus_conditional':'p_plus proportional to (X*Y*Z^4,2*X*Z,Y,2*Z^3); conditional C_plus is independent of J6 and equals covariance of the four even feature states.',
        'minus_conditional':'p_minus proportional to (Y,S,S^-1) on the three odd feature states; conditional C_minus is independent of J3 and J6.',
        'shared_field_warning':'Y is identical in both sectors and S=Z^sqrt(3); parity sectors may not be optimized independently.'
      },
      'R7P-058':{
        'probabilities':[str(x) for x in sy['probabilities']],
        'domain':sy['domain'],'reconstruction':sy['reconstruction'],
        'covariance':[[str(sp.factor(sy['covariance'][i,j])) for j in range(2)] for i in range(2)],
        'supremum':odd_supremum_proof(),
      },
      'R7P-059':dangerous_interval_data(),
    }

if __name__=='__main__':
    out=build_results();p=ROOT/'results/R7P-057_059_intraparity.json';p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
