"""R7P-017--024 independent face-certificate reconstruction from saved strict spectral intervals."""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import json, math, sys
sys.set_int_max_str_digits(100000)
import sympy as sp
from .intervals import QI, sqrt_interval, log_interval, div

ROOT=Path(__file__).resolve().parents[1]

def load_L():
    d=json.loads((ROOT/'inputs/fin_handoff_audit/results.json').read_text())
    return [QI(F(lo),F(hi)) for lo,hi in d['exact']['laplacian_intervals']]

def ev(expr, vals):
    expr=sp.sympify(expr)
    if expr.is_Rational:return QI(F(int(expr.p),int(expr.q)))
    if expr.is_Symbol:return vals[expr]
    if expr.is_Add:
        out=QI(0)
        for z in expr.args:out=out+ev(z,vals)
        return out
    if expr.is_Mul:
        out=QI(1)
        for z in expr.args:out=out*ev(z,vals)
        return out
    if expr.is_Pow and expr.exp.is_Integer:return ev(expr.base,vals)**int(expr.exp)
    raise ValueError(expr)

def bernstein_expr(poly,var,left,right):
    z=sp.symbols('z')
    P=sp.Poly(sp.expand(poly.subs(var,left+(right-left)*z)),z)
    n=P.degree()
    return [sp.factor(sum(P.nth(k)*sp.binomial(i,k)/sp.binomial(n,k) for k in range(i+1))) for i in range(n+1)]

def cert_data():
    L=load_L(); l3,l4,l5,l6=L[3],L[4],L[5],L[6]
    sigma=(2*l3*(l4+l5)-l4*l5)/(24*l3); aI=l3/6; cI=l6/(3*sigma)
    a,s,c,r,q=sp.symbols('a s c r q')
    vals={a:aI,s:sigma,c:cI}
    d0=s*(1+r)-a*r**2
    D=sp.expand((1+r)**2*d0)
    N=sp.expand(D-a*r*(1-r)*(1+r)**2-c*r*d0)
    face=[]
    for lo,hi in [(sp.Rational(0),sp.Rational(1,2)),(sp.Rational(1,2),sp.Rational(1))]:
        face.append([ev(x,vals) for x in bernstein_expr(N,r,lo,hi)])
    eta=1-c*q*(1-q)
    R=lambda x2: sp.expand(s*eta-a*q*eta-a*x2*(c*(1-q)-1))
    qmin=1-s/(2*a);qcrit=1-1/c
    low=[ev(x,vals) for x in bernstein_expr(R(2*q-1),q,qmin,qcrit)]
    upper_poly=sp.factor(R(1-s/a))
    assert sp.factor(upper_poly.subs(q,1))==0
    quotient=sp.factor(sp.cancel(upper_poly/(1-q)))
    upper=[ev(x,vals) for x in bernstein_expr(quotient,q,qcrit,1)]
    # corrected resolvent derivative, cancel only r+1>0
    Q=sp.factor(sp.diff(N,r)*D-N*sp.diff(D,r))
    assert sp.rem(sp.Poly(Q,r),sp.Poly(r+1,r))==0
    Q5=sp.factor(Q/(r+1))
    return dict(L=L,sigma=sigma,a=aI,c=cI,vals=vals,r=r,q=q,N=N,D=D,Q5=Q5,
                face=face,qmin=qmin,qcrit=qcrit,low=low,upper=upper,upper_quotient=quotient)

def interval_strings(rows):
    return [[[str(x.lo),str(x.hi)] for x in row] for row in rows]

def poly_bern_interval(poly,var,lo,hi,vals):
    return [ev(x,vals) for x in bernstein_expr(poly,var,sp.Rational(lo.numerator,lo.denominator),sp.Rational(hi.numerator,hi.denominator))]

def root_certificate():
    D=cert_data(); r=D['r'];Q5=D['Q5'];vals=D['vals']
    lo=F(6608385,10_000_000);hi=F(6608387,10_000_000)
    left=poly_bern_interval(Q5,r,F(0),lo,vals)
    right=poly_bern_interval(Q5,r,hi,F(1),vals)
    qlo=ev(Q5.subs(r,sp.Rational(lo.numerator,lo.denominator)),vals)
    qhi=ev(Q5.subs(r,sp.Rational(hi.numerator,hi.denominator)),vals)
    dQ=sp.diff(Q5,r)
    deriv=poly_bern_interval(dQ,r,lo,hi,vals)
    ok=(all(x.hi<0 for x in left) and all(x.lo>0 for x in right)
        and qlo.hi<0 and qhi.lo>0 and all(x.lo>0 for x in deriv))
    return {'ok':ok,'root_box':[str(lo),str(hi)],
            'left_bernstein':[[str(x.lo),str(x.hi)] for x in left],
            'right_bernstein':[[str(x.lo),str(x.hi)] for x in right],
            'endpoint_Q5':[[str(qlo.lo),str(qlo.hi)],[str(qhi.lo),str(qhi.hi)]],
            'derivative_bernstein':[[str(x.lo),str(x.hi)] for x in deriv]}

def location_certificate():
    D=cert_data(); rc=root_certificate(); lo,hi=map(F,rc['root_box']); RI=QI(lo,hi)
    aI=D['a'];sigma=D['sigma'];cI=D['c']
    # q=1/(1+r), t=sqrt(1-r^2)
    qI=QI(1)/(QI(1)+RI)
    tI=sqrt_interval(QI(1)-RI**2)
    # s3 = acosh(1/r)/sqrt(a) = log(1/r + sqrt(1/r^2-1))/sqrt(a)
    invr=QI(1)/RI
    acosh_arg=invr+sqrt_interval(invr**2-QI(1))
    acosh=log_interval(acosh_arg)
    s3I=acosh/sqrt_interval(aI)
    # S=N/D interval on narrow root box using direct expression via polynomial evaluation (dependency conservative)
    r= D['r']; vals=dict(D['vals']); vals[r]=RI
    NI=ev(D['N'],vals); DI=ev(D['D'],vals); SI=NI/DI
    return {'r':rc['root_box'],'q':[str(qI.lo),str(qI.hi)],'t':[str(tI.lo),str(tI.hi)],
            's3':[str(s3I.lo),str(s3I.hi)],'S':[str(SI.lo),str(SI.hi)]}

def build():
    D=cert_data(); root=root_certificate();loc=location_certificate()
    result={
      'id':'R7P-017-021-face',
      'claim_id':'CLM-FACE-STRICT',
      'domain':'corrected 1D resolvent r in [0,1] and extreme face s4=s5=0, s3,s6>=0',
      'quantifiers':'for every strict spectral tuple in the accepted outward intervals and every point in the stated face domains',
      'assumptions':['accepted strict spectral intervals','rank-seven C4 normalization','r=sech(sqrt(lambda3/6)*s3)'],
      'proof_type':'exact rational interval + Bernstein positivity + interval root isolation',
      'inputs':['inputs/fin_handoff_audit/results.json'],
      'conclusion':'corrected 1D resolvent is positive with a unique isolated minimum; extreme-face second-curvature bound holds with equality only at the compactified stated boundary',
      'global_pass':False,
      'face_resolvent':{'bernstein':interval_strings(D['face']),
          'all_positive':all(x.lo>0 for row in D['face'] for x in row),
          'denominator_endpoint_lower':[str(D['sigma'].lo),str((2*D['sigma']-D['a']).lo)]},
      'extreme_face':{'low_bernstein':interval_strings([D['low']])[0],
          'upper_quotient_bernstein':interval_strings([D['upper']])[0],
          'all_positive':all(x.lo>0 for x in D['low']+D['upper']),
          'exact_endpoint_factor':'1-q'},
      'resolvent_derivative':{'reduced_degree':int(sp.degree(D['Q5'],D['r'])),'factor_removed':'r+1 (>0 on [0,1])','symbolic':str(D['Q5'])},
      'unique_minimum':root,'location':loc,
      'equality_locus':{'finite_equality':False,'compactified':'q=1 and x=t_star on s4=s5=0; equivalent s6->+infinity at threshold balance'},
    }
    (ROOT/'certificates/R7P-017_021_face_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':
    out=build();print(json.dumps({'face_positive':out['face_resolvent']['all_positive'],'extreme_face_positive':out['extreme_face']['all_positive'],'root':out['unique_minimum'],'location':out['location']},indent=2))
