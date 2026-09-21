"""R7P-060--063 dominant-mass and intraparity closure.

Exact reductions preserve the shared fields between parity sectors.  Numerical
optimization is labelled numerical; proof decisions use rational interval
Bernstein certificates and strict spectral intervals.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import json, math
import sympy as sp
from scipy.optimize import minimize_scalar

from intervals import QI
from boundary_cover import eval_qi
import boundary_ising as bi

ROOT=Path(__file__).resolve().parents[1]


def qi_pair(I): return [str(I.lo),str(I.hi)]


def spectral_symbols():
    u,l3,l4,l5=sp.symbols('u l3 l4 l5', positive=True)
    sigma=sp.factor((2*l3*(l4+l5)-l4*l5)/(24*l3))
    N=sp.factor((8*sigma-l5*u)*(8*sigma-3*l4*u*(1-u)))
    D=sp.factor(l5*(3*l4*(1-u)-8*sigma))
    dmin2=sp.factor(N/D)
    return u,l3,l4,l5,sigma,N,D,dmin2


def bernstein_interval(expr,var,left,right,values):
    z=sp.Symbol('_z')
    P=sp.Poly(sp.expand(expr.subs(var,left+(right-left)*z)),z)
    n=P.degree()
    out=[]
    for i in range(n+1):
        c=sp.factor(sum(P.nth(j)*sp.binomial(i,j)/sp.binomial(n,j)
                        for j in range(i+1)))
        out.append(eval_qi(c,values))
    return out


def positive_cover(expr,var,left,right,values,max_depth=18):
    """Adaptive exact-rational 1D Bernstein proof of expr>0."""
    stack=[(F(left),F(right),0)]; leaves=[]
    while stack:
        a,b,depth=stack.pop()
        coeff=bernstein_interval(expr,var,a,b,values)
        margin=min(x.lo for x in coeff)
        if margin>0:
            leaves.append({'interval':[str(a),str(b)],'depth':depth,
                           'min_bernstein_lo':str(margin)})
            continue
        if depth>=max_depth:
            raise RuntimeError(f'unresolved Bernstein interval {a}..{b}: {margin}')
        m=(a+b)/2
        stack.append((m,b,depth+1)); stack.append((a,m,depth+1))
    leaves.sort(key=lambda x:F(x['interval'][0]))
    return leaves


def dmin_float(u,l3,l4,l5):
    sigma=(2*l3*(l4+l5)-l4*l5)/(24*l3)
    return ((8*sigma-l5*u)*(8*sigma-3*l4*u*(1-u)) /
            (l5*(3*l4*(1-u)-8*sigma)))


def even_weights_from_ud(u,d,X=1.0):
    """Same-field C+ weights reconstructed from the odd (u,d) variables."""
    R=(u+d)/(u-d)
    Z=R**(1/(2*math.sqrt(3)))
    Y=2*(1-u)/math.sqrt(u*u-d*d)
    w=[X*Y*Z**4,2*X*Z,Y,2*Z**3]
    s=sum(w); p=[x/s for x in w]
    return Y,Z,p


def dominant_reduction():
    """R7P-060: exact reduction plus separately labelled numerical minimum."""
    L=bi.strict_intervals()
    l3,l4,l5=[float((L[k].lo+L[k].hi)/2) for k in (3,4,5)]
    sigma=(2*l3*(l4+l5)-l4*l5)/(24*l3)
    disc=1-32*sigma/(3*l4+l5)
    ulo=(1-math.sqrt(disc))/2; uhi=(1+math.sqrt(disc))/2

    # A small proof-only bound Y^2>4 on the whole dangerous interval.  It is
    # intentionally weaker than R7P-061's sharp coarse certificate.
    u,sl3,sl4,sl5,ssigma,N,D,_=spectral_symbols()
    H4=sp.cancel((4*(1-u)**2*D-4*(u**2*D-N))*(24*sl3)**2)
    vals={sl3:L[3],sl4:L[4],sl5:L[5]}
    cover4=positive_cover(H4,u,F(2,5),F(3,5),vals)

    def objective(x):
        y=dmin_float(x,l3,l4,l5)
        if y<0 or y>=x*x: return 10.0
        d=math.sqrt(y)
        return even_weights_from_ud(x,d,1.0)[2][0]
    res=minimize_scalar(objective,bounds=(ulo,uhi),method='bounded',
                        options={'xatol':1e-14,'maxiter':1000})
    us=float(res.x); ds=math.sqrt(dmin_float(us,l3,l4,l5));Y,Z,p=even_weights_from_ud(us,ds,1.0)
    return {
      'status':'EXACT_REDUCTION_PLUS_NUMERICAL_MINIMUM',
      'same_field_formulas':{
        'R':'(u+d)/(u-d)',
        'Z':'R^(1/(2*sqrt(3)))',
        'Y':'2*(1-u)/sqrt(u^2-d^2)',
        'Cplus_weights':'(X*Y*Z^4, 2*X*Z, Y, 2*Z^3), X=exp(2J3)>=1'
      },
      'monotonicity':[
        'p1=X*a/(X*(a+b)+c) is strictly increasing in X for c>0, so the minimum has X=1 (J3=0).',
        'For fixed u, both Y and Z increase with d; p1=1/(1+2/(Y Z^3)+Z^-4+2/(Y Z)) therefore increases with d.',
        'Dangerous Cminus has d^2>=dmin^2(u), hence the minimum at each u has d=dmin(u).',
        'The interval proof Y^2>4 on the dangerous u-range makes p1 the largest Cplus mass: p1/p2=Y Z^3/2>1, p1/p3=X Z^4>=1, p1/p4=X Y Z/2>1.'
      ],
      'Y2_gt_4_cover':cover4,
      'one_dimensional_expression':'alpha(u)=1/(1+2/(Y(u) Z(u)^3)+Z(u)^(-4)+2/(Y(u) Z(u))), with d(u)=dmin(u), X=1',
      'numerical_candidate':{
        'spectral_choice':'midpoints of accepted strict intervals',
        'u':us,'d':ds,'Y':Y,'Z':Z,'p':p,'p_dom':p[0],
        'optimizer_success':bool(res.success)
      },
      'proof_scope':'The reduction/monotonicity is exact; the displayed minimizing u and p_dom are numerical, not a global interval certificate.'
    }


def coarse_dominant_certificate():
    """R7P-061 strict coarse bounds and a rational p_dom floor."""
    L=bi.strict_intervals(); u,l3,l4,l5,sigma,N,D,_=spectral_symbols()
    vals={l3:L[3],l4:L[4],l5:L[5]}
    k=sp.Rational(1703,2000)       # 0.8515
    y2=sp.Rational(229,20)         # 11.45
    # D>0 on [2/5,3/5].  The dangerous u interval is strictly inside it.
    Dpoly=sp.cancel(D*24*l3)
    Dcover=positive_cover(Dpoly,u,F(2,5),F(3,5),vals)
    Hratio=sp.cancel((N-k**2*u**2*D)*(24*l3)**2)
    HY=sp.cancel((4*(1-u)**2*D-y2*(u**2*D-N))*(24*l3)**2)
    ratio_cover=positive_cover(Hratio,u,F(2,5),F(3,5),vals)
    y_cover=positive_cover(HY,u,F(2,5),F(3,5),vals)

    # Confirm certified dangerous interval subset of the rational cover interval.
    dd=__import__('intraparity').dangerous_interval_data()
    ulo=QI(F(dd['u_minus_interval'][0]),F(dd['u_minus_interval'][1]))
    uhi=QI(F(dd['u_plus_interval'][0]),F(dd['u_plus_interval'][1]))
    assert ulo.lo>F(2,5) and uhi.hi<F(3,5)

    # Convert the two coarse field bounds into an entirely rational mass floor.
    # d/u>=k -> R=(u+d)/(u-d)>=(1+k)/(1-k)=3703/297.
    R0=(F(1)+F(k))/(F(1)-F(k))
    z0=F(2071,1000)
    y0=F(3383,1000)
    # 2 sqrt(3) < 693/200, because its square exceeds 12.  Therefore
    # 1/(2sqrt3)>200/693.  Exact integer/rational exponent comparison gives
    # R0^(200/693) >= z0 iff R0^200 >= z0^693.
    assert F(693,200)**2>12
    assert R0**200 > z0**693
    assert y0*y0 < F(y2)
    p0=1/(F(1)+F(2)/(y0*z0**3)+F(1)/(z0**4)+F(2)/(y0*z0))
    assert p0>F(711,1000)

    return {
      'status':'INTERVAL_CERTIFIED',
      'dangerous_u_containment':{'outer_rational_interval':['2/5','3/5'],
                                 'u_minus_interval':dd['u_minus_interval'],'u_plus_interval':dd['u_plus_interval']},
      'denominator_positive_cover':Dcover,
      'd_over_u_bound':{'bound':'d/u >= 1703/2000','squared_polynomial_cover':ratio_cover},
      'Y2_bound':{'bound':'Y^2 >= 229/20','polynomial_cover':y_cover},
      'field_corollaries':{
        'R_lower':str(R0),'Y_lower':str(y0),'Z_lower':str(z0),
        'Z_proof':'2*sqrt(3)<693/200 and (3703/297)^200>(2071/1000)^693'
      },
      'dominant_mass_exact_rational_floor':str(p0),
      'dominant_mass_decimal_floor':float(p0),
      'exported_simple_bound':'p_dom > 711/1000',
      'note':'The rational 711/1000 bound is freshly recomputed from the certified coarse d/u and Y^2 bounds; the imported decimal 0.7112098557 is not used as a proof input.'
    }


def covariance_envelope_certificate():
    """R7P-062 exact tetrahedral envelope reduction + interval upper bound."""
    L=bi.strict_intervals(); l3,l4,l5=sp.symbols('l3 l4 l5', positive=True)
    vals={l3:L[3],l4:L[4],l5:L[5]}
    alpha=sp.Rational(711,1000); x=sp.symbols('x', nonnegative=True)
    A=(l3*l4+l3*l5+l4*l5)/4
    B=(4*l3*l4+4*l3*l5+l4*l5)/16
    # If the dominant vertex is 1 or 3, fixed x=(p2+p4)/2 and AM-GM
    # gives p2=p4=x at the maximum.  Vertices 2 or 4 give the swapped branch.
    E13=sp.expand(2*A*alpha*x*(1-alpha-2*x)+B*x**2*(1-2*x))
    E24=sp.expand(A*x**2*(1-2*x)+2*B*alpha*x*(1-alpha-2*x))
    xmax=(1-alpha)/2
    a=sp.Rational(511,2000)
    target=a*a
    cov13=positive_cover(sp.expand(target-E13),x,F(0),F(xmax),vals)
    cov24=positive_cover(sp.expand(target-E24),x,F(0),F(xmax),vals)
    return {
      'status':'INTERVAL_CERTIFIED',
      'alpha':'711/1000',
      'e2_exact':'A*(p1*p2*p3+p1*p3*p4)+B*(p1*p2*p4+p2*p3*p4), A=(l3*l4+l3*l5+l4*l5)/4, B=(4l3l4+4l3l5+l4l5)/16',
      'branch_13':{'parameterization':'p1=alpha (or p3=alpha), p2=p4=x at the maximizing split, p3=1-alpha-2x',
                   'envelope':str(E13),'cover':cov13},
      'branch_24':{'parameterization':'p2=alpha (or p4=alpha), p1=p3=x at the maximizing split, p4=1-alpha-2x',
                   'envelope':str(E24),'cover':cov24},
      'monotonicity':'For alpha>1/2 the partial derivative of either branch with respect to alpha at fixed x is 2*C*x*(1-2alpha-2x)<0; hence p_i>=alpha is maximized at p_i=alpha.',
      'certified_e2_bound':'e2 < (511/2000)^2',
      'lambda2_bound':'lambda2(Cplus) < 511/2000',
      'justification_lambda2':'Cplus is PSD; with eigenvalues mu1>=mu2>=mu3>=0, e2=mu1*mu2+mu1*mu3+mu2*mu3 >= mu2^2.'
    }


def intraparity_weyl_certificate():
    """R7P-063 complete two-case W_par theorem on the nonnegative field domain."""
    L=bi.strict_intervals();l3,l4,l5=L[3],L[4],L[5]
    sigma=(2*l3*(l4+l5)-l4*l5)/(24*l3)
    cminus_sup=(3*l4+l5)/32
    cplus_bound=QI(F(511,2000))
    gap=2*sigma-cminus_sup-cplus_bound
    assert gap.lo>0
    return {
      'status':'INTERVAL_CERTIFIED',
      'domain':'shared nonnegative fields J3,J4,J5,J6>=0 in the four-amplitude parity decomposition',
      'sigma_interval':qi_pair(sigma),
      'Cminus_sup_interval':qi_pair(cminus_sup),
      'dangerous_Cplus_bound':'lambda2(Cplus)<511/2000',
      'q_bound':'q_even>=1/2',
      'strict_two_sigma_gap_interval':qi_pair(gap),
      'case_A':'If lambda1(Cminus)<=sigma, R7P-055 gives lambda2(Cplus)<=sigma; Weyl gives lambda2(Wpar)<=q*sigma+(1-q)*sigma=sigma.',
      'case_B':'If lambda1(Cminus)>=sigma, R7P-059 makes the odd sector dangerous; R7P-061 gives p_dom>711/1000 and R7P-062 gives lambda2(Cplus)<511/2000. R7P-058 gives lambda1(Cminus)<=(3lambda4+lambda5)/32. Since q>=1/2 and the Cminus bound exceeds 511/2000, the affine Weyl RHS is maximized at q=1/2, where the certified two-sigma gap is positive.',
      'conclusion':'lambda2(W_par)<=sigma_* throughout the stated shared-field nonnegative domain.',
      'limitations':['This is an intraparity W_par theorem, not yet the full M=W_par+b b^T off-face theorem.','No claim is made for negative J6 or independently optimized parity sectors.']
    }


def build_results():
    return {'R7P-060':dominant_reduction(),
            'R7P-061':coarse_dominant_certificate(),
            'R7P-062':covariance_envelope_certificate(),
            'R7P-063':intraparity_weyl_certificate()}


if __name__=='__main__':
    out=build_results(); p=ROOT/'results/R7P-060_063_intraparity_closure.json'
    p.write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
