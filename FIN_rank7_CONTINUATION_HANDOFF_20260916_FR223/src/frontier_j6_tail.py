"""FR5: durable global large-J6 tail certificate for the physical four-amplitude model.

The proof combines:
1. a shifted boundary-Ising Bernstein cover at theta=sigma_*-10^-6 outside a
   tiny dyadic equality box contained in the already-certified FR13 box;
2. the exact physical parity relation, which gives odd mass e<=y/(1+y) for
   y=exp(-2J6), because the zero-J6 even conditional mass is at least 1/2;
3. a covariance contamination bound M4 <= C_+ + e D^2 I with D^2<5.

At y<=1/5,000,000 the perturbation is <10^-6, so the shifted boundary reserve
closes the exterior.  The quarantined equality box is covered directly by FR13.
"""
from __future__ import annotations
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import json
import sympy as sp

import boundary_cover as bc
import boundary_ising as bi
from intervals import QI, sqrt_interval

ROOT=Path(__file__).resolve().parents[1]
DELTA=F(1,10**6)
YMAX=F(1,5_000_000)
# A dyadic box around r_* contained in FR13's projection.
LOCAL_BOX=((F(3293,8192),F(3294,8192)),(F(8191,8192),F(1)),(F(8191,8192),F(1)))


def shifted_root(delta: F=DELTA):
    inv=bi.symbolic_invariants(); l3,l4,l5=inv['l']; pvars=inv['p']
    r,s,t=sp.symbols('r s t', nonnegative=True)
    d=1+2*s*t**3+r*t**4+2*r*s*t
    probs=[1/d,2*s*t**3/d,r*t**4/d,2*r*s*t/d]
    sub=dict(zip(pvars,probs))
    e1=sp.factor(inv['e1'].subs(sub)); e2=sp.factor(inv['e2'].subs(sub)); e3=sp.factor(inv['e3'].subs(sub))
    sigma=sp.factor((2*l3*(l4+l5)-l4*l5)/(24*l3))
    theta=sigma-sp.Rational(delta.numerator,delta.denominator)
    P=sp.factor(theta**3-e1*theta**2+e2*theta-e3)
    Pp=sp.factor(3*theta**2-2*e1*theta+e2)
    A=sp.factor(sp.cancel(P*d**4*(24*l3)**3))
    B=sp.factor(sp.cancel(Pp*d**3*(24*l3)**2))
    L=bi.strict_intervals(); sv={l3:L[3],l4:L[4],l5:L[5]}
    degA,powA=bc._power_coeff_intervals(sp.expand(A),(r,s,t),sv)
    degB,powB=bc._power_coeff_intervals(sp.expand(B),(r,s,t),sv)
    return {
      'degrees_A':degA,'degrees_B':degB,
      'A':bc.power_to_bernstein(powA,degA),'B':bc.power_to_bernstein(powB,degB),
      'power_terms_A':len(powA),'power_terms_B':len(powB),
    }


def shifted_boundary_cover(max_leaves=200000,max_depth=50):
    root=shifted_root(); unit=((F(0),F(1)),)*3
    stack=[(unit,root['A'],root['B'],0,'')]; leaves=[]; splits=0
    while stack:
        box,A,B,depth,path=stack.pop()
        a=bc.bernstein_bounds(A); b=bc.bernstein_bounds(B)
        if a.hi<=0:
            reason='SAFE_A_NONPOS'; safe=True; margin=str(-a.hi)
        elif b.lo>=0:
            reason='SAFE_B_NONNEG'; safe=True; margin=str(b.lo)
        elif bc._inside(box,LOCAL_BOX):
            reason='LOCAL_FR13'; safe=True; margin='FR13 projection + y-tail parity bound'
        elif any(lo==hi==0 for lo,hi in box):
            reason='BOUNDARY_RANK_LEMMA'; safe=True; margin='rank(Cov)<=1'
        else:
            reason='UNRESOLVED'; safe=False; margin=None
        if safe:
            leaves.append({'path':path,'depth':depth,'box':bc._box_json(box),'reason':reason,
                           'A_bounds':bc.qi_json(a),'B_bounds':bc.qi_json(b),'margin':margin})
            continue
        if len(leaves)+len(stack)+1>=max_leaves:
            leaves.append({'path':path,'depth':depth,'box':bc._box_json(box),'reason':'UNRESOLVED_BUDGET',
                           'A_bounds':bc.qi_json(a),'B_bounds':bc.qi_json(b)})
            continue
        if depth>=max_depth:
            leaves.append({'path':path,'depth':depth,'box':bc._box_json(box),'reason':'UNRESOLVED_DEPTH',
                           'A_bounds':bc.qi_json(a),'B_bounds':bc.qi_json(b)})
            continue
        axis=bc._choose_axis(box,root['degrees_A'],root['degrees_B'])
        boxL,boxR=bc.split_box(box,axis)
        AL,AR=bc.split_bernstein(A,root['degrees_A'],axis)
        BL,BR=bc.split_bernstein(B,root['degrees_B'],axis)
        stack.append((boxR,AR,BR,depth+1,path+f'{axis}R'))
        stack.append((boxL,AL,BL,depth+1,path+f'{axis}L'))
        splits+=1
    counts=Counter(x['reason'] for x in leaves)
    unresolved=[x for x in leaves if x['reason'].startswith('UNRESOLVED')]
    return {
      'threshold':'sigma_*-1/1000000',
      'local_box':bc._box_json(LOCAL_BOX),
      'splits':splits,'terminal_leaves':len(leaves),'max_depth_reached':max(x['depth'] for x in leaves),
      'reason_counts':dict(counts),'unresolved_count':len(unresolved),
      'global_exterior_pass':len(unresolved)==0,
      'polynomials':{'A_degrees':list(root['degrees_A']),'B_degrees':list(root['degrees_B']),
                     'A_power_terms':root['power_terms_A'],'B_power_terms':root['power_terms_B']},
      'leaves':leaves,
    }


def diameter_bound():
    """Rigorous finite enumeration of the C4 feature diameter; only D^2<5 is used."""
    L=bi.strict_intervals(); rt3=sqrt_interval(QI(3),50)
    # By cyclicity, pair distances depend only on d=1,...,6.  Coefficients of
    # lambda3,lambda4,lambda5,lambda6 in ||C4_j-C4_0||^2.
    coeffs={
      1:(QI(F(1,6)),QI(F(3,8)),QI(F(7,24))+rt3*F(1,6),QI(F(1,3))),
      2:(QI(F(2,3)),QI(F(3,8)),QI(F(1,24)),QI(0)),
      3:(QI(F(1,6)),QI(0),QI(F(1,6)),QI(F(1,3))),
      4:(QI(0),QI(F(3,8)),QI(F(3,8)),QI(0)),
      5:(QI(F(1,6)),QI(F(3,8)),QI(F(7,24))-rt3*F(1,6),QI(F(1,3))),
      6:(QI(F(2,3)),QI(0),QI(F(2,3)),QI(0)),
    }
    rows=[]; mx=QI(0)
    for d,c in coeffs.items():
        z=c[0]*L[3]+c[1]*L[4]+c[2]*L[5]+c[3]*L[6]
        rows.append({'cyclic_distance':d,'distance2_interval':bc.qi_json(z)})
        if z.hi>mx.hi: mx=z
    if not mx.hi<F(5): raise AssertionError('feature diameter bound D^2<5 failed')
    return {'rows':rows,'diameter2_upper':str(mx.hi),'coarse_bound':'D^2<5'}


def local_box_checks():
    signs=bi.interval_signs(); T=QI(signs['t'][0],signs['t'][1]); R=(QI(1)-T)/(QI(1)+T)
    rx=F(1,6200); ru=rv=F(1,6800)
    rlo,rhi=LOCAL_BOX[0]; slo,shi=LOCAL_BOX[1]; tlo,thi=LOCAL_BOX[2]
    inside=(rlo>=R.hi-rx and rhi<=R.lo+rx and 1-slo<=ru and 1-tlo<=rv)
    emax=YMAX/(1+YMAX)
    return {
      'rstar_interval':bc.qi_json(R),'FR13_radii':{'x':'1/6200','u':'1/6800','v':'1/6800','e':'1/400'},
      'local_box_inside_FR13_projection':inside,
      'odd_mass_bound_at_ymax':str(emax),
      'odd_mass_inside_FR13_e_radius':emax<=F(1,400),
    }


def build_record(include_leaves=True):
    cov=shifted_boundary_cover(); dia=diameter_bound(); loc=local_box_checks()
    signs=bi.interval_signs(); trace_margin_lo=F(signs['trace_margin'][0])
    shifted_z2_margin=trace_margin_lo-3*DELTA
    emax=YMAX/(1+YMAX); perturb=F(5)*emax
    ok=(cov['global_exterior_pass'] and loc['local_box_inside_FR13_projection'] and
        loc['odd_mass_inside_FR13_e_radius'] and shifted_z2_margin>0 and perturb<DELTA)
    if not include_leaves: cov.pop('leaves',None)
    return {
      'id':'FR5-global-large-J6-tail','status':'INTERVAL_CERTIFIED_REPLAYED' if ok else 'FAILED',
      'domain':'J3,J4,J5,J6>=0 with y=exp(-2J6)<=1/5,000,000 and the exact physical shared-field parity law',
      'conclusion':'lambda2(M4)<=sigma_* throughout the declared global large-J6 tail.',
      'boundary_shift':{
        'delta':'1/1000000','shifted_z2_margin_lower':str(shifted_z2_margin),
        'criterion':'On the boundary C_+, lambda2<=sigma_*-delta outside LOCAL_BOX by the shifted P/Pprime Bernstein disjunction.',
        'cover':cov,
      },
      'parity_mass':{
        'input':'R7P-057 gives q0=Z_even/(Z_even+Z_odd)>=1/2 at J6=0, hence Z_odd<=Z_even.',
        'identity':'e=1-q_even=y Z_odd/(Z_even+y Z_odd)<=y/(1+y)',
        'ymax':'1/5000000','emax':str(emax),
      },
      'covariance_perturbation':{
        'diameter':dia,
        'loewner_bound':'M4-C_+ <= e D^2 I <= 5e I',
        'derivation':'For every unit z, the positive change is <=e[Var_-(zX)+(m_--m_+)^2]=e E_-[(zX-m_+)^2]<=e D^2.',
        'perturbation_upper':str(perturb),'strictly_below_delta':perturb<DELTA,
      },
      'local_quarantine':loc,
      'dependencies':['R7P-055 boundary theorem machinery','R7P-057 physical q0>=1/2','FR13 local anisotropic box'],
      'scope':'Physical four-amplitude positive orthant only; no full-seven-coordinate transfer.',
    }


def write_result(path=None):
    rec=build_record(True)
    if path is None: path=ROOT/'results/FR5_global_large_J6_tail.json'
    Path(path).write_text(json.dumps(rec,indent=2)+'\n')
    return rec

if __name__=='__main__':
    r=write_result(); print(json.dumps({'status':r['status'],'leaves':r['boundary_shift']['cover']['terminal_leaves'],
      'counts':r['boundary_shift']['cover']['reason_counts'],'perturbation_upper':r['covariance_perturbation']['perturbation_upper']},indent=2))
