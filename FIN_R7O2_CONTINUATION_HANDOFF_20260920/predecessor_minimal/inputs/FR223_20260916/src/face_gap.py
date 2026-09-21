"""R7P-022 quantitative extreme-face gap and exact local remainder identities."""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import json
from intervals import QI, sqrt_interval

ROOT=Path(__file__).resolve().parents[1]

def load_L():
    d=json.loads((ROOT/'inputs/fin_handoff_audit/results.json').read_text())
    return [QI(F(lo),F(hi)) for lo,hi in d['exact']['laplacian_intervals']]

def constants():
    L=load_L(); l3,l4,l5,l6=L[3],L[4],L[5],L[6]
    sigma=(2*l3*(l4+l5)-l4*l5)/(24*l3); a=l3/6; c=l6/(3*sigma)
    t2=QI(1)-sigma/a
    return L,sigma,a,c,t2

def qone(qlo,qhi):
    # q(1-q) decreases on [1/2,1].
    return QI(qhi*(1-qhi),qlo*(1-qlo))

def scalar_f_interval(qlo,qhi,xlo,xhi):
    _,_,a,c,_=constants(); q=QI(qlo,qhi); omq=QI(1-qhi,1-qlo)
    eta=QI(1)-c*qone(qlo,qhi); x2=QI(xlo*xlo,xhi*xhi)
    return a*(q+x2*(c*omq-QI(1))/eta)

def block_plus_interval(xlo,xhi):
    L,_,_,_,_=constants();l4,l5=L[4],L[5];x2=QI(xlo*xlo,xhi*xhi)
    disc=(l5-l4)**2+4*l4*l5*x2
    return (l4+l5+sqrt_interval(disc,40))/24

def block_minus_global_upper():
    L,_,_,_,_=constants();return L[4]/12  # lambda4<lambda5 on accepted intervals

def local_coefficients():
    L,sigma,a,c,t2=constants();l4,l5=L[4],L[5]
    # For e<=1/100 and y=x^2-t_*^2>=0,
    # G=(sigma-f)eta = a*y*(1-c e)+e*(a-c*e*(2a-sigma)+a*c*e^2).
    emax=F(1,100)
    cy=a*(QI(1)-c*emax)
    ce=a-c*emax*(2*a-sigma)  # drop the positive +a*c*e^2 term
    # For x<=t*, rationalization gives
    # sigma-B(x)=l4*l5*(t2-x2)/(6*(sqrt(D*)+sqrt(Dx)))
    # and sqrt(D)<=l4+l5 for x in [0,1].
    kb=l4*l5/(12*(l4+l5))
    assert cy.lo>F(3,10) and ce.lo>F(3,10) and kb.lo>F(9,100)
    assert (sigma-block_minus_global_upper()).lo>F(1,100)
    return {
      'tstar2_interval':[str(t2.lo),str(t2.hi)],
      'scalar_y_coefficient_interval':[str(cy.lo),str(cy.hi)],
      'scalar_e_coefficient_conservative_interval':[str(ce.lo),str(ce.hi)],
      'block_rationalized_coefficient_interval':[str(kb.lo),str(kb.hi)],
      'certified_simple_bounds':{'scalar':'sigma-f >= 3/10*((1-q)+(x^2-tstar^2)) when q>=99/100 and x^2>=tstar^2',
                                 'block':'sigma-Bplus >= 9/100*(tstar^2-x^2) when x^2<=tstar^2'},
      'exact_scalar_numerator_identity':'(sigma-f)*eta = a*c*e^3-2*a*c*e^2-a*c*e*y+a*e+a*y+c*e^2*sigma, e=1-q, y=x^2-tstar^2',
    }

GAP=F(1,10000); NQ=F(99,100); NXLO=F(21,50); NXHI=F(11,25)

def _physical_none(box):
    q0,q1,x0,x1=box; return x0*x0>2*q1-1

def _in_neighborhood(box):
    q0,q1,x0,x1=box; return q0>=NQ and x0>=NXLO and x1<=NXHI

def _classify(box):
    _,sigma,_,_,_=constants()
    if _physical_none(box): return 'OUTSIDE_PHYSICAL'
    if _in_neighborhood(box): return 'EQUALITY_NEIGHBORHOOD'
    q0,q1,x0,x1=box; target=sigma.lo-GAP
    if block_plus_interval(x0,x1).hi<=target:return 'BLOCK_GAP'
    if scalar_f_interval(q0,q1,x0,x1).hi<=target:return 'SCALAR_GAP'
    return None

def complement_cover(maxdepth=40):
    stack=[((F(1,2),F(1),F(0),F(1)),0)];leaves=[];un=[]
    while stack:
        box,d=stack.pop(); reason=_classify(box)
        if reason:
            leaves.append({'box':[str(z) for z in box],'reason':reason});continue
        if d>=maxdepth:
            un.append([str(z) for z in box]);continue
        q0,q1,x0,x1=box
        cuts=[]
        if q0<NQ<q1: cuts.append(('q',NQ))
        if x0<NXLO<x1: cuts.append(('x',NXLO))
        if x0<NXHI<x1: cuts.append(('x',NXHI))
        if cuts: dim,mid=cuts[0]
        elif q1-q0>=x1-x0: dim,mid='q',(q0+q1)/2
        else: dim,mid='x',(x0+x1)/2
        if dim=='q':
            stack.extend([((mid,q1,x0,x1),d+1),((q0,mid,x0,x1),d+1)])
        else:
            stack.extend([((q0,q1,mid,x1),d+1),((q0,q1,x0,mid),d+1)])
    counts={}
    for z in leaves:counts[z['reason']]=counts.get(z['reason'],0)+1
    return {'gap':str(GAP),'excluded_neighborhood':{'q':['99/100','1'],'x':['21/50','11/25']},
            'leaf_count':len(leaves),'reason_counts':counts,'unresolved_count':len(un),'unresolved':un,'leaves':leaves}

def build():
    loc=local_coefficients();cov=complement_cover();assert cov['unresolved_count']==0
    out={'id':'R7P-022-face-gap','status':'INTERVAL_CERTIFIED','domain':'extreme face J4=J5=0 in physical (q,x): 1/2<=q<=1, x>=0, x^2<=2q-1',
         'away_from_equality':cov,'local_exact_remainder':loc,
         'conclusion':'Outside q>=99/100, 21/50<=x<=11/25, lambda2<=sigma-1/10000. Inside that neighborhood the exact one-sided identities give a zero only at q=1,x=tstar and quantitative first-order control on both sides of x=tstar.',
         'scope':'extreme face only; transverse J4,J5 control remains R7P-068.'}
    p=ROOT/'certificates/R7P-022_face_gap.json';p.write_text(json.dumps(out,indent=2)+'\n')
    return out
if __name__=='__main__':
    o=build();print(json.dumps({k:o[k] for k in ['status','conclusion']},indent=2));print(o['away_from_equality']['reason_counts'])
