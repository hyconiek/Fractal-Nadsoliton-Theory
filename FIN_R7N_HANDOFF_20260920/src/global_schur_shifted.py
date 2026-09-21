"""Rigorous off-centre rectangular characteristic checker.

Unlike frontier_local_boxes.py, this module does not expand only at the double
root.  It certifies arbitrary small rectangles in local coordinates

    x = r-r_*,  u = 1-s,  v = 1-t,  e = 1-q_even,

using second-order interval AD and a midpoint Taylor enclosure for the shifted
characteristic polynomials P, P1, c2 at the fixed threshold sigma_*.

For a PSD 3x3 Schur-reduced matrix, R7P-044 gives lambda2<=sigma_* whenever
c2>0 and (P<=0 or P1>=0).  Thus a whole rectangle is certified if its Taylor
enclosure has c2.lo>0 and either P.hi<=0 or P1.lo>=0.

The irrational odd-law relation z=t^sqrt(3) is rigorously enclosed.  Its value
uses rational exponent brackets 19/11 < sqrt(3) < 26/15; its first/second
v-derivatives use simple rational interval bounds valid for 0<t<=1.
"""
from __future__ import annotations
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import math

from intervals import QI, sqrt_interval
import boundary_ising as bi
import off_face_local as loc

ROOT=Path(__file__).resolve().parents[1]
N=4
Jet4=loc.Jet4
AL=(19,11); AH=(26,15)

@lru_cache(maxsize=200000)
def pow_frac_point(x:F,p:int,q:int,digits:int=18):
    """Rigorous decimal-rational enclosure of x^(p/q), x>=0."""
    x=F(x)
    if x<0: raise ValueError('negative base')
    if x==0:
        if p>0:return QI(0)
        raise ValueError('zero negative power')
    if p<0:
        I=pow_frac_point(x,-p,q,digits); return QI(1/I.hi,1/I.lo)
    Q=10**digits; Nn,D=x.numerator,x.denominator
    targetN=pow(Nn,p)*pow(Q,q); targetD=pow(D,p)
    approx=float(x)**(p/q); n=max(0,int(math.floor(approx*Q)))
    def le(k): return pow(k,q)*targetD <= targetN
    while n>0 and not le(n): n-=1
    while le(n+1): n+=1
    if pow(n,q)*targetD==targetN:return QI(F(n,Q))
    return QI(F(n,Q),F(n+1,Q))


def _constants(a=None):
    L=bi.strict_intervals(); l3,l4,l5,l6=L[3],L[4],L[5],L[6]
    sigma=(2*l3*(l4+l5)-l4*l5)/(24*l3)
    threshold=sigma if a is None else QI(a)
    c=l6/(3*threshold)
    tau=sqrt_interval(((2*l3-l4)*(2*l3-l5))/(4*l3*l3),40)
    rstar=(QI(1)-tau)/(QI(1)+tau)
    a3=sqrt_interval(l3/6,40); a4=sqrt_interval(l4/6,40); a5=sqrt_interval(l5/6,40)
    rt3=sqrt_interval(QI(3),40)
    Vp=[[a3,a4,a5],[a3,-a4/2,-a5/2],[-a3,a4,-a5],[-a3,-a4/2,a5/2]]
    Vm=[[QI(0),a4,QI(0)],[QI(0),-a4/2,rt3*a5/2],[QI(0),-a4/2,-rt3*a5/2]]
    return threshold,c,rstar,rt3,Vp,Vm


def _z_of_v(v:Jet4,vlo:F,vhi:F,rt3:QI):
    """Jet enclosure for z=(1-v)^sqrt(3) on vlo<=v<=vhi."""
    tlo=F(1)-F(vhi); thi=F(1)-F(vlo)
    if not (0<tlo<=thi<=1): raise ValueError('v box outside 0<=v<1')
    # For 0<t<=1, t^alpha decreases with alpha.  AL<sqrt3<AH.
    zv=QI(pow_frac_point(tlo,*AH).lo, pow_frac_point(thi,*AL).hi)
    # d z/dv = -alpha t^(alpha-1).  Since 0<alpha-1<1,
    # t <= t^(alpha-1) <= 1.
    d1=QI(rt3.lo*tlo,rt3.hi)
    g=[QI(0) for _ in range(N)]; g[2]=-d1
    # d2 z/dv2 = alpha(alpha-1)t^(alpha-2), with -1<alpha-2<0,
    # hence 1 <= t^(alpha-2) <= 1/tlo.
    af=rt3*(rt3-QI(1))
    H=[[QI(0) for _ in range(N)] for __ in range(N)]
    H[2][2]=af*QI(1,F(1)/tlo)
    return Jet4(zv,g,H)


def _qe_exact_jet(e,elo:F,ehi:F):
    # f(e)=e(1-e), exact range and derivatives on 0<=e<=1.
    def f(x): return x*(1-x)
    vals=[f(elo),f(ehi)]
    if elo<=F(1,2)<=ehi: vals.append(F(1,4))
    v=QI(min(vals),max(vals))
    g=[QI(0) for _ in range(N)]; g[3]=QI(1-2*ehi,1-2*elo)
    H=[[QI(0) for _ in range(N)] for __ in range(N)]; H[3][3]=QI(-2)
    return Jet4(v,g,H)

def model_jet(box,a=None):
    """Return interval jets (P,P1,c2) on an arbitrary global chart rectangle."""
    box=tuple((F(lo),F(hi)) for lo,hi in box)
    if len(box)!=4 or any(lo>hi for lo,hi in box): raise ValueError('bad box')
    (xlo,xhi),(ulo,uhi),(vlo,vhi),(elo,ehi)=box
    if ulo<0 or uhi>=1 or vlo<0 or vhi>=1 or elo<0 or ehi>=1: raise ValueError('outside global physical chart')
    threshold,c,rstar,rt3,Vp,Vm=_constants(a)
    x=Jet4.var(QI(xlo,xhi),0); u=Jet4.var(QI(ulo,uhi),1)
    v=Jet4.var(QI(vlo,vhi),2); e=Jet4.var(QI(elo,ehi),3)
    r=Jet4(rstar)+x; s=Jet4(1)-u; t=Jet4(1)-v
    z=_z_of_v(v,vlo,vhi,rt3)

    wp=[Jet4(1),2*s*(t**3),r*(t**4),2*r*s*t]
    wm=[z,s,s*(z**2)]
    mp,Cp=loc._weighted_stats(wp,Vp); mm,Cm=loc._weighted_stats(wm,Vm)
    d=loc._vsub(mp,mm)
    q=Jet4(1)-e; qe=_qe_exact_jet(e,elo,ehi); eta=Jet4(1)-qe*c
    M=loc._madd(loc._mscale(q,Cp),loc._mscale(e,Cm))
    M=loc._madd(M,loc._mscale(qe/eta,loc._outer(d)))
    K=[[(Jet4(threshold)-M[i][j] if i==j else -M[i][j]) for j in range(3)] for i in range(3)]
    a,b,cc=K[0]; dd=K[1][1]; ee=K[1][2]; ff=K[2][2]
    P=a*dd*ff+2*b*cc*ee-a*ee*ee-dd*cc*cc-ff*b*b
    P1=(a*dd-b*b)+(a*ff-cc*cc)+(dd*ff-ee*ee)
    c2=a+dd+ff
    return rstar,P,P1,c2


def _center_box(box):
    mids=[]; rads=[]
    for lo,hi in box:
        lo,hi=F(lo),F(hi); m=(lo+hi)/2
        mids.append(m); rads.append((hi-lo)/2)
    return tuple((m,m) for m in mids),rads


def _absmax(I): return max(abs(I.lo),abs(I.hi))


def _range_from_jets(full,cen,rads):
    lin=cen.v
    for i,r in enumerate(rads): lin=lin+cen.g[i]*QI(-r,r)
    rem=F(0)
    for i in range(N):
        for j in range(N): rem += F(1,2)*_absmax(full.H[i][j])*rads[i]*rads[j]
    return QI(lin.lo-rem,lin.hi+rem)


def ranges(box,a=None):
    box=tuple((F(lo),F(hi)) for lo,hi in box)
    full=model_jet(box,a); cb,rads=_center_box(box); cen=model_jet(cb,a)
    return {
      'rstar':full[0],
      'P':_range_from_jets(full[1],cen[1],rads),
      'P1':_range_from_jets(full[2],cen[2],rads),
      'c2':_range_from_jets(full[3],cen[3],rads),
    }


def raw_shifted_box(box,a=None):
    R=ranges(box,a)
    reason='P_NONPOS' if R['P'].hi<=0 else ('P1_NONNEG' if R['P1'].lo>=0 else None)
    ok=R['c2'].lo>0 and reason is not None
    return {'status':'INTERVAL_CERTIFIED' if ok else 'FAILED','reason':reason,**R}
