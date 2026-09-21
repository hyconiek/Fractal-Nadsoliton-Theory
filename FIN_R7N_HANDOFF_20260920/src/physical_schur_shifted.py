from pathlib import Path
from fractions import Fraction as F
import sys
ROOT=Path('/mnt/data/fin_rank7_next_campaign');H=ROOT/'inputs/FR223_20260916';sys.path[:0]=[str(ROOT/'src'),str(H/'src')]
from intervals import QI,sqrt_interval
import off_face_local as loc
import global_schur_shifted as gs
N=4; Jet4=loc.Jet4

def model_jet(compact_cell,a=F(67,250)):
    # Input exact compact (r,s,t,y). Internal chart variables are A=sqrt(r), u=1-s, v=1-t, y.
    (rlo,rhi),(slo,shi),(tlo,thi),(ylo,yhi)=tuple((F(x),F(z)) for x,z in compact_cell)
    if not (0<rlo<=rhi<=1 and 0<slo<=shi<=1 and 0<tlo<=thi<=1 and 0<=ylo<=yhi<=1): raise ValueError('bad compact cell')
    ar=sqrt_interval(QI(rlo,rhi),60); alo,ahi=ar.lo,ar.hi
    ulo,uhi=1-shi,1-slo; vlo,vhi=1-thi,1-tlo
    threshold,c,_,rt3,Vp,Vm=gs._constants(a)
    A=Jet4.var(QI(alo,ahi),0); u=Jet4.var(QI(ulo,uhi),1); v=Jet4.var(QI(vlo,vhi),2); y=Jet4.var(QI(ylo,yhi),3)
    s=Jet4(1)-u; t=Jet4(1)-v; z=gs._z_of_v(v,vlo,vhi,rt3)
    r=A*A
    wp=[Jet4(1),2*s*(t**3),r*(t**4),2*r*s*t]
    wm=[z,s,s*(z**2)]
    mp,Cp=loc._weighted_stats(wp,Vp); mm,Cm=loc._weighted_stats(wm,Vm); d=loc._vsub(mp,mm)
    Aeven=sum(wp,Jet4(0))
    # Bodd before multiplication by physical y; common factor 2*A*t^(2-sqrt3)=2*A*t^2/z.
    Bodd=2*A*(t**2)*(z+s+s*(z**2))/z
    odd=y*Bodd; e=odd/(Aeven+odd); q=Jet4(1)-e; qe=q*e
    # Analytic physical hull: 0<=q(1-q)<=1/4. Keep derivative/Hessian enclosures, tighten only value.
    qe=Jet4(QI(max(F(0),qe.v.lo),min(F(1,4),qe.v.hi)),qe.g,qe.H)
    eta=Jet4(1)-qe*c
    M=loc._madd(loc._mscale(q,Cp),loc._mscale(e,Cm)); M=loc._madd(M,loc._mscale(qe/eta,loc._outer(d)))
    K=[[(Jet4(threshold)-M[i][j] if i==j else -M[i][j]) for j in range(3)] for i in range(3)]
    aa,b,cc=K[0];dd=K[1][1];ee=K[1][2];ff=K[2][2]
    P=aa*dd*ff+2*b*cc*ee-aa*ee*ee-dd*cc*cc-ff*b*b
    P1=(aa*dd-b*b)+(aa*ff-cc*cc)+(dd*ff-ee*ee); c2=aa+dd+ff
    return P,P1,c2,e

def center_cell(cell):
    mids=[];rads=[]
    for lo,hi in cell:
      lo,hi=F(lo),F(hi);m=(lo+hi)/2;mids.append(m);rads.append((hi-lo)/2)
    return tuple((m,m) for m in mids),rads

def absmax(I):return max(abs(I.lo),abs(I.hi))
def range_from_jets(full,cen,rads_internal):
    lin=cen.v
    for i,r in enumerate(rads_internal):lin=lin+cen.g[i]*QI(-r,r)
    rem=F(0)
    for i in range(N):
      for j in range(N): rem += F(1,2)*absmax(full.H[i][j])*rads_internal[i]*rads_internal[j]
    return QI(lin.lo-rem,lin.hi+rem)

def ranges(cell,a=F(67,250)):
    cell=tuple((F(lo),F(hi)) for lo,hi in cell)
    # Need radii in the internal (sqrt r,u,v,y) variables, not compact r,s,t,y.
    (rlo,rhi),(slo,shi),(tlo,thi),(ylo,yhi)=cell
    ar=sqrt_interval(QI(rlo,rhi),60); internal_rads=[(ar.hi-ar.lo)/2,(shi-slo)/2,(thi-tlo)/2,(yhi-ylo)/2]
    full=model_jet(cell,a); cc,_=center_cell(cell); cen=model_jet(cc,a)
    return {'P':range_from_jets(full[0],cen[0],internal_rads),'P1':range_from_jets(full[1],cen[1],internal_rads),'c2':range_from_jets(full[2],cen[2],internal_rads),'e':full[3].v}

def certify(cell,a=F(67,250)):
    R=ranges(cell,a); reason='P_NONPOS' if R['P'].hi<=0 else ('P1_NONNEG' if R['P1'].lo>=0 else None); ok=R['c2'].lo>0 and reason is not None
    return {'ok':ok,'status':'INTERVAL_CERTIFIED' if ok else 'FAILED','reason':reason,'P_hi':float(R['P'].hi),'P1_lo':float(R['P1'].lo),'c2_lo':float(R['c2'].lo),'e':[float(R['e'].lo),float(R['e'].hi)]}
