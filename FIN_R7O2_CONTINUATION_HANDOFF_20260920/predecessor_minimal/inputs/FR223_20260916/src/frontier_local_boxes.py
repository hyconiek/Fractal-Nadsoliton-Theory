"""Post-handoff anisotropic local-box certificates around the R7P-068 equality point.

This module generalizes the original R7P-068 interval-AD checker from a single
common radius rho to independent bounds in the four local coordinates

    x = r-r_*,   u = 1-s,   v = 1-t,   e = 1-q,

with r=e^{-2J3}, s=e^{-3J4/2}, t=e^{-J5/2}.

The proof logic is unchanged: for the 3x3 Schur-reduced covariance Mtilde, with
P=det(sigma I-Mtilde), P1=dP/dsigma and c2=P''/2, R7P-044 gives
lambda2<=sigma whenever c2>0 and (P<=0 or P1>=0).  Second-order rational
interval AD and two cone-Schur endpoint tests certify that disjunction.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import json, math

import off_face_local as loc
from intervals import QI, sqrt_interval
import boundary_ising as bi

ROOT=Path(__file__).resolve().parents[1]
N=4
Jet4=loc.Jet4


def _model_box(rx,ru,rv,re):
    rx,ru,rv,re=F(rx),F(ru),F(rv),F(re)
    if min(rx,ru,rv,re)<=0: raise ValueError('all radii must be positive')
    if max(rx,ru,rv,re)>=F(1,10): raise ValueError('outside local regime')

    L=bi.strict_intervals(); l3,l4,l5,l6=L[3],L[4],L[5],L[6]
    sigma=(2*l3*(l4+l5)-l4*l5)/(24*l3)
    c=l6/(3*sigma)
    tau=sqrt_interval(((2*l3-l4)*(2*l3-l5))/(4*l3*l3),40)
    rstar=(QI(1)-tau)/(QI(1)+tau)
    a3=sqrt_interval(l3/6,40); a4=sqrt_interval(l4/6,40); a5=sqrt_interval(l5/6,40)
    rt3=sqrt_interval(QI(3),40)

    Vp=[[a3,a4,a5],[a3,-a4/2,-a5/2],[-a3,a4,-a5],[-a3,-a4/2,a5/2]]
    Vm=[[QI(0),a4,QI(0)],[QI(0),-a4/2,rt3*a5/2],[QI(0),-a4/2,-rt3*a5/2]]

    x=Jet4.var(QI(-rx,rx),0)
    u=Jet4.var(QI(0,ru),1)
    v=Jet4.var(QI(0,rv),2)
    e=Jet4.var(QI(0,re),3)
    r=Jet4(rstar)+x; s=Jet4(1)-u; t=Jet4(1)-v

    # Same rigorous physical enclosure for z=t^sqrt(3) as R7P-068, but with
    # the independent v-radius rv.
    tlo=F(1)-rv
    zv=QI(tlo*tlo,1)
    zg=[QI(0) for _ in range(N)]; zg[2]=QI(-rt3.hi,-rt3.lo*tlo)
    zH=[[QI(0) for _ in range(N)] for __ in range(N)]
    af=rt3*(rt3-QI(1)); zH[2][2]=af*QI(1,F(1)/tlo)
    z=Jet4(zv,zg,zH)

    wp=[Jet4(1),2*s*(t**3),r*(t**4),2*r*s*t]
    wm=[z,s,s*(z**2)]
    mp,Cp=loc._weighted_stats(wp,Vp); mm,Cm=loc._weighted_stats(wm,Vm)
    d=loc._vsub(mp,mm)
    q=Jet4(1)-e; qe=q*e; eta=Jet4(1)-qe*c
    M=loc._madd(loc._mscale(q,Cp),loc._mscale(e,Cm))
    M=loc._madd(M,loc._mscale(qe/eta,loc._outer(d)))

    K=[[(Jet4(sigma)-M[i][j] if i==j else -M[i][j]) for j in range(3)] for i in range(3)]
    a,b,cc=K[0]; dd=K[1][1]; ee=K[1][2]; ff=K[2][2]
    P=a*dd*ff+2*b*cc*ee-a*ee*ee-dd*cc*cc-ff*b*b
    P1=(a*dd-b*b)+(a*ff-cc*cc)+(dd*ff-ee*ee)
    c2=a+dd+ff
    return rstar,P,P1,c2


def raw_box(rx,ru,rv,re,margin=F(1,10000)):
    rx,ru,rv,re,margin=F(rx),F(ru),F(rv),F(re),F(margin)
    rstar,P,P1,c2=_model_box(rx,ru,rv,re)
    g=P1.g; H=P.H
    signs={
      'P1_x_negative':g[0].hi<0,
      'P1_u_negative':g[1].hi<0,
      'P1_v_negative':g[2].hi<0,
      'P1_e_positive':g[3].lo>0,
      'P_ee_positive':H[3][3].lo>0,
      'c2_positive':c2.v.lo>0,
    }
    alpha=[(-g[i])/g[3] for i in range(3)]
    HB=[row[:3] for row in H[:3]]
    HE=[]
    for i in range(3):
        row=[]
        for j in range(3):
            row.append(H[i][j]+H[i][3]*alpha[j]+alpha[i]*H[3][j]
                       +alpha[i]*alpha[j]*H[3][3])
        HE.append(row)
    okB,sB=loc._schur_cone(HB,margin)
    okE,sE=loc._schur_cone(HE,margin)
    ok=all(signs.values()) and okB and okE
    return {
      'status':'INTERVAL_CERTIFIED' if ok else 'FAILED',
      'radii':{'x':rx,'u':ru,'v':rv,'e':re},
      'margin':margin,'rstar':rstar,'P':P,'P1':P1,'c2':c2,'signs':signs,
      'alpha':alpha,'boundary_ok':okB,'endpoint_ok':okE,
      'boundary_schur':sB,'endpoint_schur':sE,
    }


def _short(I,digits=15):
    return [format(float(I.lo),f'.{digits}g'),format(float(I.hi),f'.{digits}g')]


def fr10_record():
    good=raw_box(F(1,8192),F(1,4800),F(1,4800),F(1,3072))
    bad_e=raw_box(F(1,8192),F(1,4800),F(1,4800),F(1,3052))
    bad_uv=raw_box(F(1,8192),F(1,4700),F(1,4700),F(1,3072))
    assert good['status']=='INTERVAL_CERTIFIED'
    assert bad_e['status']=='FAILED' and not bad_e['boundary_ok']
    assert bad_uv['status']=='FAILED' and not bad_uv['boundary_ok']
    return {
      'id':'FR10-diagonal-uv-local-box',
      'status':'INTERVAL_CERTIFIED',
      'proof_type':'anisotropic second-order rational interval AD + R7P-044 characteristic-inertia disjunction',
      'domain':{
        'abs_r_minus_rstar':'<=1/8192',
        'u=1-exp(-3J4/2)':'<=1/4800',
        'v=1-exp(-J5/2)':'<=1/4800',
        'e=1-q_even':'<=1/3072',
        'J_fields':'J3,J4,J5,J6>=0 with the physical shared-field odd law',
      },
      'conclusion':'lambda2(Mtilde)<=sigma_* and hence lambda2(M4)<=sigma_* throughout the declared anisotropic box.',
      'strict_checks':{
        'signs':good['signs'],
        'boundary_schur_pass':good['boundary_ok'],
        'P1_zero_endpoint_schur_pass':good['endpoint_ok'],
        'boundary_Nvv_interval':_short(good['boundary_schur'][3]),
        'endpoint_Nvv_interval':_short(good['endpoint_schur'][3]),
        'c2_interval':_short(good['c2'].v),
      },
      'negative_controls':{
        'larger_e_1_over_3052':{
          'status':bad_e['status'],
          'boundary_Nvv_interval':_short(bad_e['boundary_schur'][3]),
          'interpretation':'checker failure only; not a physical counterexample'},
        'larger_uv_1_over_4700':{
          'status':bad_uv['status'],
          'boundary_Nvv_interval':_short(bad_uv['boundary_schur'][3]),
          'interpretation':'checker failure only; not a physical counterexample'},
      },
      'scope':'Local equality-neighborhood theorem in the physical four-amplitude model only. No full positive-orthant or full-seven-coordinate transfer.',
    }


def write_fr10(path=None):
    rec=fr10_record()
    if path is None: path=ROOT/'results/FR10_diagonal_uv_local_box.json'
    Path(path).write_text(json.dumps(rec,indent=2)+'\n')
    return rec

if __name__=='__main__':
    print(json.dumps(write_fr10(),indent=2))
