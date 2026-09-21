from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
import sys, numpy as np

ROOT=Path(__file__).resolve().parents[1]
R7N=Path('/mnt/data/r7n_full_restore/FIN_R7N_HANDOFF_20260920')
H=R7N/'inputs/FR223_20260916'
sys.path[:0]=[str(R7N/'src'),str(H/'src')]

from intervals import QI,sqrt_interval
import off_face
import off_face_local as loc
import generic_threshold_shifted as gs
import compression_interval_probe as old

TAU=F(67,250)
N=4
Jet4=loc.Jet4

def _absmax(I): return max(abs(I.lo),abs(I.hi))

def _range_from_jets(full,cen,rads):
    lin=cen.v
    for i,r in enumerate(rads): lin=lin+cen.g[i]*QI(-r,r)
    rem=F(0)
    for i in range(N):
        for j in range(N):
            rem += F(1,2)*_absmax(full.H[i][j])*rads[i]*rads[j]
    return QI(lin.lo-rem,lin.hi+rem)

def _internal_bounds(cell):
    cell=tuple((F(lo),F(hi)) for lo,hi in cell)
    (rlo,rhi),(slo,shi),(tlo,thi),(ylo,yhi)=cell
    A=sqrt_interval(QI(rlo,rhi),60)
    return ((A.lo,A.hi),(F(1)-shi,F(1)-slo),(F(1)-thi,F(1)-tlo),(ylo,yhi))

def _midpoint_bounds(bounds):
    mids=[]; rads=[]
    for lo,hi in bounds:
        m=(lo+hi)/2
        mids.append((m,m)); rads.append((hi-lo)/2)
    return tuple(mids),rads

def _basis_center(cell,den=100000):
    Bn,den,eigs,rows,det=old.center_basis(tuple((F(a),F(b)) for a,b in cell),den)
    B=[[F(int(Bn[i,j]),den) for j in range(3)] for i in range(4)]
    x=np.array([float((F(a)+F(b))/2) for a,b in cell])
    p=off_face.p_from_aligned_compact(*x)
    pa=np.array([p[0],p[4]+p[8],p[6],p[2]+p[10],p[3]+p[9],p[5]+p[7],p[1]+p[11]])
    zm=np.zeros((7,3))
    for i in range(7):
        for k in range(3):
            z=0.0
            for r in range(4):
                z += float(B[r][k])*float((old.OBS[i][r].lo+old.OBS[i][r].hi)/2)
            zm[i,k]=z
    mu=pa@zm
    c=[F(format(float(v),'.16g')) for v in mu]
    return B,Bn,den,eigs,rows,det,c

def _weights_jet(bounds):
    (Alo,Ahi),(ulo,uhi),(vlo,vhi),(ylo,yhi)=bounds
    A=Jet4.var(QI(Alo,Ahi),0)
    u=Jet4.var(QI(ulo,uhi),1)
    v=Jet4.var(QI(vlo,vhi),2)
    y=Jet4.var(QI(ylo,yhi),3)
    s=Jet4(1)-u; t=Jet4(1)-v
    rt3=sqrt_interval(QI(3),40)
    z=gs._z_of_v(v,vlo,vhi,rt3)
    r=A*A
    return [
        Jet4(1),
        2*s*(t**3),
        r*(t**4),
        2*r*s*t,
        2*A*(t**2)*y,
        2*A*s*(t**2)*y/z,
        2*A*s*(t**2)*z*y,
    ]

def _moment_jets(bounds,B,c):
    w=_weights_jet(bounds)
    D=sum(w,Jet4(0)); invD=D.inv()
    Z=[]
    for i in range(7):
        row=[]
        for k in range(3):
            zz=QI(0)
            for r in range(4): zz += B[r][k]*old.OBS[i][r]
            row.append(zz-QI(c[k]))
        Z.append(row)
    E=[[None]*3 for _ in range(3)]
    for a in range(3):
        for b in range(a,3):
            Nn=Jet4(0)
            for i in range(7): Nn += w[i]*(Z[i][a]*Z[i][b])
            E[a][b]=Nn*invD; E[b][a]=E[a][b]
    return E

def _pd_result(K):
    d1=K[0][0]
    d2=K[0][0]*K[1][1]-K[0][1]*K[0][1]
    d3=(K[0][0]*K[1][1]*K[2][2]+2*K[0][1]*K[0][2]*K[1][2]
        -K[0][0]*K[1][2]*K[1][2]-K[1][1]*K[0][2]*K[0][2]-K[2][2]*K[0][1]*K[0][1])
    sylv=d1.lo>0 and d2.lo>0 and d3.lo>0
    g=[]
    for a in range(3):
        offsum=sum(max(abs(K[a][b].lo),abs(K[a][b].hi)) for b in range(3) if b!=a)
        g.append(K[a][a].lo-offsum)
    gersh=min(g)>0
    return sylv,gersh,d1,d2,d3,g

def certify(cell,den=100000):
    cell=tuple((F(a),F(b)) for a,b in cell)
    B,Bn,den,eigs,rows,det,c=_basis_center(cell,den)
    base={'rank_rows':rows,'rank_det':det,'basis_den':den,'basis_num':Bn.tolist(),'center_eigs':eigs.tolist(),'center_c':[str(x) for x in c]}
    if det==0: return {'ok':False,'reason':'RANK',**base}
    bounds=_internal_bounds(cell); cb,rads=_midpoint_bounds(bounds)
    full=_moment_jets(bounds,B,c); cen=_moment_jets(cb,B,c)
    Er=[[None]*3 for _ in range(3)]
    for a in range(3):
        for b in range(a,3):
            Er[a][b]=_range_from_jets(full[a][b],cen[a][b],rads); Er[b][a]=Er[a][b]
    Gram=[[sum(B[r][a]*B[r][b] for r in range(4)) for b in range(3)] for a in range(3)]
    K=[[QI(TAU*Gram[a][b])-Er[a][b] for b in range(3)] for a in range(3)]
    sylv,gersh,d1,d2,d3,g=_pd_result(K)
    ok=sylv or gersh
    return {'ok':bool(ok),'reason':'PHYSICAL_CENTERED_SYLVESTER_PD' if sylv else ('PHYSICAL_CENTERED_GERSHGORIN_PD' if gersh else 'FAILED'),
            **base,'gersh_lower':float(min(g)),
            'd1':[float(d1.lo),float(d1.hi)],'d2':[float(d2.lo),float(d2.hi)],'d3':[float(d3.lo),float(d3.hi)],
            'entry_width_max':max(float(Er[a][b].hi-Er[a][b].lo) for a in range(3) for b in range(3))}
