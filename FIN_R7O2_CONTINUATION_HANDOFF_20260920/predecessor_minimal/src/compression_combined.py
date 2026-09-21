from pathlib import Path
from fractions import Fraction as F
import sys,itertools,numpy as np
ROOT=Path('/mnt/data/fin_rank7_next_campaign'); H=ROOT/'inputs/FR223_20260916'; sys.path[:0]=[str(ROOT/'src'),str(H/'src')]
from intervals import QI
import target_p_trace_cover as tc
import target_p_trace_tight_rounded as tr
import compression_interval_probe as old
import off_face
TAU=F(67,250)

def _center_basis_data(cell,den=100000):
    Bn,den,eigs,rows,det=old.center_basis(cell,den)
    B=[[F(int(Bn[i,j]),den) for j in range(3)] for i in range(4)]
    x=np.array([float((a+b)/2) for a,b in cell]); p=off_face.p_from_aligned_compact(*x)
    pa=np.array([p[0],p[4]+p[8],p[6],p[2]+p[10],p[3]+p[9],p[5]+p[7],p[1]+p[11]])
    zm=np.zeros((7,3))
    for i in range(7):
      for k in range(3):
        z=0.0
        for r in range(4): z += float(B[r][k])*float((old.OBS[i][r].lo+old.OBS[i][r].hi)/2)
        zm[i,k]=z
    mu=pa@zm
    c=[F(format(float(v),'.16g')) for v in mu]
    return B,Bn,den,eigs,rows,det,c

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
    B,Bn,den,eigs,rows,det,c=_center_basis_data(cell,den)
    base={'rank_rows':rows,'rank_det':det,'basis_den':den,'basis_num':Bn.tolist(),'center_eigs':eigs.tolist()}
    if det==0:return {'ok':False,'reason':'RANK',**base}
    raw=tc.weight_bounds(cell)
    wb=[(tr.floorq(lo),tr.ceilq(hi)) for lo,hi in raw]
    D=QI(sum(a for a,b in wb),sum(b for a,b in wb))
    Z=[]
    for i in range(7):
      row=[]
      for k in range(3):
        z=QI(0)
        for r in range(4): z += B[r][k]*old.OBS[i][r]
        row.append(z-QI(c[k]))
      Z.append(row)
    E=[[QI(0) for _ in range(3)] for __ in range(3)]
    for a in range(3):
      for b in range(a,3):
        N=QI(0)
        for i in range(7): N += QI(wb[i][0],wb[i][1])*Z[i][a]*Z[i][b]
        E[a][b]=N/D; E[b][a]=E[a][b]
    Gram=[[sum(B[r][a]*B[r][b] for r in range(4)) for b in range(3)] for a in range(3)]
    K=[[QI(TAU*Gram[a][b])-E[a][b] for b in range(3)] for a in range(3)]
    sylv,gersh,d1,d2,d3,g=_pd_result(K)
    if sylv or gersh:
      return {'ok':True,'stage':'CENTERED','reason':'CENTERED_SYLVESTER_PD' if sylv else 'CENTERED_GERSHGORIN_PD',
              **base,'center_c':[str(x) for x in c],'gersh_lower':float(min(g)),
              'd1':[float(d1.lo),float(d1.hi)],'d2':[float(d2.lo),float(d2.hi)],'d3':[float(d3.lo),float(d3.hi)]}
    # Fall through to historical pairwise covariance compression, reusing B and raw weight bounds.
    K=[[QI(TAU*Gram[a][b]) for b in range(3)] for a in range(3)]
    for i,j in itertools.combinations(range(7),2):
        q=old.pair_prob_bounds(raw,i,j)
        dz=[]
        for k in range(3):
            z=QI(0)
            for r in range(4): z += B[r][k]*(old.OBS[i][r]-old.OBS[j][r])
            dz.append(z)
        for a in range(3):
            for b in range(a,3):
                K[a][b]=K[a][b]-q*dz[a]*dz[b]
                if a!=b:K[b][a]=K[a][b]
    sylv,gersh,d1,d2,d3,g=_pd_result(K)
    return {'ok':bool(sylv or gersh),'stage':'PAIRWISE','reason':'SYLVESTER_PD' if sylv else ('GERSHGORIN_PD' if gersh else 'FAILED'),
            **base,'gersh_lower':float(min(g)),'d1':[float(d1.lo),float(d1.hi)],'d2':[float(d2.lo),float(d2.hi)],'d3':[float(d3.lo),float(d3.hi)]}
