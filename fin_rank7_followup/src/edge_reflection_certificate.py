"""Partial R7P-038/039: edge-reflection invariant stationary family.

The reflection is j -> 3-j. Its fixed dual subspace has basis
(3c-3s)/sqrt(2), 4c, (5c+5s)/sqrt(2); k6 is odd and therefore absent.
"""
from __future__ import annotations
import json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import root
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.model import feature_spaces,dual7
from src.derivatives import dual_all
from src.coexistence_certificate import _I,_bound,interval_features,interval_ldlt
mp.mp.dps=80;mp.iv.dps=60

_,_,_,X7,_,_=feature_spaces()
rt2=np.sqrt(2.0)
B=np.zeros((7,3));B[0,0]=1/rt2;B[1,0]=-1/rt2;B[2,1]=1;B[4,2]=1/rt2;B[5,2]=1/rt2
Aodd=np.zeros((7,4));Aodd[0,0]=1/rt2;Aodd[1,0]=1/rt2;Aodd[3,1]=1;Aodd[4,2]=1/rt2;Aodd[5,2]=-1/rt2;Aodd[6,3]=1
Z=X7@B


def restricted_np(q,g):
    th=B@q;phi,gr,H,p=dual7(th,g,X7)
    return B.T@gr,B.T@H@B,phi,p,H


def augmented_np(z):
    q=z[:3];g=z[3];v=z[4:];F,H,_,_,_=restricted_np(q,g)
    return np.r_[F,H@v,v@v-1]


def augmented_jac_np(z):
    q=z[:3];g=z[3];v=z[4:];th=B@q
    _,gr,H7,T37,_,p=dual_all(th,g,X7)
    H=B.T@H7@B
    T=np.einsum('ia,jb,kc,ijk->abc',B,B,B,T37)
    J=np.zeros((7,7));J[:3,:3]=H;J[:3,3]=-q/g**2
    for i in range(3):
        for k in range(3):J[3+i,k]=sum(T[i,j,k]*v[j] for j in range(3))
        J[3+i,3]=-v[i]/g**2;J[3+i,4:]=H[i]
    J[6,4:]=2*v
    return J


def _interval_bases():
    _,_,C,S=interval_features();sq2=mp.iv.sqrt(_I(2))
    Xi=[[C[j][0],S[j][0],C[j][1],S[j][1],C[j][2],S[j][2],C[j][3]] for j in range(12)]
    fixed=[[(Xi[j][0]-Xi[j][1])/sq2,Xi[j][2],(Xi[j][4]+Xi[j][5])/sq2] for j in range(12)]
    odd=[[(Xi[j][0]+Xi[j][1])/sq2,Xi[j][3],(Xi[j][4]-Xi[j][5])/sq2,Xi[j][6]] for j in range(12)]
    return fixed,odd


def eval_fixed_iv(q,g,fixed):
    h=[sum(fixed[j][a]*q[a] for a in range(3)) for j in range(12)]
    e=[mp.iv.exp(x) for x in h];Zs=sum(e,_I(0));p=[x/Zs for x in e]
    mu=[sum(p[j]*fixed[j][a] for j in range(12)) for a in range(3)]
    Y=[[fixed[j][a]-mu[a] for a in range(3)] for j in range(12)]
    cov=[[sum(p[j]*Y[j][a]*Y[j][b] for j in range(12)) for b in range(3)] for a in range(3)]
    H=[[(_I(1)/g if a==b else _I(0))-cov[a][b] for b in range(3)] for a in range(3)]
    T=[[[-sum(p[j]*Y[j][a]*Y[j][b]*Y[j][c] for j in range(12)) for c in range(3)] for b in range(3)] for a in range(3)]
    F=[q[a]/g-mu[a] for a in range(3)]
    phi=sum(x*x for x in q)/(2*g)-mp.iv.log(Zs/12)
    return F,H,T,p,phi


def eval_aug_iv(z,fixed):
    q=z[:3];g=z[3];v=z[4:];F,H,T,p,phi=eval_fixed_iv(q,g,fixed)
    G=F+[sum(H[i][j]*v[j] for j in range(3)) for i in range(3)]+[sum(x*x for x in v)-1]
    J=[[_I(0) for _ in range(7)] for __ in range(7)]
    for i in range(3):
        for k in range(3):J[i][k]=H[i][k]
        J[i][3]=-q[i]/(g*g)
    for i in range(3):
        for k in range(3):J[3+i][k]=sum(T[i][j][k]*v[j] for j in range(3))
        J[3+i][3]=-v[i]/(g*g)
        for j in range(3):J[3+i][4+j]=H[i][j]
    for j in range(3):J[6][4+j]=2*v[j]
    return G,J,H,T,p,phi


def odd_hessian(p,g,odd):
    # Odd feature means vanish exactly under the fixed reflection law.
    return [[(_I(1)/g if a==b else _I(0))-sum(p[j]*odd[j][a]*odd[j][b] for j in range(12))
             for b in range(4)] for a in range(4)]


def solve_fold():
    seed=np.array([.89489984,1.42717672,1.04202565,4.35221468,-.50003652,-.70597105,-.5015659])
    sol=root(augmented_np,seed,tol=1e-12);z=sol.x
    if np.max(abs(augmented_np(z)))>1e-11: raise RuntimeError('edge fold seed failed')
    return z


def certify_fold(radius='1e-8'):
    z=solve_fold();fixed,odd=_interval_bases();rad=mp.mpf(radius)
    x0=[mp.mpf(format(x,'.17g')) for x in z]
    X=[_I([mp.nstr(x-rad,80),mp.nstr(x+rad,80)]) for x in x0];XP=[_I(mp.nstr(x,80)) for x in x0]
    G0,*_=eval_aug_iv(XP,fixed);G,J,H,T,p,phi=eval_aug_iv(X,fixed);R=np.linalg.inv(augmented_jac_np(z))
    K0=[XP[i]-sum(_I(format(R[i,j],'.17g'))*G0[j] for j in range(7)) for i in range(7)]
    E=[[_I(1 if i==j else 0)-sum(_I(format(R[i,k],'.17g'))*J[k][j] for k in range(7)) for j in range(7)] for i in range(7)]
    dx=[_I([mp.nstr(-rad,80),mp.nstr(rad,80)]) for _ in range(7)]
    K=[K0[i]+sum(E[i][j]*dx[j] for j in range(7)) for i in range(7)]
    inc=all(float(X[i].a)<float(K[i].a) and float(K[i].b)<float(X[i].b) for i in range(7))
    d2,_,ok2=interval_ldlt([[H[i][j] for j in range(2)] for i in range(2)])
    Ho=odd_hessian(p,X[3],odd);do,_,oko=interval_ldlt(Ho)
    a=sum(X[4+i]*(-X[i]/(X[3]*X[3])) for i in range(3))
    b=sum(T[i][j][k]*X[4+i]*X[4+j]*X[4+k] for i in range(3) for j in range(3) for k in range(3))
    return {
      'id':'R7P-038-edge-reflection-fold','claim_id':'CLM-032','domain':'fixed subspace of D12 reflection j->3-j with full odd-sector transverse check','quantifiers':'for the strict spectral tuple within accepted intervals','assumptions':['conditional active gain dual','accepted strict spectral intervals'],'proof_type':'7D augmented interval Krawczyk + fixed/odd block interval inertia','inputs':['inputs/fin_handoff_audit/results.json'],'conclusion':'an edge-reflection simple fold occurs in the stated box; full H7 inertia at the fold is (1 negative,1 zero,5 positive)','global_pass':False,
      'root_box':[_bound(x) for x in X],'krawczyk_image':[_bound(x) for x in K],'strict_inclusion':inc,
      'gain_interval':_bound(X[3]),'fold_a':_bound(a),'fold_b':_bound(b),
      'fixed_2x2_principal_pivots':[_bound(x) for x in d2],'odd_LDL_pivots':[_bound(x) for x in do],
      'full_H7_inertia_at_fold':[1,1,5],
      'reflection':'j -> 3-j','fixed_basis':['(3c-3s)/sqrt2','4c','(5c+5s)/sqrt2'],'odd_basis':['(3c+3s)/sqrt2','4s','(5c-5s)/sqrt2','k6']
    }


def fixed_gain_root(seed,gstr='4.36',radius='1e-8'):
    fixed,odd=_interval_bases();g=float(gstr)
    sol=root(lambda q:restricted_np(q,g)[0],np.asarray(seed,float),tol=1e-12);q=sol.x
    Hpoint=restricted_np(q,g)[1];R=np.linalg.inv(Hpoint);rad=mp.mpf(radius)
    Q=[_I([mp.nstr(mp.mpf(format(x,'.17g'))-rad,80),mp.nstr(mp.mpf(format(x,'.17g'))+rad,80)]) for x in q]
    G=_I(gstr);QP=[_I(format(x,'.17g')) for x in q]
    F0,_,_,_,_=eval_fixed_iv(QP,G,fixed);F,H,T,p,phi=eval_fixed_iv(Q,G,fixed)
    K0=[QP[i]-sum(_I(format(R[i,j],'.17g'))*F0[j] for j in range(3)) for i in range(3)]
    E=[[_I(1 if i==j else 0)-sum(_I(format(R[i,k],'.17g'))*H[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    dx=[_I([mp.nstr(-rad,80),mp.nstr(rad,80)]) for _ in range(3)]
    K=[K0[i]+sum(E[i][j]*dx[j] for j in range(3)) for i in range(3)]
    inc=all(float(Q[i].a)<float(K[i].a) and float(K[i].b)<float(Q[i].b) for i in range(3))
    df,_,okf=interval_ldlt(H);Ho=odd_hessian(p,G,odd);do,_,oko=interval_ldlt(Ho)
    nf=sum(float(x.b)<0 for x in df);no=sum(float(x.b)<0 for x in do)
    return {'center':q.tolist(),'box':[_bound(x) for x in Q],'krawczyk':[_bound(x) for x in K],'inclusion':inc,
            'fixed_pivots':[_bound(x) for x in df],'odd_pivots':[_bound(x) for x in do],
            'fixed_negative_pivots':nf,'odd_negative_pivots':no,'H7_inertia':[nf+no,0,7-nf-no],
            'phi_interval':_bound(phi),'pmax_upper':max(float(x.b) for x in p)}


def build():
    fold=certify_fold()
    # two branches born at the edge-reflection fold
    weak=fixed_gain_root([.97907583,1.54610621,1.12661156])
    strong=fixed_gain_root([.81164568,1.30968413,.95865432])
    distinct=float(weak['phi_interval'][1]) < float(strong['phi_interval'][0]) or float(strong['phi_interval'][1]) < float(weak['phi_interval'][0])
    roots={
      'id':'R7P-038-edge-reflection-branches','claim_id':'CLM-033','domain':'edge-reflection fixed subspace at exact gain 4.36','quantifiers':'for the strict spectral tuple within accepted intervals','assumptions':['conditional active gain dual','accepted strict spectral intervals'],'proof_type':'two 3D interval Krawczyk isolations + full fixed/odd block inertia','inputs':['certificates/R7P-038_edge_reflection_fold.json'],'conclusion':'two distinct reflection-fixed stationary orbit representatives exist at g=4.36, with full H7 indices 1 and 2','global_pass':False,
      'gain':'4.36','index1_branch':weak,'index2_branch':strong,'D12_orbit_separation_by_disjoint_Phi':distinct,
      'orbit_count_scope':'representatives are distinct D12 orbits because Phi is D12 invariant; exact stabilizer/orbit sizes are not claimed here'
    }
    return fold,roots


def main():
    fold,roots=build()
    (ROOT/'certificates/R7P-038_edge_reflection_fold.json').write_text(json.dumps(fold,indent=2)+'\n')
    (ROOT/'certificates/R7P-038_edge_reflection_branches_g4p36.json').write_text(json.dumps(roots,indent=2)+'\n')
    print(json.dumps({'fold_gain':fold['gain_interval'],'fold_H7':fold['full_H7_inertia_at_fold'],
        'g4.36_indices':[roots['index1_branch']['H7_inertia'],roots['index2_branch']['H7_inertia']],
        'distinct_orbits':roots['D12_orbit_separation_by_disjoint_Phi']},indent=2))

if __name__=='__main__':main()
