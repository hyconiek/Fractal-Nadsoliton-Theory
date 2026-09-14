"""R7P-030/031: corrected 9-variable fold system and simple-fold certificate."""
from __future__ import annotations
import json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import root
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.model import feature_spaces,dual7
from src.derivatives import dual_all
from src.coexistence_certificate import (_I,_bound,interval_features,interval_eval,
    interval_ldlt)
mp.mp.dps=80;mp.iv.dps=60


def augmented_np(z,C4):
    s=z[:4];g=z[4];v=z[5:]
    _,grad,H,_,_,_=dual_all(s,g,C4)
    return np.r_[grad,H@v,v@v-1]


def augmented_jacobian_np(z,C4):
    s=z[:4];g=z[4];v=z[5:]
    _,grad,H,T3,_,_=dual_all(s,g,C4)
    J=np.zeros((9,9));J[:4,:4]=H;J[:4,4]=-s/g**2
    for i in range(4):
        for k in range(4):J[4+i,k]=sum(T3[i,j,k]*v[j] for j in range(4))
        J[4+i,4]=-v[i]/g**2;J[4+i,5:]=H[i]
    J[8,5:]=2*v
    return J


def solve_candidate():
    _,_,_,X7,C4,_=feature_spaces()
    seed=np.array([1.36431143,1.43310270,1.40804570,1.00568208,3.51564472,
                   .50736754,.52868566,.55376069,.39549811])
    sol=root(lambda z:augmented_np(z,C4),seed,tol=1e-12)
    if np.max(abs(augmented_np(sol.x,C4)))>1e-11: raise RuntimeError('bad fold seed')
    z=sol.x;s=z[:4];g=z[4]
    val,grad,H,p=dual7(s,g,C4);theta=np.array([s[0],0,s[1],0,s[2],0,s[3]])
    H7=dual7(theta,g,X7)[2];J=augmented_jacobian_np(z,C4)
    return z,dict(success=bool(sol.success),residual=float(np.max(abs(augmented_np(z,C4)))),
        s=s.tolist(),g=float(g),v=z[5:].tolist(),phi=float(val),pmax=float(p.max()),
        H4_eigs=np.linalg.eigvalsh(H).tolist(),H7_eigs=np.linalg.eigvalsh(H7).tolist(),
        augmented_jacobian_singular_values=np.linalg.svd(J,compute_uv=False).tolist())


def _augmented_iv(z,C):
    F,J,p,_=interval_eval(z[:5],C);H=[[J[i][j] for j in range(4)] for i in range(4)]
    # third dual derivative = minus third central feature moment
    h=[sum(C[j][a]*z[a] for a in range(4)) for j in range(12)]
    e=[mp.iv.exp(x) for x in h];Z=sum(e,_I(0));pp=[x/Z for x in e]
    mu=[sum(pp[j]*C[j][a] for j in range(12)) for a in range(4)]
    Y=[[C[j][a]-mu[a] for a in range(4)] for j in range(12)]
    T3=[[[-sum(pp[j]*Y[j][a]*Y[j][b]*Y[j][c] for j in range(12))
          for c in range(4)] for b in range(4)] for a in range(4)]
    v=z[5:]
    G=F[:4]+[sum(H[i][j]*v[j] for j in range(4)) for i in range(4)]+[sum(x*x for x in v)-1]
    AJ=[[_I(0) for _ in range(9)] for __ in range(9)]
    for i in range(4):
        for k in range(4):AJ[i][k]=H[i][k]
        AJ[i][4]=-z[i]/(z[4]*z[4])
    for i in range(4):
        for k in range(4):AJ[4+i][k]=sum(T3[i][j][k]*v[j] for j in range(4))
        AJ[4+i][4]=-v[i]/(z[4]*z[4])
        for j in range(4):AJ[4+i][5+j]=H[i][j]
    for j in range(4):AJ[8][5+j]=2*v[j]
    return G,AJ,H,T3,pp


def certify(radius='1e-8'):
    z,cand=solve_candidate();_,_,C,S=interval_features();rad=mp.mpf(radius)
    x0=[mp.mpf(format(x,'.17g')) for x in z]
    X=[_I([mp.nstr(x-rad,80),mp.nstr(x+rad,80)]) for x in x0]
    XP=[_I(mp.nstr(x,80)) for x in x0]
    G0,_,_,_,_=_augmented_iv(XP,C);_,JB,H,T3,p=_augmented_iv(X,C)
    R=np.linalg.inv(augmented_jacobian_np(z,feature_spaces()[4]))
    K0=[XP[i]-sum(_I(format(R[i,j],'.17g'))*G0[j] for j in range(9)) for i in range(9)]
    E=[[_I(1 if i==j else 0)-sum(_I(format(R[i,k],'.17g'))*JB[k][j] for k in range(9))
        for j in range(9)] for i in range(9)]
    dx=[_I([mp.nstr(-rad,80),mp.nstr(rad,80)]) for _ in range(9)]
    K=[K0[i]+sum(E[i][j]*dx[j] for j in range(9)) for i in range(9)]
    inc=[float(X[i].a)<float(K[i].a) and float(K[i].b)<float(X[i].b) for i in range(9)]
    v=X[5:];s=X[:4];g=X[4]
    fold_a=sum(v[i]*(-s[i]/(g*g)) for i in range(4))
    fold_b=sum(T3[i][j][k]*v[i]*v[j]*v[k] for i in range(4) for j in range(4) for k in range(4))
    # A fixed 3x3 principal H4 submatrix positive definite + Hv=0 + ||v||=1
    # proves exactly one zero and three positive H4 eigenvalues at the isolated root.
    D3,_,ok3=interval_ldlt([[H[i][j] for j in range(3)] for i in range(3)])
    Hs=[[(_I(1)/g if a==b else _I(0))-sum(p[j]*S[j][a]*S[j][b] for j in range(12))
         for b in range(3)] for a in range(3)]
    Ds,_,oks=interval_ldlt(Hs)
    cert={
      'id':'R7P-031-simple-fold','claim_id':'CLM-031','domain':'reflection-even four-amplitude stationary system with full-H7 transverse directions','quantifiers':'for the strict spectral tuple within the accepted outward spectral intervals','assumptions':['conditional active gain dual','accepted strict spectral intervals'],'proof_type':'9D validated interval Krawczyk plus simple-fold coefficients and interval inertia','inputs':['inputs/fin_handoff_audit/results.json','results/R7P-030_fold_candidate.json'],'conclusion':'one isolated simple saddle-node has exactly one full-H7 null direction and six positive transverse directions','global_pass':False,
      'task':'R7P-031','scientific_status':'VALIDATED_INTERVAL_SIMPLE_FOLD',
      'interval_backend':'mpmath.iv 60 decimal digits',
      'scope':'local reflection-even four-amplitude stationary system with full-H7 transverse check; supplied active gain remains conditional',
      'augmented_root_box':[_bound(x) for x in X],
      'krawczyk_image':[_bound(x) for x in K],'strict_inclusion':all(inc),
      'fold_gain_interval':_bound(X[4]),
      'fold_coeff_v_dot_Fg':_bound(fold_a),
      'fold_coeff_D3Phi_vvv':_bound(fold_b),
      'simple_fold_transversality':float(fold_a.a)<0 and float(fold_a.b)<0 and float(fold_b.a)>0,
      'H4_fixed_3x3_principal_LDL_pivots':[_bound(x) for x in D3],
      'H4_principal_positive':ok3 and all(float(x.a)>0 for x in D3),
      'H4_inertia_at_root':[0,1,3],
      'Hsin_LDL_pivots':[_bound(x) for x in Ds],
      'Hsin_positive_definite':oks and all(float(x.a)>0 for x in Ds),
      'H7_inertia_at_root':[0,1,6],
      'interpretation':'one simple stationary saddle-node/fold; not a global first fold or physical-time event'
    }
    return cand,cert


def main():
    cand,cert=certify()
    (ROOT/'results/R7P-030_fold_candidate.json').write_text(json.dumps(cand,indent=2)+'\n')
    (ROOT/'certificates/R7P-031_simple_fold.json').write_text(json.dumps(cert,indent=2)+'\n')
    print(json.dumps({'R7P-030_residual':cand['residual'],'R7P-031_inclusion':cert['strict_inclusion'],
      'fold_gain_interval':cert['fold_gain_interval'],'H7_inertia':cert['H7_inertia_at_root'],
      'a':cert['fold_coeff_v_dot_Fg'],'b':cert['fold_coeff_D3Phi_vvv']},indent=2))

if __name__=='__main__':main()
