#!/usr/bin/env python3
"""Kinetic-invariant numerical DIAGNOSTICS; proofs are separate. No trajectory simulations.
Uses a hash-preserved source finite-model implementation only for model reconstruction.
"""
from __future__ import annotations
import sys, json, math
from pathlib import Path
import numpy as np
from scipy.linalg import null_space, solve_continuous_lyapunov
base=Path(__file__).resolve().parent
sys.path.insert(0,str(base/'source_context'))
import theory_core as tc
Q=12

def rates(p,g,kin):
    field=g*(tc.A@p)
    delta=field[None,:]-field[:,None]
    if kin=='HB':
        v=tc.softmax_vec(field); return np.tile(v,(Q,1))
    if kin=='MET':return np.exp(np.minimum(0.,delta))/12.
    if kin=='BAR':return (2./12.)*np.exp(-np.logaddexp(0.,-delta))
    raise ValueError(kin)

def finite_rates(n,g,kin):
    N=int(n.sum());out=np.zeros((Q,Q))
    for i in range(Q):
        if n[i]<=0:continue
        fld=(g/N)*(tc.A@(n-np.eye(Q,dtype=int)[i]))
        deltas=fld-fld[i]
        if kin=='HB':probs=tc.softmax_vec(fld)
        elif kin=='MET':probs=np.exp(np.minimum(0.,deltas))/12.
        else:probs=(2./12.)*np.exp(-np.logaddexp(0.,-deltas))
        out[i]=n[i]*probs
        out[i,i]=0
    return out

def main():
    f=tc.solve_fold();g=f['g'];p=f['p'];v=f['v7'];a=f['a'];b=f['b'];A=tc.A;X=tc.X
    Sigma=np.diag(p)-np.outer(p,p)
    r=Sigma@X@v
    T=null_space(np.ones((1,Q)))
    HF=np.diag(1/p)-g*A
    Hess=T.T@HF@T
    eigenh=np.linalg.eigvalsh(Hess)
    aF=float(r@(-A@p)); bF=float(-np.sum(r**3/p**2))
    thetaHB=g*X@v
    lap={}
    for kin in ['HB','MET','BAR']:
        R=rates(p,g,kin)
        c=p[:,None]*R
        K=np.zeros((Q,Q))
        for i in range(Q):
            for j in range(i+1,Q):
                e=np.eye(Q)[j]-np.eye(Q)[i];K += c[i,j]*np.outer(e,e)
        KT=T.T@K@T
        xi=np.linalg.solve(KT,T.T@r)
        rinvt=float(r@(T@xi))
        m=1/rinvt
        ell=m*(T@xi)
        J=-K@HF
        E=null_space(np.vstack([np.ones(Q),ell]))
        B=E.T@J@E
        P=solve_continuous_lyapunov(B.T,-np.eye(10))
        DD=np.zeros((Q,Q))
        for i in range(Q):
            for j in range(Q):
                if i!=j:
                    e=np.eye(Q)[j]-np.eye(Q)[i]
                    DD+=c[i,j]*np.outer(e,e)
        ln2=np.asarray(f['theta7']);
        out={
            'mobility_m_nominal':float(m),
            'mobility_to_HB':float(m/g),
            'noise_variance_soft':float(ell@DD@ell),
            'predicted_soft_noise_variance':float(2*m),
            'tangent_K_eigen_min':float(np.min(np.linalg.eigvalsh(KT))),
            'tangent_Hessian_eigen_min_second':eigenh[:2].tolist(),
            'linearization_soft_residual':float(np.linalg.norm(J@r)),
            'left_soft_residual':float(np.linalg.norm(ell@J@T)),
            'left_right_identity':float(ell@r),
            'fast_gap_nominal':float(-np.max(np.real(np.linalg.eigvals(B)))),
            'lyapunov_P_min_eig_nominal':float(np.min(np.linalg.eigvalsh(P))),
            'lyapunov_residual':float(np.linalg.norm(B.T@P+P@B+np.eye(10))),
            'max_detailed_balance_c_defect':float(np.max(np.abs(c-c.T))),
            'D_equals_2K_error':float(np.max(np.abs(DD-2*K))),
            'J_eq_minus_K_H_error':float(np.linalg.norm(J+K@HF)),
            'ell_diff_from_thermal_norm':float(np.linalg.norm(ell-thetaHB)),
            'ell':ell.tolist(),
            'canonical_d_nominal':float((4/b)**(2/3)/(2*m))
        }
        lap[kin]=out
    # finite count detailed balance test for selected states, N=36, adjacent occupancy swap;
    n=np.array([3]*Q,dtype=int);n[0]+=1;n[1]-=1;N=n.sum();fac=math.factorial(N)
    def logwt(z):
        return math.lgamma(N+1)-sum(math.lgamma(k+1) for k in z)+(g/(2*N))*float(z@A@z)
    checks=[]
    for kin in ['HB','MET','BAR']:
        Rn=finite_rates(n,g,kin)
        for i in range(Q):
            for j in range(Q):
                if i==j or n[i]==0:continue
                nn=n.copy();nn[i]-=1;nn[j]+=1
                rr=finite_rates(nn,g,kin)[j,i]
                res=float((logwt(n)+math.log(Rn[i,j]))-(logwt(nn)+math.log(rr)))
                checks.append({'kin':kin,'i':i,'j':j,'residual':res})
    result={'input_model':'exact-decimal-formula evaluated in binary64',
      'fold_g_nominal':float(g),'fold_a_nominal':float(a),'fold_b_nominal':float(b),
      'a_F_nominal':aF,'b_F_nominal':bF,'a_delta':aF-a,'b_delta':bF-b,
      'Hess_F_tangent_eigenvalues_nominal':eigenh.tolist(),
      'r_norm':float(np.linalg.norm(r)),'models':lap,
      'finite_detailed_balance_max_log_defect_nominal':max(abs(z['residual']) for z in checks),
      'finite_detailed_balance_cases':len(checks),
      'g0_distinct_rate_per_copy':{'HB':1/12,'MET':1/12,'BAR':1/12},
      'arithmetic_status':'binary64 numerical diagnostic, not outward enclosure or exact fold certificate'}
    (base/'KI-002'/'DERIVED_SOFT_CLOCKS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'g_fold':g,'a_delta':aF-a,'b_delta':bF-b,'mobilities':{k:v['mobility_m_nominal'] for k,v in lap.items()},'log_DB_max':result['finite_detailed_balance_max_log_defect_nominal']},indent=2))
if __name__=='__main__':main()
