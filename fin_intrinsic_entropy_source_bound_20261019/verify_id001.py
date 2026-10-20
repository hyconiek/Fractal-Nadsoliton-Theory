#!/usr/bin/env python3
"""ID-001: small independent controls; NEVER a mechanized proof."""
import argparse, json, math
from pathlib import Path
import numpy as np

Q=12

def get_a7():
    # Independent reconstruction of the published PHYS-002 safe_seed.py Fourier rank-7 A.
    w=np.zeros(Q)
    for d in range(1,Q):
        x=min(d,Q-d)
        w[d]=math.cos(.18575*x+.1625)/(1+x**1.8)
    W=np.array([np.roll(w,j) for j in range(Q)])
    L=np.diag(W.sum(axis=1))-W
    lams=np.fft.fft(L[0]).real
    idx=np.arange(Q)
    row_formula=np.array([sum(lams[k]/6*math.cos(2*math.pi*k*d/Q) for k in (3,4,5)) + lams[6]/12*(-1)**d for d in idx])
    cols=[]
    for k in (3,4,5):
        cols.extend([np.sqrt(lams[k]/6)*np.cos(2*np.pi*k*idx/Q),np.sqrt(lams[k]/6)*np.sin(2*np.pi*k*idx/Q)])
    cols.append(np.sqrt(lams[6]/12)*(-1.)**idx)
    X=np.column_stack(cols)
    A=X@X.T
    defect=float(np.max(np.abs(A[0]-row_formula)))
    assert defect<1e-12
    return A,lams,defect

def entropy_nats(p):
    p=np.asarray(p)
    p=p[p>0]
    return float(-np.dot(p,np.log(p)))

def h2(e):
    return entropy_nats(np.array([e,1-e]))

def fano(e):
    return h2(e)+e*math.log(Q-1)

def diagnostics(g=3.,H=1.):
    A,lams,defect=get_a7()
    x=np.exp((g/2)*(A[0]-np.max(A[0])))
    p=x/x.sum()
    pi2=float(np.dot(p,p))
    rbar=2*(1-pi2)
    cs=float(np.sum(p*(1-p)*np.log(p)))
    c=rbar*(1-math.log(2))-2*cs
    rows=[]
    for m in (32,128,512):
        delta=1/m; t=-math.expm1(-2*delta)
        P=(1-t)*np.eye(Q)+t*np.outer(np.ones(Q),p)
        assert np.max(np.abs(P.sum(axis=1)-1))<1e-12
        assert np.max(np.abs(p@P-p))<1e-12
        h=-float(np.sum(p[:,None]*P*np.log(P)))
        # Independent decomposition: diagonal no-change + off-diagonal conditional weights.
        hd=0.
        for i in range(Q):
            for j in range(Q):
                pij=1-t*(1-p[i]) if i==j else t*p[j]
                hd-=p[i]*pij*math.log(pij)
        assert abs(h-hd)<1e-12
        exact_row=-float(np.sum(p*(1-t*(1-p))*np.log(1-t*(1-p))))-t*(1-pi2)*math.log(t)-t*cs
        assert abs(h-exact_row)<1e-12
        approx=rbar*delta*math.log(1/delta)+c*delta
        residual=h-approx
        # Exact finite-delta inverse Fano lower error if h-delta H>0.
        target=max(0.,h-delta*H)
        low,high=0.,(Q-1)/Q
        if target:
            for _ in range(80):
                mid=(low+high)/2
                if fano(mid)<target:low=mid
                else:high=mid
        emin=(low+high)/2 if target else 0.
        assert h <= delta*H+fano(emin)+1e-12
        rows.append({'m':m,'delta':delta,'h_exact':h,'h_over_delta':h/delta,
                     'leading_rbar_log_plus_c':approx/delta,'residual':residual,
                     'residual_scaled_delta2_log':residual/(delta**2*abs(math.log(delta))),
                     'finite_H_example':H,'fano_lower_e':emin,'fano_lower_e_over_delta':emin/delta})
    # Independent exact intensity: stationary rate leaving i: 2(1-p_i), average rbar.
    mean_jump_rate=float(np.sum(p*2*(1-p)))
    assert abs(mean_jump_rate-rbar)<1e-14
    # Stationary law alone does not fix time evolution.
    assert np.isclose(np.dot(p, np.ones(Q)), 1.)
    # Finite sampling Bernoulli clock target uses H=h_delta/delta, diverges at refinement.
    return {'label':'NUMERICAL_DIAGNOSTIC_ONLY','n':Q,'g':g,'p':p.tolist(),
            'a7_row':A[0].tolist(),'reconstructed_a7_defect':defect,
            'positive_p_min':float(p.min()),'rbar':rbar,'c_p':c,
            'independent_mean_jump_rate':mean_jump_rate,
            'spectrum_laplacian':lams[:7].tolist(), 'rows':rows,
            'notes':['No sampling simulation or empirical validation',
                     'Entropy/joining theorem is analytical, not checked by this script',
                     'For fixed g=3 only; no unrestricted parameter scan']}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',default='results.json')
    args=parser.parse_args()
    out=diagnostics()
    Path(args.output).write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'passed':True,'rbar':out['rbar'],'c':out['c_p'],
        'min_p':out['positive_p_min'],'m':[a['m'] for a in out['rows']]},indent=2))
if __name__=='__main__':main()
