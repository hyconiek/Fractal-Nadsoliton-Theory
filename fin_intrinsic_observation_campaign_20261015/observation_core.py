"""FIN OBS: exact, bounded, *model-only* tagged-label generator and observation checks.
Reads the archived theory_core.py, which is a numerical reimplementation of the supplied
exact-decimal kernel. Floating quantities are diagnostics, never interval certificates.
"""
from __future__ import annotations
import sys, math
from pathlib import Path
import numpy as np
from scipy.special import gammaln
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'source_context'))
from theory_core import A, Q, softmax_vec, lap
SIGMA=np.array([(-1)**j for j in range(Q)], dtype=int)
CHI=float(lap[6]/12)
GC=float(1/CHI)

def label_probs(counts,old,g):
    n=np.asarray(counts,int); N=int(n.sum())
    if N<=0 or n[old]<1:raise ValueError('invalid departure')
    n=n.copy(); n[old]-=1
    return softmax_vec((g/N)*(A@n))

def parity_refresh(counts, old, g):
    return float(SIGMA@label_probs(counts,old,g))

def generator_drift_m(counts,g):
    n=np.asarray(counts,int); N=n.sum()
    return float(sum(n[i]*(parity_refresh(n,i,g)-SIGMA[i]) for i in range(Q) if n[i]>0)/N)

def generator_bracket_m(counts,g):
    n=np.asarray(counts,int); N=n.sum()
    v=0.
    for i in range(Q):
        if not n[i]:continue
        s=label_probs(n,i,g)
        v+=float(n[i]*sum(s[j]*(SIGMA[j]-SIGMA[i])**2 for j in range(Q)))
    return v/(N*N)

def labelled_stats(labels, k, g):
    lab=np.asarray(labels,int)
    N=len(lab); n=np.bincount(lab,minlength=12); k=int(k)
    if not 1<=k<=N:raise ValueError('bad k')
    signs=SIGMA[lab]
    m=float(signs.mean()); mb=float(signs[:k].mean()); v=mb-m
    c=np.array([1/k-1/N if a<k else -1/N for a in range(N)])
    rates=np.array([label_probs(n,int(i),g) for i in lab])
    q=rates@SIGMA
    h=2*(1-signs*q)
    lv=float(np.sum(c*(q-signs)))
    lm=float(np.mean(q-signs))
    lmb=float(np.mean(q[:k]-signs[:k]))
    return dict(N=N,k=k,m=m,mB=mb,v=v,Lv=lv,Lm=lm,LmB=lmb,
        forcing=lv+v, bracket_v=float(c@np.diag(h)@c),
        bracket_m=float(sum(h)/N**2),bracket_cross=float(sum(h*c)/N),
        ideal_stationary_private_bracket=2*(1/k-1/N))

def hypergeom_pmf(N,M,k,h):
    if any(type(x) is not int for x in (N,M,k,h)): raise TypeError('integer inputs')
    if not (0<=M<=N and 0<=k<=N):raise ValueError('range')
    if h<max(0,k-(N-M)) or h>min(k,M):return 0.
    return math.comb(M,h)*math.comb(N-M,k-h)/math.comb(N,k)

def hypergeom_moments(N,M,k):
    mean_plus=k*M/N
    var_h=k*(M/N)*(1-M/N)*(N-k)/(N-1) if N>1 else 0
    return dict(mean_h=mean_plus,var_h=var_h,
                mean_mB=2*M/N-1, var_mB=4*var_h/k**2)

def two_blocks_cov(N,M,k1,k2):
    if k1+k2>N:raise ValueError('overlap')
    return -4*(M/N)*(1-M/N)/(N-1)

def joint_hypergeom_pmf(N,M,k1,k2,h1,h2):
    if any(type(x) is not int for x in (N,M,k1,k2,h1,h2)):raise TypeError
    if k1+k2>N:raise ValueError
    if not (0<=h1<=k1 and 0<=h2<=k2):return 0.
    return hypergeom_pmf(N,M,k1,h1)*hypergeom_pmf(N-k1,M-h1,k2,h2)

def quartic_density_unnorm(z,chi_lam=0.):
    return np.exp(chi_lam*z*z/2-z**4/12)

if __name__=='__main__':
    print('q',Q,'chi',CHI,'gc',GC)
