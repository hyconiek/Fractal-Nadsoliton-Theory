#!/usr/bin/env python3
"""Bounded KO spectral, conditional hypergeometric and Gaussian-event checks.
No Monte Carlo, no imports of FIN theory_core.py. Exact rational projector
identities and finite hypergeometric enumeration; interval eigen calculations.
"""
from __future__ import annotations
from fractions import Fraction as F
from itertools import combinations
from math import cos, pi, sqrt, factorial
from pathlib import Path
import json
import mpmath as mp
import numpy as np

Q = 12
A0 = F(0)

def zero(): return [[F(0) for _ in range(Q)] for _ in range(Q)]

def mm(A,B): return [[sum(A[i][k]*B[k][j] for k in range(Q)) for j in range(Q)] for i in range(Q)]

def add(A,B): return [[A[i][j]+B[i][j] for j in range(Q)] for i in range(Q)]

def cos3(d): return (F(1), F(0), F(-1), F(0))[d%4]

def cos4(d): return (F(1), F(-1,2), F(-1,2))[d%3]

def proj(cf): return [[F(1,6)*cf(i-j) for j in range(Q)] for i in range(Q)]

def rational_rank(matrix):
    M=[list(x) for x in matrix]; m=len(M); n=len(M[0]);r=0
    for j in range(n):
        pivot=next((i for i in range(r,m) if M[i][j]),None)
        if pivot is None: continue
        M[r],M[pivot]=M[pivot],M[r]
        d=M[r][j];M[r]=[x/d for x in M[r]]
        for i in range(m):
            if i!=r:
                d=M[i][j]
                if d:M[i]=[v-d*w for v,w in zip(M[i],M[r])]
        r+=1
        if r==m:break
    return r

def is_same(A,B):return all(a==b for u,v in zip(A,B) for a,b in zip(u,v))

mp.mp.dps=80;mp.iv.dps=65
iv=mp.iv
OMEGA=iv.mpf('0.18575');PHI=iv.mpf('0.16250');ETA=iv.mpf(9)/5
weights=[iv.mpf(0)] + [iv.cos(OMEGA*d+PHI)/(1+iv.mpf(d)**ETA) for d in range(1,7)]
def eig_iv(k):
    return 2*sum((weights[d]*(1-iv.cos(iv.pi*k*d/6)) for d in range(1,6)),iv.mpf(0)) +weights[6]*(1-iv.cos(iv.pi*k))

def interval_bounds_str(z):
    # str of mpi endpoints yields decimal enclosing values; expose endpoint objects separately.
    return {'lower':str(z.a),'upper':str(z.b)}

def ivcmp_gt(a,b):
    a=iv.mpf(a);b=iv.mpf(b)
    return float(a.a)>float(b.b)

def erf_iv_series(x,terms=35):
    # Alternating convergent series at x<1, enclosing tail by first omitted term.
    t=x
    s=iv.mpf(0)
    for n in range(terms):
        if n: t=t*x*x/n
        s+= ((-1)**n)*t/(2*n+1)
    tnext=t*x*x/terms
    err=abs(tnext)/(2*terms+1)
    erf=2/iv.sqrt(iv.pi)*s
    errb=2/iv.sqrt(iv.pi)*err
    return erf + iv.mpf([-1,1])*errb

def event_iv(gamma):
    variance=iv.mpf(1)+1/gamma
    x=1/iv.sqrt(2*variance)
    return 1-erf_iv_series(x)

def independently_reconstruct_A(lap):
    p={}
    for k in (3,4,5):
        p[k]=np.array([[2/12*cos(2*pi*k*(i-j)/12) for j in range(Q)] for i in range(Q)])
    sig=np.array([(-1)**i for i in range(Q)])
    p[6]=np.outer(sig,sig)/12
    return sum(float(lap[k])*p[k] for k in (3,4,5,6))

def run():
    assertions=0
    def ck(c,name):
        nonlocal assertions
        assert bool(c),name
        assertions+=1
    P3,P4=proj(cos3),proj(cos4)
    P34=[[P3[i][j]-P4[i][j] for j in range(Q)] for i in range(Q)]
    ck(is_same(mm(P3,P3),P3),'P3^2=P3')
    ck(is_same(mm(P4,P4),P4),'P4^2=P4')
    ck(is_same(mm(P3,P4),zero()),'P3P4=0')
    ck(is_same(mm(P4,P3),zero()),'P4P3=0')
    ck(rational_rank(P3)==2,'P3 rank2')
    ck(rational_rank(P4)==2,'P4 rank2')
    ck(rational_rank(P34)==4,'P3-P4 rank4')
    ck(all(P34[i][i]==0 for i in range(Q)),'diagonal fixed')
    ck(sum(P34[i][i] for i in range(Q))==0,'trace fixed')
    ck(all(sum(P34[i])==0 for i in range(Q)),'centered')
    sg=[F((-1)**i) for i in range(Q)]
    ck(all(sum(P34[i][j]*sg[j] for j in range(Q))==0 for i in range(Q)),'parity fixed')
    ck(all(P34[i][j]==P34[(i+1)%Q][(j+1)%Q] for i in range(Q) for j in range(Q)),'circulant')
    ck(all(P34[i][j]==P34[j][i] for i in range(Q) for j in range(Q)),'symmetric')
    lap={k:eig_iv(k) for k in (3,4,5,6)}
    lo={k:float(lap[k].a) for k in lap};hi={k:float(lap[k].b) for k in lap}
    for k in lap:ck(ivcmp_gt(lap[k],0),f'lambda{k}>0')
    ck(all(ivcmp_gt(lap[6],lap[k]) for k in (3,4,5)),'lambda6 unique top')
    ck(ivcmp_gt(lap[6],lap[5]+iv.mpf('0.01')),'gap remains for epsilon .01')
    gain=12/lap[6]
    ck(ivcmp_gt(gain,iv.mpf(5)) and ivcmp_gt(iv.mpf(6),gain),'5<gc<6')
    gammas={k:1-lap[k]/lap[6] for k in (3,4,5)}
    ck(all(ivcmp_gt(x,0) for x in gammas.values()),'positive gap')
    gammaplus=1-(lap[3]+iv.mpf('0.01'))/lap[6]
    ck(ivcmp_gt(gammas[3],gammaplus) and ivcmp_gt(gammaplus,0),'gamma3 decrease')
    ck(ivcmp_gt(gammaplus,gammas[5]),'max fast stays mode5')
    p0=event_iv(gammas[3]);pe=event_iv(gammaplus)
    ck(ivcmp_gt(pe,p0),'strict bounded event contrast')
    diff=pe-p0
    risk=(1-diff)/2
    ck(ivcmp_gt(risk,0) and ivcmp_gt(iv.mpf('0.5'),risk),'strict ideal advantage')
    ck(float(pe.a)>.7 and float(pe.b)<.8,'event range')
    # Exact uniform population n=(1,...,1), N12, k6; multivariate hypergeometric all 924 combinations equiprobable.
    N=12;k=6
    ck(len(list(combinations(range(N),k)))==924,'sample count')
    # enumerate with ordinary float; exact combinatorial probabilities = 1/924.
    features=np.array([[sqrt(2)*cos(2*pi*r*j/Q) for r in (1,3,4)] for j in range(Q)])
    mu=np.mean(features,axis=0)
    vals=[]
    for subset in combinations(range(N),k):
        vals.append(sqrt(N)*(np.mean(features[list(subset)],axis=0)-mu))
    vals=np.asarray(vals)
    ck(float(np.abs(vals.mean(axis=0)).max())<1e-12,'zero sample mean')
    cv=(vals.T@vals)/len(vals)
    predicted=np.eye(3)*N*(N-k)/(k*(N-1))
    ck(float(np.max(np.abs(cv-predicted)))<1e-11,'hypergeom finite correction covariance')
    # Independent exact rational enumeration of cos r=3,4 products (sqrt(2) factor cancels into 2).
    fr=[[cos3(j),cos4(j)] for j in range(12)]
    exact_cov=[[F(0) for _ in range(2)] for _ in range(2)]
    for subset in combinations(range(12),6):
        sums=[sum(fr[j][r] for j in subset) for r in range(2)]
        for r in range(2):
            for s in range(2):
                exact_cov[r][s]+=F(2*12,6*6)*sums[r]*sums[s]/924
    ck(exact_cov[0][0]==F(12,11),'exact q12 hypergeom r3 variance')
    ck(exact_cov[1][1]==F(12,11),'exact q12 hypergeom r4 variance')
    ck(exact_cov[0][1]==0 and exact_cov[1][0]==0,'exact q12 hypergeom cross covariance')
    # independent binary64 reconstruction, checksum against the original source's 7 columns made independently.
    nomlap={i:float(lap[i].mid) for i in lap}
    A=independently_reconstruct_A(nomlap)
    ck(np.max(np.abs(A-A.T))<1e-13,'symmetric independent A')
    ck(np.max(abs(A@np.ones(Q)))<1e-12,'center independent A')
    ck(np.max(abs(A@np.array([(-1)**i for i in range(Q)])-nomlap[6]*np.array([(-1)**i for i in range(Q)])))<1e-12,'parity eigen')
    P34n=np.array([[float(a) for a in row] for row in P34]); Aplus=A+.01*P34n
    ck(np.max(abs(np.diag(Aplus)-np.diag(A)))<1e-14,'constant diagonal')
    ck(np.linalg.matrix_rank(Aplus,tol=1e-10)==7,'rank7')
    ck(np.max(abs((Aplus-A)@np.array([(-1)**i for i in range(Q)])))<1e-12,'parity eigen unchanged')
    data={'scientific_status':'CONDITIONAL_ANALYTICAL_STATIC_KO_COMPARISON','assertions_passed':assertions,
      'interval_method':'mpmath.iv dps65; 35-term alternating erf series with bounded first omitted term; no external theorem mechanization',
      'source_parameters':{'omega':'0.18575','phi':'0.16250','eta':'1.8','q':12,'epsilon':'0 or +0.01','block_fraction':'1/2','event':'abs(V_3,N)>1'},
      'spectral_intervals':{str(k):interval_bounds_str(v) for k,v in lap.items()},
      'g_c_interval':interval_bounds_str(gain),
      'relative_fast_gap_lower_comparator':float(min(float(x.a) for x in [gammas[5],gammas[4]+iv.mpf('0.01'),gammaplus])),
      'gamma3_A7':interval_bounds_str(gammas[3]),'gamma3_comparator':interval_bounds_str(gammaplus),
      'probability_A7_interval':interval_bounds_str(p0),'probability_comparator_interval':interval_bounds_str(pe),
      'event_probability_gap_interval':interval_bounds_str(diff),'ideal_equal_prior_risk_interval':interval_bounds_str(risk),
      'finite_hypergeom_N12_k6':{'configurations':924,'covariance':cv.tolist(),'theoretical_diagonal':12/11,'max_residual':float(np.max(abs(cv-predicted)))},
      'nominal_A7_diag':float(A[0,0]),'nominal_A7_Mmax':float(np.max(abs(A))), 'notes':[
         'Exact projector statements use Fraction; eigendata use interval arithmetic.',
         'Finite hypergeometric check verifies one small combinatorial example, not its growing-N CLT.',
         'Static comparator/local-window theorem uses inherited MT/UF analytical premise and local extension proof.',
         'Event probability does not supply finite-N rate or calibrated physical sample size.']}
    out=Path(__file__).with_name('INDEPENDENT_VERIFICATION.json');out.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
    print('PASSED',assertions,'checks')
    print('lambda', {k:nomlap[k] for k in nomlap})
    print('g',float(gain.mid));print('p0',float(p0.mid),'pplus',float(pe.mid),'gap',float(diff.mid),'risk',float(risk.mid))
if __name__=='__main__':run()
