#!/usr/bin/env python3
"""RG-F01 reusable exactly-positive N2 difference CTMC and classical instruments.

Physics warning: these are SUPPLIED finite stochastic laws, not physical rates.
The output likelihood is compositional for fixed emission matrices; state-changing
measurement requires an independently specified joint successor/outcome kernel.
"""
import math
import numpy as np
from scipy.linalg import expm
from scipy.special import softmax
from source_model import all_models

Q=12
KINETICS=('heat_bath','metropolis','barker')

def difference_generator(A,g=3.,kind='heat_bath'):
    """Positive 12-state quotient built without averaging phase rows.
    Reference phase (a,b)=(0,d); all other phases equal only under symmetry.
    Proposal site attempts: heat bath 1; metropolis 11/12; barker 22/12.
    """
    assert A.shape==(Q,Q)
    B=np.zeros((Q,Q))
    for d in range(Q):
        a,b=0,d
        for site in (0,1):
            old,other=(a,b) if site==0 else (b,a)
            probs=softmax((g/2)*A[:,other]) if kind=='heat_bath' else None
            for new in range(Q):
                if new==old:continue
                if kind=='heat_bath':rate=float(probs[new])
                elif kind in ('metropolis','barker'):
                    ratio=math.exp((g/2)*(A[new,other]-A[old,other]))
                    accept=min(1.,ratio) if kind=='metropolis' else ratio/(1+ratio)
                    attempt={'metropolis':11/12,'barker':22/12}[kind]
                    rate=attempt/(Q-1)*accept
                else:raise ValueError(kind)
                dest=((b-new)%Q if site==0 else (new-a)%Q)
                B[d,dest]+=rate
                B[d,d]-=rate
    return B

def stationary_difference(A,g=3.):
    return softmax(g*A[0]/2)

def propagation(B,t,kind=None,pstat=None):
    if t<0:raise ValueError('nonnegative time only')
    if kind=='heat_bath':
        if pstat is None:raise ValueError('heat bath reset requires pstat')
        return math.exp(-2*t)*np.eye(Q)+(1-math.exp(-2*t))*np.ones((Q,1))*np.asarray(pstat)[None,:]
    return expm(float(t)*B)

def emission(invalid=.0,cross_talk=.0):
    """13x12 noninvasive stochastic emission; parameters are SYNTHETIC ONLY."""
    assert 0<=invalid<=1 and 0<=cross_talk<=1-invalid
    K=np.zeros((Q+1,Q));
    for d in range(Q):
        K[d,d]=1-invalid-cross_talk
        K[(d+1)%Q,d]=cross_talk
        K[Q,d]=invalid
    return K

def sequence_probability(p0,B,times,outcomes,K,kind=None,pstat=None):
    """Probability of a particular observation record; all outcomes kept.
    Noninvasive emission is diagonal in d; no normalized intermediate filter.
    """
    p=np.array(p0,dtype=float).copy()
    assert len(times)==len(outcomes)
    last=0.
    for t,o in zip(times,outcomes):
        assert t>=last
        p=p@propagation(B,t-last,kind,pstat)
        p*=K[int(o)]
        last=t
    return float(p.sum())

def full_joint(p0,B,times,K,kind=None,pstat=None):
    assert len(times)==2 and times[1]>=times[0]>=0
    P1=propagation(B,times[0],kind,pstat)
    P2=propagation(B,times[1]-times[0],kind,pstat)
    p1=p0@P1
    # unnormalized instrument branch, then propagate, then second emission
    return np.array([( (p1*K[i])@P2)@K.T for i in range(K.shape[0])])

def uniformization(B,t,rate=4.,m=90):
    """Positive stochastic polynomial. Poisson truncation loses only tail mass;
    no renormalization. Caller must separately bound input and floating errors.
    """
    assert rate>=max(-np.diag(B)) and t>=0
    P=np.eye(Q)+B/rate
    z=rate*t
    v=np.eye(Q); weight=math.exp(-z); T=weight*v.copy()
    for j in range(1,m+1):
        v=v@P
        weight*=z/j
        T+=weight*v
    return T

if __name__=='__main__':
    A=all_models()['FIN_A7'];B=difference_generator(A)
    print(np.max(np.abs(B-(2*(np.ones((Q,1))*stationary_difference(A)[None,:]-np.eye(Q))))))
