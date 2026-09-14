"""R7P-012 stable dual derivatives through fourth order."""
from __future__ import annotations
import numpy as np
from scipy.special import softmax, logsumexp
import math

def moments(features,theta):
    h=features@theta; p=softmax(h); mu=p@features; Y=features-mu
    cov=Y.T@(p[:,None]*Y)
    return p,mu,Y,cov

def dual_all(theta,g,features):
    if g<=0: raise ValueError('g must be positive')
    p,mu,Y,cov=moments(features,theta)
    val=float(theta@theta/(2*g)-logsumexp(features@theta)+math.log(len(p)))
    grad=theta/g-mu
    H=np.eye(len(theta))/g-cov
    T3=-np.einsum('i,ia,ib,ic->abc',p,Y,Y,Y)
    raw4=np.einsum('i,ia,ib,ic,id->abcd',p,Y,Y,Y,Y)
    pair=(np.einsum('ab,cd->abcd',cov,cov)+
          np.einsum('ac,bd->abcd',cov,cov)+
          np.einsum('ad,bc->abcd',cov,cov))
    T4=-(raw4-pair)
    return val,grad,H,T3,T4,p

def directional_34(theta,g,features,v):
    *_,T3,T4,p=dual_all(theta,g,features)
    d3=float(np.einsum('abc,a,b,c',T3,v,v,v))
    d4=float(np.einsum('abcd,a,b,c,d',T4,v,v,v,v))
    return d3,d4
