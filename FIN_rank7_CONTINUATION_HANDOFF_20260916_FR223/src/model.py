"""R7P-009/010/011 core coordinate and D12 machinery."""
from __future__ import annotations
import math
import numpy as np
from scipy.special import logsumexp, softmax

N=12

def strict_kernel():
    return np.array([[0.0 if i==j else
        math.cos(0.18575*min(abs(i-j),N-abs(i-j))+0.1625)/
        (1+min(abs(i-j),N-abs(i-j))**1.8)
        for j in range(N)] for i in range(N)],dtype=float)

def laplacian(W):
    return np.diag(W.sum(axis=1))-W

def sector_eigenvalues(A):
    # circulant first row convention used by the accepted audit
    return np.fft.fft(A[0]).real[:7]

def feature_spaces():
    W=strict_kernel(); A=laplacian(W); L=sector_eigenvalues(A); j=np.arange(N)
    cols=[]
    for k in (3,4,5):
        cols += [np.sqrt(L[k]/6)*np.cos(2*np.pi*k*j/N),
                 np.sqrt(L[k]/6)*np.sin(2*np.pi*k*j/N)]
    cols += [np.sqrt(L[6]/12)*(-1.0)**j]
    X7=np.column_stack(cols)
    C4=X7[:,[0,2,4,6]]
    A7=X7@X7.T
    return W,A,L,X7,C4,A7

def canonical_feature_formula(L):
    j=np.arange(N); out=[]
    for k in (3,4,5):
        out.extend([np.sqrt(L[k]/6)*np.cos(2*np.pi*k*j/N),
                    np.sqrt(L[k]/6)*np.sin(2*np.pi*k*j/N)])
    out.append(np.sqrt(L[6]/12)*(-1.0)**j)
    return np.column_stack(out)

def d12_permutation(a:int, eps:int):
    """Permutation P with (P f)[j] = f[(eps*j+a) mod 12], eps=+/-1."""
    if eps not in (-1,1): raise ValueError('eps must be +/-1')
    P=np.zeros((N,N),dtype=int)
    for j in range(N): P[j,(eps*j+a)%N]=1
    return P

def d12_actions(X7):
    gram=X7.T@X7
    inv=np.linalg.inv(gram)
    actions={}
    for eps in (1,-1):
        for a in range(N):
            P=d12_permutation(a,eps)
            T=inv@X7.T@P@X7
            actions[(a,eps)]=(P,T)
    return actions

def dual7(theta,g,X7):
    if g<=0: raise ValueError('dual formula requires g>0')
    h=X7@theta
    p=softmax(h)
    phi=float(theta@theta/(2*g)-logsumexp(h)+math.log(len(h)))
    mu=p@X7
    grad=theta/g-mu
    Y=X7-mu
    cov=Y.T@(p[:,None]*Y)
    H=np.eye(X7.shape[1])/g-cov
    return phi,grad,H,p

def primal(p,g,X7):
    p=np.asarray(p,float)
    if np.any(p<0) or not np.isclose(p.sum(),1): raise ValueError('p must lie on simplex')
    ent=float(np.sum(np.where(p>0,p*np.log(12*p),0.0)))
    q=X7.T@p
    return ent-g*float(q@q)/2
