#!/usr/bin/env python3
"""Independent local mathematical reimplementation of public safe_seed definitions.
FIN formulas originate in PHYS-002/safe_seed.py; source bytes not available here.
This is supplied-model arithmetic, not a physical data generator.
"""
import math
import numpy as np
from scipy.special import softmax
Q=12
G_M3=5.145228719489142
G_N2=3.0

def fin_matrix():
    w=np.zeros(Q)
    for d in range(1,Q):
        r=min(d,Q-d)
        w[d]=math.cos(0.18575*r+0.1625)/(1+r**1.8)
    lam=np.array([sum(w[d]*(1-math.cos(2*math.pi*k*d/Q)) for d in range(1,Q)) for k in range(Q)])
    row=np.array([sum(lam[k]/6*math.cos(2*math.pi*k*d/Q) for k in (3,4,5))+lam[6]/12*(-1)**d for d in range(Q)])
    return np.array([np.roll(row,j) for j in range(Q)])

def all_models():
    A=fin_matrix();tr=np.trace(A);q=Q
    potts=(np.eye(q)-np.ones((q,q))/q)*(tr/(q-1))
    lam=np.zeros(q);lam[[3,4,5,6,7,8,9]]=tr/7
    flat=np.array([np.roll(np.fft.ifft(lam).real,j) for j in range(q)])
    full_lam=np.fft.fft(A[0]).real
    perturb=full_lam.copy();v=.1*min(full_lam[3],full_lam[4]);perturb[[3,9]]+=v;perturb[[4,8]]-=v
    pert=np.array([np.roll(np.fft.ifft(perturb).real,j) for j in range(q)])
    return {'FIN_A7':A, 'FULL_POTTS_TRACE':potts,'FLAT_P7_TRACE':flat,'PERT_34_10P':pert}

def pair_difference(A,g=G_N2):
    return softmax((g/2)*A[0])

def m3_visible(A,g=G_M3):
    i,j,k=np.indices((Q,Q,Q))
    # n^T A n = 3*diag + 2(Aij+Aik+Ajk), const diag
    F=(g/3)*(A[i,j]+A[i,k]+A[j,k])
    from scipy.special import logsumexp
    p=np.exp(F-logsumexp(F))
    return p.sum(axis=2)

def conditional_R(A,b,c,g=G_M3):
    return softmax((g/3)*(A[:,b]+A[:,c]))

def tv(p,q):
    return float(np.sum(np.abs(np.asarray(p)-np.asarray(q)))/2)

def build_n2_generator(A,g=G_N2,kind='heat_bath'):
    """L rows x->x' in ordered 144 states; fixed g0 nonzero-mode relaxation=1.
    Metropolis/Barker propose uniformly among 11 other labels; attempts 11/12 or 22/12.
    All three therefore coincide with heat bath at g=0.
    """
    from scipy.sparse import csr_matrix
    rows=[];cols=[];data=[]
    proposal_attempt={'heat_bath':1., 'metropolis':11/12, 'barker':22/12}[kind]
    diag=np.zeros(Q*Q)
    for a in range(Q):
      for b in range(Q):
        ix=Q*a+b
        for site in [0,1]:
          old,other=(a,b) if site==0 else (b,a)
          probabilities=softmax((g/2)*A[:,other]) if kind=='heat_bath' else None
          for new in range(Q):
            if new==old: continue
            if kind=='heat_bath': rate=probabilities[new]
            else:
              r=math.exp((g/2)*(A[new,other]-A[old,other]))
              acc=min(1.,r) if kind=='metropolis' else r/(1+r)
              rate=proposal_attempt/(Q-1)*acc
            jj=Q*new+b if site==0 else Q*a+new
            rows.append(ix);cols.append(jj);data.append(rate);diag[ix]-=rate
    rows.extend(range(Q*Q));cols.extend(range(Q*Q));data.extend(diag.tolist())
    return csr_matrix((data,(rows,cols)),shape=(Q*Q,Q*Q))

def pi_n2(A,g=G_N2):
    from scipy.special import logsumexp
    x=np.arange(Q)
    f=(g/2)*A[np.ix_(x,x)]
    return np.exp(f-logsumexp(f)).ravel()

def static_mode_shape(A):
    # centered A row log-odds ratios independent of g; repeated eigenpairs (k,12-k)
    spectrum=np.fft.fft(A[0]).real
    active=[int(k) for k in range(1,7) if spectrum[k]>1e-7]
    return {'positive_modes_k1_to_k6':active,'lambdas_k1_to_k6':[float(v) for v in spectrum[1:7]],
            'active_ratios_to_first':[float(spectrum[k]/spectrum[active[0]]) for k in active] if active else [],
            'trace':float(np.trace(A))}
