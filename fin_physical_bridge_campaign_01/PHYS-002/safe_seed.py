#!/usr/bin/env python3
"""Small-N, fail-closed FIN bridge audit module.

This is an independent reconstruction from versioned finite-N formulas used in
reports 61/138 and the A7 construction present at repository commit
97c4231f33632800fd817fe2294555bb8bcb041f. It is NOT a silent replacement
for the untracked historical `FIN son` files.
"""
from __future__ import annotations
import itertools, math
import numpy as np
import scipy.sparse as sp
from scipy.special import softmax, logsumexp

Q = 12
G_FROZEN = 5.145228719489142
OMEGA = 0.18575
PHI = 0.1625
ETA = 1.8


def build_rank7():
    j=np.arange(Q)
    W=np.zeros((Q,Q),float)
    for i in range(Q):
        for k in range(Q):
            if i==k: continue
            d=min(abs(i-k),Q-abs(i-k))
            W[i,k]=math.cos(OMEGA*d+PHI)/(1+d**ETA)
    L=np.diag(W.sum(1))-W
    lap_lam=np.fft.fft(L[0]).real[:7]
    cols=[]
    for k in (3,4,5):
        cols += [np.sqrt(lap_lam[k]/6)*np.cos(2*np.pi*k*j/Q),
                 np.sqrt(lap_lam[k]/6)*np.sin(2*np.pi*k*j/Q)]
    cols += [np.sqrt(lap_lam[6]/12)*(-1.)**j]
    X=np.column_stack(cols)
    A=X@X.T
    return W,L,lap_lam,X,A


def compositions(N:int,q:int=Q):
    out=[]
    for bars in itertools.combinations(range(N+q-1),q-1):
        prev=-1; vals=[]
        for b in bars+(N+q-1,):
            vals.append(b-prev-1); prev=b
        out.append(tuple(vals))
    return out


def count_logweights(states,g,A,theta=0.0):
    N=sum(states[0]); out=[]
    for s in states:
        n=np.asarray(s,float)
        logmult=math.lgamma(N+1)-sum(math.lgamma(float(x)+1) for x in n)
        out.append(logmult+(g/(2*N))*float(n@A@n)+theta*n[0])
    return np.array(out)


def count_pi(states,g,A,theta=0.0):
    lw=count_logweights(states,g,A,theta); lw-=logsumexp(lw)
    return np.exp(lw)


def count_generator(states,g,A,kinetic='heat_bath'):
    """Continuous-time count generator; each labelled copy has attempt rate 1.

    Supported kinetics share the same Gibbs stationary law for symmetric,
    constant-diagonal A. Metropolis/Barker use uniform 12-state proposals.
    """
    N=sum(states[0]); M=len(states); idx={s:i for i,s in enumerate(states)}
    rows=[]; cols=[]; data=[]; diag=np.zeros(M)
    for ii,s in enumerate(states):
        n=np.array(s,int)
        for old in np.flatnonzero(n):
            m=n.copy(); m[old]-=1
            field=(g/N)*(A@m)
            for new in range(Q):
                if new==old: continue
                if kinetic=='heat_bath':
                    rate=n[old]*softmax(field)[new]
                else:
                    logratio=float(field[new]-field[old])
                    if kinetic=='metropolis':
                        acc=min(1.0, math.exp(logratio))
                    elif kinetic=='barker':
                        acc=1/(1+math.exp(-logratio))
                    else:
                        raise ValueError('unknown kinetic')
                    rate=n[old]*(1/Q)*acc
                n2=n.copy(); n2[old]-=1; n2[new]+=1
                jj=idx[tuple(map(int,n2))]
                rows.append(ii); cols.append(jj); data.append(rate); diag[ii]-=rate
    rows.extend(range(M)); cols.extend(range(M)); data.extend(diag.tolist())
    return sp.csr_matrix((data,(rows,cols)),shape=(M,M))


def generator_checks(states,Qmat,pi,fail_tol=1e-10):
    A=Qmat.toarray()
    row=float(np.max(np.abs(A.sum(1))))
    off=A.copy(); np.fill_diagonal(off,0)
    minoff=float(off.min())
    station=float(np.max(np.abs(pi@A)))
    db=float(np.max(np.abs(pi[:,None]*A-(pi[:,None]*A).T)))
    ok=(row<fail_tol and minoff>=-fail_tol and station<fail_tol and db<fail_tol)
    return {'row_sum_max_abs':row,'min_offdiag':minoff,'stationarity_max_abs':station,
            'detailed_balance_max_abs':db,'passed':bool(ok)}


def rotation_permutation(states,shift=1):
    idx={s:i for i,s in enumerate(states)}
    return np.array([idx[tuple(np.roll(np.asarray(s),shift))] for s in states],int)


def sector_basis(states,k):
    """Orthonormal C12 Fourier-orbit basis; robust to degenerate generator eigenvalues."""
    perm=rotation_permutation(states,1); seen=set(); vec=[]
    z=np.exp(2j*np.pi*k/Q); M=len(states)
    for start in range(M):
        if start in seen: continue
        orb=[]; cur=start
        while cur not in orb:
            orb.append(cur); seen.add(cur); cur=int(perm[cur])
        L=len(orb)
        if abs(z**L-1)>1e-8: continue
        v=np.zeros(M,complex)
        for r,ix in enumerate(orb): v[ix]=z**(-r)/math.sqrt(L)
        vec.append(v)
    return np.column_stack(vec) if vec else np.zeros((M,0),complex)


def sector_spectra(states,Qmat,pi):
    q=Qmat.toarray(); d=np.sqrt(pi)
    S=(d[:,None]*q)/d[None,:]
    rev_def=float(np.max(np.abs(S-S.T)))
    out={}
    for k in range(Q):
        B=sector_basis(states,k)
        H=B.conj().T@S@B
        herm=float(np.max(np.abs(H-H.conj().T))) if H.size else 0.0
        vals,vecs=np.linalg.eigh((H+H.conj().T)/2)
        maxres=0.0
        for val,v in zip(vals,vecs.T):
            f=B@v; maxres=max(maxres,float(np.linalg.norm(S@f-val*f)))
        out[str(k)]={'dimension':int(B.shape[1]),'eigenvalues':[float(x) for x in vals],
                     'max_eigenpair_residual':maxres,'hermiticity_defect':herm}
    return rev_def,out


def static_S(states,pi,N,k):
    phase=np.exp(2j*np.pi*k*np.arange(Q)/Q)
    vals=[]
    for s in states:
        z=np.dot(np.asarray(s,float),phase)
        vals.append((abs(z)**2)/N)
    return float(np.dot(pi,np.asarray(vals)))


def pair_difference_distribution(A,g):
    # N=2 conditional j|i=0; constant diagonal makes this also pair-difference law.
    return softmax((g/2)*A[0])


def build_countermodels(A):
    tr=float(np.trace(A)); lam=np.fft.fft(A[0]).real
    def from_lam(lv):
        row=np.fft.ifft(lv).real
        return np.array([np.roll(row,i) for i in range(Q)])
    P0=np.eye(Q)-np.ones((Q,Q))/Q
    models={'FIN_A7':A.copy(),'FULL_POTTS_TRACE':P0*(tr/11)}
    lv=np.zeros(Q); lbar=tr/7
    for k in (3,4,5,7,8,9): lv[k]=lbar
    lv[6]=lbar; models['FLAT_P7_TRACE']=from_lam(lv)
    b={3:lam[3],4:lam[4],5:lam[5],6:lam[6]}
    specs=[]
    s=.1*min(b[3],b[4]); specs.append(('PERT_34_10P',{3:s,4:-s}))
    s=.1*min(b[3],b[5]); specs.append(('PERT_35_10P',{3:s,5:-s}))
    s=.05*min(b[3],b[6]/2); specs.append(('PERT_36_5P',{3:s,6:-2*s}))
    for name,dd in specs:
        x=lam.copy()
        for k,dv in dd.items():
            x[k]+=dv
            if k!=6: x[Q-k]+=dv
        models[name]=from_lam(x)
    f=.02; x=lam.copy()*(1-f); leak=f*tr/4
    for k in (1,2,10,11): x[k]=leak
    models['LEAK_K12_2P_TRACE']=from_lam(x)
    return models

if __name__ == '__main__':
    # Deliberately small, explicit smoke test; never launches large N on import.
    import argparse, json
    ap=argparse.ArgumentParser(); ap.add_argument('--N',type=int,default=2,choices=[1,2,3,4]); ap.add_argument('--g',type=float,default=G_FROZEN)
    args=ap.parse_args()
    _,_,_,_,A=build_rank7(); st=compositions(args.N); pi=count_pi(st,args.g,A); Qm=count_generator(st,args.g,A)
    print(json.dumps({'N':args.N,'g':args.g,'states':len(st),'checks':generator_checks(st,Qm,pi)},indent=2))
