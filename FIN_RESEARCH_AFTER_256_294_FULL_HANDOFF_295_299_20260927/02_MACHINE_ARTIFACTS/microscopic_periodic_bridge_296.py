#!/usr/bin/env python3
"""FIN report 296: microscopic periodic bridge replay.

Builds the exact labelled finite-N leave-one-out Gibbs heat-bath skeleton for N=2,3,
its exact occupation-count quotient for N=2..7, checks stationarity/detailed balance,
computes slow spectral gaps, and evaluates finite cyclic-bridge cylinder TV errors
for the full labelled N=2 chain.

Numerics are double precision; the analytic bridge/convergence theorem is in the report.
"""
import itertools, json, math
import numpy as np
from scipy.special import softmax, gammaln
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import eigsh

NLAB=12
LAM={3:1.96140686197644,4:2.19956884933321,5:2.2986062720790903,6:2.3421820411463}
j=np.arange(NLAB); cc=[]
for k in (3,4,5):
    cc += [np.sqrt(LAM[k]/6)*np.cos(2*np.pi*k*j/NLAB),
           np.sqrt(LAM[k]/6)*np.sin(2*np.pi*k*j/NLAB)]
cc += [np.sqrt(LAM[6]/12)*(-1.0)**j]
X=np.column_stack(cc); A=X@X.T
G=5.145228719489144

def compositions(N,k=12):
    for bars in itertools.combinations(range(N+k-1),k-1):
        prev=-1; a=[]
        for b in bars+(N+k-1,):
            a.append(b-prev-1); prev=b
        yield tuple(a)

def labelled_chain(N):
    states=list(itertools.product(range(12),repeat=N)); idx={s:i for i,s in enumerate(states)}
    rows=[]; cols=[]; vals=[]; logw=np.empty(len(states)); fresh=0.0
    for u,x in enumerate(states):
        arr=np.array(x,int)
        logw[u]=(G/(2*N))*float(A[np.ix_(arr,arr)].sum())-N*math.log(12)
        acc={}
        for a in range(N):
            m=np.bincount(np.delete(arr,a),minlength=12)
            q=softmax((G/N)*(A@m))
            for z in range(12):
                y=list(x); y[a]=z; v=idx[tuple(y)]
                acc[v]=acc.get(v,0.0)+(1/N)*q[z]
        # P=I+Q/(2N)=(I+K)/2; K is one chosen-label Gibbs refresh event.
        acc[u]=acc.get(u,0.0)+1.0
        for v,p in acc.items(): rows.append(u); cols.append(v); vals.append(0.5*p)
    P=coo_matrix((vals,(rows,cols)),shape=(len(states),len(states))).tocsr()
    z=logw-logw.max(); pi=np.exp(z); pi/=pi.sum()
    # Average H(new symbol | current full microstate, chosen label).
    for w,x in zip(pi,states):
        arr=np.array(x,int)
        for a in range(N):
            m=np.bincount(np.delete(arr,a),minlength=12)
            q=softmax((G/N)*(A@m))
            fresh += w*(1/N)*float(-np.sum(q*np.log2(q)))
    return states,P,pi,fresh

def count_chain(N):
    states=list(compositions(N)); idx={s:i for i,s in enumerate(states)}
    rows=[]; cols=[]; vals=[]; logw=np.empty(len(states))
    for u,nt in enumerate(states):
        n=np.array(nt,int)
        logw[u]=gammaln(N+1)-np.sum(gammaln(n+1))-N*math.log(12)+(G/(2*N))*float(n@A@n)
        acc={}
        for i in np.flatnonzero(n):
            m=n.copy(); m[i]-=1
            q=softmax((G/N)*(A@m)); wi=n[i]/N
            for z in range(12):
                nn=m.copy(); nn[z]+=1; v=idx[tuple(nn)]
                acc[v]=acc.get(v,0.0)+wi*q[z]
        acc[u]=acc.get(u,0.0)+1.0
        for v,p in acc.items(): rows.append(u); cols.append(v); vals.append(0.5*p)
    P=coo_matrix((vals,(rows,cols)),shape=(len(states),len(states))).tocsr()
    z=logw-logw.max(); pi=np.exp(z); pi/=pi.sum()
    return states,P,pi

def spectral(P,pi,N,k=6):
    s=np.sqrt(pi); S=diags(s)@P@diags(1/s)
    vals=eigsh(S,k=min(k,P.shape[0]-1),which='LA',return_eigenvectors=False,tol=2e-10,maxiter=200000)
    vals=np.sort(vals)[::-1]; gp=float(1-vals[1]); gq=float(2*N*gp)
    db=(P.multiply(pi[:,None])-P.T.multiply(pi[None,:])).tocoo()
    dbmax=float(np.max(np.abs(db.data))) if db.nnz else 0.0
    return dict(top_eigenvalues=[float(x) for x in vals],gap_P=gp,tau_ticks=1/gp,
                gap_Q=gq,tau_continuous=1/gq,detailed_balance_max=dbmax,
                row_sum_max_error=float(np.max(np.abs(np.asarray(P.sum(axis=1)).ravel()-1))))

def direct_cylinder_tv(P,pi):
    D=P.toarray(); out={}
    for r in (0,1,2,3):
        Pr=np.linalg.matrix_power(D,r); rows=[]
        for L in (4,8,16,32,64,128,256,512):
            if L<=r: continue
            Pm=np.linalg.matrix_power(D,L-r); PL=np.linalg.matrix_power(D,L)
            Z=float(np.trace(PL))
            ratio=Pm.T/(pi[:,None]*Z)
            joint=pi[:,None]*Pr
            tv=0.5*float(np.sum(joint*np.abs(ratio-1)))
            rows.append(dict(L=L,closing_steps=L-r,TV=tv,Z=Z))
        out[str(r)]=rows
    return out

out={'g':G,'uniformization_rate':'2N','analytic_g0_gap_Q':1.0,'labelled':{},'count_quotient':{}}
for N in (2,3):
    st,P,pi,h=labelled_chain(N); rec=spectral(P,pi,N); rec['states']=len(st); rec['fresh_entropy_bits_per_conditioned_update']=h
    out['labelled'][str(N)]=rec
    if N==2: out['labelled_N2_cylinder_TV']=direct_cylinder_tv(P,pi)
for N in (2,3,4,5,6,7):
    st,P,pi=count_chain(N); rec=spectral(P,pi,N,k=4); rec['states']=len(st); rec['pi_min']=float(pi.min())
    out['count_quotient'][str(N)]=rec
# slow gap is in count sector for N=2,3 within numerical tolerance
for N in (2,3):
    assert abs(out['labelled'][str(N)]['gap_P']-out['count_quotient'][str(N)]['gap_P'])<2e-11
assert max(out['labelled'][str(N)]['detailed_balance_max'] for N in (2,3))<1e-12
assert max(out['count_quotient'][str(N)]['detailed_balance_max'] for N in (2,3,4,5,6,7))<1e-12
json.dump(out,open('/mnt/data/fin296/MICROSCOPIC_PERIODIC_BRIDGE_296.json','w'),indent=2)
print(json.dumps(out,indent=2))
