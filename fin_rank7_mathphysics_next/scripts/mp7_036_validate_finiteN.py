#!/usr/bin/env python3
import os
import itertools, json, math, pathlib, sys
import numpy as np
from scipy.special import logsumexp
ROOT=pathlib.Path(os.environ.get('R7P_ROOT','/mnt/data/r7p_source/unpacked/fin_rank7_followup'))
sys.path.insert(0,str(ROOT))
from src.model import feature_spaces

W,A0,L,X,C4,A=feature_spaces()
# all nonnegative compositions of N into k bins
def comps(N,k,prefix=()):
    if k==1:
        yield prefix+(N,); return
    for a in range(N+1): yield from comps(N-a,k-1,prefix+(a,))

def log_multinomial(n):
    N=sum(n); return math.lgamma(N+1)-sum(math.lgamma(x+1) for x in n)

def occ_terms(N,g,f=None):
    if f is None: f=np.zeros(7)
    out=[]
    for n in comps(N,12):
        n=np.array(n,dtype=float); p=n/N; mu=X.T@p
        logw=log_multinomial(n)-N*math.log(12)+(N*g/2)*float(p@A@p)+N*float(f@mu)
        out.append((logw,mu,n))
    return out

def Z_occ(N,g): return math.exp(logsumexp([z for z,_,__ in occ_terms(N,g)]))

def Z_labels(N,g):
    vals=[]
    for xs in itertools.product(range(12), repeat=N):
        S=X[list(xs)].sum(axis=0)
        vals.append(-N*math.log(12)+g*float(S@S)/(2*N))
    return math.exp(logsumexp(vals))

def mean_cov_mu(N,g,f):
    ts=occ_terms(N,g,f); logs=np.array([z for z,_,__ in ts]); ws=np.exp(logs-logsumexp(logs)); mus=np.array([mu for _,mu,__ in ts])
    mean=ws@mus; Y=mus-mean
    cov=np.einsum('i,ia,ib->ab',ws,Y,Y)
    return mean,cov

def d12_perm_counts(n,a,eps):
    # new occupation q_j = n_{eps*j+a}; preserves weight by circulant/reflection invariance
    q=np.zeros(12,int)
    for j in range(12): q[j]=n[(eps*j+a)%12]
    return q

out={'task':'MP7-036','scientific_state':'NUMERICAL_EVIDENCE_VALIDATING_ANALYTIC_IDENTITIES','backend':{'python':sys.version,'numpy':np.__version__},'tests':{}}
out['tests']['feature_structure']={
 'X_column_sum_max_abs':float(np.abs(X.sum(axis=0)).max()),
 'A_diag_spread':float(np.ptp(np.diag(A))),
 'A_symmetry_max_abs':float(np.abs(A-A.T).max())}
# label/occupation exact combinatorial representation checks
small=[]
for N in [1,2,3,4]:
    for g in [0.2,1.1]:
        zo=Z_occ(N,g); zl=Z_labels(N,g)
        small.append({'N':N,'g':g,'Z_labels':zl,'Z_occupations':zo,'abs_diff':abs(zl-zo),'rel_diff':abs(zl-zo)/zo})
out['tests']['label_vs_occupation']=small
# N=1 closed form constant diagonal
n1=[]
for g in [0.0,0.2,1.1,3.7]:
    z=Z_occ(1,g); closed=math.exp(g*A[0,0]/2)
    n1.append({'g':g,'Z':z,'closed':closed,'abs_diff':abs(z-closed)})
out['tests']['N1_closed']=n1
# D12 occupation-weight invariance for deterministic nontrivial examples
sym=[]
for N,n_tuple in [(4,(1,0,1,0,0,1,0,0,0,0,1,0)),(7,(2,0,1,0,1,0,0,1,0,1,1,0))]:
    n=np.array(n_tuple,int); p=n/N; base=(N*3.2/2)*float(p@A@p)+log_multinomial(n)
    errs=[]
    for eps in [1,-1]:
        for a in range(12):
            q=d12_perm_counts(n,a,eps); pq=q/N
            val=(N*3.2/2)*float(pq@A@pq)+log_multinomial(q)
            errs.append(abs(val-base))
    sym.append({'N':N,'max_abs_logweight_diff':max(errs)})
out['tests']['D12_weight_invariance']=sym
# determinant identity on deterministic interior p family
B0=np.column_stack([np.eye(12)[:,i]-np.eye(12)[:,11] for i in range(11)])
dets=[]
for t in np.linspace(0.03,0.21,12):
    q=np.arange(1,13,dtype=float); q=(q/q.sum())
    p=(1-t)*np.ones(12)/12+t*q
    Sig=np.diag(p)-np.outer(p,p); M=X.T@Sig@X; g=1.7
    Hc=B0.T@(np.diag(1/p)-g*A)@B0
    lhs=np.linalg.det(Hc)*np.prod(p); rhs=np.linalg.det(np.eye(7)-g*M)
    dets.append({'t':float(t),'lhs':float(lhs),'rhs':float(rhs),'abs_diff':float(abs(lhs-rhs))})
out['tests']['determinant_identity']=dets
# exact finite-N fluctuation response versus central finite difference of expectation
resp=[]
for N in [2,3,4]:
    g=0.9; f=np.array([0.013,-0.009,0.007,-0.004,0.006,-0.005,0.003]); mean,cov=mean_cov_mu(N,g,f)
    h=2e-5; J=np.zeros((7,7))
    for k in range(7):
        e=np.zeros(7); e[k]=h
        mp1,_=mean_cov_mu(N,g,f+e); mm1,_=mean_cov_mu(N,g,f-e); J[:,k]=(mp1-mm1)/(2*h)
    err=J-N*cov
    resp.append({'N':N,'max_abs_J_minus_NCov':float(np.abs(err).max()),'fro_error':float(np.linalg.norm(err))})
out['tests']['finite_N_response']=resp
# total covariance identity diagnostic from mixture conditional, numerical assembly
mix=[]
for N in [2,4]:
    g=0.9; mean,cov=mean_cov_mu(N,g,np.zeros(7)); covtheta=(g/N)*np.eye(7)+g*g*cov
    # compare with law-of-total-covariance independently assembled as within + between conditional theta means
    ts=occ_terms(N,g); logs=np.array([z for z,_,__ in ts]); ws=np.exp(logs-logsumexp(logs)); mus=np.array([mu for _,mu,__ in ts]); means_theta=g*mus; mt=ws@means_theta; Y=means_theta-mt
    assembled=(g/N)*np.eye(7)+np.einsum('i,ia,ib->ab',ws,Y,Y)
    mix.append({'N':N,'max_abs_diff':float(np.abs(covtheta-assembled).max())})
out['tests']['auxiliary_covariance_identity']=mix
# gates
mx_label=max(x['rel_diff'] for x in small); mx_det=max(x['abs_diff'] for x in dets); mx_resp=max(x['max_abs_J_minus_NCov'] for x in resp); mx_sym=max(x['max_abs_logweight_diff'] for x in sym)
out['gates']={
 'label_vs_occupation_rel_lt_1e-12':mx_label<1e-12,
 'D12_logweight_lt_1e-12':mx_sym<1e-12,
 'det_identity_abs_lt_1e-11':mx_det<1e-11,
 'response_fd_abs_lt_5e-8':mx_resp<5e-8,
 'aux_cov_exact_float_lt_1e-14':max(x['max_abs_diff'] for x in mix)<1e-14}
out['pass']=all(out['gates'].values())
out['scope']='Deterministic implementation diagnostics only; analytic identities are proved in MP7-031--033. No sampling or physical interpretation.'
print(json.dumps(out,indent=2))
