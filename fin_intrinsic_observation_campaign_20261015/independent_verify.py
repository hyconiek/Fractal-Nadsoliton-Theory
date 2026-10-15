"""Independent OBS checks: rebuild A directly; no import of observation_core or theory_core.
All finite checks are exact fractions or local float diagnostics, not asymptotic proof.
"""
from __future__ import annotations
import math,json,hashlib
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from scipy.integrate import quad
from scipy.special import ndtr
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,condition,details=None):
    checks.append({'check':name,'pass':bool(condition),'detail':details})
    if not condition:print('FAIL',name,details)

def H(N,M,k,h):
    if not max(0,k-(N-M))<=h<=min(k,M):return Fraction(0)
    return Fraction(math.comb(M,h)*math.comb(N-M,k-h),math.comb(N,k))

def H2(N,M,k1,k2,h1,h2):return H(N,M,k1,h1)*H(N-k1,M-h1,k2,h2)

# Exact counts and bivariate hypergeom, no float combinatoric shortcuts.
for N,M,k in [(12,6,4),(24,14,6),(36,15,7),(120,60,13)]:
    p=[H(N,M,k,h) for h in range(k+1)]
    mean=sum(Fraction(h)*x for h,x in enumerate(p))
    var=sum((h-mean)**2*x for h,x in enumerate(p))
    check(f'combinatorics_{N}_normal',sum(p)==1)
    check(f'combinatorics_{N}_mean',mean==Fraction(k*M,N))
    check(f'combinatorics_{N}_var',var==Fraction(k*M*(N-M)*(N-k),N*N*(N-1)))
for N,M,k1,k2 in [(12,6,2,3),(24,13,4,5),(120,61,12,11)]:
    pj=[[H2(N,M,k1,k2,h1,h2) for h2 in range(k2+1)] for h1 in range(k1+1)]
    check(f'joint_{N}_total',sum(map(sum,pj))==1)
    cov=sum((Fraction(h1)-Fraction(k1*M,N))*(Fraction(h2)-Fraction(k2*M,N))*pj[h1][h2] for h1 in range(k1+1) for h2 in range(k2+1))
    check(f'joint_{N}_cov',cov==-Fraction(k1*k2*M*(N-M),N*N*(N-1)))
    for h1 in range(k1+1):
        check(f'joint_{N}_marginal_{h1}',sum(pj[h1])==H(N,M,k1,h1))

# Independent Laplacian spectrum via symmetric cos kernel, no fft.
q=12
w=np.zeros((q,q))
for i in range(q):
    for j in range(q):
        if i==j:continue
        d=min((i-j)%q,(j-i)%q)
        w[i,j]=math.cos(.18575*d+.16250)/(1+d**1.8)
Lap=np.diag(w.sum(1))-w
lam=np.linalg.eigvalsh(Lap)
# Positive laplacian nonzero largest k6, but reconstruct features with orthonormal Fourier.
six=np.array([(-1)**j for j in range(q)])
chi=float(six@Lap@six/(q*q))
gc=1/chi
modes=[]
for k in [3,4,5]:
    eig=float(np.cos(2*np.pi*k*np.arange(q)/q)@Lap@np.cos(2*np.pi*k*np.arange(q)/q)/6)
    modes.extend([np.sqrt(eig/6)*np.cos(2*np.pi*k*np.arange(q)/q),np.sqrt(eig/6)*np.sin(2*np.pi*k*np.arange(q)/q)])
eig6=float(six@Lap@six/q)
modes.append(np.sqrt(eig6/q)*six)
X=np.column_stack(modes); A=X@X.T
check('q12_rank7',np.linalg.matrix_rank(A,1e-10)==7)
check('A_circulant',max(abs(A[i,j]-A[0,(j-i)%q]) for i in range(q) for j in range(q))<2e-14)
check('critical_gain',abs(gc-5.123427551398615)<2e-11)
check('equal_diag',np.ptp(np.diag(A))<2e-14)

def probabilities(n,lab,g):
    o=n.copy();o[lab]-=1
    logits=g*(A@o)/int(n.sum());z=logits-logsumexp(logits)
    return np.exp(z)

def separate_label_stats(lab,k,g):
    N=len(lab);n=np.bincount(lab,minlength=q);sig=six[lab]
    co=np.array([(1/k-1/N) if i<k else -1/N for i in range(N)])
    Lm=Lw=Bm=Bw=Bcross=0.
    for a,i in enumerate(lab):
        P=probabilities(n,i,g)
        for j in range(q):
            rate=P[j]
            dm=(six[j]-six[i])/N
            dv=co[a]*(six[j]-six[i])
            Lm+=rate*dm
            Lw+=rate*dv
            Bm+=rate*dm*dm
            Bw+=rate*dv*dv
            Bcross+=rate*dm*dv
    m=float(np.mean(sig));v=float(np.mean(sig[:k])-m)
    return Lm,Lw,v,Bm,Bw,Bcross

rng=np.random.default_rng(422)
for N,k in [(12,1),(24,2),(120,3),(120,12)]:
    lab=rng.integers(0,q,size=N)
    z=separate_label_stats(lab,k,gc)
    lab2=lab.copy(); # permutation that changes the tagged group only after checking count drift
    rng.shuffle(lab2)
    z2=separate_label_stats(lab2,k,gc)
    check(f'count_generator_perm_{N}_{k}',abs(z[0]-z2[0])<1e-12)
    check(f'count_bracket_perm_{N}_{k}',abs(z[3]-z2[3])<1e-12)
    check(f'tag_contrast_bound_{N}_{k}',abs(z[1]+z[2])*N<50) # strong diagnostic only
    check(f'bracket_psd_{N}_{k}',z[3]>=0 and z[4]>=0 and z[3]*z[4]>=z[5]**2-1e-12)
    check(f'cov_common_O_Ninv_{N}_{k}',abs(z[5])*N<10)

# At uniform n_i N/12, parity total bracket nΓ(m) close to 2.
for N in (120,1200):
    lab=np.repeat(np.arange(q),N//q)
    z=separate_label_stats(lab,1,gc)
    check(f'uniform_global_bracket_{N}',abs(N*z[3]-2)<3/N)
    check(f'uniform_tag_bracket_{N}',abs(z[4]-2*(1-1/N))<3/N)

norm=quad(lambda z:math.exp(-z**4/12),-np.inf,np.inf,epsabs=3e-12)[0]
risk=2*quad(lambda z:math.exp(-z**4/12)*ndtr(-z),0,np.inf,epsabs=3e-12)[0]/norm
source=json.loads((HERE/'NUMERICAL_RESULTS.json').read_text())
check('independent_quartic_risk',abs(risk-source['bounded_sign_risk_nominal']['alpha']['1.0'])<2e-12)
check('universal_noise_variance',abs(source['critical_private_ou_increment_variance_alpha1_delta_modeltime1']-2*(1-math.exp(-1)))<1e-14)
check('quarter_root_coefficient',abs(source['gc_nominal']-gc)<1e-11)
# Verify saved finite formulas across all selected fixtures
for t in source['exact_hypergeom_checks']:
    check(f'saved_hg_{t["N"]}',abs(t['sum']-1)<1e-13 and abs(t['variance_error'])<1e-12)
for t in source['two_block_checks']:
    check(f'saved_joint_{t["N"]}',abs(t['total']-1)<1e-13 and abs(t['cov_error'])<1e-12)
for t in source['labelled_generator_checks']:
    check(f'saved_projection_{t["N"]}_{t["k"]}',abs(t['count_projection_error'])<2e-13 and abs(t['count_bracket_error'])<2e-13)

out={'passed':sum(c['pass'] for c in checks),'total':len(checks),'independent_of_main_module':True,
     'scope':'exact Fraction hypergeometric plus separate binary64 generator construction; not stochastic theorem proof',
     'checks':checks}
(HERE/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n')
print('checks',out['passed'],'/',out['total'],'gc',gc,'risk',risk)
if out['passed']!=out['total']:raise SystemExit(2)
