#!/usr/bin/env python3
import numpy as np, math
from scipy.special import softmax
from scipy.optimize import root

Nlabel=12

def strict_kernel():
    return np.array([[0.0 if i==j else
        math.cos(0.18575*min(abs(i-j),Nlabel-abs(i-j))+0.1625)/
        (1+min(abs(i-j),Nlabel-abs(i-j))**1.8)
        for j in range(Nlabel)] for i in range(Nlabel)],float)

W=strict_kernel()
A=np.diag(W.sum(axis=1))-W
lam=np.fft.fft(A[0]).real[:7]
j=np.arange(Nlabel)
cols=[]
for k in (3,4,5):
    cols += [
        np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*j/Nlabel),
        np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*j/Nlabel)
    ]
cols += [np.sqrt(lam[6]/12)*(-1.)**j]
X=np.column_stack(cols)
A7=X@X.T
C4=X[:,[0,2,4,6]]

# circulant diagonal
assert np.ptp(np.diag(A7))<1e-14

sf=np.array([1.3643114281782411,1.4331027014285429,
             1.4080457012550378,1.0056820773201627])
gf=3.5156447168395917
v=np.array([0.50736753967724968,0.52868565531381106,
            0.5537606940687495,0.39549810524420642])

alo=-0.212571640084013874
ahi=-0.212571626020459008
blo=0.118982527450498775
bhi=0.118982667103223664

def coeff(a,b):
    return (2*b/3)*((-2*a/b)**1.5)

corners=[coeff(a,b) for a in (alo,ahi) for b in (blo,bhi)]
Clo,Chi=min(corners),max(corners)
assert 0.5357594<Clo<Chi<0.5357599

def m4(s):
    return C4.T@softmax(C4@s)

def F4(s,g):
    return s/g-m4(s)

def Phi(s,g):
    h=C4@s
    M=h.max()
    return .5*s@s/g-(M+np.log(np.exp(h-M).mean()))

amid=.5*(alo+ahi); bmid=.5*(blo+bhi)
ratios=[]
for delta in (1e-6,1e-5,1e-4,1e-3):
    amp=math.sqrt(-2*amid*delta/bmid)
    sols=[]
    for sign in (1,-1):
        z=root(lambda ss:F4(ss,gf+delta),sf+sign*amp*v,tol=1e-12)
        assert z.success or np.linalg.norm(F4(z.x,gf+delta))<1e-10
        sols.append(z.x)
    ph=sorted(Phi(s,gf+delta) for s in sols)
    ratios.append((ph[1]-ph[0])/delta**1.5)

assert abs(ratios[0]-.5*(Clo+Chi))<2e-6

# Exact Gibbs conditional detailed-balance test on microstates for small N.
g=3.6
Nc=5
rng=np.random.default_rng(20260926)

def logpi(x):
    x=np.asarray(x,dtype=int)
    S=0.0
    for a in range(Nc):
        for b in range(Nc):
            S+=A7[x[a],x[b]]
    return -Nc*math.log(12)+g*S/(2*Nc)

for _ in range(100):
    x=rng.integers(0,Nlabel,size=Nc)
    a=int(rng.integers(0,Nc))
    new=int(rng.integers(0,Nlabel))
    xp=x.copy(); xp[a]=new

    others=np.delete(x,a)
    counts=np.bincount(others,minlength=Nlabel)
    qc=softmax((g/Nc)*(A7@counts))

    # same "others" conditional for both directions
    lhs=math.exp(logpi(x))*qc[new]
    rhs=math.exp(logpi(xp))*qc[x[a]]
    assert abs(lhs-rhs) <= 2e-12*max(lhs,rhs,1.0)

# First 1/N drift correction.
p=rng.uniform(.2,2,size=Nlabel); p/=p.sum()
q=softmax(g*A7@p)
Sq=np.diag(q)-np.outer(q,q)
pred=-g*Sq@A7@p
errs=[]
for NN in (100,200,400,800):
    qb=np.zeros(Nlabel)
    for i in range(Nlabel):
        qi=softmax(g*A7@p-(g/NN)*A7[:,i])
        qb+=p[i]*qi
    scaled=NN*(qb-q)
    errs.append(np.linalg.norm(scaled-pred))
assert errs[-1] < errs[0]/6

print("PASS")
print("fold_barrier_coefficient_interval",[Clo,Chi])
print("crossover_coefficient_interval",
      [(1/Chi)**(2/3),(1/Clo)**(2/3)])
print("barrier_replay_ratios",ratios)
print("finite_N_microstate_DB_trials",100)
print("drift_correction_errors",errs)
