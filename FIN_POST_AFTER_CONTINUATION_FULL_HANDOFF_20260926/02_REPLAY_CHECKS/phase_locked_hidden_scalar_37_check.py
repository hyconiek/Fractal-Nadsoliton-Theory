#!/usr/bin/env python3
import numpy as np, math

# Real H4 representation: (k1c,k1s,k2c,k2s)
def R2(theta):
    return np.array([[math.cos(theta),-math.sin(theta)],
                     [math.sin(theta), math.cos(theta)]])

alpha=2*math.pi/12
R=np.block([[R2(alpha),np.zeros((2,2))],
            [np.zeros((2,2)),R2(2*alpha)]])
F=np.diag([1,-1,1,-1])  # reflection in real Fourier basis

# No invariant linear covector/vector.
M=np.vstack([R.T-np.eye(4),F.T-np.eye(4)])
sv=np.linalg.svd(M,compute_uv=False)
nullity=4-np.linalg.matrix_rank(M,tol=1e-12)
assert nullity==0

# Quadratic invariant symmetric matrices: solve R^T Q R=Q, F^T Q F=Q.
# Parametrize symmetric Q by 10 coordinates.
pairs=[(i,j) for i in range(4) for j in range(i,4)]
def Q_from(v):
    Q=np.zeros((4,4))
    for x,(i,j) in zip(v,pairs):
        Q[i,j]=Q[j,i]=x
    return Q

constraints=[]
for g in (R,F):
    for row in range(4):
        for col in range(row,4):
            coeff=[]
            for p in range(len(pairs)):
                e=np.zeros(len(pairs));e[p]=1
                Q=Q_from(e)
                D=g.T@Q@g-Q
                coeff.append(D[row,col])
            constraints.append(coeff)
C=np.array(constraints)
rank=np.linalg.matrix_rank(C,tol=1e-11)
qnull=len(pairs)-rank
assert qnull==2

# State-relative invariants.
rng=np.random.default_rng(20260926)
maxerr=0.0
for _ in range(1000):
    z1=rng.normal()+1j*rng.normal()
    z2=rng.normal()+1j*rng.normal()
    th=rng.uniform(-math.pi,math.pi)
    chi=np.exp(1j*th)
    s=np.array([
        np.real(z1*np.conj(chi)),
        np.real(z2*np.conj(chi)**2)
    ])

    # random D12 rotation
    a=int(rng.integers(0,12))
    ang=2*math.pi*a/12
    z1r=np.exp(1j*ang)*z1
    z2r=np.exp(2j*ang)*z2
    chir=np.exp(1j*ang)*chi
    sr=np.array([
        np.real(z1r*np.conj(chir)),
        np.real(z2r*np.conj(chir)**2)
    ])
    maxerr=max(maxerr,float(np.max(np.abs(s-sr))))

    # reflection
    z1f=np.conj(z1)
    z2f=np.conj(z2)
    chif=np.conj(chi)
    sf=np.array([
        np.real(z1f*np.conj(chif)),
        np.real(z2f*np.conj(chif)**2)
    ])
    maxerr=max(maxerr,float(np.max(np.abs(s-sf))))

assert maxerr<1e-12

print("PASS")
print("linear_invariant_nullity",nullity)
print("quadratic_invariant_dimension",qnull)
print("state_relative_invariance_error",maxerr)
print("quadratic_basis_form","span{diag(1,1,0,0), diag(0,0,1,1)}")
