#!/usr/bin/env python3
import math, numpy as np

N=12
def strict_kernel():
    return np.array([[0.0 if i==j else
        math.cos(0.18575*min(abs(i-j),N-abs(i-j))+0.1625)/
        (1+min(abs(i-j),N-abs(i-j))**1.8)
        for j in range(N)] for i in range(N)],dtype=float)

W=strict_kernel()
A=np.diag(W.sum(axis=1))-W
L=np.fft.fft(A[0]).real[:7]
j=np.arange(N)

cols=[]
for k in (3,4,5):
    cols += [np.sqrt(L[k]/6)*np.cos(2*np.pi*k*j/N),
             np.sqrt(L[k]/6)*np.sin(2*np.pi*k*j/N)]
cols += [np.sqrt(L[6]/12)*(-1.)**j]
X=np.column_stack(cols)

hidden=[]
for k in (1,2):
    hidden += [np.sqrt(2/N)*np.cos(2*np.pi*k*j/N),
               np.sqrt(2/N)*np.sin(2*np.pi*k*j/N)]
Y=np.column_stack(hidden)

T=[X.T@np.diag(Y[:,a])@X for a in range(4)]
norms=np.array([np.linalg.norm(t,'fro') for t in T])
assert np.all(norms>1.4)

# Quadratic witness f=x^T T_a x: O_a f = 2||T_a||_F^2.
witness=2*norms**2

# OU covariance identities: C(t)=Sigma+e^-2t(C0-Sigma);
# Cov(t,s)=e^{-(t-s)}C(s). Verify against closed-form propagation.
Sigma=np.eye(4)/12
C0=np.diag([.2,.1,.05,.15])
m0=np.array([.3,-.2,.1,.4])

for s,t in [(0.,.3),(.2,.8),(.7,1.1)]:
    Cs=Sigma+math.exp(-2*s)*(C0-Sigma)
    Ct=Sigma+math.exp(-2*t)*(C0-Sigma)
    cross=math.exp(-(t-s))*Cs
    # PSD block covariance consistency.
    block=np.block([[Cs,cross.T],[cross,Ct]])
    assert np.linalg.eigvalsh(block).min()>-1e-12
    mt=math.exp(-t)*m0
    assert np.allclose(mt, np.exp(-t)*m0)

print("PASS")
print("T_frobenius_norms", norms.tolist())
print("quadratic_witness_coefficients", witness.tolist())
print("uniform_Sigma_diag", np.diag(Sigma).tolist())
