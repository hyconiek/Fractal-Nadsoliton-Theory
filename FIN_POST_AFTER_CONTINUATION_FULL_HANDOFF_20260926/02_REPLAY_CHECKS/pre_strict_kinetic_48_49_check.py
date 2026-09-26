#!/usr/bin/env python3
import numpy as np, math

N=12

# S12 theorem numerical fixture.
gamma=1.7
u=np.ones(N)/N
Qreset=gamma*(np.ones((N,1))*u[None,:]-np.eye(N))
assert np.linalg.norm(Qreset@np.ones(N))<1e-14
e=np.linalg.eigvalsh(Qreset)
assert np.sum(np.isclose(e,0,atol=1e-12))==1
assert np.sum(np.isclose(e,-gamma,atol=1e-12))==N-1

# A nonuniform p cannot be invariant under all permutations.
p=np.arange(1,N+1,dtype=float)
p/=p.sum()
assert not np.allclose(p,np.ones(N)/N)

# Strict kinetic splitting.
def strict_kernel():
    return np.array([[0.0 if i==j else
        math.cos(0.18575*min(abs(i-j),N-abs(i-j))+0.1625)/
        (1+min(abs(i-j),N-abs(i-j))**1.8)
        for j in range(N)] for i in range(N)],dtype=float)

W=strict_kernel()
A=np.diag(W.sum(axis=1))-W
lam=np.fft.fft(A[0]).real[:7]

eta=.43
PC=np.eye(N)-np.ones((N,N))/N
Q=-gamma*PC-eta*A

off=Q.copy()
np.fill_diagonal(off,np.nan)
assert np.nanmin(off)>0
assert np.linalg.norm(Q@np.ones(N))<1e-12
assert np.linalg.norm(Q-Q.T)<1e-12

ev=np.sort(np.linalg.eigvalsh(Q))
expected=np.sort(np.r_[0,-(gamma+eta*lam[1]),-(gamma+eta*lam[1]),
                         -(gamma+eta*lam[2]),-(gamma+eta*lam[2]),
                         -(gamma+eta*lam[3]),-(gamma+eta*lam[3]),
                         -(gamma+eta*lam[4]),-(gamma+eta*lam[4]),
                         -(gamma+eta*lam[5]),-(gamma+eta*lam[5]),
                         -(gamma+eta*lam[6])])
assert np.max(np.abs(ev-expected))<1e-12

rates=gamma+eta*lam[1:]
assert np.all(np.diff(rates)>0)

print("PASS")
print("S12_reset_centered_rate",gamma)
print("strict_lambdas",lam[1:].tolist())
print("strict_sector_rates",rates.tolist())
print("minimum_offdiag_rate",float(np.nanmin(off)))
print("surviving_ratio_eta_over_gamma",eta/gamma)
