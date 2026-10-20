"""Bounded independent checks for the FIN Z12 conditional memory / RP004 algebra.
Not a formal PDE proof nor a live-repository custody check."""
import json, math
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.linalg import expm

BASE = Path(__file__).resolve().parent
k = np.array([0, 0.46998567264502017, 0.19204355169010282,
              0.09142861427792497, 0.04702916874565040,
              0.02413122336363006, 0.011070817321442113])
assert all(k[1:] > 0)
s = 2*sum(k[1:6]) + k[6]
W = np.zeros((12,12))
for i in range(12):
    for j in range(12):
        if i != j: W[i,j] = k[min(abs(i-j), 12-abs(i-j))]
A = s*np.eye(12)-W
D, C = A[::2,::2], A[::2,1::2]
assert np.max(np.abs(D-A[1::2,1::2])) < 1e-14
assert np.max(np.abs(C-A[1::2,::2].T)) < 1e-14
assert np.max(np.abs(A@np.ones(12))) < 1e-13
assert np.linalg.eigvalsh(A)[0] > -1e-13
assert np.linalg.eigvalsh(D)[0] > 1
u3 = np.array([1.,-1.,1.,-1.,1.,-1.])/math.sqrt(6)
assert np.linalg.norm(C.T@u3) < 1e-13
assert abs(float(u3@D@u3)-1.9614068619764466) < 1e-12
# Algebraic memory identity by direct block ODE discretization at three times;
# M(t) PSD and zero exactly on alternating Fourier channel.
lam, amp = [], []
for j in range(6):
    v = np.exp(2j*np.pi*j*np.arange(6)/6)/math.sqrt(6)
    lam.append(float(np.vdot(v,D@v).real))
    amp.append(float(max(0.0,np.vdot(v,(C@C.T)@v).real)))
for t in (0, 0.2, 1.0, 3.0):
    ld, V = np.linalg.eigh(D)
    expD = (V*np.exp(-t*ld))@V.T
    M = C@expD@C.T
    assert np.linalg.eigvalsh(M).min() > -1e-12
    assert np.linalg.norm(M@u3) < 1e-12
    v0 = np.ones(6)/math.sqrt(6)
    assert abs(float(v0@M@v0)-amp[0]*math.exp(-lam[0]*t)) < 1e-12
# At lambda=0, static Schur of the retained block has a uniform null vector.
Q = D-C@np.linalg.solve(D,C.T)
assert np.linalg.norm(Q@np.ones(6)) < 1e-12
assert np.linalg.eigvalsh(Q).min() > -1e-12
# Non-equivalence of small-time diffusive and coherent endpoint readout.
i, j = 0, 1
lag=1e-4
P=expm(-lag*A)
U=expm(-1j*lag*A)
assert abs(P[i,j]/lag-W[i,j]) < 1e-3
assert abs(abs(U[i,j])**2/lag**2-W[i,j]**2) < 1e-5
assert np.max(abs(np.sum(abs(U)**2,axis=0)-1)) < 1e-12
# Check original frozen polynomial energy and 2nd-harmonic forcing symbolically.
phi, psi, aa, gg4, gg6 = sp.symbols('phi psi a g4 g6', real=True)
V=sp.Rational(1,2)*(phi**2-1)**2 + 3*phi**2*psi**2 + aa*psi**2/2+gg4*psi**4/4+gg6*psi**6/6
assert sp.simplify(sp.diff(V,phi)-(2*phi*(phi**2-1)+6*phi*psi**2)) == 0
assert sp.simplify(sp.diff(V,psi)-((aa+6*phi**2)*psi+gg4*psi**3+gg6*psi**5)) == 0
KK, ff, theta = sp.symbols('K f theta',real=True)
assert sp.simplify(6*KK*ff**2*sp.cos(theta)**2-3*KK*ff**2*(1+sp.cos(2*theta)))==0
# Optional independent frozen RP reference diagnostics.
a_vals = (-.75, -.5, -.25)
rows=[]
for a in a_vals:
    omega=math.sqrt(a+2)
    Omega=2*omega
    q=2*math.sqrt(a+1)
    vg=q/Omega
    F4=math.pi*q*(q*q+4)/(6*math.sinh(math.pi*q/2))
    gamma=81*q*(q*q+4)*F4**2/(2048*(1+q*q))
    assert vg > 1/3
    assert gamma > 125/73728
    rows.append(dict(a=a,omega=omega,Omega=Omega,q=q,group_velocity=vg,Gamma=gamma,linear_total_energy_growth_rate=Omega*gamma/2))
# Independent constant comparison of the claimed tiny flux coefficient.
claimed_constant = 625/12386304
assert claimed_constant>0
result={
    'source_weights':k[1:].tolist(), 'strict_s':float(s),
    'memory_partition':'even nodes vs odd nodes of Z12, for explicitly chosen diffusive P(t)=exp(-tA)',
    'D_eigenvalues_by_Z6_Fourier':lam,
    'memory_amplitudes_by_Z6_Fourier':amp,
    'decoupled_channel_index':3,
    'decoupled_channel_vector':u3.tolist(),
    'static_schur_min_eig':float(np.linalg.eigvalsh(Q)[0]),
    'radiation_band':rows,
    'claimed_flux_floor_coefficient_arithmetic_only':claimed_constant,
    'replay_checks':'matrix identities, mode decoupling, positive memory, phase velocity and coefficient checks passed',
    'limits':'not a current-repo verification, not a dispersion estimate, not an independently audited moving-cone theorem, not a source law'
}
(BASE/'RESULTS.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
print(json.dumps(result,indent=2))
