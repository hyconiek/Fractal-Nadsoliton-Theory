"""fin_core.py -- seed reference implementation for the FIN <-> known-physics bridge campaign.

Independent re-implementation (written from the repository's model definition, NOT copied
from its result files).  Reproduces: V_loc(G), B4, rho_N (k=4 slow mode), R_k(N) affine law.

Model (repo convention):
  labels j in Z12; strict kernel W_d = cos(0.18575 d + 0.1625)/(1 + d^1.8) (d = ring distance, W_0 = 0)
  A = Laplacian(W) (circulant);  L_k = Re FFT(A[0])_k
  A7 = X7 X7^T = projection of A onto Fourier modes k = 3,4,5,6 (eigenvalue L_k on mode k)
  Gibbs law on counts n (N copies):  pi(n) ~ N!/prod n_j! * exp[(g/(2N)) n^T A7 n]
  leave-one-out heat bath: pick a copy at rate 1, remove it, resample label ~ softmax_j[(g/N) (A7 n')_j]
  mean-field functional: V_g(p) = sum p log(12 p) - (g/2) p^T A7 p
"""
import itertools
import math

import numpy as np
import scipy.sparse as sp
from scipy.optimize import brentq, minimize
from scipy.sparse.linalg import eigsh

Q = 12
G_FROZEN = 5.145228719489142  # repo's frozen working gain (origin not documented in the reports read)


# ----------------------------------------------------------------------------- coupling matrix
def strict_L(omega=0.18575, phi=0.1625, eta=1.8, beta=1.0):
    d = lambda i, k: min(abs(i - k), Q - abs(i - k))
    W = np.array([[0.0 if i == k else math.cos(omega * d(i, k) + phi) / (1 + beta * d(i, k) ** eta)
                   for k in range(Q)] for i in range(Q)])
    A = np.diag(W.sum(1)) - W
    return np.fft.fft(A[0]).real  # L_k, k = 0..11


def A_from_lams(lam):
    """Circulant matrix with eigenvalue lam[k] on Fourier mode k (k in 1..6), 0 elsewhere."""
    first = np.zeros(Q)
    for k, l in lam.items():
        w = 1 if k == 6 else 2
        first += w * l * np.array([math.cos(2 * math.pi * k * d / Q) for d in range(Q)]) / Q
    return np.array([[first[(a - b) % Q] for b in range(Q)] for a in range(Q)])


def A7_default():
    L = strict_L()
    return A_from_lams({3: L[3], 4: L[4], 5: L[5], 6: L[6]}), L


def projector_countermodel():
    """Same cutoff (k>=3) but all four eigenvalues equal to their mean: kills the strict-kernel tilt."""
    _, L = A7_default()
    lbar = float(np.mean([L[3], L[4], L[5], L[6]]))
    return A_from_lams({k: lbar for k in (3, 4, 5, 6)}), lbar


# ----------------------------------------------------------------------------- mean-field statics
def _softmax(e):
    e = e - e.max()
    x = np.exp(e)
    return x / x.sum()


def _Vg(eta, g, A):
    p = _softmax(eta)
    pc = np.clip(p, 1e-300, None)
    Ap = A @ p
    v = float(np.sum(pc * np.log(Q * pc)) - 0.5 * g * p @ Ap)
    dVdp = np.log(Q * pc) + 1 - g * Ap
    return v, p * (dVdp - p @ dVdp)


def V_localized(g, A, starts=3):
    """Minimum of V_g over NON-uniform local minima found from delta-like starts (None if none)."""
    best = None
    for s in range(starts):
        eta0 = np.full(Q, -5.0)
        eta0[0] = 0.0
        if s:
            eta0 += 0.2 * np.random.default_rng(s).standard_normal(Q)
        r = minimize(_Vg, eta0, args=(g, A), jac=True, method="L-BFGS-B",
                     options={"ftol": 1e-15, "gtol": 1e-11, "maxiter": 2000})
        p = _softmax(r.x)
        if np.abs(p - 1 / Q).sum() > 2e-2 and (best is None or r.fun < best[0]):
            best = (r.fun, p)
    return best


def g_eq(A, lo=0.5, hi=40.0, n=80):
    """Coexistence gain: V_localized = 0 (uniform has V = 0).  Robust to uniform-minimum collapse."""
    f = lambda g: (lambda b: 1e-2 if b is None else b[0])(V_localized(g, A))
    gs = np.linspace(lo, hi, n)
    vs = [f(g) for g in gs]
    for a, b, va, vb in zip(gs[:-1], gs[1:], vs[:-1], vs[1:]):
        if va > 0 and vb <= 0:
            return brentq(f, a, b, xtol=1e-10)
    return None


def V_d4_saddle(g, A):
    """Minimum of V_g on the reflection-symmetric subspace p_j = p_{4-j}: the 'd4' communication saddle."""
    orbits, seen = [], set()
    for a in range(Q):
        if a in seen:
            continue
        orb = sorted({a, (4 - a) % Q})
        orbits.append(orb)
        seen |= set(orb)

    def expand(eta):
        s = np.exp(eta)
        p = np.zeros(Q)
        for si, orb in zip(s, orbits):
            for a in orb:
                p[a] = si
        return p / p.sum()

    def f(eta):
        p = expand(eta)
        pc = np.clip(p, 1e-300, None)
        return float(np.sum(pc * np.log(Q * pc)) - 0.5 * g * p @ A @ p)

    eta0 = np.full(len(orbits), -6.0)
    for i, orb in enumerate(orbits):
        if 0 in orb:
            eta0[i] = 0.0
    r = minimize(f, eta0, method="BFGS", options={"gtol": 1e-12})
    return r.fun, expand(r.x)


# ----------------------------------------------------------------------------- finite-N chain
def states(N):
    S = []
    for c in itertools.combinations_with_replacement(range(Q), N):
        n = [0] * Q
        for x in c:
            n[x] += 1
        S.append(tuple(n))
    return S


def build_generator(N, g, A, kind="heatbath"):
    """Exact count-space generator.  kind: heatbath (FIN convention) | metropolis | barker.
    Each copy is selected at rate 1; time unit = one update per copy."""
    S = states(N)
    idx = {s: i for i, s in enumerate(S)}
    M = len(S)
    rows, cols, vals = [], [], []
    diag = np.zeros(M)
    for i, s in enumerate(S):
        n = np.array(s, dtype=float)
        for a in range(Q):
            if s[a] == 0:
                continue
            nm = n.copy()
            nm[a] -= 1
            h = (g / N) * (A @ nm)
            if kind == "heatbath":
                hh = h - h.max()
                q = np.exp(hh)
                q /= q.sum()
            else:
                d = h - h[a]
                acc = np.minimum(1.0, np.exp(d)) if kind == "metropolis" else 1 / (1 + np.exp(-d))
                q = acc / (Q - 1)
            for j in range(Q):
                if j == a:
                    continue
                t = list(s)
                t[a] -= 1
                t[j] += 1
                r = s[a] * q[j]
                rows.append(i)
                cols.append(idx[tuple(t)])
                vals.append(r)
                diag[i] -= r
    Qm = sp.csr_matrix((vals, (rows, cols)), shape=(M, M)) + sp.diags(diag)
    lw = np.array([math.lgamma(N + 1) - sum(math.lgamma(x + 1) for x in s)
                   + (g / (2 * N)) * float(np.array(s) @ A @ np.array(s)) for s in S])
    lw -= lw.max()
    pi = np.exp(lw)
    pi /= pi.sum()
    return S, idx, Qm, pi


def slow_modes(N, g, A, kind="heatbath"):
    """Returns (rho, {k: R_k}, detailed_balance_asymmetry, start_of_fast_band).
    rho = slow eigenvalue of the Fourier-k=4 sector; R_k = eigenvalue_k / rho (R_4 = 1 by definition).
    k is assigned by the eigenvalue of the label-rotation operator on each degenerate slow multiplet."""
    S, idx, Qm, pi = build_generator(N, g, A, kind)
    Dh, Dih = sp.diags(np.sqrt(pi)), sp.diags(1 / np.sqrt(pi))
    Ss = Dh @ Qm @ Dih
    Ss = (Ss + Ss.T) / 2
    asym = abs(sp.diags(pi) @ Qm - (sp.diags(pi) @ Qm).T).max()
    w, v = eigsh(-Ss, k=13, sigma=-1e-4, which="LM")
    o = np.argsort(w)
    w, v = w[o], v[:, o]
    perm = np.array([idx[tuple(s[(j - 1) % Q] for j in range(Q))] for s in S])
    groups, i = [], 1
    while i < 12:
        j = i + 1
        while j < 12 and abs(w[j] - w[i]) < 1e-7 * max(1, abs(w[i])):
            j += 1
        groups.append((i, j))
        i = j
    res = {}
    for a, b in groups:
        V = v[:, a:b]
        TV = np.empty_like(V)
        TV[perm] = V
        ev = np.linalg.eigvals(V.T @ TV)
        k = int(round(abs(np.angle(ev[0])) * 12 / (2 * math.pi)))
        res[k] = float(w[a])
    rho = res[4]
    return rho, {k: res[k] / rho for k in sorted(res)}, float(asym), float(w[12])


# ----------------------------------------------------------------------------- exact static Fourier variances
def static_S(N, g, A):
    """Exact equal-time S_k = E|sum_j n_j exp(2 pi i k j/12)|^2 / N under the finite-N Gibbs law (S_k = 1 at g = 0)."""
    S = np.zeros((math.comb(N + Q - 1, Q - 1), Q), dtype=np.int16)
    for i, c in enumerate(itertools.combinations_with_replacement(range(Q), N)):
        for x in c:
            S[i, x] += 1
    n = S.astype(float)
    lgf = np.array([math.lgamma(x + 1) for x in range(N + 1)])
    lw = math.lgamma(N + 1) - lgf[S].sum(axis=1) + (g / (2 * N)) * np.einsum("si,ij,sj->s", n, A, n)
    lw -= lw.max()
    w = np.exp(lw)
    w /= w.sum()
    out = {}
    for k in range(1, 7):
        F = n @ np.exp(2j * math.pi * k * np.arange(Q) / Q)
        out[k] = float(w @ (np.abs(F) ** 2) / N)
    return out
