#!/usr/bin/env python3
"""ST132: construct and certify the ST108 fold center without reading it.

The discovery stage starts from a coarse localized ansatz written below and
uses damped Newton iteration.  The proof stage independently rebuilds the
transcendental strict-kernel intervals and applies the same interval Krawczyk
test as ST120.  The prior ST108/ST120 certificate is never opened.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "FIN_ST132_Independent_Center_Isolation_Certificate.json"
N = 12


def iv(value):
    if isinstance(value, tuple):
        return mp.iv.mpf([str(value[0]), str(value[1])])
    if isinstance(value, (float, np.floating)):
        return mp.iv.mpf([
            str(float(np.nextafter(value, -np.inf))),
            str(float(np.nextafter(value, np.inf))),
        ])
    return mp.iv.mpf(str(value))


def bounds(x):
    return (
        float(np.nextafter(float(x.a), -np.inf)),
        float(np.nextafter(float(x.b), np.inf)),
    )


def prod(a0, a1, b0, b1):
    vals = [a0*b0, a0*b1, a1*b0, a1*b1]
    return float(np.nextafter(min(vals), -np.inf)), float(np.nextafter(max(vals), np.inf))


def matvec(alo, ahi, xlo, xhi):
    lo = np.zeros(alo.shape[0]); hi = np.zeros(alo.shape[0])
    for i in range(alo.shape[0]):
        for j in range(alo.shape[1]):
            q0, q1 = prod(alo[i, j], ahi[i, j], xlo[j], xhi[j])
            lo[i] = np.nextafter(lo[i] + q0, -np.inf)
            hi[i] = np.nextafter(hi[i] + q1, np.inf)
    return lo, hi


def left_product(r, alo, ahi):
    lo = np.zeros((r.shape[0], alo.shape[1])); hi = np.zeros_like(lo)
    for i in range(r.shape[0]):
        for j in range(alo.shape[1]):
            for k in range(r.shape[1]):
                q0, q1 = prod(r[i, k], r[i, k], alo[k, j], ahi[k, j])
                lo[i, j] = np.nextafter(lo[i, j] + q0, -np.inf)
                hi[i, j] = np.nextafter(hi[i, j] + q1, np.inf)
    return lo, hi


def strict_float_matrix():
    omega, phi, eta = 0.18575, 0.16250, 9/5
    weights = {d: math.cos(omega*d + phi)/(1+d**eta) for d in range(1, 7)}
    s = 2*sum(weights[d] for d in range(1, 6)) + weights[6]
    a = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            a[i, j] = s if i == j else -weights[min((i-j) % N, (j-i) % N)]
    return a


def strict_interval_matrix():
    omega, phi, eta = iv("0.18575"), iv("0.16250"), iv(9)/iv(5)
    weights = {d: mp.iv.cos(omega*d + phi)/(1+iv(d)**eta) for d in range(1, 7)}
    s = 2*sum((weights[d] for d in range(1, 6)), iv(0)) + weights[6]
    a = [[iv(0) for _ in range(N)] for _ in range(N)]
    for i in range(N):
        for j in range(N):
            a[i][j] = s if i == j else -weights[min((i-j) % N, (j-i) % N)]
    return a, weights, s


def float_system_jacobian(x, a):
    q, kappa, v = x[:7], x[7], x[8:]
    idx = np.array([i if i <= 6 else N-i for i in range(N)])
    u = q[idx]; w = v[idx]
    au = a @ u; aw = a @ w
    h = np.zeros(N); rdiag = np.zeros(N); drdu = np.zeros(N)
    for i, item in enumerate(u):
        rho = item*item; den = 1 + rho/2
        qfun = rho/den; qp = den**-2; qpp = -den**-3; qppp = 1.5*den**-4
        hh = -qfun*qp + .075; hp = -(qp*qp + qfun*qpp); hpp = -(3*qp*qpp + qfun*qppp)
        h[i] = hh; rdiag[i] = 2*hh + 4*rho*hp; drdu[i] = 2*item*(6*hp + 4*rho*hpp)
    g = np.r_[kappa*au[:7] + 2*u[:7]*h[:7], kappa*aw[:7] + rdiag[:7]*w[:7], .5*(v@v-1)]
    jac = np.zeros((15, 15))
    for i in range(7):
        for col in range(7):
            total = kappa*sum(a[i, site] for site in range(N) if idx[site] == col)
            if i == col: total += rdiag[i]
            jac[i, col] = total; jac[7+i, 8+col] = total
        jac[i, 7] = au[i]; jac[7+i, i] = drdu[i]*w[i]
        jac[7+i, 7] = aw[i]; jac[14, 8+i] = v[i]
    return g, jac


def discover_center():
    # Qualitative localized seed: one high central amplitude and one dominant
    # odd tangent component.  It contains no digits from ST108/ST120 output.
    seed = np.array([
        2.5, .2, .08, .03, .015, .008, .004, .02,
        0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ])
    a = strict_float_matrix(); x = seed.copy(); history = []
    for iteration in range(80):
        g, jac = float_system_jacobian(x, a); norm = float(np.linalg.norm(g, ord=np.inf))
        history.append({"iteration": iteration, "residual_inf": norm})
        if norm < 2e-14: break
        step = np.linalg.solve(jac, -g)
        damping = 1.0
        while damping > 2**-20:
            trial = x + damping*step
            trial_norm = float(np.linalg.norm(float_system_jacobian(trial, a)[0], ord=np.inf))
            if trial_norm < norm: x = trial; break
            damping /= 2
        else:
            raise RuntimeError("Damped Newton failed to decrease the residual")
    return seed, x, history


def interval_system(center, radius, a):
    q = [iv((center[i]-radius, center[i]+radius)) for i in range(7)]
    kappa = iv((center[7]-radius, center[7]+radius))
    v = [iv((center[8+i]-radius, center[8+i]+radius)) for i in range(7)]
    u, w = [], []
    for site in range(N):
        index = site if site <= 6 else N-site
        u.append(q[index]); w.append(v[index])
    au = [sum((a[i][j]*u[j] for j in range(N)), iv(0)) for i in range(N)]
    aw = [sum((a[i][j]*w[j] for j in range(N)), iv(0)) for i in range(N)]
    h, rdiag, drdu = [], [], []
    for item in u:
        rho = item**2; den = 1+rho/2; qfun = rho/den
        qp = den**-2; qpp = -den**-3; qppp = iv("1.5")*den**-4
        hh = -qfun*qp+iv("0.075"); hp = -(qp**2+qfun*qpp); hpp = -(3*qp*qpp+qfun*qppp)
        h.append(hh); rdiag.append(2*hh+4*rho*hp); drdu.append(2*item*(6*hp+4*rho*hpp))
    g = [kappa*au[i]+2*u[i]*h[i] for i in range(7)]
    g += [kappa*aw[i]+rdiag[i]*w[i] for i in range(7)]
    g += [iv("0.5")*(sum((item**2 for item in v), iv(0))-1)]
    jac = [[iv(0) for _ in range(15)] for _ in range(15)]
    for i in range(7):
        for col in range(7):
            total = iv(0)
            for site in range(N):
                if (site if site <= 6 else N-site) == col: total += kappa*a[i][site]
            if i == col: total += rdiag[i]
            jac[i][col] = total; jac[7+i][8+col] = total
        jac[i][7] = au[i]; jac[7+i][i] = drdu[i]*w[i]
        jac[7+i][7] = aw[i]; jac[14][8+i] = v[i]
    glo = np.array([bounds(x)[0] for x in g]); ghi = np.array([bounds(x)[1] for x in g])
    jlo = np.array([[bounds(x)[0] for x in row] for row in jac])
    jhi = np.array([[bounds(x)[1] for x in row] for row in jac])
    return glo, ghi, jlo, jhi


def run(write=True):
    mp.iv.dps = 70
    seed, center, history = discover_center()
    a, weights, s = strict_interval_matrix()
    g0lo, g0hi, j0lo, j0hi = interval_system(center, 0.0, a)
    preconditioner = np.linalg.inv(.5*(j0lo+j0hi))
    attempts = []; accepted = None
    for radius in [1e-8, 3e-9, 1e-9, 3e-10]:
        _, _, jlo, jhi = interval_system(center, radius, a)
        cglo, cghi = matvec(preconditioner, preconditioner, g0lo, g0hi)
        ylo = np.nextafter(center-cghi, -np.inf); yhi = np.nextafter(center-cglo, np.inf)
        cjlo, cjhi = left_product(preconditioner, jlo, jhi); mlo, mhi = -cjhi, -cjlo
        for i in range(15):
            mlo[i, i] = np.nextafter(mlo[i, i]+1, -np.inf)
            mhi[i, i] = np.nextafter(mhi[i, i]+1, np.inf)
        mdlo, mdhi = matvec(mlo, mhi, np.full(15, -radius), np.full(15, radius))
        klo = np.nextafter(ylo+mdlo, -np.inf); khi = np.nextafter(yhi+mdhi, np.inf)
        margin = float(min(np.min(klo-(center-radius)), np.min((center+radius)-khi)))
        row = {"radius": radius, "minimum_strict_inclusion_margin": margin,
               "maximum_Krawczyk_half_width": float(np.max((khi-klo)/2)), "included": margin > 0}
        attempts.append(row)
        if accepted is None and margin > 0: accepted = row
    result = {
        "prior_certificate_opened": False,
        "dependencies": ["Python standard library", "NumPy", "mpmath"],
        "coarse_seed": seed.tolist(), "discovered_center": center.tolist(),
        "Newton_iterations": len(history)-1, "Newton_history": history,
        "final_point_residual_inf": float(np.linalg.norm(float_system_jacobian(center, strict_float_matrix())[0], ord=np.inf)),
        "attempts": attempts, "accepted": accepted,
        "row_sum_interval": list(bounds(s)),
        "kernel_widths": {str(d): bounds(weights[d])[1]-bounds(weights[d])[0] for d in weights},
        "scope": (
            "The center is generated from a coarse localized seed without opening ST108/ST120. "
            "Krawczyk inclusion proves one unique root in the accepted local box for every strict operator in the declared transcendental interval family."
        ),
        "boundary": "This is local basin isolation, not an exhaustive proof that no other fold root exists outside the certified box.",
    }
    if write: OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    return result


if __name__ == "__main__":
    print(json.dumps(run(True), indent=2, sort_keys=True))
