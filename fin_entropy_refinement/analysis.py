"""ST8547: independently checked entropy/refinement rigidity and correlation defect.

Natural logs; positive mass coordinates; row-sum-zero Laplacians.
This is a conditional mathematical construction, not a physical FIN source.
Run: python3 analysis.py --output results.json
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm


def strict_weights():
    n = 12
    W = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = min(j - i, n - j + i)
            W[i, j] = W[j, i] = math.cos(0.18575*d + 0.16250)/(1+d**(9/5))
    return W


def laplacian(W):
    return np.diag(W.sum(axis=1)) - W


def log_mean(a, b):
    """Cancellation-safe logarithmic mean, including its continuous diagonal."""
    if a <= 0 or b <= 0:
        raise ValueError("The implemented differential identity is interior-only.")
    h, ell = max(a, b), min(a, b)
    z = (h - ell) / ell
    if z == 0:
        return a
    if z < 0.1:
        return ell * z / math.log1p(z)
    return (h - ell) / (math.log(h) - math.log(ell))


OMEGA = 2*math.pi/math.log(2)
EPSILON = 0.25


def periodic_phi(x, eps=EPSILON):
    t = OMEGA*np.log(x)
    return x*np.log(x) - eps*x*(np.cos(t)+OMEGA*np.sin(t))/(OMEGA*(1+OMEGA**2))


def periodic_dphi(x, eps=EPSILON):
    return np.log(x)+1-eps*np.cos(OMEGA*np.log(x))/OMEGA


def periodic_ddphi(x, eps=EPSILON):
    return (1+eps*np.sin(OMEGA*np.log(x)))/x


def periodic_mean(a, b, eps=EPSILON):
    if a <= 0 or b <= 0:
        raise ValueError("Positive masses required")
    if a == b:
        return 1/float(periodic_ddphi(a, eps))
    # Average curvature over the log interval, evaluated without subtracting cosines.
    la, lb = math.log(a), math.log(b)
    dl = la-lb
    if abs(dl) < 1e-7:
        mid = (la+lb)/2
        # (cos(omega la)-cos(omega lb))/dl using sin difference.
        correction = 2*eps*math.sin(OMEGA*mid)*math.sin(OMEGA*dl/2)/(OMEGA*dl)
        return log_mean(a, b)/(1+correction)
    denom = dl-eps*(math.cos(OMEGA*la)-math.cos(OMEGA*lb))/OMEGA
    return (a-b)/denom


def mobility(W, p, mean=log_mean):
    p = np.asarray(p, dtype=float)
    if np.any(p <= 0):
        raise ValueError("Mobility requires interior state")
    n = len(p)
    if W.shape != (n, n):
        raise ValueError("State and conductance dimensions disagree")
    C = np.zeros_like(W, dtype=float)
    for i in range(n):
        for j in range(i+1, n):
            C[i, j] = C[j, i] = W[i, j]*mean(p[i], p[j])
    return laplacian(C)


def product_weights(W, B):
    return np.kron(W, np.eye(len(B))) + np.kron(np.eye(len(W)), B)


def collapse(n, m):
    """Sums probabilities over fibers. This is NOT the Hilbert isometry R."""
    return np.kron(np.eye(n), np.ones((1, m)))


def dissipation(W, p):
    p = np.asarray(p)
    return sum(W[i, j]*(p[i]-p[j])*(math.log(p[i])-math.log(p[j]))
               for i in range(len(p)) for j in range(i+1, len(p)))


def kl(a, b):
    return float(np.dot(a, np.log(a/b)))


def hidden_horizontal_dissipation(W, P):
    """Exact analytic RHS of the horizontal entropy-production gap identity."""
    p = P.sum(axis=1)
    cond = P/p[:, None]
    return sum(W[i, j]*(p[i]*kl(cond[i], cond[j])+p[j]*kl(cond[j], cond[i]))
               for i in range(len(p)) for j in range(i+1, len(p)))


def hidden_horizontal_direct(W, P):
    p = P.sum(axis=1)
    return sum(dissipation(W, P[:, b]) for b in range(P.shape[1]))-dissipation(W, p)


def correlation_case(W, p, epsilon):
    v = np.cos(2*np.pi*np.arange(len(p))/len(p))
    v -= p@v
    v /= np.max(np.abs(v))
    cond = np.column_stack(((1+epsilon*v)/2, (1-epsilon*v)/2))
    return p[:, None]*cond


def projected_defect(W, P, B):
    C = collapse(len(W), len(B))
    return mobility(W, P.sum(axis=1))-C@mobility(product_weights(W, B), P.ravel())@C.T


def artifact_audit(repo):
    """Count payload-only JSON in the previous range, never mistake it for proof."""
    counts = {"program_records": 0, "empty_payload_files": 0}
    for lo in range(8397, 8547, 15):
        file = repo/f"FIN_ST{lo}_ST{lo+14}_Results.json"
        if not file.exists():
            continue
        results = json.loads(file.read_text())
        for row in results.values():
            counts["program_records"] += 1
            payload = json.loads((repo/row["packet_file"]).read_text())
            counts["empty_payload_files"] += (payload == {})
    return counts


def run():
    W = strict_weights()
    A = laplacian(W)
    n = len(W)
    p = np.arange(1, n+1, dtype=float)
    p /= p.sum()
    baseK = mobility(W, p)
    B = np.array([[0., .2], [.2, 0.]])
    Wf = product_weights(W, B)
    C = collapse(n, 2)
    rng = np.random.default_rng(8547)
    # These tests recompute quantities; they do not check stored pass flags.
    identity_errors = []
    pair_records = []
    for _ in range(80):
        pp = np.exp(rng.uniform(-7, 2, n))
        pp /= pp.sum()
        identity_errors.append(float(np.max(np.abs(mobility(W, pp)@np.log(pp)-A@pp))))
    for a, b in [(0.07, .19), (.11, .31), (.03, .12), (.2, .2)]:
        M = periodic_mean(a, b)
        pair_records.append({
            "a": a, "b": b,
            "binary_relative_defect": abs(2*periodic_mean(a/2, b/2)/M-1),
            "ternary_relative_defect": abs(3*periodic_mean(a/3, b/3)/M-1),
            "biased_0_3_relative_defect":
                abs((periodic_mean(.3*a, .3*b)+periodic_mean(.7*a, .7*b))/M-1),
        })
    product_errors = []
    for r in [np.array([.5, .5]), np.array([.3, .7]), np.array([.05, .95])]:
        P = p[:, None]*r[None, :]
        product_errors.append(float(np.linalg.norm(projected_defect(W, P, B), np.inf)))
    P = correlation_case(W, p, .6)
    D = projected_defect(W, P, B)
    lam = np.linalg.eigvalsh(D)
    q0, q1 = hidden_horizontal_direct(W, P), hidden_horizontal_dissipation(W, P)
    time = .4
    coarse_error = float(np.max(np.abs(C@expm(-time*laplacian(Wf))@P.ravel()-expm(-time*A)@p)))
    sweeps = []
    for eps in [0., .1, .3, .6, .9]:
        PP = correlation_case(W, p, eps)
        DD = projected_defect(W, PP, B)
        sweeps.append({"epsilon": eps, "defect_trace": float(np.trace(DD)),
                       "entropy_gap": hidden_horizontal_direct(W, PP)})
    # Direct counterexample to global metric uniqueness from its value on one gradient.
    v = np.arange(n, dtype=float)
    g = np.log(p)-np.mean(np.log(p))
    v -= v.mean()
    v -= g*(g@v)/(g@g)
    v /= np.linalg.norm(v)
    scale = .1*min(-baseK[i, j] for i in range(n) for j in range(i+1, n))
    altered = baseK+scale*np.outer(v, v)
    # Reference-density counterexample: quadratic relative entropy, not Shannon.
    pi = np.full(n, 1/n)
    r = np.array([.3, .7])
    cw = pi[:, None]*W
    relative_coarse = laplacian(cw)  # M_phi=1 for phi(u)=u^2/2
    relative_fine = laplacian(np.kron(cw, np.diag(r)))
    reference_error = float(np.linalg.norm(C@relative_fine@C.T-relative_coarse, np.inf))
    result = {
        "schema": "FIN-ENTROPY-REFINEMENT-RESEARCH-v1",
        "programme": "ST8547",
        "scientific_status": "Analytic conditional theorems with numerical checks; no physical source claim",
        "strict": {"s": float(A[0, 0]), "min_edge": float(W[0, 6]),
                   "max_heat_identity_error_80_states": max(identity_errors)},
        "dyadic_counterexample": {
            "epsilon": EPSILON, "omega_log": OMEGA,
            "curvature_lower_bound_x_phi_second": 1-EPSILON,
            "diagonal_mobility_second_derivative_at_one":
                2*EPSILON**2*OMEGA**2-EPSILON*OMEGA,
            "pairs": pair_records},
        "log_mean_refinement": {
            "product_defects_inf": product_errors,
            "correlated_defect_eigenvalues": lam.tolist(),
            "correlated_horizontal_entropy_gap": q0,
            "KL_formula_gap": q1, "identity_error": abs(q0-q1),
            "coarse_heat_identity_error": coarse_error, "sweep": sweeps},
        "nonlocal_metric_counterexample": {
            "operator_change_norm": float(np.linalg.norm(altered-baseK)),
            "heat_identity_error": float(np.max(np.abs(altered@np.log(p)-A@p))),
            "max_offdiagonal": float(np.max(altered-np.diag(np.diag(altered))))},
        "reference_density_quadratic_split_error": reference_error,
        "scope": {
            "Shannon_rigidity_requires": "heat matching + dyadic mobility compatibility + diagonal concavity, in mass coordinates",
            "additional_alternative": "binary plus ternary equal split compatibility",
            "binary_refinement_alone_selects_Shannon": False,
            "known_log_mean_heat_identity": "Maas 2011",
            "microscopic_rates_derived": False, "selector_sourced": False,
            "physical_units_sourced": False},
    }
    return result


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, default=Path("results.json"))
    args = ap.parse_args()
    result = run()
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print(json.dumps({k: result[k] for k in ("strict", "dyadic_counterexample", "log_mean_refinement")}, indent=2))
