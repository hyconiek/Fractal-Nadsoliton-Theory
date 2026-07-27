#!/usr/bin/env python3
"""Corrected post-Programs-41--50 audit and mirror-coupling tests.

All objects are finite and dimensionless.  The script does not export physical
units, a strict selector, a completed legacy-to-strict bridge, role transfer,
or a Theory-of-Everything claim.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "FIN_Post_41_50_Corrected_Results.json"
FIG = ROOT / "FIN_Post_41_50_Correction_Figures"
FIG.mkdir(exist_ok=True)

N = 12
ALPHA_GEO = 4.0 * math.log(2.0)
OMEGA_L = math.pi / 4.0
PHI_L = math.pi / 6.0
BETA_TORS = 0.01
OMEGA_S = 0.18575
PHI_S = 0.16250
BETA_S = 1.0
ETA_S = 1.8


def k_family(d, amplitude: float, beta: float):
    d = np.asarray(d, dtype=float)
    return amplitude * np.cos(OMEGA_L * d + PHI_L) / (1.0 + beta * d)


def k_legacy(d):
    return k_family(d, ALPHA_GEO, BETA_TORS)


def k_strict(d):
    d = np.asarray(d, dtype=float)
    return np.cos(OMEGA_S * d + PHI_S) / (1.0 + BETA_S * d**ETA_S)


def cycle_distance(i: int, j: int, n: int) -> int:
    return min((i - j) % n, (j - i) % n)


def reflection(n: int) -> np.ndarray:
    r = np.zeros((n, n))
    for i in range(n):
        r[i, (-i) % n] = 1.0
    return r


def undirected_matrix(kernel, n: int, diagonal: bool = True) -> np.ndarray:
    w = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            if not diagonal and i == j:
                continue
            w[i, j] = float(kernel(cycle_distance(i, j, n)))
    return w


def directed_matrix(kernel, n: int, diagonal: bool = False) -> np.ndarray:
    w = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            if not diagonal and i == j:
                continue
            w[i, j] = float(kernel((j - i) % n))
    return w


def laplacian(weights: np.ndarray) -> np.ndarray:
    return np.diag(weights.sum(axis=1)) - weights


def beta_threshold(n: int) -> dict:
    def f(beta: float) -> float:
        w = undirected_matrix(lambda d: k_family(d, 1.0, beta), n, diagonal=True)
        return float(np.linalg.eigvalsh(w).min())

    grid = np.logspace(-4, 3, 3000)
    values = [f(float(x)) for x in grid]
    bracket = None
    for left, right, fl, fr in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
        if fl < 0.0 <= fr:
            bracket = [float(left), float(right)]
            break
    if bracket is None:
        return {"n": n, "threshold": None}
    left, right = bracket
    for _ in range(100):
        mid = 0.5 * (left + right)
        if f(mid) < 0.0:
            left = mid
        else:
            right = mid
    root = 0.5 * (left + right)
    return {"n": n, "threshold": root, "root_residual": f(root)}


def metric_audit(beta: float, normalized: bool) -> dict:
    w0 = abs(float(k_family(0, 1.0, beta)))
    dmat = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            value = abs(float(k_family(cycle_distance(i, j, N), 1.0, beta)))
            ratio = value / w0 if normalized else value
            dmat[i, j] = -math.log(max(ratio, 1e-300))
    violations = []
    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            for k in range(N):
                if k in (i, j):
                    continue
                excess = dmat[i, k] - dmat[i, j] - dmat[j, k]
                if excess > 1e-12:
                    violations.append(float(excess))
    diagonal_zero = bool(np.max(np.abs(np.diag(dmat))) <= 1e-12)
    nonnegative = bool(dmat.min() >= -1e-12)
    separates = bool(
        all(dmat[i, j] > 1e-12 for i in range(N) for j in range(N) if i != j)
    )
    return {
        "beta": beta,
        "normalized_by_K0": normalized,
        "diagonal_zero": diagonal_zero,
        "nonnegative": nonnegative,
        "separates_points": separates,
        "triangle_violations": len(violations),
        "max_triangle_excess": max(violations) if violations else 0.0,
        "is_metric": bool(
            diagonal_zero and nonnegative and separates and not violations
        ),
    }


def corrected_hazard_fit() -> dict:
    d = np.arange(1.0, 7.0)
    star = np.abs(k_family(d, 2.9, 0.05))
    target = np.abs(k_strict(d))
    train = np.arange(4)

    def objective(x):
        c, eta = x
        pred = star[train] * np.exp(-c * d[train] ** eta)
        return float(np.sum((pred - target[train]) ** 2))

    fit = minimize(
        objective,
        np.array([0.47, 1.9]),
        method="L-BFGS-B",
        bounds=[(1e-9, 10.0), (0.1, 5.0)],
        options={"ftol": 1e-15, "gtol": 1e-12, "maxiter": 10000},
    )
    c, eta = map(float, fit.x)
    pred = star * np.exp(-c * d**eta)
    return {
        "training_distances": [1, 2, 3, 4],
        "held_out_distances": [5, 6],
        "c": c,
        "eta": eta,
        "train_relative_l2": float(
            np.linalg.norm(pred[:4] - target[:4]) / np.linalg.norm(target[:4])
        ),
        "holdout_relative_l2": float(
            np.linalg.norm(pred[4:] - target[4:]) / np.linalg.norm(target[4:])
        ),
        "optimizer_success": bool(fit.success),
    }


def phase_tautology() -> dict:
    a = OMEGA_S / OMEGA_L
    b = PHI_S - a * PHI_L
    d = np.arange(1.0, 49.0)
    residual = np.linalg.norm(
        a * (OMEGA_L * d + PHI_L) + b - (OMEGA_S * d + PHI_S)
    )
    return {
        "a_exact_definition": "omega_strict/omega_legacy",
        "a": a,
        "b_exact_definition": "phi_strict-a*phi_legacy",
        "b": b,
        "residual_d1_48": float(residual),
        "target_inserted": True,
        "predictive_content": False,
    }


def mirror_quantum_audit() -> dict:
    r = reflection(N)
    h0 = undirected_matrix(k_legacy, N, diagonal=False)
    directed = directed_matrix(k_legacy, N, diagonal=False)
    c_raw = 0.5j * (directed - directed.T)
    c_norm = c_raw / np.linalg.norm(c_raw, "fro")
    lam = 0.4
    hp = h0 + lam * c_norm
    hm = h0 - lam * c_norm
    symmetry_residual = float(np.linalg.norm(r @ h0 @ r.T - h0, "fro"))
    odd_residual = float(np.linalg.norm(r @ c_norm @ r.T + c_norm, "fro"))
    covariance = float(np.linalg.norm(r @ hp @ r.T - hm, "fro"))
    isospectral = float(
        np.linalg.norm(np.linalg.eigvalsh(hp) - np.linalg.eigvalsh(hm))
    )
    commutator = float(np.linalg.norm(hp @ r - r @ hp, "fro"))

    psi = np.zeros(N, dtype=complex)
    psi[2] = 1.0
    t = 0.7
    pp = np.abs(expm(-1j * t * hp) @ psi) ** 2
    pm = np.abs(expm(-1j * t * hm) @ psi) ** 2
    pm_reflected_input = np.abs(expm(-1j * t * hm) @ (r @ psi)) ** 2

    eye = np.eye(N)
    g = 0.3
    block_plus = np.block([[h0, g * eye], [g * eye, h0]])
    block_minus = np.block([[h0, -g * eye], [-g * eye, h0]])
    gauge = np.block([[eye, np.zeros_like(eye)], [np.zeros_like(eye), -eye]])
    block_gauge_residual = float(
        np.linalg.norm(gauge @ block_plus @ gauge - block_minus, "fro")
    )
    sector_swap = np.block([[np.zeros_like(r), r], [r, np.zeros_like(r)]])
    block_mirror_commutator = float(
        np.linalg.norm(block_plus @ sector_swap - sector_swap @ block_plus, "fro")
    )
    return {
        "radial_legacy_inversion_residual": symmetry_residual,
        "directed_odd_carrier_fro_norm_before_normalization": float(
            np.linalg.norm(c_raw, "fro")
        ),
        "odd_carrier_inversion_residual": odd_residual,
        "fixed_branch_reflection_commutator": commutator,
        "R_Hplus_R_minus_Hminus_residual": covariance,
        "plus_minus_isospectral_residual": isospectral,
        "same_site_probability_difference": float(np.max(np.abs(pp - pm))),
        "mirror_related_probability_residual": float(
            np.max(np.abs(r @ pp - pm_reflected_input))
        ),
        "two_copy_sign_gauge_residual": block_gauge_residual,
        "two_copy_mirror_commutator": block_mirror_commutator,
        "conclusion": (
            "A passive mirror copy and a symmetric two-copy coupling preserve a Z2 "
            "symmetry. An inversion-odd Hermitian carrier C breaks reflection for "
            "fixed nonzero lambda, but reflection exchanges the isospectral +/- "
            "branches. K_legacy supplies neither the directed chart nor lambda's sign."
        ),
    }


def chiral_markov_audit(lam: float = 0.4) -> dict:
    base = np.abs(k_legacy(np.arange(1, 7, dtype=float)))
    q = np.zeros((N, N))
    for i in range(N):
        for distance in range(1, 6):
            q[i, (i + distance) % N] = base[distance - 1] * math.exp(lam)
            q[i, (i - distance) % N] = base[distance - 1] * math.exp(-lam)
        q[i, (i + 6) % N] = base[5]
    np.fill_diagonal(q, -q.sum(axis=1))
    qm = np.zeros_like(q)
    for i in range(N):
        for distance in range(1, 6):
            qm[i, (i + distance) % N] = base[distance - 1] * math.exp(-lam)
            qm[i, (i - distance) % N] = base[distance - 1] * math.exp(lam)
        qm[i, (i + 6) % N] = base[5]
    np.fill_diagonal(qm, -qm.sum(axis=1))
    r = reflection(N)
    uniform = np.full(N, 1.0 / N)
    stationary_residual = float(np.linalg.norm(uniform @ q))
    mirror_covariance = float(np.linalg.norm(r @ q @ r.T - qm, "fro"))
    j1 = float(
        uniform[0] * q[0, 1] - uniform[1] * q[1, 0]
    )
    entropy_production = 0.0
    for i in range(N):
        for j in range(i + 1, N):
            if q[i, j] > 0 and q[j, i] > 0:
                current = uniform[i] * q[i, j] - uniform[j] * q[j, i]
                entropy_production += current * math.log(q[i, j] / q[j, i])
    return {
        "lambda": lam,
        "uniform_stationary_residual": stationary_residual,
        "mirror_covariance_residual": mirror_covariance,
        "nearest_edge_stationary_current": j1,
        "entropy_production_rate": float(entropy_production),
        "current_is_odd_under_lambda": True,
        "entropy_production_is_even_under_lambda": True,
        "imports_orientation_and_lambda_sign": True,
    }


def kernel_policy() -> dict:
    d = np.arange(1.0, 13.0)
    canonical = k_legacy(d)
    star_same = k_family(d, ALPHA_GEO, BETA_TORS)
    hist = k_family(d, 2.9, 0.05)
    z12 = k_family(d, 1.0, 1.0)
    return {
        "Kstar_family_contains_canonical_legacy_max_error": float(
            np.max(np.abs(canonical - star_same))
        ),
        "canonical_legacy_parameters": {
            "A": ALPHA_GEO,
            "beta": BETA_TORS,
            "omega": "pi/4",
            "phi": "pi/6",
        },
        "historical_freeze_relative_l2_to_canonical": float(
            np.linalg.norm(hist - canonical) / np.linalg.norm(canonical)
        ),
        "Z12_freeze_relative_l2_to_canonical": float(
            np.linalg.norm(z12 - canonical) / np.linalg.norm(canonical)
        ),
        "decision": {
            "primary_strict_work": "K_strict_gate",
            "only_canonical_non_strict_bridge_kernel": "K_legacy_ont with A=4 ln 2 and beta_tors=0.01",
            "Kstar_status": "parameter family containing K_legacy_ont, not a third canonical kernel",
            "A2.9_beta0.05_status": "historical stress-test freeze only",
            "A1_beta1_status": "finite operator benchmark only; not canonical",
            "K_rejected_status": "retired",
            "mirror_C_status": "conditional sector/operator candidate, not a kernel replacement",
            "Markov_policy": (
                "Use strict for primary diffusion. For legacy sensitivity report "
                "both W+=max(W,0) and |W| because positivity repair is not unique."
            ),
        },
    }


def figures(results: dict) -> list[str]:
    d = np.arange(1.0, 13.0)
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    ax.plot(d, k_legacy(d), "o-", label="canonical legacy")
    ax.plot(d, k_family(d, 2.9, 0.05), "s--", label="A=2.9, beta=0.05 benchmark")
    ax.plot(d, k_family(d, 1.0, 1.0), "^--", label="A=1, beta=1 benchmark")
    ax.plot(d, k_strict(d), "d-", label="strict")
    ax.axhline(0, color="black", lw=0.7)
    ax.set_xlabel("distance d")
    ax.set_ylabel("kernel value")
    ax.set_title("Canonical kernel policy: one strict and one legacy bridge object")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    p1 = FIG / "kernel_policy_profiles.png"
    fig.savefig(p1, dpi=180)
    plt.close(fig)

    thresholds = results["psd_threshold_support_scan"]
    fig, ax = plt.subplots(figsize=(7.8, 4.4))
    ax.plot(
        [row["n"] for row in thresholds],
        [row["threshold"] for row in thresholds],
        "o-",
    )
    ax.axhline(BETA_TORS, color="#a61b1b", ls="--", label="canonical beta=0.01")
    ax.set_xlabel("cycle size N")
    ax.set_ylabel("Gram PSD threshold beta*")
    ax.set_title("The Gram-PSD threshold is support dependent")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    p2 = FIG / "psd_threshold_support_scan.png"
    fig.savefig(p2, dpi=180)
    plt.close(fig)

    lam = np.linspace(-1.0, 1.0, 201)
    base1 = abs(float(k_legacy(1)))
    current = base1 * (np.exp(lam) - np.exp(-lam)) / N
    entropy = 2.0 * N * current * lam
    fig, ax = plt.subplots(figsize=(7.8, 4.4))
    ax.plot(lam, current, label="oriented current (odd)")
    ax.plot(lam, entropy, label="entropy-production proxy (even)")
    ax.axvline(0, color="black", lw=0.7)
    ax.set_xlabel("externally signed mirror bias lambda")
    ax.set_title("Mirror bias creates an arrow but does not select its sign")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    p3 = FIG / "mirror_bias_current_entropy.png"
    fig.savefig(p3, dpi=180)
    plt.close(fig)
    return [str(p.relative_to(ROOT)) for p in (p1, p2, p3)]


def main() -> None:
    results = {
        "metadata": {
            "date": "2026-07-27",
            "scope": "corrected finite dimensionless audit after Programs 41-50",
            "no_firecrawl": True,
        },
        "kernel_policy": kernel_policy(),
        "psd_threshold_support_scan": [
            beta_threshold(n) for n in [8, 10, 12, 14, 16, 20, 24, 32]
        ],
        "metric_audits": [
            metric_audit(beta, normalized)
            for beta in [0.01, 1.0, 2.0]
            for normalized in [False, True]
        ],
        "phase_map_correction": phase_tautology(),
        "hazard_fit_correction": corrected_hazard_fit(),
        "mirror_quantum_audit": mirror_quantum_audit(),
        "mirror_markov_audit": chiral_markov_audit(),
        "methodology_corrections": {
            "program42": "target-inserted affine identity, not reconstruction",
            "program43": "training objective restricted to d=1..4",
            "program45": "chosen dilation witness, not a kernel-derived environment",
            "program46": "actual mirror-odd Hamiltonian coupling now used",
            "program48": "closed-loop work replaces open-spiral proxy",
            "program49": "two-semigroup discrimination, not full process tensor",
            "program50": "strict-generated synthetic self-recovery, not validation",
            "P3172": "full metric axioms and bisection PSD threshold",
        },
    }
    results["figures"] = figures(results)
    OUT.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(OUT), "figures": results["figures"]}, indent=2))


if __name__ == "__main__":
    main()
