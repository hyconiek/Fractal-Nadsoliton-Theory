#!/usr/bin/env python3
"""
QW-1741: Constrained global derivation of (beta_tors, omega, phi) from signed dynamics.

Key upgrades vs QW-1739:
1) Multi-run global fit (shared omega,phi,beta; per-run amplitudes).
2) Robust objective + physically motivated constraints (not canonical priors).
3) Bootstrap uncertainty from run-resampling.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1741_constrained_global_derivation.json"
OUT_MD = ROOT / "RAPORT_QW1741_CONSTRAINED_GLOBAL_DERIVATION.md"


@dataclass
class SimCfg:
    n_nodes: int
    seed: int
    xi_local: float
    p_short: float
    eta_short: float
    tau_torsion: float
    rho_antisym: float
    dyn_c2: float
    dyn_damp: float
    n_steps: int = 280
    burn: int = 120


def ring_dist(i: int, j: int, n: int) -> int:
    d = abs(i - j)
    return min(d, n - d)


def smooth_periodic(x: np.ndarray, iters: int = 5) -> np.ndarray:
    y = x.copy()
    for _ in range(iters):
        y = 0.25 * np.roll(y, -1) + 0.5 * y + 0.25 * np.roll(y, 1)
    return y


def micro_state(n: int, rng: np.random.Generator) -> Tuple[np.ndarray, np.ndarray]:
    raw = rng.normal(size=n)
    th = math.pi * np.tanh(smooth_periodic(raw, iters=8) / max(np.std(raw), 1e-6))
    q = rng.integers(-2, 3, size=n).astype(float)
    return th, q


def build_signed_matrix(cfg: SimCfg, theta: np.ndarray, q: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    n = cfg.n_nodes
    w = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(i + 1, n):
            d = ring_dist(i, j, n)
            amp = math.exp(-d / cfg.xi_local)
            if rng.random() < cfg.p_short / (d ** cfg.eta_short):
                amp += 0.6 * (d ** -0.85) * (1.0 + 0.12 * rng.normal())

            dt = abs(theta[i] - theta[j]) / math.pi
            dq = abs(q[i] - q[j])
            tors = dq + 0.3 * dt

            sym_sign = 1.0 if (math.cos(theta[i] - theta[j]) + 0.18 * (q[i] - q[j])) >= 0 else -1.0
            sym = sym_sign * amp * math.exp(-cfg.tau_torsion * tors)
            anti = cfg.rho_antisym * math.sin(theta[i] - theta[j]) / (d ** 1.2)
            w[i, j] = sym + anti
            w[j, i] = sym - anti

    smax = float(np.linalg.norm(w, 2))
    if smax > 1e-12:
        w = w / smax
    return w


def simulate_profile(cfg: SimCfg) -> np.ndarray:
    rng = np.random.default_rng(cfg.seed)
    theta, q = micro_state(cfg.n_nodes, rng)
    w = build_signed_matrix(cfg, theta, q, rng)

    n = cfg.n_nodes
    x_prev = 0.03 * rng.normal(size=n)
    x_curr = 0.02 * rng.normal(size=n)
    x_curr[0] += 1.0

    lap = w - np.eye(n)
    snaps: List[np.ndarray] = []
    for _ in range(cfg.n_steps):
        x_next = 2.0 * x_curr - x_prev + cfg.dyn_c2 * (lap @ x_curr) - cfg.dyn_damp * (x_curr - x_prev)
        x_next = np.tanh(1.3 * x_next)
        x_prev, x_curr = x_curr, x_next
        snaps.append(x_curr.copy())

    arr = np.array(snaps[cfg.burn :], dtype=float)
    dmax = min(12, n // 2)
    prof = []
    for d in range(1, dmax + 1):
        c = arr * np.roll(arr, -d, axis=1)
        prof.append(float(np.mean(c)))
    y = np.array(prof, dtype=float)
    if abs(y[0]) > 1e-12:
        y = y / abs(y[0])
    return y


def robust_loss(res: np.ndarray, delta: float = 0.08) -> np.ndarray:
    a = np.abs(res)
    quad = 0.5 * (a ** 2)
    lin = delta * (a - 0.5 * delta)
    return np.where(a <= delta, quad, lin)


def predict_basis(d: np.ndarray, omega: float, phi: float, beta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * d)


def global_objective(theta: Tuple[float, float, float], y: np.ndarray, weights: np.ndarray) -> float:
    omega, phi, beta = theta
    if not (0.15 <= omega <= 1.6 and -math.pi <= phi <= math.pi and 0.001 <= beta <= 0.25):
        return float("inf")
    d = np.arange(1, y.shape[1] + 1, dtype=float)
    b = predict_basis(d, omega, phi, beta)
    bb = float(np.dot(b, b))
    if bb <= 1e-12:
        return float("inf")

    # per-run nuisance amplitudes
    a = (y @ b) / bb
    pred = a[:, None] * b[None, :]
    res = y - pred
    core = np.sum(weights[:, None] * robust_loss(res))

    # Physical regularizers (weak, non-canonical):
    env = 1.0 / (1.0 + beta * d)
    # envelope should be non-increasing
    pen_env = float(np.sum(np.maximum(0.0, np.diff(env))))
    # prevent near-zero frequency collapse from prior failures
    pen_omega = float(np.maximum(0.0, 0.25 - omega) ** 2)
    # avoid extreme beta boundary behavior
    pen_beta = float(np.maximum(0.0, beta - 0.18) ** 2)

    return float(core + 8.0 * pen_env + 2.5 * pen_omega + 4.0 * pen_beta)


def refine(theta0: Tuple[float, float, float], y: np.ndarray, w: np.ndarray) -> Tuple[Tuple[float, float, float], float]:
    cur = (float(theta0[0]), float(theta0[1]), float(theta0[2]))
    fcur = global_objective(cur, y, w)

    steps = [
        (0.20, 0.80, 0.05),
        (0.08, 0.30, 0.02),
        (0.03, 0.12, 0.008),
        (0.012, 0.05, 0.003),
    ]
    for so, sp, sb in steps:
        improved = True
        while improved:
            improved = False

            # omega
            cand = np.linspace(max(0.15, cur[0] - so), min(1.6, cur[0] + so), 9)
            best = (cur, fcur)
            for om in cand:
                th = (float(om), cur[1], cur[2])
                f = global_objective(th, y, w)
                if f < best[1]:
                    best = (th, f)
            if best[1] < fcur:
                cur, fcur = best
                improved = True

            # phi
            cand = np.linspace(cur[1] - sp, cur[1] + sp, 9)
            best = (cur, fcur)
            for ph in cand:
                ph = float((ph + math.pi) % (2.0 * math.pi) - math.pi)
                th = (cur[0], ph, cur[2])
                f = global_objective(th, y, w)
                if f < best[1]:
                    best = (th, f)
            if best[1] < fcur:
                cur, fcur = best
                improved = True

            # beta
            cand = np.linspace(max(0.001, cur[2] - sb), min(0.25, cur[2] + sb), 9)
            best = (cur, fcur)
            for be in cand:
                th = (cur[0], cur[1], float(be))
                f = global_objective(th, y, w)
                if f < best[1]:
                    best = (th, f)
            if best[1] < fcur:
                cur, fcur = best
                improved = True

    return cur, fcur


def bootstrap_fit(y: np.ndarray, w: np.ndarray, theta_init: Tuple[float, float, float], n_boot: int = 120, seed: int = 1741) -> List[Tuple[float, float, float]]:
    rng = np.random.default_rng(seed)
    n = y.shape[0]
    out = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yb = y[idx]
        wb = w[idx]
        th, _ = refine(theta_init, yb, wb)
        out.append(th)
    return out


def circular_mean(x: np.ndarray) -> float:
    return float(np.arctan2(np.mean(np.sin(x)), np.mean(np.cos(x))))


def circular_std(x: np.ndarray) -> float:
    r = np.sqrt(np.mean(np.sin(x)) ** 2 + np.mean(np.cos(x)) ** 2)
    if r <= 1e-12:
        return float("inf")
    return float(np.sqrt(-2.0 * np.log(r)))


def main() -> None:
    n_grid = [72, 96, 120]
    seeds_per_n = 10

    cfgs: List[SimCfg] = []
    for n in n_grid:
        for k in range(seeds_per_n):
            seed = 174100 + 100 * n + k
            cfgs.append(
                SimCfg(
                    n_nodes=n,
                    seed=seed,
                    xi_local=1.7 + 0.15 * ((k % 3) - 1),
                    p_short=0.08 + 0.02 * (k % 2),
                    eta_short=1.55 + 0.20 * (k % 3),
                    tau_torsion=0.55 + 0.12 * ((k + 1) % 3),
                    rho_antisym=0.14 + 0.06 * (k % 2),
                    dyn_c2=0.30 + 0.05 * (k % 3),
                    dyn_damp=0.045 + 0.015 * (k % 2),
                )
            )

    profiles = np.array([simulate_profile(c) for c in cfgs], dtype=float)
    var = np.var(profiles, axis=1)
    weights = 1.0 / np.clip(var, 1e-5, None)

    starts = [
        (0.35, 0.2, 0.02),
        (0.55, 0.8, 0.04),
        (0.75, -0.3, 0.06),
        (0.95, 0.5, 0.08),
        (1.15, -0.9, 0.03),
        (0.45, -1.2, 0.12),
    ]

    sols = [refine(st, profiles, weights) for st in starts]
    sols = sorted(sols, key=lambda x: x[1])
    theta_best, f_best = sols[0]

    # Derived run-wise amplitudes and fit metrics.
    d = np.arange(1, profiles.shape[1] + 1, dtype=float)
    b = predict_basis(d, theta_best[0], theta_best[1], theta_best[2])
    bb = float(np.dot(b, b))
    a = (profiles @ b) / bb
    pred = a[:, None] * b[None, :]
    rmse_runs = np.sqrt(np.mean((profiles - pred) ** 2, axis=1))
    ss_tot = np.sum((profiles - np.mean(profiles, axis=1, keepdims=True)) ** 2, axis=1)
    ss_res = np.sum((profiles - pred) ** 2, axis=1)
    r2_runs = np.where(ss_tot > 1e-12, 1.0 - ss_res / ss_tot, 1.0)

    # Bootstrap uncertainty.
    boots = bootstrap_fit(profiles, weights, theta_best, n_boot=120, seed=1741)
    bo = np.array([x[0] for x in boots], dtype=float)
    bp = np.array([x[1] for x in boots], dtype=float)
    bbeta = np.array([x[2] for x in boots], dtype=float)

    omega_med = float(np.median(bo))
    omega_ci = [float(np.quantile(bo, 0.025)), float(np.quantile(bo, 0.975))]
    phi_mean = circular_mean(bp)
    phi_std = circular_std(bp)
    phi_ci = [float(np.quantile(bp, 0.025)), float(np.quantile(bp, 0.975))]
    beta_med = float(np.median(bbeta))
    beta_ci = [float(np.quantile(bbeta, 0.025)), float(np.quantile(bbeta, 0.975))]

    # Not used in fitting; diagnostic only.
    omega_ref = math.pi / 4.0
    phi_ref = math.pi / 6.0
    beta_ref = 0.01
    dev = {
        "omega_rel_vs_pi_over_4": abs(omega_med - omega_ref) / omega_ref,
        "phi_abs_over_pi_vs_pi_over_6": abs(((phi_mean - phi_ref + math.pi) % (2 * math.pi)) - math.pi) / math.pi,
        "beta_rel_vs_001": abs(beta_med - beta_ref) / beta_ref,
    }

    pass_fit = float(np.median(r2_runs)) >= 0.70
    pass_ident = (omega_ci[1] - omega_ci[0]) <= 0.22 and (beta_ci[1] - beta_ci[0]) <= 0.08 and phi_std <= 0.8
    pass_nonboundary = 0.005 <= beta_med <= 0.18 and 0.25 <= omega_med <= 1.4

    if pass_fit and pass_ident and pass_nonboundary:
        verdict = "CONSTRAINED_GLOBAL_DERIVATION_IDENTIFIABLE"
    elif pass_fit and pass_nonboundary:
        verdict = "CONSTRAINED_GLOBAL_DERIVATION_PARTIAL"
    else:
        verdict = "CONSTRAINED_GLOBAL_DERIVATION_NOT_CLOSED"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "n_grid": n_grid,
            "seeds_per_n": seeds_per_n,
            "n_total_runs": len(cfgs),
            "n_bootstrap": 120,
            "objective": "robust global fit with weak physical constraints",
        },
        "optimum": {
            "omega": float(theta_best[0]),
            "phi": float(theta_best[1]),
            "beta": float(theta_best[2]),
            "objective_value": float(f_best),
        },
        "run_metrics": {
            "median_r2": float(np.median(r2_runs)),
            "q10_r2": float(np.quantile(r2_runs, 0.1)),
            "median_rmse": float(np.median(rmse_runs)),
        },
        "bootstrap_summary": {
            "omega_median": omega_med,
            "omega_ci95": omega_ci,
            "phi_circular_mean": phi_mean,
            "phi_circular_std": phi_std,
            "phi_quantile_ci95_raw": phi_ci,
            "beta_median": beta_med,
            "beta_ci95": beta_ci,
        },
        "reference_deviation": dev,
        "pass_flags": {
            "fit_quality": pass_fit,
            "identifiability": pass_ident,
            "nonboundary_solution": pass_nonboundary,
        },
        "verdict": verdict,
        "solutions_from_starts": [
            {"theta": [float(v) for v in th], "objective": float(f)} for th, f in sols
        ],
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1741: CONSTRAINED GLOBAL DERIVATION",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Runs: {len(cfgs)}",
        f"- Optimum: omega={theta_best[0]:.6f}, phi={theta_best[1]:.6f}, beta={theta_best[2]:.6f}",
        f"- Run median R2={np.median(r2_runs):.4f}, median RMSE={np.median(rmse_runs):.4f}",
        f"- Bootstrap medians: omega={omega_med:.6f}, phi_mean={phi_mean:.6f}, beta={beta_med:.6f}",
        f"- Dev vs refs: domega_rel={dev['omega_rel_vs_pi_over_4']:.3f}, dphi/pi={dev['phi_abs_over_pi_vs_pi_over_6']:.3f}, dbeta_rel={dev['beta_rel_vs_001']:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- fit_quality: {pass_fit}",
        f"- identifiability: {pass_ident}",
        f"- nonboundary_solution: {pass_nonboundary}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1741] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1741] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
