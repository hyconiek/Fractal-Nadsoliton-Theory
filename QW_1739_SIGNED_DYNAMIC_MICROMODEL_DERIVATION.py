#!/usr/bin/env python3
"""
QW-1739: Signed dynamic micromodel -> derivation attempt for (beta_tors, omega, phi).

Design principles:
1) No direct ansatz fixing of omega, phi, beta.
2) Microscopic signed couplings: local, shortcut, torsion damping, antisymmetric drift.
3) Effective kernel profile emerges from dynamic path summation.
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
OUT_JSON = ROOT / "report_qw1739_signed_dynamic_micromodel_derivation.json"
OUT_MD = ROOT / "RAPORT_QW1739_SIGNED_DYNAMIC_MICROMODEL_DERIVATION.md"


@dataclass
class MicroCfg:
    n_nodes: int
    seed: int
    local_cut: int = 4
    xi_local: float = 1.75
    p_short: float = 0.12
    eta_short: float = 1.7
    mu_short: float = 0.85
    short_scale: float = 0.80
    tau_torsion: float = 0.65
    charge_bias: float = 0.22
    rho_antisym: float = 0.18
    dyn_lambda: float = 0.77
    dyn_kmax: int = 7


def smooth_periodic(x: np.ndarray, n_iter: int = 5) -> np.ndarray:
    y = x.copy()
    for _ in range(n_iter):
        y = 0.25 * np.roll(y, -1) + 0.5 * y + 0.25 * np.roll(y, 1)
    return y


def build_micro_state(cfg: MicroCfg, rng: np.random.Generator) -> Tuple[np.ndarray, np.ndarray]:
    theta_raw = rng.normal(0.0, 1.0, size=cfg.n_nodes)
    theta_sm = smooth_periodic(theta_raw, n_iter=7)
    theta = math.pi * np.tanh(theta_sm / np.std(theta_sm))
    q = rng.integers(-2, 3, size=cfg.n_nodes)  # topological charges in [-2,2]
    return theta, q.astype(float)


def ring_dist(i: int, j: int, n: int) -> int:
    d = abs(i - j)
    return min(d, n - d)


def build_signed_dynamic_matrix(cfg: MicroCfg, theta: np.ndarray, q: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    n = cfg.n_nodes
    w = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(i + 1, n):
            d = ring_dist(i, j, n)
            amp = 0.0

            if d <= cfg.local_cut:
                amp += math.exp(-d / cfg.xi_local)

            if rng.random() < cfg.p_short / (d ** cfg.eta_short):
                amp += cfg.short_scale * (d ** (-cfg.mu_short)) * (1.0 + 0.12 * rng.normal())

            if amp <= 0.0:
                continue

            dt = abs(theta[i] - theta[j]) / math.pi
            dq = abs(q[i] - q[j])
            tors = dq + 0.35 * dt

            signed_scalar = math.cos(theta[i] - theta[j]) + cfg.charge_bias * (q[i] - q[j])
            sgn = 1.0 if signed_scalar >= 0.0 else -1.0

            sym = sgn * amp * math.exp(-cfg.tau_torsion * tors)
            anti = cfg.rho_antisym * math.sin(theta[i] - theta[j]) / (d ** 1.15)

            w[i, j] = sym + anti
            w[j, i] = sym - anti

    # Normalize by spectral norm to keep dynamics stable.
    smax = float(np.linalg.norm(w, 2))
    if smax > 1e-12:
        w = w / smax
    return w


def effective_path_sum(w: np.ndarray, lam: float, kmax: int) -> np.ndarray:
    n = w.shape[0]
    g = np.zeros((n, n), dtype=float)
    p = w.copy()
    coeff = 1.0
    for _ in range(1, kmax + 1):
        g += coeff * p
        coeff *= lam
        p = p @ w
    return g


def profile_from_matrix(g: np.ndarray, dmax: int = 12) -> np.ndarray:
    n = g.shape[0]
    dmax = min(dmax, n // 2)
    prof = []
    for d in range(1, dmax + 1):
        vals = []
        for i in range(n):
            j1 = (i + d) % n
            j2 = (i - d) % n
            vals.append(0.5 * (g[i, j1] + g[i, j2]))
        prof.append(float(np.mean(vals)))
    y = np.array(prof, dtype=float)
    if abs(y[0]) > 1e-12:
        y = y / abs(y[0])
    return y


def estimate_beta_from_envelope(y: np.ndarray) -> float:
    d = np.arange(1, len(y) + 1, dtype=float)
    a = np.abs(y)
    q = float(np.quantile(a, 0.35))
    mask = a > max(q, 1e-4)
    x = d[mask]
    z = 1.0 / np.clip(a[mask], 1e-6, None)

    beta = np.nan
    if len(x) >= 3:
        m, c = np.polyfit(x, z, deg=1)
        if c > 0 and m > 0:
            beta = float(m / c)

    if not np.isfinite(beta) or beta <= 0:
        beta_grid = np.linspace(0.001, 0.30, 1800)
        best_rmse = float("inf")
        best_beta = 0.05
        for b in beta_grid:
            env = 1.0 / (1.0 + b * d)
            # Best scale for abs(y) against env.
            scale = float(np.dot(a, env) / max(np.dot(env, env), 1e-12))
            pred = scale * env
            rmse = float(np.sqrt(np.mean((a - pred) ** 2)))
            if rmse < best_rmse:
                best_rmse = rmse
                best_beta = float(b)
        beta = best_beta

    return float(np.clip(beta, 0.001, 0.30))


def fit_phase_with_fixed_omega(y_flat: np.ndarray, d: np.ndarray, omega: float) -> Tuple[float, float, float]:
    c = np.cos(omega * d)
    s = np.sin(omega * d)
    b = np.column_stack([c, s])
    coeff, *_ = np.linalg.lstsq(b, y_flat, rcond=None)
    c0, s0 = float(coeff[0]), float(coeff[1])
    amp = float(np.hypot(c0, s0))
    phi = float(np.arctan2(-s0, c0))
    pred = amp * np.cos(omega * d + phi)
    sse = float(np.sum((y_flat - pred) ** 2))
    return amp, phi, sse


def derive_omega_phi_beta(y: np.ndarray) -> Dict[str, float]:
    d = np.arange(1, len(y) + 1, dtype=float)
    beta = estimate_beta_from_envelope(y)
    y_flat = y * (1.0 + beta * d)
    y_ctr = y_flat - np.mean(y_flat)

    # Initial omega from dominant spatial frequency.
    n_fft = 256
    spec = np.fft.rfft(y_ctr, n=n_fft)
    freqs = np.fft.rfftfreq(n_fft, d=1.0)  # cycles/sample
    power = np.abs(spec) ** 2
    if len(power) > 1:
        idx = int(np.argmax(power[1:]) + 1)
        omega0 = float(2.0 * math.pi * freqs[idx])
    else:
        omega0 = 0.5
    omega0 = float(np.clip(omega0, 0.05, 1.8))

    # Refine omega by SSE scan.
    omega_grid = np.linspace(max(0.05, omega0 - 0.25), min(1.8, omega0 + 0.25), 701)
    best = {"omega": omega0, "amp": 0.0, "phi": 0.0, "sse": float("inf")}
    for om in omega_grid:
        amp, phi, sse = fit_phase_with_fixed_omega(y_flat, d, float(om))
        if sse < best["sse"]:
            best = {"omega": float(om), "amp": amp, "phi": float(phi), "sse": sse}

    omega = best["omega"]
    amp = best["amp"]
    phi = best["phi"]
    y_hat = amp * np.cos(omega * d + phi) / (1.0 + beta * d)

    ss_res = float(np.sum((y - y_hat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 if ss_tot <= 1e-15 else float(1.0 - ss_res / ss_tot)
    rmse = float(np.sqrt(np.mean((y - y_hat) ** 2)))

    return {
        "omega_hat": omega,
        "phi_hat": phi,
        "beta_hat": beta,
        "amp_hat": amp,
        "rmse": rmse,
        "r2": r2,
    }


def circular_mean(x: np.ndarray) -> float:
    return float(np.arctan2(np.mean(np.sin(x)), np.mean(np.cos(x))))


def circular_std(x: np.ndarray) -> float:
    r = np.sqrt(np.mean(np.sin(x)) ** 2 + np.mean(np.cos(x)) ** 2)
    if r <= 1e-12:
        return float("inf")
    return float(np.sqrt(-2.0 * np.log(r)))


def main() -> None:
    n_grid = [64, 96]
    seeds_per_n = 14

    rows: List[Dict[str, object]] = []
    for n in n_grid:
        for sidx in range(seeds_per_n):
            seed = 173900 + 100 * n + sidx
            cfg = MicroCfg(n_nodes=n, seed=seed)
            rng = np.random.default_rng(seed)
            theta, q = build_micro_state(cfg, rng)
            w = build_signed_dynamic_matrix(cfg, theta, q, rng)
            g = effective_path_sum(w, lam=cfg.dyn_lambda, kmax=cfg.dyn_kmax)
            y = profile_from_matrix(g, dmax=12)
            est = derive_omega_phi_beta(y)

            rows.append(
                {
                    "n_nodes": n,
                    "seed": seed,
                    "omega_hat": est["omega_hat"],
                    "phi_hat": est["phi_hat"],
                    "beta_hat": est["beta_hat"],
                    "amp_hat": est["amp_hat"],
                    "rmse": est["rmse"],
                    "r2": est["r2"],
                    "profile": {str(i + 1): float(y[i]) for i in range(len(y))},
                }
            )

    omega = np.array([float(r["omega_hat"]) for r in rows], dtype=float)
    phi = np.array([float(r["phi_hat"]) for r in rows], dtype=float)
    beta = np.array([float(r["beta_hat"]) for r in rows], dtype=float)
    r2 = np.array([float(r["r2"]) for r in rows], dtype=float)
    rmse = np.array([float(r["rmse"]) for r in rows], dtype=float)

    omega_med = float(np.median(omega))
    omega_iqr = float(np.quantile(omega, 0.75) - np.quantile(omega, 0.25))
    phi_mean = circular_mean(phi)
    phi_std = circular_std(phi)
    beta_med = float(np.median(beta))
    beta_iqr = float(np.quantile(beta, 0.75) - np.quantile(beta, 0.25))
    r2_med = float(np.median(r2))
    rmse_med = float(np.median(rmse))

    # Reference values are for comparison only (not used in derivation).
    omega_ref = math.pi / 4.0
    phi_ref = math.pi / 6.0
    beta_ref = 0.01
    dev = {
        "omega_rel_vs_pi_over_4": abs(omega_med - omega_ref) / omega_ref,
        "phi_abs_over_pi_vs_pi_over_6": abs(((phi_mean - phi_ref + math.pi) % (2 * math.pi)) - math.pi) / math.pi,
        "beta_rel_vs_001": abs(beta_med - beta_ref) / beta_ref,
    }

    pass_fit = r2_med >= 0.65
    pass_stability = omega_iqr <= 0.22 and beta_iqr <= 0.08 and phi_std <= 1.0
    pass_reference = dev["omega_rel_vs_pi_over_4"] <= 0.20 and dev["phi_abs_over_pi_vs_pi_over_6"] <= 0.20

    if pass_fit and pass_stability and pass_reference:
        verdict = "SIGNED_MICROMODEL_DERIVATION_SUPPORTED_NEAR_CANONICAL"
    elif pass_fit and pass_stability:
        verdict = "SIGNED_MICROMODEL_DERIVATION_SUPPORTED_BUT_SHIFTED"
    elif pass_fit:
        verdict = "SIGNED_MICROMODEL_PARTIAL_LOW_IDENTIFIABILITY"
    else:
        verdict = "SIGNED_MICROMODEL_DERIVATION_NOT_CLOSED"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "n_grid": n_grid,
            "seeds_per_n": seeds_per_n,
            "n_total_runs": len(rows),
            "derivation_note": "omega/phi from spectrum+phase, beta from envelope; no direct parameter fixing.",
        },
        "summary": {
            "omega_median": omega_med,
            "omega_iqr": omega_iqr,
            "phi_circular_mean": phi_mean,
            "phi_circular_std": phi_std,
            "beta_median": beta_med,
            "beta_iqr": beta_iqr,
            "median_r2": r2_med,
            "median_rmse": rmse_med,
        },
        "reference_deviation": dev,
        "pass_flags": {
            "fit_quality": pass_fit,
            "stability": pass_stability,
            "near_reference": pass_reference,
        },
        "verdict": verdict,
        "rows": rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1739: SIGNED DYNAMIC MICROMODEL DERIVATION",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Runs: {len(rows)}",
        f"- omega median/IQR: {omega_med:.6f} / {omega_iqr:.6f}",
        f"- phi circular mean/std: {phi_mean:.6f} / {phi_std:.6f}",
        f"- beta median/IQR: {beta_med:.6f} / {beta_iqr:.6f}",
        f"- median R2={r2_med:.6f}, median RMSE={rmse_med:.6f}",
        f"- deviation vs refs: domega_rel={dev['omega_rel_vs_pi_over_4']:.3f}, dphi/pi={dev['phi_abs_over_pi_vs_pi_over_6']:.3f}, dbeta_rel={dev['beta_rel_vs_001']:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- fit_quality: {pass_fit}",
        f"- stability: {pass_stability}",
        f"- near_reference: {pass_reference}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1739] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1739] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
