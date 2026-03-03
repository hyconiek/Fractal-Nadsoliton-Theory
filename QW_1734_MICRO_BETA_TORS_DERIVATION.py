#!/usr/bin/env python3
"""
QW-1734: Microscopic derivation test for beta_tors.

Goal:
1) Generate effective inter-octave coupling from a micro topological model.
2) Fit hyperbolic beta in K_eff(d) ~ A/(1+beta*d) without using target beta=0.01.
3) Evaluate out-of-sample stability across network size/seed.
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
OUT_JSON = ROOT / "report_qw1734_micro_beta_tors_derivation.json"
OUT_MD = ROOT / "RAPORT_QW1734_MICRO_BETA_TORS_DERIVATION.md"


@dataclass
class Config:
    n_nodes: int
    p_shortcut: float
    eta_shortcut: float
    tau_torsion: float
    gamma_prop: float
    seed: int


def ring_distance(i: int, j: int, n: int) -> int:
    d = abs(i - j)
    return min(d, n - d)


def build_weight_matrix(cfg: Config, rng: np.random.Generator) -> np.ndarray:
    n = cfg.n_nodes
    w = np.zeros((n, n), dtype=float)

    # Local lattice (nearest + next-nearest) is fixed and non-fitted.
    for i in range(n):
        j1 = (i + 1) % n
        j2 = (i - 1) % n
        w[i, j1] += 1.0
        w[i, j2] += 1.0
        k1 = (i + 2) % n
        k2 = (i - 2) % n
        w[i, k1] += 0.35
        w[i, k2] += 0.35

    # Fractal shortcuts with torsion-attenuated transfer.
    for i in range(n):
        for j in range(i + 3, n):
            d = ring_distance(i, j, n)
            if d < 3:
                continue
            p = cfg.p_shortcut / (d ** cfg.eta_shortcut)
            if rng.random() < p:
                theta = rng.uniform(-math.pi, math.pi)
                torsion_factor = 1.0 + 0.25 * math.cos(theta)
                dist_factor = math.exp(-cfg.tau_torsion * (d / n) * 4.0)
                val = max(0.0, 0.45 * torsion_factor * dist_factor)
                if val > 0:
                    w[i, j] += val
                    w[j, i] += val

    # Stable normalization by largest singular value.
    smax = float(np.linalg.norm(w, 2))
    if smax > 1e-12:
        w = w / smax
    return w


def propagator_series(w: np.ndarray, gamma: float, kmax: int = 8) -> np.ndarray:
    n = w.shape[0]
    g = np.zeros((n, n), dtype=float)
    term = np.eye(n, dtype=float)
    gw = gamma * w
    for _ in range(kmax):
        term = term @ gw
        g += term
    return g


def distance_profile(g: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    n = g.shape[0]
    max_d = min(12, n // 2)
    dvals = np.arange(1, max_d + 1, dtype=float)
    y = []
    for d in range(1, max_d + 1):
        vals = []
        for i in range(n):
            j = (i + d) % n
            vals.append(abs(g[i, j]))
            j2 = (i - d) % n
            vals.append(abs(g[i, j2]))
        y.append(float(np.mean(vals)))
    y = np.array(y, dtype=float)
    if y[0] <= 1e-12:
        y = y + 1e-12
    y = y / y[0]
    return dvals, y


def fit_beta_hyperbolic(dvals: np.ndarray, y: np.ndarray) -> Dict[str, float]:
    beta_grid = np.linspace(0.001, 0.2, 2000)
    pred = 1.0 / (1.0 + beta_grid[:, None] * (dvals[None, :] - 1.0))
    mse = np.mean((pred - y[None, :]) ** 2, axis=1)
    idx = int(np.argmin(mse))
    beta = float(beta_grid[idx])
    yhat = pred[idx]
    rmse = float(np.sqrt(mse[idx]))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    if ss_tot <= 1e-15:
        r2 = 1.0
    else:
        r2 = float(1.0 - np.sum((y - yhat) ** 2) / ss_tot)
    return {"beta_hat": beta, "rmse": rmse, "r2": r2}


def main() -> None:
    configs: List[Config] = []
    n_nodes_grid = [48, 72]
    p_short_grid = [0.03, 0.06, 0.1]
    eta_grid = [1.4, 1.8]
    tau_grid = [0.25, 0.5, 0.9]
    gamma_grid = [0.18, 0.24]
    seeds = [173401 + i for i in range(5)]

    for n in n_nodes_grid:
        for p in p_short_grid:
            for eta in eta_grid:
                for tau in tau_grid:
                    for gamma in gamma_grid:
                        for seed in seeds:
                            configs.append(
                                Config(
                                    n_nodes=n,
                                    p_shortcut=p,
                                    eta_shortcut=eta,
                                    tau_torsion=tau,
                                    gamma_prop=gamma,
                                    seed=seed,
                                )
                            )

    rows: List[Dict[str, float]] = []
    for cfg in configs:
        rng = np.random.default_rng(cfg.seed)
        w = build_weight_matrix(cfg, rng)
        g = propagator_series(w, cfg.gamma_prop, kmax=8)
        dvals, y = distance_profile(g)
        fit = fit_beta_hyperbolic(dvals, y)
        rows.append(
            {
                "n_nodes": cfg.n_nodes,
                "p_shortcut": cfg.p_shortcut,
                "eta_shortcut": cfg.eta_shortcut,
                "tau_torsion": cfg.tau_torsion,
                "gamma_prop": cfg.gamma_prop,
                "seed": cfg.seed,
                "beta_hat": fit["beta_hat"],
                "rmse": fit["rmse"],
                "r2": fit["r2"],
            }
        )

    beta_all = np.array([r["beta_hat"] for r in rows], dtype=float)
    r2_all = np.array([r["r2"] for r in rows], dtype=float)
    rmse_all = np.array([r["rmse"] for r in rows], dtype=float)

    by_n: Dict[int, np.ndarray] = {}
    for n in n_nodes_grid:
        by_n[n] = np.array([r["beta_hat"] for r in rows if r["n_nodes"] == n], dtype=float)

    median_beta = float(np.median(beta_all))
    q10, q90 = float(np.quantile(beta_all, 0.1)), float(np.quantile(beta_all, 0.9))
    iqr = float(np.quantile(beta_all, 0.75) - np.quantile(beta_all, 0.25))
    median_r2 = float(np.median(r2_all))
    median_rmse = float(np.median(rmse_all))

    target_beta = 0.01
    rel_err = abs(median_beta - target_beta) / target_beta
    oos_shift = float(abs(np.median(by_n[n_nodes_grid[0]]) - np.median(by_n[n_nodes_grid[1]])))

    pass_beta = rel_err <= 0.10
    pass_stability = (iqr <= 0.015) and (oos_shift <= 0.01)
    pass_fit = median_r2 >= 0.95

    if pass_beta and pass_stability and pass_fit:
        verdict = "BETA_TORS_MICRO_DERIVATION_SUPPORTED"
    elif pass_fit and pass_stability:
        verdict = "BETA_TORS_MICRO_DERIVATION_PARTIAL_NOT_TARGETED"
    else:
        verdict = "BETA_TORS_MICRO_DERIVATION_NOT_CLOSED"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "n_nodes_grid": n_nodes_grid,
            "p_shortcut_grid": p_short_grid,
            "eta_shortcut_grid": eta_grid,
            "tau_torsion_grid": tau_grid,
            "gamma_grid": gamma_grid,
            "n_seeds_per_config": len(seeds),
            "n_total_runs": len(rows),
        },
        "summary": {
            "beta_hat_median": median_beta,
            "beta_hat_q10": q10,
            "beta_hat_q90": q90,
            "beta_hat_iqr": iqr,
            "median_r2": median_r2,
            "median_rmse": median_rmse,
            "relative_error_to_target_beta_001": rel_err,
            "oos_median_shift_between_sizes": oos_shift,
        },
        "by_size": {
            str(n): {
                "median": float(np.median(arr)),
                "q10": float(np.quantile(arr, 0.1)),
                "q90": float(np.quantile(arr, 0.9)),
            }
            for n, arr in by_n.items()
        },
        "pass_flags": {
            "target_beta_match_pm10pct": pass_beta,
            "oos_stability": pass_stability,
            "fit_quality": pass_fit,
        },
        "verdict": verdict,
        "rows_head_120": rows[:120],
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1734: MICRO BETA_TORS DERIVATION",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Runs: {len(rows)}",
        f"- beta_hat median: {median_beta:.6f}",
        f"- beta_hat q10..q90: [{q10:.6f}, {q90:.6f}]",
        f"- median R2: {median_r2:.4f}",
        f"- rel.error vs 0.01: {rel_err:.3f}",
        f"- OOS size-shift: {oos_shift:.6f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- target_beta_match_pm10pct: {pass_beta}",
        f"- oos_stability: {pass_stability}",
        f"- fit_quality: {pass_fit}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1734] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1734] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
