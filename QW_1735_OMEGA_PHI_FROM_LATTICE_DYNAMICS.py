#!/usr/bin/env python3
"""
QW-1735: Derive (omega, phi) from lattice dynamics (not direct ansatz fixing).

Goal:
1) Simulate second-order wave dynamics on a topological lattice.
2) Extract distance-correlation profile C(d).
3) Fit C(d) with A*cos(omega*d+phi)/(1+beta*d) and estimate (omega, phi).
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1735_omega_phi_from_lattice_dynamics.json"
OUT_MD = ROOT / "RAPORT_QW1735_OMEGA_PHI_FROM_LATTICE_DYNAMICS.md"


def load_beta_prior(default: float = 0.01) -> float:
    p = ROOT / "report_qw1734_micro_beta_tors_derivation.json"
    if not p.exists():
        return default
    try:
        d = json.loads(p.read_text(encoding="utf-8"))
        x = float(d["summary"]["beta_hat_median"])
        if np.isfinite(x) and x > 0:
            return x
    except Exception:
        pass
    return default


def build_lattice(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    w = np.zeros((n, n), dtype=float)

    # Local and mid-range couplings.
    for i in range(n):
        for d, val in [(1, 1.0), (2, 0.45), (3, 0.18)]:
            j1 = (i + d) % n
            j2 = (i - d) % n
            w[i, j1] += val
            w[i, j2] += val

    # Sparse topological shortcuts.
    for i in range(n):
        for j in range(i + 4, n):
            dist = min(abs(i - j), n - abs(i - j))
            p = 0.05 / (dist ** 1.6)
            if rng.random() < p:
                theta = rng.uniform(-math.pi, math.pi)
                val = max(0.0, 0.38 * (1.0 + 0.2 * math.cos(theta)) * math.exp(-0.7 * dist / n))
                w[i, j] += val
                w[j, i] += val

    smax = float(np.linalg.norm(w, 2))
    if smax > 1e-12:
        w = w / smax
    return w


def simulate_wave_profile(
    w: np.ndarray,
    n_steps: int = 220,
    burn_in: int = 80,
    c2: float = 0.34,
    damp: float = 0.04,
    max_d: int = 12,
    seed: int = 0,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n = w.shape[0]
    x_prev = 0.02 * rng.normal(size=n)
    x_curr = np.zeros(n, dtype=float)
    x_curr[0] = 1.0
    x_curr += 0.01 * rng.normal(size=n)

    # Dynamic operator: second-order damped wave over graph Laplacian.
    lap = w - np.eye(n, dtype=float)
    snapshots: List[np.ndarray] = []

    for _ in range(n_steps):
        x_next = (
            2.0 * x_curr
            - x_prev
            + c2 * (lap @ x_curr)
            - damp * (x_curr - x_prev)
        )
        # Mild bounded nonlinearity for stability in long runs.
        x_next = np.tanh(1.2 * x_next)
        x_prev, x_curr = x_curr, x_next
        snapshots.append(x_curr.copy())

    arr = np.array(snapshots[burn_in:], dtype=float)  # [T, N]
    max_d = min(max_d, n // 2)
    prof = []
    for d in range(1, max_d + 1):
        vals = arr * np.roll(arr, -d, axis=1)
        prof.append(float(np.mean(vals)))
    y = np.array(prof, dtype=float)
    if abs(y[0]) < 1e-12:
        y[0] = 1e-12
    y = y / abs(y[0])
    return y


def fit_omega_phi(y: np.ndarray, beta: float) -> Dict[str, float]:
    d = np.arange(1, len(y) + 1, dtype=float)
    denom = 1.0 / (1.0 + beta * d)

    omega_grid_coarse = np.linspace(0.25, 1.35, 111)
    phi_grid_coarse = np.linspace(-math.pi, math.pi, 315)

    best = {
        "omega": float("nan"),
        "phi": float("nan"),
        "amp": float("nan"),
        "rmse": float("inf"),
        "r2": float("-inf"),
    }
    y_mean = float(np.mean(y))
    ss_tot = float(np.sum((y - y_mean) ** 2))

    def scan(omega_grid: np.ndarray, phi_grid: np.ndarray, best_now: Dict[str, float]) -> Dict[str, float]:
        best_local = dict(best_now)
        for omega in omega_grid:
            wd = omega * d
            for phi in phi_grid:
                basis = np.cos(wd + phi) * denom
                norm = float(np.dot(basis, basis))
                if norm < 1e-14:
                    continue
                amp = float(np.dot(y, basis) / norm)
                pred = amp * basis
                rmse = float(np.sqrt(np.mean((y - pred) ** 2)))
                if rmse < best_local["rmse"]:
                    if ss_tot <= 1e-15:
                        r2 = 1.0
                    else:
                        r2 = float(1.0 - np.sum((y - pred) ** 2) / ss_tot)
                    best_local = {
                        "omega": float(omega),
                        "phi": float(phi),
                        "amp": amp,
                        "rmse": rmse,
                        "r2": r2,
                    }
        return best_local

    best = scan(omega_grid_coarse, phi_grid_coarse, best)

    # Local refinement around coarse optimum.
    omega_center = best["omega"]
    phi_center = best["phi"]
    omega_grid_refine = np.linspace(max(0.2, omega_center - 0.06), min(1.4, omega_center + 0.06), 121)
    phi_grid_refine = np.linspace(phi_center - 0.35, phi_center + 0.35, 121)
    best = scan(omega_grid_refine, phi_grid_refine, best)

    return best


def local_minima_abs(y: np.ndarray) -> List[int]:
    a = np.abs(y)
    mins: List[int] = []
    for i in range(1, len(a) - 1):
        if a[i] <= a[i - 1] and a[i] <= a[i + 1]:
            mins.append(i + 1)  # convert to 1-based d
    return mins


def circular_mean(phases: np.ndarray) -> float:
    s = np.mean(np.sin(phases))
    c = np.mean(np.cos(phases))
    return float(np.arctan2(s, c))


def circular_std(phases: np.ndarray) -> float:
    r = np.sqrt(np.mean(np.sin(phases)) ** 2 + np.mean(np.cos(phases)) ** 2)
    if r <= 1e-12:
        return float("inf")
    return float(np.sqrt(-2.0 * np.log(r)))


def main() -> None:
    beta_prior = load_beta_prior(default=0.01)

    n = 84
    seeds = [173500 + i for i in range(16)]
    rows: List[Dict[str, object]] = []

    for seed in seeds:
        w = build_lattice(n, seed)
        y = simulate_wave_profile(w, n_steps=220, burn_in=80, c2=0.34, damp=0.04, max_d=12, seed=seed + 19)
        fit = fit_omega_phi(y, beta=beta_prior)
        mins = local_minima_abs(y)
        rows.append(
            {
                "seed": seed,
                "beta_used": beta_prior,
                "omega_hat": fit["omega"],
                "phi_hat": fit["phi"],
                "amp_hat": fit["amp"],
                "rmse": fit["rmse"],
                "r2": fit["r2"],
                "minima_abs_positions": mins,
                "profile": {str(i + 1): float(y[i]) for i in range(len(y))},
            }
        )

    omega = np.array([float(r["omega_hat"]) for r in rows], dtype=float)
    phi = np.array([float(r["phi_hat"]) for r in rows], dtype=float)
    rmse = np.array([float(r["rmse"]) for r in rows], dtype=float)
    r2 = np.array([float(r["r2"]) for r in rows], dtype=float)

    omega_med = float(np.median(omega))
    omega_std = float(np.std(omega))
    phi_cmean = float(circular_mean(phi))
    phi_cstd = float(circular_std(phi))
    rmse_med = float(np.median(rmse))
    r2_med = float(np.median(r2))

    omega_ref = math.pi / 4.0
    phi_ref = math.pi / 6.0
    omega_rel_dev = abs(omega_med - omega_ref) / omega_ref
    phi_dev = abs(((phi_cmean - phi_ref + math.pi) % (2 * math.pi)) - math.pi)
    phi_rel_pi = phi_dev / math.pi

    stable_pair = (omega_std <= 0.08) and (phi_cstd <= 0.35) and (r2_med >= 0.85)
    near_reference = (omega_rel_dev <= 0.10) and (phi_rel_pi <= 0.10)

    if stable_pair and near_reference:
        verdict = "OMEGA_PHI_DERIVATION_SUPPORTED_NEAR_CANONICAL"
    elif stable_pair:
        verdict = "OMEGA_PHI_DERIVATION_SUPPORTED_BUT_SHIFTED"
    else:
        verdict = "OMEGA_PHI_DERIVATION_NOT_STABLE"

    minima_pool: List[int] = []
    for r in rows:
        minima_pool.extend([int(x) for x in r["minima_abs_positions"]])
    minima_hist: Dict[str, int] = {}
    for m in sorted(minima_pool):
        k = str(m)
        minima_hist[k] = minima_hist.get(k, 0) + 1

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "n_nodes": n,
            "n_runs": len(rows),
            "seeds": seeds,
            "beta_used_for_fit": beta_prior,
            "wave_params": {"n_steps": 220, "burn_in": 80, "c2": 0.34, "damp": 0.04},
        },
        "summary": {
            "omega_median": omega_med,
            "omega_std": omega_std,
            "phi_circular_mean": phi_cmean,
            "phi_circular_std": phi_cstd,
            "median_rmse": rmse_med,
            "median_r2": r2_med,
            "omega_rel_dev_vs_pi_over_4": omega_rel_dev,
            "phi_abs_dev_over_pi_vs_pi_over_6": phi_rel_pi,
        },
        "reference": {"omega_pi_over_4": omega_ref, "phi_pi_over_6": phi_ref},
        "minima_histogram_abs_profile": minima_hist,
        "pass_flags": {
            "stable_pair": stable_pair,
            "near_reference": near_reference,
        },
        "verdict": verdict,
        "rows": rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1735: OMEGA/PHI FROM LATTICE DYNAMICS",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Runs: {len(rows)}",
        f"- beta used: {beta_prior:.6f}",
        f"- omega median/std: {omega_med:.6f} / {omega_std:.6f}",
        f"- phi circular mean/std: {phi_cmean:.6f} / {phi_cstd:.6f}",
        f"- median RMSE={rmse_med:.6f}, median R2={r2_med:.6f}",
        f"- dev vs refs: d_omega_rel={omega_rel_dev:.3f}, d_phi/pi={phi_rel_pi:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- stable_pair: {stable_pair}",
        f"- near_reference: {near_reference}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1735] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1735] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
