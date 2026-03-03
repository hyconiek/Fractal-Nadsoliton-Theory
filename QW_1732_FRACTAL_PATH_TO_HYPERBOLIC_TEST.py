#!/usr/bin/env python3
"""
QW-1732: Test robustness of exp/path-sum -> hyperbolic denominator transition.

Purpose:
1) Probe if simplified fractal path model naturally yields 1/(1+beta*d) shape.
2) Quantify whether fit is robust or requires narrow tuning.
3) Provide a stricter empirical basis for the denominator derivation claim.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1732_fractal_path_to_hyperbolic_test.json"
OUT_MD = ROOT / "RAPORT_QW1732_FRACTAL_PATH_TO_HYPERBOLIC_TEST.md"


def fractal_path_shape(d: np.ndarray, nu: float, lam: float) -> np.ndarray:
    # N_paths(d) ~ d^nu and per-path decay ~ (1+d)^(-lam)
    y = np.power(d, nu) / np.power(1.0 + d, lam)
    y = y / y[0]  # normalize at d=1
    return y


def hyperbolic_shape(d: np.ndarray, beta: float) -> np.ndarray:
    y = 1.0 / (1.0 + beta * (d - 1.0))
    y = y / y[0]  # also normalized at d=1
    return y


def r2_score(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    ss_res = float(np.sum((y_true - y_pred) ** 2))
    ss_tot = float(np.sum((y_true - np.mean(y_true)) ** 2))
    if ss_tot <= 1e-15:
        return 1.0
    return 1.0 - ss_res / ss_tot


def best_beta_fit_from_precomputed(
    target: np.ndarray, beta_grid: np.ndarray, hyper_grid_matrix: np.ndarray
) -> Dict[str, float]:
    # hyper_grid_matrix shape: [n_beta, n_d]
    residual = hyper_grid_matrix - target[None, :]
    mse = np.mean(residual * residual, axis=1)
    idx = int(np.argmin(mse))
    pred = hyper_grid_matrix[idx]
    rmse = float(np.sqrt(mse[idx]))
    return {"beta": float(beta_grid[idx]), "rmse": rmse, "r2": r2_score(target, pred)}


def main() -> None:
    d = np.arange(1.0, 13.0, dtype=float)
    nu_grid = np.linspace(0.0, 3.0, 91)
    lam_grid = np.linspace(0.0, 4.0, 121)
    beta_grid = np.linspace(0.001, 0.2, 800)

    # Precompute all candidate hyperbolic shapes once for speed.
    hyper_grid_matrix = np.array([hyperbolic_shape(d, float(b)) for b in beta_grid], dtype=float)

    rows: List[Dict[str, float]] = []
    for nu in nu_grid:
        for lam in lam_grid:
            target = fractal_path_shape(d, float(nu), float(lam))
            fit = best_beta_fit_from_precomputed(target, beta_grid, hyper_grid_matrix)
            rows.append(
                {
                    "nu": float(nu),
                    "lam": float(lam),
                    "best_beta": fit["beta"],
                    "rmse": fit["rmse"],
                    "r2": fit["r2"],
                    "delta_lam_minus_nu": float(lam - nu),
                }
            )

    rows_sorted = sorted(rows, key=lambda r: r["rmse"])
    best = rows_sorted[0]

    r2_hi = [r for r in rows if r["r2"] >= 0.995]
    r2_mid = [r for r in rows if r["r2"] >= 0.98]
    total = len(rows)

    frac_hi = len(r2_hi) / total
    frac_mid = len(r2_mid) / total

    if r2_hi:
        delta_vals = np.array([r["delta_lam_minus_nu"] for r in r2_hi], dtype=float)
        beta_vals = np.array([r["best_beta"] for r in r2_hi], dtype=float)
        ridge_stats = {
            "count": len(r2_hi),
            "delta_lam_minus_nu_mean": float(np.mean(delta_vals)),
            "delta_lam_minus_nu_std": float(np.std(delta_vals)),
            "best_beta_mean": float(np.mean(beta_vals)),
            "best_beta_std": float(np.std(beta_vals)),
        }
    else:
        ridge_stats = {
            "count": 0,
            "delta_lam_minus_nu_mean": None,
            "delta_lam_minus_nu_std": None,
            "best_beta_mean": None,
            "best_beta_std": None,
        }

    tuning_required = frac_hi < 0.1
    if best["r2"] < 0.95:
        verdict = "HYPERBOLIC_REDUCTION_NOT_SUPPORTED"
    elif tuning_required:
        verdict = "HYPERBOLIC_REDUCTION_PLAUSIBLE_BUT_TUNED"
    else:
        verdict = "HYPERBOLIC_REDUCTION_ROBUST"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "grid": {
            "nu_range": [float(nu_grid[0]), float(nu_grid[-1]), len(nu_grid)],
            "lam_range": [float(lam_grid[0]), float(lam_grid[-1]), len(lam_grid)],
            "beta_range": [float(beta_grid[0]), float(beta_grid[-1]), len(beta_grid)],
            "total_models": total,
        },
        "best_global_fit": best,
        "fit_coverage": {
            "fraction_r2_ge_0_995": frac_hi,
            "fraction_r2_ge_0_98": frac_mid,
        },
        "high_quality_ridge_stats": ridge_stats,
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1732: FRACTAL PATH -> HYPERBOLIC TEST",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Total models scanned: {total}",
        f"- Best fit: nu={best['nu']:.3f}, lam={best['lam']:.3f}, beta={best['best_beta']:.5f}",
        f"- Best RMSE={best['rmse']:.6e}, R2={best['r2']:.6f}",
        f"- Coverage R2>=0.995: {frac_hi:.4f}",
        f"- Coverage R2>=0.98: {frac_mid:.4f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Ridge Statistics (R2>=0.995)",
        f"- count: {ridge_stats['count']}",
        f"- mean(lam-nu): {ridge_stats['delta_lam_minus_nu_mean']}",
        f"- std(lam-nu): {ridge_stats['delta_lam_minus_nu_std']}",
        f"- mean(beta): {ridge_stats['best_beta_mean']}",
        f"- std(beta): {ridge_stats['best_beta_std']}",
        "",
        "## Interpretation",
        "- If high-quality region is narrow, path->hyperbolic claim is mechanistically plausible but not uniquely derived.",
        "- If high-quality region is broad, reduction is robust under model perturbations.",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1732] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1732] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
