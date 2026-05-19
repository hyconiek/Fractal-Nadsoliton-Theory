#!/usr/bin/env python3
"""Scratch verification: identifiability/stability of bridge-fit parameters.

Strictly exploratory, non-theorem, guardrail-compliant.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import scipy.optimize as so
import scipy.stats as ss

HERE = Path(__file__).resolve().parent
OUT = HERE / "bridge_verification_report.json"


def k_legacy(d: np.ndarray, alpha_geo: float, omega: float, phi: float, beta_tors: float) -> np.ndarray:
    return alpha_geo * np.cos(omega * d + phi) / (1.0 + beta_tors * d)


def k_strict(d: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * (d ** eta))


def fit_once(d: np.ndarray, y: np.ndarray, omega: float, phi: float, alpha_geo: float, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    x0 = np.array([
        float(rng.uniform(0.1, 4.0)),
        float(rng.uniform(0.05, 2.2)),
        float(rng.uniform(-np.pi, np.pi)),
        float(rng.uniform(0.01, 8.0)),
    ], dtype=float)

    def obj(x: np.ndarray) -> float:
        lam, a, dphi, bt = x
        y_leg = k_legacy(lam * (d**a), alpha_geo, omega, phi + dphi, bt)
        c = float(np.dot(y_leg, y) / max(np.dot(y_leg, y_leg), 1e-15))
        return float(np.mean((c * y_leg - y) ** 2))

    res = so.minimize(
        obj,
        x0=x0,
        bounds=[(0.1, 5.0), (0.01, 2.5), (-np.pi, np.pi), (1e-4, 10.0)],
        method="L-BFGS-B",
    )
    lam, a, dphi, bt = [float(v) for v in res.x]
    y_leg = k_legacy(lam * (d**a), alpha_geo, omega, phi + dphi, bt)
    c = float(np.dot(y_leg, y) / max(np.dot(y_leg, y_leg), 1e-15))
    y_fit = c * y_leg
    ss_res = float(np.sum((y_fit - y) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = float(1.0 - ss_res / max(ss_tot, 1e-15))
    return {"lambda": lam, "alpha": a, "delta_phi": dphi, "beta_tors_eff": bt, "c": c, "r2": r2}


def main() -> None:
    d = np.linspace(1.0, 11.0, 400)
    omega, phi, beta, eta = 0.18575, 0.16250, 1.0, 1.8
    alpha_geo = float(4.0 * np.log(2.0))
    y = k_strict(d, omega, phi, beta, eta)

    rows = [fit_once(d, y, omega, phi, alpha_geo, 20260519 + i) for i in range(64)]
    a_vals = np.array([r["alpha"] for r in rows], dtype=float)
    lam_vals = np.array([r["lambda"] for r in rows], dtype=float)
    r2_vals = np.array([r["r2"] for r in rows], dtype=float)

    alpha_span = float(np.max(a_vals) - np.min(a_vals))
    lambda_span = float(np.max(lam_vals) - np.min(lam_vals))
    corr_alpha_lambda = float(ss.pearsonr(a_vals, lam_vals).statistic)
    r2_q = [float(x) for x in np.quantile(r2_vals, [0.05, 0.5, 0.95])]

    identifiability_flag = bool(alpha_span < 0.25 and lambda_span < 1.0)

    report = {
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "result_kind": "SCRATCH_BRIDGE_IDENTIFIABILITY_CHECK__NOT_A_THEOREM",
        "num_multistarts": 64,
        "fit_rows": rows[:12],  # preview only
        "summary": {
            "alpha_span": alpha_span,
            "lambda_span": lambda_span,
            "corr_alpha_lambda": corr_alpha_lambda,
            "r2_q05_q50_q95": r2_q,
            "identifiability_flag_strict_threshold": identifiability_flag,
        },
        "hard_limits": [
            "Parameter stability is necessary but not sufficient for a bridge theorem.",
            "No transfer of legacy physical-role claims is made.",
            "Kernel split remains open without explicit theorem.",
        ],
        "next_honest_step": "Promote this to theorem-prep only if identifiability remains stable under changed d-ranges and noise models; otherwise classify as non-bridge fit degeneracy.",
    }
    OUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()

