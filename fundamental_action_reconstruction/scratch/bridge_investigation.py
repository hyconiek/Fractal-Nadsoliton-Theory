#!/usr/bin/env python3
"""Scratch investigation: legacy->strict bridge hypothesis stress-test.

This is intentionally scratch-scoped and does not alter canonical artifacts.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import scipy.optimize as so
import scipy.stats as ss
import sympy as sp

HERE = Path(__file__).resolve().parent
OUT = HERE / "bridge_analysis_report.json"


def k_legacy(d: np.ndarray, alpha_geo: float, omega: float, phi: float, beta_tors: float) -> np.ndarray:
    return alpha_geo * np.cos(omega * d + phi) / (1.0 + beta_tors * d)


def k_strict(d: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * (d ** eta))


def main() -> None:
    d = np.linspace(1.0, 11.0, 400)
    omega, phi, beta, eta = 0.18575, 0.16250, 1.0, 1.8
    alpha_geo = float(4.0 * np.log(2.0))
    y_strict = k_strict(d, omega, phi, beta, eta)

    # H1: simple fractal-distance + dephasing + torsion fit
    def obj_h1(x: np.ndarray) -> float:
        lam, a, dphi, bt = x
        y_leg = k_legacy(lam * (d**a), alpha_geo, omega, phi + dphi, bt)
        c = float(np.dot(y_leg, y_strict) / max(np.dot(y_leg, y_leg), 1e-15))
        return float(np.mean((c * y_leg - y_strict) ** 2))

    h1 = so.minimize(obj_h1, x0=np.array([1.0, 1.0, 0.0, 1.0]), bounds=[(0.1, 5.0), (0.01, 2.5), (-np.pi, np.pi), (1e-4, 10.0)], method="L-BFGS-B")
    lam, a, dphi, bt = [float(v) for v in h1.x]
    y1 = k_legacy(lam * (d**a), alpha_geo, omega, phi + dphi, bt)
    c1 = float(np.dot(y1, y_strict) / max(np.dot(y1, y1), 1e-15))
    y1f = c1 * y1
    ss_res = float(np.sum((y1f - y_strict) ** 2))
    ss_tot = float(np.sum((y_strict - np.mean(y_strict)) ** 2))
    r2_h1 = float(1.0 - ss_res / max(ss_tot, 1e-15))
    max_gap_h1 = float(np.max(np.abs(y1f - y_strict)))

    # H2: force alpha=1 (no fractal exponent), fit only lambda/dephase/beta
    def obj_h2(x: np.ndarray) -> float:
        lam2, dphi2, bt2 = x
        y_leg = k_legacy(lam2 * d, alpha_geo, omega, phi + dphi2, bt2)
        c = float(np.dot(y_leg, y_strict) / max(np.dot(y_leg, y_leg), 1e-15))
        return float(np.mean((c * y_leg - y_strict) ** 2))

    h2 = so.minimize(obj_h2, x0=np.array([1.0, 0.0, 1.0]), bounds=[(0.1, 5.0), (-np.pi, np.pi), (1e-4, 10.0)], method="L-BFGS-B")
    lam2, dphi2, bt2 = [float(v) for v in h2.x]
    y2 = k_legacy(lam2 * d, alpha_geo, omega, phi + dphi2, bt2)
    c2 = float(np.dot(y2, y_strict) / max(np.dot(y2, y2), 1e-15))
    y2f = c2 * y2
    ss_res2 = float(np.sum((y2f - y_strict) ** 2))
    r2_h2 = float(1.0 - ss_res2 / max(ss_tot, 1e-15))

    # symbolic non-identity check
    ds, b1, b2 = sp.symbols("d b1 b2", positive=True, real=True)
    expr = sp.simplify((1 + b1 * ds) - (1 + b2 * ds ** sp.Rational(9, 5)))

    # simple bootstrap around residual for H1 fit quality uncertainty
    rng = np.random.default_rng(20260519)
    resid = y_strict - y1f
    r2_boot = []
    for _ in range(300):
        yb = y1f + rng.choice(resid, size=resid.size, replace=True)
        ssr = float(np.sum((yb - y1f) ** 2))
        sst = float(np.sum((yb - np.mean(yb)) ** 2))
        r2_boot.append(float(1.0 - ssr / max(sst, 1e-15)))
    q05, q50, q95 = [float(x) for x in np.quantile(np.array(r2_boot), [0.05, 0.5, 0.95])]

    # one-sided test whether H1 materially exceeds H2
    delta = np.array(r2_boot) - r2_h2
    p_nonworse = float(np.mean(delta > 0))
    pval = float(ss.wilcoxon(delta, alternative="greater", zero_method="zsplit", method="approx").pvalue) if np.any(np.abs(delta) > 1e-15) else 1.0

    report = {
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "result_kind": "SCRATCH_BRIDGE_HYPOTHESIS_TEST__NOT_A_THEOREM",
        "hypothesis_h1_fractal_map": {"lambda": lam, "alpha": a, "delta_phi": dphi, "beta_tors_eff": bt, "amplitude_renorm_c": c1, "r2": r2_h1, "max_abs_gap": max_gap_h1},
        "hypothesis_h2_no_fractal_exponent": {"lambda": lam2, "delta_phi": dphi2, "beta_tors_eff": bt2, "amplitude_renorm_c": c2, "r2": r2_h2},
        "model_comparison": {"delta_r2_h1_minus_h2": float(r2_h1 - r2_h2), "bootstrap_r2_h1_q05_q50_q95": [q05, q50, q95], "probability_h1_better_than_h2": p_nonworse, "wilcoxon_pvalue_h1_gt_h2": pval},
        "symbolic_nonidentity": {"expr_linear_minus_fractional": sp.sstr(expr), "exact_zero_identity": bool(sp.simplify(expr) == 0)},
        "hard_limits": [
            "High R2 cannot prove a bridge theorem.",
            "No legacy-role transfer to strict is inferred.",
            "Strict selector closure remains open."
        ],
    }
    OUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()

