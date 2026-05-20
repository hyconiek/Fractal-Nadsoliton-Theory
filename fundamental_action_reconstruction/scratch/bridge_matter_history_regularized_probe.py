#!/usr/bin/env python3
"""Scratch probe: regularized/tighter-prior memory hypothesis check.

Purpose:
- test whether mild physically-motivated priors reduce parameter degeneracy,
- keep conclusion non-theorem and guardrail-compliant.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import scipy.optimize as so
import sympy as sp

HERE = Path(__file__).resolve().parent
OUT = HERE / "bridge_matter_history_regularized_report.json"


@dataclass
class Scenario:
    name: str
    d_min: float
    d_max: float
    npts: int
    phase_shift: float
    eta_shift: float


def k_legacy(d: np.ndarray, alpha_geo: float, omega: float, phi: float, beta_tors: float) -> np.ndarray:
    return alpha_geo * np.cos(omega * d + phi) / (1.0 + beta_tors * d)


def k_strict(d: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * (d**eta))


def exp_memory_filter(x: np.ndarray, rho: float) -> np.ndarray:
    y = np.empty_like(x)
    y[0] = x[0]
    a = 1.0 - rho
    for i in range(1, x.size):
        y[i] = a * x[i] + rho * y[i - 1]
    return y


def fit_regularized(d: np.ndarray, y: np.ndarray, alpha_geo: float, omega: float, phi: float, seed: int, starts: int = 10) -> dict:
    rng = np.random.default_rng(seed)

    # Tighter prior window chosen as methodological regularization, not physical theorem.
    rho_center, rho_scale = 0.22, 0.08
    gamma_center, gamma_scale = 1.0, 1.2
    lam = 3e-3

    def obj(x: np.ndarray) -> float:
        dphi, bt, rho, gamma = x
        raw = k_legacy(d, alpha_geo, omega, phi + dphi, bt)
        mem = exp_memory_filter(raw, rho)
        feat = np.vstack([raw, gamma * mem]).T
        coef, *_ = np.linalg.lstsq(feat, y, rcond=None)
        yhat = feat @ coef
        mse = float(np.mean((yhat - y) ** 2))
        prior = ((rho - rho_center) / rho_scale) ** 2 + ((gamma - gamma_center) / gamma_scale) ** 2
        return mse + lam * prior

    best = None
    for _ in range(starts):
        x0 = np.array([
            rng.uniform(-0.4, 0.4),
            rng.uniform(1e-3, 30.0),
            rng.uniform(0.05, 0.45),
            rng.uniform(-1.5, 3.5),
        ])
        res = so.minimize(
            obj,
            x0=x0,
            bounds=[(-np.pi, np.pi), (1e-5, 40.0), (0.01, 0.60), (-2.0, 4.0)],
            method="L-BFGS-B",
        )
        if best is None or res.fun < best.fun:
            best = res

    dphi, bt, rho, gamma = [float(v) for v in best.x]
    raw = k_legacy(d, alpha_geo, omega, phi + dphi, bt)
    mem = exp_memory_filter(raw, rho)
    feat = np.vstack([raw, gamma * mem]).T
    coef, *_ = np.linalg.lstsq(feat, y, rcond=None)
    yhat = feat @ coef
    resid = yhat - y
    sse = float(np.sum(resid * resid))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = float(1.0 - sse / max(ss_tot, 1e-15))

    return {
        "delta_phi": dphi,
        "beta_tors": bt,
        "rho_memory": rho,
        "gamma_memory": gamma,
        "coef_raw": float(coef[0]),
        "coef_memory": float(coef[1]),
        "objective": float(best.fun),
        "r2": r2,
    }


def main() -> None:
    alpha_geo = float(4.0 * np.log(2.0))
    omega0, phi0, beta0, eta0 = 0.18575, 0.16250, 1.0, 1.8

    scenarios = [
        Scenario("base_window", 1.0, 11.0, 700, 0.00, 0.00),
        Scenario("short_window", 1.0, 8.0, 640, 0.00, 0.00),
        Scenario("long_window", 1.0, 14.0, 800, 0.00, 0.00),
        Scenario("phase_eta_perturbed_plus", 1.0, 11.0, 700, +0.04, +0.10),
        Scenario("phase_eta_perturbed_minus", 1.0, 11.0, 700, -0.04, -0.10),
    ]

    seeds = [20260519 + i for i in range(8)]
    rows = []
    for sc in scenarios:
        d = np.linspace(sc.d_min, sc.d_max, sc.npts)
        y = k_strict(d, omega0, phi0 + sc.phase_shift, beta0, eta0 + sc.eta_shift)
        fits = [fit_regularized(d, y, alpha_geo, omega0, phi0 + sc.phase_shift, sd, starts=8) for sd in seeds]
        rho = np.array([f["rho_memory"] for f in fits], dtype=float)
        gamma = np.array([f["gamma_memory"] for f in fits], dtype=float)
        r2 = np.array([f["r2"] for f in fits], dtype=float)
        rows.append({
            "scenario": asdict(sc),
            "seed_count": len(seeds),
            "rho_q05_q50_q95": [float(v) for v in np.quantile(rho, [0.05, 0.5, 0.95])],
            "gamma_q05_q50_q95": [float(v) for v in np.quantile(gamma, [0.05, 0.5, 0.95])],
            "r2_q05_q50_q95": [float(v) for v in np.quantile(r2, [0.05, 0.5, 0.95])],
            "rho_cv": float(np.std(rho) / max(abs(np.mean(rho)), 1e-12)),
            "gamma_cv": float(np.std(gamma) / max(abs(np.mean(gamma)), 1e-12)),
        })

    rho_cv_all = np.array([r["rho_cv"] for r in rows], dtype=float)
    gamma_cv_all = np.array([r["gamma_cv"] for r in rows], dtype=float)
    stable_flag = bool(np.max(rho_cv_all) < 0.35 and np.max(gamma_cv_all) < 0.60)

    # symbolic prior expression
    rho, gamma, rho0, g0, sr, sg, lam = sp.symbols("rho gamma rho0 g0 sr sg lam", positive=True, real=True)
    reg_expr = sp.expand(lam * (((rho - rho0) / sr) ** 2 + ((gamma - g0) / sg) ** 2))

    report = {
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "result_kind": "SCRATCH_MEMORY_REGULARIZED_STABILITY__EVIDENCE_NOT_THEOREM",
        "scenarios": rows,
        "global_summary": {
            "max_rho_cv": float(np.max(rho_cv_all)),
            "max_gamma_cv": float(np.max(gamma_cv_all)),
            "stability_flag_strict_threshold": stable_flag,
        },
        "regularization_note": {
            "type": "quadratic_prior_penalty",
            "intent": "degeneracy control only, not physical proof",
        },
        "symbolic_trace": {
            "regularization_term": sp.sstr(reg_expr),
        },
        "hard_limits": [
            "Regularization can suppress variance without proving mechanism truth.",
            "No theorem-level bridge is claimed.",
            "Legacy physical-role transfer to strict remains blocked.",
        ],
        "next_honest_step": "Compare regularized vs unregularized out-of-window prediction error; if regularization only shrinks variance without predictive gain, classify as non-bridge fit stabilization artifact.",
    }

    OUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
