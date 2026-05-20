#!/usr/bin/env python3
"""Scratch probe: out-of-window predictive check (regularized vs unregularized memory model)."""
from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import scipy.optimize as so
import scipy.stats as ss
import sympy as sp

HERE = Path(__file__).resolve().parent
OUT = HERE / "bridge_matter_history_out_of_window_report.json"


@dataclass
class WindowCase:
    name: str
    train_min: float
    train_max: float
    hold_min: float
    hold_max: float
    n_train: int
    n_hold: int
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


def fit_unregularized(train_d: np.ndarray, y_train: np.ndarray, alpha_geo: float, omega: float, phi: float, seed: int) -> dict:
    rng = np.random.default_rng(seed)

    def loss(x: np.ndarray) -> float:
        dphi, bt, rho, gamma = x
        raw = k_legacy(train_d, alpha_geo, omega, phi + dphi, bt)
        mem = exp_memory_filter(raw, rho)
        feat = np.vstack([raw, gamma * mem]).T
        coef, *_ = np.linalg.lstsq(feat, y_train, rcond=None)
        pred = feat @ coef
        return float(np.mean((pred - y_train) ** 2))

    best = None
    for _ in range(10):
        x0 = np.array([rng.uniform(-np.pi, np.pi), rng.uniform(1e-4, 40.0), rng.uniform(1e-5, 0.95), rng.uniform(-8, 8)])
        res = so.minimize(loss, x0=x0, bounds=[(-np.pi, np.pi), (1e-5, 60.0), (1e-5, 0.999), (-10.0, 10.0)], method="L-BFGS-B")
        if best is None or res.fun < best.fun:
            best = res
    dphi, bt, rho, gamma = [float(v) for v in best.x]
    return {"delta_phi": dphi, "beta_tors": bt, "rho": rho, "gamma": gamma}


def fit_regularized(train_d: np.ndarray, y_train: np.ndarray, alpha_geo: float, omega: float, phi: float, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    rho_center, rho_scale = 0.22, 0.08
    gamma_center, gamma_scale = 1.0, 1.2
    lam = 3e-3

    def loss(x: np.ndarray) -> float:
        dphi, bt, rho, gamma = x
        raw = k_legacy(train_d, alpha_geo, omega, phi + dphi, bt)
        mem = exp_memory_filter(raw, rho)
        feat = np.vstack([raw, gamma * mem]).T
        coef, *_ = np.linalg.lstsq(feat, y_train, rcond=None)
        pred = feat @ coef
        mse = float(np.mean((pred - y_train) ** 2))
        prior = ((rho - rho_center) / rho_scale) ** 2 + ((gamma - gamma_center) / gamma_scale) ** 2
        return mse + lam * prior

    best = None
    for _ in range(10):
        x0 = np.array([rng.uniform(-0.4, 0.4), rng.uniform(1e-4, 30.0), rng.uniform(0.05, 0.45), rng.uniform(-1.5, 3.5)])
        res = so.minimize(loss, x0=x0, bounds=[(-np.pi, np.pi), (1e-5, 40.0), (0.01, 0.60), (-2.0, 4.0)], method="L-BFGS-B")
        if best is None or res.fun < best.fun:
            best = res
    dphi, bt, rho, gamma = [float(v) for v in best.x]
    return {"delta_phi": dphi, "beta_tors": bt, "rho": rho, "gamma": gamma}


def predict(d: np.ndarray, params: dict, alpha_geo: float, omega: float, phi: float) -> np.ndarray:
    raw = k_legacy(d, alpha_geo, omega, phi + params["delta_phi"], params["beta_tors"])
    mem = exp_memory_filter(raw, params["rho"])
    feat = np.vstack([raw, params["gamma"] * mem]).T
    return feat @ np.linalg.lstsq(feat, np.zeros_like(d), rcond=None)[0]


def mse_with_refit_linear(d: np.ndarray, y: np.ndarray, params: dict, alpha_geo: float, omega: float, phi: float) -> float:
    raw = k_legacy(d, alpha_geo, omega, phi + params["delta_phi"], params["beta_tors"])
    mem = exp_memory_filter(raw, params["rho"])
    feat = np.vstack([raw, params["gamma"] * mem]).T
    coef, *_ = np.linalg.lstsq(feat, y, rcond=None)
    pred = feat @ coef
    return float(np.mean((pred - y) ** 2))


def main() -> None:
    alpha_geo = float(4.0 * np.log(2.0))
    omega0, phi0, beta0, eta0 = 0.18575, 0.16250, 1.0, 1.8

    cases = [
        WindowCase("expand_right", 1.0, 8.0, 8.0, 14.0, 500, 400, 0.00, 0.00),
        WindowCase("expand_left", 4.0, 11.0, 1.0, 4.0, 500, 300, 0.00, 0.00),
        WindowCase("phase_eta_plus", 1.0, 8.0, 8.0, 14.0, 500, 400, +0.04, +0.10),
        WindowCase("phase_eta_minus", 4.0, 11.0, 1.0, 4.0, 500, 300, -0.04, -0.10),
    ]

    seeds = [20260519 + i for i in range(8)]
    rows = []
    hold_gain = []
    for c in cases:
        d_tr = np.linspace(c.train_min, c.train_max, c.n_train)
        d_ho = np.linspace(c.hold_min, c.hold_max, c.n_hold)
        y_tr = k_strict(d_tr, omega0, phi0 + c.phase_shift, beta0, eta0 + c.eta_shift)
        y_ho = k_strict(d_ho, omega0, phi0 + c.phase_shift, beta0, eta0 + c.eta_shift)

        case_metrics = []
        for s in seeds:
            pu = fit_unregularized(d_tr, y_tr, alpha_geo, omega0, phi0 + c.phase_shift, s)
            pr = fit_regularized(d_tr, y_tr, alpha_geo, omega0, phi0 + c.phase_shift, s)
            mse_u = mse_with_refit_linear(d_ho, y_ho, pu, alpha_geo, omega0, phi0 + c.phase_shift)
            mse_r = mse_with_refit_linear(d_ho, y_ho, pr, alpha_geo, omega0, phi0 + c.phase_shift)
            case_metrics.append({"seed": s, "holdout_mse_unregularized": mse_u, "holdout_mse_regularized": mse_r, "delta_u_minus_r": mse_u - mse_r})
            hold_gain.append(mse_u - mse_r)

        deltas = np.array([m["delta_u_minus_r"] for m in case_metrics], dtype=float)
        rows.append({
            "case": asdict(c),
            "seed_count": len(seeds),
            "delta_holdout_mse_u_minus_r_q05_q50_q95": [float(v) for v in np.quantile(deltas, [0.05, 0.5, 0.95])],
            "regularized_better_freq": float(np.mean(deltas > 0.0)),
            "metrics": case_metrics,
        })

    hold_gain = np.array(hold_gain, dtype=float)
    p_w = float(ss.wilcoxon(hold_gain, alternative="greater", zero_method="zsplit", method="approx").pvalue) if np.any(np.abs(hold_gain) > 1e-15) else 1.0

    # symbolic ridge-prior style term
    x, x0, lam = sp.symbols("x x0 lam", real=True)
    pen = sp.expand(lam * (x - x0) ** 2)

    report = {
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "result_kind": "SCRATCH_MEMORY_OUT_OF_WINDOW_COMPARISON__EVIDENCE_NOT_THEOREM",
        "cases": rows,
        "global": {
            "delta_holdout_mse_u_minus_r_q05_q50_q95": [float(v) for v in np.quantile(hold_gain, [0.05, 0.5, 0.95])],
            "regularized_better_freq": float(np.mean(hold_gain > 0.0)),
            "wilcoxon_pvalue_regularized_better": p_w,
            "predictive_gain_flag": bool(float(np.mean(hold_gain > 0.0)) > 0.70 and p_w < 0.01),
        },
        "symbolic_trace": {
            "quadratic_penalty_scalar_form": sp.sstr(pen),
        },
        "hard_limits": [
            "Holdout gain strengthens practical plausibility but is not a bridge theorem.",
            "No legacy->strict role transfer is implied.",
            "Mechanism identification remains open without explicit theorem and identifiability export.",
        ],
        "next_honest_step": "If predictive gain is weak/inconsistent, classify as non-bridge stabilization artifact; if robust across additional independent observables, prepare assumption-explicit theorem candidate draft.",
    }

    OUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
