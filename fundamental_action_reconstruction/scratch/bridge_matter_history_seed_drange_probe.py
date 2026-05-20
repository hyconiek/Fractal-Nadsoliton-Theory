#!/usr/bin/env python3
"""Scratch probe: seed-ensemble + d-range perturbation for memory hypothesis.

Goal: next honest step after cross-channel check.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import scipy.optimize as so
import scipy.stats as ss
import sympy as sp

HERE = Path(__file__).resolve().parent
OUT = HERE / "bridge_matter_history_seed_drange_report.json"


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


def fit_m1_multistart(d: np.ndarray, y: np.ndarray, alpha_geo: float, omega: float, phi: float, seed: int, starts: int = 24) -> dict:
    rng = np.random.default_rng(seed)

    def mse(x: np.ndarray) -> float:
        dphi, bt, rho, gamma = x
        raw = k_legacy(d, alpha_geo, omega, phi + dphi, bt)
        mem = exp_memory_filter(raw, rho)
        feat = np.vstack([raw, gamma * mem]).T
        coef, *_ = np.linalg.lstsq(feat, y, rcond=None)
        yhat = feat @ coef
        return float(np.mean((yhat - y) ** 2))

    best = None
    for _ in range(starts):
        x0 = np.array([
            rng.uniform(-np.pi, np.pi),
            rng.uniform(1e-4, 60.0),
            rng.uniform(1e-4, 0.98),
            rng.uniform(-8.0, 8.0),
        ])
        res = so.minimize(mse, x0=x0, bounds=[(-np.pi, np.pi), (1e-5, 60.0), (1e-5, 0.999), (-10.0, 10.0)], method="L-BFGS-B")
        if (best is None) or (res.fun < best.fun):
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
        "mse": float(best.fun),
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

    seed_list = [20260519 + i for i in range(8)]
    rows = []
    for sc in scenarios:
        d = np.linspace(sc.d_min, sc.d_max, sc.npts)
        y = k_strict(d, omega0, phi0 + sc.phase_shift, beta0, eta0 + sc.eta_shift)
        scenario_fits = []
        for sd in seed_list:
            scenario_fits.append(fit_m1_multistart(d, y, alpha_geo, omega0, phi0 + sc.phase_shift, seed=sd, starts=8))
        rho_vals = np.array([r["rho_memory"] for r in scenario_fits], dtype=float)
        gamma_vals = np.array([r["gamma_memory"] for r in scenario_fits], dtype=float)
        r2_vals = np.array([r["r2"] for r in scenario_fits], dtype=float)
        rows.append({
            "scenario": asdict(sc),
            "seed_count": len(seed_list),
            "rho_q05_q50_q95": [float(v) for v in np.quantile(rho_vals, [0.05, 0.5, 0.95])],
            "gamma_q05_q50_q95": [float(v) for v in np.quantile(gamma_vals, [0.05, 0.5, 0.95])],
            "r2_q05_q50_q95": [float(v) for v in np.quantile(r2_vals, [0.05, 0.5, 0.95])],
            "rho_cv": float(np.std(rho_vals) / max(abs(np.mean(rho_vals)), 1e-12)),
            "gamma_cv": float(np.std(gamma_vals) / max(abs(np.mean(gamma_vals)), 1e-12)),
        })

    rho_cv_all = np.array([r["rho_cv"] for r in rows], dtype=float)
    gamma_cv_all = np.array([r["gamma_cv"] for r in rows], dtype=float)
    stable = bool(np.max(rho_cv_all) < 0.35 and np.max(gamma_cv_all) < 0.60)

    # Sympy check: memory recurrence remains first-order causal operator.
    y_t, y_tm1, x_t, rho = sp.symbols("y_t y_tm1 x_t rho", real=True)
    recurrence = sp.Eq(y_t - rho * y_tm1, (1 - rho) * x_t)
    z = sp.symbols("z", complex=True)
    transfer = sp.simplify((1 - rho) / (1 - rho / z))

    report = {
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "result_kind": "SCRATCH_MEMORY_SEED_DRANGE_STABILITY__EVIDENCE_NOT_THEOREM",
        "scenarios": rows,
        "global_summary": {
            "max_rho_cv": float(np.max(rho_cv_all)),
            "max_gamma_cv": float(np.max(gamma_cv_all)),
            "stability_flag_strict_threshold": stable,
        },
        "symbolic_trace": {
            "memory_recurrence_rearranged": sp.sstr(recurrence),
            "z_transfer_form": sp.sstr(transfer),
        },
        "hard_limits": [
            "Even seed-range stability does not identify a unique physical mechanism.",
            "No theorem-level legacy->strict bridge is claimed.",
            "Physical-role transfer remains blocked without explicit bridge theorem.",
        ],
        "next_honest_step": "If instability persists, classify memory hypothesis as non-bridge fit family; if stabilized under tighter priors and independent observables, draft assumption-explicit theorem candidate.",
    }

    OUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
