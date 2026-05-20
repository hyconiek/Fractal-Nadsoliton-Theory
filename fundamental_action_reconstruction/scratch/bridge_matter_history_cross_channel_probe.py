#!/usr/bin/env python3
"""Scratch probe: cross-channel stability for matter-history memory hypothesis.

This is a strict, non-theorem follow-up:
- evaluates whether memory parameters remain stable across channel variants,
- keeps OPEN_OBSTRUCTION status unless stability is demonstrated.
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
OUT = HERE / "bridge_matter_history_cross_channel_report.json"


@dataclass
class ChannelSpec:
    name: str
    omega: float
    phi: float
    beta: float
    eta: float


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


def fit_models(d: np.ndarray, y: np.ndarray, alpha_geo: float, omega: float, phi: float) -> dict:
    n = len(d)

    def mse_m0(x: np.ndarray) -> float:
        dphi, bt = x
        raw = k_legacy(d, alpha_geo, omega, phi + dphi, bt)
        c = float(np.dot(raw, y) / max(np.dot(raw, raw), 1e-15))
        return float(np.mean((c * raw - y) ** 2))

    m0 = so.minimize(mse_m0, x0=np.array([0.0, 1.0]), bounds=[(-np.pi, np.pi), (1e-5, 60.0)], method="L-BFGS-B")
    dphi0, bt0 = [float(v) for v in m0.x]
    raw0 = k_legacy(d, alpha_geo, omega, phi + dphi0, bt0)
    c0 = float(np.dot(raw0, y) / max(np.dot(raw0, raw0), 1e-15))
    y0 = c0 * raw0
    r0 = y0 - y
    sse0 = float(np.sum(r0 * r0))

    def mse_m1(x: np.ndarray) -> float:
        dphi, bt, rho, gamma = x
        raw = k_legacy(d, alpha_geo, omega, phi + dphi, bt)
        mem = exp_memory_filter(raw, rho)
        feat = np.vstack([raw, gamma * mem]).T
        coef, *_ = np.linalg.lstsq(feat, y, rcond=None)
        yhat = feat @ coef
        return float(np.mean((yhat - y) ** 2))

    m1 = so.minimize(
        mse_m1,
        x0=np.array([0.0, 1.0, 0.25, 0.5]),
        bounds=[(-np.pi, np.pi), (1e-5, 60.0), (1e-5, 0.999), (-8.0, 8.0)],
        method="L-BFGS-B",
    )
    dphi1, bt1, rho1, gamma1 = [float(v) for v in m1.x]
    raw1 = k_legacy(d, alpha_geo, omega, phi + dphi1, bt1)
    mem1 = exp_memory_filter(raw1, rho1)
    feat1 = np.vstack([raw1, gamma1 * mem1]).T
    coef1, *_ = np.linalg.lstsq(feat1, y, rcond=None)
    y1 = feat1 @ coef1
    r1 = y1 - y
    sse1 = float(np.sum(r1 * r1))

    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2_0 = float(1.0 - sse0 / max(ss_tot, 1e-15))
    r2_1 = float(1.0 - sse1 / max(ss_tot, 1e-15))

    k0, k1 = 3, 6
    aic0 = float(n * np.log(max(sse0 / n, 1e-15)) + 2 * k0)
    aic1 = float(n * np.log(max(sse1 / n, 1e-15)) + 2 * k1)

    bic0 = float(n * np.log(max(sse0 / n, 1e-15)) + k0 * np.log(n))
    bic1 = float(n * np.log(max(sse1 / n, 1e-15)) + k1 * np.log(n))

    return {
        "M0": {"delta_phi": dphi0, "beta_tors": bt0, "r2": r2_0, "sse": sse0, "aic": aic0, "bic": bic0},
        "M1": {
            "delta_phi": dphi1,
            "beta_tors": bt1,
            "rho_memory": rho1,
            "gamma_memory": gamma1,
            "coef_raw": float(coef1[0]),
            "coef_memory": float(coef1[1]),
            "r2": r2_1,
            "sse": sse1,
            "aic": aic1,
            "bic": bic1,
        },
        "delta": {
            "delta_r2_M1_minus_M0": r2_1 - r2_0,
            "delta_aic_M0_minus_M1": aic0 - aic1,
            "delta_bic_M0_minus_M1": bic0 - bic1,
        },
    }


def main() -> None:
    alpha_geo = float(4.0 * np.log(2.0))
    d = np.linspace(1.0, 11.0, 700)

    channels = [
        ChannelSpec("strict_baseline", 0.18575, 0.16250, 1.0, 1.8),
        ChannelSpec("strict_phase_plus", 0.18575, 0.21250, 1.0, 1.8),
        ChannelSpec("strict_phase_minus", 0.18575, 0.11250, 1.0, 1.8),
        ChannelSpec("strict_eta_low", 0.18575, 0.16250, 1.0, 1.65),
        ChannelSpec("strict_eta_high", 0.18575, 0.16250, 1.0, 1.95),
    ]

    rows = []
    for ch in channels:
        y = k_strict(d, ch.omega, ch.phi, ch.beta, ch.eta)
        fit = fit_models(d, y, alpha_geo, ch.omega, ch.phi)
        rows.append({"channel": asdict(ch), **fit})

    rho_vals = np.array([r["M1"]["rho_memory"] for r in rows], dtype=float)
    gamma_vals = np.array([r["M1"]["gamma_memory"] for r in rows], dtype=float)
    daic_vals = np.array([r["delta"]["delta_aic_M0_minus_M1"] for r in rows], dtype=float)

    rho_cv = float(np.std(rho_vals) / max(abs(np.mean(rho_vals)), 1e-12))
    gamma_cv = float(np.std(gamma_vals) / max(abs(np.mean(gamma_vals)), 1e-12))

    stability_flag = bool((rho_cv < 0.35) and (gamma_cv < 0.60) and np.all(daic_vals > 10.0))

    # symbolic sanity: memory extension is nonlocal in index (depends on y_{t-1})
    y_t, y_tm1, x_t, rho = sp.symbols("y_t y_tm1 x_t rho", real=True)
    recurrence = sp.Eq(y_t, (1 - rho) * x_t + rho * y_tm1)

    report = {
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "result_kind": "SCRATCH_CROSS_CHANNEL_MEMORY_STABILITY__EVIDENCE_NOT_THEOREM",
        "channels_evaluated": rows,
        "cross_channel_summary": {
            "rho_memory_values": [float(v) for v in rho_vals],
            "gamma_memory_values": [float(v) for v in gamma_vals],
            "delta_aic_values": [float(v) for v in daic_vals],
            "rho_cv": rho_cv,
            "gamma_cv": gamma_cv,
            "stability_flag_strict_threshold": stability_flag,
        },
        "symbolic_trace": {
            "memory_recurrence": sp.sstr(recurrence),
            "note": "Memory hypothesis introduces a causal recurrence, not just a static algebraic re-scaling.",
        },
        "hard_limits": [
            "Cross-channel stability improves plausibility but is still not a bridge theorem.",
            "No automatic legacy->strict physical-role transfer is licensed.",
            "Without identifiability + exported mechanism theorem, status remains open.",
        ],
        "next_honest_step": "Add seed-ensemble + d-range perturbation per channel and test whether rho/gamma remain stable jointly; only then draft explicit assumption-level theorem candidate.",
    }

    OUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
