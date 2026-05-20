#!/usr/bin/env python3
"""Scratch probe: does matter-history memory improve legacy->strict bridge fit?

Guardrail discipline:
- exploratory model-selection only,
- no theorem/closure claim,
- no role-transfer claim.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import scipy.optimize as so
import scipy.stats as ss
import sympy as sp

HERE = Path(__file__).resolve().parent
OUT = HERE / "bridge_matter_history_memory_report.json"


def k_legacy(d: np.ndarray, alpha_geo: float, omega: float, phi: float, beta_tors: float) -> np.ndarray:
    return alpha_geo * np.cos(omega * d + phi) / (1.0 + beta_tors * d)


def k_strict(d: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * (d**eta))


def exp_memory_filter(x: np.ndarray, rho: float) -> np.ndarray:
    """Causal memory accumulator: y_t = (1-rho) x_t + rho y_{t-1}."""
    y = np.empty_like(x)
    y[0] = x[0]
    a = 1.0 - rho
    for i in range(1, x.size):
        y[i] = a * x[i] + rho * y[i - 1]
    return y


def main() -> None:
    d = np.linspace(1.0, 11.0, 700)
    omega, phi, beta, eta = 0.18575, 0.16250, 1.0, 1.8
    alpha_geo = float(4.0 * np.log(2.0))
    y = k_strict(d, omega, phi, beta, eta)
    n = len(d)

    # M0: no memory, but with phase/beta_tors adjustments
    def mse_m0(x: np.ndarray) -> float:
        dphi, bt = x
        y0 = k_legacy(d, alpha_geo, omega, phi + dphi, bt)
        c = float(np.dot(y0, y) / max(np.dot(y0, y0), 1e-15))
        return float(np.mean((c * y0 - y) ** 2))

    fit0 = so.minimize(mse_m0, x0=np.array([0.0, 1.0]), bounds=[(-np.pi, np.pi), (1e-4, 30.0)], method="L-BFGS-B")
    dphi0, bt0 = [float(v) for v in fit0.x]
    raw0 = k_legacy(d, alpha_geo, omega, phi + dphi0, bt0)
    c0 = float(np.dot(raw0, y) / max(np.dot(raw0, raw0), 1e-15))
    pred0 = c0 * raw0

    # M1: memory-augmented legacy with one causal history channel
    # pred = c * raw + gamma * memory(raw; rho)
    def mse_m1(x: np.ndarray) -> float:
        dphi, bt, rho, gamma = x
        raw = k_legacy(d, alpha_geo, omega, phi + dphi, bt)
        mem = exp_memory_filter(raw, rho)
        feat = np.vstack([raw, gamma * mem]).T
        coef, *_ = np.linalg.lstsq(feat, y, rcond=None)
        yhat = feat @ coef
        return float(np.mean((yhat - y) ** 2))

    fit1 = so.minimize(
        mse_m1,
        x0=np.array([0.0, 1.0, 0.35, 0.5]),
        bounds=[(-np.pi, np.pi), (1e-4, 30.0), (1e-5, 0.999), (-5.0, 5.0)],
        method="L-BFGS-B",
    )
    dphi1, bt1, rho1, gamma1 = [float(v) for v in fit1.x]
    raw1 = k_legacy(d, alpha_geo, omega, phi + dphi1, bt1)
    mem1 = exp_memory_filter(raw1, rho1)
    feat1 = np.vstack([raw1, gamma1 * mem1]).T
    coef1, *_ = np.linalg.lstsq(feat1, y, rcond=None)
    pred1 = feat1 @ coef1

    r0 = pred0 - y
    r1 = pred1 - y
    sse0 = float(np.sum(r0 * r0))
    sse1 = float(np.sum(r1 * r1))

    # effective parameter counts: M0=2 + one linear scale, M1=4 + two linear scales
    k0, k1 = 3, 6
    aic0 = float(n * np.log(max(sse0 / n, 1e-15)) + 2 * k0)
    aic1 = float(n * np.log(max(sse1 / n, 1e-15)) + 2 * k1)
    bic0 = float(n * np.log(max(sse0 / n, 1e-15)) + k0 * np.log(n))
    bic1 = float(n * np.log(max(sse1 / n, 1e-15)) + k1 * np.log(n))

    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2_0 = float(1.0 - sse0 / max(ss_tot, 1e-15))
    r2_1 = float(1.0 - sse1 / max(ss_tot, 1e-15))

    rng = np.random.default_rng(20260519)
    dsse = []
    for _ in range(400):
        idx = rng.integers(0, n, size=n)
        dsse.append(float(np.sum(r0[idx] * r0[idx]) - np.sum(r1[idx] * r1[idx])))
    dsse = np.array(dsse, dtype=float)
    p_sup = float(np.mean(dsse > 0.0))
    p_w = float(ss.wilcoxon(dsse, alternative="greater", zero_method="zsplit", method="approx").pvalue) if np.any(np.abs(dsse) > 1e-15) else 1.0

    # Symbolic: memory term is a true extension unless gamma=0 or rho=0 degenerate case
    g, r, x_t, y_tm1 = sp.symbols("g r x_t y_tm1", real=True)
    mem_update = sp.expand((1 - r) * x_t + r * y_tm1)
    extra_term = sp.expand(g * mem_update)

    report = {
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "result_kind": "SCRATCH_MATTER_HISTORY_MEMORY_PROBE__EVIDENCE_NOT_THEOREM",
        "models": {
            "M0_no_memory": {
                "delta_phi": dphi0,
                "beta_tors": bt0,
                "r2": r2_0,
                "sse": sse0,
                "aic": aic0,
                "bic": bic0,
            },
            "M1_memory_augmented": {
                "delta_phi": dphi1,
                "beta_tors": bt1,
                "rho_memory": rho1,
                "gamma_memory": gamma1,
                "linear_coef_raw": float(coef1[0]),
                "linear_coef_memory_channel": float(coef1[1]),
                "r2": r2_1,
                "sse": sse1,
                "aic": aic1,
                "bic": bic1,
            },
        },
        "comparison": {
            "delta_r2_M1_minus_M0": r2_1 - r2_0,
            "delta_aic_M0_minus_M1": aic0 - aic1,
            "delta_bic_M0_minus_M1": bic0 - bic1,
            "bootstrap_prob_M1_better": p_sup,
            "wilcoxon_pvalue_M1_better": p_w,
        },
        "symbolic_trace": {
            "memory_update": sp.sstr(mem_update),
            "extra_term": sp.sstr(extra_term),
            "note": "Memory extension collapses to baseline only on degenerate submanifolds (e.g., gamma=0).",
        },
        "hard_limits": [
            "A better fit from memory channel does not prove physical causality.",
            "No legacy->strict role transfer is licensed by this probe.",
            "Bridge theorem remains open unless identifiability and mechanism export are proven.",
        ],
        "next_honest_step": "Run cross-observable validation: check if (rho_memory, gamma_memory) are stable across independent channels; unstable values classify this as fit-level artifact.",
    }

    OUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
