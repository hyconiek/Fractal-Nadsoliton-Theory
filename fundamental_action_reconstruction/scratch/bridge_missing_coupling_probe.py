#!/usr/bin/env python3
"""Scratch probe: does strict behavior suggest a missing legacy coupling/property?

Hypothesis tested:
- Legacy may miss one nadsoliton property (e.g., fractional/fractal transport coupling).
- Add one extra coupling exponent to legacy denominator and compare against strict target.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import scipy.optimize as so
import scipy.stats as ss
import sympy as sp

HERE = Path(__file__).resolve().parent
OUT = HERE / "bridge_missing_coupling_report.json"


def k_legacy(d: np.ndarray, alpha_geo: float, omega: float, phi: float, beta_tors: float) -> np.ndarray:
    return alpha_geo * np.cos(omega * d + phi) / (1.0 + beta_tors * d)


def k_legacy_ext(d: np.ndarray, alpha_geo: float, omega: float, phi: float, beta_tors: float, eta_eff: float) -> np.ndarray:
    return alpha_geo * np.cos(omega * d + phi) / (1.0 + beta_tors * (d ** eta_eff))


def k_strict(d: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * (d ** eta))


def main() -> None:
    d = np.linspace(1.0, 11.0, 500)
    omega, phi, beta, eta = 0.18575, 0.16250, 1.0, 1.8
    alpha_geo = float(4.0 * np.log(2.0))
    y = k_strict(d, omega, phi, beta, eta)
    n = len(d)

    # M0: legacy-only (no extra coupling exponent)
    def mse_m0(x: np.ndarray) -> float:
        dphi, bt = x
        y0 = k_legacy(d, alpha_geo, omega, phi + dphi, bt)
        c = float(np.dot(y0, y) / max(np.dot(y0, y0), 1e-15))
        return float(np.mean((c * y0 - y) ** 2))

    m0 = so.minimize(mse_m0, x0=np.array([0.0, 1.0]), bounds=[(-np.pi, np.pi), (1e-4, 20.0)], method="L-BFGS-B")
    dphi0, bt0 = [float(v) for v in m0.x]
    y0 = k_legacy(d, alpha_geo, omega, phi + dphi0, bt0)
    c0 = float(np.dot(y0, y) / max(np.dot(y0, y0), 1e-15))
    r0 = c0 * y0 - y
    sse0 = float(np.sum(r0 * r0))

    # M1: legacy + one additional property/coupling (eta_eff)
    def mse_m1(x: np.ndarray) -> float:
        dphi, bt, eta_eff = x
        y1 = k_legacy_ext(d, alpha_geo, omega, phi + dphi, bt, eta_eff)
        c = float(np.dot(y1, y) / max(np.dot(y1, y1), 1e-15))
        return float(np.mean((c * y1 - y) ** 2))

    m1 = so.minimize(mse_m1, x0=np.array([0.0, 1.0, 1.0]), bounds=[(-np.pi, np.pi), (1e-4, 20.0), (0.2, 3.0)], method="L-BFGS-B")
    dphi1, bt1, eta1 = [float(v) for v in m1.x]
    y1 = k_legacy_ext(d, alpha_geo, omega, phi + dphi1, bt1, eta1)
    c1 = float(np.dot(y1, y) / max(np.dot(y1, y1), 1e-15))
    r1 = c1 * y1 - y
    sse1 = float(np.sum(r1 * r1))

    # information criteria (k excludes scale c because projected analytically both models)
    k0, k1 = 2, 3
    aic0 = float(n * np.log(max(sse0 / n, 1e-15)) + 2 * k0)
    aic1 = float(n * np.log(max(sse1 / n, 1e-15)) + 2 * k1)
    bic0 = float(n * np.log(max(sse0 / n, 1e-15)) + k0 * np.log(n))
    bic1 = float(n * np.log(max(sse1 / n, 1e-15)) + k1 * np.log(n))
    delta_aic = aic0 - aic1
    delta_bic = bic0 - bic1

    # bootstrap support for delta_sse > 0
    rng = np.random.default_rng(20260519)
    dsse = []
    for _ in range(300):
        idx = rng.integers(0, n, size=n)
        dsse.append(float(np.sum(r0[idx] * r0[idx]) - np.sum(r1[idx] * r1[idx])))
    dsse_arr = np.array(dsse, dtype=float)
    p_superior = float(np.mean(dsse_arr > 0.0))
    pval = float(ss.wilcoxon(dsse_arr, alternative="greater", zero_method="zsplit", method="approx").pvalue) if np.any(np.abs(dsse_arr) > 1e-15) else 1.0

    # symbolic statement: adding eta_eff creates distinct model class unless eta_eff=1
    d_sym, b_sym, eta_sym = sp.symbols("d b eta", positive=True, real=True)
    denom_diff = sp.simplify((1 + b_sym * d_sym) - (1 + b_sym * d_sym**eta_sym))

    report = {
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "result_kind": "SCRATCH_MISSING_COUPLING_HYPOTHESIS_TEST__NOT_A_THEOREM",
        "models": {
            "M0_legacy": {"delta_phi": dphi0, "beta_tors": bt0, "scale_c": c0, "sse": sse0, "aic": aic0, "bic": bic0},
            "M1_legacy_plus_eta_eff": {"delta_phi": dphi1, "beta_tors": bt1, "eta_eff": eta1, "scale_c": c1, "sse": sse1, "aic": aic1, "bic": bic1},
        },
        "comparison": {
            "delta_aic_M0_minus_M1": delta_aic,
            "delta_bic_M0_minus_M1": delta_bic,
            "bootstrap_prob_M1_better": p_superior,
            "wilcoxon_pvalue_M1_better": pval,
        },
        "symbolic_class_separation": {
            "denominator_difference_symbolic": sp.sstr(denom_diff),
            "note": "Model classes coincide only on constrained submanifolds (e.g., eta_eff=1).",
        },
        "hard_limits": [
            "Evidence for missing coupling is model-selection evidence, not a bridge theorem.",
            "No role-transfer from legacy ontology to strict kernel is inferred.",
            "Kernel split remains open pending explicit theorem/non-bridge theorem.",
        ],
        "next_honest_step": "Test whether eta_eff remains stable across independent observables/channels; if unstable, classify as phenomenological fit artifact rather than missing fundamental coupling.",
    }
    OUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()

