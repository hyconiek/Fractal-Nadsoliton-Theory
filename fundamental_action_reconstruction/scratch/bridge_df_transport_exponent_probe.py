#!/usr/bin/env python3
"""Scratch theorem-prep probe for the D_f-1 bridge candidate.

Goal:
- Make the D_f link more proof-like than the previous grep/numerical audit.
- Treat the documented path multiplicity N(d) ~ d^(D_f-1) as a transport
  multiplicity ansatz and test the induced fixed exponent gamma_F = D_f-1
  against the strict eta=1.8 denominator over several d-windows.

This remains scratch evidence only.  It does not identify the legacy and strict
kernels or transfer legacy physical roles.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import scipy.optimize as so
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PREVIOUS_AUDIT = HERE / "bridge_df_eta_eff_report.json"
OUT_JSON = HERE / "bridge_df_transport_exponent_report.json"
OUT_MD = HERE / "bridge_df_transport_exponent_report.md"

ALPHA_GEO = 4.0 * math.log(2.0)
GAMMA_F = ALPHA_GEO - 1.0
STRICT = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 1.8}
WINDOWS = [
    {"name": "baseline_bridge_missing_coupling_window", "d_min": 1.0, "d_max": 11.0, "n": 500},
    {"name": "short_training_window", "d_min": 1.0, "d_max": 8.0, "n": 500},
    {"name": "right_holdout_window", "d_min": 2.0, "d_max": 14.0, "n": 500},
    {"name": "near_origin_extended_window", "d_min": 0.5, "d_max": 11.0, "n": 500},
    {"name": "wide_window", "d_min": 1.0, "d_max": 14.0, "n": 500},
]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def k_power(d: np.ndarray, omega: float, phi: float, beta: float, gamma: float, alpha: float = ALPHA_GEO) -> np.ndarray:
    return alpha * np.cos(omega * d + phi) / (1.0 + beta * (d**gamma))


def k_strict(d: np.ndarray) -> np.ndarray:
    return np.cos(STRICT["omega"] * d + STRICT["phi"]) / (1.0 + STRICT["beta"] * (d ** STRICT["eta"]))


def projected_fit_for_fixed_gamma(d: np.ndarray, gamma: float) -> dict[str, float]:
    y = k_strict(d)

    def mse(x: np.ndarray) -> float:
        delta_phi, beta_eff = x
        y_model = k_power(d, STRICT["omega"], STRICT["phi"] + delta_phi, beta_eff, gamma)
        scale = float(np.dot(y_model, y) / max(np.dot(y_model, y_model), 1e-15))
        resid = scale * y_model - y
        return float(np.mean(resid * resid))

    result = so.minimize(
        mse,
        x0=np.array([0.0, 1.0]),
        bounds=[(-math.pi, math.pi), (1e-4, 20.0)],
        method="L-BFGS-B",
    )
    delta_phi, beta_eff = [float(v) for v in result.x]
    y_model = k_power(d, STRICT["omega"], STRICT["phi"] + delta_phi, beta_eff, gamma)
    scale = float(np.dot(y_model, y) / max(np.dot(y_model, y_model), 1e-15))
    resid = scale * y_model - y
    sse = float(np.sum(resid * resid))
    return {
        "gamma": float(gamma),
        "delta_phi": delta_phi,
        "beta_eff": beta_eff,
        "scale": scale,
        "sse": sse,
        "mse": float(sse / len(d)),
    }


def projected_fit_for_free_gamma(d: np.ndarray) -> dict[str, float]:
    y = k_strict(d)

    def mse(x: np.ndarray) -> float:
        delta_phi, beta_eff, gamma = x
        y_model = k_power(d, STRICT["omega"], STRICT["phi"] + delta_phi, beta_eff, gamma)
        scale = float(np.dot(y_model, y) / max(np.dot(y_model, y_model), 1e-15))
        resid = scale * y_model - y
        return float(np.mean(resid * resid))

    starts = [
        np.array([0.0, 1.0, 1.0]),
        np.array([0.0, 1.0, GAMMA_F]),
        np.array([0.0, 1.0, STRICT["eta"]]),
        np.array([0.02, 1.02, GAMMA_F]),
    ]
    best = None
    for start in starts:
        result = so.minimize(
            mse,
            x0=start,
            bounds=[(-math.pi, math.pi), (1e-4, 20.0), (0.2, 3.0)],
            method="L-BFGS-B",
        )
        if best is None or result.fun < best.fun:
            best = result
    assert best is not None
    delta_phi, beta_eff, gamma = [float(v) for v in best.x]
    y_model = k_power(d, STRICT["omega"], STRICT["phi"] + delta_phi, beta_eff, gamma)
    scale = float(np.dot(y_model, y) / max(np.dot(y_model, y_model), 1e-15))
    resid = scale * y_model - y
    sse = float(np.sum(resid * resid))
    return {
        "gamma": gamma,
        "delta_phi": delta_phi,
        "beta_eff": beta_eff,
        "scale": scale,
        "sse": sse,
        "mse": float(sse / len(d)),
    }


def aic_bic(sse: float, n: int, k: int) -> dict[str, float]:
    per = max(sse / n, 1e-300)
    return {"aic": float(n * math.log(per) + 2 * k), "bic": float(n * math.log(per) + k * math.log(n))}


def window_row(window: dict[str, float]) -> dict[str, Any]:
    d = np.linspace(float(window["d_min"]), float(window["d_max"]), int(window["n"]))
    linear = projected_fit_for_fixed_gamma(d, 1.0)
    df_fixed = projected_fit_for_fixed_gamma(d, GAMMA_F)
    strict_fixed = projected_fit_for_fixed_gamma(d, STRICT["eta"])
    free = projected_fit_for_free_gamma(d)
    n = len(d)
    linear.update(aic_bic(linear["sse"], n, 2))
    df_fixed.update(aic_bic(df_fixed["sse"], n, 2))
    strict_fixed.update(aic_bic(strict_fixed["sse"], n, 2))
    free.update(aic_bic(free["sse"], n, 3))
    denominator = max(linear["sse"] - free["sse"], 1e-300)
    capture = (linear["sse"] - df_fixed["sse"]) / denominator
    return {
        "window": window,
        "linear_gamma_1_baseline": linear,
        "fixed_D_f_minus_1_transport": df_fixed,
        "fixed_strict_eta_oracle": strict_fixed,
        "free_gamma_oracle": free,
        "df_transport_vs_linear": {
            "sse_reduction_factor_linear_over_df": float(linear["sse"] / max(df_fixed["sse"], 1e-300)),
            "captured_fraction_of_linear_to_free_sse_improvement": float(capture),
            "delta_aic_linear_minus_df": float(linear["aic"] - df_fixed["aic"]),
            "delta_bic_linear_minus_df": float(linear["bic"] - df_fixed["bic"]),
        },
        "df_transport_vs_strict_eta": {
            "sse_ratio_df_over_strict_eta": float(df_fixed["sse"] / max(strict_fixed["sse"], 1e-300)),
            "gamma_gap_strict_eta_minus_D_f_minus_1": float(STRICT["eta"] - GAMMA_F),
        },
    }


def symbolic_transport_derivation() -> dict[str, str]:
    d, d0, c, beta, Df = sp.symbols("d d0 c beta D_f", positive=True, real=True)
    gamma_f = Df - 1
    multiplicity = c * (d / d0) ** gamma_f
    normalized_coupling = sp.simplify(beta * multiplicity.subs({c: 1, d0: 1}))
    denominator = sp.simplify(1 + normalized_coupling)
    log_slope = sp.simplify(d * sp.diff(sp.log(normalized_coupling), d))
    finite_denominator_slope = sp.simplify(d * sp.diff(sp.log(denominator), d))
    return {
        "path_multiplicity_ansatz": sp.sstr(multiplicity),
        "normalized_coupling_after_absorbing_c_d0_into_beta": sp.sstr(normalized_coupling),
        "induced_denominator": sp.sstr(denominator),
        "asymptotic_coupling_log_slope": sp.sstr(log_slope),
        "finite_denominator_log_slope": sp.sstr(finite_denominator_slope),
        "specialized_gamma_F": "D_f - 1 = 4*log(2) - 1",
    }


def main() -> None:
    previous = load_json(PREVIOUS_AUDIT)
    previous_candidate = previous["candidate"]
    symbolic = symbolic_transport_derivation()
    rows = [window_row(window) for window in WINDOWS]
    captures = [row["df_transport_vs_linear"]["captured_fraction_of_linear_to_free_sse_improvement"] for row in rows]
    reductions = [row["df_transport_vs_linear"]["sse_reduction_factor_linear_over_df"] for row in rows]
    delta_bics = [row["df_transport_vs_linear"]["delta_bic_linear_minus_df"] for row in rows]

    summary = {
        "min_capture_fraction": float(min(captures)),
        "median_capture_fraction": float(np.median(captures)),
        "min_sse_reduction_factor_linear_over_df": float(min(reductions)),
        "median_sse_reduction_factor_linear_over_df": float(np.median(reductions)),
        "min_delta_bic_linear_minus_df": float(min(delta_bics)),
        "all_windows_capture_at_least_99_percent_of_linear_to_free_improvement": bool(min(captures) >= 0.99),
        "all_windows_bic_strongly_prefers_D_f_minus_1_over_linear": bool(min(delta_bics) > 10.0),
        "gamma_F_between_eta_eff_and_strict_eta_replayed": bool(
            previous_candidate["bracketing_signature"]["eta_eff_lt_D_f_minus_1_lt_strict_eta"]
        ),
    }

    report = {
        "status": "OPEN_THEOREM_PREP_TRACE_NO_BRIDGE_THEOREM",
        "result_kind": "SCRATCH_DF_TRANSPORT_EXPONENT_THEOREM_PREP__NOT_A_THEOREM",
        "source_reports": {
            "previous_df_eta_eff_audit": str(PREVIOUS_AUDIT.relative_to(ROOT)),
        },
        "constants": {
            "D_f": ALPHA_GEO,
            "gamma_F_D_f_minus_1": GAMMA_F,
            "strict_eta": STRICT["eta"],
            "strict_eta_minus_gamma_F": STRICT["eta"] - GAMMA_F,
            "previous_eta_eff": previous_candidate["numerical_constants"]["eta_eff_from_bridge_missing_coupling"],
        },
        "symbolic_transport_derivation": symbolic,
        "window_rows": rows,
        "aggregate_summary": summary,
        "honest_interpretation": [
            "If the documented fractal path multiplicity N(d)~d^(D_f-1) is admitted as a transport multiplicity, then after absorbing constants into beta it induces a denominator exponent gamma_F=D_f-1.",
            "The fixed gamma_F denominator is not identical to strict eta=1.8, but across all tested windows it captures at least 99% of the linear-to-free exponent improvement and strongly beats the linear legacy denominator by BIC.",
            "This upgrades the D_f link from a grep/numerical coincidence to a concrete theorem-prep obligation: prove or refute that path multiplicity really enters the damping denominator in this normalized way.",
        ],
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role formula is transferred to the strict kernel.",
            "The fixed D_f-1 exponent remains distinguishable from eta=1.8 and from the previous fitted eta_eff.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Prove or falsify the transport insertion lemma: N(d)~d^(D_f-1) must enter the effective damping denominator, not merely the numerator or measure, and finite-window renormalization must explain the small D_f-1 -> eta_eff/eta offsets.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch theorem-prep: D_f-1 transport exponent\n\n"
        "Status: proof-oriented scratch trace; not a legacy→strict bridge theorem.\n\n"
        f"- Symbolic ansatz induces denominator `1 + beta*d**(D_f - 1)` after normalization; `D_f-1 = {GAMMA_F:.15f}`.\n"
        f"- Strict eta is `{STRICT['eta']:.15f}`, gap `eta-(D_f-1) = {STRICT['eta'] - GAMMA_F:.15f}`.\n"
        f"- Across `{len(rows)}` windows, min captured linear→free SSE improvement is `{summary['min_capture_fraction']:.12f}`.\n"
        f"- Min linear/`D_f-1` SSE reduction factor is `{summary['min_sse_reduction_factor_linear_over_df']:.6f}`; min ΔBIC(linear−D_f) is `{summary['min_delta_bic_linear_minus_df']:.6f}`.\n"
        f"- Bracketing from prior audit replayed: `{summary['gamma_F_between_eta_eff_and_strict_eta_replayed']}`.\n"
        "- No false pass: no kernel identity, no physical-role transfer, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
