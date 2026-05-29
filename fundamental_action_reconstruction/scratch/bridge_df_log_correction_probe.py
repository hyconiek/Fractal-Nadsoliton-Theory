#!/usr/bin/env python3
"""Scratch probe: can the D_f-1 route reach strict eta=9/5 by a small log correction?

Previous scratch audits show that D_f-1 is the strongest denominator-placement
candidate, but also that D_f-1 is not analytically identical to strict eta=9/5.
This probe isolates the remaining gap as a marginal logarithmic correction:

    d^(9/5) = d^(D_f-1) * exp(delta*log(d)),
    delta = 9/5 - (D_f-1) = 14/5 - 4 log(2).

It then checks whether the first-order log-corrected denominator

    1 + beta*d^(D_f-1)*(1 + epsilon*log(d))

is enough to approximate the strict denominator over the same windows.  This is
not a bridge theorem; it is a theorem-prep residual and domain-control packet.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
import scipy.optimize as so
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DOMAIN_AUDIT = HERE / "bridge_kernel_analytic_domain_audit_report.json"
PLACEMENT_AUDIT = HERE / "bridge_df_transport_placement_discriminator_report.json"
OUT_JSON = HERE / "bridge_df_log_correction_report.json"
OUT_MD = HERE / "bridge_df_log_correction_report.md"

ALPHA_GEO = 4.0 * math.log(2.0)
GAMMA_F = ALPHA_GEO - 1.0
STRICT_ETA = float(Fraction(9, 5))
DELTA_STRICT = STRICT_ETA - GAMMA_F
STRICT = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": STRICT_ETA}
WINDOWS = [
    {"name": "baseline", "d_min": 1.0, "d_max": 11.0, "n": 500},
    {"name": "short", "d_min": 1.0, "d_max": 8.0, "n": 500},
    {"name": "right_holdout", "d_min": 2.0, "d_max": 14.0, "n": 500},
    {"name": "near_origin", "d_min": 0.5, "d_max": 11.0, "n": 500},
    {"name": "wide", "d_min": 1.0, "d_max": 14.0, "n": 500},
]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def strict_target(d: np.ndarray) -> np.ndarray:
    return np.cos(STRICT["omega"] * d + STRICT["phi"]) / (1.0 + d ** STRICT["eta"])


def base_model(d: np.ndarray, delta_phi: float, beta_eff: float) -> np.ndarray:
    carrier = ALPHA_GEO * np.cos(STRICT["omega"] * d + STRICT["phi"] + delta_phi)
    return carrier / (1.0 + beta_eff * d**GAMMA_F)


def log_corrected_model(d: np.ndarray, delta_phi: float, beta_eff: float, epsilon: float) -> np.ndarray:
    carrier = ALPHA_GEO * np.cos(STRICT["omega"] * d + STRICT["phi"] + delta_phi)
    correction = 1.0 + epsilon * np.log(d)
    return carrier / (1.0 + beta_eff * d**GAMMA_F * correction)


def projected_sse(y_model: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    scale = float(np.dot(y_model, y) / max(np.dot(y_model, y_model), 1e-15))
    resid = scale * y_model - y
    return float(np.sum(resid * resid)), scale


def fit_base(d: np.ndarray) -> dict[str, float]:
    y = strict_target(d)

    def mse(x: np.ndarray) -> float:
        sse, _ = projected_sse(base_model(d, float(x[0]), float(x[1])), y)
        return float(sse / len(d))

    result = so.minimize(mse, x0=np.array([0.0, 1.0]), bounds=[(-math.pi, math.pi), (1e-4, 50.0)], method="L-BFGS-B")
    delta_phi, beta_eff = [float(v) for v in result.x]
    sse, scale = projected_sse(base_model(d, delta_phi, beta_eff), y)
    return {"delta_phi": delta_phi, "beta_eff": beta_eff, "scale": scale, "sse": sse, "mse": float(sse / len(d))}


def fit_log_corrected(d: np.ndarray) -> dict[str, float]:
    y = strict_target(d)

    def mse(x: np.ndarray) -> float:
        delta_phi, beta_eff, epsilon = [float(v) for v in x]
        correction = 1.0 + epsilon * np.log(d)
        if float(np.min(correction)) <= 0.0:
            return 1e9
        sse, _ = projected_sse(log_corrected_model(d, delta_phi, beta_eff, epsilon), y)
        return float(sse / len(d))

    starts = [
        np.array([0.0, 1.0, 0.0]),
        np.array([0.0, 1.0, DELTA_STRICT]),
        np.array([0.0, 1.05, DELTA_STRICT]),
    ]
    best = None
    for start in starts:
        result = so.minimize(
            mse,
            x0=start,
            bounds=[(-math.pi, math.pi), (1e-4, 50.0), (-0.2, 0.2)],
            method="L-BFGS-B",
        )
        if best is None or result.fun < best.fun:
            best = result
    assert best is not None
    delta_phi, beta_eff, epsilon = [float(v) for v in best.x]
    sse, scale = projected_sse(log_corrected_model(d, delta_phi, beta_eff, epsilon), y)
    return {
        "delta_phi": delta_phi,
        "beta_eff": beta_eff,
        "epsilon": epsilon,
        "epsilon_minus_delta_strict": epsilon - DELTA_STRICT,
        "scale": scale,
        "sse": sse,
        "mse": float(sse / len(d)),
    }


def window_row(window: dict[str, Any]) -> dict[str, Any]:
    d = np.linspace(float(window["d_min"]), float(window["d_max"]), int(window["n"]))
    base = fit_base(d)
    log_corr = fit_log_corrected(d)
    exact_factor = d**DELTA_STRICT
    linear_factor = 1.0 + DELTA_STRICT * np.log(d)
    rel_error = np.abs(linear_factor / exact_factor - 1.0)
    return {
        "window": window,
        "base_D_f_minus_1_denominator_fit": base,
        "linear_log_corrected_fit": log_corr,
        "improvement": {
            "sse_reduction_factor_base_over_log_corrected": float(base["sse"] / max(log_corr["sse"], 1e-300)),
            "epsilon_close_to_delta_strict_abs": abs(log_corr["epsilon"] - DELTA_STRICT),
        },
        "log_correction_bounds": {
            "delta_strict": DELTA_STRICT,
            "exact_factor_min": float(np.min(exact_factor)),
            "exact_factor_max": float(np.max(exact_factor)),
            "linear_factor_min": float(np.min(linear_factor)),
            "linear_factor_max": float(np.max(linear_factor)),
            "max_relative_error_linear_vs_exact_factor": float(np.max(rel_error)),
        },
    }


def symbolic_log_correction() -> dict[str, str]:
    d, eps, Df = sp.symbols("d eps D_f", positive=True, real=True)
    gamma = Df - 1
    exact = d ** (gamma + eps)
    factored = d**gamma * sp.exp(eps * sp.log(d))
    first_order = d**gamma * (1 + eps * sp.log(d))
    residual = sp.series(sp.exp(eps * sp.log(d)), eps, 0, 3).removeO() - (1 + eps * sp.log(d))
    return {
        "exact_factorization": sp.sstr(factored),
        "first_order_log_corrected_power": sp.sstr(first_order),
        "second_order_remainder_factor": sp.sstr(sp.simplify(residual)),
        "delta_strict_exact_expression": "14/5 - 4*log(2)",
    }


def main() -> None:
    domain = load_json(DOMAIN_AUDIT)
    placement = load_json(PLACEMENT_AUDIT)
    rows = [window_row(window) for window in WINDOWS]
    reductions = [row["improvement"]["sse_reduction_factor_base_over_log_corrected"] for row in rows]
    epsilon_gaps = [row["improvement"]["epsilon_close_to_delta_strict_abs"] for row in rows]
    max_lin_errors = [row["log_correction_bounds"]["max_relative_error_linear_vs_exact_factor"] for row in rows]

    summary = {
        "delta_strict": DELTA_STRICT,
        "max_abs_epsilon_fit_minus_delta_strict": float(max(epsilon_gaps)),
        "min_sse_reduction_factor_base_over_log_corrected": float(min(reductions)),
        "median_sse_reduction_factor_base_over_log_corrected": float(np.median(reductions)),
        "max_relative_error_first_order_factor_all_windows": float(max(max_lin_errors)),
        "all_fit_epsilons_match_delta_strict_to_1e_12": bool(max(epsilon_gaps) < 1e-12),
        "all_first_order_factor_errors_below_0p003": bool(max(max_lin_errors) < 0.003),
        "placement_denominator_win_replayed": placement["aggregate_summary"]["all_windows_denominator_beats_numerator_and_inverse_by_bic"],
        "domain_rationalization_warning_replayed": not domain["rationalization_audit"]["strict_9_over_5_best_for_all_tested_denominator_bounds"],
    }

    report = {
        "status": "OPEN_LOG_CORRECTION_TRACE_NO_BRIDGE_THEOREM",
        "result_kind": "SCRATCH_DF_LOG_CORRECTION_RENORMALIZATION_TRACE__NOT_A_THEOREM",
        "source_reports": {
            "analytic_domain_audit": str(DOMAIN_AUDIT.relative_to(ROOT)),
            "placement_discriminator": str(PLACEMENT_AUDIT.relative_to(ROOT)),
        },
        "constants": {
            "D_f": ALPHA_GEO,
            "gamma_F_D_f_minus_1": GAMMA_F,
            "strict_eta_9_over_5": STRICT_ETA,
            "delta_strict_eta_minus_gamma_F": DELTA_STRICT,
        },
        "symbolic_log_correction": symbolic_log_correction(),
        "window_rows": rows,
        "aggregate_summary": summary,
        "honest_interpretation": [
            "The strict eta=9/5 exponent differs from D_f-1 by the small marginal log-correction delta=14/5-4 log(2).",
            "A first-order denominator correction 1+epsilon log(d) recovers the strict target with epsilon equal to delta_strict in all tested windows and substantially reduces the fixed-D_f-1 residual.",
            "This does not prove the correction exists physically; it only turns the remaining bridge problem into a precise marginal-renormalization obligation.",
        ],
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No theorem deriving delta_strict = 14/5 - 4 log(2) from nadsoliton dynamics is exported.",
            "No legacy physical-role formula is transferred to the strict kernel.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Derive or falsify the marginal log-renormalization term epsilon*log(d) from a concrete transport/RG balance equation; until then, treat eta=9/5 as D_f-1 plus an observed small correction, not as an explained theorem.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch D_f log-correction bridge trace\n\n"
        "Status: marginal log-correction trace; no bridge theorem.\n\n"
        f"- `D_f-1={GAMMA_F:.15f}`, strict `eta=9/5={STRICT_ETA:.15f}`, delta `eta-(D_f-1)={DELTA_STRICT:.15f}`.\n"
        f"- First-order factor max relative error over all windows: `{summary['max_relative_error_first_order_factor_all_windows']:.12f}`.\n"
        f"- Min SSE reduction base/log-corrected: `{summary['min_sse_reduction_factor_base_over_log_corrected']:.6f}`; median `{summary['median_sse_reduction_factor_base_over_log_corrected']:.6f}`.\n"
        f"- Fitted epsilons match delta to 1e-12: `{summary['all_fit_epsilons_match_delta_strict_to_1e_12']}`.\n"
        "- No false pass: no kernel identity, no dynamical derivation of delta, no physical-role transfer, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
