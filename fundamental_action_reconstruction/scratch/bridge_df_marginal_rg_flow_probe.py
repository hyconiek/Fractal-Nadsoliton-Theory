#!/usr/bin/env python3
"""Scratch theorem-prep probe: marginal RG flow from D_f-1 to strict eta=9/5.

Previous packet isolated the gap

    delta = 9/5 - (D_f - 1) = 14/5 - 4 log(2)

as a small log correction.  This packet records the exact ODE form behind that
correction:

    ell = log(d),    dB/dell = delta * B,
    denominator = 1 + B(ell) * exp((D_f-1)*ell).

The solution B(ell)=B0 exp(delta ell) gives exactly

    1 + B0*d^(D_f-1+delta) = 1 + B0*d^(9/5).

This is still not a bridge theorem: the ODE is a candidate marginal RG law, not
a derived nadsoliton dynamic.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LOG_CORRECTION = HERE / "bridge_df_log_correction_report.json"
DOMAIN_AUDIT = HERE / "bridge_kernel_analytic_domain_audit_report.json"
OUT_JSON = HERE / "bridge_df_marginal_rg_flow_report.json"
OUT_MD = HERE / "bridge_df_marginal_rg_flow_report.md"

ALPHA_GEO = 4.0 * math.log(2.0)
GAMMA_F = ALPHA_GEO - 1.0
STRICT_ETA = float(Fraction(9, 5))
DELTA = STRICT_ETA - GAMMA_F
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


def exact_rg_model(d: np.ndarray) -> np.ndarray:
    # Includes alpha_geo amplitude; projected scale should be exactly 1/alpha_geo.
    return ALPHA_GEO * np.cos(STRICT["omega"] * d + STRICT["phi"]) / (
        1.0 + d**GAMMA_F * np.exp(DELTA * np.log(d))
    )


def projected_residual(y_model: np.ndarray, y: np.ndarray) -> dict[str, float]:
    scale = float(np.dot(y_model, y) / max(np.dot(y_model, y_model), 1e-15))
    resid = scale * y_model - y
    return {
        "scale": scale,
        "sse": float(np.sum(resid * resid)),
        "max_abs_residual": float(np.max(np.abs(resid))),
    }


def window_row(window: dict[str, Any]) -> dict[str, Any]:
    d = np.linspace(float(window["d_min"]), float(window["d_max"]), int(window["n"]))
    ell = np.log(d)
    y = strict_target(d)
    residual = projected_residual(exact_rg_model(d), y)
    x = DELTA * ell
    exp_factor = np.exp(x)
    first_order = 1.0 + x
    second_order = 1.0 + x + 0.5 * x * x
    return {
        "window": window,
        "ell_bounds": {"ell_min": float(np.min(ell)), "ell_max": float(np.max(ell))},
        "delta_ell_bounds": {"min": float(np.min(x)), "max": float(np.max(x)), "max_abs": float(np.max(np.abs(x)))},
        "exact_rg_residual_vs_strict": residual,
        "expected_projected_scale_1_over_alpha_geo": 1.0 / ALPHA_GEO,
        "first_order_factor_error": {
            "max_abs_exp_minus_linear": float(np.max(np.abs(exp_factor - first_order))),
            "max_rel_linear_vs_exp": float(np.max(np.abs(first_order / exp_factor - 1.0))),
        },
        "second_order_factor_error": {
            "max_abs_exp_minus_second_order": float(np.max(np.abs(exp_factor - second_order))),
            "max_rel_second_order_vs_exp": float(np.max(np.abs(second_order / exp_factor - 1.0))),
        },
    }


def symbolic_rg_flow() -> dict[str, str]:
    ell, beta0 = sp.symbols("ell beta0", real=True)
    gamma_f = 4 * sp.log(2) - 1
    delta = sp.Rational(14, 5) - 4 * sp.log(2)
    beta_running = beta0 * sp.exp(delta * ell)
    denominator = 1 + beta_running * sp.exp(gamma_f * ell)
    strict_denominator = 1 + beta0 * sp.exp(sp.Rational(9, 5) * ell)
    residual = sp.simplify(denominator - strict_denominator)
    ode_residual = sp.simplify(sp.diff(beta_running, ell) - delta * beta_running)
    first_order = sp.series(sp.exp(delta * ell), delta, 0, 3).removeO()
    return {
        "gamma_F": sp.sstr(gamma_f),
        "delta": sp.sstr(delta),
        "gamma_F_plus_delta": sp.sstr(sp.simplify(gamma_f + delta)),
        "running_beta_solution": sp.sstr(beta_running),
        "rg_ode_residual": sp.sstr(ode_residual),
        "denominator_after_flow": sp.sstr(denominator),
        "strict_denominator": sp.sstr(strict_denominator),
        "denominator_residual": sp.sstr(residual),
        "first_order_running_beta_factor": sp.sstr(first_order),
    }


def main() -> None:
    log_correction = load_json(LOG_CORRECTION)
    domain = load_json(DOMAIN_AUDIT)
    rows = [window_row(window) for window in WINDOWS]
    max_exact_sse = max(row["exact_rg_residual_vs_strict"]["sse"] for row in rows)
    max_exact_abs = max(row["exact_rg_residual_vs_strict"]["max_abs_residual"] for row in rows)
    max_delta_ell = max(row["delta_ell_bounds"]["max_abs"] for row in rows)
    max_first_rel = max(row["first_order_factor_error"]["max_rel_linear_vs_exp"] for row in rows)
    max_second_rel = max(row["second_order_factor_error"]["max_rel_second_order_vs_exp"] for row in rows)

    summary = {
        "delta": DELTA,
        "max_abs_delta_ell_all_windows": float(max_delta_ell),
        "max_exact_rg_sse_vs_strict": float(max_exact_sse),
        "max_exact_rg_abs_residual_vs_strict": float(max_exact_abs),
        "max_first_order_relative_factor_error": float(max_first_rel),
        "max_second_order_relative_factor_error": float(max_second_rel),
        "exact_rg_matches_strict_to_numeric_precision": bool(max_exact_abs < 1e-14),
        "second_order_improves_first_order_factor_error": bool(max_second_rel < max_first_rel),
        "previous_first_order_epsilon_match_replayed": log_correction["aggregate_summary"]["all_fit_epsilons_match_delta_strict_to_1e_12"],
        "domain_warning_replayed": not domain["rationalization_audit"]["strict_9_over_5_best_for_all_tested_denominator_bounds"],
    }

    report = {
        "status": "OPEN_MARGINAL_RG_FLOW_TRACE_NO_BRIDGE_THEOREM",
        "result_kind": "SCRATCH_DF_MARGINAL_RG_FLOW_TRACE__NOT_A_THEOREM",
        "source_reports": {
            "log_correction_trace": str(LOG_CORRECTION.relative_to(ROOT)),
            "analytic_domain_audit": str(DOMAIN_AUDIT.relative_to(ROOT)),
        },
        "constants": {
            "D_f": ALPHA_GEO,
            "gamma_F_D_f_minus_1": GAMMA_F,
            "strict_eta_9_over_5": STRICT_ETA,
            "delta_14_over_5_minus_4log2": DELTA,
        },
        "symbolic_rg_flow": symbolic_rg_flow(),
        "window_rows": rows,
        "aggregate_summary": summary,
        "honest_interpretation": [
            "An exact marginal flow dB/dlog(d)=delta*B converts the D_f-1 denominator into the strict 9/5 denominator with zero symbolic residual.",
            "The previous first-order log correction is the small-delta Euler/Taylor shadow of this flow; second-order terms reduce the factor error further.",
            "This is a sharper bridge target, not a bridge theorem: the missing ingredient is a derivation of this marginal RG flow and of delta=14/5-4log(2) from strict-side/nadsoliton dynamics.",
        ],
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No derivation of the marginal RG law dB/dlog(d)=delta*B is exported.",
            "No derivation of delta=14/5-4log(2) from strict-side dynamics is exported.",
            "No legacy physical-role formula is transferred to the strict kernel.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Search for or construct a strict-side transport/RG balance whose linearized beta function is dB/dlog(d)=(14/5-4log(2))*B; absent that, keep this as an exact reparameterization target only.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch D_f marginal RG-flow bridge trace\n\n"
        "Status: exact RG-flow target recorded; no bridge theorem.\n\n"
        f"- `D_f-1={GAMMA_F:.15f}`, strict `eta=9/5={STRICT_ETA:.15f}`, delta `{DELTA:.15f}`.\n"
        f"- Symbolic flow gives denominator residual `{report['symbolic_rg_flow']['denominator_residual']}` and ODE residual `{report['symbolic_rg_flow']['rg_ode_residual']}`.\n"
        f"- Max exact RG residual vs strict over windows: `{summary['max_exact_rg_abs_residual_vs_strict']:.3e}`.\n"
        f"- Max first/second-order factor relative errors: `{summary['max_first_order_relative_factor_error']:.12f}` / `{summary['max_second_order_relative_factor_error']:.12f}`.\n"
        "- No false pass: no kernel identity, no derived RG law, no derived delta, no physical-role transfer, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
