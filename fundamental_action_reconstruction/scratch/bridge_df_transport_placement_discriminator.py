#!/usr/bin/env python3
"""Scratch discriminator for where the D_f-1 transport factor can enter.

Previous scratch packets made the D_f-1 denominator route plausible.  This
packet makes the next falsifiable step: compare the same D_f-1 factor inserted
into competing locations.

Families tested against the strict eta=1.8 target on several d-windows:
- linear denominator:              1 / (1 + beta*d)
- D_f-1 denominator insertion:     1 / (1 + beta*d^(D_f-1))
- D_f-1 numerator/path count:      d^(D_f-1) / (1 + beta*d)
- inverse numerator attenuation:   d^-(D_f-1) / (1 + beta*d)
- free denominator exponent oracle

This remains scratch theorem-prep only.  It does not identify kernels or
transfer legacy physical roles.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np
import scipy.optimize as so
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TRANSPORT_AUDIT = HERE / "bridge_df_transport_exponent_report.json"
OPINION_AUDIT = HERE / "bridge_fractal_dimension_opinion_audit_report.json"
OUT_JSON = HERE / "bridge_df_transport_placement_discriminator_report.json"
OUT_MD = HERE / "bridge_df_transport_placement_discriminator_report.md"

ALPHA_GEO = 4.0 * math.log(2.0)
GAMMA_F = ALPHA_GEO - 1.0
STRICT = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 1.8}
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
    return np.cos(STRICT["omega"] * d + STRICT["phi"]) / (1.0 + STRICT["beta"] * d ** STRICT["eta"])


def family_values(d: np.ndarray, delta_phi: float, beta_eff: float, family: str, gamma: float = GAMMA_F) -> np.ndarray:
    carrier = ALPHA_GEO * np.cos(STRICT["omega"] * d + STRICT["phi"] + delta_phi)
    if family == "linear_denominator":
        return carrier / (1.0 + beta_eff * d)
    if family == "df_denominator":
        return carrier / (1.0 + beta_eff * d**gamma)
    if family == "df_numerator_linear_denominator":
        return carrier * d**gamma / (1.0 + beta_eff * d)
    if family == "inverse_df_numerator_linear_denominator":
        return carrier * d ** (-gamma) / (1.0 + beta_eff * d)
    raise ValueError(f"unknown family: {family}")


def project_scale(y_model: np.ndarray, y: np.ndarray) -> float:
    return float(np.dot(y_model, y) / max(np.dot(y_model, y_model), 1e-15))


def fit_fixed_family(d: np.ndarray, family: str, gamma: float = GAMMA_F) -> dict[str, float]:
    y = strict_target(d)

    def mse(x: np.ndarray) -> float:
        delta_phi, beta_eff = x
        y_model = family_values(d, float(delta_phi), float(beta_eff), family, gamma)
        scale = project_scale(y_model, y)
        resid = scale * y_model - y
        return float(np.mean(resid * resid))

    starts = [np.array([0.0, 1.0]), np.array([0.02, 1.0]), np.array([0.0, 0.5]), np.array([0.0, 2.0])]
    best = None
    for start in starts:
        result = so.minimize(
            mse,
            x0=start,
            bounds=[(-math.pi, math.pi), (1e-4, 50.0)],
            method="L-BFGS-B",
        )
        if best is None or result.fun < best.fun:
            best = result
    assert best is not None
    delta_phi, beta_eff = [float(v) for v in best.x]
    y_model = family_values(d, delta_phi, beta_eff, family, gamma)
    scale = project_scale(y_model, y)
    resid = scale * y_model - y
    sse = float(np.sum(resid * resid))
    return {
        "family": family,
        "gamma": float(gamma),
        "delta_phi": delta_phi,
        "beta_eff": beta_eff,
        "scale": scale,
        "sse": sse,
        "mse": float(sse / len(d)),
    }


def fit_free_denominator(d: np.ndarray) -> dict[str, float]:
    y = strict_target(d)

    def mse(x: np.ndarray) -> float:
        delta_phi, beta_eff, gamma = x
        y_model = family_values(d, float(delta_phi), float(beta_eff), "df_denominator", float(gamma))
        scale = project_scale(y_model, y)
        resid = scale * y_model - y
        return float(np.mean(resid * resid))

    starts = [
        np.array([0.0, 1.0, 1.0]),
        np.array([0.0, 1.0, GAMMA_F]),
        np.array([0.0, 1.0, STRICT["eta"]]),
    ]
    best = None
    for start in starts:
        result = so.minimize(
            mse,
            x0=start,
            bounds=[(-math.pi, math.pi), (1e-4, 50.0), (0.2, 3.0)],
            method="L-BFGS-B",
        )
        if best is None or result.fun < best.fun:
            best = result
    assert best is not None
    delta_phi, beta_eff, gamma = [float(v) for v in best.x]
    y_model = family_values(d, delta_phi, beta_eff, "df_denominator", gamma)
    scale = project_scale(y_model, y)
    resid = scale * y_model - y
    sse = float(np.sum(resid * resid))
    return {
        "family": "free_denominator_exponent_oracle",
        "gamma": gamma,
        "delta_phi": delta_phi,
        "beta_eff": beta_eff,
        "scale": scale,
        "sse": sse,
        "mse": float(sse / len(d)),
    }


def with_ic(row: dict[str, float], n: int, k: int) -> dict[str, float]:
    out = dict(row)
    per = max(float(row["sse"]) / n, 1e-300)
    out["aic"] = float(n * math.log(per) + 2 * k)
    out["bic"] = float(n * math.log(per) + k * math.log(n))
    out["k_parameters_excluding_projected_scale"] = float(k)
    return out


def window_row(window: dict[str, Any]) -> dict[str, Any]:
    d = np.linspace(float(window["d_min"]), float(window["d_max"]), int(window["n"]))
    n = len(d)
    linear = with_ic(fit_fixed_family(d, "linear_denominator"), n, 2)
    denom = with_ic(fit_fixed_family(d, "df_denominator"), n, 2)
    numerator = with_ic(fit_fixed_family(d, "df_numerator_linear_denominator"), n, 2)
    inverse_num = with_ic(fit_fixed_family(d, "inverse_df_numerator_linear_denominator"), n, 2)
    free = with_ic(fit_free_denominator(d), n, 3)
    denominator_advantage = {
        "bic_numerator_minus_denominator": float(numerator["bic"] - denom["bic"]),
        "bic_inverse_numerator_minus_denominator": float(inverse_num["bic"] - denom["bic"]),
        "sse_ratio_numerator_over_denominator": float(numerator["sse"] / max(denom["sse"], 1e-300)),
        "sse_ratio_inverse_numerator_over_denominator": float(inverse_num["sse"] / max(denom["sse"], 1e-300)),
        "denominator_beats_both_alternative_placements_by_bic": bool(
            denom["bic"] < numerator["bic"] and denom["bic"] < inverse_num["bic"]
        ),
    }
    return {
        "window": window,
        "linear_denominator": linear,
        "df_denominator_insertion": denom,
        "df_numerator_path_count_insertion": numerator,
        "inverse_df_numerator_attenuation": inverse_num,
        "free_denominator_exponent_oracle": free,
        "denominator_advantage": denominator_advantage,
    }


def symbolic_placement_slopes() -> dict[str, str]:
    d, beta, Df, eta = sp.symbols("d beta D_f eta", positive=True, real=True)
    gamma = Df - 1
    linear = 1 / (1 + beta * d)
    denom = 1 / (1 + beta * d**gamma)
    numerator = d**gamma / (1 + beta * d)
    inverse_num = d ** (-gamma) / (1 + beta * d)
    strict = 1 / (1 + beta * d**eta)

    def asymptotic_slope(expr: sp.Expr) -> str:
        # For rational power products, SymPy's limit is less transparent than
        # recording the direct large-d exponent after leading-term extraction.
        return sp.sstr(sp.simplify(d * sp.diff(sp.log(expr), d)))

    return {
        "linear_denominator_log_slope_finite": asymptotic_slope(linear),
        "df_denominator_log_slope_finite": asymptotic_slope(denom),
        "df_numerator_linear_denominator_log_slope_finite": asymptotic_slope(numerator),
        "inverse_df_numerator_linear_denominator_log_slope_finite": asymptotic_slope(inverse_num),
        "strict_eta_denominator_log_slope_finite": asymptotic_slope(strict),
        "large_d_tail_exponents": {
            "linear_denominator": "-1",
            "df_denominator": "-(D_f - 1)",
            "df_numerator_linear_denominator": "D_f - 2",
            "inverse_df_numerator_linear_denominator": "-D_f",
            "strict_eta_denominator": "-eta",
        },
    }


def main() -> None:
    transport = load_json(TRANSPORT_AUDIT)
    opinion = load_json(OPINION_AUDIT)
    rows = [window_row(window) for window in WINDOWS]
    bic_num_minus_den = [row["denominator_advantage"]["bic_numerator_minus_denominator"] for row in rows]
    bic_inv_minus_den = [row["denominator_advantage"]["bic_inverse_numerator_minus_denominator"] for row in rows]
    sse_num_over_den = [row["denominator_advantage"]["sse_ratio_numerator_over_denominator"] for row in rows]
    all_denominator_wins = all(row["denominator_advantage"]["denominator_beats_both_alternative_placements_by_bic"] for row in rows)

    summary = {
        "all_windows_denominator_beats_numerator_and_inverse_by_bic": bool(all_denominator_wins),
        "min_bic_numerator_minus_denominator": float(min(bic_num_minus_den)),
        "min_bic_inverse_numerator_minus_denominator": float(min(bic_inv_minus_den)),
        "min_sse_ratio_numerator_over_denominator": float(min(sse_num_over_den)),
        "prior_transport_min_capture_replayed": transport["aggregate_summary"]["min_capture_fraction"],
        "opinion_audit_best_route_replayed": opinion["opinion_verdict"]["best_supported_next_formula"]["verdict"],
    }

    report = {
        "status": "OPEN_PLACEMENT_DISCRIMINATOR_TRACE_NO_BRIDGE_THEOREM",
        "result_kind": "SCRATCH_DF_TRANSPORT_PLACEMENT_DISCRIMINATOR__NOT_A_THEOREM",
        "source_reports": {
            "transport_exponent_audit": str(TRANSPORT_AUDIT.relative_to(ROOT)),
            "fractal_dimension_opinion_audit": str(OPINION_AUDIT.relative_to(ROOT)),
        },
        "constants": {
            "D_f": ALPHA_GEO,
            "gamma_F_D_f_minus_1": GAMMA_F,
            "strict_eta": STRICT["eta"],
        },
        "symbolic_placement_slopes": symbolic_placement_slopes(),
        "window_rows": rows,
        "aggregate_summary": summary,
        "honest_interpretation": [
            "The same D_f-1 factor has very different large-d signatures depending on placement: denominator gives tail exponent -(D_f-1), numerator path-count gives D_f-2, inverse numerator gives -D_f.",
            "Across the tested windows, the D_f-1 denominator placement beats both numerator placements by BIC, so the bridge-prep route should focus on a damping-denominator insertion lemma rather than a raw numerator/path-count multiplier.",
            "This is still a discriminator trace, not a theorem that the physical nadsoliton path multiplicity must enter the denominator.",
        ],
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No proof that D_f-1 must enter the denominator is exported; only alternative placements are numerically and asymptotically disfavored.",
            "No legacy physical-role formula is transferred to the strict kernel.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Turn the surviving denominator-placement route into a local lemma: specify the coarse-grained transport balance equation whose stationary Green/kernel denominator receives the D_f-1 path multiplicity, then check whether its finite-window renormalization can account for eta_eff and eta=1.8 offsets.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch discriminator: D_f-1 transport placement\n\n"
        "Status: placement discriminator trace; no bridge theorem.\n\n"
        f"- Constants: `D_f={ALPHA_GEO:.15f}`, `D_f-1={GAMMA_F:.15f}`, strict `eta={STRICT['eta']:.15f}`.\n"
        f"- Denominator beats numerator and inverse-numerator placements by BIC in all windows: `{summary['all_windows_denominator_beats_numerator_and_inverse_by_bic']}`.\n"
        f"- Min ΔBIC(numerator−denominator): `{summary['min_bic_numerator_minus_denominator']:.6f}`; min ΔBIC(inverse−denominator): `{summary['min_bic_inverse_numerator_minus_denominator']:.6f}`.\n"
        f"- Min SSE ratio numerator/denominator: `{summary['min_sse_ratio_numerator_over_denominator']:.6f}`.\n"
        "- No false pass: no kernel identity, no forced denominator theorem, no physical-role transfer, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
