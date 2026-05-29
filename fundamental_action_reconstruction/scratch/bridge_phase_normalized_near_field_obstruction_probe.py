#!/usr/bin/env python3
"""Scratch theorem-prep probe for phase-normalized near-field bridge claims.

Previous robustness audit found that the raw full-kernel limit
K_legacy/K_strict at d->0 is not alpha_geo for the repo phases.  This probe
asks the next narrower question: can a phase-normalized near-field distance map
make the legacy damping picture and the strict eta=9/5 damping picture agree?

It records two honest facts:

1. Exact denominator matching forces x(d)=(beta/beta_tors)*d^eta.  For eta>1
   this kills the linear phase slope at the origin, while the strict phase has a
   nonzero linear slope.  Thus the exact-denominator map is not a full
   phase-normalized near-field bridge.
2. A local two-term map x(d)=a*d+c*d^eta can match the first ordinary phase
   slope and the first fractional damping coefficient, but it is only a local
   asymptotic witness, not a kernel identity or role-transfer theorem.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_phase_normalized_near_field_obstruction_report.json"
OUT_MD = HERE / "bridge_phase_normalized_near_field_obstruction_report.md"
ROBUSTNESS_AUDIT = HERE / "bridge_strict_kernel_robustness_opinion_audit_report.json"
RG_FLOW = HERE / "bridge_df_marginal_rg_flow_report.json"

ALPHA_GEO = 4.0 * math.log(2.0)
STRICT_ETA = 9.0 / 5.0
LEGACY = {"omega": math.pi / 4.0, "phi": math.pi / 6.0, "beta_tors": 0.01}
STRICT = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": STRICT_ETA}
SMALL_D = np.array([1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 1e-1], dtype=float)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def legacy_normalized(x: np.ndarray) -> np.ndarray:
    return (
        np.cos(LEGACY["omega"] * x + LEGACY["phi"])
        / math.cos(LEGACY["phi"])
        / (1.0 + LEGACY["beta_tors"] * x)
    )


def strict_normalized(d: np.ndarray) -> np.ndarray:
    return (
        np.cos(STRICT["omega"] * d + STRICT["phi"])
        / math.cos(STRICT["phi"])
        / (1.0 + STRICT["beta"] * d ** STRICT["eta"])
    )


def symbolic_coefficients() -> dict[str, Any]:
    omega_l, phi_l, beta_t, omega_s, phi_s, beta, eta = sp.symbols(
        "omega_l phi_l beta_t omega_s phi_s beta eta", positive=True, real=True
    )
    legacy_linear_in_x = -(omega_l * sp.tan(phi_l) + beta_t)
    strict_linear_in_d = -omega_s * sp.tan(phi_s)
    exact_denominator_x = beta / beta_t
    exact_denominator_phase_slope = sp.Integer(0)
    local_a = sp.simplify(strict_linear_in_d / legacy_linear_in_x)
    local_c = sp.simplify(-beta / legacy_linear_in_x)
    return {
        "legacy_normalized_linear_coefficient_in_x": sp.sstr(legacy_linear_in_x),
        "strict_normalized_linear_coefficient_in_d": sp.sstr(strict_linear_in_d),
        "exact_denominator_map_x_of_d": "(beta/beta_tors)*d**eta",
        "exact_denominator_prefactor_beta_over_beta_tors": sp.sstr(exact_denominator_x),
        "exact_denominator_map_linear_phase_slope_at_origin_for_eta_gt_1": sp.sstr(exact_denominator_phase_slope),
        "linear_phase_slope_residual_under_exact_denominator_map": sp.sstr(-strict_linear_in_d),
        "local_two_term_map": "x(d)=a*d+c*d**eta",
        "local_two_term_a_matching_linear_slope": sp.sstr(local_a),
        "local_two_term_c_matching_fractional_damping_coefficient": sp.sstr(local_c),
    }


def numeric_coefficients() -> dict[str, float]:
    legacy_linear = -(LEGACY["omega"] * math.tan(LEGACY["phi"]) + LEGACY["beta_tors"])
    strict_linear = -STRICT["omega"] * math.tan(STRICT["phi"])
    a = strict_linear / legacy_linear
    c = -STRICT["beta"] / legacy_linear
    return {
        "legacy_normalized_linear_coefficient_in_x": legacy_linear,
        "strict_normalized_linear_coefficient_in_d": strict_linear,
        "exact_denominator_prefactor_beta_over_beta_tors": STRICT["beta"] / LEGACY["beta_tors"],
        "linear_phase_slope_residual_under_exact_denominator_map": -strict_linear,
        "local_two_term_a": a,
        "local_two_term_c": c,
    }


def residual_rows() -> list[dict[str, float]]:
    coeffs = numeric_coefficients()
    d = SMALL_D
    target = strict_normalized(d)
    x_exact_denominator = coeffs["exact_denominator_prefactor_beta_over_beta_tors"] * d ** STRICT_ETA
    x_two_term = coeffs["local_two_term_a"] * d + coeffs["local_two_term_c"] * d ** STRICT_ETA
    exact_resid = legacy_normalized(x_exact_denominator) - target
    two_term_resid = legacy_normalized(x_two_term) - target
    rows = []
    for idx, value in enumerate(d):
        rows.append(
            {
                "d": float(value),
                "strict_normalized": float(target[idx]),
                "exact_denominator_x": float(x_exact_denominator[idx]),
                "exact_denominator_residual": float(exact_resid[idx]),
                "exact_denominator_abs_residual_over_d": float(abs(exact_resid[idx]) / value),
                "two_term_x": float(x_two_term[idx]),
                "two_term_residual": float(two_term_resid[idx]),
                "two_term_abs_residual_over_d_eta": float(abs(two_term_resid[idx]) / (value ** STRICT_ETA)),
            }
        )
    return rows


def main() -> None:
    robustness = load_json(ROBUSTNESS_AUDIT)
    rg_flow = load_json(RG_FLOW)
    symbolic = symbolic_coefficients()
    numeric = numeric_coefficients()
    rows = residual_rows()

    report = {
        "status": "OPEN_PHASE_NORMALIZED_NEAR_FIELD_PROBE_NO_BRIDGE_THEOREM",
        "result_kind": "SCRATCH_PHASE_NORMALIZED_NEAR_FIELD_OBSTRUCTION_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "strict_kernel_robustness_opinion_audit": str(ROBUSTNESS_AUDIT.relative_to(HERE.parents[1])),
            "df_marginal_rg_flow": str(RG_FLOW.relative_to(HERE.parents[1])),
        },
        "constants": {
            "alpha_geo": ALPHA_GEO,
            "strict_eta": STRICT_ETA,
            "legacy": LEGACY,
            "strict": STRICT,
        },
        "symbolic_coefficients": symbolic,
        "numeric_coefficients": numeric,
        "restricted_obstruction": {
            "claim": "Exact denominator matching x(d)=(beta/beta_tors)*d**eta cannot also match the nonzero strict linear phase slope at d=0 when eta>1.",
            "obstruction_residual_linear_coefficient": numeric["linear_phase_slope_residual_under_exact_denominator_map"],
            "is_zero": abs(numeric["linear_phase_slope_residual_under_exact_denominator_map"]) < 1e-14,
            "interpretation": "This blocks the strongest naive near-field story: exact fractal damping replacement plus a single shared distance variable is not a full normalized-kernel bridge.",
        },
        "local_two_term_witness": {
            "claim": "A two-term distance map x(d)=a*d+c*d**eta matches the first ordinary and first fractional coefficients of the normalized kernels.",
            "a": numeric["local_two_term_a"],
            "c": numeric["local_two_term_c"],
            "status": "ASYMPTOTIC_WITNESS_ONLY_NOT_GLOBAL_IDENTITY",
            "interpretation": "The bridge candidate would need a phase-carrying linear channel plus a fractal damping channel; denominator-only D_f is insufficient by itself.",
        },
        "sample_residual_rows": rows,
        "upstream_replay": {
            "robustness_support": robustness["support_score"],
            "rg_flow_status": rg_flow["status"],
            "rg_denominator_residual": rg_flow["symbolic_rg_flow"]["denominator_residual"],
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No physical-role transfer from legacy alpha_geo/beta_tors formulas is licensed.",
            "No theorem derives the two-term distance map from nadsoliton dynamics.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "If continuing this route, derive or falsify a two-channel near-field transport law: one analytic phase channel proportional to d and one fractal damping/RG channel proportional to d**eta.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch phase-normalized near-field obstruction probe\n\n"
        "Status: restricted obstruction plus local two-term witness; no bridge theorem.\n\n"
        f"- Exact denominator map `x=(beta/beta_tors)*d^eta` has prefactor `{numeric['exact_denominator_prefactor_beta_over_beta_tors']:.12f}`.\n"
        f"- It leaves nonzero strict linear phase residual `{numeric['linear_phase_slope_residual_under_exact_denominator_map']:.12f}`, so exact denominator matching is not a full near-field kernel bridge.\n"
        f"- Local two-term witness: `x(d)=a*d+c*d^eta` with `a={numeric['local_two_term_a']:.12f}`, `c={numeric['local_two_term_c']:.12f}`.\n"
        f"- Smallest sampled two-term residual at `d={rows[0]['d']}` is `{rows[0]['two_term_residual']:.12e}`; largest sampled exact-denominator residual over d is `{max(row['exact_denominator_abs_residual_over_d'] for row in rows):.12e}`.\n"
        "- No false pass: no kernel identity, no physical-role transfer, no derived two-channel law, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
