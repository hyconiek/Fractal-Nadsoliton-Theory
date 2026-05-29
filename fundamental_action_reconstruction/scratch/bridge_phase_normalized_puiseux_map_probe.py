#!/usr/bin/env python3
"""Scratch theorem-prep probe: phase-normalized Puiseux near-field map.

The previous near-field obstruction showed that exact denominator matching
x=(beta/beta_tors)*d**eta cannot also carry the strict linear phase.  This
probe asks the constructive follow-up: if a near-field map is allowed to have
both an analytic phase channel and a fractional damping channel,

    x(d) = a*d + c*d**eta + b*d**2 + ...,

can the normalized legacy kernel match the normalized strict kernel as a local
Puiseux expansion?  We match coefficients through ordinary order d**2 for
eta=9/5 and record the first remaining residual scale.  This is still only a
formal local witness, not a bridge theorem or a derived transport law.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_phase_normalized_puiseux_map_report.json"
OUT_MD = HERE / "bridge_phase_normalized_puiseux_map_report.md"
OBSTRUCTION = HERE / "bridge_phase_normalized_near_field_obstruction_report.json"
ROBUSTNESS = HERE / "bridge_strict_kernel_robustness_opinion_audit_report.json"

ALPHA_GEO = 4.0 * math.log(2.0)
ETA = 9.0 / 5.0
LEGACY = {"omega": math.pi / 4.0, "phi": math.pi / 6.0, "beta_tors": 0.01}
STRICT = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": ETA}
SMALL_D = np.array([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 5e-2], dtype=float)


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


def coefficient_data() -> dict[str, Any]:
    tan_l = math.tan(LEGACY["phi"])
    tan_s = math.tan(STRICT["phi"])
    omega_l = LEGACY["omega"]
    omega_s = STRICT["omega"]
    beta_t = LEGACY["beta_tors"]
    beta = STRICT["beta"]

    legacy_l1 = -(omega_l * tan_l + beta_t)
    legacy_l2 = beta_t**2 + beta_t * omega_l * tan_l - 0.5 * omega_l**2
    strict_s1 = -omega_s * tan_s
    strict_s_eta = -beta
    strict_s2 = -0.5 * omega_s**2

    a = strict_s1 / legacy_l1
    c = strict_s_eta / legacy_l1
    b = (strict_s2 - legacy_l2 * a**2) / legacy_l1

    matched = {
        "d^1_residual": legacy_l1 * a - strict_s1,
        "d^eta_residual": legacy_l1 * c - strict_s_eta,
        "d^2_residual": legacy_l1 * b + legacy_l2 * a**2 - strict_s2,
    }
    return {
        "legacy_coefficients": {
            "L1_x": legacy_l1,
            "L2_x2": legacy_l2,
            "formula_L1_x": "-(omega_l*tan(phi_l)+beta_tors)",
            "formula_L2_x2": "beta_tors**2 + beta_tors*omega_l*tan(phi_l) - omega_l**2/2",
        },
        "strict_coefficients": {
            "S1_d": strict_s1,
            "S_eta_d_eta": strict_s_eta,
            "S2_d2": strict_s2,
            "formula_S1_d": "-omega_s*tan(phi_s)",
            "formula_S_eta_d_eta": "-beta",
            "formula_S2_d2": "-omega_s**2/2",
        },
        "puiseux_map_coefficients": {
            "a": a,
            "c": c,
            "b": b,
            "map": "x(d)=a*d+c*d**eta+b*d**2+O(d**(eta+1))",
        },
        "matched_coefficient_residuals": matched,
        "first_unmatched_expected_order": ETA + 1.0,
    }


def residual_rows(coeffs: dict[str, Any]) -> list[dict[str, float]]:
    a = coeffs["puiseux_map_coefficients"]["a"]
    c = coeffs["puiseux_map_coefficients"]["c"]
    b = coeffs["puiseux_map_coefficients"]["b"]
    d = SMALL_D
    x_two = a * d + c * d**ETA
    x_three = x_two + b * d**2
    target = strict_normalized(d)
    two_resid = legacy_normalized(x_two) - target
    three_resid = legacy_normalized(x_three) - target
    rows = []
    for idx, value in enumerate(d):
        rows.append(
            {
                "d": float(value),
                "two_term_residual": float(two_resid[idx]),
                "three_term_puiseux_residual": float(three_resid[idx]),
                "abs_three_term_residual_over_d_eta_plus_1": float(abs(three_resid[idx]) / value ** (ETA + 1.0)),
                "improvement_abs_two_over_three": float(abs(two_resid[idx]) / max(abs(three_resid[idx]), 1e-300)),
                "x_three_term": float(x_three[idx]),
            }
        )
    return rows


def main() -> None:
    obstruction = load_json(OBSTRUCTION)
    robustness = load_json(ROBUSTNESS)
    coeffs = coefficient_data()
    rows = residual_rows(coeffs)
    max_match_abs = max(abs(value) for value in coeffs["matched_coefficient_residuals"].values())

    report = {
        "status": "OPEN_PHASE_NORMALIZED_PUISEUX_MAP_WITNESS_NO_BRIDGE_THEOREM",
        "result_kind": "SCRATCH_PHASE_NORMALIZED_PUISEUX_MAP_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "near_field_obstruction": str(OBSTRUCTION.relative_to(HERE.parents[1])),
            "strict_kernel_robustness": str(ROBUSTNESS.relative_to(HERE.parents[1])),
        },
        "constants": {
            "alpha_geo": ALPHA_GEO,
            "eta": ETA,
            "legacy": LEGACY,
            "strict": STRICT,
        },
        "coefficient_match": coeffs,
        "coefficient_match_passes_to_order_2": max_match_abs < 1e-14,
        "sample_residual_rows": rows,
        "upstream_replay": {
            "previous_obstruction_is_zero": obstruction["restricted_obstruction"]["is_zero"],
            "previous_two_term_status": obstruction["local_two_term_witness"]["status"],
            "robustness_support": robustness["support_score"],
        },
        "honest_interpretation": [
            "The exact-denominator one-channel map remains obstructed by phase slope.",
            "A formal local Puiseux map with analytic d and fractional d**eta channels can match normalized kernels through d**2.",
            "This is a local asymptotic construction only; no repo artifact derives this map from strict nadsoliton dynamics.",
        ],
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No global distance map on Z12 is derived.",
            "No physical-role transfer from legacy electroweak/gravity formulas is licensed.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Try to derive the Puiseux coefficients a,c,b from a strict-side two-channel transport/RG equation, or show that such a derivation is impossible under the current axioms.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch phase-normalized Puiseux near-field map probe\n\n"
        "Status: local coefficient witness through d^2; no bridge theorem.\n\n"
        f"- Coefficient residual max through `d`, `d^eta`, `d^2`: `{max_match_abs:.3e}`.\n"
        f"- Puiseux map: `x(d)=a*d+c*d^eta+b*d^2+O(d^(eta+1))` with `a={coeffs['puiseux_map_coefficients']['a']:.12f}`, `c={coeffs['puiseux_map_coefficients']['c']:.12f}`, `b={coeffs['puiseux_map_coefficients']['b']:.12f}`.\n"
        f"- First unmatched expected order: `d^{ETA + 1.0:.1f}`; stable sampled residual/d^(eta+1) at `d={rows[2]['d']}` is `{rows[2]['abs_three_term_residual_over_d_eta_plus_1']:.6e}`.\n"
        "- Honest read: one-channel exact denominator map is obstructed; two-channel local Puiseux matching is possible but not derived.\n"
        "- No false pass: no kernel identity, no global Z12 map, no physical-role transfer, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
