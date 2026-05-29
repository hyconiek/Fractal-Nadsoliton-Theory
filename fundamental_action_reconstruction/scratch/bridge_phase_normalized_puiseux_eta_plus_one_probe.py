#!/usr/bin/env python3
"""Scratch theorem-prep probe: extend Puiseux matching through d^(eta+1).

The previous Puiseux packet matched the normalized legacy and strict kernels at
orders d, d^eta, and d^2 using

    x(d) = a*d + c*d^eta + b*d^2 + ... .

This packet takes one more honest local step.  It adds q*d^(eta+1) and matches
the mixed phase-damping order d^(eta+1).  For eta=9/5, the first unmatched order
then moves to d^3.  This is still only a formal local asymptotic witness: no
strict-side transport law derives the map and no global Z12 bridge is exported.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_phase_normalized_puiseux_eta_plus_one_report.json"
OUT_MD = HERE / "bridge_phase_normalized_puiseux_eta_plus_one_report.md"
PUISEUX_PREV = HERE / "bridge_phase_normalized_puiseux_map_report.json"
OBSTRUCTION = HERE / "bridge_phase_normalized_near_field_obstruction_report.json"

ETA = 9.0 / 5.0
ALPHA_GEO = 4.0 * math.log(2.0)
LEGACY = {"omega": math.pi / 4.0, "phi": math.pi / 6.0, "beta_tors": 0.01}
STRICT = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": ETA}
SMALL_D = np.array([1e-5, 1e-4, 1e-3, 1e-2, 2e-2, 5e-2], dtype=float)


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
        / (1.0 + STRICT["beta"] * d ** ETA)
    )


def coefficients() -> dict[str, Any]:
    tan_l = math.tan(LEGACY["phi"])
    tan_s = math.tan(STRICT["phi"])
    wl = LEGACY["omega"]
    ws = STRICT["omega"]
    bt = LEGACY["beta_tors"]
    beta = STRICT["beta"]

    # Normalized legacy L(x)=1+L1*x+L2*x^2+O(x^3).
    l1 = -(wl * tan_l + bt)
    l2 = bt**2 + bt * wl * tan_l - 0.5 * wl**2
    # Normalized strict S(d)=1+S1*d+S2*d^2+Seta*d^eta+Se1*d^(eta+1)+...
    s1 = -ws * tan_s
    s2 = -0.5 * ws**2
    s_eta = -beta
    s_eta_plus_1 = -s1 * beta

    a = s1 / l1
    c = s_eta / l1
    b = (s2 - l2 * a**2) / l1
    q = (s_eta_plus_1 - 2.0 * l2 * a * c) / l1

    residuals = {
        "d^1_residual": l1 * a - s1,
        "d^eta_residual": l1 * c - s_eta,
        "d^2_residual": l1 * b + l2 * a**2 - s2,
        "d^(eta+1)_residual": l1 * q + 2.0 * l2 * a * c - s_eta_plus_1,
    }
    return {
        "legacy_expansion": {
            "L1": l1,
            "L2": l2,
            "L_of_x_to_needed_order": "1 + L1*x + L2*x**2 + O(x**3)",
        },
        "strict_expansion": {
            "S1": s1,
            "S2": s2,
            "S_eta": s_eta,
            "S_eta_plus_1": s_eta_plus_1,
            "S_of_d_to_needed_order": "1 + S1*d + S_eta*d**eta + S2*d**2 + S_eta_plus_1*d**(eta+1) + O(d**3)",
        },
        "puiseux_coefficients": {
            "a": a,
            "c": c,
            "b": b,
            "q": q,
            "map": "x(d)=a*d+c*d**eta+b*d**2+q*d**(eta+1)+O(d**3)",
        },
        "matched_coefficient_residuals": residuals,
        "first_unmatched_expected_order": 3.0,
    }


def residual_rows(coeffs: dict[str, Any]) -> list[dict[str, float]]:
    pc = coeffs["puiseux_coefficients"]
    d = SMALL_D
    x3 = pc["a"] * d + pc["c"] * d**ETA + pc["b"] * d**2
    x4 = x3 + pc["q"] * d ** (ETA + 1.0)
    target = strict_normalized(d)
    resid3 = legacy_normalized(x3) - target
    resid4 = legacy_normalized(x4) - target
    rows = []
    for idx, value in enumerate(d):
        rows.append(
            {
                "d": float(value),
                "three_term_residual": float(resid3[idx]),
                "eta_plus_one_residual": float(resid4[idx]),
                "abs_eta_plus_one_residual_over_d3": float(abs(resid4[idx]) / value**3),
                "improvement_abs_three_over_eta_plus_one": float(abs(resid3[idx]) / max(abs(resid4[idx]), 1e-300)),
                "x_eta_plus_one": float(x4[idx]),
            }
        )
    return rows


def main() -> None:
    previous = load_json(PUISEUX_PREV)
    obstruction = load_json(OBSTRUCTION)
    coeffs = coefficients()
    rows = residual_rows(coeffs)
    max_match_abs = max(abs(value) for value in coeffs["matched_coefficient_residuals"].values())

    report = {
        "status": "OPEN_PHASE_NORMALIZED_PUISEUX_ETA_PLUS_ONE_WITNESS_NO_BRIDGE_THEOREM",
        "result_kind": "SCRATCH_PHASE_NORMALIZED_PUISEUX_ETA_PLUS_ONE_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "previous_puiseux_map": str(PUISEUX_PREV.relative_to(HERE.parents[1])),
            "near_field_obstruction": str(OBSTRUCTION.relative_to(HERE.parents[1])),
        },
        "constants": {"alpha_geo": ALPHA_GEO, "eta": ETA, "legacy": LEGACY, "strict": STRICT},
        "coefficient_match": coeffs,
        "coefficient_match_passes_through_eta_plus_one": max_match_abs < 1e-14,
        "sample_residual_rows": rows,
        "upstream_replay": {
            "previous_order_2_match": previous["coefficient_match_passes_to_order_2"],
            "previous_obstruction_is_zero": obstruction["restricted_obstruction"]["is_zero"],
            "previous_first_unmatched_expected_order": previous["coefficient_match"]["first_unmatched_expected_order"],
        },
        "honest_interpretation": [
            "The local two-channel Puiseux map can be extended one more mixed order through d^(eta+1).",
            "For eta=9/5 this pushes the first unmatched formal order to d^3.",
            "This strengthens the local asymptotic bridge candidate but still does not derive a transport law or global Z12 map.",
        ],
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No global distance map on Z12 is derived.",
            "No strict-side derivation of the Puiseux coefficients is exported.",
            "No legacy physical-role transfer is licensed.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Derive the coefficient recurrence for all Puiseux orders, or identify an obstruction at finite order before trying any global Z12 extension.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch phase-normalized Puiseux eta+1 probe\n\n"
        "Status: local coefficient witness through d^(eta+1); no bridge theorem.\n\n"
        f"- Coefficient residual max through `d`, `d^eta`, `d^2`, `d^(eta+1)`: `{max_match_abs:.3e}`.\n"
        f"- Extended map coefficient `q={coeffs['puiseux_coefficients']['q']:.12f}` in `x=a*d+c*d^eta+b*d^2+q*d^(eta+1)+O(d^3)`.\n"
        f"- First unmatched expected order is now `d^{coeffs['first_unmatched_expected_order']:.1f}`; stable residual/d^3 at `d={rows[2]['d']}` is `{rows[2]['abs_eta_plus_one_residual_over_d3']:.6e}`.\n"
        "- Honest read: stronger local two-channel asymptotic witness, still no derived transport law or global Z12 bridge.\n"
        "- No false pass: no kernel identity, no global Z12 map, no strict-side coefficient derivation, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
