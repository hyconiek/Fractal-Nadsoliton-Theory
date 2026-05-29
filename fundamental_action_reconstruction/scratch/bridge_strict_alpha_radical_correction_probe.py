#!/usr/bin/env python3
"""Scratch probe: exact radical correction route for strict alpha -> eta.

The Diophantine slot probe ruled out an exact eta=9/5 theorem from rational
finite-slot exclusions B/(B-k).  This packet records the complementary exact
route.  If the strict alpha scale is 4 and the support count is 12, then exact
eta=9/5 is equivalent to an algebraic correction

    q_rad = 4^(9/5) / 12 = 2^(8/5) / 3,
    243*q_rad^5 - 256 = 0,
    ln(12*q_rad) / ln(4) = 9/5.

This is an exact algebraic target, not a derivation.  It sharpens the next proof
obligation: derive the fifth-root correction from strict nadsoliton geometry, or
admit that eta=9/5 remains gate-selected with a precise radical shadow.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DIOPHANTINE = HERE / "bridge_strict_alpha_diophantine_slot_obstruction_report.json"
FINITE_SIZE = HERE / "bridge_strict_alpha_finite_size_correction_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_radical_correction_report.json"
OUT_MD = HERE / "bridge_strict_alpha_radical_correction_report.md"

NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0
STRICT_TARGET_ETA = 9.0 / 5.0
RATIONAL_SHADOW = 96 / 95


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def eta_from_correction(correction: float) -> float:
    return math.log(NAD12_SUPPORT_SIZE * correction) / math.log(ALPHA_SCALE)


def main() -> None:
    diophantine_report = load_json(DIOPHANTINE)
    finite_size_report = load_json(FINITE_SIZE)
    q_radical = ALPHA_SCALE**STRICT_TARGET_ETA / NAD12_SUPPORT_SIZE
    effective_count = NAD12_SUPPORT_SIZE * q_radical
    eta_exact = eta_from_correction(q_radical)
    polynomial_residual = 243.0 * q_radical**5 - 256.0
    rational_eta = eta_from_correction(RATIONAL_SHADOW)
    rational_eta_residual = rational_eta - STRICT_TARGET_ETA
    best_bounded = diophantine_report["bounded_integer_search"]["best_row"]

    report = {
        "status": "OPEN_STRICT_ALPHA_RADICAL_CORRECTION_TARGET_NO_DERIVATION_THEOREM",
        "result_kind": "SCRATCH_STRICT_ALPHA_RADICAL_CORRECTION_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "diophantine_slot_obstruction": str(DIOPHANTINE.relative_to(ROOT)),
            "finite_size_correction": str(FINITE_SIZE.relative_to(ROOT)),
        },
        "exact_radical_target": {
            "q_radical": q_radical,
            "closed_form": "q_radical = 4^(9/5)/12 = 2^(8/5)/3 = (256/243)^(1/5)",
            "minimal_polynomial_over_q_candidate": "243*q^5 - 256",
            "polynomial_residual_243_q5_minus_256": polynomial_residual,
            "effective_count_12_q": effective_count,
            "effective_count_closed_form": "12*q_radical = 4^(9/5) = 2^(18/5)",
            "eta_from_radical": eta_exact,
            "eta_residual_vs_9_5": eta_exact - STRICT_TARGET_ETA,
        },
        "comparison_to_rational_shadows": {
            "natural_96_over_95": {
                "q": RATIONAL_SHADOW,
                "eta": rational_eta,
                "eta_residual_vs_9_5": rational_eta_residual,
                "q_residual_vs_radical": RATIONAL_SHADOW - q_radical,
            },
            "best_bounded_diophantine_row": best_bounded,
            "diophantine_exact_rational_slot_possible": diophantine_report["exact_obstruction"]["exact_rational_slot_exclusion_possible"],
            "finite_size_best_candidate": finite_size_report["bounded_rational_search"]["best_candidate"],
        },
        "candidate_interpretation": {
            "supported_by_this_probe": bool(
                abs(polynomial_residual) < 1e-12
                and abs(eta_exact - STRICT_TARGET_ETA) < 1e-12
                and diophantine_report["exact_obstruction"]["exact_rational_slot_exclusion_possible"] is False
            ),
            "content": "Exact eta=9/5 is equivalent to a fifth-root strict-alpha correction q=(256/243)^(1/5), not to any rational slot exclusion.",
            "why_this_is_more_proof_like": "It gives an exact algebraic target with a polynomial identity and separates it from the rational approximant family.",
            "why_this_is_not_enough": "No strict-side geometric, Shannon, or selector theorem derives the fifth-root correction from nadsoliton structure.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "No theorem derives q_radical=(256/243)^(1/5) from strict nadsoliton geometry.",
            "No theorem derives eta=9/5 from alpha_geo_strict_derived_v1 without adopting the radical correction as an extra premise.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Try to derive the fifth-root correction from a strict five-branch compression/selector mechanism; if that fails, classify eta=9/5 as gate-selected with an exact radical shadow rather than a derived terminal formula.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha radical correction probe\n\n"
        "Status: exact radical target for eta=9/5; no derivation theorem.\n\n"
        f"- Exact correction: `q=(256/243)^(1/5)={q_radical:.15f}`, with polynomial residual `243*q^5-256={polynomial_residual:.3e}`.\n"
        f"- Effective count: `12*q={effective_count:.12f}=4^(9/5)`, giving eta residual `{eta_exact-STRICT_TARGET_ETA:.3e}`.\n"
        f"- Rational shadow replay: `96/95` eta residual `{rational_eta_residual:.3e}`; best bounded rational residual `{best_bounded['eta_signed_residual_vs_9_5']:.3e}`.\n"
        "- Honest read: exact eta needs a fifth-root correction; rational slots are approximants unless a new strict derivation appears.\n"
        "- No false pass: no kernel identity, no legacy role transfer, no radical-correction theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
