#!/usr/bin/env python3
"""Scratch probe: finite-size correction search for the strict-alpha eta shadow.

The preceding strict-alpha fractal-dimension probe found

    eta_alpha_z12 = ln(12)/ln(4) = 1.7924812503605783,

which is close to, but not equal to, the strict gate exponent eta=9/5.  This
packet makes the next proof obligation more concrete: if the strict Shannon
scale 4 and nad12 support count 12 are the right primitives, then exact eta=9/5
requires a small effective-count correction

    q_target = 4^(9/5) / 12 = 1.010477711006932.

We search simple rational finite-size corrections q=p/q_den with bounded
denominators and record the best candidates.  This is a discriminator only: a
small rational correction such as 96/95 nearly closes the residual, but no
strict-side theorem derives that correction.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ALPHA_DIMENSION = HERE / "bridge_strict_alpha_fractal_dimension_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_finite_size_correction_report.json"
OUT_MD = HERE / "bridge_strict_alpha_finite_size_correction_report.md"

NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0
STRICT_TARGET_ETA = 9.0 / 5.0
MAX_DENOMINATOR = 144
TOP_K = 12


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def eta_from_count_correction(correction: float) -> float:
    return math.log(NAD12_SUPPORT_SIZE * correction) / math.log(ALPHA_SCALE)


def rational_search(target: float, max_denominator: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for denominator in range(1, max_denominator + 1):
        lo = max(1, math.floor(target * denominator) - 3)
        hi = math.ceil(target * denominator) + 4
        for numerator in range(lo, hi + 1):
            correction = numerator / denominator
            if correction <= 0:
                continue
            eta_candidate = eta_from_count_correction(correction)
            rows.append(
                {
                    "numerator": numerator,
                    "denominator": denominator,
                    "correction": correction,
                    "correction_signed_residual": correction - target,
                    "correction_abs_residual": abs(correction - target),
                    "eta_candidate": eta_candidate,
                    "eta_signed_residual_vs_9_5": eta_candidate - STRICT_TARGET_ETA,
                    "eta_abs_residual_vs_9_5": abs(eta_candidate - STRICT_TARGET_ETA),
                }
            )
    rows.sort(key=lambda row: (row["eta_abs_residual_vs_9_5"], row["denominator"], row["numerator"]))
    deduped: list[dict[str, Any]] = []
    seen: set[tuple[int, int]] = set()
    for row in rows:
        fraction = Fraction(int(row["numerator"]), int(row["denominator"]))
        key = (fraction.numerator, fraction.denominator)
        if key in seen:
            continue
        seen.add(key)
        row["reduced_numerator"] = fraction.numerator
        row["reduced_denominator"] = fraction.denominator
        row["reduced_label"] = f"{fraction.numerator}/{fraction.denominator}"
        deduped.append(row)
        if len(deduped) >= TOP_K:
            break
    return deduped


def main() -> None:
    alpha_dimension_report = load_json(ALPHA_DIMENSION)
    eta_alpha = alpha_dimension_report["fractal_dimension_candidate"]["eta_alpha_z12"]
    baseline_residual = eta_alpha - STRICT_TARGET_ETA
    q_target = ALPHA_SCALE**STRICT_TARGET_ETA / NAD12_SUPPORT_SIZE
    exact_count = NAD12_SUPPORT_SIZE * q_target
    rational_rows = rational_search(q_target, MAX_DENOMINATOR)
    best = rational_rows[0]
    nearest_fraction = Fraction(q_target).limit_denominator(MAX_DENOMINATOR)
    nearest_eta = eta_from_count_correction(float(nearest_fraction))
    improvement_factor = abs(baseline_residual) / best["eta_abs_residual_vs_9_5"]

    report = {
        "status": "OPEN_STRICT_ALPHA_FINITE_SIZE_CORRECTION_SEARCH_NO_ETA_THEOREM",
        "result_kind": "SCRATCH_STRICT_ALPHA_FINITE_SIZE_CORRECTION_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "strict_alpha_fractal_dimension": str(ALPHA_DIMENSION.relative_to(ROOT)),
        },
        "problem_statement": {
            "baseline_formula": "eta_alpha_z12 = ln(12)/ln(4)",
            "baseline_eta_alpha_z12": eta_alpha,
            "strict_target_eta_9_5": STRICT_TARGET_ETA,
            "baseline_signed_residual": baseline_residual,
            "need": "Find a strict-side reason for a small effective-count correction q so ln(12*q)/ln(4)=9/5.",
        },
        "exact_correction_required": {
            "q_target": q_target,
            "effective_count_target_12_times_q": exact_count,
            "effective_count_residual_vs_12": exact_count - NAD12_SUPPORT_SIZE,
            "relative_count_correction": q_target - 1.0,
            "closed_form": "q_target = 4^(9/5)/12 = 2^(8/5)/3",
        },
        "bounded_rational_search": {
            "max_denominator": MAX_DENOMINATOR,
            "best_candidate": best,
            "top_candidates": rational_rows,
            "nearest_fraction_limit_denominator": {
                "label": f"{nearest_fraction.numerator}/{nearest_fraction.denominator}",
                "correction": float(nearest_fraction),
                "eta_candidate": nearest_eta,
                "eta_signed_residual_vs_9_5": nearest_eta - STRICT_TARGET_ETA,
                "eta_abs_residual_vs_9_5": abs(nearest_eta - STRICT_TARGET_ETA),
            },
            "baseline_to_best_eta_residual_improvement_factor": improvement_factor,
        },
        "candidate_interpretation": {
            "supported_by_this_probe": bool(
                best["reduced_label"] == "96/95"
                and best["eta_abs_residual_vs_9_5"] < 5e-5
                and improvement_factor > 100.0
            ),
            "content": "A tiny effective-count correction q≈96/95 turns the strict-alpha dimension shadow into eta≈1.8000346966, reducing the residual by more than two orders of magnitude.",
            "why_this_is_more_proof_like": "The open correction is no longer vague: exact eta fixes q_target, and a bounded rational search isolates a concrete small finite-size candidate.",
            "why_this_is_not_enough": "No strict-side combinatorial or Shannon theorem derives q=96/95 or q_target; this remains a numerically sharp proof obligation.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "No theorem derives q_target or 96/95 from strict nadsoliton combinatorics.",
            "No theorem derives eta=9/5 exactly from alpha_geo_strict_derived_v1 and nad12 support.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Try to derive q_target=2^(8/5)/3 or the rational shadow 96/95 from a strict finite-size orbit/stabilizer correction; otherwise keep eta=9/5 as gate-selected with a near strict-alpha dimension shadow.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha finite-size correction probe\n\n"
        "Status: finite-size correction search; no eta theorem.\n\n"
        f"- Exact correction required: `q_target={q_target:.15f}` so `ln(12*q)/ln(4)=9/5`; effective count `{exact_count:.12f}`.\n"
        f"- Best bounded rational (`den<=144`): `{best['reduced_label']}` with eta `{best['eta_candidate']:.12f}` and residual `{best['eta_signed_residual_vs_9_5']:.3e}`.\n"
        f"- Baseline residual `{baseline_residual:.3e}` improves by factor `{improvement_factor:.3f}` under the best rational correction.\n"
        "- Honest read: the correction target is now concrete, but no strict combinatorial theorem derives it.\n"
        "- No false pass: no kernel identity, no legacy role transfer, no exact eta theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
