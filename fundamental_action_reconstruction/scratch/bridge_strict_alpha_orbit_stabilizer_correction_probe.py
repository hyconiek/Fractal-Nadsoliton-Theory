#!/usr/bin/env python3
"""Scratch probe: orbit/stabilizer-style source for the eta finite-size correction.

The finite-size correction probe isolated a very sharp rational shadow:

    q ~= 96/95,      eta = ln(12*q)/ln(4) ~= 1.8000346966.

This packet asks whether 96 has a strict-side combinatorial interpretation that
is more than numerology.  Using only strict-side counts already in the current
route, the natural slot count

    B = nad12_support_size * sqrt(exp(alpha_geo_strict_derived_v1)) * 2
      = 12 * 4 * 2
      = 96

can be read as a 12-site support crossed with the Shannon half-microstate / two-
orientation phase budget.  A one-slot exclusion then gives B/(B-1)=96/95.  This
is still not a theorem: no repo packet proves that exactly one such slot is
forbidden by strict nadsoliton combinatorics or by QW-2191 selector resolution.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FINITE_SIZE = HERE / "bridge_strict_alpha_finite_size_correction_report.json"
ALPHA_PACKET = ROOT / "fundamental_action_reconstruction" / "generated" / "alpha_geo_strict_derived_v1.json"
OUT_JSON = HERE / "bridge_strict_alpha_orbit_stabilizer_correction_report.json"
OUT_MD = HERE / "bridge_strict_alpha_orbit_stabilizer_correction_report.md"

NAD12_SUPPORT_SIZE = 12
ALPHA_MICROSTATES = 16
ALPHA_SCALE = 4.0
ORIENTATION_FACTOR = 2
STRICT_TARGET_ETA = 9.0 / 5.0


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def eta_from_correction(correction: float) -> float:
    return math.log(NAD12_SUPPORT_SIZE * correction) / math.log(ALPHA_SCALE)


def one_slot_exclusion_row(label: str, slot_count: int) -> dict[str, Any]:
    correction = slot_count / (slot_count - 1)
    eta = eta_from_correction(correction)
    return {
        "label": label,
        "slot_count_B": slot_count,
        "excluded_slots": 1,
        "correction": correction,
        "correction_label": f"{slot_count}/{slot_count - 1}",
        "eta_candidate": eta,
        "eta_signed_residual_vs_9_5": eta - STRICT_TARGET_ETA,
        "eta_abs_residual_vs_9_5": abs(eta - STRICT_TARGET_ETA),
    }


def candidate_slot_counts() -> list[tuple[str, int]]:
    alpha_half_microstates = ALPHA_MICROSTATES // 2
    z12_pairs = NAD12_SUPPORT_SIZE // 2
    nonzero_z12_distances = NAD12_SUPPORT_SIZE - 1
    return [
        ("nad12_support", NAD12_SUPPORT_SIZE),
        ("nad12_support_times_orientation", NAD12_SUPPORT_SIZE * ORIENTATION_FACTOR),
        ("alpha_microstates_times_orientation", ALPHA_MICROSTATES * ORIENTATION_FACTOR),
        ("nad12_support_times_alpha_scale", NAD12_SUPPORT_SIZE * int(ALPHA_SCALE)),
        ("alpha_microstates_times_alpha_scale", ALPHA_MICROSTATES * int(ALPHA_SCALE)),
        ("nad12_pairs_times_alpha_microstates", z12_pairs * ALPHA_MICROSTATES),
        ("nad12_support_times_alpha_half_microstates", NAD12_SUPPORT_SIZE * alpha_half_microstates),
        ("nad12_support_times_alpha_scale_times_orientation", NAD12_SUPPORT_SIZE * int(ALPHA_SCALE) * ORIENTATION_FACTOR),
        ("nad12_support_times_nonzero_distances", NAD12_SUPPORT_SIZE * nonzero_z12_distances),
        ("nad12_support_squared", NAD12_SUPPORT_SIZE**2),
        ("nad12_support_times_alpha_microstates", NAD12_SUPPORT_SIZE * ALPHA_MICROSTATES),
    ]


def main() -> None:
    finite_size_report = load_json(FINITE_SIZE)
    alpha_packet = load_json(ALPHA_PACKET)
    exact_q = finite_size_report["exact_correction_required"]["q_target"]
    exact_one_slot_B = exact_q / (exact_q - 1.0)
    natural_B = NAD12_SUPPORT_SIZE * int(ALPHA_SCALE) * ORIENTATION_FACTOR
    natural_q = natural_B / (natural_B - 1.0)
    natural_eta = eta_from_correction(natural_q)
    rows = [one_slot_exclusion_row(label, count) for label, count in candidate_slot_counts()]
    rows.sort(key=lambda row: row["eta_abs_residual_vs_9_5"])
    best = rows[0]

    report = {
        "status": "OPEN_STRICT_ALPHA_ORBIT_STABILIZER_CORRECTION_CANDIDATE_NO_THEOREM",
        "result_kind": "SCRATCH_STRICT_ALPHA_ORBIT_STABILIZER_CORRECTION_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "finite_size_correction": str(FINITE_SIZE.relative_to(ROOT)),
            "alpha_geo_strict_derived_v1": str(ALPHA_PACKET.relative_to(ROOT)),
        },
        "strict_counts_used": {
            "alpha_geo_strict_derived_v1": alpha_packet["value"],
            "nad12_support_size": NAD12_SUPPORT_SIZE,
            "alpha_microstates": ALPHA_MICROSTATES,
            "alpha_scale_sqrt_exp_alpha": ALPHA_SCALE,
            "orientation_factor": ORIENTATION_FACTOR,
            "natural_slot_count_B": natural_B,
            "natural_slot_count_formula": "12 * sqrt(exp(4 ln 2)) * 2 = 12 * 4 * 2 = 96",
        },
        "exact_vs_natural_one_slot_exclusion": {
            "exact_q_target": exact_q,
            "exact_B_if_q_equals_B_over_B_minus_1": exact_one_slot_B,
            "natural_B": natural_B,
            "natural_q_B_over_B_minus_1": natural_q,
            "natural_q_signed_residual_vs_exact": natural_q - exact_q,
            "natural_eta": natural_eta,
            "natural_eta_signed_residual_vs_9_5": natural_eta - STRICT_TARGET_ETA,
            "natural_eta_abs_residual_vs_9_5": abs(natural_eta - STRICT_TARGET_ETA),
        },
        "slot_count_discriminator": {
            "one_slot_exclusion_rows_sorted": rows,
            "best_natural_row": best,
            "natural_96_is_best_in_candidate_library": best["slot_count_B"] == natural_B,
        },
        "candidate_interpretation": {
            "supported_by_this_probe": bool(
                best["slot_count_B"] == natural_B
                and best["correction_label"] == "96/95"
                and best["eta_abs_residual_vs_9_5"] < 5e-5
                and abs(exact_one_slot_B - natural_B) < 0.5
            ),
            "content": "The rational shadow 96/95 has a concrete strict-side slot-count reading: one excluded slot out of 12*4*2=96 phase/support slots.",
            "why_this_is_more_proof_like": "It turns the previous rational correction into an orbit/stabilizer-style finite-slot candidate and checks competing natural slot counts.",
            "why_this_is_not_enough": "No theorem identifies the 96 slots as the correct strict compression state space, and no theorem derives exactly one forbidden slot from QW-2191 or strict nadsoliton symmetry.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "No theorem proves that 12*4*2 is the intrinsic strict compression slot space.",
            "No theorem proves exactly one forbidden slot, so 96/95 remains a candidate correction.",
            "No theorem derives eta=9/5 exactly from this correction.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Try to derive the one-forbidden-slot rule from strict orbit/stabilizer combinatorics or from a selector premise; without that derivation, keep 96/95 as a sharp rational shadow, not an eta theorem.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha orbit/stabilizer correction probe\n\n"
        "Status: one-slot exclusion candidate for the strict-alpha eta correction; no theorem.\n\n"
        f"- Natural slot count: `B=12*4*2={natural_B}`, giving correction `{natural_B}/{natural_B-1}` and eta `{natural_eta:.12f}`.\n"
        f"- Exact one-slot B from q_target: `{exact_one_slot_B:.12f}`; natural residual `{natural_B-exact_one_slot_B:.3e}`.\n"
        f"- Candidate library best: `{best['label']}` with correction `{best['correction_label']}` and eta residual `{best['eta_signed_residual_vs_9_5']:.3e}`.\n"
        "- Honest read: 96/95 now has a concrete strict slot-count reading, but the one-forbidden-slot rule is not derived.\n"
        "- No false pass: no kernel identity, no legacy role transfer, no exact eta theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
