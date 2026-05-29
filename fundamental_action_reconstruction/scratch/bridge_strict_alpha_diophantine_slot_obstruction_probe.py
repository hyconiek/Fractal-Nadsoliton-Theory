#!/usr/bin/env python3
"""Scratch probe: Diophantine obstruction for exact slot-exclusion eta.

The orbit/stabilizer correction probe gave the sharp natural shadow
96/95.  This packet checks the exactness question: can a rational finite-slot
exclusion correction

    q = B / (B-k)       with integers B>k>0

produce the exact strict eta=9/5 when the alpha scale is 4 and the support count
is 12?  The required correction is

    q_target = 4^(9/5)/12 = 2^(8/5)/3,

which is irrational.  Therefore no rational finite-slot exclusion can be an
exact eta theorem.  We still scan bounded integer slots to record the best
approximants and to keep the natural 96/95 shadow in context.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ORBIT_STABILIZER = HERE / "bridge_strict_alpha_orbit_stabilizer_correction_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_diophantine_slot_obstruction_report.json"
OUT_MD = HERE / "bridge_strict_alpha_diophantine_slot_obstruction_report.md"

NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0
STRICT_TARGET_ETA = 9.0 / 5.0
MAX_SLOT_COUNT = 4096
MAX_EXCLUDED_SLOTS = 16
TOP_K = 16


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def eta_from_correction(correction: float) -> float:
    return math.log(NAD12_SUPPORT_SIZE * correction) / math.log(ALPHA_SCALE)


def slot_row(slot_count: int, excluded_slots: int, q_target: float) -> dict[str, Any]:
    correction = slot_count / (slot_count - excluded_slots)
    eta = eta_from_correction(correction)
    return {
        "slot_count_B": slot_count,
        "excluded_slots_k": excluded_slots,
        "correction": correction,
        "correction_label": f"{slot_count}/{slot_count - excluded_slots}",
        "correction_signed_residual_vs_q_target": correction - q_target,
        "correction_abs_residual_vs_q_target": abs(correction - q_target),
        "eta_candidate": eta,
        "eta_signed_residual_vs_9_5": eta - STRICT_TARGET_ETA,
        "eta_abs_residual_vs_9_5": abs(eta - STRICT_TARGET_ETA),
        "complexity_B_plus_k": slot_count + excluded_slots,
    }


def bounded_slot_search(q_target: float) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for excluded_slots in range(1, MAX_EXCLUDED_SLOTS + 1):
        center = q_target * excluded_slots / (q_target - 1.0)
        lo = max(excluded_slots + 1, int(math.floor(center)) - 5)
        hi = min(MAX_SLOT_COUNT, int(math.ceil(center)) + 5)
        for slot_count in range(lo, hi + 1):
            rows.append(slot_row(slot_count, excluded_slots, q_target))
    rows.sort(key=lambda row: (row["eta_abs_residual_vs_9_5"], row["complexity_B_plus_k"]))
    return rows[:TOP_K]


def main() -> None:
    orbit_report = load_json(ORBIT_STABILIZER)
    q_target = ALPHA_SCALE**STRICT_TARGET_ETA / NAD12_SUPPORT_SIZE
    exact_ratio_B_over_k = q_target / (q_target - 1.0)
    natural_B = orbit_report["strict_counts_used"]["natural_slot_count_B"]
    natural_row = slot_row(int(natural_B), 1, q_target)
    top_rows = bounded_slot_search(q_target)
    best = top_rows[0]

    report_exact_possible = False

    report = {
        "status": "OPEN_STRICT_ALPHA_DIOPHANTINE_SLOT_OBSTRUCTION_NO_ETA_THEOREM",
        "result_kind": "SCRATCH_STRICT_ALPHA_DIOPHANTINE_SLOT_OBSTRUCTION_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "orbit_stabilizer_correction": str(ORBIT_STABILIZER.relative_to(ROOT)),
        },
        "exact_obstruction": {
            "q_target": q_target,
            "closed_form": "q_target = 4^(9/5)/12 = 2^(8/5)/3",
            "irrationality_reason": "2^(1/5) is irrational, hence 2^(8/5)/3 is irrational.",
            "rational_slot_exclusion_form": "B/(B-k) is rational for integer B,k.",
            "exact_rational_slot_exclusion_possible": report_exact_possible,
            "exact_B_over_k_if_formally_solving_B_over_B_minus_k": exact_ratio_B_over_k,
        },
        "natural_96_shadow_replay": {
            "natural_slot_count_B": natural_B,
            "natural_one_slot_row": natural_row,
            "orbit_probe_supported": orbit_report["candidate_interpretation"]["supported_by_this_probe"],
        },
        "bounded_integer_search": {
            "max_slot_count_B": MAX_SLOT_COUNT,
            "max_excluded_slots_k": MAX_EXCLUDED_SLOTS,
            "best_row": best,
            "top_rows": top_rows,
            "best_is_natural_96_one_slot": bool(best["slot_count_B"] == natural_B and best["excluded_slots_k"] == 1),
            "natural_96_rank_message": "96/95 is the simplest strict slot-count shadow, not the best bounded rational approximant once larger B,k are allowed.",
        },
        "candidate_interpretation": {
            "supported_by_this_probe": bool(
                report_exact_possible is False
                and natural_row["eta_abs_residual_vs_9_5"] < 5e-5
                and best["eta_abs_residual_vs_9_5"] < 1e-6
            ),
            "content": "The exact rational slot-exclusion theorem is obstructed; 96/95 remains a simple strict slot-count shadow, while larger rational slots approximate eta even better.",
            "why_this_is_more_proof_like": "It separates an exact impossibility statement from bounded approximation evidence, preventing a false exact eta derivation from rational slot counts.",
            "why_this_is_not_enough": "It does not derive the exact irrational correction q_target, nor does it prove which approximate rational shadow should be physically selected.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "No exact eta=9/5 theorem follows from any rational finite-slot exclusion B/(B-k).",
            "No theorem selects 96/95 over other rational approximants.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "If exact eta=9/5 is required, search for a strict derivation of the irrational correction q_target=2^(8/5)/3; if only rational slot shadows are allowed, explicitly downgrade them to approximants.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Diophantine slot obstruction probe\n\n"
        "Status: exact rational slot-exclusion eta theorem is obstructed; rational shadows remain approximants.\n\n"
        f"- Exact target: `q_target={q_target:.15f}=2^(8/5)/3`, irrational, so no integer `B/(B-k)` can equal it exactly.\n"
        f"- Natural 96 shadow: `{natural_row['correction_label']}` gives eta residual `{natural_row['eta_signed_residual_vs_9_5']:.3e}`.\n"
        f"- Best bounded search (`B<={MAX_SLOT_COUNT}`, `k<={MAX_EXCLUDED_SLOTS}`): `B={best['slot_count_B']}`, `k={best['excluded_slots_k']}`, eta residual `{best['eta_signed_residual_vs_9_5']:.3e}`.\n"
        "- Honest read: 96/95 is a simple strict slot-count shadow, not an exact theorem; exact eta needs an irrational strict correction.\n"
        "- No false pass: no kernel identity, no legacy role transfer, no exact rational-slot eta theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
