#!/usr/bin/env python3
"""P1493 S4.43: export proof skeleton and explicit falsifier set for QW-2191."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1491 = GEN / "p1491_s441_qw2191_kappa_robustness_sweep_summary.json"
P1492 = GEN / "p1492_s442_qw2191_selector_theorem_candidate_summary.json"

SUMMARY = GEN / "p1493_s443_qw2191_proof_skeleton_and_falsifier_summary.json"


def main() -> None:
    s1491 = json.loads(P1491.read_text(encoding="utf-8"))
    s1492 = json.loads(P1492.read_text(encoding="utf-8"))

    rows = s1491["rows"]
    safe_rows = [r for r in rows if r["safe"]]

    l1 = all(abs(r["delta_sb"]) <= s1491["safety_margin"] for r in safe_rows)
    l2 = all(r["gap_after"] < r["gap_before"] for r in safe_rows)
    l3 = s1492["theorem_candidate"]["assumptions"]["A3_orientation_stable"]

    c1 = l1 and l2 and l3

    falsifier_rows = [
        {
            "kappa": r["kappa"],
            "reason": "safe_but_no_gap_reduction_or_orientation_flip"
        }
        for r in safe_rows
        if not (r["gap_after"] < r["gap_before"])
    ]

    summary = {
        "packet": "P1493",
        "status": "PASS_PROOF_SKELETON_LOCAL_ONLY" if c1 else "FAIL_PROOF_SKELETON_LOCAL_ONLY",
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "lemmas": {
            "L1_safe_margin_holds": l1,
            "L2_gap_reduction_holds": l2,
            "L3_orientation_stable": l3,
        },
        "conclusion_C1_local_selector_consistency_region": c1,
        "falsifier_dataset": falsifier_rows,
        "falsifier_rule": "Any safe-row violating L2 (or orientation stability failure) rejects C1.",
        "qw2191_closed": False,
        "next_step_recommendation": "S4.44: promote C1 into theorem draft with explicit quantifiers and attempt contradiction proof for uniqueness obstruction branch.",
        "layman_explanation": "To jak lista warunków bezpieczeństwa w samolocie: jeśli wszystkie są spełnione, lot lokalnie jest stabilny. Ale pełna certyfikacja (zamknięcie problemu) to jeszcze osobny etap.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1493] status={summary['status']} falsifier_count={len(falsifier_rows)}")


if __name__ == "__main__":
    main()
