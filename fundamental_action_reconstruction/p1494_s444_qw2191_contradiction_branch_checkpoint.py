#!/usr/bin/env python3
"""P1494 S4.44: contradiction-branch check for QW-2191 uniqueness obstruction."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1493 = GEN / "p1493_s443_qw2191_proof_skeleton_and_falsifier_summary.json"
P1492 = GEN / "p1492_s442_qw2191_selector_theorem_candidate_summary.json"

SUMMARY = GEN / "p1494_s444_qw2191_contradiction_branch_summary.json"


def main() -> None:
    s1493 = json.loads(P1493.read_text(encoding="utf-8"))
    s1492 = json.loads(P1492.read_text(encoding="utf-8"))

    lemmas = s1493["lemmas"]
    robust_range = s1492["theorem_candidate"]["robust_kappa_range"]

    h0_nonexistence_selector_source = True

    evidence_selector_consistency = bool(
        lemmas["L1_safe_margin_holds"]
        and lemmas["L2_gap_reduction_holds"]
        and lemmas["L3_orientation_stable"]
    )

    contradiction_found = bool(h0_nonexistence_selector_source and evidence_selector_consistency)

    summary = {
        "packet": "P1494",
        "status": "PASS_CONTRADICTION_BRANCH_LOCAL_ONLY" if contradiction_found else "FAIL_CONTRADICTION_BRANCH_LOCAL_ONLY",
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "hypothesis_H0": "No consistent selector source exists on robust kappa range.",
        "robust_kappa_range": robust_range,
        "evidence_selector_consistency": evidence_selector_consistency,
        "contradiction_found": contradiction_found,
        "qw2191_closed": False,
        "next_step_recommendation": "S4.45: convert contradiction branch into explicit theorem text with quantified assumptions and independent verification script.",
        "layman_explanation": "Sprawdzamy hipotezę 'nie ma mechanizmu wyboru'. Dane pokazują stabilny mechanizm w bezpiecznym zakresie, więc ta hipoteza lokalnie przeczy faktom.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1494] status={summary['status']} contradiction_found={contradiction_found}")


if __name__ == "__main__":
    main()
