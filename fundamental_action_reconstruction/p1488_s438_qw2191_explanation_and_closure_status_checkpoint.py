#!/usr/bin/env python3
"""P1488 S4.38: transparent explanation artifact for QW-2191 closure status."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1487 = GEN / "p1487_s437_qw2191_selector_injection_test_summary.json"
P1486 = GEN / "p1486_s436_qw2191_symmetry_breaking_premise_probe_summary.json"
P1485 = GEN / "p1485_s435_qw2191_strict_closure_criterion_summary.json"

SUMMARY = GEN / "p1488_s438_qw2191_explanation_and_closure_status_summary.json"


def main() -> None:
    s1487 = json.loads(P1487.read_text(encoding="utf-8"))
    s1486 = json.loads(P1486.read_text(encoding="utf-8"))
    s1485 = json.loads(P1485.read_text(encoding="utf-8"))

    g0 = float(s1487["gaps"]["baseline_g0"])
    g1 = float(s1487["gaps"]["injected_g1"])
    delta_sb = float(s1487["delta_sb"])

    summary = {
        "packet": "P1488",
        "status": "EXPLAINED_LOCAL_STATUS",
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "qw2191_closed": False,
        "numbers_explained": {
            "W_SM_minus_W_GR_gap_before": g0,
            "W_SM_minus_W_GR_gap_after": g1,
            "delta_sb": delta_sb,
            "gap_reduction": g0 - g1,
        },
        "physical_reading": {
            "premise_is_controlled": bool(s1486["premise_verified_local"]),
            "closure_criterion_satisfied": bool(s1485["qw2191_closed"]),
            "selector_progress_local": g1 < g0,
        },
        "hard_truth": "Current repo evidence supports local selector progress, not theorem-level final closure of QW-2191.",
        "next_step_recommendation": "S4.39: perform kappa-window robustness sweep and require monotonic gap-reduction region before attempting strict closure theorem draft.",
        "layman_explanation": "To jak zmniejszanie luzu w mechanizmie: luz się zmniejszył, ale nie udowodniliśmy jeszcze, że mechanizm zawsze działa perfekcyjnie w całym zakresie ustawień.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("[P1488] status=EXPLAINED_LOCAL_STATUS")


if __name__ == "__main__":
    main()
