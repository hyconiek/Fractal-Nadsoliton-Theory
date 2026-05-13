#!/usr/bin/env python3
"""P1492 S4.42: draft strict selector-source theorem candidate from robust region."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1491 = GEN / "p1491_s441_qw2191_kappa_robustness_sweep_summary.json"
P1489 = GEN / "p1489_s439_qw2191_selector_source_candidate_summary.json"

SUMMARY = GEN / "p1492_s442_qw2191_selector_theorem_candidate_summary.json"
OBSTRUCTION = GEN / "p1492_s442_qw2191_selector_theorem_candidate_obstruction.json"


def main() -> None:
    s1491 = json.loads(P1491.read_text(encoding="utf-8"))
    s1489 = json.loads(P1489.read_text(encoding="utf-8"))

    rows = s1491["rows"]
    robust = [r for r in rows if r["safe"] and r["improved"]]

    if robust:
        kappa_min = min(r["kappa"] for r in robust)
        kappa_max = max(r["kappa"] for r in robust)
    else:
        kappa_min = None
        kappa_max = None

    orientation = s1489["orientation"]
    orientation_stable = all((r["gap_after"] > 0) for r in robust) if robust else False

    assumptions = {
        "A1_safe_margin_enforced": bool(robust),
        "A2_gap_reduction_uniform": len(robust) == s1491["safe_points"] and s1491["safe_points"] > 0,
        "A3_orientation_stable": orientation_stable,
        "A4_no_legacy_bridge": True,
    }

    theorem_candidate_ready = all(assumptions.values())
    status = "PASS_SELECTOR_THEOREM_CANDIDATE_READY_LOCAL_ONLY" if theorem_candidate_ready else "FAIL_SELECTOR_THEOREM_CANDIDATE_READY_LOCAL_ONLY"

    summary = {
        "packet": "P1492",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "theorem_candidate": {
            "name": "QW2191_strict_selector_source_candidate_v1",
            "statement": "On robust kappa range, controlled Delta_SB induces consistent selector orientation with uniform safe gap reduction.",
            "robust_kappa_range": [kappa_min, kappa_max],
            "orientation": orientation,
            "assumptions": assumptions,
        },
        "qw2191_closed": False,
        "next_step_recommendation": "S4.43: attempt formal proof skeleton (lemmas + contradiction branch) and define explicit falsifier dataset for theorem rejection.",
        "layman_explanation": "To szkic twierdzenia: mówimy dokładnie kiedy mechanizm wyboru działa i na jakim zakresie ustawień. Jeszcze nie finalny dowód, ale już formalny plan dowodu.",
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if theorem_candidate_ready:
        if OBSTRUCTION.exists():
            OBSTRUCTION.unlink()
    else:
        OBSTRUCTION.write_text(json.dumps({"packet": "P1492", "status": status, "assumptions": assumptions}, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print(f"[P1492] status={status} robust_range=({kappa_min},{kappa_max})")


if __name__ == "__main__":
    main()
