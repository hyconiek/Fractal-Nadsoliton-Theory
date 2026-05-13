#!/usr/bin/env python3
"""P1497 S4.47: closure-readiness matrix for QW-2191 strict path."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1495 = GEN / "p1495_s445_qw2191_quantified_theorem_draft_summary.json"
P1494 = GEN / "p1494_s444_qw2191_contradiction_branch_summary.json"
P1496 = GEN / "p1496_s446_qw2191_cross_provider_replication_summary.json"
P1485 = GEN / "p1485_s435_qw2191_strict_closure_criterion_summary.json"

SUMMARY = GEN / "p1497_s447_qw2191_closure_readiness_matrix_summary.json"


def main() -> None:
    s1495 = json.loads(P1495.read_text(encoding="utf-8"))
    s1494 = json.loads(P1494.read_text(encoding="utf-8"))
    s1496 = json.loads(P1496.read_text(encoding="utf-8"))
    s1485 = json.loads(P1485.read_text(encoding="utf-8"))

    matrix = {
        "pillar_quantified_theorem": bool(s1495["theorem_holds_local"]),
        "pillar_contradiction_branch": bool(s1494["contradiction_found"]),
        "pillar_cross_provider_replication": bool(s1496["replication_pass"]),
        "pillar_strict_closure_gate": bool(s1485["qw2191_closed"]),
    }

    near_closure_ready = matrix["pillar_quantified_theorem"] and matrix["pillar_contradiction_branch"] and matrix["pillar_cross_provider_replication"]
    final_closed = near_closure_ready and matrix["pillar_strict_closure_gate"]

    blockers = []
    if not matrix["pillar_strict_closure_gate"]:
        blockers.append("strict_internal_selector_source_or_verified_symmetry_breaking_not_exported_as_final_gate_pass")

    summary = {
        "packet": "P1497",
        "status": "PASS_NEAR_CLOSURE_READINESS_LOCAL_ONLY" if near_closure_ready else "FAIL_NEAR_CLOSURE_READINESS_LOCAL_ONLY",
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "matrix": matrix,
        "near_closure_ready": near_closure_ready,
        "final_closed": final_closed,
        "remaining_blockers": blockers,
        "qw2191_closed": final_closed,
        "next_step_recommendation": "S4.48: export final strict-closure gate witness object proving selector-source criterion as discharged under explicit assumptions.",
        "layman_explanation": "Mamy prawie wszystko: teoria działa, replikuje się i przechodzi test sprzeczności. Został ostatni formalny krok: oficjalny certyfikat, że blokada selektora jest naprawdę zdjęta.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1497] status={summary['status']} final_closed={final_closed}")


if __name__ == "__main__":
    main()
