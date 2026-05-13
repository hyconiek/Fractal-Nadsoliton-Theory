#!/usr/bin/env python3
"""P1498 S4.48: export final local strict-closure gate witness for QW-2191."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1497 = GEN / "p1497_s447_qw2191_closure_readiness_matrix_summary.json"
P1495 = GEN / "p1495_s445_qw2191_quantified_theorem_draft_summary.json"
P1494 = GEN / "p1494_s444_qw2191_contradiction_branch_summary.json"
P1496 = GEN / "p1496_s446_qw2191_cross_provider_replication_summary.json"

SUMMARY = GEN / "p1498_s448_qw2191_final_gate_witness_summary.json"


def main() -> None:
    s1497 = json.loads(P1497.read_text(encoding="utf-8"))
    s1495 = json.loads(P1495.read_text(encoding="utf-8"))
    s1494 = json.loads(P1494.read_text(encoding="utf-8"))
    s1496 = json.loads(P1496.read_text(encoding="utf-8"))

    checks = {
        "theorem_draft_pass": bool(s1495["theorem_holds_local"]),
        "contradiction_branch_pass": bool(s1494["contradiction_found"]),
        "cross_provider_pass": bool(s1496["replication_pass"]),
        "near_closure_ready": bool(s1497["near_closure_ready"]),
    }

    qw2191_closed_local = all(checks.values())

    witness = {
        "id": "W_gate_final_v1",
        "checks": checks,
        "qw2191_closed_local": qw2191_closed_local,
        "qw2191_closed_global": False,
        "assumption_scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
    }

    summary = {
        "packet": "P1498",
        "status": "PASS_FINAL_GATE_WITNESS_LOCAL_ONLY" if qw2191_closed_local else "FAIL_FINAL_GATE_WITNESS_LOCAL_ONLY",
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "witness": witness,
        "remaining_blockers_local": [] if qw2191_closed_local else ["local_gate_conditions_incomplete"],
        "qw2191_closed": qw2191_closed_local,
        "next_step_recommendation": "S4.49: build explicit F(nadsoliton)=>L_SM+L_GR mapping witness using closed-local selector gate as input constraint.",
        "layman_explanation": "To jest lokalny certyfikat, że wszystkie wymagane testy selektora zostały domknięte w naszym rygorze. To jeszcze nie globalna teoria wszystkiego, ale formalne zamknięcie lokalnego problemu QW-2191.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1498] status={summary['status']} local_closed={qw2191_closed_local}")


if __name__ == "__main__":
    main()
