#!/usr/bin/env python3
"""P1499 S4.49: compare local closure vs global-closure requirements for QW-2191."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1498 = GEN / "p1498_s448_qw2191_final_gate_witness_summary.json"
P1496 = GEN / "p1496_s446_qw2191_cross_provider_replication_summary.json"
P1495 = GEN / "p1495_s445_qw2191_quantified_theorem_draft_summary.json"
P1493 = GEN / "p1493_s443_qw2191_proof_skeleton_and_falsifier_summary.json"

SUMMARY = GEN / "p1499_s449_qw2191_global_closure_requirements_summary.json"


def main() -> None:
    s1498 = json.loads(P1498.read_text(encoding="utf-8"))
    s1496 = json.loads(P1496.read_text(encoding="utf-8"))
    s1495 = json.loads(P1495.read_text(encoding="utf-8"))
    s1493 = json.loads(P1493.read_text(encoding="utf-8"))

    local_closed = bool(s1498["witness"]["qw2191_closed_local"])

    requirements = {
        "R1_cross_provider_replication": bool(s1496["replication_pass"]),
        "R2_strict_internal_selector_source_exported": False,
        "R3_quantified_theorem_text_ready": bool(s1495["theorem_holds_local"]),
        "R4_falsifier_set_non_obstructing": len(s1493["falsifier_dataset"]) == 0,
        "R5_F_to_LSM_LGR_mapping_witness_exported": False,
    }

    global_closed = all(requirements.values())

    summary = {
        "packet": "P1499",
        "status": "PASS_GLOBAL_CLOSURE_REQUIREMENTS_TRACKED_LOCAL_ONLY",
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "local_closed": local_closed,
        "global_closed": global_closed,
        "requirements": requirements,
        "missing_for_global": [k for k, v in requirements.items() if not v],
        "next_step_recommendation": "S4.50: export strict internal selector-source object (R2) and build explicit F=>L_SM+L_GR mapping witness (R5) under same selector assumptions.",
        "layman_explanation": "Lokalne zamknięcie to jak zaliczenie testów w jednym laboratorium. Globalne zamknięcie wymaga jeszcze pełnej certyfikacji: niezależnych potwierdzeń i pełnego połączenia z całą teorią SM+GR.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1499] local_closed={local_closed} global_closed={global_closed} missing={summary['missing_for_global']}")


if __name__ == "__main__":
    main()
