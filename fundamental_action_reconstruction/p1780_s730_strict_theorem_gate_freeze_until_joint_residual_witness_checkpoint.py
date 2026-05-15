#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1767 = GEN / "p1767_s717_strict_bianchi_ward_to_brst_cutkosky_gate_sequencing_checkpoint.json"
IN1778 = GEN / "p1778_s728_strict_componentwise_h1_metric_joint_execution_lock_checkpoint.json"
IN1779 = GEN / "p1779_s729_strict_current_priority_success_condition_tracker_checkpoint.json"
OUT = GEN / "p1780_s730_strict_theorem_gate_freeze_until_joint_residual_witness_checkpoint.json"


def main() -> None:
    p1767 = json.loads(IN1767.read_text(encoding="utf-8"))
    p1778 = json.loads(IN1778.read_text(encoding="utf-8"))
    p1779 = json.loads(IN1779.read_text(encoding="utf-8"))

    joint_allowed = p1778.get("gate_decision", {}).get("joint_run_allowed_now", False)

    payload = {
        "checkpoint": "P1780_S730",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_anchors": ["p1767_s717", "p1778_s728", "p1779_s729"],
        "theorem_gate_freeze": {
            "freeze_active": True,
            "freeze_scope": [
                "G_BW_promotion",
                "G_BRST_promotion",
                "G_CUT_promotion",
                "renormalization_theorem_gate_update",
                "background_independence_theorem_gate_update",
            ],
            "unfreeze_conditions": [
                "W1 accepted as FULL_EXPORT under P1774 contract",
                "joint H1+metric componentwise run executed on same background and index convention",
                "joint run published with PASS_ZERO or OBSTRUCTION_WITH_DIVERGENCE_TRACE and explicit witness trace",
            ],
            "current_blocker": p1778.get("gate_decision", {}),
        },
        "consistency_link": {
            "bw_brst_cut_sequence_ref": p1767.get("gate_chain", {}),
            "priority_success_tracker_ref": p1779.get("success_condition_tracker", {}),
        },
        "proofs_and_nonproofs": {
            "proven_now": [
                "A formal freeze contract is exported to prevent theorem-gate drift before joint residual witness publication.",
            ],
            "still_open": [
                "No new residual theorem witness was generated in this freeze step.",
                "QG theorem gates remain blocked pending joint witness execution.",
            ],
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Dowieźć W1 FULL_EXPORT, wykonać joint H1+metric run i dopiero po publikacji witness trace rozważać odblokowanie BW->BRST->CUT oraz dalsze theorem-gates.",
        "lay_summary": "Wprowadzono twardą blokadę formalną: dopóki nie ma wspólnego, jawnego wyniku kluczowych testów, nie wolno ogłaszać postępu w bramkach końcowych.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
