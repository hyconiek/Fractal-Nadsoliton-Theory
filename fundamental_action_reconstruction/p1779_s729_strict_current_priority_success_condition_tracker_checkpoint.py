#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1764 = GEN / "p1764_s714_strict_covariant_nonproxy_ea_eh_explicit_export_checkpoint.json"
IN1765 = GEN / "p1765_s715_strict_nonproxy_metric_elg_explicit_export_checkpoint.json"
IN1778 = GEN / "p1778_s728_strict_componentwise_h1_metric_joint_execution_lock_checkpoint.json"
OUT = GEN / "p1779_s729_strict_current_priority_success_condition_tracker_checkpoint.json"


def main() -> None:
    p1764 = json.loads(IN1764.read_text(encoding="utf-8"))
    p1765 = json.loads(IN1765.read_text(encoding="utf-8"))
    p1778 = json.loads(IN1778.read_text(encoding="utf-8"))

    tracker = {
        "explicit_covariant_nonproxy_E_A_mu": "EXPORTED_OPERATOR_LEVEL",
        "explicit_covariant_nonproxy_E_H": "EXPORTED_OPERATOR_LEVEL",
        "metric_EL_g_export": "EXPORTED_OPERATOR_LEVEL",
        "boundary_term_control": "FINALIZED_CONTRACT_PRESENT",
        "H1_4D_weak_form_readiness": p1764.get("readiness_update", {}).get("H1_4D_strict_local_readiness", "OPEN"),
        "Bianchi_Ward_consistency": "BLOCKED_PENDING_W1_AND_COMPONENTWISE_DIVERGENCE",
        "BRST_Cutkosky_theorem_gates": "BLOCKED_PENDING_BW_PASS_AND_THEOREM_WITNESSES",
    }

    payload = {
        "checkpoint": "P1779_S729",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_anchors": ["p1764_s714", "p1765_s715", "p1778_s728"],
        "current_priority_tracker": tracker,
        "success_condition_tracker": {
            "K_strict_to_full_variational_structure": "OPEN_COMPONENTWISE_REQUIRED",
            "full_variational_structure_to_covariant_EOM_bundle": "PARTIAL_OPERATOR_LEVEL_DONE",
            "covariant_EOM_bundle_to_quantum_consistency_constraints": "OPEN_THEOREM_LEVEL_REQUIRED",
        },
        "blocking_fact": p1778.get("gate_decision", {}),
        "proofs_and_nonproofs": {
            "proven_now": [
                "Priority tracker aligned to currently exported strict nonproxy operator-level artifacts.",
                "Success condition decomposition exported with explicit OPEN/PARTIAL states (no solved claim).",
            ],
            "still_open": [
                "Componentwise H1 and EL_g-E_munu residual witnesses.",
                "Bianchi/Ward divergence closure at theorem gate quality.",
                "BRST nilpotency and Cutkosky unitarity theorem witnesses.",
            ],
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Zamknąć W1 FULL_EXPORT i uruchomić wspólny componentwise run H1+metric; dopiero po publikacji PASS_ZERO/OBSTRUCTION aktualizować BW, następnie BRST/Cutkosky.",
        "lay_summary": "Mamy mapę postępu do celu końcowego: część operatorowa jest gotowa, ale kluczowe dowody komponentowe i kwantowe są nadal przed nami.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
