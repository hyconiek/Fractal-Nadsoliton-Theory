#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1716 = GEN / "p1716_s666_strict_metric_index_convention_normalization_audit_checkpoint.json"
IN1760 = GEN / "p1760_s710_strict_shared_background_family_contract_delivery_checkpoint.json"
IN1775 = GEN / "p1775_s725_strict_w1_hr2_full_export_delivery_attempt_checkpoint.json"
IN1777 = GEN / "p1777_s727_strict_forward_reverse_state_vector_sync_with_priority_blocker_matrix_checkpoint.json"
OUT = GEN / "p1778_s728_strict_componentwise_h1_metric_joint_execution_lock_checkpoint.json"


def main() -> None:
    p1716 = json.loads(IN1716.read_text(encoding="utf-8"))
    p1760 = json.loads(IN1760.read_text(encoding="utf-8"))
    p1775 = json.loads(IN1775.read_text(encoding="utf-8"))
    p1777 = json.loads(IN1777.read_text(encoding="utf-8"))

    w1_ok = p1775.get("acceptance_verdict") == "PASS_ZERO"
    h1_status = p1777.get("synced_state_vector", {}).get("reverse_chain", {}).get("H1_cross_variation", "OPEN")

    payload = {
        "checkpoint": "P1778_S728",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_anchors": ["p1716_s666", "p1760_s710", "p1775_s725", "p1777_s727"],
        "joint_execution_lock": {
            "scope": "R1_phase_1_H1_and_phase_2_metric_ELg_minus_Emunu",
            "required_same_background_family": True,
            "required_same_index_sign_convention": True,
            "required_same_residual_basis": ["B1", "B2", "B3", "C1", "C2"],
            "admissible_outcomes": ["PASS_ZERO", "OBSTRUCTION_WITH_DIVERGENCE_TRACE"],
            "forbidden": [
                "mixed-background evaluation across H1 and metric residual",
                "index/sign drift between H1 and metric residual runs",
                "declaring PASS without explicit residual witness",
            ],
        },
        "readiness_snapshot": {
            "index_convention_audit_anchor": p1716.get("checkpoint", "P1716_S666"),
            "shared_background_contract_anchor": p1760.get("checkpoint", "P1760_S710"),
            "w1_full_export_gate": p1775.get("acceptance_verdict", "OPEN"),
            "h1_status": h1_status,
        },
        "gate_decision": {
            "joint_run_allowed_now": w1_ok,
            "reason": "W1 not FULL_EXPORT yet" if not w1_ok else "All preconditions met for strict joint run",
            "next_blocker": "W1_B2_B3_projection_and_HR2_divergence_normalization" if not w1_ok else "EXECUTE_JOINT_RUN_NOW",
        },
        "proofs_and_nonproofs": {
            "proven_now": [
                "A single-lock contract for same-background/same-convention H1+metric componentwise execution is exported.",
                "Current gate remains blocked because W1 is not accepted as FULL_EXPORT.",
            ],
            "still_open": [
                "No H1 componentwise residual witness yet.",
                "No EL_g-E_munu componentwise residual witness yet.",
                "No Bianchi/Ward theorem-level divergence witness yet.",
            ],
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Domknąć W1 do FULL_EXPORT i natychmiast uruchomić joint H1+metric run na tej samej rodzinie teł i tej samej konwencji indeksowej; opublikować wyłącznie PASS_ZERO albo OBSTRUCTION_WITH_DIVERGENCE_TRACE.",
        "lay_summary": "Przygotowano rygorystyczny plan jednego wspólnego uruchomienia dwóch kluczowych testów. Jednak start jest zablokowany, dopóki W1 nie będzie formalnie domknięte.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
