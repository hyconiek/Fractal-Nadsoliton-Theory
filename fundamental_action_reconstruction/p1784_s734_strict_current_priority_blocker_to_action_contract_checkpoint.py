#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1782 = GEN / "p1782_s732_strict_priority_closure_gap_matrix_checkpoint.json"
IN1783 = GEN / "p1783_s733_strict_priority_gate_exit_criteria_matrix_checkpoint.json"
OUT = GEN / "p1784_s734_strict_current_priority_blocker_to_action_contract_checkpoint.json"


def main() -> None:
    p1782 = json.loads(IN1782.read_text(encoding="utf-8"))
    p1783 = json.loads(IN1783.read_text(encoding="utf-8"))

    blockers = p1782.get("priority_closure_gap_matrix", {})
    exits = p1783.get("priority_gate_exit_criteria_matrix", {})

    action_contract = {
        "A1_W1_full_export_close": {
            "trigger": exits.get("W1_FULL_EXPORT", {}).get("met", False) is False,
            "required_delivery": [
                "B2_B3_projection_coefficients_explicit",
                "HR2_divergence_normalization_map_explicit",
            ],
            "success_flag": "W1_FULL_EXPORT=True",
        },
        "A2_joint_H1_metric_componentwise_run": {
            "trigger": exits.get("JOINT_H1_METRIC_COMPONENTWISE_RUN", {}).get("met", False) is False,
            "required_delivery": [
                "shared_background_family_locked",
                "shared_index_sign_convention_locked",
                "shared_basis_B1_B2_B3_C1_C2_locked",
            ],
            "result_policy": ["PASS_ZERO", "OBSTRUCTION_WITH_DIVERGENCE_TRACE"],
        },
        "A3_bw_brst_cut_progression": {
            "trigger": exits.get("THEOREM_GATE_FREEZE_RELEASE", {}).get("met", False) is False,
            "required_delivery": [
                "joint_witness_trace_published",
                "formal_freeze_release_conditions_met",
            ],
            "allowed_after": ["A1", "A2"],
        },
    }

    payload = {
        "checkpoint": "P1784_S734",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_anchors": ["p1782_s732", "p1783_s733"],
        "current_priority_blockers": blockers,
        "blocker_to_action_contract": action_contract,
        "global_exit_ready": p1783.get("global_exit_ready", False),
        "proofs_and_nonproofs": {
            "proven_now": [
                "Each current blocker is mapped to an explicit next-action contract with trigger and success flag semantics.",
            ],
            "still_open": [
                "No new residual witness produced in this mapping step.",
            ],
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Wykonać kolejno A1->A2->A3 zgodnie z kontraktem, publikując jawne witnessy i utrzymując klasyfikację PASS_ZERO/OBSTRUCTION_WITH_DIVERGENCE_TRACE.",
        "lay_summary": "To plan działań przypięty do realnych blokad: każdy problem ma przypisany konkretny następny krok i warunek sukcesu.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
