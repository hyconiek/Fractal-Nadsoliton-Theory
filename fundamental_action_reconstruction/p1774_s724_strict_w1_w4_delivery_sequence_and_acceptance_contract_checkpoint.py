#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1772 = GEN / "p1772_s722_strict_gbw_tensor_closure_workplan_checkpoint.json"
IN1773 = GEN / "p1773_s723_strict_reverse_gate_state_vector_sync_after_tensor_workplan_checkpoint.json"
OUT = GEN / "p1774_s724_strict_w1_w4_delivery_sequence_and_acceptance_contract_checkpoint.json"


def main() -> None:
    p1772 = json.loads(IN1772.read_text(encoding="utf-8"))
    p1773 = json.loads(IN1773.read_text(encoding="utf-8"))

    payload = {
        "checkpoint": "P1774_S724",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_anchors": ["p1772_s722", "p1773_s723"],
        "w1_w4_delivery_sequence": [
            "W1_H_R2_componentwise",
            "W2_H_Ric2_componentwise",
            "W3_H_Riem2_componentwise",
            "W4_CT_tensor_basis",
        ],
        "acceptance_contract": {
            "for_each_Wi_require": [
                "FULL_EXPORT_CLASSIFICATION",
                "componentwise_symbol_list",
                "index_sign_convention_lock_proof",
                "residual_basis_projection_map_to_B1_B2_B3_C1_C2",
            ],
            "global_require_before_GBW_rerun": [
                "All W1..W4 accepted as FULL_EXPORT",
                "No change in background family",
                "No change in residual basis B1/B2/B3/C1/C2",
            ],
            "allowed_final_verdicts": ["PASS_ZERO", "OBSTRUCTION_WITH_DIVERGENCE_TRACE"],
        },
        "gate_sync_snapshot": p1773.get("updated_reverse_state_vector", {}),
        "status": "READY_FOR_W1_W4_DELIVERY_PHASE",
        "no_false_pass_claim": True,
        "next_honest_step": "Dostarczyć W1 jako FULL_EXPORT z mapą projekcji do B1/B2/B3/C1/C2 i dopiero potem przejść do W2.",
        "lay_summary": "Ustalono jasne kryteria odbioru dla czterech brakujących bloków. Bez ich pełnego odbioru nie wolno uruchomić finałowego testu.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
