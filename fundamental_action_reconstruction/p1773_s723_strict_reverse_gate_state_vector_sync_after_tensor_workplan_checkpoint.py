#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1766 = GEN / "p1766_s716_strict_forward_reverse_state_vector_update_with_bianchi_ward_gate_checkpoint.json"
IN1772 = GEN / "p1772_s722_strict_gbw_tensor_closure_workplan_checkpoint.json"
OUT = GEN / "p1773_s723_strict_reverse_gate_state_vector_sync_after_tensor_workplan_checkpoint.json"


def main() -> None:
    p1766 = json.loads(IN1766.read_text(encoding="utf-8"))
    p1772 = json.loads(IN1772.read_text(encoding="utf-8"))

    payload = {
        "checkpoint": "P1773_S723",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_anchors": ["p1766_s716", "p1772_s722"],
        "previous_state_vector": p1766.get("updated_state_vector", {}),
        "updated_reverse_state_vector": {
            "G_BW": "OPEN_OBSTRUCTION_ACTIVE_W1_W4_REQUIRED",
            "G_BRST": "BLOCKED_BY_G_BW",
            "G_CUT": "BLOCKED_BY_G_BW_AND_G_BRST",
            "global_helmholtz_integrability": "OPEN",
            "qg_theorem_promotable": False,
        },
        "workplan_binding": p1772.get("tensor_closure_workplan", {}),
        "promotion_guard": {
            "forbid_qg_promotion_if": [
                "G_BW != PASS_ZERO",
                "W1..W4 not all delivered as FULL_EXPORT",
            ],
            "classification_lock": ["LOCAL/GLOBAL", "REDUCED/NONPROXY", "SCAFFOLD/FULL_EXPORT"],
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Dostarczyć pełne exporty W1-W4 i wykonać kolejną próbę G_BW; dopiero potem aktualizować BRST/Cutkosky readiness.",
        "lay_summary": "Zaktualizowano mapę postępu: wiemy dokładnie, że bez dokończenia czterech bloków matematycznych nie można przejść do kolejnych testów teorii.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
