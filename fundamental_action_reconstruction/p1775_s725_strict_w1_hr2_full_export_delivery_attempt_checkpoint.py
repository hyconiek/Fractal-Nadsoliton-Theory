#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1774 = GEN / "p1774_s724_strict_w1_w4_delivery_sequence_and_acceptance_contract_checkpoint.json"
OUT = GEN / "p1775_s725_strict_w1_hr2_full_export_delivery_attempt_checkpoint.json"


def main() -> None:
    p1774 = json.loads(IN1774.read_text(encoding="utf-8"))

    payload = {
        "checkpoint": "P1775_S725",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_anchor": "p1774_s724",
        "delivery_target": "W1_H_R2_componentwise",
        "w1_delivery_attempt": {
            "classification": "PARTIAL_NOT_FULL_EXPORT",
            "componentwise_symbol_list_exported": True,
            "index_sign_convention_lock_proof": True,
            "residual_basis_projection_map_to_B1_B2_B3_C1_C2": "PARTIAL",
            "missing_items": [
                "explicit projection coefficients for B2/B3 in mixed-derivative terms",
                "fully normalized divergence-side contraction map for H_R2 branch",
            ],
        },
        "acceptance_verdict": "OBSTRUCTION_W1_NOT_YET_FULL_EXPORT",
        "gate_policy_effect": {
            "G_BW_rerun_allowed": False,
            "reason": "W1 not accepted as FULL_EXPORT under P1774 acceptance contract",
            "G_BRST": "BLOCKED",
            "G_CUT": "BLOCKED",
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Dowieźć brakujące współczynniki projekcji B2/B3 i normalizację kontrakcji dywergencji dla W1, następnie ponowić odbiór W1 jako FULL_EXPORT.",
        "lay_summary": "Pierwsza dostawa W1 jest blisko, ale jeszcze niepełna. Trzeba dopracować dwa brakujące elementy mapy projekcji, zanim wolno iść dalej.",
        "acceptance_contract_ref": p1774.get("acceptance_contract", {}),
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
