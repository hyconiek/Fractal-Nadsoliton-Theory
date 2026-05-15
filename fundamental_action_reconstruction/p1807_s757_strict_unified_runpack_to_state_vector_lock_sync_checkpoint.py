#!/usr/bin/env python3
"""P1807: sync P1806 gate verdict into strict state-vector lock fields."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
INP = GEN / "p1806_s756_strict_unified_nonproxy_ea_eh_elg_residual_runpack_contract_checkpoint.json"
OUT = GEN / "p1807_s757_strict_unified_runpack_to_state_vector_lock_sync_checkpoint.json"


def main() -> None:
    p1806 = json.loads(INP.read_text(encoding="utf-8"))
    tg1 = p1806.get("gate_vector", {}).get("TG1_BW", "OPEN_LOCKED_BY_UNIFIED_NONPROXY_RESIDUAL")

    if tg1 == "PASS_ZERO":
        sv6 = "PASS_ZERO"
        sv7 = "OPEN_LOCKED_BY_NILPOTENCY_WITNESS"
        sv8 = "OPEN_LOCKED_BY_TG2"
    else:
        sv6 = "OPEN_LOCKED_BY_UNIFIED_NONPROXY_RESIDUAL"
        sv7 = "OPEN_LOCKED_BY_TG1"
        sv8 = "OPEN_LOCKED_BY_TG2"

    checkpoint = {
        "packet_id": "P1807",
        "stage_id": "S757",
        "status": "PASS_ZERO_SYNC_APPLIED",
        "source_checkpoint": INP.name,
        "source_tg1_bw": tg1,
        "state_vector_lock_sync": {
            "SV6_BW_global": sv6,
            "SV7_BRST_global": sv7,
            "SV8_Cutkosky_global": sv8,
        },
        "anti_false_pass": {
            "manual_pass_promotion_forbidden": True,
            "brst_cut_require_explicit_witness_pack": True,
        },
    }

    OUT.write_text(json.dumps(checkpoint, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(str(OUT))


if __name__ == "__main__":
    main()
