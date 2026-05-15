#!/usr/bin/env python3
"""P1806: validate unified nonproxy residual run-pack and emit checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
INPUT = GEN / "p1806_s756_unified_nonproxy_ea_eh_elg_runpack_input.json"
OUT = GEN / "p1806_s756_strict_unified_nonproxy_ea_eh_elg_residual_runpack_contract_checkpoint.json"


def verdict_zero(section: dict) -> bool:
    return section.get("verdict") == "PASS_ZERO"


def main() -> None:
    payload = json.loads(INPUT.read_text(encoding="utf-8"))

    required = ["shared_background_family_id", "manifest_id", "boundary_control_clause_id", "weak_form_h1_4d_ready", "ea_mu_residual", "eh_residual", "elg_residual"]
    missing = [k for k in required if k not in payload]

    tg1 = "OPEN_LOCKED_BY_UNIFIED_NONPROXY_RESIDUAL"
    if not missing and payload.get("weak_form_h1_4d_ready") is True:
        if all(verdict_zero(payload[k]) for k in ["ea_mu_residual", "eh_residual", "elg_residual"]):
            tg1 = "PASS_ZERO"

    checkpoint = {
        "packet_id": "P1806",
        "stage_id": "S756",
        "status": "PASS_ZERO" if tg1 == "PASS_ZERO" else "OPEN_OBSTRUCTION_WITH_TRACE",
        "missing_required_fields": missing,
        "gate_vector": {
            "TG1_BW": tg1,
            "TG2_BRST": "UNLOCKED" if tg1 == "PASS_ZERO" else "LOCKED_BY_TG1",
            "TG3_CUT": "LOCKED_BY_TG2",
        },
        "checks": {
            "weak_form_h1_4d_ready": payload.get("weak_form_h1_4d_ready") is True,
            "ea_mu_zero": verdict_zero(payload.get("ea_mu_residual", {})),
            "eh_zero": verdict_zero(payload.get("eh_residual", {})),
            "elg_zero": verdict_zero(payload.get("elg_residual", {})),
        },
    }

    OUT.write_text(json.dumps(checkpoint, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(str(OUT))


if __name__ == "__main__":
    main()
