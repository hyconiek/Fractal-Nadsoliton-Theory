#!/usr/bin/env python3
"""P1794 S744 strict SV1->SV5 verdict-to-state update protocol checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"


def load(name: str) -> dict:
    p = GENERATED / name
    if not p.exists():
        return {"_missing": name}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    p1789 = load("p1789_s739_strict_current_priority_bidirectional_closure_state_vector_checkpoint.json")
    p1788 = load("p1788_s738_strict_theorem_gate_prerequisite_lock_matrix_checkpoint.json")
    p1793 = load("p1793_s743_sv1_sv5_evidence_intake_validation_result.json")

    verdict = p1793.get("verdict", "OPEN_OBSTRUCTION_WITH_TRACE")
    obstruction_tags = p1793.get("obstruction_tags", [])

    if verdict == "PASS_ZERO":
        sv15_status = "PASS_STRICT_LOCAL_EVIDENCE_ACCEPTED"
    else:
        sv15_status = "OPEN_OBSTRUCTION_WITH_TRACE"

    base_sv = p1789.get("state_vector", {})
    theorem_locks = p1788.get("theorem_gate_locks", {})

    out = {
        "checkpoint_id": "P1794_S744",
        "title": "STRICT_SV1_TO_SV5_VERDICT_TO_STATE_UPDATE_PROTOCOL",
        "input_reuse": [
            "p1789_s739_strict_current_priority_bidirectional_closure_state_vector_checkpoint.json",
            "p1788_s738_strict_theorem_gate_prerequisite_lock_matrix_checkpoint.json",
            "p1793_s743_sv1_sv5_evidence_intake_validation_result.json",
        ],
        "source_verdict": verdict,
        "source_obstruction_tags": obstruction_tags,
        "updated_state_vector": {
            "SV1_E_A_mu_nonproxy_covariant_explicit": sv15_status,
            "SV2_E_H_nonproxy_covariant_explicit": sv15_status,
            "SV3_EL_g_nonproxy_explicit": sv15_status,
            "SV4_boundary_term_control": sv15_status,
            "SV5_H1_4D_weak_form": sv15_status,
            "SV6_Bianchi_Ward_global": theorem_locks.get("TG1_BIANCHI_WARD_GLOBAL", {}).get("status", base_sv.get("SV6_Bianchi_Ward_global", "OPEN_LOCKED")),
            "SV7_BRST_global": theorem_locks.get("TG2_BRST_GLOBAL_NILPOTENCY", {}).get("status", base_sv.get("SV7_BRST_global", "OPEN_LOCKED")),
            "SV8_Cutkosky_global": theorem_locks.get("TG3_CUTKOSKY_GLOBAL_UNITARITY", {}).get("status", base_sv.get("SV8_Cutkosky_global", "OPEN_LOCKED")),
        },
        "global_promotion_policy": "NO_AUTOPROMOTION_OF_SV6_TO_SV8",
        "next_honest_step": "Use this deterministic state update as input for subsequent BW/BRST/Cutkosky witness runs, without theorem-level claims from local acceptance alone.",
        "status": "STATE_UPDATE_PROTOCOL_APPLIED",
    }

    out_path = GENERATED / "p1794_s744_strict_sv1_to_sv5_verdict_to_state_update_protocol_checkpoint.json"
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(out_path)


if __name__ == "__main__":
    main()
