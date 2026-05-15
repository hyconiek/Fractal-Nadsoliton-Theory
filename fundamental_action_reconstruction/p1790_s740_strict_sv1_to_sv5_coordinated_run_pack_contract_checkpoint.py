#!/usr/bin/env python3
"""P1790 S740 strict SV1->SV5 coordinated run-pack contract checkpoint."""

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

    sv = p1789.get("state_vector", {})
    locks = p1788.get("theorem_gate_locks", {})

    out = {
        "checkpoint_id": "P1790_S740",
        "title": "STRICT_SV1_TO_SV5_COORDINATED_RUN_PACK_CONTRACT",
        "input_reuse": [
            "p1789_s739_strict_current_priority_bidirectional_closure_state_vector_checkpoint.json",
            "p1788_s738_strict_theorem_gate_prerequisite_lock_matrix_checkpoint.json",
        ],
        "execution_contract": {
            "strict_only": True,
            "nonproxy_only": True,
            "single_background_family": True,
            "single_index_convention": True,
            "single_boundary_control_clause": True,
            "decision_policy": ["PASS_ZERO", "OPEN_OBSTRUCTION_WITH_TRACE"],
        },
        "run_order_no_skip": ["SV1", "SV2", "SV3", "SV4", "SV5"],
        "current_sv_snapshot": {
            "SV1": sv.get("SV1_E_A_mu_nonproxy_covariant_explicit", "OPEN"),
            "SV2": sv.get("SV2_E_H_nonproxy_covariant_explicit", "OPEN"),
            "SV3": sv.get("SV3_EL_g_nonproxy_explicit", "OPEN"),
            "SV4": sv.get("SV4_boundary_term_control", "OPEN"),
            "SV5": sv.get("SV5_H1_4D_weak_form", "OPEN"),
        },
        "pass_requirements": [
            "componentwise_ledger_E_A_mu",
            "componentwise_ledger_E_H",
            "componentwise_ledger_EL_g",
            "componentwise_H1_zero_vector_ledger",
            "boundary_control_confirmed_on_same_freeze",
        ],
        "theorem_gate_promotion_after_run_pack": {
            "allowed": False,
            "reason": "BW_BRST_CUT_still_locked_until_global_witnesses",
            "locks": {
                "TG1": locks.get("TG1_BIANCHI_WARD_GLOBAL", {}).get("status", "OPEN_LOCKED"),
                "TG2": locks.get("TG2_BRST_GLOBAL_NILPOTENCY", {}).get("status", "OPEN_LOCKED"),
                "TG3": locks.get("TG3_CUTKOSKY_GLOBAL_UNITARITY", {}).get("status", "OPEN_LOCKED"),
            },
        },
        "next_honest_step": "Execute SV1->SV5 coordinated run-pack and update only SV1..SV5 from explicit ledgers, keeping SV6..SV8 unchanged until global theorem witnesses are exported.",
        "status": "RUN_PACK_CONTRACT_READY_EXECUTION_REQUIRED",
    }

    GENERATED.mkdir(parents=True, exist_ok=True)
    out_path = GENERATED / "p1790_s740_strict_sv1_to_sv5_coordinated_run_pack_contract_checkpoint.json"
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(out_path)


if __name__ == "__main__":
    main()
