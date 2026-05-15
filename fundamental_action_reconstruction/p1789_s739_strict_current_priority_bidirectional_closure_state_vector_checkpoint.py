#!/usr/bin/env python3
"""P1789 S739 strict current-priority bidirectional closure state vector checkpoint."""

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
    p1764 = load("p1764_s714_strict_covariant_nonproxy_ea_eh_explicit_export_checkpoint.json")
    p1765 = load("p1765_s715_strict_nonproxy_metric_elg_explicit_export_checkpoint.json")
    p1787 = load("p1787_s737_strict_unified_ea_eh_componentwise_export_and_h1_run_contract_checkpoint.json")
    p1788 = load("p1788_s738_strict_theorem_gate_prerequisite_lock_matrix_checkpoint.json")

    r1764 = p1764.get("readiness_update", {})
    theorem_locks = p1788.get("theorem_gate_locks", {})

    out = {
        "checkpoint_id": "P1789_S739",
        "title": "STRICT_CURRENT_PRIORITY_BIDIRECTIONAL_CLOSURE_STATE_VECTOR",
        "input_reuse": [
            "p1764_s714_strict_covariant_nonproxy_ea_eh_explicit_export_checkpoint.json",
            "p1765_s715_strict_nonproxy_metric_elg_explicit_export_checkpoint.json",
            "p1787_s737_strict_unified_ea_eh_componentwise_export_and_h1_run_contract_checkpoint.json",
            "p1788_s738_strict_theorem_gate_prerequisite_lock_matrix_checkpoint.json",
        ],
        "state_vector": {
            "SV1_E_A_mu_nonproxy_covariant_explicit": "OPEN_COMPONENTWISE_REQUIRED",
            "SV2_E_H_nonproxy_covariant_explicit": "OPEN_COMPONENTWISE_REQUIRED",
            "SV3_EL_g_nonproxy_explicit": "OPEN_COMPONENTWISE_REQUIRED",
            "SV4_boundary_term_control": "PARTIAL_OPEN",
            "SV5_H1_4D_weak_form": r1764.get("H1_4D_strict_local_readiness", "OPEN"),
            "SV6_Bianchi_Ward_global": theorem_locks.get("TG1_BIANCHI_WARD_GLOBAL", {}).get("status", "OPEN_LOCKED"),
            "SV7_BRST_global": theorem_locks.get("TG2_BRST_GLOBAL_NILPOTENCY", {}).get("status", "OPEN_LOCKED"),
            "SV8_Cutkosky_global": theorem_locks.get("TG3_CUTKOSKY_GLOBAL_UNITARITY", {}).get("status", "OPEN_LOCKED"),
        },
        "interpretation_discipline": {
            "local_pass_implies_global_pass": False,
            "scaffold_implies_full_export": False,
            "missing_residual_or_witness_means_open": True,
        },
        "minimal_bidirectional_closure_condition": [
            "single_background_family_componentwise_run_for_EA_EH_ELg",
            "explicit_H1_ledger_and_ELg_minus_E_munu_ledger",
            "strict_theorem_gate_sequence_BW_then_BRST_then_Cutkosky",
        ],
        "dependency_snapshot": {
            "p1764_status": p1764.get("status", "MISSING"),
            "p1765_status": p1765.get("status", "MISSING"),
            "p1787_status": p1787.get("status", "MISSING"),
            "p1788_status": p1788.get("status", "MISSING"),
        },
        "next_honest_step": "Execute coordinated run-pack for SV1..SV5 and publish only PASS_ZERO or OPEN_OBSTRUCTION_WITH_TRACE, then evaluate theorem gates without premature promotion.",
        "status": "CURRENT_PRIORITY_VECTOR_SYNCHRONIZED_GLOBAL_STILL_OPEN",
    }

    GENERATED.mkdir(parents=True, exist_ok=True)
    out_path = GENERATED / "p1789_s739_strict_current_priority_bidirectional_closure_state_vector_checkpoint.json"
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(out_path)


if __name__ == "__main__":
    main()
