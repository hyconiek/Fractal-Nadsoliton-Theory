#!/usr/bin/env python3
"""P1787 S737 strict unified EA/EH componentwise export + H1 run-contract checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"


def load(name: str) -> dict:
    path = GENERATED / name
    if not path.exists():
        return {"_missing": name}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p1760 = load("p1760_s710_strict_shared_background_family_contract_delivery_checkpoint.json")
    p1761 = load("p1761_s711_strict_boundary_term_control_clause_delivery_checkpoint.json")
    p1762 = load("p1762_s712_strict_boundary_control_contract_finalization_checkpoint.json")
    p1764 = load("p1764_s714_strict_covariant_nonproxy_ea_eh_explicit_export_checkpoint.json")
    p1786 = load("p1786_s736_strict_current_priority_nonproxy_closure_gate_matrix_checkpoint.json")

    gates = p1786.get("gate_matrix", {})

    out = {
        "checkpoint_id": "P1787_S737",
        "title": "STRICT_UNIFIED_EA_EH_COMPONENTWISE_EXPORT_AND_H1_RUN_CONTRACT",
        "scope": "strict-only local execution contract for G1+G2+G5",
        "input_reuse": [
            "p1760_s710_strict_shared_background_family_contract_delivery_checkpoint.json",
            "p1761_s711_strict_boundary_term_control_clause_delivery_checkpoint.json",
            "p1762_s712_strict_boundary_control_contract_finalization_checkpoint.json",
            "p1764_s714_strict_covariant_nonproxy_ea_eh_explicit_export_checkpoint.json",
            "p1786_s736_strict_current_priority_nonproxy_closure_gate_matrix_checkpoint.json",
        ],
        "strict_contract": {
            "route_policy": "STRICT_ONLY_NO_LEGACY_BRIDGE",
            "artifact_policy": "NONPROXY_EXPORTS_ONLY",
            "background_policy": "SINGLE_SHARED_BACKGROUND_FAMILY_REQUIRED",
            "index_policy": "SINGLE_INDEX_CONVENTION_REQUIRED",
            "pass_policy": "PASS_REQUIRES_EXPLICIT_COMPONENTWISE_RESIDUAL_ZERO_AND_WITNESS_LEDGER",
        },
        "run_order": ["G1_EA_NONPROXY_EXPLICIT_EXPORT", "G2_EH_NONPROXY_EXPLICIT_EXPORT", "G5_H1_4D_WEAK_FORM"],
        "gate_status_snapshot": {
            "G1": gates.get("G1_EA_NONPROXY_EXPLICIT_EXPORT", "OPEN"),
            "G2": gates.get("G2_EH_NONPROXY_EXPLICIT_EXPORT", "OPEN"),
            "G5": gates.get("G5_H1_4D_WEAK_FORM", "OPEN"),
            "boundary_control": gates.get("G4_BOUNDARY_TERM_CONTROL", "PARTIAL_OPEN"),
        },
        "dependency_snapshot": {
            "shared_background_contract_status": p1760.get("status", "MISSING"),
            "boundary_clause_status": p1761.get("status", "MISSING"),
            "boundary_finalization_status": p1762.get("status", "MISSING"),
            "ea_eh_export_status": p1764.get("status", "MISSING"),
        },
        "status_discipline": {
            "if_componentwise_missing": "OPEN_COMPONENTWISE_REQUIRED",
            "if_nonzero_residual": "OPEN_OBSTRUCTION_WITH_TRACE",
            "if_residual_zero_with_witness": "PASS_STRICT_LOCAL_H1_ZERO",
            "theorem_level_qg": "KEEP_OPEN_UNTIL_BRST_CUTKOSKY_BIANCHI_WARD_CLOSURE",
        },
        "next_honest_step": "Execute unified componentwise EA/EH export under one background+index freeze and publish H1 residual ledger without PASS claim unless explicit zero is produced.",
        "status": "UNIFIED_CONTRACT_READY_GATES_STILL_OPEN",
    }

    GENERATED.mkdir(parents=True, exist_ok=True)
    out_path = GENERATED / "p1787_s737_strict_unified_ea_eh_componentwise_export_and_h1_run_contract_checkpoint.json"
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(out_path)


if __name__ == "__main__":
    main()
