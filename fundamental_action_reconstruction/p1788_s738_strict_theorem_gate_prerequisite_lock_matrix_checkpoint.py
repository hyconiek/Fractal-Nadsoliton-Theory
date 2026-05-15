#!/usr/bin/env python3
"""P1788 S738 strict theorem-gate prerequisite lock matrix checkpoint."""

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
    p1767 = load("p1767_s717_strict_bianchi_ward_to_brst_cutkosky_gate_sequencing_checkpoint.json")
    p1787 = load("p1787_s737_strict_unified_ea_eh_componentwise_export_and_h1_run_contract_checkpoint.json")

    readiness_1764 = p1764.get("readiness_update", {})
    readiness_1765 = p1765.get("readiness_update", {})

    out = {
        "checkpoint_id": "P1788_S738",
        "title": "STRICT_THEOREM_GATE_PREREQUISITE_LOCK_MATRIX",
        "input_reuse": [
            "p1764_s714_strict_covariant_nonproxy_ea_eh_explicit_export_checkpoint.json",
            "p1765_s715_strict_nonproxy_metric_elg_explicit_export_checkpoint.json",
            "p1767_s717_strict_bianchi_ward_to_brst_cutkosky_gate_sequencing_checkpoint.json",
            "p1787_s737_strict_unified_ea_eh_componentwise_export_and_h1_run_contract_checkpoint.json",
        ],
        "theorem_gate_locks": {
            "TG1_BIANCHI_WARD_GLOBAL": {
                "required": "componentwise_divergence_trace_on_shared_background_family",
                "status": "OPEN_LOCKED_BY_COMPONENTWISE_DIVERGENCE",
            },
            "TG2_BRST_GLOBAL_NILPOTENCY": {
                "required": "explicit_nilpotency_witness_on_current_nonproxy_bundle",
                "status": "OPEN_LOCKED_BY_NILPOTENCY_WITNESS",
            },
            "TG3_CUTKOSKY_GLOBAL_UNITARITY": {
                "required": "explicit_cut_compatibility_witness_without_ghost_pole_contradiction",
                "status": "OPEN_LOCKED_BY_UNITARITY_WITNESS",
            },
        },
        "local_to_global_promotion_policy": {
            "H1_local_pass_sufficient_for_global": False,
            "ELg_local_pass_sufficient_for_global": False,
            "global_promotion_requires": ["TG1", "TG2", "TG3"],
        },
        "evidence_snapshot": {
            "p1764_H1_4D": readiness_1764.get("H1_4D_strict_local_readiness", "OPEN"),
            "p1765_metric_local": readiness_1765.get("metric_residual_gate", "OPEN"),
            "p1764_BW": readiness_1764.get("Bianchi_Ward_consistency_gate", "OPEN"),
            "p1764_BRST": readiness_1764.get("BRST_gate", "OPEN"),
            "p1764_CUT": readiness_1764.get("Cutkosky_gate", "OPEN"),
            "p1787_status": p1787.get("status", "MISSING"),
            "p1767_status": p1767.get("status", "MISSING"),
        },
        "next_honest_step": "After unified G1/G2/G5 execution, run global gate-pack in strict order: BW divergence witness -> BRST nilpotency witness -> Cutkosky unitarity witness, with OPEN maintained until all three are closed.",
        "status": "THEOREM_GATES_LOCKED_KEEP_GLOBAL_OPEN",
    }

    GENERATED.mkdir(parents=True, exist_ok=True)
    out_path = GENERATED / "p1788_s738_strict_theorem_gate_prerequisite_lock_matrix_checkpoint.json"
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(out_path)


if __name__ == "__main__":
    main()
