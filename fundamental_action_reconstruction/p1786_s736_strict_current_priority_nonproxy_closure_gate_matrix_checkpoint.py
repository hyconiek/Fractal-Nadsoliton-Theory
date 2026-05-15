#!/usr/bin/env python3
"""P1786 S736 strict current-priority nonproxy closure gate matrix checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"


def load_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p1764 = load_json(GENERATED / "p1764_s714_strict_covariant_nonproxy_ea_eh_explicit_export_checkpoint.json")
    p1765 = load_json(GENERATED / "p1765_s715_strict_nonproxy_metric_elg_explicit_export_checkpoint.json")

    out = {
        "checkpoint_id": "P1786_S736",
        "title": "STRICT_CURRENT_PRIORITY_NONPROXY_CLOSURE_GATE_MATRIX",
        "strict_only_policy": {
            "legacy_bridge_used": False,
            "proxy_shortcut_used": False,
            "pass_requires_explicit_witness": True,
        },
        "input_reuse": [
            "p1764_s714_strict_covariant_nonproxy_ea_eh_explicit_export_checkpoint.json",
            "p1765_s715_strict_nonproxy_metric_elg_explicit_export_checkpoint.json",
        ],
        "gate_matrix": {
            "G1_EA_NONPROXY_EXPLICIT_EXPORT": "OPEN",
            "G2_EH_NONPROXY_EXPLICIT_EXPORT": "OPEN",
            "G3_ELG_NONPROXY_EXPLICIT_EXPORT": "OPEN",
            "G4_BOUNDARY_TERM_CONTROL": "PARTIAL_OPEN",
            "G5_H1_4D_WEAK_FORM": "OPEN",
            "G6_BIANCHI_WARD": "OPEN",
            "G7_BRST_NILPOTENCY": "OPEN",
            "G8_CUTKOSKY_UNITARITY": "OPEN",
        },
        "evidence_snapshot": {
            "p1764_status": p1764.get("status", "MISSING"),
            "p1765_status": p1765.get("status", "MISSING"),
            "p1764_readiness": p1764.get("readiness_update", {}),
            "p1765_readiness": p1765.get("readiness_update", {}),
        },
        "false_pass_risks": [
            "FULL_EXPORT_claim_without_componentwise_covariant_expansion",
            "PASS_claim_without_residual_or_witness",
            "LOCAL_or_REDUCED_state_promoted_to_GLOBAL_NONPROXY",
        ],
        "next_honest_step": "Run unified componentwise export for E_A^mu and E_H on shared background family, then execute H1 weak-form witness without PASS claim before explicit residual=0.",
        "status": "STRICT_PRIORITY_GATES_OPEN_KEEP_THEOREM_LEVEL_OPEN",
    }

    GENERATED.mkdir(parents=True, exist_ok=True)
    out_path = GENERATED / "p1786_s736_strict_current_priority_nonproxy_closure_gate_matrix_checkpoint.json"
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(out_path)


if __name__ == "__main__":
    main()
