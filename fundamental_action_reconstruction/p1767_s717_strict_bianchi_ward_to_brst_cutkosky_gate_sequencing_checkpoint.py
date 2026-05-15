#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1766 = GEN / "p1766_s716_strict_forward_reverse_state_vector_update_with_bianchi_ward_gate_checkpoint.json"
OUT = GEN / "p1767_s717_strict_bianchi_ward_to_brst_cutkosky_gate_sequencing_checkpoint.json"


def main() -> None:
    p1766 = json.loads(IN1766.read_text(encoding="utf-8"))

    sequencing = {
        "G_BW": {
            "name": "Bianchi_Ward_divergence_gate",
            "required_input": "componentwise_EL_g_minus_E_munu_B1_B2_B3_C1_C2",
            "admissible_outputs": ["PASS_ZERO", "OBSTRUCTION_WITH_DIVERGENCE_TRACE"],
            "status": "OPEN",
        },
        "G_BRST": {
            "name": "BRST_nilpotency_gate",
            "requires": ["G_BW:PASS_ZERO", "ghost_sector_nonproxy_export", "BV_BRST_operator_map"],
            "admissible_outputs": ["PASS_NILPOTENT", "OBSTRUCTION_WITH_NILPOTENCY_TRACE"],
            "status": "OPEN_BLOCKED_BY_G_BW",
        },
        "G_CUT": {
            "name": "Cutkosky_unitarity_gate",
            "requires": ["G_BW:PASS_ZERO", "G_BRST:PASS_NILPOTENT", "on_shell_cut_channel_pack"],
            "admissible_outputs": ["PASS_UNITARITY", "OBSTRUCTION_WITH_CUT_TRACE"],
            "status": "OPEN_BLOCKED_BY_G_BW_AND_G_BRST",
        },
    }

    payload = {
        "checkpoint": "P1767_S717",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_anchor": "p1766_s716",
        "state_vector_ref": p1766.get("updated_state_vector", {}),
        "gate_sequencing_contract": sequencing,
        "policy": {
            "no_false_pass": True,
            "promotion_rule": "QG theorem gates cannot be promoted if G_BW is not PASS_ZERO",
            "classification_required": ["LOCAL_vs_GLOBAL", "REDUCED_vs_NONPROXY", "SCAFFOLD_vs_FULL_EXPORT"],
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "next_honest_step": "Wykonać G_BW na tej samej rodzinie teł i wydać PASS_ZERO albo OBSTRUCTION_WITH_DIVERGENCE_TRACE; bez tego nie uruchamiać formalnie G_BRST ani G_CUT.",
        "lay_summary": "Ustalono twardą kolejność testów: najpierw zgodność Bianchi/Ward, potem BRST, a dopiero na końcu unitarność Cutkosky. Bez pierwszego kroku nie wolno ogłaszać sukcesu dalej.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
