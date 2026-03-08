#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f82_first_exported_preobserver_source_object_second_clause_typing_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f81 = load_json(
        "fundamental_action_reconstruction/generated/f81_first_additive_preobserver_strict_core_source_object_export_packet_summary.json"
    )
    state = f81["state"]["closed_form"]

    summary = {
        "stage": "F82",
        "lane": "second_admissibility_clause_typing_only",
        "goal": "freeze_one_minimal_typed_carrier_structure_sufficient_to_make_later_E_orient_export_meaningful",
        "status": "F82_EXECUTED_FIRST_EXPORTED_PREOBSERVER_SOURCE_OBJECT_SECOND_CLAUSE_TYPING_PACKET_NO_FALSE_PASS",
        "exported_object": "S_preLM_strict_core_source_object_v1",
        "carrier": {
            "typed_sum": ["V_topo", "L_int", "M_int"],
            "basis": ["u_T", "u_L", "u_M"],
            "topological_seed_slot": "u_T",
            "light_transport_slot": "u_L",
            "matter_encoding_slot": "u_M",
        },
        "state_support": {
            "coeff_u_T": state[0],
            "coeff_u_L": state[1],
            "coeff_u_M": state[2],
            "nonzero_light_support": state[1] != 0.0,
            "nonzero_matter_support": state[2] != 0.0,
        },
        "future_orientation_export_target_frame": ["u_T", "u_L"],
        "E_orient_exported": False,
        "full_admissibility_pass": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
