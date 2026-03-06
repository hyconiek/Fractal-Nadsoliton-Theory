#!/usr/bin/env python3
import json
from pathlib import Path

root = Path(__file__).resolve().parent
out = root / "generated"
out.mkdir(exist_ok=True)

packet = {
    "stage": "O4",
    "object_id": "A_1_ext",
    "pair": "pair1",
    "basis_order": ["c_1", "s_1"],
    "population_routes": {
        "Route_P1": {
            "kind": "direct_matrix_entry_export",
            "required_fields": [
                "object_id=A_1_ext",
                "pair=pair1",
                "basis_order=[c_1,s_1]",
                "entry=(i,j)",
                "value=v_x",
                "operator_origin=exported_composite_A_1",
                "lane=hypothesis_extension_only",
                "base_kernel_contains_obs=false",
                "strict_core_reinterpretation=false",
                "selector_smuggling=false"
            ]
        },
        "Route_P2": {
            "kind": "basis_scalar_export",
            "required_fields": [
                "x=<basis_i,A_1_ext basis_j>",
                "value=v_x",
                "operator_origin=exported_composite_A_1",
                "lane=hypothesis_extension_only",
                "base_kernel_contains_obs=false",
                "strict_core_reinterpretation=false",
                "selector_smuggling=false"
            ]
        }
    },
    "slot_map": {
        "a_1": {
            "entry": [0, 0],
            "scalar_form": "<c_1, A_1_ext c_1>"
        },
        "b_1": {
            "entry": [0, 1],
            "scalar_form": "<c_1, A_1_ext s_1>"
        },
        "d_1": {
            "entry": [1, 1],
            "scalar_form": "<s_1, A_1_ext s_1>"
        }
    },
    "actual_population_present_now": False,
    "no_false_pass": True
}

summary = {
    "stage": "O4",
    "status": "PASS_PARTIAL_A1_EXT_ENTRY_POPULATION_RULE_DEFINED_NO_ACTUAL_POPULATED_ENTRIES",
    "frontier": [
        "O4_B1",
        "O3_B1",
        "O2_B1",
        "H28_B1",
        "T12_B1",
        "T2_B1",
        "C32_B2"
    ],
    "theorem_level_pass": False,
    "full_closure_pass": False
}

(out / "o4_a1_ext_entry_population_rule.json").write_text(
    json.dumps(packet, indent=2) + "\n", encoding="utf-8"
)
(out / "o4_a1_ext_entry_population_rule_summary.json").write_text(
    json.dumps(summary, indent=2) + "\n", encoding="utf-8"
)
