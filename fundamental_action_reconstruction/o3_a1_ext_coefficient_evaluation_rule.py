#!/usr/bin/env python3
import json
from pathlib import Path

root = Path(__file__).resolve().parent
out = root / "generated"
out.mkdir(exist_ok=True)

packet = {
    "stage": "O3",
    "object_id": "A_1_ext",
    "pair": "pair1",
    "basis_order": ["c_1", "s_1"],
    "carrier": "V_1 = span{c_1, s_1}",
    "input_matrix_form": [["a_1", "b_1"], ["b_1", "d_1"]],
    "evaluation_rule": {
        "a_1": "read_entry(A_1_ext,0,0)",
        "b_1": "read_entry(A_1_ext,0,1)",
        "d_1": "read_entry(A_1_ext,1,1)",
        "trace_A_1": "read_entry(A_1_ext,0,0) + read_entry(A_1_ext,1,1)",
        "Delta_1": [
            "read_entry(A_1_ext,0,0) - read_entry(A_1_ext,1,1)",
            "read_entry(A_1_ext,0,1)"
        ]
    },
    "current_entry_state": {
        "a_1": "SYMBOLIC_PLACEHOLDER",
        "b_1": "SYMBOLIC_PLACEHOLDER",
        "d_1": "SYMBOLIC_PLACEHOLDER"
    },
    "lane": "hypothesis_extension_only",
    "computes_values_now": False,
    "no_false_pass": True
}

summary = {
    "stage": "O3",
    "status": "PASS_PARTIAL_A1_EXT_COEFFICIENT_EVALUATION_RULE_DEFINED_VALUES_STILL_SYMBOLIC",
    "frontier": [
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

(out / "o3_a1_ext_coefficient_evaluation_rule.json").write_text(
    json.dumps(packet, indent=2) + "\n", encoding="utf-8"
)
(out / "o3_a1_ext_coefficient_evaluation_rule_summary.json").write_text(
    json.dumps(summary, indent=2) + "\n", encoding="utf-8"
)
