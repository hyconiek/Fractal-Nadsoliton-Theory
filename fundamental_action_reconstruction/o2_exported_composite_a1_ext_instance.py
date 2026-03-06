#!/usr/bin/env python3
import json
from pathlib import Path

root = Path(__file__).resolve().parent
out = root / "generated"
out.mkdir(exist_ok=True)

instance = {
    "stage": "O2",
    "object_id": "A_1_ext",
    "mode": "exported_composite_A_1",
    "pair": "pair1",
    "carrier": "V_1 = span{c_1, s_1}",
    "basis_order": ["c_1", "s_1"],
    "matrix_form": [["a_1", "b_1"], ["b_1", "d_1"]],
    "coefficient_semantics": {
        "a_1": "<c_1, A_1_ext c_1>",
        "b_1": "<c_1, A_1_ext s_1>",
        "d_1": "<s_1, A_1_ext s_1>",
    },
    "provenance": {
        "lane": "hypothesis_extension_only",
        "base_kernel_contains_obs": False,
        "strict_core_reinterpretation": False,
        "selector_smuggling": False,
        "operator_origin": "exported_composite_A_1",
        "current_composite_export_witness": "A_1_cand",
    },
    "coefficient_state": {
        "a_1": "UNRESOLVED_VALUE",
        "b_1": "UNRESOLVED_VALUE",
        "d_1": "UNRESOLVED_VALUE",
    },
    "computable_now": False,
    "no_false_pass": True,
}

summary = {
    "stage": "O2",
    "status": "PASS_PARTIAL_EXPORTED_COMPOSITE_A1_EXT_INSTANCE_CREATED_VALUES_UNRESOLVED",
    "frontier": [
        "O2_B1",
        "H28_B1",
        "T12_B1",
        "T2_B1",
        "C32_B2",
    ],
    "theorem_level_pass": False,
    "full_closure_pass": False,
}

(out / "o2_exported_composite_a1_ext_instance.json").write_text(
    json.dumps(instance, indent=2) + "\n", encoding="utf-8"
)
(out / "o2_exported_composite_a1_ext_instance_summary.json").write_text(
    json.dumps(summary, indent=2) + "\n", encoding="utf-8"
)
