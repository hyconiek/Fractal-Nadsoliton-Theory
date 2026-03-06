#!/usr/bin/env python3
import json
from pathlib import Path

root = Path(__file__).resolve().parent
out = root / "generated"
out.mkdir(exist_ok=True)

artifact = {
    "stage": "O1",
    "target_object": "A_1_ext",
    "carrier": "V_1 = span{c_1, s_1}",
    "lane": "hypothesis_extension_only",
    "admissible_modes": [
        "exported_composite_A_1",
        "P_1 E_1^* G_light^* R_mat^* O_readout R_mat G_light E_1 P_1",
    ],
    "legacy_alias": {"O_readout": "O_obs"},
    "computable_targets": ["a_1", "b_1", "d_1", "trace_A_1", "Delta_1"],
    "persisted_instance_present": False,
    "no_false_pass": True,
}

summary = {
    "stage": "O1",
    "status": "PASS_PARTIAL_EXPLICIT_A1_OPERATOR_DEFINITION_SPEC_WRITTEN_INSTANCE_ABSENT",
    "frontier": [
        "O1_B1",
        "H28_B1",
        "T12_B1",
        "T2_B1",
        "C32_B2",
    ],
    "theorem_level_pass": False,
    "full_closure_pass": False,
}

(out / "o1_explicit_a1_operator_definition_spec.json").write_text(
    json.dumps(artifact, indent=2) + "\n", encoding="utf-8"
)
(out / "o1_explicit_a1_operator_definition_spec_summary.json").write_text(
    json.dumps(summary, indent=2) + "\n", encoding="utf-8"
)
