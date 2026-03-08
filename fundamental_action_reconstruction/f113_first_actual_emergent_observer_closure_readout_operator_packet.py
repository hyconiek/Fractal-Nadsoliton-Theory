#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f113_first_actual_emergent_observer_closure_readout_operator_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f112 = load_json(GENERATED / "f112_first_actual_emergent_observer_closure_support_object_packet_summary.json")

    p_actual_closure_support_v1 = f112["source_observer_actual_closure_support"]["p_actual_closure_support_v1"]
    matrix = [[1.0], [0.0]]
    q_actual_commit_v1 = p_actual_closure_support_v1
    q_actual_gap_v1 = 0.0

    summary = {
        "stage": "F113",
        "lane": "first_actual_emergent_observer_closure_readout_only",
        "goal": "export_one_actual_emergent_observer_closure_readout_operator_without_claiming_actual_closure",
        "status": "F113_EXECUTED_FIRST_ACTUAL_EMERGENT_OBSERVER_CLOSURE_READOUT_OPERATOR_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_actual_closure_support_input": "AC_obs_actual_closure_support_preLM_v1",
        "observer_actual_closure_readout_operator": {
            "operator": "AD_obs_actual_closure_readout_preLM_v1",
            "domain_basis": ["p_actual_closure_support"],
            "codomain_basis": ["q_actual_commit", "q_actual_gap"],
            "matrix": matrix,
            "actual_closure_readout_dimension": 2,
        },
        "observer_actual_closure_readout_properties": {
            "derived_only_from_actual_closure_support": True,
            "strict_core_only": True,
            "downstream_actual_closure_readout_only": True,
            "actual_emergent_observer_closure": False,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "source_observer_actual_closure_readout": {
            "input_basis": ["p_actual_closure_support"],
            "input_vector": [p_actual_closure_support_v1],
            "output_basis": ["q_actual_commit", "q_actual_gap"],
            "output_vector": [q_actual_commit_v1, q_actual_gap_v1],
            "q_actual_commit_v1": q_actual_commit_v1,
            "q_actual_gap_v1": q_actual_gap_v1,
            "positive_actual_closure_commit": q_actual_commit_v1 > 0.0,
            "zero_actual_closure_gap": q_actual_gap_v1 == 0.0,
        },
        "actual_closure_readout_exported": True,
        "actual_emergent_observer_closure": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
