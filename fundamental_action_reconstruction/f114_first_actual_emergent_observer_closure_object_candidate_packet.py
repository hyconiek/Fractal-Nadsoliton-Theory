#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f114_first_actual_emergent_observer_closure_object_candidate_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f113 = load_json(GENERATED / "f113_first_actual_emergent_observer_closure_readout_operator_packet_summary.json")

    q_actual_commit_v1 = f113["source_observer_actual_closure_readout"]["q_actual_commit_v1"]
    matrix = [[1.0, 0.0]]
    r_actual_closure_obj_v1 = q_actual_commit_v1

    summary = {
        "stage": "F114",
        "lane": "first_actual_emergent_observer_closure_object_candidate_only",
        "goal": "export_one_actual_emergent_observer_closure_object_candidate_without_claiming_actual_closure",
        "status": "F114_EXECUTED_FIRST_ACTUAL_EMERGENT_OBSERVER_CLOSURE_OBJECT_CANDIDATE_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_actual_closure_readout_input": "AD_obs_actual_closure_readout_preLM_v1",
        "observer_actual_closure_object_candidate_operator": {
            "operator": "AE_obs_actual_closure_object_candidate_preLM_v1",
            "domain_basis": ["q_actual_commit", "q_actual_gap"],
            "codomain_basis": ["r_actual_closure_obj"],
            "matrix": matrix,
            "actual_closure_object_candidate_dimension": 1,
        },
        "observer_actual_closure_object_candidate_properties": {
            "derived_only_from_actual_closure_readout": True,
            "strict_core_only": True,
            "downstream_actual_closure_object_candidate_only": True,
            "actual_emergent_observer_closure": False,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "source_observer_actual_closure_object_candidate": {
            "input_basis": ["q_actual_commit", "q_actual_gap"],
            "input_vector": [
                q_actual_commit_v1,
                f113["source_observer_actual_closure_readout"]["q_actual_gap_v1"],
            ],
            "output_basis": ["r_actual_closure_obj"],
            "output_vector": [r_actual_closure_obj_v1],
            "r_actual_closure_obj_v1": r_actual_closure_obj_v1,
            "positive_actual_closure_object_candidate_amplitude": r_actual_closure_obj_v1 > 0.0,
            "actual_gap_annihilated": True,
        },
        "actual_closure_object_candidate_exported": True,
        "actual_emergent_observer_closure": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
