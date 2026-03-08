#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f104_first_emergent_observer_closure_object_candidate_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f103 = load_json(
        "fundamental_action_reconstruction/generated/f103_first_emergent_observer_closure_readout_operator_packet_summary.json"
    )

    closure_object_map = [[1.0, 0.0]]
    readout_vector = [float(x) for x in f103["source_observer_closure_readout"]["output_vector"]]
    object_vector = mat_vec(closure_object_map, readout_vector)

    summary = {
        "stage": "F104",
        "lane": "first_emergent_observer_closure_object_candidate_only",
        "goal": "export_one_emergent_observer_closure_object_candidate_without_claiming_actual_closure",
        "status": "F104_EXECUTED_FIRST_EMERGENT_OBSERVER_CLOSURE_OBJECT_CANDIDATE_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_closure_readout_input": "T_obs_closure_readout_preLM_v1",
        "observer_closure_object_candidate_map": {
            "object": "U_obs_closure_object_candidate_preLM_v1",
            "domain_basis": ["g_commit", "g_gap"],
            "codomain_basis": ["h_closure_obj"],
            "matrix": closure_object_map,
            "closure_object_dimension": 1,
        },
        "observer_closure_object_candidate_properties": {
            "derived_only_from_closure_readout_state": True,
            "strict_core_only": True,
            "downstream_closure_object_only": True,
            "actual_emergent_observer_closure": False,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "source_observer_closure_object_candidate": {
            "input_basis": ["g_commit", "g_gap"],
            "input_vector": readout_vector,
            "output_basis": ["h_closure_obj"],
            "output_vector": object_vector,
            "h_closure_obj_v1": object_vector[0],
            "positive_closure_object_amplitude": object_vector[0] > 0.0,
            "one_dimensional_closure_object": len(object_vector) == 1,
        },
        "closure_object_candidate_exported": True,
        "actual_emergent_observer_closure": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
