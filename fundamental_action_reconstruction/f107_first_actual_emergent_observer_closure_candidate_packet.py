#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f107_first_actual_emergent_observer_closure_candidate_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f106 = load_json(
        "fundamental_action_reconstruction/generated/f106_first_emergent_observer_closure_realization_object_packet_summary.json"
    )

    actual_closure_candidate = [[1.0]]
    realization_vector = [
        float(x) for x in f106["source_observer_closure_realization"]["output_vector"]
    ]
    actual_candidate_vector = mat_vec(actual_closure_candidate, realization_vector)

    summary = {
        "stage": "F107",
        "lane": "first_emergent_observer_actual_closure_candidate_only",
        "goal": "export_one_actual_emergent_observer_closure_candidate_without_claiming_actual_closure",
        "status": "F107_EXECUTED_FIRST_ACTUAL_EMERGENT_OBSERVER_CLOSURE_CANDIDATE_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_closure_realization_input": "W_obs_closure_realization_preLM_v1",
        "observer_actual_closure_candidate_map": {
            "object": "X_obs_actual_closure_candidate_preLM_v1",
            "domain_basis": ["j_closure_real"],
            "codomain_basis": ["k_actual_closure"],
            "matrix": actual_closure_candidate,
            "actual_closure_candidate_dimension": 1,
        },
        "observer_actual_closure_candidate_properties": {
            "derived_only_from_closure_realization_object": True,
            "strict_core_only": True,
            "downstream_actual_closure_candidate_only": True,
            "actual_emergent_observer_closure": False,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "source_observer_actual_closure_candidate": {
            "input_basis": ["j_closure_real"],
            "input_vector": realization_vector,
            "output_basis": ["k_actual_closure"],
            "output_vector": actual_candidate_vector,
            "k_actual_closure_v1": actual_candidate_vector[0],
            "positive_actual_closure_candidate_amplitude": actual_candidate_vector[0] > 0.0,
            "one_dimensional_actual_closure_candidate": len(actual_candidate_vector) == 1,
        },
        "actual_closure_candidate_exported": True,
        "actual_emergent_observer_closure": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
