#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f116_first_actual_emergent_observer_closure_realization_map_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f115 = load_json(GENERATED / "f115_first_actual_emergent_observer_closure_commit_map_packet_summary.json")

    actual_closure_realization = [[1.0, 0.0]]
    commit_vector = [
        float(x) for x in f115["source_observer_actual_closure_commit"]["output_vector"]
    ]
    realization_vector = mat_vec(actual_closure_realization, commit_vector)

    summary = {
        "stage": "F116",
        "lane": "first_actual_emergent_observer_closure_realization_only",
        "goal": "export_one_actual_emergent_observer_closure_realization_map_without_claiming_actual_closure",
        "status": "F116_EXECUTED_FIRST_ACTUAL_EMERGENT_OBSERVER_CLOSURE_REALIZATION_MAP_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_actual_closure_commit_input": "AF_obs_actual_closure_commit_preLM_v1",
        "observer_actual_closure_realization_operator": {
            "operator": "AG_obs_actual_closure_realization_preLM_v1",
            "domain_basis": ["s_actual_commit", "s_actual_residual"],
            "codomain_basis": ["t_actual_closure_real"],
            "matrix": actual_closure_realization,
            "actual_closure_realization_dimension": 1,
        },
        "observer_actual_closure_realization_properties": {
            "derived_only_from_actual_closure_commit_state": True,
            "strict_core_only": True,
            "downstream_actual_closure_realization_only": True,
            "actual_emergent_observer_closure": False,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "source_observer_actual_closure_realization": {
            "input_basis": ["s_actual_commit", "s_actual_residual"],
            "input_vector": commit_vector,
            "output_basis": ["t_actual_closure_real"],
            "output_vector": realization_vector,
            "t_actual_closure_real_v1": realization_vector[0],
            "positive_actual_closure_realization_amplitude": realization_vector[0] > 0.0,
            "actual_residual_annihilated": commit_vector[1] == 0.0,
        },
        "actual_closure_realization_exported": True,
        "actual_emergent_observer_closure": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
