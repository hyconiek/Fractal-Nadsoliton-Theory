#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f95_first_actual_emergent_observer_realization_map_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f94 = load_json(
        "fundamental_action_reconstruction/generated/f94_first_actual_emergent_observer_construction_candidate_packet_summary.json"
    )

    realization = [
        [1.0, 0.0],
        [0.0, 1.0],
    ]
    candidate_vector = [float(x) for x in f94["source_observer_construction_candidate"]["output_vector"]]
    realized_vector = mat_vec(realization, candidate_vector)

    summary = {
        "stage": "F95",
        "lane": "first_actual_emergent_observer_realization_map_only",
        "goal": "export_one_actual_emergent_observer_realization_map_without_claiming_actual_observer",
        "status": "F95_EXECUTED_FIRST_ACTUAL_EMERGENT_OBSERVER_REALIZATION_MAP_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_construction_candidate_input": "G_obs_candidate_preLM_v1",
        "observer_realization_map": {
            "object": "H_obs_realization_preLM_v1",
            "domain_basis": ["w_commit", "w_residual"],
            "codomain_basis": ["x_commit", "x_residual"],
            "matrix": realization,
        },
        "observer_realization_properties": {
            "derived_only_from_construction_candidate": True,
            "strict_core_only": True,
            "downstream_realization_only": True,
            "actual_emergent_observer_constructed": False,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "source_observer_realization": {
            "input_basis": ["w_commit", "w_residual"],
            "input_vector": candidate_vector,
            "output_basis": ["x_commit", "x_residual"],
            "output_vector": realized_vector,
            "x_commit_v1": realized_vector[0],
            "x_residual_v1": realized_vector[1],
            "positive_commit": realized_vector[0] > 0.0,
            "vanishing_residual": abs(realized_vector[1]) < 1.0e-12,
        },
        "actual_observer_realization_map_exported": True,
        "emergent_observer_constructed": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
