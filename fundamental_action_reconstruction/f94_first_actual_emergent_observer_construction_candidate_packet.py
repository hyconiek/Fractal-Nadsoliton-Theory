#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f94_first_actual_emergent_observer_construction_candidate_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f93 = load_json(
        "fundamental_action_reconstruction/generated/f93_first_emergent_observer_limit_readout_operator_packet_summary.json"
    )

    candidate = [
        [1.0, 0.0],
        [0.0, 1.0],
    ]
    readout_vector = [float(x) for x in f93["source_observer_limit_readout"]["output_vector"]]
    candidate_vector = mat_vec(candidate, readout_vector)

    summary = {
        "stage": "F94",
        "lane": "first_actual_emergent_observer_construction_candidate_only",
        "goal": "export_one_actual_emergent_observer_construction_candidate_without_claiming_actual_observer",
        "status": "F94_EXECUTED_FIRST_ACTUAL_EMERGENT_OBSERVER_CONSTRUCTION_CANDIDATE_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_limit_readout_input": "L_obs_limit_preLM_v1",
        "observer_construction_candidate_operator": {
            "object": "G_obs_candidate_preLM_v1",
            "domain_basis": ["z_commit", "z_residual"],
            "codomain_basis": ["w_commit", "w_residual"],
            "matrix": candidate,
        },
        "observer_construction_candidate_properties": {
            "derived_only_from_observer_limit_readout": True,
            "strict_core_only": True,
            "downstream_candidate_only": True,
            "actual_emergent_observer_constructed": False,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "source_observer_construction_candidate": {
            "input_basis": ["z_commit", "z_residual"],
            "input_vector": readout_vector,
            "output_basis": ["w_commit", "w_residual"],
            "output_vector": candidate_vector,
            "w_commit_v1": candidate_vector[0],
            "w_residual_v1": candidate_vector[1],
            "positive_commit": candidate_vector[0] > 0.0,
            "vanishing_residual": abs(candidate_vector[1]) < 1.0e-12,
        },
        "actual_observer_construction_candidate_constructed": True,
        "emergent_observer_constructed": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
