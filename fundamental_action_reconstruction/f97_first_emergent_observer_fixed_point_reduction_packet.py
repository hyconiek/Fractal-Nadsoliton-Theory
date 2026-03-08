#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f97_first_emergent_observer_fixed_point_reduction_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f96 = load_json(
        "fundamental_action_reconstruction/generated/f96_first_emergent_observer_self_consistency_operator_packet_summary.json"
    )

    reduction = [
        [1.0, 0.0],
    ]
    self_consistency_vector = [float(x) for x in f96["source_observer_self_consistency"]["output_vector"]]
    fixed_point_vector = mat_vec(reduction, self_consistency_vector)

    summary = {
        "stage": "F97",
        "lane": "first_emergent_observer_fixed_point_reduction_only",
        "goal": "export_one_emergent_observer_fixed_point_reduction_without_claiming_actual_observer",
        "status": "F97_EXECUTED_FIRST_EMERGENT_OBSERVER_FIXED_POINT_REDUCTION_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_self_consistency_input": "J_obs_self_consistency_preLM_v1",
        "observer_fixed_point_reduction_operator": {
            "object": "K_obs_fixed_point_preLM_v1",
            "domain_basis": ["u_commit", "u_residual"],
            "codomain_basis": ["f_commit"],
            "matrix": reduction,
            "fixed_point_support_dimension": 1,
        },
        "observer_fixed_point_reduction_properties": {
            "derived_only_from_self_consistency_state": True,
            "strict_core_only": True,
            "downstream_fixed_point_only": True,
            "actual_emergent_observer_constructed": False,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "source_observer_fixed_point": {
            "input_basis": ["u_commit", "u_residual"],
            "input_vector": self_consistency_vector,
            "output_basis": ["f_commit"],
            "output_vector": fixed_point_vector,
            "f_commit_v1": fixed_point_vector[0],
            "positive_fixed_point_amplitude": fixed_point_vector[0] > 0.0,
            "source_state_in_fixed_point_support": abs(self_consistency_vector[1]) < 1.0e-12,
        },
        "fixed_point_object_candidate_exported": True,
        "emergent_observer_constructed": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
