#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f106_first_emergent_observer_closure_realization_object_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f105 = load_json(
        "fundamental_action_reconstruction/generated/f105_first_emergent_observer_closure_commit_map_packet_summary.json"
    )

    closure_realization = [[1.0, 0.0]]
    commit_vector = [float(x) for x in f105["source_observer_closure_commit"]["output_vector"]]
    realization_vector = mat_vec(closure_realization, commit_vector)

    summary = {
        "stage": "F106",
        "lane": "first_emergent_observer_closure_realization_only",
        "goal": "export_one_emergent_observer_closure_realization_object_without_claiming_actual_closure",
        "status": "F106_EXECUTED_FIRST_EMERGENT_OBSERVER_CLOSURE_REALIZATION_OBJECT_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_closure_commit_input": "V_obs_closure_commit_preLM_v1",
        "observer_closure_realization_map": {
            "object": "W_obs_closure_realization_preLM_v1",
            "domain_basis": ["i_commit", "i_residual"],
            "codomain_basis": ["j_closure_real"],
            "matrix": closure_realization,
            "closure_realization_dimension": 1,
        },
        "observer_closure_realization_properties": {
            "derived_only_from_closure_commit_state": True,
            "strict_core_only": True,
            "downstream_closure_realization_only": True,
            "actual_emergent_observer_closure": False,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "source_observer_closure_realization": {
            "input_basis": ["i_commit", "i_residual"],
            "input_vector": commit_vector,
            "output_basis": ["j_closure_real"],
            "output_vector": realization_vector,
            "j_closure_real_v1": realization_vector[0],
            "positive_closure_realization_amplitude": realization_vector[0] > 0.0,
            "one_dimensional_closure_realization": len(realization_vector) == 1,
        },
        "closure_realization_exported": True,
        "actual_emergent_observer_closure": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
