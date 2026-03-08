#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f100_first_emergent_observer_closure_realization_map_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f99 = load_json(
        "fundamental_action_reconstruction/generated/f99_first_emergent_observer_closure_candidate_packet_summary.json"
    )

    realization_map = [[1.0]]
    closure_candidate_vector = [float(x) for x in f99["source_observer_closure_candidate"]["output_vector"]]
    closure_real_vector = mat_vec(realization_map, closure_candidate_vector)

    summary = {
        "stage": "F100",
        "lane": "first_emergent_observer_closure_realization_only",
        "goal": "export_one_emergent_observer_closure_realization_without_claiming_actual_closure",
        "status": "F100_EXECUTED_FIRST_EMERGENT_OBSERVER_CLOSURE_REALIZATION_MAP_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_closure_candidate_input": "N_obs_closure_candidate_preLM_v1",
        "observer_closure_realization_map": {
            "object": "Q_obs_closure_realization_preLM_v1",
            "domain_basis": ["c_closure"],
            "codomain_basis": ["d_closure"],
            "matrix": realization_map,
            "closure_realization_dimension": 1,
        },
        "observer_closure_realization_properties": {
            "derived_only_from_closure_candidate_state": True,
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
            "input_basis": ["c_closure"],
            "input_vector": closure_candidate_vector,
            "output_basis": ["d_closure"],
            "output_vector": closure_real_vector,
            "d_closure_v1": closure_real_vector[0],
            "positive_closure_realization_amplitude": closure_real_vector[0] > 0.0,
            "one_dimensional_closure_realization": len(closure_real_vector) == 1,
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
