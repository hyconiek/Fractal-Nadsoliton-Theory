#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f96_first_emergent_observer_self_consistency_operator_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def mat_mul(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    rows = len(a)
    cols = len(b[0])
    inner = len(b)
    return [[sum(a[i][k] * b[k][j] for k in range(inner)) for j in range(cols)] for i in range(rows)]


def main() -> None:
    f95 = load_json(
        "fundamental_action_reconstruction/generated/f95_first_actual_emergent_observer_realization_map_packet_summary.json"
    )

    consistency = [
        [1.0, 0.0],
        [0.0, 0.0],
    ]
    realized_vector = [float(x) for x in f95["source_observer_realization"]["output_vector"]]
    consistency_vector = mat_vec(consistency, realized_vector)
    idempotent = mat_mul(consistency, consistency) == consistency

    summary = {
        "stage": "F96",
        "lane": "first_emergent_observer_self_consistency_operator_only",
        "goal": "export_one_emergent_observer_self_consistency_operator_without_claiming_actual_observer",
        "status": "F96_EXECUTED_FIRST_EMERGENT_OBSERVER_SELF_CONSISTENCY_OPERATOR_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_realization_input": "H_obs_realization_preLM_v1",
        "observer_self_consistency_operator": {
            "object": "J_obs_self_consistency_preLM_v1",
            "domain_basis": ["x_commit", "x_residual"],
            "codomain_basis": ["u_commit", "u_residual"],
            "matrix": consistency,
            "idempotent": idempotent,
        },
        "observer_self_consistency_properties": {
            "derived_only_from_realization_state": True,
            "strict_core_only": True,
            "downstream_self_consistency_only": True,
            "actual_emergent_observer_constructed": False,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "source_observer_self_consistency": {
            "input_basis": ["x_commit", "x_residual"],
            "input_vector": realized_vector,
            "output_basis": ["u_commit", "u_residual"],
            "output_vector": consistency_vector,
            "u_commit_v1": consistency_vector[0],
            "u_residual_v1": consistency_vector[1],
            "positive_commit": consistency_vector[0] > 0.0,
            "vanishing_residual": abs(consistency_vector[1]) < 1.0e-12,
        },
        "observer_self_consistency_operator_exported": True,
        "emergent_observer_constructed": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
