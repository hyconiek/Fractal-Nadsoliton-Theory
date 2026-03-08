#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f103_first_emergent_observer_closure_readout_operator_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f102 = load_json(
        "fundamental_action_reconstruction/generated/f102_first_emergent_observer_closure_support_object_packet_summary.json"
    )

    closure_readout = [
        [1.0],
        [0.0],
    ]
    support_vector = [float(x) for x in f102["source_observer_closure_support"]["output_vector"]]
    readout_vector = mat_vec(closure_readout, support_vector)

    summary = {
        "stage": "F103",
        "lane": "first_emergent_observer_closure_readout_only",
        "goal": "export_one_emergent_observer_closure_readout_without_claiming_actual_closure",
        "status": "F103_EXECUTED_FIRST_EMERGENT_OBSERVER_CLOSURE_READOUT_OPERATOR_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_closure_support_input": "S_obs_closure_support_preLM_v1",
        "observer_closure_readout_operator": {
            "object": "T_obs_closure_readout_preLM_v1",
            "domain_basis": ["f_closure_support"],
            "codomain_basis": ["g_commit", "g_gap"],
            "matrix": closure_readout,
            "closure_readout_dimension": 2,
        },
        "observer_closure_readout_properties": {
            "derived_only_from_closure_support_state": True,
            "strict_core_only": True,
            "downstream_closure_readout_only": True,
            "actual_emergent_observer_closure": False,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "source_observer_closure_readout": {
            "input_basis": ["f_closure_support"],
            "input_vector": support_vector,
            "output_basis": ["g_commit", "g_gap"],
            "output_vector": readout_vector,
            "g_commit_v1": readout_vector[0],
            "g_gap_v1": readout_vector[1],
            "positive_commit_amplitude": readout_vector[0] > 0.0,
            "zero_gap_channel": abs(readout_vector[1]) < 1e-12,
        },
        "closure_readout_exported": True,
        "actual_emergent_observer_closure": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
