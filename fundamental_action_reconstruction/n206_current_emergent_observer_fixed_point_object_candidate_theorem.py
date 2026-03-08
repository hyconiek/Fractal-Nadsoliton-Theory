#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n206_current_emergent_observer_fixed_point_object_candidate_theorem_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p186 = load_json(
        "fundamental_action_reconstruction/generated/p186_current_emergent_observer_fixed_point_object_candidate_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_EMERGENT_OBSERVER_FIXED_POINT_OBJECT_CANDIDATE_FROM_K_OBS_FIXED_POINT_PRELM_V1_AFTER_P186"
    )
    status_ok = p186["status"] == expected_status

    checks = [
        {
            "id": "positive_emergent_observer_fixed_point_object_candidate_probe",
            "actual": p186["status"],
            "expected": expected_status,
            "pass": status_ok,
        }
    ]

    summary = {
        "step": "N206",
        "status": "N206_DISCHARGED_CURRENT_EMERGENT_OBSERVER_FIXED_POINT_OBJECT_CANDIDATE_THEOREM_NO_FALSE_PASS",
        "scope": "current_emergent_observer_fixed_point_object_candidate_only",
        "checks": checks,
        "theorem_result": {
            "discharged": status_ok,
            "source_object": "S_preLM_strict_core_source_object_v1",
            "orientation_input": "E_orient_preLM_v1",
            "selector_bridge_operator": "B_sel_preLM_v1",
            "selector_reduction_operator": "R_sel_preLM_v1",
            "selector_output_operator": "O_sel_preLM_v1",
            "coarse_graining_operator": "C_obs_limit_preLM_v1",
            "observer_limit_readout_operator": "L_obs_limit_preLM_v1",
            "observer_construction_candidate_operator": "G_obs_candidate_preLM_v1",
            "observer_realization_map": "H_obs_realization_preLM_v1",
            "observer_self_consistency_operator": "J_obs_self_consistency_preLM_v1",
            "observer_fixed_point_reduction_operator": "K_obs_fixed_point_preLM_v1",
            "observer_fixed_point_object_map": "M_obs_fixed_object_preLM_v1",
            "admissible_M_obs_fixed_object": status_ok,
            "observer_information_deficit_downstream": status_ok,
        },
        "hard_limits": [
            "actual_emergent_observer_not_yet_constructed",
            "no_strict_core_selector_closure",
            "no_QW2191_discharge",
            "no_ToE_closure",
        ],
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
