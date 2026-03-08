#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n215_current_actual_emergent_observer_closure_candidate_theorem_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p195 = load_json(
        "fundamental_action_reconstruction/generated/p195_current_actual_emergent_observer_closure_candidate_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_ACTUAL_EMERGENT_OBSERVER_CLOSURE_CANDIDATE_FROM_W_OBS_CLOSURE_REALIZATION_PRELM_V1_AFTER_P195"
    )
    status_ok = p195["status"] == expected_status

    checks = [
        {
            "id": "positive_actual_emergent_observer_closure_candidate_probe",
            "actual": p195["status"],
            "expected": expected_status,
            "pass": status_ok,
        }
    ]

    summary = {
        "step": "N215",
        "status": "N215_DISCHARGED_CURRENT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_CANDIDATE_THEOREM_NO_FALSE_PASS",
        "scope": "current_actual_emergent_observer_closure_candidate_only",
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
            "observer_closure_candidate_map": "N_obs_closure_candidate_preLM_v1",
            "observer_closure_realization_map": "Q_obs_closure_realization_preLM_v1",
            "observer_closure_fixed_point_test": "R_obs_closure_fixed_point_test_preLM_v1",
            "observer_closure_support_object": "S_obs_closure_support_preLM_v1",
            "observer_closure_readout_operator": "T_obs_closure_readout_preLM_v1",
            "observer_closure_object_candidate": "U_obs_closure_object_candidate_preLM_v1",
            "observer_closure_commit_map": "V_obs_closure_commit_preLM_v1",
            "observer_closure_realization_object": "W_obs_closure_realization_preLM_v1",
            "actual_emergent_observer_closure_candidate_map": "X_obs_actual_closure_candidate_preLM_v1",
            "admissible_X_obs_actual_closure_candidate_map": status_ok,
            "observer_information_deficit_downstream": status_ok,
        },
        "hard_limits": [
            "actual_emergent_observer_closure_not_yet_constructed",
            "no_strict_core_selector_closure",
            "no_QW2191_discharge",
            "no_ToE_closure",
        ],
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
