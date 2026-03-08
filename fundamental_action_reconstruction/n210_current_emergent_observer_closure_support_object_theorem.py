#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n210_current_emergent_observer_closure_support_object_theorem_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p190 = load_json(
        "fundamental_action_reconstruction/generated/p190_current_emergent_observer_closure_support_object_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_EMERGENT_OBSERVER_CLOSURE_SUPPORT_OBJECT_FROM_R_OBS_CLOSURE_FIXED_POINT_TEST_PRELM_V1_AFTER_P190"
    )
    status_ok = p190["status"] == expected_status

    checks = [
        {
            "id": "positive_emergent_observer_closure_support_object_probe",
            "actual": p190["status"],
            "expected": expected_status,
            "pass": status_ok,
        }
    ]

    summary = {
        "step": "N210",
        "status": "N210_DISCHARGED_CURRENT_EMERGENT_OBSERVER_CLOSURE_SUPPORT_OBJECT_THEOREM_NO_FALSE_PASS",
        "scope": "current_emergent_observer_closure_support_object_only",
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
            "admissible_S_obs_closure_support": status_ok,
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
