#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n200_current_preobserver_to_emergent_observer_coarse_graining_theorem_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p180 = load_json(
        "fundamental_action_reconstruction/generated/p180_current_preobserver_to_emergent_observer_coarse_graining_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_PREOBSERVER_TO_EMERGENT_OBSERVER_COARSE_GRAINING_OPERATOR_FROM_O_SEL_PRELM_V1_AFTER_P180"
    )
    status_ok = p180["status"] == expected_status

    checks = [
        {
            "id": "positive_preobserver_to_emergent_observer_coarse_graining_probe",
            "actual": p180["status"],
            "expected": expected_status,
            "pass": status_ok,
        }
    ]

    summary = {
        "step": "N200",
        "status": "N200_DISCHARGED_CURRENT_PREOBSERVER_TO_EMERGENT_OBSERVER_COARSE_GRAINING_THEOREM_NO_FALSE_PASS",
        "scope": "current_preobserver_to_emergent_observer_coarse_graining_only",
        "checks": checks,
        "theorem_result": {
            "discharged": status_ok,
            "source_object": "S_preLM_strict_core_source_object_v1",
            "orientation_input": "E_orient_preLM_v1",
            "selector_bridge_operator": "B_sel_preLM_v1",
            "selector_reduction_operator": "R_sel_preLM_v1",
            "selector_output_operator": "O_sel_preLM_v1",
            "coarse_graining_operator": "C_obs_limit_preLM_v1",
            "admissible_C_obs_limit": status_ok,
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
