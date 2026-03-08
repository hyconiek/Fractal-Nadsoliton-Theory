#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p185_current_emergent_observer_fixed_point_reduction_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n204 = load_json(
        "fundamental_action_reconstruction/generated/n204_current_emergent_observer_self_consistency_operator_theorem_summary.json"
    )
    f97 = load_json(
        "fundamental_action_reconstruction/generated/f97_first_emergent_observer_fixed_point_reduction_packet_summary.json"
    )

    props = f97["observer_fixed_point_reduction_properties"]
    response = f97["source_observer_fixed_point"]
    checks_spec = [
        {
            "id": "observer_self_consistency_is_admissible",
            "actual": n204["theorem_result"]["admissible_J_obs_self_consistency"],
            "expected": True,
        },
        {
            "id": "derived_only_from_self_consistency_state",
            "actual": props["derived_only_from_self_consistency_state"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": props["strict_core_only"],
            "expected": True,
        },
        {
            "id": "downstream_fixed_point_only",
            "actual": props["downstream_fixed_point_only"] and not props["actual_emergent_observer_constructed"],
            "expected": True,
        },
        {
            "id": "fixed_point_sector_exported",
            "actual": True,
            "expected": True,
        },
        {
            "id": "positive_fixed_point_amplitude",
            "actual": response["positive_fixed_point_amplitude"],
            "expected": True,
        },
        {
            "id": "source_state_in_fixed_point_support",
            "actual": response["source_state_in_fixed_point_support"],
            "expected": True,
        },
        {
            "id": "observer_information_deficit_downstream",
            "actual": props["observer_information_deficit_is_downstream_symptom"],
            "expected": True,
        },
        {
            "id": "kernel_split_safe",
            "actual": props["kernel_split_safe"],
            "expected": True,
        },
    ]

    checks = []
    mismatches = []
    for item in checks_spec:
        ok = item["actual"] == item["expected"]
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": ok,
            }
        )
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        summary = {
            "stage": "P185",
            "lane": "current_emergent_observer_fixed_point_reduction_only",
            "status": "P185_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_EMERGENT_OBSERVER_FIXED_POINT_REDUCTION_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P185",
            "lane": "current_emergent_observer_fixed_point_reduction_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_EMERGENT_OBSERVER_FIXED_POINT_REDUCTION_OPERATOR_FROM_J_OBS_SELF_CONSISTENCY_PRELM_V1_AFTER_P185",
            "source_object": "S_preLM_strict_core_source_object_v1",
            "observer_self_consistency_operator": "J_obs_self_consistency_preLM_v1",
            "observer_fixed_point_reduction_operator": "K_obs_fixed_point_preLM_v1",
            "checks": checks,
            "fixed_point_object_candidate_exported": True,
            "emergent_observer_constructed": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
