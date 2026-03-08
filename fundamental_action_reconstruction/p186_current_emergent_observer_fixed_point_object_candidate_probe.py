#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p186_current_emergent_observer_fixed_point_object_candidate_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n205 = load_json(
        "fundamental_action_reconstruction/generated/n205_current_emergent_observer_fixed_point_reduction_theorem_summary.json"
    )
    f98 = load_json(
        "fundamental_action_reconstruction/generated/f98_first_emergent_observer_fixed_point_object_candidate_packet_summary.json"
    )

    props = f98["observer_fixed_point_object_properties"]
    response = f98["source_observer_fixed_point_object"]
    checks_spec = [
        {
            "id": "observer_fixed_point_reduction_is_admissible",
            "actual": n205["theorem_result"]["admissible_K_obs_fixed_point"],
            "expected": True,
        },
        {
            "id": "derived_only_from_fixed_point_state",
            "actual": props["derived_only_from_fixed_point_state"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": props["strict_core_only"],
            "expected": True,
        },
        {
            "id": "downstream_fixed_point_object_only",
            "actual": props["downstream_fixed_point_object_only"]
            and not props["actual_emergent_observer_constructed"],
            "expected": True,
        },
        {
            "id": "one_dimensional_fixed_point_object_exported",
            "actual": response["one_dimensional_fixed_point_object"],
            "expected": True,
        },
        {
            "id": "positive_fixed_point_object_amplitude",
            "actual": response["positive_fixed_point_object_amplitude"],
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
            "stage": "P186",
            "lane": "current_emergent_observer_fixed_point_object_candidate_only",
            "status": "P186_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_EMERGENT_OBSERVER_FIXED_POINT_OBJECT_CANDIDATE_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P186",
            "lane": "current_emergent_observer_fixed_point_object_candidate_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_EMERGENT_OBSERVER_FIXED_POINT_OBJECT_CANDIDATE_FROM_K_OBS_FIXED_POINT_PRELM_V1_AFTER_P186",
            "source_object": "S_preLM_strict_core_source_object_v1",
            "observer_fixed_point_reduction_operator": "K_obs_fixed_point_preLM_v1",
            "observer_fixed_point_object_map": "M_obs_fixed_object_preLM_v1",
            "checks": checks,
            "fixed_point_object_candidate_exported": True,
            "actual_emergent_observer_constructed": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
