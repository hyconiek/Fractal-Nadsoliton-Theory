#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p187_current_emergent_observer_closure_candidate_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n206 = load_json(
        "fundamental_action_reconstruction/generated/n206_current_emergent_observer_fixed_point_object_candidate_theorem_summary.json"
    )
    f99 = load_json(
        "fundamental_action_reconstruction/generated/f99_first_emergent_observer_closure_candidate_packet_summary.json"
    )

    props = f99["observer_closure_candidate_properties"]
    response = f99["source_observer_closure_candidate"]
    checks_spec = [
        {
            "id": "observer_fixed_point_object_is_admissible",
            "actual": n206["theorem_result"]["admissible_M_obs_fixed_object"],
            "expected": True,
        },
        {
            "id": "derived_only_from_fixed_point_object_state",
            "actual": props["derived_only_from_fixed_point_object_state"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": props["strict_core_only"],
            "expected": True,
        },
        {
            "id": "downstream_closure_candidate_only",
            "actual": props["downstream_closure_candidate_only"]
            and not props["actual_emergent_observer_closure"],
            "expected": True,
        },
        {
            "id": "one_dimensional_closure_candidate_exported",
            "actual": response["one_dimensional_closure_candidate"],
            "expected": True,
        },
        {
            "id": "positive_closure_candidate_amplitude",
            "actual": response["positive_closure_candidate_amplitude"],
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
            "stage": "P187",
            "lane": "current_emergent_observer_closure_candidate_only",
            "status": "P187_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_EMERGENT_OBSERVER_CLOSURE_CANDIDATE_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P187",
            "lane": "current_emergent_observer_closure_candidate_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_EMERGENT_OBSERVER_CLOSURE_CANDIDATE_FROM_M_OBS_FIXED_OBJECT_PRELM_V1_AFTER_P187",
            "source_object": "S_preLM_strict_core_source_object_v1",
            "observer_fixed_point_object_map": "M_obs_fixed_object_preLM_v1",
            "observer_closure_candidate_map": "N_obs_closure_candidate_preLM_v1",
            "checks": checks,
            "closure_candidate_exported": True,
            "actual_emergent_observer_closure": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
