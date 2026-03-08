#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p196_current_actual_emergent_observer_closure_object_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n215 = load_json(
        "fundamental_action_reconstruction/generated/n215_current_actual_emergent_observer_closure_candidate_theorem_summary.json"
    )
    f108 = load_json(
        "fundamental_action_reconstruction/generated/f108_first_actual_emergent_observer_closure_object_packet_summary.json"
    )

    props = f108["observer_actual_closure_object_properties"]
    response = f108["source_observer_actual_closure_object"]
    checks_spec = [
        {
            "id": "observer_actual_closure_candidate_is_admissible",
            "actual": n215["theorem_result"]["admissible_X_obs_actual_closure_candidate_map"],
            "expected": True,
        },
        {
            "id": "derived_only_from_actual_closure_candidate_state",
            "actual": props["derived_only_from_actual_closure_candidate_state"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": props["strict_core_only"],
            "expected": True,
        },
        {
            "id": "downstream_actual_closure_object_only",
            "actual": props["downstream_actual_closure_object_only"]
            and not props["actual_emergent_observer_closure"],
            "expected": True,
        },
        {
            "id": "positive_actual_closure_object_amplitude",
            "actual": response["positive_actual_closure_object_amplitude"],
            "expected": True,
        },
        {
            "id": "one_dimensional_actual_closure_object",
            "actual": response["one_dimensional_actual_closure_object"],
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
            "stage": "P196",
            "lane": "current_actual_emergent_observer_closure_object_only",
            "status": "P196_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_OBJECT_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P196",
            "lane": "current_actual_emergent_observer_closure_object_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_ACTUAL_EMERGENT_OBSERVER_CLOSURE_OBJECT_FROM_X_OBS_ACTUAL_CLOSURE_CANDIDATE_PRELM_V1_AFTER_P196",
            "source_object": "S_preLM_strict_core_source_object_v1",
            "observer_actual_closure_candidate_map": "X_obs_actual_closure_candidate_preLM_v1",
            "observer_actual_closure_object_map": "Y_obs_actual_closure_object_preLM_v1",
            "checks": checks,
            "actual_closure_object_exported": True,
            "actual_emergent_observer_closure": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
