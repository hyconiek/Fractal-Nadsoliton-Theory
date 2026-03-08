#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p194_current_emergent_observer_closure_realization_object_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n213 = load_json(
        "fundamental_action_reconstruction/generated/n213_current_emergent_observer_closure_commit_map_theorem_summary.json"
    )
    f106 = load_json(
        "fundamental_action_reconstruction/generated/f106_first_emergent_observer_closure_realization_object_packet_summary.json"
    )

    props = f106["observer_closure_realization_properties"]
    response = f106["source_observer_closure_realization"]
    checks_spec = [
        {
            "id": "observer_closure_commit_is_admissible",
            "actual": n213["theorem_result"]["admissible_V_obs_closure_commit_map"],
            "expected": True,
        },
        {
            "id": "derived_only_from_closure_commit_state",
            "actual": props["derived_only_from_closure_commit_state"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": props["strict_core_only"],
            "expected": True,
        },
        {
            "id": "downstream_closure_realization_only",
            "actual": props["downstream_closure_realization_only"]
            and not props["actual_emergent_observer_closure"],
            "expected": True,
        },
        {
            "id": "positive_closure_realization_amplitude",
            "actual": response["positive_closure_realization_amplitude"],
            "expected": True,
        },
        {
            "id": "one_dimensional_closure_realization",
            "actual": response["one_dimensional_closure_realization"],
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
            "stage": "P194",
            "lane": "current_emergent_observer_closure_realization_only",
            "status": "P194_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_EMERGENT_OBSERVER_CLOSURE_REALIZATION_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P194",
            "lane": "current_emergent_observer_closure_realization_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_EMERGENT_OBSERVER_CLOSURE_REALIZATION_OBJECT_FROM_V_OBS_CLOSURE_COMMIT_PRELM_V1_AFTER_P194",
            "source_object": "S_preLM_strict_core_source_object_v1",
            "observer_closure_commit_map": "V_obs_closure_commit_preLM_v1",
            "observer_closure_realization_map": "W_obs_closure_realization_preLM_v1",
            "checks": checks,
            "closure_realization_exported": True,
            "actual_emergent_observer_closure": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
