#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p210_current_next_actual_emergent_observer_closure_realization_map_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    n229 = load_json(
        GENERATED / "n229_current_next_actual_emergent_observer_closure_commit_map_theorem_summary.json"
    )
    f122 = load_json(
        GENERATED / "f122_next_actual_emergent_observer_closure_realization_map_packet_summary.json"
    )

    props = f122["observer_actual_closure_realization_properties"]
    response = f122["source_observer_actual_closure_realization"]
    checks_spec = [
        {
            "id": "actual_closure_commit_is_admissible",
            "actual": n229["theorem_result"]["admissible_AL_obs_actual_closure_commit"],
            "expected": True,
        },
        {
            "id": "derived_only_from_actual_closure_commit_state",
            "actual": props["derived_only_from_actual_closure_commit_state"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": props["strict_core_only"],
            "expected": True,
        },
        {
            "id": "downstream_actual_closure_realization_only",
            "actual": props["downstream_actual_closure_realization_only"]
            and not props["actual_emergent_observer_closure"],
            "expected": True,
        },
        {
            "id": "positive_actual_closure_realization_amplitude",
            "actual": response["positive_actual_closure_realization_amplitude"],
            "expected": True,
        },
        {
            "id": "actual_closure_residual_annihilated",
            "actual": response["actual_closure_residual_annihilated"],
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
            "stage": "P210",
            "lane": "current_next_actual_emergent_observer_closure_realization_only",
            "status": "P210_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_REALIZATION_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P210",
            "lane": "current_next_actual_emergent_observer_closure_realization_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_REALIZATION_MAP_FROM_AL_OBS_ACTUAL_CLOSURE_COMMIT_PRELM_V1_AFTER_P210",
            "source_object": "S_preLM_strict_core_source_object_v1",
            "observer_actual_closure_commit_operator": "AL_obs_actual_closure_commit_preLM_v1",
            "observer_actual_closure_realization_operator": "AM_obs_actual_closure_realization_preLM_v1",
            "checks": checks,
            "actual_closure_realization_exported": True,
            "actual_emergent_observer_closure": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
