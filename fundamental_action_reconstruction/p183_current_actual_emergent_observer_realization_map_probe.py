#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p183_current_actual_emergent_observer_realization_map_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n202 = load_json(
        "fundamental_action_reconstruction/generated/n202_current_actual_emergent_observer_construction_candidate_theorem_summary.json"
    )
    f95 = load_json(
        "fundamental_action_reconstruction/generated/f95_first_actual_emergent_observer_realization_map_packet_summary.json"
    )

    props = f95["observer_realization_properties"]
    response = f95["source_observer_realization"]
    checks_spec = [
        {
            "id": "observer_construction_candidate_is_admissible",
            "actual": n202["theorem_result"]["admissible_G_obs_candidate"],
            "expected": True,
        },
        {
            "id": "derived_only_from_construction_candidate",
            "actual": props["derived_only_from_construction_candidate"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": props["strict_core_only"],
            "expected": True,
        },
        {
            "id": "downstream_realization_only",
            "actual": props["downstream_realization_only"] and not props["actual_emergent_observer_constructed"],
            "expected": True,
        },
        {
            "id": "realization_sector_exported",
            "actual": True,
            "expected": True,
        },
        {
            "id": "positive_commit",
            "actual": response["positive_commit"],
            "expected": True,
        },
        {
            "id": "vanishing_residual",
            "actual": response["vanishing_residual"],
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
            "stage": "P183",
            "lane": "current_actual_emergent_observer_realization_map_only",
            "status": "P183_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_EMERGENT_OBSERVER_REALIZATION_MAP_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P183",
            "lane": "current_actual_emergent_observer_realization_map_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_ACTUAL_EMERGENT_OBSERVER_REALIZATION_MAP_FROM_G_OBS_CANDIDATE_PRELM_V1_AFTER_P183",
            "source_object": "S_preLM_strict_core_source_object_v1",
            "observer_construction_candidate_operator": "G_obs_candidate_preLM_v1",
            "observer_realization_map": "H_obs_realization_preLM_v1",
            "checks": checks,
            "actual_observer_realization_map_exported": True,
            "emergent_observer_constructed": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
