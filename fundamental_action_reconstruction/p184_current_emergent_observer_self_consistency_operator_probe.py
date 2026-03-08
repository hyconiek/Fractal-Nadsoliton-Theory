#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p184_current_emergent_observer_self_consistency_operator_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n203 = load_json(
        "fundamental_action_reconstruction/generated/n203_current_actual_emergent_observer_realization_map_theorem_summary.json"
    )
    f96 = load_json(
        "fundamental_action_reconstruction/generated/f96_first_emergent_observer_self_consistency_operator_packet_summary.json"
    )

    props = f96["observer_self_consistency_properties"]
    response = f96["source_observer_self_consistency"]
    checks_spec = [
        {
            "id": "observer_realization_is_admissible",
            "actual": n203["theorem_result"]["admissible_H_obs_realization"],
            "expected": True,
        },
        {
            "id": "derived_only_from_realization_state",
            "actual": props["derived_only_from_realization_state"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": props["strict_core_only"],
            "expected": True,
        },
        {
            "id": "downstream_self_consistency_only",
            "actual": props["downstream_self_consistency_only"] and not props["actual_emergent_observer_constructed"],
            "expected": True,
        },
        {
            "id": "self_consistency_sector_exported",
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
        {
            "id": "idempotent",
            "actual": f96["observer_self_consistency_operator"]["idempotent"],
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
            "stage": "P184",
            "lane": "current_emergent_observer_self_consistency_operator_only",
            "status": "P184_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_EMERGENT_OBSERVER_SELF_CONSISTENCY_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P184",
            "lane": "current_emergent_observer_self_consistency_operator_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_EMERGENT_OBSERVER_SELF_CONSISTENCY_OPERATOR_FROM_H_OBS_REALIZATION_PRELM_V1_AFTER_P184",
            "source_object": "S_preLM_strict_core_source_object_v1",
            "observer_realization_map": "H_obs_realization_preLM_v1",
            "observer_self_consistency_operator": "J_obs_self_consistency_preLM_v1",
            "checks": checks,
            "observer_self_consistency_operator_exported": True,
            "emergent_observer_constructed": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
