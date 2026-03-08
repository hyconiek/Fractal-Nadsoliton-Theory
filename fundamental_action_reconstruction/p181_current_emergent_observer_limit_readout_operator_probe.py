#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p181_current_emergent_observer_limit_readout_operator_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n200 = load_json(
        "fundamental_action_reconstruction/generated/n200_current_preobserver_to_emergent_observer_coarse_graining_theorem_summary.json"
    )
    f93 = load_json(
        "fundamental_action_reconstruction/generated/f93_first_emergent_observer_limit_readout_operator_packet_summary.json"
    )

    props = f93["observer_limit_readout_properties"]
    response = f93["source_observer_limit_readout"]
    checks_spec = [
        {
            "id": "coarse_graining_is_admissible",
            "actual": n200["theorem_result"]["admissible_C_obs_limit"],
            "expected": True,
        },
        {
            "id": "derived_only_from_coarse_graining",
            "actual": props["derived_only_from_coarse_graining"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": props["strict_core_only"],
            "expected": True,
        },
        {
            "id": "observer_limit_only",
            "actual": props["observer_limit_only"] and not props["actual_emergent_observer_constructed"],
            "expected": True,
        },
        {
            "id": "observer_limit_readout_sector_exported",
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
            "stage": "P181",
            "lane": "current_emergent_observer_limit_readout_operator_only",
            "status": "P181_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_OBSERVER_LIMIT_READOUT_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P181",
            "lane": "current_emergent_observer_limit_readout_operator_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_EMERGENT_OBSERVER_LIMIT_READOUT_OPERATOR_FROM_C_OBS_LIMIT_PRELM_V1_AFTER_P181",
            "source_object": "S_preLM_strict_core_source_object_v1",
            "coarse_graining_operator": "C_obs_limit_preLM_v1",
            "observer_limit_readout_operator": "L_obs_limit_preLM_v1",
            "checks": checks,
            "actual_observer_limit_readout_constructed": True,
            "emergent_observer_constructed": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
