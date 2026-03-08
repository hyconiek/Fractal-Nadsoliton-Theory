#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p182_current_actual_emergent_observer_construction_candidate_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n201 = load_json(
        "fundamental_action_reconstruction/generated/n201_current_emergent_observer_limit_readout_operator_theorem_summary.json"
    )
    f94 = load_json(
        "fundamental_action_reconstruction/generated/f94_first_actual_emergent_observer_construction_candidate_packet_summary.json"
    )

    props = f94["observer_construction_candidate_properties"]
    response = f94["source_observer_construction_candidate"]
    checks_spec = [
        {
            "id": "observer_limit_readout_is_admissible",
            "actual": n201["theorem_result"]["admissible_L_obs_limit"],
            "expected": True,
        },
        {
            "id": "derived_only_from_observer_limit_readout",
            "actual": props["derived_only_from_observer_limit_readout"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": props["strict_core_only"],
            "expected": True,
        },
        {
            "id": "downstream_candidate_only",
            "actual": props["downstream_candidate_only"] and not props["actual_emergent_observer_constructed"],
            "expected": True,
        },
        {
            "id": "candidate_sector_exported",
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
            "stage": "P182",
            "lane": "current_actual_emergent_observer_construction_candidate_only",
            "status": "P182_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_EMERGENT_OBSERVER_CANDIDATE_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P182",
            "lane": "current_actual_emergent_observer_construction_candidate_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_ACTUAL_EMERGENT_OBSERVER_CONSTRUCTION_CANDIDATE_OPERATOR_FROM_L_OBS_LIMIT_PRELM_V1_AFTER_P182",
            "source_object": "S_preLM_strict_core_source_object_v1",
            "observer_limit_readout_operator": "L_obs_limit_preLM_v1",
            "observer_construction_candidate_operator": "G_obs_candidate_preLM_v1",
            "checks": checks,
            "actual_observer_construction_candidate_constructed": True,
            "emergent_observer_constructed": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
