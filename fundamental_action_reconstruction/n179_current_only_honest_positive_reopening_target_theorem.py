#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n179_current_only_honest_positive_reopening_target_theorem_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n163 = load_json(
        "fundamental_action_reconstruction/generated/n163_current_observer_information_deficit_downstream_symptom_theorem_summary.json"
    )
    n164 = load_json(
        "fundamental_action_reconstruction/generated/n164_current_selector_construction_stop_condition_theorem_summary.json"
    )
    n166 = load_json(
        "fundamental_action_reconstruction/generated/n166_current_only_honest_positive_work_contract_theorem_summary.json"
    )
    n178 = load_json(
        "fundamental_action_reconstruction/generated/n178_current_fixed_first_contract_compliant_additive_attempt_nonreopening_theorem_summary.json"
    )
    p160 = load_json(
        "fundamental_action_reconstruction/generated/p160_current_first_preobserver_light_matter_source_provider_class_target_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "n163_observer_information_deficit_downstream",
            "actual": n163["theorem_result"]["observer_information_deficit_downstream_symptom_on_current_repo_state"],
            "expected": True,
            "meaning": "N163 already keeps observer information deficit downstream",
        },
        {
            "id": "n164_lane_stopped",
            "actual": n164["theorem_result"]["current_selector_construction_lane_should_be_treated_as_stopped"],
            "expected": True,
            "meaning": "N164 already stops the current selector-construction lane",
        },
        {
            "id": "n166_additive_upstream_contract",
            "actual": [
                n166["theorem_result"]["remaining_positive_work_must_be_genuinely_additive"],
                n166["theorem_result"]["remaining_positive_work_must_stay_upstream_of_observer"],
                n166["theorem_result"]["remaining_positive_work_must_be_kernel_split_safe"],
                n166["theorem_result"]["remaining_positive_work_must_be_source_object_first"],
            ],
            "expected": [True, True, True, True],
            "meaning": "N166 already fixes the contract for honest positive work",
        },
        {
            "id": "n178_fixed_attempt_nonreopening",
            "actual": n178["theorem_result"]["current_fixed_first_contract_compliant_additive_attempt_does_not_reopen_constructive_selector_frontier"],
            "expected": True,
            "meaning": "N178 already excludes the fixed attempt as a reopening route",
        },
        {
            "id": "p160_provider_class_target",
            "actual": p160["provider_class_target"],
            "expected": "preobserver_light_matter_source_provider_class_v1",
            "meaning": "P160 already reduces the reopening target to one explicit provider class",
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
                "meaning": item["meaning"],
            }
        )
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        summary = {
            "step": "N179",
            "status": "N179_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ONLY_POSITIVE_REOPENING_TARGET_STATE",
            "scope": "current_only_honest_positive_reopening_target_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N179",
            "status": "N179_DISCHARGED_CURRENT_ONLY_HONEST_POSITIVE_REOPENING_TARGET_THEOREM_NO_FALSE_PASS",
            "scope": "current_only_honest_positive_reopening_target_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "only_honest_positive_reopening_target": "preobserver_light_matter_source_provider_class_v1",
                "observer_information_deficit_remains_downstream": True,
                "stopped_selector_construction_lane_remains_stopped": True,
                "fixed_first_contract_compliant_additive_attempt_remains_nonreopening": True,
                "full_closure_pass": False,
            },
            "hard_limits": [
                "constructed_source_object_not_yet_exported",
                "admissible_S_sel_int_not_yet_constructed",
                "admissible_E_orient_not_yet_constructed",
                "downstream_chain_not_yet_constructed",
                "no_strict_core_selector_closure",
                "no_QW2191_discharge",
                "no_ToE_closure",
            ],
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
