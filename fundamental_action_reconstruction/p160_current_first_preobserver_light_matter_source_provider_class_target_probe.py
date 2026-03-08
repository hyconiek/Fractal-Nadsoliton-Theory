#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p160_current_first_preobserver_light_matter_source_provider_class_target_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f73 = load_json(
        "fundamental_action_reconstruction/generated/f73_first_preobserver_light_matter_source_provider_class_target_packet_summary.json"
    )
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

    checks_spec = [
        {
            "id": "f73_provider_class_target",
            "actual": f73["provider_class_target"],
            "expected": "preobserver_light_matter_source_provider_class_v1",
            "meaning": "F73 already freezes one explicit provider class target",
        },
        {
            "id": "n163_observer_downstream",
            "actual": n163["theorem_result"]["observer_information_deficit_downstream_symptom_on_current_repo_state"],
            "expected": True,
            "meaning": "N163 already keeps observer information deficit downstream",
        },
        {
            "id": "n164_lane_stopped",
            "actual": n164["theorem_result"]["current_selector_construction_lane_should_be_treated_as_stopped"],
            "expected": True,
            "meaning": "N164 already stops the old selector-construction lane",
        },
        {
            "id": "n166_upstream_additive_contract",
            "actual": [
                n166["theorem_result"]["remaining_positive_work_must_be_genuinely_additive"],
                n166["theorem_result"]["remaining_positive_work_must_stay_upstream_of_observer"],
                n166["theorem_result"]["remaining_positive_work_must_be_source_object_first"],
            ],
            "expected": [True, True, True],
            "meaning": "N166 already fixes the contract for honest positive work",
        },
        {
            "id": "n178_fixed_attempt_nonreopening",
            "actual": n178["theorem_result"]["current_fixed_first_contract_compliant_additive_attempt_does_not_reopen_constructive_selector_frontier"],
            "expected": True,
            "meaning": "N178 already removes the fixed attempt as a reopening route",
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
            "stage": "P160",
            "lane": "future_preobserver_light_matter_source_provider_class_target_probe_only",
            "status": "P160_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_PROVIDER_CLASS_TARGET_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P160",
            "lane": "future_preobserver_light_matter_source_provider_class_target_probe_only",
            "goal": "test_whether_the_next_honest_positive_reopening_target_is_one_explicit_future_preobserver_light_matter_source_provider_class",
            "status": "CURRENT_REPO_REDUCES_THE_NEXT_HONEST_POSITIVE_REOPENING_TARGET_TO_ONE_EXPLICIT_PREOBSERVER_LIGHT_MATTER_SOURCE_PROVIDER_CLASS_AFTER_P160",
            "provider_class_target": "preobserver_light_matter_source_provider_class_v1",
            "checks": checks,
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
