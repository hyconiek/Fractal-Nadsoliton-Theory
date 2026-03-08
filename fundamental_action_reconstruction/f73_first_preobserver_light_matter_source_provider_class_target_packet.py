#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f73_first_preobserver_light_matter_source_provider_class_target_packet_summary.json"


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

    checks_spec = [
        {
            "id": "n163_observer_deficit_downstream",
            "actual": n163["theorem_result"]["observer_information_deficit_downstream_symptom_on_current_repo_state"],
            "expected": True,
            "meaning": "N163 already keeps observer information deficit downstream",
        },
        {
            "id": "n163_primary_gap_upstream",
            "actual": n163["theorem_result"]["primary_missing_selector_source_gap_upstream_of_observer"],
            "expected": True,
            "meaning": "N163 already keeps the primary selector gap upstream of observer",
        },
        {
            "id": "n164_lane_stopped",
            "actual": n164["theorem_result"]["current_selector_construction_lane_should_be_treated_as_stopped"],
            "expected": True,
            "meaning": "N164 already stops the current selector-construction lane",
        },
        {
            "id": "n166_positive_work_additive",
            "actual": n166["theorem_result"]["remaining_positive_work_must_be_genuinely_additive"],
            "expected": True,
            "meaning": "N166 already requires genuinely additive positive work",
        },
        {
            "id": "n166_positive_work_upstream_of_observer",
            "actual": n166["theorem_result"]["remaining_positive_work_must_stay_upstream_of_observer"],
            "expected": True,
            "meaning": "N166 already requires positive work to stay upstream of observer",
        },
        {
            "id": "n178_fixed_attempt_nonreopening",
            "actual": n178["theorem_result"]["current_fixed_first_contract_compliant_additive_attempt_does_not_reopen_constructive_selector_frontier"],
            "expected": True,
            "meaning": "N178 already excludes the fixed first contract-compliant additive attempt as a reopening route",
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
            "stage": "F73",
            "lane": "future_preobserver_light_matter_source_provider_class_target_only",
            "status": "F73_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_PROVIDER_CLASS_TARGET_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "F73",
            "lane": "future_preobserver_light_matter_source_provider_class_target_only",
            "goal": "freeze_one_noncyclic_provider_class_target_for_future_additive_upstream_source_work_while_keeping_observer_information_deficit_downstream",
            "status": "F73_EXECUTED_FIRST_PREOBSERVER_LIGHT_MATTER_SOURCE_PROVIDER_CLASS_TARGET_PACKET_NO_FALSE_PASS",
            "provider_class_target": "preobserver_light_matter_source_provider_class_v1",
            "ordering": [
                "nadsoliton",
                "light",
                "matter",
                "emergent_observer",
            ],
            "provider_class_contract": [
                "future_only",
                "genuinely_additive",
                "upstream_of_observer",
                "light_before_observer",
                "matter_as_encoding_not_primary_selector_source",
                "kernel_split_safe",
                "no_external_selector_import",
                "source_object_first",
            ],
            "checks": checks,
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
