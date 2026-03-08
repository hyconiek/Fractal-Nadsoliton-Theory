#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p163_current_additive_preobserver_source_object_construction_attempt_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f74 = load_json(
        "fundamental_action_reconstruction/generated/f74_first_preobserver_light_matter_asymmetric_provider_packet_summary.json"
    )
    p161 = load_json(
        "fundamental_action_reconstruction/generated/p161_current_preobserver_light_matter_asymmetric_provider_packet_consistency_probe_summary.json"
    )
    n180 = load_json(
        "fundamental_action_reconstruction/generated/n180_current_first_guardrail_consistent_preobserver_provider_packet_theorem_summary.json"
    )
    f75 = load_json(
        "fundamental_action_reconstruction/generated/f75_first_preobserver_light_matter_source_object_target_packet_summary.json"
    )
    p162 = load_json(
        "fundamental_action_reconstruction/generated/p162_current_first_preobserver_light_matter_source_object_target_probe_summary.json"
    )
    n181 = load_json(
        "fundamental_action_reconstruction/generated/n181_current_first_preobserver_light_matter_source_object_target_theorem_summary.json"
    )
    f76 = load_json(
        "fundamental_action_reconstruction/generated/f76_first_additive_preobserver_source_object_construction_attempt_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "provider_packet_still_present",
            "actual": f74["provider_packet"],
            "expected": "preobserver_light_matter_source_provider_packet_v1",
        },
        {
            "id": "provider_packet_still_consistent",
            "actual": p161["provider_packet"],
            "expected": "preobserver_light_matter_source_provider_packet_v1",
        },
        {
            "id": "provider_packet_future_only",
            "actual": [
                n180["theorem_result"]["guardrail_consistent"],
                n180["theorem_result"]["future_only"],
            ],
            "expected": [True, True],
        },
        {
            "id": "source_target_still_present",
            "actual": f75["source_object_target"],
            "expected": "preobserver_light_matter_source_object_target_v1",
        },
        {
            "id": "source_target_probe_still_positive",
            "actual": p162["source_object_target"],
            "expected": "preobserver_light_matter_source_object_target_v1",
        },
        {
            "id": "source_target_theorem_future_only",
            "actual": [
                n181["theorem_result"]["future_only"],
                n181["theorem_result"]["observer_excluded_from_target_carrier"],
            ],
            "expected": [True, True],
        },
        {
            "id": "additive_attempt_present",
            "actual": f76["construction_attempt"],
            "expected": "S_preLM_additive_candidate_v1",
        },
        {
            "id": "additive_attempt_guardrails",
            "actual": [
                f76["guardrails"]["future_only"],
                f76["guardrails"]["upstream_of_observer"],
                f76["guardrails"]["kernel_split_safe"],
                f76["guardrails"]["no_external_selector_import"],
            ],
            "expected": [True, True, True, True],
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
            "stage": "P163",
            "lane": "future_additive_preobserver_source_object_construction_attempt_probe_only",
            "status": "P163_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ADDITIVE_PREOBSERVER_SOURCE_OBJECT_ATTEMPT_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P163",
            "lane": "future_additive_preobserver_source_object_construction_attempt_probe_only",
            "goal": "test_whether_one_explicit_additive_preobserver_source_object_construction_attempt_is_now_exported_without_observer_promotion",
            "status": "CURRENT_REPO_EXPORTS_ONE_INTERNAL_GUARDRAIL_CONSISTENT_ADDITIVE_PREOBSERVER_SOURCE_OBJECT_CONSTRUCTION_ATTEMPT_AFTER_P163",
            "construction_attempt": "S_preLM_additive_candidate_v1",
            "checks": checks,
            "constructed_source_object": False,
            "admissible_S_sel_int": False,
            "admissible_E_orient": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
