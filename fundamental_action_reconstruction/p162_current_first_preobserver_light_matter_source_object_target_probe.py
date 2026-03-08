#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p162_current_first_preobserver_light_matter_source_object_target_probe_summary.json"


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

    checks_spec = [
        {
            "id": "f74_provider_packet",
            "actual": f74["provider_packet"],
            "expected": "preobserver_light_matter_source_provider_packet_v1",
            "meaning": "F74 already exports the provider packet",
        },
        {
            "id": "p161_provider_packet_consistent",
            "actual": p161["provider_packet"],
            "expected": "preobserver_light_matter_source_provider_packet_v1",
            "meaning": "P161 already confirms packet consistency",
        },
        {
            "id": "n180_future_only",
            "actual": [
                n180["theorem_result"]["guardrail_consistent"],
                n180["theorem_result"]["future_only"],
            ],
            "expected": [True, True],
            "meaning": "N180 already keeps the packet future-only and guardrail-consistent",
        },
        {
            "id": "f75_target_exported",
            "actual": f75["source_object_target"],
            "expected": "preobserver_light_matter_source_object_target_v1",
            "meaning": "F75 exports one source-object target",
        },
        {
            "id": "f75_observer_excluded",
            "actual": [
                f75["carrier"]["observer_excluded"],
                f75["guardrails"]["upstream_of_observer"],
                f75["guardrails"]["source_object_first"],
            ],
            "expected": [True, True, True],
            "meaning": "F75 keeps observer excluded and source-object-first",
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
            "stage": "P162",
            "lane": "future_preobserver_light_matter_source_object_target_probe_only",
            "status": "P162_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_SOURCE_OBJECT_TARGET_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P162",
            "lane": "future_preobserver_light_matter_source_object_target_probe_only",
            "goal": "test_whether_the_first_preobserver_light_matter_source_object_target_stays_guardrail_consistent_and_upstream_of_observer",
            "status": "CURRENT_REPO_EXPORTS_ONE_INTERNAL_GUARDRAIL_CONSISTENT_PREOBSERVER_LIGHT_MATTER_SOURCE_OBJECT_TARGET_AFTER_P162",
            "provider_packet": "preobserver_light_matter_source_provider_packet_v1",
            "source_object_target": "preobserver_light_matter_source_object_target_v1",
            "checks": checks,
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
