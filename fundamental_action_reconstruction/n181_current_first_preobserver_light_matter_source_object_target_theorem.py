#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n181_current_first_preobserver_light_matter_source_object_target_theorem_summary.json"


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

    checks_spec = [
        {
            "id": "f74_provider_packet",
            "actual": f74["provider_packet"],
            "expected": "preobserver_light_matter_source_provider_packet_v1",
            "meaning": "F74 exports the provider packet",
        },
        {
            "id": "p161_provider_packet_consistent",
            "actual": p161["provider_packet"],
            "expected": "preobserver_light_matter_source_provider_packet_v1",
            "meaning": "P161 confirms provider-packet consistency",
        },
        {
            "id": "n180_guardrail_consistent",
            "actual": [
                n180["theorem_result"]["guardrail_consistent"],
                n180["theorem_result"]["future_only"],
            ],
            "expected": [True, True],
            "meaning": "N180 keeps the packet future-only and guardrail-consistent",
        },
        {
            "id": "f75_source_target",
            "actual": f75["source_object_target"],
            "expected": "preobserver_light_matter_source_object_target_v1",
            "meaning": "F75 exports one explicit source-object target",
        },
        {
            "id": "p162_target_consistent",
            "actual": p162["source_object_target"],
            "expected": "preobserver_light_matter_source_object_target_v1",
            "meaning": "P162 confirms the source-object target",
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
            "step": "N181",
            "status": "N181_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIRST_PREOBSERVER_SOURCE_OBJECT_TARGET_STATE",
            "scope": "current_first_preobserver_light_matter_source_object_target_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N181",
            "status": "N181_DISCHARGED_CURRENT_FIRST_PREOBSERVER_LIGHT_MATTER_SOURCE_OBJECT_TARGET_THEOREM_NO_FALSE_PASS",
            "scope": "current_first_preobserver_light_matter_source_object_target_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "provider_packet": "preobserver_light_matter_source_provider_packet_v1",
                "source_object_target": "preobserver_light_matter_source_object_target_v1",
                "future_only": True,
                "upstream_of_observer": True,
                "observer_excluded_from_target_carrier": True,
                "kernel_split_safe": True,
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
