#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n180_current_first_guardrail_consistent_preobserver_provider_packet_theorem_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f73 = load_json(
        "fundamental_action_reconstruction/generated/f73_first_preobserver_light_matter_source_provider_class_target_packet_summary.json"
    )
    f74 = load_json(
        "fundamental_action_reconstruction/generated/f74_first_preobserver_light_matter_asymmetric_provider_packet_summary.json"
    )
    p160 = load_json(
        "fundamental_action_reconstruction/generated/p160_current_first_preobserver_light_matter_source_provider_class_target_probe_summary.json"
    )
    p161 = load_json(
        "fundamental_action_reconstruction/generated/p161_current_preobserver_light_matter_asymmetric_provider_packet_consistency_probe_summary.json"
    )
    n179 = load_json(
        "fundamental_action_reconstruction/generated/n179_current_only_honest_positive_reopening_target_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "f73_provider_class_target",
            "actual": f73["provider_class_target"],
            "expected": "preobserver_light_matter_source_provider_class_v1",
            "meaning": "F73 already fixes the provider class target",
        },
        {
            "id": "n179_only_reopening_target",
            "actual": n179["theorem_result"]["only_honest_positive_reopening_target"],
            "expected": "preobserver_light_matter_source_provider_class_v1",
            "meaning": "N179 already reduces reopening to this one provider class",
        },
        {
            "id": "f74_provider_packet",
            "actual": f74["provider_packet"],
            "expected": "preobserver_light_matter_source_provider_packet_v1",
            "meaning": "F74 already exports one explicit provider packet",
        },
        {
            "id": "p160_provider_class_probe",
            "actual": p160["provider_class_target"],
            "expected": "preobserver_light_matter_source_provider_class_v1",
            "meaning": "P160 already confirms the provider class target",
        },
        {
            "id": "p161_provider_packet_consistent",
            "actual": p161["provider_packet"],
            "expected": "preobserver_light_matter_source_provider_packet_v1",
            "meaning": "P161 already confirms the provider packet",
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
            "step": "N180",
            "status": "N180_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_PREOBSERVER_PROVIDER_PACKET_STATE",
            "scope": "current_first_guardrail_consistent_preobserver_provider_packet_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N180",
            "status": "N180_DISCHARGED_CURRENT_FIRST_GUARDRAIL_CONSISTENT_PREOBSERVER_PROVIDER_PACKET_THEOREM_NO_FALSE_PASS",
            "scope": "current_first_guardrail_consistent_preobserver_provider_packet_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "provider_class_target": "preobserver_light_matter_source_provider_class_v1",
                "provider_packet": "preobserver_light_matter_source_provider_packet_v1",
                "guardrail_consistent": True,
                "future_only": True,
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
