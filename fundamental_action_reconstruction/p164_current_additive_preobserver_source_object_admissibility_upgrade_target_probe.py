#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p164_current_additive_preobserver_source_object_admissibility_upgrade_target_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f76 = load_json(
        "fundamental_action_reconstruction/generated/f76_first_additive_preobserver_source_object_construction_attempt_packet_summary.json"
    )
    p163 = load_json(
        "fundamental_action_reconstruction/generated/p163_current_additive_preobserver_source_object_construction_attempt_probe_summary.json"
    )
    n182 = load_json(
        "fundamental_action_reconstruction/generated/n182_current_first_additive_preobserver_source_object_construction_attempt_theorem_summary.json"
    )
    f77 = load_json(
        "fundamental_action_reconstruction/generated/f77_first_additive_preobserver_source_object_admissibility_upgrade_target_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "additive_attempt_present",
            "actual": f76["construction_attempt"],
            "expected": "S_preLM_additive_candidate_v1",
        },
        {
            "id": "additive_attempt_probe_positive",
            "actual": p163["construction_attempt"],
            "expected": "S_preLM_additive_candidate_v1",
        },
        {
            "id": "additive_attempt_theorem_future_only",
            "actual": [
                n182["theorem_result"]["future_only"],
                n182["theorem_result"]["observer_excluded_from_carrier"],
            ],
            "expected": [True, True],
        },
        {
            "id": "upgrade_target_present",
            "actual": f77["admissibility_upgrade_target"],
            "expected": "upgrade_to_admissible_source_v1(S_preLM_additive_candidate_v1)",
        },
        {
            "id": "upgrade_target_guardrails",
            "actual": [
                f77["guardrails"]["future_only"],
                f77["guardrails"]["upstream_of_observer"],
                f77["guardrails"]["kernel_split_safe"],
                f77["guardrails"]["no_external_selector_import"],
            ],
            "expected": [True, True, True, True],
        },
        {
            "id": "upgrade_target_reuses_f34",
            "actual": f77["reused_contract_source"],
            "expected": "F34",
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
            "stage": "P164",
            "lane": "future_additive_preobserver_source_object_admissibility_upgrade_target_probe_only",
            "status": "P164_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ADDITIVE_PREOBSERVER_SOURCE_OBJECT_ADMISSIBILITY_UPGRADE_TARGET_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P164",
            "lane": "future_additive_preobserver_source_object_admissibility_upgrade_target_probe_only",
            "goal": "test_whether_one_explicit_admissibility_upgrade_target_is_now_exported_above_the_first_additive_preobserver_source_object_construction_attempt",
            "status": "CURRENT_REPO_EXPORTS_ONE_INTERNAL_GUARDRAIL_CONSISTENT_ADDITIVE_PREOBSERVER_SOURCE_OBJECT_ADMISSIBILITY_UPGRADE_TARGET_AFTER_P164",
            "construction_attempt": "S_preLM_additive_candidate_v1",
            "admissibility_upgrade_target": "upgrade_to_admissible_source_v1(S_preLM_additive_candidate_v1)",
            "checks": checks,
            "admissible_S_sel_int": False,
            "admissible_E_orient": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
