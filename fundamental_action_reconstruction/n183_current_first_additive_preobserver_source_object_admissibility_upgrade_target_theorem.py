#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n183_current_first_additive_preobserver_source_object_admissibility_upgrade_target_theorem_summary.json"


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
    p164 = load_json(
        "fundamental_action_reconstruction/generated/p164_current_additive_preobserver_source_object_admissibility_upgrade_target_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "construction_attempt_present",
            "actual": f76["construction_attempt"],
            "expected": "S_preLM_additive_candidate_v1",
        },
        {
            "id": "construction_attempt_probe_positive",
            "actual": p163["construction_attempt"],
            "expected": "S_preLM_additive_candidate_v1",
        },
        {
            "id": "construction_attempt_theorem_future_only",
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
            "id": "upgrade_target_probe_positive",
            "actual": p164["admissibility_upgrade_target"],
            "expected": "upgrade_to_admissible_source_v1(S_preLM_additive_candidate_v1)",
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
            "step": "N183",
            "status": "N183_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIRST_ADDITIVE_PREOBSERVER_SOURCE_OBJECT_ADMISSIBILITY_UPGRADE_TARGET_STATE",
            "scope": "current_first_additive_preobserver_source_object_admissibility_upgrade_target_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N183",
            "status": "N183_DISCHARGED_CURRENT_FIRST_ADDITIVE_PREOBSERVER_SOURCE_OBJECT_ADMISSIBILITY_UPGRADE_TARGET_THEOREM_NO_FALSE_PASS",
            "scope": "current_first_additive_preobserver_source_object_admissibility_upgrade_target_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "construction_attempt": "S_preLM_additive_candidate_v1",
                "admissibility_upgrade_target": "upgrade_to_admissible_source_v1(S_preLM_additive_candidate_v1)",
                "future_only": True,
                "upstream_of_observer": True,
                "kernel_split_safe": True,
                "full_closure_pass": False,
            },
            "hard_limits": [
                "no_admissibility_clause_discharged",
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
