#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f63_current_future_additive_upstream_source_work_contract_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n123 = load_json(
        "fundamental_action_reconstruction/generated/n123_current_legacy_to_strict_kernel_package_level_nonbridge_theorem_summary.json"
    )
    n163 = load_json(
        "fundamental_action_reconstruction/generated/n163_current_observer_information_deficit_downstream_symptom_theorem_summary.json"
    )
    n165 = load_json(
        "fundamental_action_reconstruction/generated/n165_current_selector_construction_handoff_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n123_package_nonbridge",
            "actual": n123["theorem_result"]["package_level_nonbridge_on_current_repo_state"],
            "expected": True,
            "meaning": "N123 keeps the package-level legacy-strict nonbridge active",
        },
        {
            "id": "n163_observer_downstream",
            "actual": n163["theorem_result"][
                "observer_information_deficit_downstream_symptom_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N163 keeps observer information deficit downstream",
        },
        {
            "id": "n165_handoff_active",
            "actual": n165["theorem_result"][
                "stopped_selector_construction_lane_handed_off"
            ],
            "expected": True,
            "meaning": "N165 already hands off the stopped selector-construction lane",
        },
        {
            "id": "n165_remaining_move_class",
            "actual": n165["theorem_result"]["remaining_positive_move_class"],
            "expected": "future_genuinely_additive_upstream_source_work_only",
            "meaning": "N165 already restricts remaining positive work to genuinely additive upstream source work only",
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
            "stage": "F63",
            "lane": "current_future_additive_upstream_source_work_contract_only",
            "status": "F63_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ADDITIVE_UPSTREAM_WORK_CONTRACT_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "F63",
            "lane": "current_future_additive_upstream_source_work_contract_only",
            "goal": "freeze_the_only_remaining_honest_positive_work_contract_after_selector_construction_handoff",
            "status": "F63_EXECUTED_CURRENT_FUTURE_ADDITIVE_UPSTREAM_SOURCE_WORK_CONTRACT_PACKET_NO_FALSE_PASS",
            "checks": checks,
            "contract": {
                "current_future_additive_upstream_source_work_contract_active": True,
                "admitted_contract": {
                    "genuinely_additive": True,
                    "upstream_of_observer": True,
                    "kernel_split_safe": True,
                    "no_external_selector_import": True,
                    "strict_core_source_object_first": True,
                },
            },
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
