#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n152_only_remaining_positive_move_class_reduced_to_one_first_additive_construction_attempt_instance_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p139 = load_json(
        "fundamental_action_reconstruction/generated/p139_first_future_genuinely_additive_source_object_construction_attempt_instance_probe_summary.json"
    )
    n151 = load_json(
        "fundamental_action_reconstruction/generated/n151_only_remaining_positive_move_class_reduced_to_one_first_additive_construction_attempt_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n151_target_reduction",
            "actual": n151["theorem_result"]["explicit_future_additive_attempt_target"],
            "expected": "S_sel_int_additive_attempt_target_v1",
            "meaning": "N151 already fixes one explicit future additive-attempt target",
        },
        {
            "id": "p139_attempt_instance",
            "actual": p139["target_state"]["explicit_future_additive_attempt_instance"],
            "expected": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "meaning": "P139 reduces the move class to one explicit construction-attempt instance",
        },
        {
            "id": "p139_single_remaining_branch",
            "actual": p139["target_state"]["remaining_open_branches"],
            "expected": [
                "future_success_or_failure_evaluation_of_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)"
            ],
            "meaning": "P139 leaves only one future evaluation branch",
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
            "step": "N152",
            "status": "N152_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FUTURE_ADDITIVE_ATTEMPT_INSTANCE_REDUCTION_STATE",
            "scope": "only_remaining_positive_move_class_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N152",
            "status": "N152_DISCHARGED_ONLY_REMAINING_POSITIVE_MOVE_CLASS_REDUCED_TO_ONE_FIRST_ADDITIVE_CONSTRUCTION_ATTEMPT_INSTANCE_THEOREM_NO_FALSE_PASS",
            "scope": "only_remaining_positive_move_class_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "only_remaining_positive_move_class": "future_genuinely_additive_new_strict_core_source_object_construction",
                "explicit_future_additive_attempt_instance": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "future_success_or_failure_evaluation_of_construct_attempt_v1(S_sel_int_additive_attempt_target_v1)"
            ],
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

    OUT.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    print(OUT)


if __name__ == "__main__":
    main()
