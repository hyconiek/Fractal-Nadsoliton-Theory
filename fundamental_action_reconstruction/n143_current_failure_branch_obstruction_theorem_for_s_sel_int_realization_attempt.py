#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n143_current_failure_branch_obstruction_theorem_for_s_sel_int_realization_attempt_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p132 = load_json(
        "fundamental_action_reconstruction/generated/p132_current_failure_verdict_discharge_probe_for_s_sel_int_realization_attempt_summary.json"
    )

    checks_spec = [
        {
            "id": "failure_branch_selected_as_first_attack",
            "actual": p132["target_state"]["failure_branch_selected_as_first_attack"],
            "expected": True,
            "meaning": "P132 stays on the failure-first branch ordering",
        },
        {
            "id": "explicit_failure_verdict_exported",
            "actual": p132["target_state"]["explicit_failure_verdict_exported"],
            "expected": False,
            "meaning": "the current repo still does not export an explicit failure verdict",
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
            "step": "N143",
            "status": "N143_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FAILURE_BRANCH_STATE",
            "scope": "current_failure_verdict_discharge_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N143",
            "status": "N143_DISCHARGED_CURRENT_FAILURE_BRANCH_OBSTRUCTION_THEOREM_FOR_S_SEL_INT_REALIZATION_ATTEMPT_NO_FALSE_PASS",
            "scope": "current_failure_verdict_discharge_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "explicit_failure_verdict_exported": False,
                "remaining_open_branches": p132["target_state"]["remaining_open_branches"],
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "future_success_verdict_discharge_for_S_sel_int_new_object_constructed_realization_attempt_v0",
                "future_admissibility_test_of_a_future_constructed_source_object",
                "future_derivation_of_admissible_E_orient_from_a_future_new_source_object",
                "future_completion_of_B_sel_R_sel_O_sel_after_new_source_object_construction",
            ],
            "hard_limits": [
                "explicit_failure_verdict_not_yet_exported",
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
