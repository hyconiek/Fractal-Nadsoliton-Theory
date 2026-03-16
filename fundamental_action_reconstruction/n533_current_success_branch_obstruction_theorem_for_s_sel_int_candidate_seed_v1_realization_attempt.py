#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n533_current_success_branch_obstruction_theorem_for_s_sel_int_candidate_seed_v1_realization_attempt_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p642 = load_json(
        "fundamental_action_reconstruction/generated/p642_current_success_verdict_discharge_probe_for_s_sel_int_candidate_seed_v1_realization_attempt_summary.json"
    )

    checks_spec = [
        {
            "id": "success_branch_selected_as_next_attack",
            "actual": p642["target_state"]["success_branch_selected_as_next_attack"],
            "expected": True,
            "meaning": "P642 stays on the remaining success branch ordering (v1)",
        },
        {
            "id": "explicit_success_verdict_exported",
            "actual": p642["target_state"]["explicit_success_verdict_exported"],
            "expected": False,
            "meaning": "the current repo still does not export an explicit success verdict (v1)",
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
            "step": "N533",
            "status": "N533_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_SUCCESS_BRANCH_STATE_FOR_SEED_V1",
            "scope": "current_success_verdict_discharge_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
            "no_false_pass": True,
        }
    else:
        summary = {
            "step": "N533",
            "status": "N533_DISCHARGED_CURRENT_SUCCESS_BRANCH_OBSTRUCTION_THEOREM_FOR_S_SEL_INT_REALIZATION_ATTEMPT_V1_NO_FALSE_PASS",
            "scope": "current_success_verdict_discharge_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "explicit_success_verdict_exported": False,
                "remaining_open_branches": p642["target_state"]["remaining_open_branches"],
                "full_closure_pass": False,
            },
            "remaining_open_branches": p642["target_state"]["remaining_open_branches"],
            "hard_limits": [
                "explicit_success_verdict_not_yet_exported",
                "constructed_source_object_not_yet_exported",
                "admissible_S_sel_int_not_yet_constructed",
                "admissible_E_orient_not_yet_constructed",
                "downstream_chain_not_yet_constructed",
                "no_strict_core_selector_closure",
                "no_QW2191_discharge",
                "no_ToE_closure",
            ],
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

