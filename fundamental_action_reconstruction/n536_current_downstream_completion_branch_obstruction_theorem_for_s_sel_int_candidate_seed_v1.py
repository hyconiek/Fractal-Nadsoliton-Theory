#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n536_current_downstream_completion_branch_obstruction_theorem_for_s_sel_int_candidate_seed_v1_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p645 = load_json(
        "fundamental_action_reconstruction/generated/p645_current_downstream_completion_branch_discharge_probe_for_s_sel_int_candidate_seed_v1_summary.json"
    )

    checks_spec = [
        {
            "id": "downstream_completion_branch_selected_as_next_attack",
            "actual": p645["target_state"]["downstream_completion_branch_selected_as_next_attack"],
            "expected": True,
            "meaning": "P645 stays on the seed-v1 downstream-completion ordering",
        },
        {
            "id": "explicit_downstream_completion_branch_discharge_exported",
            "actual": p645["target_state"]["explicit_downstream_completion_branch_discharge_exported"],
            "expected": True,
            "meaning": "the current repo exports an explicit downstream-completion branch discharge (seed v1)",
        },
        {
            "id": "remaining_open_branches_empty",
            "actual": p645["target_state"]["remaining_open_branches"],
            "expected": [],
            "meaning": "after P645 the lower-branch list is exhausted on the current repo state (seed v1)",
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
            "step": "N536",
            "status": "N536_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_DOWNSTREAM_COMPLETION_BRANCH_STATE_FOR_SEED_V1",
            "scope": "current_downstream_completion_branch_discharge_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
            "no_false_pass": True,
        }
    else:
        summary = {
            "step": "N536",
            "status": "N536_DISCHARGED_CURRENT_DOWNSTREAM_COMPLETION_BRANCH_OBSTRUCTION_THEOREM_FOR_SEED_V1_NO_FALSE_PASS",
            "scope": "current_downstream_completion_branch_discharge_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "explicit_downstream_completion_branch_discharge_exported": True,
                "remaining_open_branches": [],
                "full_closure_pass": False,
            },
            "remaining_open_branches": [],
            "hard_limits": [
                "admissible_S_sel_int_not_yet_constructed",
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
