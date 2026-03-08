#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n145_current_admissibility_branch_obstruction_theorem_for_future_constructed_source_object_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p134 = load_json(
        "fundamental_action_reconstruction/generated/p134_current_admissibility_branch_discharge_probe_for_future_constructed_source_object_summary.json"
    )

    checks_spec = [
        {
            "id": "admissibility_branch_selected_as_next_attack",
            "actual": p134["target_state"]["admissibility_branch_selected_as_next_attack"],
            "expected": True,
            "meaning": "P134 stays on the post-verdict admissibility-first ordering",
        },
        {
            "id": "explicit_admissibility_branch_discharge_exported",
            "actual": p134["target_state"]["explicit_admissibility_branch_discharge_exported"],
            "expected": False,
            "meaning": "the current repo still does not export an explicit admissibility-branch discharge",
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
            "step": "N145",
            "status": "N145_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ADMISSIBILITY_BRANCH_STATE",
            "scope": "current_admissibility_branch_discharge_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N145",
            "status": "N145_DISCHARGED_CURRENT_ADMISSIBILITY_BRANCH_OBSTRUCTION_THEOREM_FOR_FUTURE_CONSTRUCTED_SOURCE_OBJECT_NO_FALSE_PASS",
            "scope": "current_admissibility_branch_discharge_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "explicit_admissibility_branch_discharge_exported": False,
                "remaining_open_branches": p134["target_state"]["remaining_open_branches"],
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "future_derivation_of_admissible_E_orient_from_a_future_new_source_object",
                "future_completion_of_B_sel_R_sel_O_sel_after_new_source_object_construction",
            ],
            "hard_limits": [
                "explicit_admissibility_branch_discharge_not_yet_exported",
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
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    print(OUT)


if __name__ == "__main__":
    main()
