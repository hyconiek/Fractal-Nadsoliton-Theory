#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n149_current_repo_constructive_selector_frontier_exhaustion_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n123 = load_json(
        "fundamental_action_reconstruction/generated/n123_current_legacy_to_strict_kernel_package_level_nonbridge_theorem_summary.json"
    )
    n125 = load_json(
        "fundamental_action_reconstruction/generated/n125_current_selector_requirement_theory_level_acceptance_theorem_summary.json"
    )
    n126 = load_json(
        "fundamental_action_reconstruction/generated/n126_current_repo_exports_no_admissible_strict_core_internal_selector_source_object_theorem_summary.json"
    )
    n148 = load_json(
        "fundamental_action_reconstruction/generated/n148_current_post_verdict_lower_branch_full_negative_closure_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n123_package_nonbridge_discharged",
            "actual": n123["theorem_result"]["package_level_nonbridge_on_current_repo_state"],
            "expected": True,
            "meaning": "N123 already discharges package-level nonbridge on the current repo state",
        },
        {
            "id": "n125_selector_requirement_accepted_outside_strict_core",
            "actual": (
                n125["theorem_result"]["selector_requirement_accepted_at_theory_level"]
                and n125["theorem_result"]["accepted_scope"] == "axiom_augmented_only"
                and n125["theorem_result"]["strict_core_changed"] is False
            ),
            "expected": True,
            "meaning": "N125 already freezes theory-level acceptance outside strict core",
        },
        {
            "id": "n126_no_admissible_source_object_present",
            "actual": n126["theorem_result"]["admissible_strict_core_internal_selector_source_object_present"],
            "expected": False,
            "meaning": "N126 already excludes all current exports from the role of admissible strict-core selector source object",
        },
        {
            "id": "n148_lower_branch_frontier_closed",
            "actual": n148["theorem_result"]["post_verdict_lower_branch_frontier_closed_negatively_on_current_repo_state"],
            "expected": True,
            "meaning": "N148 already closes the whole post-verdict lower-branch frontier negatively",
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
            "step": "N149",
            "status": "N149_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_CONSTRUCTIVE_SELECTOR_FRONTIER_STATE",
            "scope": "current_constructive_selector_frontier_as_a_whole",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N149",
            "status": "N149_DISCHARGED_CURRENT_REPO_CONSTRUCTIVE_SELECTOR_FRONTIER_EXHAUSTION_THEOREM_NO_FALSE_PASS",
            "scope": "current_constructive_selector_frontier_as_a_whole",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "package_level_nonbridge_on_current_repo_state": True,
                "selector_requirement_accepted_at_theory_level_outside_strict_core": True,
                "admissible_strict_core_internal_selector_source_object_present": False,
                "post_verdict_lower_branch_frontier_closed_negatively_on_current_repo_state": True,
                "constructive_selector_frontier_exhausted_on_current_repo_state": True,
                "only_remaining_positive_move_class": "future_genuinely_additive_new_strict_core_source_object_construction",
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "future_genuinely_additive_new_strict_core_source_object_construction"
            ],
            "hard_limits": [
                "no_proof_that_future_additive_source_object_construction_is_impossible_forever",
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
