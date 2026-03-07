#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n116_current_legacy_physical_role_package_full_negative_closure_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n83 = load_json(
        "fundamental_action_reconstruction/generated/n83_current_legacy_weinberg_full_claim_specific_negative_closure_theorem_summary.json"
    )
    n99 = load_json(
        "fundamental_action_reconstruction/generated/n99_current_legacy_fine_structure_full_claim_specific_negative_closure_theorem_summary.json"
    )
    n115 = load_json(
        "fundamental_action_reconstruction/generated/n115_current_legacy_gravity_hierarchy_full_claim_specific_negative_closure_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n83_legacy_weinberg_closed",
            "actual": n83["theorem_result"]["legacy_weinberg_claim_specific_frontier_closed_negatively_on_current_repo_state"],
            "expected": True,
            "meaning": "N83 already closes the full legacy Weinberg claim-specific frontier negatively on the current repo state",
        },
        {
            "id": "n99_legacy_fine_structure_closed",
            "actual": n99["theorem_result"]["legacy_fine_structure_claim_specific_frontier_closed_negatively_on_current_repo_state"],
            "expected": True,
            "meaning": "N99 already closes the full legacy fine-structure claim-specific frontier negatively on the current repo state",
        },
        {
            "id": "n115_legacy_gravity_hierarchy_closed",
            "actual": n115["theorem_result"]["legacy_gravity_hierarchy_claim_specific_frontier_closed_negatively_on_current_repo_state"],
            "expected": True,
            "meaning": "N115 already closes the full legacy gravity-hierarchy claim-specific frontier negatively on the current repo state",
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
            "step": "N116",
            "status": "N116_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_LEGACY_PHYSICAL_ROLE_PACKAGE_STATE",
            "scope": "current_legacy_physical_role_package_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N116",
            "status": "N116_DISCHARGED_CURRENT_LEGACY_PHYSICAL_ROLE_PACKAGE_FULL_NEGATIVE_CLOSURE_THEOREM_NO_FALSE_PASS",
            "scope": "current_legacy_physical_role_package_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "legacy_weinberg_claim_specific_frontier_closed_negatively_on_current_repo_state": True,
                "legacy_fine_structure_claim_specific_frontier_closed_negatively_on_current_repo_state": True,
                "legacy_gravity_hierarchy_claim_specific_frontier_closed_negatively_on_current_repo_state": True,
                "legacy_physical_role_package_closed_negatively_on_current_repo_state": True,
                "full_closure_pass": False,
            },
            "remaining_open_branches": [],
            "hard_limits": [
                "no_proof_that_legacy_physical_role_transfer_is_impossible_forever",
                "no_proof_that_future_strict_side_successor_semantics_cannot_exist",
                "no_rigorous_bridge_or_nonbridge_theorem_between_k_legacy_ont_and_k_strict_gate",
                "no_selector_closure",
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
