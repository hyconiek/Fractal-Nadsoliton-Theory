#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
OUT = (
    ROOT
    / "generated"
    / "n88_current_legacy_fine_structure_semantic_transfer_subbranch_obstruction_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p83 = load_json(
        "fundamental_action_reconstruction/generated/p83_legacy_fine_structure_semantic_transfer_subbranch_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p83_status_matches_expected_boundary",
            "actual": p83["status"],
            "expected": "CURRENT_REPO_EXPORTS_NEITHER_TEXTUAL_SUCCESSOR_NOR_OBJECT_LINEAGE_UPGRADE_TRANSFER_VERDICT_FOR_THE_LEGACY_FINE_STRUCTURE_ROLE_AFTER_P83",
            "meaning": "P83 confirms that the current repo exports neither retained semantic-transfer sub-branch for the legacy fine-structure role",
        },
        {
            "id": "textual_retained_successor_verdict_absent",
            "actual": p83["textual_retained_successor_verdict_present"],
            "expected": False,
            "meaning": "the explicit textual retained-successor verdict is still absent",
        },
        {
            "id": "object_lineage_upgrade_verdict_absent",
            "actual": p83["object_lineage_upgrade_verdict_present"],
            "expected": False,
            "meaning": "the explicit object-lineage-upgrade verdict is still absent",
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
            "step": "N88",
            "status": "N88_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_LEGACY_FINE_STRUCTURE_SEMANTIC_TRANSFER_STATE",
            "scope": "current_legacy_fine_structure_retained_semantic_transfer_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N88",
            "status": "N88_DISCHARGED_CURRENT_LEGACY_FINE_STRUCTURE_SEMANTIC_TRANSFER_SUBBRANCH_OBSTRUCTION_THEOREM_NO_FALSE_PASS",
            "scope": "current_legacy_fine_structure_retained_semantic_transfer_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "textual_retained_successor_verdict_present": False,
                "object_lineage_upgrade_verdict_present": False,
                "retained_branch_fully_discharged": False,
                "full_closure_pass": False,
            },
            "remaining_open_branch": [
                "explicit_textual_retained_successor_verdict_identifying_alpha_em_inv_mz_as_the_retained_strict_side_successor_of_the_legacy_fine_structure_role",
                "explicit_object_lineage_upgrade_verdict_elevating_the_existing_alpha_em_inv_mz_candidate_chain_into_retained_strict_side_fine_structure_role_transfer",
            ],
            "hard_limits": [
                "no_proof_that_either_semantic_transfer_subbranch_is_impossible_forever",
                "no_proof_that_the_replaced_branch_is_solved",
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
