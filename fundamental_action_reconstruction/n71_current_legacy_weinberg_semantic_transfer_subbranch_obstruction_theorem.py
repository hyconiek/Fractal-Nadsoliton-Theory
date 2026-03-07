#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
OUT = (
    ROOT
    / "generated"
    / "n71_current_legacy_weinberg_semantic_transfer_subbranch_obstruction_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p68 = load_json(
        "fundamental_action_reconstruction/generated/p68_legacy_weinberg_semantic_transfer_subbranch_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p68_probe_negative",
            "actual": p68["status"],
            "expected": "CURRENT_REPO_EXPORTS_NEITHER_TEXTUAL_SUCCESSOR_NOR_LINEAGE_UPGRADE_TRANSFER_VERDICT_FOR_THE_LEGACY_WEINBERG_ROLE_AFTER_P68",
            "meaning": "P68 confirms that both narrowed semantic-transfer sub-branches remain absent",
        },
        {
            "id": "textual_successor_verdict_absent",
            "actual": p68["textual_successor_verdict_present"],
            "expected": False,
            "meaning": "textual retained-successor verdict is absent",
        },
        {
            "id": "lineage_upgrade_verdict_absent",
            "actual": p68["lineage_upgrade_verdict_present"],
            "expected": False,
            "meaning": "lineage-upgrade verdict is absent",
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
            "step": "N71",
            "status": "N71_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_LEGACY_WEINBERG_SEMANTIC_TRANSFER_SUBBRANCH_STATE",
            "scope": "current_legacy_weinberg_semantic_transfer_subbranch_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N71",
            "status": "N71_DISCHARGED_CURRENT_LEGACY_WEINBERG_SEMANTIC_TRANSFER_SUBBRANCH_OBSTRUCTION_THEOREM_NO_FALSE_PASS",
            "scope": "current_legacy_weinberg_semantic_transfer_subbranch_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "textual_successor_verdict_present": False,
                "lineage_upgrade_verdict_present": False,
                "retained_semantic_transfer_branch_fully_discharged": False,
                "full_closure_pass": False,
            },
            "remaining_open_branch": [
                "explicit_textual_retained_successor_verdict_identifying_sin2_theta_w_mz_as_the_retained_strict_side_successor_of_the_legacy_weinberg_angle_role",
                "explicit_lineage_upgrade_verdict_elevating_the_qw2093_alpha_geo_touchpoint_into_retained_strict_side_weinberg_role_transfer",
            ],
            "hard_limits": [
                "no_proof_that_either_subbranch_is_impossible_forever",
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
