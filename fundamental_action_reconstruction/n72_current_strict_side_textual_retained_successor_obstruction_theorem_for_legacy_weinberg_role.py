#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
OUT = (
    ROOT
    / "generated"
    / "n72_current_strict_side_textual_retained_successor_obstruction_theorem_for_legacy_weinberg_role_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p69 = load_json(
        "fundamental_action_reconstruction/generated/p69_strict_side_textual_retained_successor_probe_for_legacy_weinberg_role_summary.json"
    )

    checks_spec = [
        {
            "id": "p69_probe_negative",
            "actual": p69["status"],
            "expected": "CURRENT_STRICT_SIDE_SOURCES_DO_NOT_EXPORT_AN_EXPLICIT_TEXTUAL_RETAINED_SUCCESSOR_VERDICT_FOR_THE_LEGACY_WEINBERG_ROLE_AFTER_P69",
            "meaning": "P69 confirms that the current strict-side source set exports no explicit textual retained-successor verdict",
        },
        {
            "id": "textual_successor_verdict_absent",
            "actual": p69["textual_successor_verdict_present"],
            "expected": False,
            "meaning": "textual retained-successor verdict is absent on the current repo state",
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
            "step": "N72",
            "status": "N72_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_TEXTUAL_RETAINED_SUCCESSOR_STATE",
            "scope": "current_textual_retained_successor_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N72",
            "status": "N72_DISCHARGED_CURRENT_STRICT_SIDE_TEXTUAL_RETAINED_SUCCESSOR_OBSTRUCTION_THEOREM_FOR_LEGACY_WEINBERG_ROLE_NO_FALSE_PASS",
            "scope": "current_textual_retained_successor_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "textual_successor_verdict_present": False,
                "textual_successor_path_closed_on_current_repo_state": True,
                "full_closure_pass": False,
            },
            "remaining_open_branch": [
                "explicit_lineage_upgrade_verdict_elevating_the_qw2093_alpha_geo_touchpoint_into_retained_strict_side_weinberg_role_transfer"
            ],
            "hard_limits": [
                "no_proof_that_the_lineage_upgrade_subbranch_is_absent",
                "no_proof_that_the_replaced_branch_is_solved",
                "no_proof_that_textual_successor_semantics_are_impossible_forever",
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
