#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n122_current_selector_requirement_full_theory_level_decision_negative_closure_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n120 = load_json(
        "fundamental_action_reconstruction/generated/n120_current_selector_requirement_acceptance_branch_full_negative_closure_theorem_summary.json"
    )
    n121 = load_json(
        "fundamental_action_reconstruction/generated/n121_current_selector_requirement_deferral_branch_full_negative_closure_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "acceptance_branch_closed_negatively",
            "actual": n120["theorem_result"][
                "acceptance_branch_closed_negatively_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N120 already closes the acceptance branch negatively",
        },
        {
            "id": "deferral_branch_closed_negatively",
            "actual": n121["theorem_result"][
                "deferral_branch_closed_negatively_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N121 already closes the deferral branch negatively",
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
            "step": "N122",
            "status": "N122_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_SELECTOR_REQUIREMENT_DECISION_FRONTIER_STATE",
            "scope": "current_selector_requirement_theory_level_decision_frontier_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N122",
            "status": "N122_DISCHARGED_CURRENT_SELECTOR_REQUIREMENT_FULL_THEORY_LEVEL_DECISION_NEGATIVE_CLOSURE_THEOREM_NO_FALSE_PASS",
            "scope": "current_selector_requirement_theory_level_decision_frontier_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "acceptance_branch_closed_negatively_on_current_repo_state": True,
                "deferral_branch_closed_negatively_on_current_repo_state": True,
                "full_theory_level_decision_frontier_closed_negatively_on_current_repo_state": True,
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "explicit_strict_core_internal_selector_source_derivation_discharge",
                "explicit_legacy_to_strict_kernel_bridge_or_nonbridge_theorem_with_package_level_scope",
            ],
            "hard_limits": [
                "no_theory_level_selector_decision_verdict",
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
