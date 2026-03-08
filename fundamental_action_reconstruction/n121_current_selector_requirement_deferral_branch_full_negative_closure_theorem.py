#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n121_current_selector_requirement_deferral_branch_full_negative_closure_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p111 = load_json(
        "fundamental_action_reconstruction/generated/p111_strict_side_theory_level_deferral_verdict_probe_for_selector_requirement_summary.json"
    )

    checks_spec = [
        {
            "id": "p111_deferral_verdict_absent",
            "actual": p111["status"],
            "expected": "CURRENT_REPO_DOES_NOT_EXPORT_AN_EXPLICIT_THEORY_LEVEL_DEFERRAL_VERDICT_FOR_SELECTOR_OR_SYMMETRY_BREAKING_REQUIREMENT_AFTER_P111",
            "meaning": "P111 already proves that the deferral verdict is absent on the current repo state",
        },
        {
            "id": "deferral_branch_present",
            "actual": p111["deferral_branch_state"][
                "theory_level_deferral_verdict_present"
            ],
            "expected": False,
            "meaning": "the deferral branch is not exported on the current repo state",
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
            "step": "N121",
            "status": "N121_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_SELECTOR_REQUIREMENT_DEFERRAL_BRANCH_STATE",
            "scope": "current_selector_requirement_deferral_branch_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N121",
            "status": "N121_DISCHARGED_CURRENT_SELECTOR_REQUIREMENT_DEFERRAL_BRANCH_FULL_NEGATIVE_CLOSURE_THEOREM_NO_FALSE_PASS",
            "scope": "current_selector_requirement_deferral_branch_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "deferral_branch_closed_negatively_on_current_repo_state": True,
                "selector_requirement_supported_on_current_repo_state": True,
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "explicit_strict_core_internal_selector_source_derivation_discharge",
                "explicit_legacy_to_strict_kernel_bridge_or_nonbridge_theorem_with_package_level_scope",
            ],
            "hard_limits": [
                "no_theory_level_deferral_verdict",
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
