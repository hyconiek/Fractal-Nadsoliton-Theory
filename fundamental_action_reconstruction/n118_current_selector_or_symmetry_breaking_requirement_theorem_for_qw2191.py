#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n118_current_selector_or_symmetry_breaking_requirement_theorem_for_qw2191_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p108 = load_json(
        "fundamental_action_reconstruction/generated/p108_current_selector_symmetry_breaking_requirement_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p108_supports_requirement_conclusion",
            "actual": p108["status"],
            "expected": "CURRENT_REPO_SUPPORTS_THE_SELECTOR_OR_SYMMETRY_BREAKING_REQUIREMENT_CONCLUSION_FOR_THE_QW2191_UNIQUENESS_FRONTIER_AFTER_P108",
            "meaning": "P108 already supports the selector/symmetry-breaking requirement conclusion on the current repo state",
        },
        {
            "id": "q2191_kernel_alone_obstructed",
            "actual": p108["q2191_state"]["kernel_alone_obstructed"],
            "expected": True,
            "meaning": "QW-2191 keeps kernel-alone uniqueness obstructed",
        },
        {
            "id": "b2_internal_source_missing",
            "actual": p108["strict_core_source_state"]["internal_orientation_datum_status"],
            "expected": "not_found_in_strict_core",
            "meaning": "B2 keeps the strict-core internal selector source absent",
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
            "step": "N118",
            "status": "N118_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_QW2191_SELECTOR_REQUIREMENT_STATE",
            "scope": "current_qw2191_selector_requirement_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N118",
            "status": "N118_DISCHARGED_CURRENT_SELECTOR_OR_SYMMETRY_BREAKING_REQUIREMENT_THEOREM_FOR_QW2191_NO_FALSE_PASS",
            "scope": "current_qw2191_selector_requirement_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "kernel_alone_uniqueness_obstructed": True,
                "axiom_augmented_selector_route_available": True,
                "robust_selector_family_available": True,
                "strict_core_internal_selector_source_present": False,
                "selector_or_symmetry_breaking_requirement_supported_on_current_repo_state": True,
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "explicit_strict_core_internal_selector_source_derivation_discharge",
                "explicit_theory_level_acceptance_of_selector_or_symmetry_breaking_requirement_if_no_internal_source_is_derived",
                "explicit_legacy_to_strict_kernel_bridge_or_nonbridge_theorem_with_package_level_scope"
            ],
            "hard_limits": [
                "no_strict_core_selector_closure",
                "no_QW2191_discharge",
                "no_proof_that_no_future_internal_selector_source_can_exist",
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
