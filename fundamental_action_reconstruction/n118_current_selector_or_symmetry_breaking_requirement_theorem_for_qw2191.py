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


def b2_findings_map(b2: dict) -> dict[str, str]:
    return {entry["object"]: entry["status"] for entry in b2["b2"]["findings"]}


def main() -> None:
    p108 = load_json(
        "fundamental_action_reconstruction/generated/p108_current_selector_symmetry_breaking_requirement_probe_summary.json"
    )
    b2 = load_json(
        "fundamental_action_reconstruction/generated/b2_internal_orientation_datum_source_audit_summary.json"
    )
    b2_map = b2_findings_map(b2)

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
            "id": "b2_internal_orientation_datum_axis_only_found",
            "actual": b2_map.get("strict internal orientation datum"),
            "expected": "found_axis_only_residual_z2",
            "meaning": "B2 confirms an axis-only strict-core internal orientation datum exists (O(2)->Z2 cut); residual Z2 sign remains",
        },
        {
            "id": "b2_kernel_invariant_selecting_one_o2_point_missing",
            "actual": b2_map.get("kernel invariant selecting one O(2) point"),
            "expected": "not_found_in_strict_core",
            "meaning": "B2 confirms no strict-core kernel-invariant selects one unique O(2) point (full uniqueness remains obstructed)",
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
                "strict_core_internal_orientation_datum_axis_only_present": True,
                "residual_z2_sign_remaining": True,
                "kernel_invariant_unique_o2_point_present": False,
                "selector_or_symmetry_breaking_requirement_supported_on_current_repo_state": True,
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "explicit_strict_core_internal_sign_sensitive_orientation_datum_or_sign_lift_discharge",
                "explicit_strict_core_kernel_invariant_selecting_one_o2_point_or_equivalent_internal_selector_source",
                "explicit_theory_level_acceptance_of_selector_or_symmetry_breaking_requirement_if_no_sign_lift_is_derived"
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
