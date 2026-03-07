#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n119_current_selector_requirement_theory_level_decision_boundary_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p109 = load_json(
        "fundamental_action_reconstruction/generated/p109_current_selector_requirement_theory_acceptance_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p109_support_and_decision_inputs_present",
            "actual": p109["status"],
            "expected": "CURRENT_REPO_EXPORTS_SELECTOR_REQUIREMENT_SUPPORT_AND_DECISION_INPUTS_BUT_NO_EXPLICIT_THEORY_LEVEL_DECISION_VERDICT_AFTER_P109",
            "meaning": "P109 already proves that support and decision inputs are present while no explicit decision verdict is exported",
        },
        {
            "id": "acceptance_verdict_absent",
            "actual": "explicit_theory_level_acceptance_verdict_adopting_the_selector_or_symmetry_breaking_requirement_into_axiom_augmented_theory_scope_if_no_internal_source_is_derived"
            in p109["remaining_missing_objects"],
            "expected": True,
            "meaning": "the theory-level acceptance verdict remains absent",
        },
        {
            "id": "deferral_verdict_absent",
            "actual": "explicit_theory_level_deferral_verdict_keeping_the_selector_or_symmetry_breaking_requirement_as_an_active_boundary_while_internal_source_search_continues"
            in p109["remaining_missing_objects"],
            "expected": True,
            "meaning": "the theory-level deferral verdict also remains absent",
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
            "step": "N119",
            "status": "N119_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_SELECTOR_REQUIREMENT_DECISION_STATE",
            "scope": "current_selector_requirement_theory_level_decision_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N119",
            "status": "N119_DISCHARGED_CURRENT_SELECTOR_REQUIREMENT_THEORY_LEVEL_DECISION_BOUNDARY_THEOREM_NO_FALSE_PASS",
            "scope": "current_selector_requirement_theory_level_decision_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "selector_requirement_supported_on_current_repo_state": True,
                "theory_level_acceptance_verdict_present": False,
                "theory_level_deferral_verdict_present": False,
                "theory_level_decision_boundary_crossed": False,
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "explicit_theory_level_acceptance_verdict_adopting_the_selector_or_symmetry_breaking_requirement_into_axiom_augmented_theory_scope_if_no_internal_source_is_derived",
                "explicit_theory_level_deferral_verdict_keeping_the_selector_or_symmetry_breaking_requirement_as_an_active_boundary_while_internal_source_search_continues",
                "explicit_strict_core_internal_selector_source_derivation_discharge",
                "explicit_legacy_to_strict_kernel_bridge_or_nonbridge_theorem_with_package_level_scope",
            ],
            "hard_limits": [
                "no_strict_core_selector_closure",
                "no_QW2191_discharge",
                "no_accepted_theory_level_selector_requirement_verdict",
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
