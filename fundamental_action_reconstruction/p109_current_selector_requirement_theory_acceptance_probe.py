#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED / "p109_current_selector_requirement_theory_acceptance_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p109_current_selector_requirement_theory_acceptance_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f27 = load_json(
        "fundamental_action_reconstruction/generated/f27_current_selector_requirement_theory_acceptance_refinement_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "selector_requirement_supported",
            "actual": f27["support_state"]["selector_requirement_supported"],
            "expected": True,
            "meaning": "the selector requirement conclusion is already supported on the current repo state",
        },
        {
            "id": "strategic_priority_formalized",
            "actual": f27["support_state"]["strategic_priority_formalized"],
            "expected": True,
            "meaning": "the project already records this as a top-level strategic priority",
        },
        {
            "id": "theory_level_acceptance_verdict_present",
            "actual": f27["support_state"]["theory_level_acceptance_verdict_present"],
            "expected": False,
            "meaning": "no explicit theory-level acceptance verdict is exported",
        },
        {
            "id": "theory_level_deferral_verdict_present",
            "actual": f27["support_state"]["theory_level_deferral_verdict_present"],
            "expected": False,
            "meaning": "no explicit theory-level deferral verdict is exported",
        },
    ]

    checks: list[dict[str, Any]] = []
    for item in checks_spec:
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": item["actual"] == item["expected"],
                "meaning": item["meaning"],
            }
        )

    artifact = {
        "stage": "P109",
        "lane": "current_selector_requirement_theory_acceptance_probe_qw2191_followup_only",
        "goal": "test_whether_the_current_repo_already_exports_any_explicit_theory_level_decision_verdict_about_the_selector_or_symmetry_breaking_requirement_after_qw2191",
        "status": "CURRENT_REPO_EXPORTS_SELECTOR_REQUIREMENT_SUPPORT_AND_DECISION_INPUTS_BUT_NO_EXPLICIT_THEORY_LEVEL_DECISION_VERDICT_AFTER_P109",
        "reason": "the current repo already supports the selector requirement conclusion after QW-2191 and already elevates that question strategically, but still exports neither an explicit theory-level acceptance verdict nor an explicit theory-level deferral verdict; therefore the theory-level decision remains open",
        "support_state": f27["support_state"],
        "remaining_missing_objects": [
            "explicit_theory_level_acceptance_verdict_adopting_the_selector_or_symmetry_breaking_requirement_into_axiom_augmented_theory_scope_if_no_internal_source_is_derived",
            "explicit_theory_level_deferral_verdict_keeping_the_selector_or_symmetry_breaking_requirement_as_an_active_boundary_while_internal_source_search_continues",
            "explicit_strict_core_internal_selector_source_derivation_discharge",
            "explicit_legacy_to_strict_kernel_bridge_or_nonbridge_theorem_with_package_level_scope",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P109",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "support_state": artifact["support_state"],
        "remaining_missing_objects": artifact["remaining_missing_objects"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    print(OUT_JSON)


if __name__ == "__main__":
    main()
