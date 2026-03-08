#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED
    / "p111_strict_side_theory_level_deferral_verdict_probe_for_selector_requirement.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p111_strict_side_theory_level_deferral_verdict_probe_for_selector_requirement_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p109 = load_json(
        "fundamental_action_reconstruction/generated/p109_current_selector_requirement_theory_acceptance_probe_summary.json"
    )
    n120 = load_json(
        "fundamental_action_reconstruction/generated/n120_current_selector_requirement_acceptance_branch_full_negative_closure_theorem_summary.json"
    )

    deferral_missing_object = (
        "explicit_theory_level_deferral_verdict_keeping_the_selector_or_symmetry_breaking_requirement_as_an_active_boundary_while_internal_source_search_continues"
    )

    checks_spec = [
        {
            "id": "p109_decision_verdict_absent",
            "actual": p109["status"],
            "expected": "CURRENT_REPO_EXPORTS_SELECTOR_REQUIREMENT_SUPPORT_AND_DECISION_INPUTS_BUT_NO_EXPLICIT_THEORY_LEVEL_DECISION_VERDICT_AFTER_P109",
            "meaning": "P109 already proves that no explicit theory-level decision verdict is currently exported",
        },
        {
            "id": "n120_acceptance_branch_closed_negatively",
            "actual": n120["theorem_result"][
                "acceptance_branch_closed_negatively_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N120 already closes the acceptance branch negatively, leaving deferral as the only remaining decision branch",
        },
        {
            "id": "deferral_missing_object_present",
            "actual": deferral_missing_object in p109["remaining_missing_objects"],
            "expected": True,
            "meaning": "the deferral branch remains explicitly missing after P109",
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
        "stage": "P111",
        "lane": "selector_requirement_deferral_branch_current_repo_state_only",
        "goal": "test_whether_the_current_repo_already_exports_an_explicit_theory_level_deferral_verdict_for_the_selector_or_symmetry_breaking_requirement",
        "status": "CURRENT_REPO_DOES_NOT_EXPORT_AN_EXPLICIT_THEORY_LEVEL_DEFERRAL_VERDICT_FOR_SELECTOR_OR_SYMMETRY_BREAKING_REQUIREMENT_AFTER_P111",
        "reason": "P109 already keeps the full decision question open, N120 already closes the acceptance branch negatively, and the deferral branch remains an explicit missing object; therefore the current repo still exports no theory-level deferral verdict for the selector or symmetry-breaking requirement",
        "deferral_branch_state": {
            "theory_level_deferral_verdict_present": False,
        },
        "remaining_missing_objects": [
            deferral_missing_object,
            "explicit_strict_core_internal_selector_source_derivation_discharge",
            "explicit_legacy_to_strict_kernel_bridge_or_nonbridge_theorem_with_package_level_scope",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P111",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "deferral_branch_state": artifact["deferral_branch_state"],
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
