#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED / "ax15_selector_requirement_theory_level_acceptance_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "ax15_selector_requirement_theory_level_acceptance_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n118 = load_json(
        "fundamental_action_reconstruction/generated/n118_current_selector_or_symmetry_breaking_requirement_theorem_for_qw2191_summary.json"
    )
    n122 = load_json(
        "fundamental_action_reconstruction/generated/n122_current_selector_requirement_full_theory_level_decision_negative_closure_theorem_summary.json"
    )
    n124 = load_json(
        "fundamental_action_reconstruction/generated/n124_current_strict_core_internal_selector_source_derivation_full_negative_closure_theorem_summary.json"
    )

    requirement_supported = n118["theorem_result"][
        "selector_or_symmetry_breaking_requirement_supported_on_current_repo_state"
    ]
    no_theory_level_decision_yet = n122["theorem_result"][
        "full_theory_level_decision_frontier_closed_negatively_on_current_repo_state"
    ]
    no_strict_core_source_discharge = (
        n124["theorem_result"][
            "strict_core_internal_selector_source_derivation_discharge_present"
        ]
        is False
    )

    checks_spec = [
        {
            "id": "n118_requirement_supported",
            "actual": requirement_supported,
            "expected": True,
            "meaning": "N118 already supports the selector requirement after QW-2191",
        },
        {
            "id": "n122_no_decision_yet",
            "actual": no_theory_level_decision_yet,
            "expected": True,
            "meaning": "N122 already shows that no theory-level decision verdict had yet been exported",
        },
        {
            "id": "n124_no_strict_core_source_discharge",
            "actual": no_strict_core_source_discharge,
            "expected": True,
            "meaning": "N124 already shows that no strict-core internal selector source discharge is present",
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
        "stage": "AX15",
        "lane": "selector_requirement_theory_level_acceptance_axiom_augmented_only",
        "goal": "add_an_explicit_theory_level_acceptance_verdict_for_the_selector_requirement_without_claiming_strict_core_derivation",
        "status": "AX15_EXECUTED_SELECTOR_REQUIREMENT_THEORY_LEVEL_ACCEPTANCE_PACKET_NO_FALSE_PASS",
        "acceptance_state": {
            "selector_requirement_supported": requirement_supported,
            "strict_core_internal_selector_source_derivation_discharge_present": False,
            "explicit_theory_level_acceptance_verdict_present": True,
            "accepted_scope": "axiom_augmented_only",
            "strict_core_changed": False,
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "AX15",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "acceptance_state": artifact["acceptance_state"],
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
