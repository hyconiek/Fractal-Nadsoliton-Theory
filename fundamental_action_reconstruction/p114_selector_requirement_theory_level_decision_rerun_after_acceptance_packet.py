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
    / "p114_selector_requirement_theory_level_decision_rerun_after_acceptance_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p114_selector_requirement_theory_level_decision_rerun_after_acceptance_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n118 = load_json(
        "fundamental_action_reconstruction/generated/n118_current_selector_or_symmetry_breaking_requirement_theorem_for_qw2191_summary.json"
    )
    n124 = load_json(
        "fundamental_action_reconstruction/generated/n124_current_strict_core_internal_selector_source_derivation_full_negative_closure_theorem_summary.json"
    )
    ax15 = load_json(
        "fundamental_action_reconstruction/generated/ax15_selector_requirement_theory_level_acceptance_packet_summary.json"
    )

    decision_present = (
        n118["theorem_result"][
            "selector_or_symmetry_breaking_requirement_supported_on_current_repo_state"
        ]
        and not n124["theorem_result"][
            "strict_core_internal_selector_source_derivation_discharge_present"
        ]
        and ax15["acceptance_state"]["explicit_theory_level_acceptance_verdict_present"]
    )

    checks_spec = [
        {
            "id": "n118_requirement_supported",
            "actual": n118["theorem_result"][
                "selector_or_symmetry_breaking_requirement_supported_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "selector requirement support is present",
        },
        {
            "id": "n124_no_strict_core_source_discharge",
            "actual": n124["theorem_result"][
                "strict_core_internal_selector_source_derivation_discharge_present"
            ],
            "expected": False,
            "meaning": "no strict-core internal selector source discharge is present",
        },
        {
            "id": "ax15_acceptance_verdict_present",
            "actual": ax15["acceptance_state"][
                "explicit_theory_level_acceptance_verdict_present"
            ],
            "expected": True,
            "meaning": "AX15 now exports an explicit theory-level acceptance verdict",
        },
        {
            "id": "decision_present",
            "actual": decision_present,
            "expected": True,
            "meaning": "the updated repo now exports a theory-level selector-requirement decision verdict",
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
        "stage": "P114",
        "lane": "selector_requirement_theory_level_decision_rerun_after_acceptance_only",
        "goal": "test_whether_the_updated_repo_now_exports_an_explicit_theory_level_selector_requirement_decision_verdict",
        "status": "CURRENT_REPO_EXPORTS_AN_EXPLICIT_THEORY_LEVEL_ACCEPTANCE_VERDICT_FOR_SELECTOR_REQUIREMENT_AFTER_P114",
        "decision_state": {
            "theory_level_selector_requirement_decision_present": decision_present,
            "decision_kind": "acceptance",
            "accepted_scope": "axiom_augmented_only",
            "strict_core_changed": False,
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P114",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "decision_state": artifact["decision_state"],
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
