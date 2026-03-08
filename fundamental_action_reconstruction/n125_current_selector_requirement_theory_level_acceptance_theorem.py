#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n125_current_selector_requirement_theory_level_acceptance_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p114 = load_json(
        "fundamental_action_reconstruction/generated/p114_selector_requirement_theory_level_decision_rerun_after_acceptance_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "p114_decision_present",
            "actual": p114["decision_state"][
                "theory_level_selector_requirement_decision_present"
            ],
            "expected": True,
            "meaning": "P114 already proves that a theory-level selector-requirement decision is now present",
        },
        {
            "id": "decision_kind_acceptance",
            "actual": p114["decision_state"]["decision_kind"],
            "expected": "acceptance",
            "meaning": "the decision is acceptance, not deferral",
        },
        {
            "id": "accepted_scope_axiom_augmented_only",
            "actual": p114["decision_state"]["accepted_scope"],
            "expected": "axiom_augmented_only",
            "meaning": "the accepted scope stays outside current strict core",
        },
        {
            "id": "strict_core_unchanged",
            "actual": p114["decision_state"]["strict_core_changed"],
            "expected": False,
            "meaning": "strict core remains unchanged",
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
            "step": "N125",
            "status": "N125_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_SELECTOR_REQUIREMENT_ACCEPTANCE_STATE",
            "scope": "updated_repo_selector_requirement_decision_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N125",
            "status": "N125_DISCHARGED_CURRENT_SELECTOR_REQUIREMENT_THEORY_LEVEL_ACCEPTANCE_THEOREM_NO_FALSE_PASS",
            "scope": "updated_repo_selector_requirement_decision_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "selector_requirement_accepted_at_theory_level": True,
                "accepted_scope": "axiom_augmented_only",
                "strict_core_changed": False,
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "future_new_strict_core_internal_selector_source_object",
            ],
            "hard_limits": [
                "no_strict_core_selector_closure",
                "no_QW2191_discharge",
                "no_legacy_to_strict_bridge",
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
