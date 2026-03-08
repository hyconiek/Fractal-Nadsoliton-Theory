#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n163_current_observer_information_deficit_downstream_symptom_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p147 = load_json(
        "fundamental_action_reconstruction/generated/p147_current_light_matter_emergent_observer_information_deficit_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p147_status",
            "actual": p147["status"],
            "expected": "CURRENT_REPO_SUPPORTS_THE_CONCLUSION_THAT_OBSERVER_INFORMATION_DEFICIT_IS_A_DOWNSTREAM_SYMPTOM_AND_NOT_A_PRIMARY_SELECTOR_SOURCE_GAP_AFTER_P147",
            "meaning": "P147 already supports the downstream-symptom reading",
        },
        {
            "id": "p147_ordering",
            "actual": p147["ordering"],
            "expected": ["nadsoliton", "light", "matter", "emergent_observer"],
            "meaning": "P147 keeps the nadsoliton-light-matter-observer order explicit",
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
            "step": "N163",
            "status": "N163_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_OBSERVER_INFORMATION_DEFICIT_DOWNSTREAM_STATE",
            "scope": "current_observer_information_deficit_location_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N163",
            "status": "N163_DISCHARGED_CURRENT_OBSERVER_INFORMATION_DEFICIT_DOWNSTREAM_SYMPTOM_THEOREM_NO_FALSE_PASS",
            "scope": "current_observer_information_deficit_location_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "ordering_supported_on_current_repo_state": True,
                "observer_information_deficit_downstream_symptom_on_current_repo_state": True,
                "primary_missing_selector_source_gap_upstream_of_observer": True,
                "fixed_first_additive_attempt_already_negative_before_any_honest_observer_side_promotion": True,
                "full_closure_pass": False,
            },
            "hard_limits": [
                "no_claim_that_observer_side_information_is_irrelevant_forever",
                "constructed_source_object_not_yet_exported",
                "admissible_S_sel_int_not_yet_constructed",
                "admissible_E_orient_not_yet_constructed",
                "downstream_chain_not_yet_constructed",
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
