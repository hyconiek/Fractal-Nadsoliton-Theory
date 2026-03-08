#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "p147_current_light_matter_emergent_observer_information_deficit_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f60 = load_json(
        "fundamental_action_reconstruction/generated/f60_current_light_matter_emergent_observer_information_deficit_synthesis_packet_summary.json"
    )
    n118 = load_json(
        "fundamental_action_reconstruction/generated/n118_current_selector_or_symmetry_breaking_requirement_theorem_for_qw2191_summary.json"
    )
    n149 = load_json(
        "fundamental_action_reconstruction/generated/n149_current_repo_constructive_selector_frontier_exhaustion_theorem_summary.json"
    )
    n162 = load_json(
        "fundamental_action_reconstruction/generated/n162_current_fixed_first_additive_construction_attempt_full_negative_closure_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "f60_observer_deficit_downstream",
            "actual": f60["classification"]["observer_information_deficit"],
            "expected": "downstream_symptom",
            "meaning": "F60 already classifies observer information deficit as downstream",
        },
        {
            "id": "f60_primary_gap_upstream",
            "actual": f60["classification"]["primary_missing_selector_source_gap"],
            "expected": "upstream_of_observer",
            "meaning": "F60 already places the primary selector gap upstream of observer",
        },
        {
            "id": "n118_internal_source_missing",
            "actual": n118["theorem_result"]["strict_core_internal_selector_source_present"],
            "expected": False,
            "meaning": "N118 already excludes a current strict-core internal selector source",
        },
        {
            "id": "n149_constructive_frontier_exhausted",
            "actual": n149["theorem_result"][
                "constructive_selector_frontier_exhausted_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N149 already exhausts the present constructive frontier",
        },
        {
            "id": "n162_fixed_attempt_negative",
            "actual": n162["theorem_result"][
                "fixed_first_additive_construction_attempt_closed_negatively_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N162 already keeps the fixed first additive attempt negative",
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
            "stage": "P147",
            "lane": "current_light_matter_observer_information_deficit_probe_only",
            "status": "P147_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_LIGHT_MATTER_OBSERVER_INFORMATION_DEFICIT_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P147",
            "lane": "current_light_matter_observer_information_deficit_probe_only",
            "goal": "test_whether_the_current_repo_supports_classifying_observer_information_deficit_as_a_downstream_symptom_rather_than_as_the_primary_missing_selector_source_gap",
            "status": "CURRENT_REPO_SUPPORTS_THE_CONCLUSION_THAT_OBSERVER_INFORMATION_DEFICIT_IS_A_DOWNSTREAM_SYMPTOM_AND_NOT_A_PRIMARY_SELECTOR_SOURCE_GAP_AFTER_P147",
            "ordering": [
                "nadsoliton",
                "light",
                "matter",
                "emergent_observer",
            ],
            "checks": checks,
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    print(OUT)


if __name__ == "__main__":
    main()
