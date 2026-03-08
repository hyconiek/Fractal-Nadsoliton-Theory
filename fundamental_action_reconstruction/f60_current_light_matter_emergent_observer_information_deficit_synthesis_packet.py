#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "f60_current_light_matter_emergent_observer_information_deficit_synthesis_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    r3 = load_json(
        "fundamental_action_reconstruction/generated/r3_minimal_internal_light_propagator_packet_for_kobs_summary.json"
    )
    r5 = load_json(
        "fundamental_action_reconstruction/generated/r5_minimal_light_to_matter_response_packet_for_kobs_summary.json"
    )
    r6 = load_json(
        "fundamental_action_reconstruction/generated/r6_minimal_observer_readout_packet_for_kobs_summary.json"
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
            "id": "r3_light_packet_present",
            "actual": r3["result"],
            "expected": "explicit_internal_light_propagator_packet_present_but_eigenchannel_only_and_not_yet_factorized_from_existing_kernel_feedback",
            "meaning": "R3 exports an explicit light packet before matter and observer",
        },
        {
            "id": "r5_matter_packet_present",
            "actual": r5["result"],
            "expected": "explicit_light_to_matter_response_packet_present_but_current_pair_scoped_and_not_yet_factorized_from_existing_kernel_feedback",
            "meaning": "R5 exports an explicit light-to-matter packet before observer",
        },
        {
            "id": "r6_observer_packet_present",
            "actual": r6["result"],
            "expected": "explicit_observer_readout_packet_present_but_current_pair_scoped_and_not_yet_factorized_from_existing_kernel_feedback",
            "meaning": "R6 exports only a downstream observer-readout packet",
        },
        {
            "id": "n118_internal_source_missing",
            "actual": n118["theorem_result"]["strict_core_internal_selector_source_present"],
            "expected": False,
            "meaning": "N118 keeps the primary selector-source gap upstream of observer",
        },
        {
            "id": "n149_constructive_frontier_exhausted",
            "actual": n149["theorem_result"][
                "constructive_selector_frontier_exhausted_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N149 already exhausts the present constructive selector frontier",
        },
        {
            "id": "n162_fixed_additive_attempt_negative",
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
            "stage": "F60",
            "lane": "current_light_matter_observer_synthesis_only",
            "status": "F60_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_LIGHT_MATTER_OBSERVER_SYNTHESIS_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "F60",
            "lane": "current_light_matter_observer_information_deficit_synthesis_only",
            "goal": "freeze_the_current_repo_state_reading_that_observer_information_deficit_is_downstream_of_light_and_matter_and_not_the_primary_missing_selector_source",
            "status": "F60_EXECUTED_CURRENT_LIGHT_MATTER_EMERGENT_OBSERVER_INFORMATION_DEFICIT_SYNTHESIS_PACKET_NO_FALSE_PASS",
            "ordering": [
                "nadsoliton",
                "light",
                "matter",
                "emergent_observer",
            ],
            "classification": {
                "observer_information_deficit": "downstream_symptom",
                "primary_missing_selector_source_gap": "upstream_of_observer",
            },
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
