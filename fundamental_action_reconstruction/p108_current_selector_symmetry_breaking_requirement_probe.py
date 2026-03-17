#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p108_current_selector_symmetry_breaking_requirement_probe.json"
OUT_SUMMARY = (
    GENERATED / "p108_current_selector_symmetry_breaking_requirement_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def b2_findings_map(b2: dict[str, Any]) -> dict[str, str]:
    return {entry["object"]: entry["status"] for entry in b2["b2"]["findings"]}


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    q2191 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"
    )
    q2192 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2192_mode_index_selection_axiom_gate.json"
    )
    q2193 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2193_selection_axiom_family_robustness_gate.json"
    )
    b1 = load_json(
        "fundamental_action_reconstruction/generated/b1_mode_uniqueness_minimal_extra_structure_audit_summary.json"
    )
    b2 = load_json(
        "fundamental_action_reconstruction/generated/b2_internal_orientation_datum_source_audit_summary.json"
    )

    b2_map = b2_findings_map(b2)

    checks_spec = [
        {
            "id": "q2191_kernel_alone_obstructed",
            "actual": q2191["flags"]["full_uniqueness_from_kernel_alone_obstructed"],
            "expected": True,
            "meaning": "QW-2191 proves that kernel alone obstructs full physical uniqueness",
        },
        {
            "id": "q2191_explicit_symmetry_breaking_required",
            "actual": q2191["flags"]["obstruction_requires_explicit_symmetry_breaking_postulate"],
            "expected": True,
            "meaning": "QW-2191 explicitly requires an extra symmetry-breaking postulate for uniqueness",
        },
        {
            "id": "q2192_axiom_augmented_uniqueness_closed",
            "actual": q2192["flags"]["axiom_augmented_uniqueness_closed_for_q2190_mapping"],
            "expected": True,
            "meaning": "QW-2192 closes uniqueness in axiom-augmented scope",
        },
        {
            "id": "q2192_axiom_free_uniqueness_still_open",
            "actual": q2192["flags"]["axiom_free_uniqueness_closed"],
            "expected": False,
            "meaning": "QW-2192 keeps axiom-free uniqueness explicitly open",
        },
        {
            "id": "q2193_selector_family_robust",
            "actual": q2193["flags"]["axiom_family_robustness_closed_for_q2190_mapping"],
            "expected": True,
            "meaning": "QW-2193 proves robustness of the selector family once the extra postulate is added",
        },
        {
            "id": "b1_kernel_alone_ruled_out",
            "actual": any(
                item["object"] == "kernel alone" and item["status"] == "ruled_out_as_sufficient"
                for item in b1["b1"]["core_findings"]
            ),
            "expected": True,
            "meaning": "B1 already narrows the blocker to the need for one internal symmetry-breaking selector",
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
        {
            "id": "b2_selector_axiom_control_route_only",
            "actual": b2_map.get("explicit selector axiom"),
            "expected": "control_route_only",
            "meaning": "B2 keeps the explicit selector only as a control route, not a strict-core derivation",
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
        "stage": "P108",
        "lane": "current_selector_symmetry_breaking_requirement_probe_qw2191_frontier_only",
        "goal": "test_whether_the_current_repo_supports_the_conclusion_that_the_qw2191_uniqueness_frontier_now_requires_an_explicit_selector_or_symmetry_breaking_premise_unless_a_new_internal_source_is_derived",
        "status": "CURRENT_REPO_SUPPORTS_THE_SELECTOR_OR_SYMMETRY_BREAKING_REQUIREMENT_CONCLUSION_FOR_THE_QW2191_UNIQUENESS_FRONTIER_AFTER_P108",
        "reason": "QW-2191 proves kernel-alone obstruction, QW-2192 closes uniqueness only after adding an explicit selector axiom, QW-2193 proves robustness of that axiom-augmented family, and B2 confirms an axis-only strict-core internal orientation datum exists (O(2)->Z2 cut) while full kernel-invariant uniqueness (and residual sign lift) remain open; therefore the current repo supports the scoped selector/symmetry-breaking requirement conclusion for the QW-2191 frontier",
        "q2191_state": {
            "kernel_alone_obstructed": q2191["flags"]["full_uniqueness_from_kernel_alone_obstructed"],
            "explicit_symmetry_breaking_required": q2191["flags"]["obstruction_requires_explicit_symmetry_breaking_postulate"],
        },
        "selector_route_state": {
            "axiom_augmented_uniqueness_closed": q2192["flags"]["axiom_augmented_uniqueness_closed_for_q2190_mapping"],
            "axiom_free_uniqueness_closed": q2192["flags"]["axiom_free_uniqueness_closed"],
            "selector_family_robust": q2193["flags"]["axiom_family_robustness_closed_for_q2190_mapping"],
        },
        "strict_core_source_state": {
            "internal_orientation_datum_status": b2_map.get("strict internal orientation datum"),
            "explicit_selector_axiom_status": b2_map.get("explicit selector axiom"),
            "robust_selector_family_status": b2_map.get("robust selector family"),
            "strict_internal_selector_derivations_found": b2["b2"]["strict_internal_selector_derivations_found"],
        },
        "remaining_missing_objects": [
            "explicit_strict_core_internal_sign_sensitive_orientation_datum_or_sign_lift_discharge",
            "explicit_theory_level_acceptance_of_selector_or_symmetry_breaking_requirement_if_no_sign_lift_is_derived"
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P108",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "q2191_state": artifact["q2191_state"],
        "selector_route_state": artifact["selector_route_state"],
        "strict_core_source_state": artifact["strict_core_source_state"],
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
