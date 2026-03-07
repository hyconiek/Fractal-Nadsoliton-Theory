#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent


def load_json(repo_relative_path: str) -> Any:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def get_by_path(payload: Any, path: tuple[str, ...]) -> Any:
    current = payload
    for key in path:
        if not isinstance(current, dict) or key not in current:
            raise KeyError(".".join(path))
        current = current[key]
    return current


def main() -> None:
    sources = {
        "QW2191": load_json("material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"),
        "B2": load_json("fundamental_action_reconstruction/generated/b2_internal_orientation_datum_source_audit_summary.json"),
        "H30": load_json("fundamental_action_reconstruction/generated/h30_kernel_invariant_psi0_anchor_audit.json"),
        "H31": load_json("fundamental_action_reconstruction/generated/h31_psi0_to_pair1_reduction_audit.json"),
        "H33": load_json("fundamental_action_reconstruction/generated/h33_pair1_selector_target_justification_audit.json"),
        "H34": load_json("fundamental_action_reconstruction/generated/h34_basis_covariance_target_independence_audit_summary.json"),
        "H35": load_json("fundamental_action_reconstruction/generated/h35_pair1_axis_selection_audit_summary.json"),
        "H36": load_json("fundamental_action_reconstruction/generated/h36_directed_axis_orientation_audit.json"),
        "H37": load_json("fundamental_action_reconstruction/generated/h37_sign_distinction_state_audit.json"),
        "H38": load_json("fundamental_action_reconstruction/generated/h38_projective_selector_state_audit.json"),
        "H42": load_json("fundamental_action_reconstruction/generated/h42_c_based_retardation_operator_on_pair1_audit.json"),
        "P1": load_json("fundamental_action_reconstruction/generated/pair1_operator_probe_report.json"),
    }

    checks_spec = [
        {
            "id": "qw2191_kernel_uniqueness_obstructed",
            "source": "QW2191",
            "path": ("flags", "full_uniqueness_from_kernel_alone_obstructed"),
            "expected": True,
            "meaning": "kernel alone is theorem-level obstructed from full uniqueness",
        },
        {
            "id": "qw2191_extra_symmetry_breaking_required",
            "source": "QW2191",
            "path": ("flags", "obstruction_requires_explicit_symmetry_breaking_postulate"),
            "expected": True,
            "meaning": "the obstruction requires extra symmetry breaking",
        },
        {
            "id": "b2_no_internal_selector_derivation",
            "source": "B2",
            "path": ("b2", "strict_internal_selector_derivations_found"),
            "expected": 0,
            "meaning": "no internal orientation datum is derived in strict core",
        },
        {
            "id": "h30_psi0_not_selector_datum",
            "source": "H30",
            "path": ("classification", "strict_core_selector_export"),
            "expected": False,
            "meaning": "psi0 is not exported as strict-core selector datum",
        },
        {
            "id": "h31_no_selector_reduction",
            "source": "H31",
            "path": ("classification", "strict_core_selector_reduction_present"),
            "expected": False,
            "meaning": "psi0->pair1 is not a strict-core selector reduction",
        },
        {
            "id": "h33_pair1_only_local_chart",
            "source": "H33",
            "path": ("result",),
            "expected": "pair1_is_available_as_a_deterministic_local_chart_for_the_primary_psi0_lane_but_not_yet_justified_as_a_uniquely_selector_relevant_target",
            "meaning": "pair1 is only a local chart for the psi0 lane",
        },
        {
            "id": "h34_no_chart_independent_reduction",
            "source": "H34",
            "path": ("result",),
            "expected": "strict_core_supports_only_local_chart_embeddings_for_psi0_and_not_a_basis_covariant_or_target_independent_selector_reduction",
            "meaning": "no basis-covariant / target-independent reduction exists",
        },
        {
            "id": "h35_no_physical_axis",
            "source": "H35",
            "path": ("result",),
            "expected": "strict_core_supports_only_a_coordinate_level_direction_u_psi0_pair1_inside_pair1_and_not_a_strict_physical_axis_selection",
            "meaning": "no strict physical axis is selected",
        },
        {
            "id": "h36_no_directed_orientation",
            "source": "H36",
            "path": ("result",),
            "expected": "strict_core_supports_only_a_coordinate_level_undirected_axis_representative_u_psi0_pair1_inside_pair1_and_not_a_strict_directed_orientation_selection",
            "meaning": "no strict directed orientation is selected",
        },
        {
            "id": "h37_no_sign_sensitive_state",
            "source": "H37",
            "path": ("result",),
            "expected": "strict_core_contains_no_sign_sensitive_state_object_or_observable_on_pair1_and_therefore_does_not_distinguish_u_from_minus_u_as_physically_different_selector_states",
            "meaning": "no sign-sensitive selector state exists",
        },
        {
            "id": "h38_only_projective_state",
            "source": "H38",
            "path": ("result",),
            "expected": "strict_core_supports_at_most_a_local_projective_or_ray_level_selector_representative_on_pair1_and_does_not_furnish_a_physically_individuated_directed_selector_state",
            "meaning": "only a projective/ray-level representative exists",
        },
        {
            "id": "h42_psi0_plus_c_not_strict_core",
            "source": "H42",
            "path": ("hard_limits", "psi0_plus_c_is_strict_core"),
            "expected": False,
            "meaning": "the first nontrivial c-based split is not strict core",
        },
        {
            "id": "p1_split_is_anchor_imported",
            "source": "P1",
            "path": ("selected_operator", "classifier"),
            "expected": "ANCHOR_IMPORTED_SPLIT",
            "meaning": "the first computed nontrivial split is anchor-imported",
        },
        {
            "id": "p1_split_not_promoted_to_strict_core",
            "source": "P1",
            "path": ("strict_core_promotion",),
            "expected": False,
            "meaning": "the computed split is not promoted to strict core",
        },
    ]

    checks: list[dict[str, Any]] = []
    mismatches: list[str] = []
    for item in checks_spec:
        try:
            actual = get_by_path(sources[item["source"]], item["path"])
        except KeyError as exc:
            actual = f"MISSING:{exc}"
        ok = actual == item["expected"]
        checks.append(
            {
                "id": item["id"],
                "source": item["source"],
                "path": ".".join(item["path"]),
                "actual": actual,
                "expected": item["expected"],
                "pass": ok,
                "meaning": item["meaning"],
            }
        )
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        summary = {
            "step": "N5",
            "status": "N5_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FRONTIER",
            "goal": "Check whether the current strict-core psi0 route is structurally obstructed from selector closure.",
            "scope": "current_strict_core_psi0_route_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {
                "discharged": False,
                "reason": "the expected route-specific obstruction frontier has changed or a required evidence object is missing",
            },
            "required_next_step": "REVIEW_CHANGED_FRONTIER_BEFORE_CLAIMING_PSI0_ROUTE_OBSTRUCTION",
        }
    else:
        summary = {
            "step": "N5",
            "status": "N5_DISCHARGED_CURRENT_STRICT_CORE_PSI0_ROUTE_OBSTRUCTION_NO_FALSE_PASS",
            "goal": "Discharge a route-specific theorem: the current strict-core psi0 lane cannot close selector generation without extra symmetry-breaking structure.",
            "scope": "current_strict_core_psi0_route_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "kernel_alone_is_sufficient": False,
                "current_internal_orientation_datum_present": False,
                "current_strict_core_psi0_route_closes_selector": False,
                "extra_symmetry_breaking_structure_required": True,
            },
            "minimal_missing_structure_classes": [
                "internal_orientation_datum",
                "chart_independent_selector_reduction",
                "strict_physical_axis_or_equivalent_orientation_object",
                "directed_or_sign_sensitive_selector_state",
                "strict_core_operator_map_to_nontrivial_pair1_split",
            ],
            "hard_limits": [
                "no global impossibility theorem",
                "no claim that future strict-core extensions are impossible",
                "no claim that a specific future selector mechanism is already known",
                "no claim that QW-2191 is discharged",
                "no claim that ToE is closed",
            ],
            "required_next_step": "ADD_ONE_CONCRETE_STRICT_CORE_SYMMETRY_BREAKING_STRUCTURE_OR_FIND_A_NEW_ARGUMENT_THAT_GLOBALIZES_BEYOND_THIS_ROUTE_WITHOUT_FALLING_BACK_TO_T12",
        }

    out = ROOT / "generated" / "n5_current_strict_core_psi0_route_obstruction_theorem_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out)


if __name__ == "__main__":
    main()
