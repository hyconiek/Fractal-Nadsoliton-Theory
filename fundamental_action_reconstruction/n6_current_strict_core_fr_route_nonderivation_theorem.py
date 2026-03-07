#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def blocker_present(items: list[str], needle: str) -> bool:
    return needle in items


def main() -> None:
    sources = {
        "B4": load_json("fundamental_action_reconstruction/generated/b4_minimal_sigma_int_candidate_summary.json"),
        "B5": load_json("fundamental_action_reconstruction/generated/b5_sigma_int_local_stability_audit_summary.json"),
        "B6": load_json("fundamental_action_reconstruction/generated/b6_sigma_to_selector_factorized_bridge_summary.json"),
        "B7": load_json("fundamental_action_reconstruction/generated/b7_factorized_selector_mode_scaffold_compatibility_audit_summary.json"),
        "B8": load_json("fundamental_action_reconstruction/generated/b8_selector_track_anti_overclaim_audit_summary.json"),
        "T2": load_json("fundamental_action_reconstruction/generated/t2_sigma_int_to_residual_datum_bridge_theorem_spec_summary.json"),
        "C35": load_json("fundamental_action_reconstruction/generated/c35_actual_phase_source_branch_audit_summary.json"),
        "C37": load_json("fundamental_action_reconstruction/generated/c37_residual_orientation_datum_internalization_candidate_audit_summary.json"),
        "C38": load_json("fundamental_action_reconstruction/generated/c38_sigma_int_residual_datum_theorem_spec_audit_summary.json"),
    }

    b8_blockers = sources["B8"]["residual_blockers"]

    checks_spec = [
        {
            "id": "b4_candidate_exists",
            "actual": sources["B4"]["b4"]["candidate_name"],
            "expected": "sigma_int_candidate",
            "meaning": "the route starts from a candidate object",
        },
        {
            "id": "b8_no_strict_derivation_of_sigma",
            "actual": blocker_present(b8_blockers, "no_strict_derivation_of_sigma_int_candidate"),
            "expected": True,
            "meaning": "sigma_int_candidate is not strict-derived",
        },
        {
            "id": "b5_gauge_quotient_safety_open",
            "actual": sources["B5"]["b5"]["findings"][2]["status"],
            "expected": "open",
            "meaning": "full gauge-quotient safety remains open",
        },
        {
            "id": "b6_sigma_alone_to_theta_not_shown",
            "actual": sources["B6"]["findings"]["sigma_alone_selects_theta"]["status"],
            "expected": "not_shown",
            "meaning": "sigma alone does not derive theta",
        },
        {
            "id": "b7_overlay_only",
            "actual": sources["B7"]["findings"]["compatibility_with_a6_boundary"]["status"],
            "expected": "partial_control_route_only",
            "meaning": "bridge compatibility remains overlay-only",
        },
        {
            "id": "t2_map_absent",
            "actual": sources["T2"]["findings"]["strict_core_equivalence_or_export_map_present"],
            "expected": False,
            "meaning": "strict-core bridge map is absent",
        },
        {
            "id": "c37_equivalence_not_shown",
            "actual": sources["C37"]["result"]["strict_core_equivalence_bridge"],
            "expected": "not_shown",
            "meaning": "candidate-fit is not strict-core equivalence",
        },
        {
            "id": "c38_theorem_spec_absent",
            "actual": sources["C38"]["findings"]["strict_core_theorem_spec_present"],
            "expected": False,
            "meaning": "no strict-core theorem-spec exists for the bridge",
        },
        {
            "id": "c35_actual_theta_source_absent",
            "actual": sources["C35"]["result"]["strict_core_actual_phase_source"],
            "expected": "not_shown",
            "meaning": "the route does not supply actual theta-source in strict core",
        },
    ]

    checks: list[dict[str, Any]] = []
    mismatches: list[str] = []
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
            "step": "N6",
            "status": "N6_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FR_ROUTE_FRONTIER",
            "goal": "Check whether the current strict-core FR/topological route derives an internal selector source.",
            "scope": "current_strict_core_fr_route_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {
                "discharged": False,
                "reason": "the expected FR-route frontier has changed or required evidence objects are missing",
            },
            "required_next_step": "REVIEW_CHANGED_FR_ROUTE_FRONTIER_BEFORE_CLAIMING_FR_ROUTE_NONDERIVATION",
        }
    else:
        summary = {
            "step": "N6",
            "status": "N6_DISCHARGED_CURRENT_STRICT_CORE_FR_ROUTE_NONDERIVATION_NO_FALSE_PASS",
            "goal": "Discharge a route-specific theorem: the current strict-core FR/topological route does not yet derive an internal selector source.",
            "scope": "current_strict_core_fr_route_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "strict_core_fr_route_derives_residual_orientation_datum": False,
                "strict_core_fr_route_derives_actual_theta_source": False,
                "strict_core_fr_route_serves_as_internal_selector_source": False,
                "extra_strict_core_bridge_objects_required": True,
            },
            "missing_structure_classes": [
                "strict_derivation_of_sigma_int_candidate",
                "theorem_level_gauge_quotient_safety",
                "strict_core_equivalence_or_export_map_to_residual_orientation_datum",
                "strict_core_sigma_to_theta_map_or_internal_selector_family_derivation",
                "strict_core_actual_theta_source",
            ],
            "hard_limits": [
                "no global impossibility theorem",
                "no claim that future FR-based strict-core routes are impossible",
                "no claim that QW-2191 is discharged",
                "no claim that ToE is closed",
            ],
            "required_next_step": "ADD_ONE_MISSING_STRICT_CORE_FR_BRIDGE_OBJECT_AND_RERUN_P3_OR_FORMALIZE_A_STRONGER_NEGATIVE_THEOREM_IF_A_NEW_ARGUMENT_APPEARS",
        }

    out = ROOT / "generated" / "n6_current_strict_core_fr_route_nonderivation_theorem_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out)


if __name__ == "__main__":
    main()
