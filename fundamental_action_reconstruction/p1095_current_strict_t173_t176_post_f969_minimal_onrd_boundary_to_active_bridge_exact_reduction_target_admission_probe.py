#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-29"
ALLOWED_CONTRACT = "search_one_genuinely_new_non_same_lane_upgrade_route_within_exported_noncyclic_provider_split_family_or_one_genuinely_new_inversion_sensitive_source_side_provider_class"
EXPECTED_BOUNDARY = "MinimalOrientedNonreciprocalDephasingNewImportBoundary_v1"
EXPECTED_BOUNDARY_GRADE = "candidate_provider_class_seed_only"
EXPECTED_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1091 = GENERATED / "p1091_current_strict_t173_t176_post_f965_failure_map_to_exported_noncyclic_provider_split_non_same_lane_upgrade_route_decision_audit_probe_summary.json"
IN_F966 = GENERATED / "f966_current_strict_t173_t176_post_f965_failure_map_to_exported_noncyclic_provider_split_non_same_lane_upgrade_route_packet_summary.json"
IN_F969 = GENERATED / "f969_current_strict_t173_t176_minimal_oriented_nonreciprocal_dephasing_new_import_boundary_packet_summary.json"
IN_F960 = GENERATED / "f960_current_strict_t173_t176_existing_t183_residual_datum_pair12_orbit_direction_selection_bridge_target_packet_summary.json"

OUT_JSON = GENERATED / "p1095_current_strict_t173_t176_post_f969_minimal_onrd_boundary_to_active_bridge_exact_reduction_target_admission_probe.json"
OUT_SUMMARY = GENERATED / "p1095_current_strict_t173_t176_post_f969_minimal_onrd_boundary_to_active_bridge_exact_reduction_target_admission_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)
    prereq = [IN_P1091, IN_F966, IN_F969, IN_F960]
    missing = [rel(p) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P1095",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1091 = load_json(IN_P1091)
    f966 = load_json(IN_F966)
    f969 = load_json(IN_F969)
    f960 = load_json(IN_F960)

    checks: list[dict[str, Any]] = []
    blocking: list[str] = []

    def add_check(check_id: str, actual: Any, expected: Any, meaning: str) -> None:
        passed = actual == expected
        checks.append(
            {
                "id": check_id,
                "actual": actual,
                "expected": expected,
                "pass": passed,
                "meaning": meaning,
            }
        )
        if not passed:
            blocking.append(check_id)

    add_check(
        "new_inversion_sensitive_provider_class_route_is_allowed",
        (
            p1091.get("allowed_next_move_contract") == ALLOWED_CONTRACT
            and f966.get("allowed_next_move_contract") == ALLOWED_CONTRACT
        ),
        True,
        "F966/P1091 must still allow one genuinely new inversion-sensitive source-side provider-class route.",
    )
    add_check(
        "f969_boundary_seed_is_exported_only_as_seed",
        (
            f969.get("boundary_object_id") == EXPECTED_BOUNDARY
            and f969.get("boundary_grade") == EXPECTED_BOUNDARY_GRADE
            and bool(f969.get("admissible_as_candidate_provider_class_seed_only")) is True
        ),
        True,
        "F969 must export the minimal ONRD boundary only as a candidate provider-class seed.",
    )
    add_check(
        "boundary_and_bridge_share_same_frontier_context",
        (
            f969.get("frontier_context_bridge") == EXPECTED_BRIDGE
            and f960.get("target_object_id") == EXPECTED_BRIDGE
        ),
        True,
        "The new boundary and the active bridge must point to the same frontier context.",
    )
    add_check(
        "exact_reduction_is_still_missing",
        bool(f969.get("exact_reduction_to_frontier_context_bridge_exported")),
        False,
        "Exact reduction from the boundary to the active bridge must still be missing here.",
    )
    add_check(
        "same_lane_reentry_remains_disallowed",
        (
            bool(p1091.get("pair12_entry_same_lane_reentry_disallowed_as_primary_move"))
            and bool(p1091.get("pair_side_sharper_same_lane_reentry_disallowed_as_primary_move"))
            and bool(p1091.get("feeder_sharper_same_lane_reentry_disallowed_as_primary_move"))
        ),
        True,
        "Old same-lane reentry must remain blocked.",
    )
    add_check(
        "boundary_is_still_below_supplier_solution_and_orientation_export",
        (
            bool(f969.get("counts_as_lawful_supplier")) is False
            and bool(f969.get("counts_as_solution")) is False
            and bool(f969.get("counts_as_strict_physical_orientation_datum")) is False
        ),
        True,
        "The boundary must remain below supplier / solution / strict orientation export.",
    )

    discharged = len(blocking) == 0
    status = (
        "PASS_CURRENT_STRICT_T173_T176_POST_F969_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_EXACT_REDUCTION_TARGET_ADMITTED"
        if discharged
        else "P1095_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ROUTE_STATE"
    )

    artifact = {
        "stage": "P1095",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "p1091_route_decision_summary": rel(IN_P1091),
            "f966_route_packet_summary": rel(IN_F966),
            "f969_boundary_packet_summary": rel(IN_F969),
            "f960_active_bridge_target_summary": rel(IN_F960),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "admissible_target": {
            "object_id": "MinimalONRDBoundaryToActiveBridgeExactReductionTarget_v1",
            "source_boundary_object_id": EXPECTED_BOUNDARY,
            "source_boundary_grade": EXPECTED_BOUNDARY_GRADE,
            "target_bridge_object_id": EXPECTED_BRIDGE,
            "search_mode_selected": "genuinely_new_inversion_sensitive_source_side_provider_class_route",
            "within_exported_noncyclic_provider_split_family": False,
            "admissible_as_exact_reduction_target": discharged,
            "counts_as_lawful_supplier": False,
            "counts_as_solution": False,
            "counts_as_strict_physical_orientation_datum": False,
        },
        "classification": {
            "exact_reduction_already_exported": False,
            "next_honest_move": "audit_actual_realization_or_nonexport_status_of_the_new_exact_reduction_target_before_any_attempt",
        },
        "hard_limits": [
            "No exact reduction export.",
            "No lawful supplier export.",
            "No solution export.",
            "No strict physical orientation datum export.",
            "No T183 discharge.",
            "No T176 discharge.",
            "No QW-2191 discharge.",
            "No ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "target_object_id": artifact["admissible_target"]["object_id"],
        "source_boundary_object_id": artifact["admissible_target"]["source_boundary_object_id"],
        "target_bridge_object_id": artifact["admissible_target"]["target_bridge_object_id"],
        "search_mode_selected": artifact["admissible_target"]["search_mode_selected"],
        "within_exported_noncyclic_provider_split_family": artifact["admissible_target"]["within_exported_noncyclic_provider_split_family"],
        "admissible_as_exact_reduction_target": artifact["admissible_target"]["admissible_as_exact_reduction_target"],
        "counts_as_lawful_supplier": False,
        "counts_as_solution": False,
        "counts_as_strict_physical_orientation_datum": False,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
