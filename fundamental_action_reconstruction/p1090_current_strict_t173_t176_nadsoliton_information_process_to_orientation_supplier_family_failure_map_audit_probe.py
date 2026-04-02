#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-28"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F179 = ROOT / "F179_FIRST_CANONICAL_ONTOLOGY_SUPPORTED_FRACTAL_INFORMATION_NADSOLITON_PAIR_POPULATION_NONCYCLIC_ANCHOR_HYPOTHESIS_PACKET.md"
IN_T41 = ROOT / "T41_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_PREORIENTED_PRIMORDIAL_INFORMATION_BLOCKER_CUT_TARGET_SPEC.md"
IN_F196 = ROOT / "F196_FIRST_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_PRIMORDIAL_PREORIENTATION_LAW_TARGET_PACKET.md"
IN_P1011 = GENERATED / "p1011_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_provider_shift_candidate_reference_lane_admission_probe_summary.json"
IN_P1034 = GENERATED / "p1034_current_strict_qw2191_nadsoliton_neural_character_information_primary_selector_support_reference_admission_probe_summary.json"
IN_F953 = GENERATED / "f953_current_strict_qw2191_nadsoliton_neural_character_information_primary_selector_support_reference_packet_summary.json"
IN_P1047 = GENERATED / "p1047_current_strict_qw2191_post_stop_neural_bridge_lane_to_strict_int_selector_source_frontier_route_decision_audit_probe_summary.json"
IN_F956 = GENERATED / "f956_current_strict_qw2191_post_stop_neural_bridge_lane_to_strict_int_selector_source_frontier_route_decision_packet_summary.json"
IN_P1061 = GENERATED / "p1061_current_strict_t173_t176_existing_t183_residual_datum_pair12_orbit_direction_selection_bridge_existing_export_or_near_miss_candidate_audit_probe_summary.json"
IN_F960 = GENERATED / "f960_current_strict_t173_t176_existing_t183_residual_datum_pair12_orbit_direction_selection_bridge_target_packet_summary.json"
IN_P1062 = GENERATED / "p1062_current_strict_t173_t176_existing_t183_negative_3_cycle_sign_holonomy_obstruction_reference_admission_probe_summary.json"
IN_F961 = GENERATED / "f961_current_strict_t173_t176_existing_t183_negative_3_cycle_sign_holonomy_obstruction_reference_packet_summary.json"
IN_P1089 = GENERATED / "p1089_current_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_refinement_stagnation_and_stop_audit_probe_summary.json"
IN_F964 = GENERATED / "f964_current_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_refinement_stagnation_stop_packet_summary.json"
IN_P708 = GENERATED / "p708_current_strict_t173_frontier_dashboard_probe_summary.json"

OUT_JSON = GENERATED / "p1090_current_strict_t173_t176_nadsoliton_information_process_to_orientation_supplier_family_failure_map_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1090_current_strict_t173_t176_nadsoliton_information_process_to_orientation_supplier_family_failure_map_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [
        IN_F179,
        IN_T41,
        IN_F196,
        IN_P1011,
        IN_P1034,
        IN_F953,
        IN_P1047,
        IN_F956,
        IN_P1061,
        IN_F960,
        IN_P1062,
        IN_F961,
        IN_P1089,
        IN_F964,
        IN_P708,
    ]
    missing = [rel(p) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P1090",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    f179_text = IN_F179.read_text(encoding="utf-8")
    t41_text = IN_T41.read_text(encoding="utf-8")
    f196_text = IN_F196.read_text(encoding="utf-8")
    p1011 = load_json(IN_P1011)
    p1034 = load_json(IN_P1034)
    f953 = load_json(IN_F953)
    p1047 = load_json(IN_P1047)
    f956 = load_json(IN_F956)
    p1061 = load_json(IN_P1061)
    f960 = load_json(IN_F960)
    p1062 = load_json(IN_P1062)
    f961 = load_json(IN_F961)
    p1089 = load_json(IN_P1089)
    f964 = load_json(IN_F964)
    p708 = load_json(IN_P708)

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
        "f179_explicitly_treats_nadsoliton_as_primordial_information",
        "the nadsoliton is itself the primordial information of the universe" in f179_text,
        True,
        "F179 already treats nadsoliton as an informational carrier, not merely a geometric object.",
    )
    add_check(
        "t41_keeps_preorientation_on_same_informational_carrier",
        "the nadsoliton itself is the primordial information of the universe" in t41_text,
        True,
        "T41 keeps the preorientation proposal on the same informational carrier rather than a lower substrate.",
    )
    add_check(
        "f196_names_future_preorientation_law_only",
        "Pi_primordial_preorientation_law_target_v1" in f196_text,
        True,
        "F196 names one missing preorientation law but keeps it future-only.",
    )
    add_check(
        "p1011_info_primary_lane_is_reference_only",
        (
            bool(p1011.get("info_primary_scpc_like_selector_provider_shift_candidate_reference_lane_admitted"))
            and p1011.get("candidate_reference_lane_grade") == "reference_context_candidate_only"
            and bool(p1011.get("provider_class_shift_realized")) is False
            and bool(p1011.get("strict_selector_interface_exported")) is False
        ),
        True,
        "The info-primary SCPC-like lane is admitted only as reference context and does not realize provider shift or interface export.",
    )
    add_check(
        "p1034_f953_neural_character_lane_stays_support_reference_only",
        (
            bool(p1034.get("nadsoliton_neural_character_support_reference_admitted"))
            and p1034.get("support_reference_grade") == "cross_repo_support_reference_only"
            and bool(p1034.get("strict_selector_interface_exported")) is False
            and bool(p1034.get("strict_selector_source_exported")) is False
            and f953.get("support_reference_grade") == "cross_repo_support_reference_only"
            and f953.get("strict_selector_interface_status") == "blocked_nonexport"
            and f953.get("strict_selector_source_status") == "blocked_nonexport"
        ),
        True,
        "The nadsoliton neural/information lane exists, but only as support-reference with blocked strict interface/source.",
    )
    add_check(
        "p1047_f956_positive_bridge_not_justified",
        (
            bool(p1047.get("positive_bridge_branch_selection_not_justified"))
            and bool(p1047.get("strict_int_selector_source_frontier_open"))
            and p1047.get("primary_continuation_route") == "explicit_strict_internal_selector_source_derivation_frontier"
            and f956.get("primary_continuation_route") == "explicit_strict_internal_selector_source_derivation_frontier"
        ),
        True,
        "After the stopped neural bridge lane, the repo rejoined the explicit strict internal selector-source frontier instead of promoting the bridge.",
    )
    add_check(
        "p1061_f960_no_actual_t183_supplier_found",
        (
            bool(p1061.get("current_repo_already_exports_actual_t183_supplier")) is False
            and bool(p1061.get("direction_free_shannon_candidate_family_negative"))
            and f960.get("status")
            == "F960_EXECUTED_CURRENT_STRICT_T173_T176_EXISTING_T183_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_SELECTION_BRIDGE_TARGET_PACKET_NO_FALSE_PASS"
        ),
        True,
        "The active T183 bridge target still has no lawful supplier on current repo state.",
    )
    add_check(
        "p1062_f961_negative_3_cycle_is_reference_not_supplier",
        (
            bool(p1062.get("admissible_as_reference_only"))
            and bool(p1062.get("admissible_as_lawful_supplier")) is False
            and bool(p1062.get("counts_as_strict_physical_orientation_datum")) is False
            and bool(f961.get("counts_as_reference_only"))
            and bool(f961.get("counts_as_lawful_supplier")) is False
        ),
        True,
        "The negative 3-cycle/sign-holonomy object is frontier-relevant but remains obstruction-reference only.",
    )
    add_check(
        "p1089_f964_downstream_feeder_same_lane_is_stopped",
        (
            bool(p1089.get("stop_condition_triggered"))
            and bool(p1089.get("same_lane_stagnation_boundary_reached"))
            and bool(p1089.get("restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route"))
            and bool(f964.get("same_lane_stagnation_boundary_reached"))
        ),
        True,
        "The downstream feeder same-lane witness family is formally stopped and cannot be continued honestly by deeper same-lane descent.",
    )
    add_check(
        "p708_current_strict_frontier_still_requires_new_provider_class",
        (
            bool(p708.get("QW2191_kernel_alone_discharge")) is False
            and bool(p708.get("t176_global_provider_exported")) is False
            and bool(p708.get("surviving_pair12_residual_datum_carrier_remains_selector_neutral"))
            and bool(p708.get("next_honest_primary_t173_move_requires_genuinely_new_provider_class_beyond_current_exported_continuation_family"))
        ),
        True,
        "The strict frontier still lacks a lawful provider that would turn the informational process reading into strict orientation selection.",
    )

    process_reading_is_repo_real = len(blocking) == 0

    families = [
        {
            "family_id": "canonical_preorientation_and_fractal_information_roots",
            "repo_objects": [rel(IN_F179), rel(IN_T41), rel(IN_F196)],
            "role": "upstream_ontology_support_for_nadsoliton_as_information_process",
            "strongest_export_grade": "future_only_provider_hypothesis_or_law_target_only",
            "counts_as_lawful_supplier": False,
            "exact_failure_localization": (
                "no actual preorientation law, no actual pair-indexed support, no actual E_orient, "
                "no actual selector-source support"
            ),
        },
        {
            "family_id": "info_primary_scpc_like_selector_provider_shift_candidate_lane",
            "repo_objects": [rel(IN_P1011)],
            "role": "reference_context_process_to_selector_candidate_lane",
            "strongest_export_grade": p1011.get("candidate_reference_lane_grade"),
            "counts_as_lawful_supplier": False,
            "exact_failure_localization": (
                "provider_class_shift_realized = false; strict_selector_interface_exported = false"
            ),
        },
        {
            "family_id": "nadsoliton_neural_character_information_primary_selector_support_reference",
            "repo_objects": [rel(IN_P1034), rel(IN_F953)],
            "role": "cross_repo_support_reference_for_information_primary_selector_hypothesis",
            "strongest_export_grade": p1034.get("support_reference_grade"),
            "counts_as_lawful_supplier": False,
            "exact_failure_localization": (
                "strict_selector_interface_exported = false; strict_selector_source_exported = false"
            ),
        },
        {
            "family_id": "post_stop_neural_bridge_route_decision",
            "repo_objects": [rel(IN_P1047), rel(IN_F956)],
            "role": "frontier_route_decision_after_processive_bridge_attempt",
            "strongest_export_grade": "route_decision_only",
            "counts_as_lawful_supplier": False,
            "exact_failure_localization": (
                "positive_bridge_branch_selection_not_justified = true; continuation returns to explicit_strict_internal_selector_source_derivation_frontier"
            ),
        },
        {
            "family_id": "t183_residual_datum_pair12_orbit_direction_selection_bridge_search",
            "repo_objects": [rel(IN_P1061), rel(IN_F960)],
            "role": "strict_process_to_orientation_bridge_search",
            "strongest_export_grade": "target_packet_only",
            "counts_as_lawful_supplier": False,
            "exact_failure_localization": (
                "current_repo_already_exports_actual_t183_supplier = false; direction_free_shannon_candidate_family_negative = true"
            ),
        },
        {
            "family_id": "negative_3_cycle_sign_holonomy_obstruction_reference",
            "repo_objects": [rel(IN_P1062), rel(IN_F961)],
            "role": "orientation_obstruction_reference_for_supplier_search",
            "strongest_export_grade": p1062.get("reference_grade"),
            "counts_as_lawful_supplier": False,
            "exact_failure_localization": (
                "admissible_as_lawful_supplier = false; counts_as_strict_physical_orientation_datum = false"
            ),
        },
        {
            "family_id": "downstream_feeder_support_side_provider_support_same_lane_family",
            "repo_objects": [rel(IN_P1089), rel(IN_F964)],
            "role": "downstream_source_side_support_refinement_family_below_process_reading",
            "strongest_export_grade": "same_lane_family_formally_stopped",
            "counts_as_lawful_supplier": False,
            "exact_failure_localization": (
                "same_lane_stagnation_boundary_reached = true; restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route = true"
            ),
        },
    ]

    active_missing_bridge = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"
    current_supplier_status = False
    status = (
        "PASS_CURRENT_STRICT_T173_T176_NADSOLITON_INFORMATION_PROCESS_TO_ORIENTATION_SUPPLIER_FAMILY_FAILURE_MAP_AUDITED"
        if process_reading_is_repo_real
        else "P1090_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FRONTIER_STATE"
    )

    artifact = {
        "stage": "P1090",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "F179": rel(IN_F179),
            "T41": rel(IN_T41),
            "F196": rel(IN_F196),
            "P1011": rel(IN_P1011),
            "P1034": rel(IN_P1034),
            "F953": rel(IN_F953),
            "P1047": rel(IN_P1047),
            "F956": rel(IN_F956),
            "P1061": rel(IN_P1061),
            "F960": rel(IN_F960),
            "P1062": rel(IN_P1062),
            "F961": rel(IN_F961),
            "P1089": rel(IN_P1089),
            "F964": rel(IN_F964),
            "P708": rel(IN_P708),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "failure_map_object_id": "NadsolitonInformationProcessToOrientationSupplierFailureMap_global_C_v1_strict_v1",
        "process_reading_of_nadsoliton_is_repo_real": process_reading_is_repo_real,
        "current_repo_exports_lawful_process_to_orientation_supplier": current_supplier_status,
        "active_missing_bridge": active_missing_bridge,
        "tested_families": families,
        "current_strict_frontier_snapshot": {
            "QW2191_kernel_alone_discharge": bool(p708.get("QW2191_kernel_alone_discharge")),
            "t176_global_provider_exported": bool(p708.get("t176_global_provider_exported")),
            "surviving_pair12_residual_datum_carrier_remains_selector_neutral": bool(
                p708.get("surviving_pair12_residual_datum_carrier_remains_selector_neutral")
            ),
            "next_honest_primary_t173_move_requires_genuinely_new_provider_class_beyond_current_exported_continuation_family": bool(
                p708.get("next_honest_primary_t173_move_requires_genuinely_new_provider_class_beyond_current_exported_continuation_family")
            ),
        },
        "audit_conclusion": {
            "process_reading_is_present_in_repo": process_reading_is_repo_real,
            "informational_process_has_not_yet_become_lawful_orientation_supplier": True,
            "current_exact_missing_bridge": active_missing_bridge,
            "next_honest_move": (
                "search_one_genuinely_new_inversion_sensitive_source_side_provider_class_or_non_same_lane_upgrade_route_"
                "that_can_lawfully_bridge_process_reading_to_strict_orientation_selection"
            ),
        },
        "hard_limits": [
            "No lawful process-to-orientation supplier exported on current repo state.",
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
        "failure_map_object_id": artifact["failure_map_object_id"],
        "process_reading_of_nadsoliton_is_repo_real": artifact["process_reading_of_nadsoliton_is_repo_real"],
        "current_repo_exports_lawful_process_to_orientation_supplier": False,
        "active_missing_bridge": active_missing_bridge,
        "tested_family_count": len(families),
        "next_honest_move": artifact["audit_conclusion"]["next_honest_move"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
