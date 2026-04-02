#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1060 = GENERATED / "p1060_current_strict_t173_t176_post_stop_source_side_input_leg_lane_to_existing_t183_residual_datum_pair12_orbit_direction_frontier_route_decision_audit_probe_summary.json"
IN_F959 = GENERATED / "f959_current_strict_t173_t176_post_stop_source_side_input_leg_lane_to_existing_t183_residual_datum_pair12_orbit_direction_frontier_route_decision_packet_summary.json"
IN_P708 = GENERATED / "p708_current_strict_t173_frontier_dashboard_probe_summary.json"
IN_P729 = GENERATED / "p729_current_strict_t183_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"
IN_P730 = GENERATED / "p730_current_strict_t184_direction_free_shannon_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"
IN_P731 = GENERATED / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe_summary.json"

OUT_JSON = GENERATED / "p1061_current_strict_t173_t176_existing_t183_residual_datum_pair12_orbit_direction_selection_bridge_existing_export_or_near_miss_candidate_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1061_current_strict_t173_t176_existing_t183_residual_datum_pair12_orbit_direction_selection_bridge_existing_export_or_near_miss_candidate_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1060, IN_F959, IN_P708, IN_P729, IN_P730, IN_P731]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1061",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1060 = load_json(IN_P1060)
    f959 = load_json(IN_F959)
    p708 = load_json(IN_P708)
    p729 = load_json(IN_P729)
    p730 = load_json(IN_P730)
    p731 = load_json(IN_P731)

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

    post_stop_route_to_existing_t183_exported = (
        p1060.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_POST_STOP_SOURCE_SIDE_INPUT_LEG_LANE_TO_EXISTING_T183_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_FRONTIER_ROUTE_DECISION_AUDITED"
        and p1060.get("primary_continuation_route")
        == "existing_t183_residual_datum_pair12_orbit_direction_selection_frontier"
        and f959.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_POST_STOP_SOURCE_SIDE_INPUT_LEG_LANE_TO_EXISTING_T183_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_FRONTIER_ROUTE_DECISION_PACKET_EXPORTED"
        and f959.get("packet_exported_on_current_repo_state") is True
        and f959.get("primary_continuation_route")
        == "existing_t183_residual_datum_pair12_orbit_direction_selection_frontier"
    )

    t183_bridge_still_unexported_but_exactly_localized = (
        p729.get("status")
        == "PASS_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_SELECTION_BRIDGE_NONEXPORT_AUDITED"
        and p729.get("t183_target_name")
        == "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"
        and p729.get("t183_target_exported_on_current_repo_state") is False
        and p729.get("remaining_pair12_split_localized_as_opposite_orbit_directions") is True
        and p729.get("pair1_orbit_branch_kind") == "delta_k_positive_index_branch"
        and p729.get("pair2_orbit_branch_kind") == "delta_minus_k_negative_index_branch"
    )

    direction_free_shannon_candidate_family_negative = (
        p730.get("status")
        == "PASS_DIRECTION_FREE_SHANNON_PAIR12_ORBIT_DIRECTION_SELECTION_BRIDGE_NONEXPORT_AUDITED"
        and p730.get("t184_target_exported_on_current_repo_state") is False
        and p730.get("current_direction_free_shannon_lane_already_exports_pair1_pair2_o2_to_z2_cuts")
        is True
        and p730.get("current_direction_free_shannon_lane_selects_pair12_orbit_direction_branch") is False
    )

    w_break_near_miss_family_real_but_still_unlawful = (
        p731.get("status")
        == "PASS_W_BREAK_WITNESS_PAYLOAD_PAIR12_ORBIT_DIRECTION_PROMOTION_BRIDGE_NONEXPORT_AUDITED"
        and p731.get("t185_target_exported_on_current_repo_state") is False
        and p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches")
        is True
        and p708.get("current_pair12_witness_split_current_exported_continuation_family_named_members_all_real")
        is True
        and p708.get("current_pair12_witness_split_current_exported_continuation_family_named_members_all_negative")
        is True
        and p708.get("same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move")
        is True
    )

    no_current_lawful_t183_supplier_found = (
        post_stop_route_to_existing_t183_exported
        and t183_bridge_still_unexported_but_exactly_localized
        and direction_free_shannon_candidate_family_negative
        and w_break_near_miss_family_real_but_still_unlawful
        and p708.get("t183_residual_datum_pair12_orbit_direction_selection_bridge_exported") is False
    )

    add_check(
        "post_stop_route_to_existing_t183_exported",
        post_stop_route_to_existing_t183_exported,
        True,
        "The post-stop route decision to the existing T183 frontier is already exported.",
    )
    add_check(
        "t183_bridge_still_unexported_but_exactly_localized",
        t183_bridge_still_unexported_but_exactly_localized,
        True,
        "P729 already localizes the exact delta_k vs delta_-k split while keeping T183 unexported.",
    )
    add_check(
        "direction_free_shannon_candidate_family_negative",
        direction_free_shannon_candidate_family_negative,
        True,
        "The current direction-free Shannon candidate family is real but explicitly insufficient for T183.",
    )
    add_check(
        "w_break_near_miss_family_real_but_still_unlawful",
        w_break_near_miss_family_real_but_still_unlawful,
        True,
        "The current w_break family is a real inversion-sensitive near miss, but still remains below typed promotion and outside the admitted active-primary route.",
    )
    add_check(
        "no_current_lawful_t183_supplier_found",
        no_current_lawful_t183_supplier_found,
        True,
        "Therefore no current export yet lawfully supplies the exact T183 bridge.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_T183_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_SELECTION_BRIDGE_EXISTING_EXPORT_OR_NEAR_MISS_CANDIDATES_AUDITED_NO_LAWFUL_SUPPLIER_FOUND"
        if not blocking and no_current_lawful_t183_supplier_found
        else "FAIL_CURRENT_STRICT_T173_T176_EXISTING_T183_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_SELECTION_BRIDGE_EXISTING_EXPORT_OR_NEAR_MISS_CANDIDATE_AUDIT"
    )

    artifact = {
        "stage": "P1061",
        "status": status,
        "as_of": AS_OF,
        "post_stop_route_to_existing_t183_exported": post_stop_route_to_existing_t183_exported,
        "t183_target_name": "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1",
        "current_repo_already_exports_actual_t183_supplier": False,
        "t183_bridge_still_unexported_but_exactly_localized": t183_bridge_still_unexported_but_exactly_localized,
        "direction_free_shannon_candidate_family_negative": direction_free_shannon_candidate_family_negative,
        "w_break_near_miss_family_real_but_still_unlawful": w_break_near_miss_family_real_but_still_unlawful,
        "current_best_positive_near_miss_family": (
            "w_break_witness_payload_residual_datum_pair12_branch_separation_family_below_typed_promotion_and_below_admitted_active_primary_route"
        ),
        "current_negative_candidate_family": "direction_free_shannon_ord_reference_pair12_lane",
        "no_current_lawful_t183_supplier_found": no_current_lawful_t183_supplier_found,
        "checks": checks,
        "blocking_checks": blocking,
        "next_honest_move": "freeze_the_exact_missing_t183_target_as_the_active_post_stop_frontier_packet_without_promoting_t184_or_t185_by_fiat",
        "hard_limits": [
            "No T183 discharge claim.",
            "No actual source-side branch selection claim.",
            "No T176 discharge claim.",
            "No strict physical orientation datum claim.",
            "No QW-2191 discharge claim.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "t183_target_name": artifact["t183_target_name"],
        "current_repo_already_exports_actual_t183_supplier": artifact["current_repo_already_exports_actual_t183_supplier"],
        "direction_free_shannon_candidate_family_negative": artifact[
            "direction_free_shannon_candidate_family_negative"
        ],
        "w_break_near_miss_family_real_but_still_unlawful": artifact[
            "w_break_near_miss_family_real_but_still_unlawful"
        ],
        "next_honest_move": artifact["next_honest_move"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
