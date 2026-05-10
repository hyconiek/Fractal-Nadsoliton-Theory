#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1061 = GENERATED / "p1061_current_strict_t173_t176_existing_t183_residual_datum_pair12_orbit_direction_selection_bridge_existing_export_or_near_miss_candidate_audit_probe_summary.json"
IN_F959 = GENERATED / "f959_current_strict_t173_t176_post_stop_source_side_input_leg_lane_to_existing_t183_residual_datum_pair12_orbit_direction_frontier_route_decision_packet_summary.json"

OUT_JSON = GENERATED / "f960_current_strict_t173_t176_existing_t183_residual_datum_pair12_orbit_direction_selection_bridge_target_packet.json"
OUT_SUMMARY = GENERATED / "f960_current_strict_t173_t176_existing_t183_residual_datum_pair12_orbit_direction_selection_bridge_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1061, IN_F959]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "F960",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1061 = load_json(IN_P1061)
    f959 = load_json(IN_F959)

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

    p1061_audit_passed = (
        p1061.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_EXISTING_T183_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_SELECTION_BRIDGE_EXISTING_EXPORT_OR_NEAR_MISS_CANDIDATES_AUDITED_NO_LAWFUL_SUPPLIER_FOUND"
        and p1061.get("t183_target_name")
        == "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"
        and p1061.get("current_repo_already_exports_actual_t183_supplier") is False
        and p1061.get("next_honest_move")
        == "freeze_the_exact_missing_t183_target_as_the_active_post_stop_frontier_packet_without_promoting_t184_or_t185_by_fiat"
    )

    f959_route_exported = (
        f959.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_POST_STOP_SOURCE_SIDE_INPUT_LEG_LANE_TO_EXISTING_T183_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_FRONTIER_ROUTE_DECISION_PACKET_EXPORTED"
        and f959.get("packet_exported_on_current_repo_state") is True
        and f959.get("primary_continuation_route")
        == "existing_t183_residual_datum_pair12_orbit_direction_selection_frontier"
    )

    target_packet_exported_on_current_repo_state = p1061_audit_passed and f959_route_exported

    add_check(
        "p1061_audit_passed",
        p1061_audit_passed,
        True,
        "P1061 already proves that no current export lawfully supplies T183 and that the exact target should be frozen.",
    )
    add_check(
        "f959_route_exported",
        f959_route_exported,
        True,
        "F959 already exports the post-stop route decision to the existing T183 frontier.",
    )
    add_check(
        "target_packet_exported_on_current_repo_state",
        target_packet_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one honest exact target packet for the active post-stop T183 frontier.",
    )

    status = (
        "F960_EXECUTED_CURRENT_STRICT_T173_T176_EXISTING_T183_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_SELECTION_BRIDGE_TARGET_PACKET_NO_FALSE_PASS"
        if not blocking and target_packet_exported_on_current_repo_state
        else "F960_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "F960",
        "packet_name": "CurrentStrictT173T176ExistingT183ResidualDatumPair12OrbitDirectionSelectionBridgeTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p1061_existing_export_or_near_miss_candidate_audit": rel(IN_P1061),
            "f959_post_stop_route_decision_packet": rel(IN_F959),
        },
        "checks": checks,
        "blocking_checks": blocking,
        "target_object": {
            "object_id": "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1",
            "goal": "Freeze the exact missing bridge selecting between the surviving residual-datum pair1/pair2 orbit-direction branches delta_k and delta_-k on full C_v1 after the source-side input-leg same lane has already been stopped.",
            "required_fields": [
                {
                    "name": "post_stop_route_decision_ref",
                    "required": True,
                    "hard_limit": "Must keep explicit that this target is active only as the lawful post-stop continuation route.",
                },
                {
                    "name": "exact_pair12_orbit_direction_split_ref",
                    "required": True,
                    "hard_limit": "Must preserve the exact surviving split delta_k versus delta_-k on the residual-datum carrier.",
                },
                {
                    "name": "direction_free_shannon_negative_boundary_ref",
                    "required": True,
                    "hard_limit": "Must keep explicit that the current direction-free Shannon lane is insufficient.",
                },
                {
                    "name": "w_break_near_miss_nonpromotion_ref",
                    "required": True,
                    "hard_limit": "Must keep explicit that the current w_break family remains below typed promotion and below the admitted active-primary route.",
                },
                {
                    "name": "exact_t183_target_statement",
                    "required": True,
                    "hard_limit": "Must state the exact missing T183 object without claiming that it already exists.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny T183 discharge, T176 discharge, strict physical orientation export, QW-2191 discharge, and ToE closure.",
                },
            ],
        },
        "target_refs": {
            "post_stop_route_decision_ref": "existing_t183_residual_datum_pair12_orbit_direction_selection_frontier",
            "exact_pair12_orbit_direction_split_ref": {
                "pair1_branch": "delta_k_positive_index_branch",
                "pair2_branch": "delta_minus_k_negative_index_branch",
            },
            "direction_free_shannon_negative_boundary_ref": "DirectionFreeShannonResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1",
            "w_break_near_miss_nonpromotion_ref": "WBreakWitnessPayloadResidualDatumPair12OrbitDirectionPromotionBridge_global_C_v1_strict_v1",
            "exact_t183_target_statement": {
                "object_id": "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1",
                "typed_meaning": "one strict source-side bridge on full C_v1 selecting between the surviving residual-datum pair1/pair2 orbit-direction branches delta_k and delta_-k without smuggling nonadmissible witness-side or rooted-convention data into strict-core physics",
            },
        },
        "current_honest_reading": [
            "The exact surviving source-side ambiguity is already localized as delta_k versus delta_-k on the residual-datum carrier.",
            "The current direction-free Shannon lane is already proven insufficient for selecting one of those two branches.",
            "The current w_break witness payload already separates the two branches, but remains below typed promotion and below the admitted active-primary route.",
            "So the exact missing active post-stop object is now the T183 bridge itself.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any genuinely new inversion-sensitive source-side provider class beyond the current direction-free Shannon lane and beyond nonadmissible witness-side promotion can lawfully realize ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1.",
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "target_object_id": artifact["target_object"]["object_id"],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
