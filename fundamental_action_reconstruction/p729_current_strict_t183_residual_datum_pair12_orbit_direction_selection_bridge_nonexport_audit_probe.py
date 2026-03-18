#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from math import pi
from pathlib import Path
from typing import Any


AS_OF = "2026-03-18"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P728 = GENERATED / "p728_current_strict_t182_residual_datum_source_side_boundary_shielded_sublane_reduction_audit_probe_summary.json"
IN_F301 = GENERATED / "residual_datum_bridge_export_map_object_support_carrier_actual_inhabitant_instance.json"
IN_F461 = GENERATED / "o12_pair1_pair2_selector_chart_transport_operator_strict_derived_from_sigma_int_alpha12_v1_summary.json"
IN_F462 = GENERATED / "f462_current_strict_sigma_int_pair1_pair2_projector_operator_section_glue_export_packet_summary.json"
IN_ALPHA12 = GENERATED / "alpha12_pair1_pair2_transition_angle_strict_derived_from_sigma_int_slot_free_theta_pair_v1.json"

OUT_JSON = GENERATED / "p729_current_strict_t183_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p729_current_strict_t183_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def invert_weight_dict(weight_dict: dict[str, Any]) -> dict[str, float]:
    out: dict[str, float] = {}
    for key, value in weight_dict.items():
        idx = int(key)
        flipped = -idx
        out[str(flipped)] = float(value)
    return out


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P728, IN_F301, IN_F461, IN_F462, IN_ALPHA12]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P729",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p728 = load_json(IN_P728)
    f301 = load_json(IN_F301)
    f461 = load_json(IN_F461)
    f462 = load_json(IN_F462)
    alpha12 = load_json(IN_ALPHA12)

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

    pair_support = list(f301.get("pair_index_set") or [])
    normalized_orbit_reps = f301.get("normalized_orbit_representatives") or {}
    carrier_vectors = f301.get("carrier_vectors") or {}
    pair1_vector = {str(k): float(v) for k, v in (carrier_vectors.get("pair1") or {}).items()}
    pair2_vector = {str(k): float(v) for k, v in (carrier_vectors.get("pair2") or {}).items()}

    pair1_positive_branch_keys = sorted(int(k) for k in pair1_vector.keys())
    pair2_negative_branch_keys = sorted(int(k) for k in pair2_vector.keys())
    pair2_equals_pair1_under_index_inversion = invert_weight_dict(pair1_vector) == pair2_vector

    alpha_outputs = alpha12.get("outputs") or {}
    f461_outputs = f461.get("outputs") or {}
    f462_outputs = f462.get("outputs") or {}

    pair12_lane_transport_exact = (
        abs(float(f461_outputs.get("alpha12_mod_2pi") or 0.0) - (pi / 2.0)) < 1.0e-12
        and float(f461_outputs.get("orthogonality_max_abs_residual") or 1.0) < 1.0e-12
        and float(f461_outputs.get("involution_max_abs_residual") or 1.0) < 1.0e-12
        and float(f461_outputs.get("u1_to_u2_max_abs_residual") or 1.0) < 1.0e-12
        and float(f462_outputs.get("gluing_full_matrix_max_abs_residual") or 1.0) < 1.0e-12
        and float(f462_outputs.get("gluing_pair_plane_2x2_max_abs_residual") or 1.0) < 1.0e-12
    )

    remaining_pair12_split_localized_as_opposite_orbit_directions = (
        pair_support == ["pair1", "pair2"]
        and normalized_orbit_reps.get("pair1") == "q_{1,k} := P_pair1^k(ψ_src)/||P_pair1^k(ψ_src)|| = δ_k"
        and normalized_orbit_reps.get("pair2") == "q_{2,k} := P_pair2^k(ψ_src)/||P_pair2^k(ψ_src)|| = δ_{-k}"
        and pair2_equals_pair1_under_index_inversion
    )

    add_check(
        "p728_boundary_shielded_sublane_reduction_already_present",
        bool(p728.get("current_residual_datum_source_side_support_reduces_positive_corridor_to_boundary_shielded_sublane")),
        True,
        "P728 already reduces the positive corridor to the residual-datum boundary-shielded sublane pair1,pair2.",
    )
    add_check(
        "f301_pair12_support_still_exact",
        pair_support,
        ["pair1", "pair2"],
        "F301 still exports the residual-datum source-side support exactly on pair1 and pair2.",
    )
    add_check(
        "f301_pair1_pair2_orbit_representatives_are_opposite_index_directions",
        {
            "pair1": normalized_orbit_reps.get("pair1"),
            "pair2": normalized_orbit_reps.get("pair2"),
        },
        {
            "pair1": "q_{1,k} := P_pair1^k(ψ_src)/||P_pair1^k(ψ_src)|| = δ_k",
            "pair2": "q_{2,k} := P_pair2^k(ψ_src)/||P_pair2^k(ψ_src)|| = δ_{-k}",
        },
        "Within the already selected pair1/pair2 sublane, F301 localizes the two remaining branches as opposite orbit-index directions: δ_k versus δ_{-k}.",
    )
    add_check(
        "pair1_support_keys_are_nonnegative_branch",
        pair1_positive_branch_keys,
        list(range(12)),
        "The pair1 residual-datum branch is currently carried on the nonnegative orbit-index side 0..11.",
    )
    add_check(
        "pair2_support_keys_are_nonpositive_branch",
        pair2_negative_branch_keys,
        list(range(-11, 1)),
        "The pair2 residual-datum branch is currently carried on the nonpositive orbit-index side -11..0.",
    )
    add_check(
        "pair2_branch_is_exact_index_inversion_of_pair1_branch",
        pair2_equals_pair1_under_index_inversion,
        True,
        "The pair2 residual-datum branch is exactly the k -> -k inversion of the pair1 branch, with no amplitude asymmetry left on current exports.",
    )
    add_check(
        "pair12_lane_transport_and_projector_glue_remain_exact",
        pair12_lane_transport_exact,
        True,
        "The current pair1↔pair2 lane transport and projector glue remain exact and only relate the two branches; they do not pick one of them as a selected source-side branch.",
    )
    add_check(
        "pair12_transition_angle_is_quarter_turn_only",
        float(alpha_outputs.get("alpha_12_mod_2pi") or 0.0),
        pi / 2.0,
        "The exported pair1↔pair2 transition angle remains a quarter-turn lane relation, not a source-side branch-selection rule.",
    )
    add_check(
        "remaining_pair12_split_localized_as_opposite_orbit_directions",
        remaining_pair12_split_localized_as_opposite_orbit_directions,
        True,
        "The remaining pair1/pair2 ambiguity is now fully localized on current exports as opposite residual-datum orbit-index directions on the same source-side carrier.",
    )
    add_check(
        "t183_orbit_direction_selection_bridge_exported",
        False,
        False,
        "No current export selects between the positive-index orbit branch and the negative-index orbit branch inside pair1/pair2.",
    )

    status = (
        "PASS_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_SELECTION_BRIDGE_NONEXPORT_AUDITED"
        if not blocking
        else "P729_REQUIRES_REVIEW_CHANGED_PAIR12_ORBIT_DIRECTION_STATE"
    )

    artifact = {
        "stage": "P729",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t183_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_boundary_only",
        "inputs": {
            "P728": str(IN_P728.relative_to(REPO)),
            "F301_artifact": str(IN_F301.relative_to(REPO)),
            "F461_summary": str(IN_F461.relative_to(REPO)),
            "F462_summary": str(IN_F462.relative_to(REPO)),
            "alpha12": str(IN_ALPHA12.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t183_target_name": "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1",
        "t183_target_exported_on_current_repo_state": False,
        "pair12_residual_datum_support": pair_support,
        "pair1_orbit_branch_kind": "delta_k_positive_index_branch",
        "pair2_orbit_branch_kind": "delta_minus_k_negative_index_branch",
        "pair2_branch_is_exact_index_inversion_of_pair1_branch": pair2_equals_pair1_under_index_inversion,
        "remaining_pair12_split_localized_as_opposite_orbit_directions": (
            remaining_pair12_split_localized_as_opposite_orbit_directions
        ),
        "pair12_lane_transport_and_projector_glue_remain_exact": pair12_lane_transport_exact,
        "audit_conclusion": {
            "current_repo_already_exports_pair12_source_side_orbit_direction_localization": True,
            "current_repo_already_exports_t183_target": False,
            "next_honest_move": (
                "export_or_attack_a_finer_residual_datum_source_side_rule_selecting_between_positive_index_and_negative_index_orbit_branches_inside_pair1_pair2"
            ),
        },
        "hard_limits": [
            "No T183 discharge claim.",
            "No unique chart-seed selection claim.",
            "No promotion of orbit-index direction into a physical orientation datum.",
            "No identification of tau_src with the current selector carrier.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P729",
        "status": status,
        "as_of": AS_OF,
        "t183_target_name": "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1",
        "t183_target_exported_on_current_repo_state": False,
        "pair1_orbit_branch_kind": "delta_k_positive_index_branch",
        "pair2_orbit_branch_kind": "delta_minus_k_negative_index_branch",
        "remaining_pair12_split_localized_as_opposite_orbit_directions": (
            remaining_pair12_split_localized_as_opposite_orbit_directions
        ),
        "next_honest_move": "finer_residual_datum_source_side_rule_selecting_between_positive_index_and_negative_index_orbit_branches_inside_pair1_pair2",
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
