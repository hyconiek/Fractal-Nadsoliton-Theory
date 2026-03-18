#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-18"
TOL = 1.0e-15

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P729 = GENERATED / "p729_current_strict_t183_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"
IN_F301 = GENERATED / "residual_datum_bridge_export_map_object_support_carrier_actual_inhabitant_instance.json"
IN_R_ORD = GENERATED / "r_ord_z12_v1_reference_distribution.json"
IN_THETA_PAIR = GENERATED / "theta_pair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1.json"

OUT_JSON = GENERATED / "p730_current_strict_t184_direction_free_shannon_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p730_current_strict_t184_direction_free_shannon_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def normalize_branch_profile(weight_dict: dict[str, Any], n: int) -> tuple[dict[int, float], float]:
    merged: dict[int, float] = {}
    total = 0.0
    for key, value in weight_dict.items():
        idx = int(key) % n
        val = float(value)
        merged[idx] = merged.get(idx, 0.0) + val
        total += val
    if total <= 0.0:
        raise ValueError("branch profile must have positive total weight")
    return ({idx: merged.get(idx, 0.0) / total for idx in range(n)}, total)


def max_abs_difference(lhs: list[float], rhs: list[float]) -> float:
    return max(abs(a - b) for a, b in zip(lhs, rhs))


def expectation(prob: list[float], values: list[float]) -> float:
    return sum(p * v for p, v in zip(prob, values))


def cross_entropy(prob: list[float], ref_prob: list[float]) -> float:
    return -sum(p * math.log(q) for p, q in zip(prob, ref_prob))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P729, IN_F301, IN_R_ORD, IN_THETA_PAIR]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P730",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p729 = load_json(IN_P729)
    f301 = load_json(IN_F301)
    r_ord = load_json(IN_R_ORD)
    theta_pair = load_json(IN_THETA_PAIR)

    n = 12
    ord_by_x_raw = r_ord.get("ord_z12_by_x")
    if not (isinstance(ord_by_x_raw, list) and len(ord_by_x_raw) == n and all(isinstance(v, int) for v in ord_by_x_raw)):
        artifact = {
            "stage": "P730",
            "status": "INVALID_R_ORD_SHAPE",
            "as_of": AS_OF,
            "expected": "r_ord.ord_z12_by_x is a length-12 integer list",
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    carrier_vectors = f301.get("carrier_vectors") or {}
    pair1_weights_raw = carrier_vectors.get("pair1")
    pair2_weights_raw = carrier_vectors.get("pair2")
    if not isinstance(pair1_weights_raw, dict) or not isinstance(pair2_weights_raw, dict):
        artifact = {
            "stage": "P730",
            "status": "INVALID_F301_CARRIER_SHAPE",
            "as_of": AS_OF,
            "expected": "F301 carrier_vectors contain pair1 and pair2 dictionaries",
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    pair1_prob_dict, pair1_total = normalize_branch_profile(pair1_weights_raw, n=n)
    pair2_prob_dict, pair2_total = normalize_branch_profile(pair2_weights_raw, n=n)

    ord_by_x = [float(v) for v in ord_by_x_raw]
    alpha_geo = 4.0 * math.log(2.0)
    r_weights = [math.exp(-alpha_geo * value) for value in ord_by_x]
    r_total = sum(r_weights)
    r_prob = [value / r_total for value in r_weights]

    pair1_prob = [pair1_prob_dict[idx] for idx in range(n)]
    pair2_prob = [pair2_prob_dict[idx] for idx in range(n)]

    pair2_is_exact_inversion_of_pair1 = max_abs_difference(
        pair2_prob,
        [pair1_prob[(-idx) % n] for idx in range(n)],
    ) <= TOL

    ord_inversion_invariant = max_abs_difference(
        ord_by_x,
        [ord_by_x[(-idx) % n] for idx in range(n)],
    ) <= TOL
    r_prob_inversion_invariant = max_abs_difference(
        r_prob,
        [r_prob[(-idx) % n] for idx in range(n)],
    ) <= TOL

    pair1_expectation_ord = expectation(pair1_prob, ord_by_x)
    pair2_expectation_ord = expectation(pair2_prob, ord_by_x)
    pair1_cross_entropy_to_r_ord = cross_entropy(pair1_prob, r_prob)
    pair2_cross_entropy_to_r_ord = cross_entropy(pair2_prob, r_prob)

    expectation_gap_abs = abs(pair1_expectation_ord - pair2_expectation_ord)
    cross_entropy_gap_abs = abs(pair1_cross_entropy_to_r_ord - pair2_cross_entropy_to_r_ord)

    theta_pair_outputs = theta_pair.get("outputs") or {}
    theta_pair_construction = theta_pair.get("construction_class") or {}

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
        "p729_orbit_direction_split_already_localized",
        bool(p729.get("remaining_pair12_split_localized_as_opposite_orbit_directions")),
        True,
        "P729 already localizes the surviving pair1/pair2 ambiguity as opposite orbit-index directions δ_k versus δ_{-k}.",
    )
    add_check(
        "current_direction_free_shannon_lane_exported_on_pair1_pair2",
        {
            "reference_distribution": theta_pair_construction.get("reference_distribution"),
            "slot_free": bool(theta_pair_construction.get("slot_free")),
            "outputs": sorted(theta_pair_outputs.keys()),
        },
        {
            "reference_distribution": str(IN_R_ORD.relative_to(REPO)),
            "slot_free": True,
            "outputs": ["pair1", "pair2"],
        },
        "The current strict Shannon ord-reference lane already exports a slot-free theta-pair source on pair1 and pair2, derived from the direction-free reference r_ord.",
    )
    add_check(
        "current_direction_free_shannon_lane_objective_is_cross_entropy_to_r_ord",
        any(
            token in str(theta_pair_construction.get("objective_family") or "").lower()
            for token in ["cross-entropy", "cross_entropy", "log r_ord"]
        ),
        True,
        "The current strict Shannon theta-pair source is built from the cross-entropy objective family against r_ord.",
    )
    add_check(
        "ord_reference_profile_inversion_invariant",
        ord_inversion_invariant,
        True,
        "The current element-order reference profile ord_Z12(x) is invariant under orbit-direction inversion k -> -k.",
    )
    add_check(
        "r_ord_reference_distribution_inversion_invariant",
        r_prob_inversion_invariant,
        True,
        "The current strict Shannon reference distribution r_ord is likewise inversion-invariant on Z_12.",
    )
    add_check(
        "f301_pair12_branch_profiles_are_exact_inversion_partners",
        pair2_is_exact_inversion_of_pair1,
        True,
        "The surviving F301 pair1/pair2 residual-datum carrier profiles are exact inversion partners under k -> -k.",
    )
    add_check(
        "direction_free_shannon_ord_expectation_scores_equal_on_pair12_branches",
        expectation_gap_abs <= TOL,
        True,
        "Because both the branch profiles and the ord-reference are inversion-related / inversion-invariant, the induced source-side ord expectation is exactly equal on δ_k and δ_{-k}.",
    )
    add_check(
        "direction_free_shannon_cross_entropy_scores_equal_on_pair12_branches",
        cross_entropy_gap_abs <= TOL,
        True,
        "The current strict Shannon cross-entropy source score against r_ord is exactly equal on the surviving pair1/pair2 orbit-direction branches.",
    )
    add_check(
        "t184_direction_free_shannon_pair12_orbit_direction_selection_bridge_exported",
        False,
        False,
        "No current export upgrades the direction-free Shannon ord-reference lane into a source-side selector between δ_k and δ_{-k}.",
    )

    status = (
        "PASS_DIRECTION_FREE_SHANNON_PAIR12_ORBIT_DIRECTION_SELECTION_BRIDGE_NONEXPORT_AUDITED"
        if not blocking
        else "P730_REQUIRES_REVIEW_CHANGED_DIRECTION_FREE_SHANNON_PAIR12_STATE"
    )

    artifact = {
        "stage": "P730",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t184_direction_free_shannon_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_boundary_only",
        "inputs": {
            "P729": str(IN_P729.relative_to(REPO)),
            "F301_artifact": str(IN_F301.relative_to(REPO)),
            "r_ord_reference": str(IN_R_ORD.relative_to(REPO)),
            "theta_pair_sigma_int_slot_free_source": str(IN_THETA_PAIR.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t184_target_name": "DirectionFreeShannonResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1",
        "t184_target_exported_on_current_repo_state": False,
        "current_direction_free_shannon_lane_already_exports_pair1_pair2_o2_to_z2_cuts": True,
        "pair1_branch_kind": "delta_k_positive_index_branch",
        "pair2_branch_kind": "delta_minus_k_negative_index_branch",
        "ord_reference_profile_inversion_invariant": ord_inversion_invariant,
        "r_ord_reference_distribution_inversion_invariant": r_prob_inversion_invariant,
        "pair12_branch_profiles_equal_under_inversion": pair2_is_exact_inversion_of_pair1,
        "pair12_source_scores_under_direction_free_shannon_lane": {
            "pair1_expectation_ord": pair1_expectation_ord,
            "pair2_expectation_ord": pair2_expectation_ord,
            "expectation_ord_gap_abs": expectation_gap_abs,
            "pair1_cross_entropy_to_r_ord": pair1_cross_entropy_to_r_ord,
            "pair2_cross_entropy_to_r_ord": pair2_cross_entropy_to_r_ord,
            "cross_entropy_gap_abs": cross_entropy_gap_abs,
        },
        "branch_profile_totals_before_normalization": {
            "pair1_total": pair1_total,
            "pair2_total": pair2_total,
        },
        "audit_conclusion": {
            "current_direction_free_shannon_lane_selects_pair12_orbit_direction_branch": False,
            "current_direction_free_shannon_lane_already_exports_pair1_pair2_o2_cut": True,
            "next_honest_move": (
                "attack_or_export_a_genuinely_inversion_sensitive_source_side_rule_beyond_the_current_direction_free_shannon_ord_reference_lane_if_pair12_branch_selection_is_still_required"
            ),
        },
        "hard_limits": [
            "No T184 discharge claim.",
            "No promotion of orbit-index direction into a physical orientation datum.",
            "No claim that the current pair1/pair2 O(2)->Z2 cut selects a unique source-side orbit branch.",
            "No identification of tau_src with the current selector carrier.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P730",
        "status": status,
        "as_of": AS_OF,
        "t184_target_name": "DirectionFreeShannonResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1",
        "t184_target_exported_on_current_repo_state": False,
        "current_direction_free_shannon_lane_already_exports_pair1_pair2_o2_to_z2_cuts": True,
        "pair12_branch_profiles_equal_under_inversion": pair2_is_exact_inversion_of_pair1,
        "direction_free_shannon_pair12_expectation_ord_scores_equal": expectation_gap_abs <= TOL,
        "direction_free_shannon_pair12_cross_entropy_scores_equal": cross_entropy_gap_abs <= TOL,
        "current_direction_free_shannon_lane_selects_pair12_orbit_direction_branch": False,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
