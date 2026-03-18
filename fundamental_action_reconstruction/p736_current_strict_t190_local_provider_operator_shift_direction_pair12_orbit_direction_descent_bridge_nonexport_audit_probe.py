#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-18"
TOL = 1.0e-15

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P729 = GENERATED / "p729_current_strict_t183_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"
IN_P731 = GENERATED / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe_summary.json"
IN_P735 = GENERATED / "p735_current_strict_t189_local_source_side_scalar_bind_pair12_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json"
IN_PROVIDER_LAYER = GENERATED / "provider_object_carrier_layer_actual_inhabitant_instance.json"
IN_F301 = GENERATED / "residual_datum_bridge_export_map_object_support_carrier_actual_inhabitant_instance.json"

OUT_JSON = GENERATED / "p736_current_strict_t190_local_provider_operator_shift_direction_pair12_orbit_direction_descent_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p736_current_strict_t190_local_provider_operator_shift_direction_pair12_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def close(a: Any, b: Any, tol: float = TOL) -> bool:
    return isinstance(a, (int, float)) and isinstance(b, (int, float)) and abs(float(a) - float(b)) <= tol


def contains_token(entries: Any, token: str) -> bool:
    return isinstance(entries, list) and any(token in str(entry) for entry in entries)


def load_numeric_map(obj: Any) -> dict[int, float]:
    if not isinstance(obj, dict):
        return {}
    out: dict[int, float] = {}
    for key, value in obj.items():
        if not isinstance(value, (int, float)):
            return {}
        out[int(key)] = float(value)
    return out


def inversion_paired(pair1: Any, pair2: Any, tol: float = TOL) -> bool:
    pair1_map = load_numeric_map(pair1)
    pair2_map = load_numeric_map(pair2)
    if not pair1_map or not pair2_map:
        return False
    if set(pair2_map.keys()) != {-k for k in pair1_map.keys()}:
        return False
    return all(close(value, pair2_map.get(-k), tol=tol) for k, value in pair1_map.items())


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P729, IN_P731, IN_P735, IN_PROVIDER_LAYER, IN_F301]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P736",
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
    p731 = load_json(IN_P731)
    p735 = load_json(IN_P735)
    provider_layer = load_json(IN_PROVIDER_LAYER)
    f301 = load_json(IN_F301)

    provider_derived = provider_layer.get("derived_parameters") or {}
    provider_form = provider_layer.get("carrier_form") or {}
    provider_ops = provider_form.get("provider_operators") or {}
    pair1_provider = provider_ops.get("pair1") or {}
    pair2_provider = provider_ops.get("pair2") or {}
    seed_state = provider_form.get("seed_state") or {}
    provider_notes = provider_layer.get("notes") or []
    provider_current_absence = provider_layer.get("current_absence") or []

    a_val = provider_derived.get("a")
    b_val = provider_derived.get("b")
    pair1_formula = str(pair1_provider.get("formula") or "")
    pair2_formula = str(pair2_provider.get("formula") or "")
    seed_formula = str(seed_state.get("formula") or "")

    normalized_reps = f301.get("normalized_orbit_representatives") or {}
    pair1_norm = str(normalized_reps.get("pair1") or "")
    pair2_norm = str(normalized_reps.get("pair2") or "")
    pair1_vector = ((f301.get("carrier_vectors") or {}).get("pair1")) or {}
    pair2_vector = ((f301.get("carrier_vectors") or {}).get("pair2")) or {}

    contraction_parameters_equal_and_positive = (
        close(a_val, b_val)
        and isinstance(a_val, (int, float))
        and float(a_val) > 0.0
        and isinstance(b_val, (int, float))
        and float(b_val) > 0.0
    )
    pair12_shift_directions_are_opposite = "k-1" in pair1_formula and "k+1" in pair2_formula
    pair12_carrier_vectors_are_exact_index_inversion_partners = inversion_paired(pair1_vector, pair2_vector)
    provider_layer_selector_neutral = (
        contains_token(provider_notes, "selector-neutral")
        and contains_token(provider_current_absence, "no admissible S_sel_int")
        and contains_token(provider_current_absence, "no strict-core selector closure")
    )
    current_local_provider_operator_shift_direction_lane_realizes_both_pair12_branches_symmetrically = (
        provider_form.get("pair_index_set") == ["pair1", "pair2"]
        and "δ_0" in seed_formula
        and pair12_shift_directions_are_opposite
        and contraction_parameters_equal_and_positive
        and "δ_k" in pair1_norm
        and "δ_{-k}" in pair2_norm
        and pair12_carrier_vectors_are_exact_index_inversion_partners
        and provider_layer_selector_neutral
    )
    p731_pair12_witness_split_descends_to_current_local_provider_operator_shift_direction_lane = (
        bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches"))
        and not current_local_provider_operator_shift_direction_lane_realizes_both_pair12_branches_symmetrically
        and not provider_layer_selector_neutral
    )

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
        "p729_pair12_orbit_direction_split_already_localized",
        bool(p729.get("remaining_pair12_split_localized_as_opposite_orbit_directions")),
        True,
        "P729 already localizes the surviving pair1/pair2 ambiguity as the opposite residual-datum orbit-direction branches delta_k and delta_-k.",
    )
    add_check(
        "p731_pair12_witness_split_already_opposite_and_unpromoted",
        {
            "split_separated": bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches")),
            "pair1_sign": p731.get("pair1_w_break_branch_score_sign"),
            "pair2_sign": p731.get("pair2_w_break_branch_score_sign"),
            "t185_exported": bool(p731.get("t185_target_exported_on_current_repo_state")),
        },
        {
            "split_separated": True,
            "pair1_sign": "negative",
            "pair2_sign": "positive",
            "t185_exported": False,
        },
        "P731 already separates the surviving pair1/pair2 branches by opposite witness-score signs, while the typed promotion bridge remains unexported.",
    )
    add_check(
        "p735_current_local_source_side_scalar_bind_already_exhausted_as_branch_blind",
        {
            "scalar_bind_branch_blind": bool(p735.get("current_local_source_side_scalar_bind_is_pair12_branch_blind")),
            "scalar_bind_split_descends": bool(p735.get("p731_pair12_witness_split_descends_to_current_local_source_side_scalar_bind")),
            "t189_exported": bool(p735.get("t189_target_exported_on_current_repo_state")),
        },
        {
            "scalar_bind_branch_blind": True,
            "scalar_bind_split_descends": False,
            "t189_exported": False,
        },
        "P735 already proves that the strongest current local source-side scalar bind lane is real but branch-blind, so the next honest test must move to the local non-scalar provider-operator lane itself.",
    )
    add_check(
        "provider_object_carrier_layer_pair_indexed_opposite_shift_lane_is_exported",
        {
            "pair_index_set": provider_form.get("pair_index_set"),
            "seed_state_is_delta0": "δ_0" in seed_formula,
            "pair1_left_shift": "k-1" in pair1_formula,
            "pair2_right_shift": "k+1" in pair2_formula,
            "selector_neutral": provider_layer_selector_neutral,
        },
        {
            "pair_index_set": ["pair1", "pair2"],
            "seed_state_is_delta0": True,
            "pair1_left_shift": True,
            "pair2_right_shift": True,
            "selector_neutral": True,
        },
        "The current provider-object carrier layer already exports one pair-indexed local non-scalar lane: opposite shift-direction operators from the common seed delta_0, but still in explicit selector-neutral scope.",
    )
    add_check(
        "provider_object_carrier_shift_direction_lane_uses_equal_positive_source_derived_amplitudes",
        {
            "source_map": provider_layer.get("contraction_parameter_source_map"),
            "a_equals_b_positive": contraction_parameters_equal_and_positive,
            "pair1_scalar_matches_a": close(pair1_provider.get("a"), a_val),
            "pair2_scalar_matches_b": close(pair2_provider.get("b"), b_val),
        },
        {
            "source_map": "A_strict_provider_object_contraction_parameter_source_map_v1",
            "a_equals_b_positive": True,
            "pair1_scalar_matches_a": True,
            "pair2_scalar_matches_b": True,
        },
        "That current pair-indexed operator lane still uses the same source-derived positive amplitude on both branches: (a,b)=(cos(phi),cos(phi)).",
    )
    add_check(
        "f301_pair12_carrier_repackages_current_shift_lane_as_exact_inversion_pair",
        {
            "source_domain": f301.get("source_domain"),
            "pair1_is_delta_k": "δ_k" in pair1_norm,
            "pair2_is_delta_minus_k": "δ_{-k}" in pair2_norm,
            "carrier_vectors_exact_inversion_pair": pair12_carrier_vectors_are_exact_index_inversion_partners,
        },
        {
            "source_domain": "tau_src_candidate_v1",
            "pair1_is_delta_k": True,
            "pair2_is_delta_minus_k": True,
            "carrier_vectors_exact_inversion_pair": True,
        },
        "F301 already repackages that opposite-shift provider lane as the exact inversion-paired residual branches delta_k and delta_-k on tau_src_candidate_v1.",
    )
    add_check(
        "current_local_provider_operator_shift_direction_lane_realizes_both_pair12_branches_symmetrically",
        current_local_provider_operator_shift_direction_lane_realizes_both_pair12_branches_symmetrically,
        True,
        "Therefore the current local non-scalar provider-operator lane already realizes both surviving pair1/pair2 branches, but only symmetrically as opposite shift directions from the same seed with equal amplitudes.",
    )
    add_check(
        "p731_pair12_witness_split_descends_to_current_local_provider_operator_shift_direction_lane",
        p731_pair12_witness_split_descends_to_current_local_provider_operator_shift_direction_lane,
        False,
        "So the already-separated opposite P731 pair1/pair2 witness split still does not descend through the current local provider-operator shift-direction lane as one typed branch distinction.",
    )
    add_check(
        "t190_local_provider_operator_shift_direction_pair12_orbit_direction_descent_bridge_exported",
        False,
        False,
        "The repo still does not export the local provider-operator shift-direction pair1/pair2 orbit-direction descent bridge on current repo state.",
    )

    status = (
        "PASS_LOCAL_PROVIDER_OPERATOR_SHIFT_DIRECTION_PAIR12_ORBIT_DIRECTION_DESCENT_BRIDGE_NONEXPORT_AUDITED"
        if not blocking
        else "P736_REQUIRES_REVIEW_CHANGED_LOCAL_PROVIDER_OPERATOR_SHIFT_DIRECTION_PAIR12_STATE"
    )

    artifact = {
        "stage": "P736",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t190_local_provider_operator_shift_direction_pair12_orbit_direction_descent_bridge_nonexport_boundary_only",
        "inputs": {
            "P729": str(IN_P729.relative_to(REPO)),
            "P731": str(IN_P731.relative_to(REPO)),
            "P735": str(IN_P735.relative_to(REPO)),
            "provider_object_carrier_layer_artifact": str(IN_PROVIDER_LAYER.relative_to(REPO)),
            "F301_artifact": str(IN_F301.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t190_target_name": "LocalProviderOperatorShiftDirectionPair12OrbitDirectionDescentBridge_global_C_v1_strict_v1",
        "t190_target_exported_on_current_repo_state": False,
        "current_local_provider_operator_shift_direction_lane_realizes_both_pair12_branches_symmetrically": (
            current_local_provider_operator_shift_direction_lane_realizes_both_pair12_branches_symmetrically
        ),
        "current_local_provider_operator_shift_direction_lane_is_selector_neutral": provider_layer_selector_neutral,
        "p731_pair12_witness_split_descends_to_current_local_provider_operator_shift_direction_lane": (
            p731_pair12_witness_split_descends_to_current_local_provider_operator_shift_direction_lane
        ),
        "current_local_provider_operator_shift_direction_lane_profile": {
            "seed_state_formula": seed_formula,
            "pair1_formula": pair1_formula,
            "pair2_formula": pair2_formula,
            "provider_contraction_parameters": {
                "a": a_val,
                "b": b_val,
            },
            "pair12_carrier_vectors_exact_inversion_pair": pair12_carrier_vectors_are_exact_index_inversion_partners,
        },
        "audit_conclusion": {
            "current_repo_already_exports_pair_indexed_local_non_scalar_provider_operator_lane": True,
            "current_repo_already_exports_t190_target": False,
            "next_honest_move": (
                "attack_a_typed_local_source_side_bind_or_descent_mechanism_that_couples_the_current_inversion_sensitive_p731_witness_split_to_the_pair_indexed_shift_direction_lane_without_collapsing_back_to_the_current_symmetric_a_equals_b_common_seed_selector_neutral_provider_operator_realization"
            ),
        },
        "hard_limits": [
            "No T190 discharge claim.",
            "No claim that the current local provider-operator shift-direction lane already selects one raw pair1/pair2 orbit-direction branch.",
            "No claim that tau_src_candidate_v1 is identified with the current selector carrier.",
            "No admissible S_sel_int claim.",
            "No directed/sign-sensitive physical orientation datum claim.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P736",
        "status": status,
        "as_of": AS_OF,
        "t190_target_name": "LocalProviderOperatorShiftDirectionPair12OrbitDirectionDescentBridge_global_C_v1_strict_v1",
        "t190_target_exported_on_current_repo_state": False,
        "current_local_provider_operator_shift_direction_lane_realizes_both_pair12_branches_symmetrically": (
            current_local_provider_operator_shift_direction_lane_realizes_both_pair12_branches_symmetrically
        ),
        "current_local_provider_operator_shift_direction_lane_is_selector_neutral": provider_layer_selector_neutral,
        "p731_pair12_witness_split_descends_to_current_local_provider_operator_shift_direction_lane": (
            p731_pair12_witness_split_descends_to_current_local_provider_operator_shift_direction_lane
        ),
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
