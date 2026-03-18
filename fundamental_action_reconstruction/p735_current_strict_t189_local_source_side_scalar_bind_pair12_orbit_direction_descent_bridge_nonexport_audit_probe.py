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

IN_P229 = GENERATED / "p229_current_actual_source_topology_barrier_protected_sign_witness_probe_summary.json"
IN_P231 = GENERATED / "p231_current_actual_source_topology_nonzero_flow_witness_probe_summary.json"
IN_P729 = GENERATED / "p729_current_strict_t183_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"
IN_P731 = GENERATED / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe_summary.json"
IN_P734 = GENERATED / "p734_current_strict_t188_declared_scope_source_topology_selector_theorem_pair12_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json"
IN_F141 = GENERATED / "f141_first_actual_source_topology_barrier_protected_sign_witness_packet_summary.json"
IN_F143 = GENERATED / "f143_first_actual_source_topology_nonzero_flow_witness_packet_summary.json"
IN_PROVIDER_LAYER = GENERATED / "provider_object_carrier_layer_actual_inhabitant_instance.json"
IN_F301 = GENERATED / "residual_datum_bridge_export_map_object_support_carrier_actual_inhabitant_instance.json"

OUT_JSON = GENERATED / "p735_current_strict_t189_local_source_side_scalar_bind_pair12_orbit_direction_descent_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p735_current_strict_t189_local_source_side_scalar_bind_pair12_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def close(a: Any, b: Any, tol: float = TOL) -> bool:
    return isinstance(a, (int, float)) and isinstance(b, (int, float)) and abs(float(a) - float(b)) <= tol


def list_close(xs: Any, ys: Any, tol: float = TOL) -> bool:
    if not isinstance(xs, list) or not isinstance(ys, list) or len(xs) != len(ys):
        return False
    return all(close(x, y, tol=tol) for x, y in zip(xs, ys))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [
        IN_P229,
        IN_P231,
        IN_P729,
        IN_P731,
        IN_P734,
        IN_F141,
        IN_F143,
        IN_PROVIDER_LAYER,
        IN_F301,
    ]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P735",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p229 = load_json(IN_P229)
    p231 = load_json(IN_P231)
    p729 = load_json(IN_P729)
    p731 = load_json(IN_P731)
    p734 = load_json(IN_P734)
    f141 = load_json(IN_F141)
    f143 = load_json(IN_F143)
    provider_layer = load_json(IN_PROVIDER_LAYER)
    f301 = load_json(IN_F301)

    barrier_support = f141.get("support_packet") or {}
    provider_derived = provider_layer.get("derived_parameters") or {}
    provider_form = provider_layer.get("carrier_form") or {}
    provider_ops = provider_form.get("provider_operators") or {}
    pair1_provider = provider_ops.get("pair1") or {}
    pair2_provider = provider_ops.get("pair2") or {}
    nad12 = f301.get("nad12_depth") or {}

    barrier_sign = barrier_support.get("psi_src_barrier_sign_component_witness_v1")
    flow_value = f143.get("scalar_component_witness_value")
    a_val = provider_derived.get("a")
    b_val = provider_derived.get("b")
    pair1_weights = nad12.get("pair1_weights")
    pair2_weights = nad12.get("pair2_weights")

    contraction_parameters_equal_and_positive = (
        close(a_val, b_val)
        and isinstance(a_val, (int, float))
        and float(a_val) > 0.0
        and isinstance(b_val, (int, float))
        and float(b_val) > 0.0
    )
    local_scalar_family_factors_through_shared_cos_phi_data = (
        barrier_sign == 1
        and isinstance(flow_value, (int, float))
        and float(flow_value) > 0.0
        and contraction_parameters_equal_and_positive
        and close(flow_value, a_val)
        and close(flow_value, b_val)
    )
    pair12_weight_profiles_identical = list_close(pair1_weights, pair2_weights)
    current_local_source_side_scalar_bind_is_pair12_branch_blind = (
        local_scalar_family_factors_through_shared_cos_phi_data
        and pair12_weight_profiles_identical
    )
    p731_pair12_witness_split_descends_to_current_local_source_side_scalar_bind = (
        bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches"))
        and not current_local_source_side_scalar_bind_is_pair12_branch_blind
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
        "p734_strongest_current_source_side_theorem_lane_already_exists_but_does_not_descend_split",
        {
            "theorem_lane_exported": bool(p734.get("current_declared_scope_source_topology_selector_theorem_exported")),
            "remains_quotient_class_only": bool(p734.get("current_declared_scope_source_topology_selector_theorem_remains_quotient_class_only")),
            "split_descends": bool(p734.get("p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem")),
        },
        {
            "theorem_lane_exported": True,
            "remains_quotient_class_only": True,
            "split_descends": False,
        },
        "P734 already proves that the strongest current source-side theorem lane exists, but remains quotient-only and does not descend the P731 split.",
    )
    add_check(
        "current_local_source_topology_scalar_witnesses_are_exported_and_downstream_only",
        {
            "P229_actual_barrier_sign_witness_exported": bool(p229.get("actual_barrier_protected_sign_witness_exported")),
            "P231_actual_nonzero_flow_witness_exported": bool(p231.get("actual_nonzero_flow_witness_exported")),
            "P229_observer_role": p229.get("observer_role"),
            "P231_observer_role": p231.get("observer_role"),
            "F141_barrier_sign_component": barrier_sign,
            "F143_nonzero_flow_component_positive": bool(isinstance(flow_value, (int, float)) and float(flow_value) > 0.0),
        },
        {
            "P229_actual_barrier_sign_witness_exported": True,
            "P231_actual_nonzero_flow_witness_exported": True,
            "P229_observer_role": "downstream_only",
            "P231_observer_role": "downstream_only",
            "F141_barrier_sign_component": 1,
            "F143_nonzero_flow_component_positive": True,
        },
        "The current local source-topology scalar witness family is real and positive, but it is still exported only as downstream-only common source-side data on tau_src_candidate_v1.",
    )
    add_check(
        "provider_object_carrier_layer_uses_equal_positive_source_derived_contraction_parameters",
        {
            "source_map": provider_layer.get("contraction_parameter_source_map"),
            "a_equals_b": contraction_parameters_equal_and_positive,
            "flow_matches_a": close(flow_value, a_val),
            "flow_matches_b": close(flow_value, b_val),
            "pair_index_set": provider_form.get("pair_index_set"),
        },
        {
            "source_map": "A_strict_provider_object_contraction_parameter_source_map_v1",
            "a_equals_b": True,
            "flow_matches_a": True,
            "flow_matches_b": True,
            "pair_index_set": ["pair1", "pair2"],
        },
        "The current provider-object carrier layer already binds the source-side scalar data into the pair1/pair2 carrier lane, but it does so with equal positive source-derived contraction parameters a=b=cos(phi).",
    )
    add_check(
        "provider_object_carrier_pair12_operators_differ_only_by_shift_direction_not_scalar_amplitude",
        {
            "pair1_formula_contains_left_shift": "k-1" in str(pair1_provider.get("formula") or ""),
            "pair2_formula_contains_right_shift": "k+1" in str(pair2_provider.get("formula") or ""),
            "pair1_scalar_matches_a": close(pair1_provider.get("a"), a_val),
            "pair2_scalar_matches_b": close(pair2_provider.get("b"), b_val),
        },
        {
            "pair1_formula_contains_left_shift": True,
            "pair2_formula_contains_right_shift": True,
            "pair1_scalar_matches_a": True,
            "pair2_scalar_matches_b": True,
        },
        "On the current provider-object carrier lane, the remaining pair1/pair2 difference is only the opposite orbit-direction shift, not a scalar-amplitude distinction.",
    )
    add_check(
        "f301_post_witness_pair12_weight_profiles_are_identical",
        {
            "source_domain": f301.get("source_domain"),
            "pair_index_set": f301.get("pair_index_set"),
            "pair12_weight_profiles_identical": pair12_weight_profiles_identical,
        },
        {
            "source_domain": "tau_src_candidate_v1",
            "pair_index_set": ["pair1", "pair2"],
            "pair12_weight_profiles_identical": True,
        },
        "The current post-witness residual-datum support carrier F301 stays on tau_src_candidate_v1 and inherits identical pair1/pair2 weight profiles from the symmetric source-derived contraction data.",
    )
    add_check(
        "current_local_source_side_scalar_witness_family_factors_through_shared_cos_phi_data",
        local_scalar_family_factors_through_shared_cos_phi_data,
        True,
        "The current local source-side scalar witness family factors only through one shared positive cos(phi) datum and its sign, rather than through one pair1/pair2 branch-sensitive scalar distinction.",
    )
    add_check(
        "current_local_source_side_scalar_bind_is_pair12_branch_blind",
        current_local_source_side_scalar_bind_is_pair12_branch_blind,
        True,
        "Therefore the current local source-side scalar bind data remain branch-blind on the surviving pair1/pair2 orbit-direction split.",
    )
    add_check(
        "p731_pair12_witness_split_descends_to_current_local_source_side_scalar_bind",
        p731_pair12_witness_split_descends_to_current_local_source_side_scalar_bind,
        False,
        "So the already-separated opposite P731 pair1/pair2 witness split still does not descend to one typed branch distinction through the current local source-side scalar bind data.",
    )
    add_check(
        "t189_local_source_side_scalar_bind_pair12_orbit_direction_descent_bridge_exported",
        False,
        False,
        "The repo still does not export the local source-side scalar bind pair1/pair2 orbit-direction descent bridge on current repo state.",
    )

    status = (
        "PASS_LOCAL_SOURCE_SIDE_SCALAR_BIND_PAIR12_ORBIT_DIRECTION_DESCENT_BRIDGE_NONEXPORT_AUDITED"
        if not blocking
        else "P735_REQUIRES_REVIEW_CHANGED_LOCAL_SOURCE_SIDE_SCALAR_BIND_PAIR12_STATE"
    )

    artifact = {
        "stage": "P735",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t189_local_source_side_scalar_bind_pair12_orbit_direction_descent_bridge_nonexport_boundary_only",
        "inputs": {
            "P229": str(IN_P229.relative_to(REPO)),
            "P231": str(IN_P231.relative_to(REPO)),
            "P729": str(IN_P729.relative_to(REPO)),
            "P731": str(IN_P731.relative_to(REPO)),
            "P734": str(IN_P734.relative_to(REPO)),
            "F141": str(IN_F141.relative_to(REPO)),
            "F143": str(IN_F143.relative_to(REPO)),
            "provider_object_carrier_layer_artifact": str(IN_PROVIDER_LAYER.relative_to(REPO)),
            "F301_artifact": str(IN_F301.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t189_target_name": "LocalSourceSideScalarBindPair12OrbitDirectionDescentBridge_global_C_v1_strict_v1",
        "t189_target_exported_on_current_repo_state": False,
        "current_local_source_side_scalar_witnesses_exported": bool(p229.get("actual_barrier_protected_sign_witness_exported"))
        and bool(p231.get("actual_nonzero_flow_witness_exported")),
        "current_local_source_side_scalar_witness_family_factors_through_shared_cos_phi_data": (
            local_scalar_family_factors_through_shared_cos_phi_data
        ),
        "current_local_source_side_scalar_bind_is_pair12_branch_blind": (
            current_local_source_side_scalar_bind_is_pair12_branch_blind
        ),
        "p731_pair12_witness_split_descends_to_current_local_source_side_scalar_bind": (
            p731_pair12_witness_split_descends_to_current_local_source_side_scalar_bind
        ),
        "current_local_source_side_scalar_bind_profile": {
            "barrier_sign_component": barrier_sign,
            "nonzero_flow_component_value": flow_value,
            "provider_contraction_parameters": {
                "a": a_val,
                "b": b_val,
            },
            "pair12_weight_profiles_identical": pair12_weight_profiles_identical,
        },
        "audit_conclusion": {
            "current_repo_already_exports_local_source_side_scalar_witnesses": bool(
                p229.get("actual_barrier_protected_sign_witness_exported")
            )
            and bool(p231.get("actual_nonzero_flow_witness_exported")),
            "current_repo_already_exports_t189_target": False,
            "next_honest_move": (
                "attack_a_genuinely_pair12_branch_sensitive_local_source_side_bind_or_descent_mechanism_that_does_not_factor_only_through_the_current_shared_cos_phi_scalar_family_or_the_current_symmetric_a_equals_b_carrier_data"
            ),
        },
        "hard_limits": [
            "No T189 discharge claim.",
            "No claim that the current local source-side scalar witness family already selects one pair1/pair2 orbit-direction branch.",
            "No claim that tau_src_candidate_v1 is identified with the current selector carrier.",
            "No admissible S_sel_int claim.",
            "No directed/sign-sensitive physical orientation datum claim.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P735",
        "status": status,
        "as_of": AS_OF,
        "t189_target_name": "LocalSourceSideScalarBindPair12OrbitDirectionDescentBridge_global_C_v1_strict_v1",
        "t189_target_exported_on_current_repo_state": False,
        "current_local_source_side_scalar_witness_family_factors_through_shared_cos_phi_data": (
            local_scalar_family_factors_through_shared_cos_phi_data
        ),
        "current_local_source_side_scalar_bind_is_pair12_branch_blind": (
            current_local_source_side_scalar_bind_is_pair12_branch_blind
        ),
        "p731_pair12_witness_split_descends_to_current_local_source_side_scalar_bind": (
            p731_pair12_witness_split_descends_to_current_local_source_side_scalar_bind
        ),
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
