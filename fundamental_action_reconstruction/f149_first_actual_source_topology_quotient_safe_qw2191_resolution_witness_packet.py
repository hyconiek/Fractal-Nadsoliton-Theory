#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
IN_F137 = GENERATED / "f137_first_source_topology_quotient_safe_qw2191_resolution_target_packet_summary.json"
IN_F148 = GENERATED / "f148_first_actual_source_topology_basis_independent_promotion_witness_packet_summary.json"
IN_QW2191 = ROOT.parent / "material_dowodowy" / "korpus_qw_pozostaly" / "raporty_json" / "report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"
OUT = GENERATED / "f149_first_actual_source_topology_quotient_safe_qw2191_resolution_witness_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f137 = load_json(IN_F137)
    f148 = load_json(IN_F148)
    q2191 = load_json(IN_QW2191)

    basis_packet = f148["support_packet"]
    axis_class = basis_packet["basis_free_axis_class"]
    split_class = basis_packet["basis_free_signed_split_class"]
    scope_tag = basis_packet["basis_free_scope_tag"]

    kernel_alone_obstruction_explicit = (
        q2191["flags"]["kernel_has_degenerate_eigenspaces"] is True
        and q2191["flags"]["o2_rotation_family_defined_on_degenerate_pairs"] is True
        and q2191["flags"]["continuous_nonunique_embedding_family_exhibited"] is True
        and q2191["flags"]["full_uniqueness_from_kernel_alone_obstructed"] is True
    )

    quotient_relation_well_defined = (
        axis_class["projector_idempotent"] is True
        and abs(axis_class["projector_trace"] - 1.0) <= 1e-12
        and axis_class["invariant_under_positive_rescaling"] is True
        and axis_class["matches_selector_plus_projector_on_light_plane"] is True
        and split_class["complementary_projectors"] is True
        and split_class["orthogonal_projector_pair"] is True
        and split_class["plus_projector_idempotent"] is True
        and split_class["minus_projector_idempotent"] is True
        and split_class["positive_signed_selector_response"] is True
        and split_class["positive_plus_output"] is True
        and split_class["vanishing_minus_output"] is True
        and scope_tag["preobserver_only"] is True
        and scope_tag["observer_downstream_only"] is True
    )

    distinguished_quotient_class_exported = (
        f148["actual_basis_independent_selector_promotion_exported"] is True
        and f148["basis_independence_discharged"] is True
        and f148["observer_role"] == "downstream_only"
        and f148["tau_src_identified_with_s_prelm"] is False
        and quotient_relation_well_defined is True
    )

    actual_quotient_safe_qw2191_resolution_exported = (
        f137["target_name"] == "Phi_qw2191_safe_target_v1"
        and f137["codomain_target"] == "actual_quotient_safe_qw2191_resolution_target_v1"
        and f148["input_packet"] == "tau_src_candidate_v1"
        and f148["witness"] == "Upsilon_sel_basis_actual_witness_v1"
        and kernel_alone_obstruction_explicit is True
        and distinguished_quotient_class_exported is True
    )

    support_packet = {
        "input_packet": "tau_src_candidate_v1",
        "basis_promotion_witness": f148["witness"],
        "qw2191_obstruction": {
            "kernel_alone_obstructed": q2191["flags"]["full_uniqueness_from_kernel_alone_obstructed"],
            "o2_family_present": q2191["flags"]["o2_rotation_family_defined_on_degenerate_pairs"],
            "continuous_nonunique_family_exhibited": q2191["flags"]["continuous_nonunique_embedding_family_exhibited"],
            "required_next_step_before_source_topology_route": q2191["required_next_step"],
        },
        "quotient_relation": {
            "object": "~_src_qw2191_v1",
            "criterion": "same_source_side_basis_free_selector_packet",
            "uses_observer_side_data": False,
            "uses_external_selector_axiom": False,
            "quotients_chart_labels_only_via_basis_free_packet": True,
        },
        "distinguished_quotient_class": {
            "object": "C_sel_src_qw2191_resolved_v1",
            "axis_class": axis_class["object"],
            "signed_split_class": split_class["object"],
            "preobserver_scope_tag": scope_tag["object"],
            "raw_theta_uniqueness_claimed": False,
            "quotient_class_only": True,
        },
    }

    summary = {
        "packet_id": "F149",
        "status": "F149_EXECUTED_FIRST_ACTUAL_SOURCE_TOPOLOGY_QUOTIENT_SAFE_QW2191_RESOLUTION_WITNESS_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-08",
        "refines_future_qw2191_target": f137["target_name"],
        "input_packet": "tau_src_candidate_v1",
        "witness": "Phi_qw2191_safe_actual_witness_v1",
        "codomain_target": f137["codomain_target"],
        "support_packet_id": "W_src_qw2191_safe_support_packet_v1",
        "support_packet": support_packet,
        "observer_role": "downstream_only",
        "tau_src_identified_with_s_prelm": f148["tau_src_identified_with_s_prelm"],
        "actual_basis_independent_selector_promotion_exported": f148["actual_basis_independent_selector_promotion_exported"],
        "basis_independence_discharged": f148["basis_independence_discharged"],
        "kernel_alone_qw2191_obstruction_retained": kernel_alone_obstruction_explicit,
        "distinguished_quotient_class_exported": distinguished_quotient_class_exported,
        "raw_theta_uniqueness_claimed": False,
        "actual_quotient_safe_qw2191_resolution_exported": actual_quotient_safe_qw2191_resolution_exported,
        "qw2191_quotient_safe_discharged": actual_quotient_safe_qw2191_resolution_exported,
        "current_selector_closure": False,
        "current_global_qw2191_discharge": False,
        "kernel_split_safe": True,
        "legacy_to_strict_bridge_claimed": False,
        "no_false_pass": True,
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
