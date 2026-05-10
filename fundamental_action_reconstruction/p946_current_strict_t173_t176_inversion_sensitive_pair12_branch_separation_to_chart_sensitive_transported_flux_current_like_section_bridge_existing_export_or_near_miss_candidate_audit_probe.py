#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-21"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P945 = GENERATED / "p945_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_nonexport_audit_probe.json"
IN_F945 = GENERATED / "f945_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_target_packet.json"
IN_P712 = GENERATED / "p712_current_strict_t176_existing_global_directed_sign_coherence_provider_nonexport_audit_probe_summary.json"
IN_P718 = GENERATED / "p718_current_strict_t176_single_mixed_linear_weight_span_provider_insufficiency_audit_probe_summary.json"
IN_P719 = GENERATED / "p719_current_strict_t176_low_complexity_odd_polynomial_two_readout_provider_class_audit_probe_summary.json"
IN_P720 = GENERATED / "p720_current_strict_t176_observer_facing_signed_output_channel_projection_provider_class_audit_probe_summary.json"
IN_P721 = GENERATED / "p721_current_strict_t176_source_topology_basis_free_qw2191_safe_provider_nonupgrade_audit_probe_summary.json"
IN_P722 = GENERATED / "p722_current_strict_t177_chart_sensitive_transported_flux_current_like_section_nonexport_audit_probe_summary.json"
IN_P734 = GENERATED / "p734_current_strict_t188_declared_scope_source_topology_selector_theorem_pair12_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json"
IN_P742 = GENERATED / "p742_current_strict_t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P747 = GENERATED / "p747_current_strict_t201_actual_source_topology_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge_nonexport_audit_probe_summary.json"
IN_P763 = GENERATED / "p763_current_strict_t217_pair12_source_side_branch_selection_provider_actual_realization_attempt_immediate_missing_interface_nonexport_audit_probe_summary.json"
IN_P764 = GENERATED / "p764_current_strict_t218_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_probe_summary.json"
IN_P765 = GENERATED / "p765_current_strict_t219_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_nonexport_audit_probe_summary.json"
IN_P766 = GENERATED / "p766_current_strict_t220_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_probe_summary.json"
IN_P767 = GENERATED / "p767_current_strict_t221_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_nonexport_audit_probe_summary.json"
IN_P768 = GENERATED / "p768_current_strict_t222_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_target_probe_summary.json"
IN_P769 = GENERATED / "p769_current_strict_t223_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_nonexport_audit_probe_summary.json"
IN_P770 = GENERATED / "p770_current_strict_t224_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_attempt_probe_summary.json"
IN_P683 = GENERATED / "p683_current_strict_t173_rooted_transport_sign_lift_from_s_sel_int_w_break_audit_probe_summary.json"
IN_P684 = GENERATED / "p684_current_strict_t173_w_break_rooted_directed_state_lift_audit_probe_summary.json"
IN_F684 = GENERATED / "f684_first_exported_selector_state_global_c_v1_directed_rooted_transport_from_s_sel_int_w_break_packet_summary.json"
IN_F685 = GENERATED / "f685_first_exported_selector_closure_global_c_v1_directed_rooted_transport_from_s_sel_int_w_break_packet_summary.json"

OUT_JSON = GENERATED / "p946_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_existing_export_or_near_miss_candidate_audit_probe.json"
OUT_SUMMARY = GENERATED / "p946_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_existing_export_or_near_miss_candidate_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [
        IN_P945,
        IN_F945,
        IN_P712,
        IN_P718,
        IN_P719,
        IN_P720,
        IN_P721,
        IN_P722,
        IN_P734,
        IN_P742,
        IN_P747,
        IN_P763,
        IN_P764,
        IN_P765,
        IN_P766,
        IN_P767,
        IN_P768,
        IN_P769,
        IN_P770,
        IN_P683,
        IN_P684,
        IN_F684,
        IN_F685,
    ]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P946",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p945 = load_json(IN_P945)
    f945 = load_json(IN_F945)
    p712 = load_json(IN_P712)
    p718 = load_json(IN_P718)
    p719 = load_json(IN_P719)
    p720 = load_json(IN_P720)
    p721 = load_json(IN_P721)
    p722 = load_json(IN_P722)
    p734 = load_json(IN_P734)
    p742 = load_json(IN_P742)
    p747 = load_json(IN_P747)
    p763 = load_json(IN_P763)
    p764 = load_json(IN_P764)
    p765 = load_json(IN_P765)
    p766 = load_json(IN_P766)
    p767 = load_json(IN_P767)
    p768 = load_json(IN_P768)
    p769 = load_json(IN_P769)
    p770 = load_json(IN_P770)
    p683 = load_json(IN_P683)
    p684 = load_json(IN_P684)
    f684 = load_json(IN_F684)
    f685 = load_json(IN_F685)

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

    f945_target = f945.get("target_object") or {}
    f945_schema = f945.get("bridge_output_schema") or {}
    p945_conclusion = p945.get("audit_conclusion") or {}

    rooted_transport_family_is_real_but_nonlawful = (
        bool(p683.get("all_pairs_supported"))
        and bool(p683.get("directed_output_sign_consistent_if_all_pairs_supported"))
        and bool(p684.get("descends_to_projective"))
        and not bool(p684.get("counts_as_strict_physical_orientation_datum"))
        and bool(f685.get("glued"))
        and not bool(f684.get("counts_as_strict_physical_orientation_datum"))
    )

    static_provider_classes_still_have_no_exact_candidate = (
        not bool(p718.get("single_mixed_linear_weight_span_exact_root_independent_section_exists"))
        and int(p719.get("exact_candidates_found", -1)) == 0
        and int(p720.get("exact_candidates_found", -1)) == 0
        and not bool(p712.get("t176_target_exported_on_current_repo_state"))
    )

    source_topology_bridge_family_still_below_required_level = (
        not bool(p734.get("t188_target_exported_on_current_repo_state"))
        and bool(p734.get("current_declared_scope_source_topology_selector_theorem_remains_quotient_class_only"))
        and not bool(p742.get("t196_target_exported_on_current_repo_state"))
        and bool(p742.get("current_actual_selector_witness_codomain_continuation_remains_basis_free_not_pair12_typed"))
        and not bool(p747.get("t201_target_exported_on_current_repo_state"))
        and bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe"))
        and bool(p747.get("current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas"))
    )

    exact_t216_t224_near_miss_family_is_frozen_but_still_nonexport = (
        bool(p763.get("current_t216_attempt_immediate_missing_interface_is_still_unexported"))
        and bool(p764.get("current_t218_target_is_future_route_only"))
        and not bool(p765.get("t219_target_exported_on_current_repo_state"))
        and bool(p766.get("t220_attempt_exported_on_current_repo_state"))
        and bool(p767.get("current_repo_still_does_not_export_actual_realization_of_t220_immediate_missing_subinterface"))
        and bool(p768.get("current_t222_target_is_future_route_only"))
        and not bool(p769.get("t223_target_exported_on_current_repo_state"))
        and bool(p770.get("t224_attempt_exported_on_current_repo_state"))
    )

    no_current_lawful_bridge_supplier_found = (
        rooted_transport_family_is_real_but_nonlawful
        and static_provider_classes_still_have_no_exact_candidate
        and source_topology_bridge_family_still_below_required_level
        and exact_t216_t224_near_miss_family_is_frozen_but_still_nonexport
        and not bool(p722.get("t177_target_exported_on_current_repo_state"))
        and not bool(p945.get("bridge_target_exported_on_current_repo_state"))
    )

    add_check(
        "f945_still_freezes_exact_missing_bridge_output_schema_under_t176_t177_nonexport",
        {
            "p945_bridge_target_exported": bool(p945.get("bridge_target_exported_on_current_repo_state")),
            "p945_pair12_to_t176_bridge_exported": bool(
                p945_conclusion.get("current_repo_already_exports_inversion_sensitive_pair12_to_t176_bridge")
            ),
            "f945_target_object_id": f945_target.get("object_id"),
            "f945_bridge_output_schema_minimum_properties": f945_schema.get("minimum_properties"),
            "p722_t177_exported": bool(p722.get("t177_target_exported_on_current_repo_state")),
        },
        {
            "p945_bridge_target_exported": False,
            "p945_pair12_to_t176_bridge_exported": False,
            "f945_target_object_id": "inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_global_C_v1_strict_v1",
            "f945_bridge_output_schema_minimum_properties": [
                "full_C_v1_scope",
                "pair12_branch_sensitive",
                "chart_sensitive_transported_section_level",
                "nonconvention_nonprojective_nonpremise_smuggled",
            ],
            "p722_t177_exported": False,
        },
        "F945 still freezes one exact missing bridge object together with its minimum schema, and T177 remains unexported.",
    )
    add_check(
        "rooted_transport_family_is_global_but_convention_projective_and_nonphysical",
        {
            "all_pairs_supported": bool(p683.get("all_pairs_supported")),
            "directed_output_sign_consistent": bool(p683.get("directed_output_sign_consistent_if_all_pairs_supported")),
            "descends_to_projective": bool(p684.get("descends_to_projective")),
            "p684_counts_as_strict_physical_orientation_datum": bool(
                p684.get("counts_as_strict_physical_orientation_datum")
            ),
            "f685_glued": bool(f685.get("glued")),
            "f684_counts_as_strict_physical_orientation_datum": bool(
                f684.get("counts_as_strict_physical_orientation_datum")
            ),
            "can_lawfully_supply_f945_bridge": False,
        },
        {
            "all_pairs_supported": True,
            "directed_output_sign_consistent": True,
            "descends_to_projective": True,
            "p684_counts_as_strict_physical_orientation_datum": False,
            "f685_glued": True,
            "f684_counts_as_strict_physical_orientation_datum": False,
            "can_lawfully_supply_f945_bridge": False,
        },
        "The rooted transport family is real and global, but it remains projective/convention-layer and therefore cannot lawfully supply the F945 bridge.",
    )
    add_check(
        "static_t176_provider_classes_still_have_zero_exact_candidates",
        {
            "p718_exact_exists": bool(p718.get("single_mixed_linear_weight_span_exact_root_independent_section_exists")),
            "p719_exact_candidates_found": int(p719.get("exact_candidates_found", -1)),
            "p720_exact_candidates_found": int(p720.get("exact_candidates_found", -1)),
            "p712_t176_exported": bool(p712.get("t176_target_exported_on_current_repo_state")),
        },
        {
            "p718_exact_exists": False,
            "p719_exact_candidates_found": 0,
            "p720_exact_candidates_found": 0,
            "p712_t176_exported": False,
        },
        "The currently scanned static provider classes still contain no exact T176 provider candidate.",
    )
    add_check(
        "source_topology_theorem_witness_bridge_family_remains_quotient_or_projector_only",
        {
            "p721_source_topology_lane_upgrades_to_exact_t176_provider": bool(
                p721.get("source_topology_lane_upgrades_to_exact_t176_provider")
            ),
            "p734_t188_exported": bool(p734.get("t188_target_exported_on_current_repo_state")),
            "p734_quotient_only": bool(p734.get("current_declared_scope_source_topology_selector_theorem_remains_quotient_class_only")),
            "p742_t196_exported": bool(p742.get("t196_target_exported_on_current_repo_state")),
            "p742_pair12_typed_continuation_exported": bool(
                p742.get("current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation")
            ),
            "p747_t201_exported": bool(p747.get("t201_target_exported_on_current_repo_state")),
            "p747_local_atlas_projector_only": bool(
                p747.get("current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe")
            ),
        },
        {
            "p721_source_topology_lane_upgrades_to_exact_t176_provider": False,
            "p734_t188_exported": False,
            "p734_quotient_only": True,
            "p742_t196_exported": False,
            "p742_pair12_typed_continuation_exported": False,
            "p747_t201_exported": False,
            "p747_local_atlas_projector_only": True,
        },
        "The source-topology theorem/witness family remains quotient-only or projector-only and still does not reach the required chart-sensitive transported-section bridge level.",
    )
    add_check(
        "t216_t224_positive_attempt_family_is_closest_near_miss_but_still_not_actual_bridge_export",
        {
            "p763_exact_interface_unexported": bool(
                p763.get("current_t216_attempt_immediate_missing_interface_is_still_unexported")
            ),
            "p764_future_only_target": bool(p764.get("current_t218_target_is_future_route_only")),
            "p765_actual_interface_exported": bool(p765.get("t219_target_exported_on_current_repo_state")),
            "p766_attempt_exported": bool(p766.get("t220_attempt_exported_on_current_repo_state")),
            "p767_exact_subinterface_unexported": bool(
                p767.get("current_repo_still_does_not_export_actual_realization_of_t220_immediate_missing_subinterface")
            ),
            "p768_future_only_subinterface_target": bool(p768.get("current_t222_target_is_future_route_only")),
            "p769_actual_subinterface_exported": bool(p769.get("t223_target_exported_on_current_repo_state")),
            "p770_attempt_exported": bool(p770.get("t224_attempt_exported_on_current_repo_state")),
        },
        {
            "p763_exact_interface_unexported": True,
            "p764_future_only_target": True,
            "p765_actual_interface_exported": False,
            "p766_attempt_exported": True,
            "p767_exact_subinterface_unexported": True,
            "p768_future_only_subinterface_target": True,
            "p769_actual_subinterface_exported": False,
            "p770_attempt_exported": True,
        },
        "The positive T216/T218/T220/T222 family is the closest current near-miss, but it still remains target-only or attempt-only and not an actual bridge export.",
    )
    add_check(
        "no_current_export_lawfully_supplies_f945_bridge_output_schema",
        no_current_lawful_bridge_supplier_found,
        True,
        "Therefore no current export yet lawfully supplies the bridge_output_schema frozen by F945.",
    )

    status = (
        "PASS_T173_T176_BRIDGE_EXISTING_EXPORT_OR_NEAR_MISS_CANDIDATES_AUDITED_NO_LAWFUL_SUPPLIER_FOUND"
        if not blocking
        else "P946_REQUIRES_REVIEW_CHANGED_CANDIDATE_STATE"
    )

    candidate_families = [
        {
            "candidate_family_id": "global_rooted_transport_convention_family",
            "evidence_refs": [rel(IN_P683), rel(IN_P684), rel(IN_F684), rel(IN_F685)],
            "current_state": "real_exported_global_rooted_transport_and_closure_family",
            "schema_fit": {
                "full_C_v1_scope": True,
                "pair12_branch_sensitive": False,
                "chart_sensitive_transported_section_level": False,
                "nonconvention_nonprojective_nonpremise_smuggled": False,
            },
            "why_not_lawful_supplier": [
                "The family is rooted-transport/convention based rather than an exported chart-sensitive pair12-typed descent bridge.",
                "It explicitly remains below any strict physical orientation datum.",
                "It explicitly descends only to projective data.",
            ],
        },
        {
            "candidate_family_id": "static_t176_provider_scan_families",
            "evidence_refs": [rel(IN_P712), rel(IN_P718), rel(IN_P719), rel(IN_P720)],
            "current_state": "multiple_provider_classes_scanned_but_zero_exact_candidates",
            "schema_fit": {
                "full_C_v1_scope": False,
                "pair12_branch_sensitive": False,
                "chart_sensitive_transported_section_level": False,
                "nonconvention_nonprojective_nonpremise_smuggled": False,
            },
            "why_not_lawful_supplier": [
                "No exact candidate is exported in the currently scanned classes.",
                "The surviving exact next direction already shifts away from these static classes toward more physical flux-like providers.",
            ],
        },
        {
            "candidate_family_id": "source_topology_pair12_typed_descent_family_t188_to_t224",
            "evidence_refs": [
                rel(IN_P734),
                rel(IN_P742),
                rel(IN_P747),
                rel(IN_P763),
                rel(IN_P764),
                rel(IN_P765),
                rel(IN_P766),
                rel(IN_P767),
                rel(IN_P768),
                rel(IN_P769),
                rel(IN_P770),
            ],
            "current_state": "closest_current_near_miss_family",
            "schema_fit": {
                "full_C_v1_scope": False,
                "pair12_branch_sensitive": True,
                "chart_sensitive_transported_section_level": False,
                "nonconvention_nonprojective_nonpremise_smuggled": False,
            },
            "why_not_lawful_supplier": [
                "The upper theorem/witness layers remain quotient-only or projector-only.",
                "The exact pair12-typed interface and its deeper subinterface remain unexported.",
                "The deepest current positive objects are future-only targets or open attempts, not actual exported bridge objects on full C_v1.",
            ],
        },
    ]

    narrowest_honest_next_probe_question = (
        "whether the current T220/T222 route already exports one actual chart_label_retaining_pair12_typed_seed_or_subinterface "
        "from Sigma_sel_src_target_v1 toward the surviving F301 pair12 carrier that avoids both Q_basis_sel_v1 terminal collapse "
        "and projector_only_local_pair12_atlas collapse, and can therefore lawfully serve as the source_side_input_leg "
        "of inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_global_C_v1_strict_v1 "
        "without using rooted_w_break_transport_convention_as_physical_orientation_data"
    )

    artifact = {
        "stage": "P946",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t173_t176_bridge_output_schema_existing_export_or_near_miss_candidate_audit_only",
        "inputs": {
            "P945": rel(IN_P945),
            "F945": rel(IN_F945),
            "P712": rel(IN_P712),
            "P718": rel(IN_P718),
            "P719": rel(IN_P719),
            "P720": rel(IN_P720),
            "P721": rel(IN_P721),
            "P722": rel(IN_P722),
            "P734": rel(IN_P734),
            "P742": rel(IN_P742),
            "P747": rel(IN_P747),
            "P763": rel(IN_P763),
            "P764": rel(IN_P764),
            "P765": rel(IN_P765),
            "P766": rel(IN_P766),
            "P767": rel(IN_P767),
            "P768": rel(IN_P768),
            "P769": rel(IN_P769),
            "P770": rel(IN_P770),
            "P683": rel(IN_P683),
            "P684": rel(IN_P684),
            "F684": rel(IN_F684),
            "F685": rel(IN_F685),
        },
        "bridge_target_object_id": f945_target.get("object_id"),
        "bridge_output_schema_minimum_properties": f945_schema.get("minimum_properties"),
        "checks": checks,
        "blocking_mismatches": blocking,
        "candidate_families": candidate_families,
        "no_current_lawful_bridge_supplier_found": no_current_lawful_bridge_supplier_found,
        "closest_current_near_miss_family": {
            "family_id": "source_topology_pair12_typed_descent_family_t188_to_t224",
            "deepest_current_exported_attempt": p770.get("t224_attempt_name"),
            "deepest_current_missing_subinterface": p767.get("exact_named_missing_subinterface"),
            "why_this_is_closest": [
                "It is already pair12 typed and explicitly chart-label retaining.",
                "It explicitly avoids both Q_basis terminal collapse and projector-only atlas collapse.",
                "It is the deepest current positive route still alive below the frozen F945 bridge target.",
            ],
            "why_it_still_fails_f945_now": [
                "It does not yet export the targeted interface or subinterface itself.",
                "It remains below actual provider export and below global T176 discharge.",
                "Its deepest current objects are attempts, not an exported chart-sensitive transported flux/current-like section on full C_v1.",
            ],
        },
        "audit_conclusion": {
            "rooted_transport_convention_family_is_not_lawful_bridge_supplier": True,
            "static_t176_provider_scans_supply_no_exact_bridge_supplier": True,
            "source_topology_pair12_typed_descent_family_is_closest_near_miss_only": True,
            "narrowest_honest_next_probe_question": narrowest_honest_next_probe_question,
        },
        "hard_limits": [
            "No T176 discharge claim.",
            "No T177 discharge claim.",
            "No T185 discharge claim.",
            "No promotion of F647 to admissible S_sel_int.",
            "No promotion of rooted w_break transport into a strict physical orientation datum.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P946",
        "status": status,
        "as_of": AS_OF,
        "bridge_target_object_id": artifact["bridge_target_object_id"],
        "no_current_lawful_bridge_supplier_found": no_current_lawful_bridge_supplier_found,
        "closest_current_near_miss_family": artifact["closest_current_near_miss_family"]["family_id"],
        "deepest_current_exported_attempt": artifact["closest_current_near_miss_family"][
            "deepest_current_exported_attempt"
        ],
        "narrowest_honest_next_probe_question": narrowest_honest_next_probe_question,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
