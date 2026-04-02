#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-21"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F945 = GENERATED / "f945_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_target_packet.json"
IN_P945 = GENERATED / "p945_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_nonexport_audit_probe_summary.json"
IN_P721 = GENERATED / "p721_current_strict_t176_source_topology_basis_free_qw2191_safe_provider_nonupgrade_audit_probe_summary.json"
IN_P722 = GENERATED / "p722_current_strict_t177_chart_sensitive_transported_flux_current_like_section_nonexport_audit_probe_summary.json"
IN_P729 = GENERATED / "p729_current_strict_t183_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"
IN_P731 = GENERATED / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe_summary.json"
IN_P742 = GENERATED / "p742_current_strict_t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P747 = GENERATED / "p747_current_strict_t201_actual_source_topology_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge_nonexport_audit_probe_summary.json"
IN_F647 = GENERATED / "f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json"
IN_P684 = GENERATED / "p684_current_strict_t173_w_break_rooted_directed_state_lift_audit_probe.json"

OUT = GENERATED / "p946_current_strict_t173_t176_output_schema_class_candidate_supported_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_output_schema_blocked.json"
OUT_SUMMARY = GENERATED / "p946_current_strict_t173_t176_output_schema_class_candidate_supported_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_output_schema_blocked_summary.json"

BRIDGE_TARGET_NAME = (
    "inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_global_C_v1_strict_v1"
)
BRIDGE_OUTPUT_SCHEMA_TARGET_NAME = (
    "inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_output_schema_target_v1"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def has_required_field(obj: dict[str, Any], name: str) -> bool:
    return any(
        isinstance(item, dict) and item.get("name") == name
        for item in (obj.get("required_fields") or [])
    )


def existing_export_mentions(target: str) -> list[str]:
    matches: list[str] = []
    for path in GENERATED.glob("*.json"):
        if path in {OUT, OUT_SUMMARY}:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            text = path.read_text(encoding="latin-1")
        if target in text:
            matches.append(rel(path))
    return sorted(matches)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [
        IN_F945,
        IN_P945,
        IN_P721,
        IN_P722,
        IN_P729,
        IN_P731,
        IN_P742,
        IN_P747,
        IN_F647,
        IN_P684,
    ]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P946",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f945 = load_json(IN_F945)
    p945 = load_json(IN_P945)
    p721 = load_json(IN_P721)
    p722 = load_json(IN_P722)
    p729 = load_json(IN_P729)
    p731 = load_json(IN_P731)
    p742 = load_json(IN_P742)
    p747 = load_json(IN_P747)
    f647 = load_json(IN_F647)
    p684 = load_json(IN_P684)

    f945_target = f945.get("target_object") or {}
    existing_output_schema_mentions = existing_export_mentions(BRIDGE_OUTPUT_SCHEMA_TARGET_NAME)

    checks = [
        {
            "id": "f945_already_names_exact_missing_bridge_output_schema_field",
            "pass": (
                f945.get("status")
                == "F945_EXECUTED_CURRENT_STRICT_T173_T176_INVERSION_SENSITIVE_PAIR12_BRANCH_SEPARATION_TO_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_BRIDGE_TARGET_PACKET_NO_FALSE_PASS"
                and f945_target.get("object_id") == BRIDGE_TARGET_NAME
                and has_required_field(f945_target, "bridge_output_schema")
            ),
            "details": "F945 already isolates bridge_output_schema as one exact missing field of the frozen T173/T176 bridge target.",
        },
        {
            "id": "p945_already_keeps_bridge_unexported_while_pair12_branch_separation_present",
            "pass": (
                p945.get("status")
                == "PASS_T173_T176_INVERSION_SENSITIVE_PAIR12_TO_CHART_SENSITIVE_PROVIDER_BRIDGE_NONEXPORT_AUDITED"
                and p945.get("bridge_target_name") == "InversionSensitivePair12BranchSeparationToChartSensitiveTransportedFluxCurrentLikeSectionBridge_global_C_v1_strict_v1"
                and p945.get("bridge_target_exported_on_current_repo_state") is False
                and p945.get("current_repo_already_exports_inversion_sensitive_pair12_branch_separation")
                is True
                and p945.get("current_repo_already_exports_inversion_sensitive_pair12_to_t176_bridge")
                is False
            ),
            "details": "P945 already keeps the frozen bridge itself unexported while preserving that inversion-sensitive pair1/pair2 branch separation already exists on the current repo state.",
        },
        {
            "id": "p721_and_p722_preserve_upstream_provider_class_and_chart_sensitive_target_class",
            "pass": (
                p721.get("status")
                == "PASS_SOURCE_TOPOLOGY_BASIS_FREE_QW2191_SAFE_PROVIDER_NONUPGRADE_AUDITED"
                and p721.get("source_topology_physically_interpretable_strict_ingredients_present")
                is True
                and p721.get("source_topology_lane_upgrades_to_exact_t176_provider") is False
                and p721.get("current_best_source_topology_output_is_basis_free_quotient_class_only")
                is True
                and p721.get("preferred_next_direction")
                == "chart_sensitive_transported_flux_or_current_like_provider_class_over_further_static_basis_free_or_output_axis_classes"
                and p722.get("status")
                == "PASS_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_NONEXPORT_AUDITED"
                and p722.get("t177_target_name")
                == "ChartSensitiveTransportedFluxCurrentLikeSection_global_C_v1_source_topology_bridge_v1"
                and p722.get("t177_target_exported_on_current_repo_state") is False
                and p722.get("current_source_topology_lane_is_physics_facing_but_chart_blind")
                is True
            ),
            "details": "P721/P722 already preserve a physics-facing source-side provider class and the exact chart-sensitive transported flux/current-like section target class, but still below exact T176 export.",
        },
        {
            "id": "p742_and_p747_preserve_neighboring_partial_bridge_classes",
            "pass": (
                p742.get("status")
                == "PASS_ACTUAL_SOURCE_TOPOLOGY_SELECTOR_WITNESS_TO_RESIDUAL_DATUM_PAIR12_TYPED_CARRIER_BRIDGE_NONEXPORT_AUDITED"
                and p742.get("t196_target_name")
                == "ActualSourceTopologySelectorWitnessToResidualDatumPair12TypedCarrierBridge_global_C_v1_strict_v1"
                and p742.get("current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation")
                is True
                and p742.get("current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation")
                is False
                and p742.get("current_actual_selector_witness_codomain_continuation_remains_basis_free_not_pair12_typed")
                is True
                and p747.get("status")
                == "PASS_ACTUAL_SOURCE_TOPOLOGY_SELECTOR_WITNESS_TARGET_TO_LOCAL_PAIR12_CHART_SENSITIVE_ATLAS_BRIDGE_NONEXPORT_AUDITED"
                and p747.get("t201_target_name")
                == "ActualSourceTopologySelectorWitnessTargetToLocalPair12ChartSensitiveAtlasBridge_sigma_int_corridor_strict_v1"
                and p747.get("current_local_pair12_chart_sensitive_atlas_lane_exported") is True
                and p747.get("current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe")
                is True
                and p747.get("current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge")
                is False
            ),
            "details": "P742/P747 already preserve neighboring partial bridge classes on the pair1/pair2 typed-carrier side and on the local chart-sensitive atlas side, but neither one supplies the frozen full-C_v1 bridge output schema.",
        },
        {
            "id": "p729_and_p731_preserve_pair12_branch_sensitivity_support",
            "pass": (
                p729.get("status")
                == "PASS_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_SELECTION_BRIDGE_NONEXPORT_AUDITED"
                and p729.get("t183_target_name")
                == "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"
                and p729.get("remaining_pair12_split_localized_as_opposite_orbit_directions")
                is True
                and p729.get("pair1_orbit_branch_kind") == "delta_k_positive_index_branch"
                and p729.get("pair2_orbit_branch_kind") == "delta_minus_k_negative_index_branch"
                and p731.get("status")
                == "PASS_W_BREAK_WITNESS_PAYLOAD_PAIR12_ORBIT_DIRECTION_PROMOTION_BRIDGE_NONEXPORT_AUDITED"
                and p731.get("t185_target_name")
                == "WBreakWitnessPayloadResidualDatumPair12OrbitDirectionPromotionBridge_global_C_v1_strict_v1"
                and p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches")
                is True
                and p731.get("pair1_w_break_branch_score_sign") == "negative"
                and p731.get("pair2_w_break_branch_score_sign") == "positive"
                and p731.get("w_break_pair12_branch_scores_are_antisymmetric") is True
            ),
            "details": "P729/P731 already preserve exact pair1/pair2 branch sensitivity together with antisymmetric witness-side separation, so branch-sensitive support is real but still below the frozen bridge output schema.",
        },
        {
            "id": "f647_and_p684_preserve_nonpromotion_limits",
            "pass": (
                f647.get("status")
                == "F647_EXECUTED_STRICT_WITNESS_PROVIDER_EXPORT_PACKET_FOR_SEED_V1_REALIZATION_ATTEMPT_NO_FALSE_PASS"
                and f647.get("constructed_source_object_exported") is True
                and f647.get("admissible_S_sel_int") is False
                and ((f647.get("strict_core_export_properties") or {}).get("kernel_split_safe")) is True
                and ((p684.get("sign_lift") or {}).get("counts_as_strict_physical_orientation_datum"))
                is False
            ),
            "details": "F647/P684 keep the current lane explicitly below admissible S_sel_int and below any strict physical orientation datum, so no hidden promotion can supply the missing bridge output schema.",
        },
        {
            "id": "no_exact_new_bridge_output_schema_export_present",
            "pass": (
                len(existing_output_schema_mentions) == 0
                and p722.get("t177_target_name") != BRIDGE_OUTPUT_SCHEMA_TARGET_NAME
                and p742.get("t196_target_name") != BRIDGE_OUTPUT_SCHEMA_TARGET_NAME
                and p747.get("t201_target_name") != BRIDGE_OUTPUT_SCHEMA_TARGET_NAME
                and p729.get("t183_target_name") != BRIDGE_OUTPUT_SCHEMA_TARGET_NAME
                and p731.get("t185_target_name") != BRIDGE_OUTPUT_SCHEMA_TARGET_NAME
            ),
            "details": "No current generated export names one exact bridge_output_schema object for the frozen inversion-sensitive pair1/pair2 to chart-sensitive transported section bridge.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "bridge_output_schema_class_candidate_supported_now": (
            checks[2]["pass"] and checks[3]["pass"] and checks[4]["pass"] and checks[5]["pass"]
        ),
        "bridge_output_schema_exported_now": False if all_pass else None,
        "sharp_blocker_field": "bridge_output_schema" if all_pass else None,
        "target_bridge_output_schema_object_id": BRIDGE_OUTPUT_SCHEMA_TARGET_NAME if all_pass else None,
        "next_honest_move_is_freeze_exact_bridge_output_schema_target": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P946_CURRENT_STRICT_T173_T176_OUTPUT_SCHEMA_CLASS_CANDIDATE_SUPPORTED_INVERSION_SENSITIVE_PAIR12_BRANCH_SEPARATION_TO_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_BRIDGE_OUTPUT_SCHEMA_BLOCKED"
        if all_pass
        else "P946_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P946",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f945_bridge_target_packet": rel(IN_F945),
            "p945_bridge_nonexport_audit_probe": rel(IN_P945),
            "p721_t176_source_topology_nonupgrade_audit": rel(IN_P721),
            "p722_t177_chart_sensitive_section_nonexport_audit": rel(IN_P722),
            "p729_pair12_orbit_direction_localization_audit": rel(IN_P729),
            "p731_w_break_pair12_branch_separation_audit": rel(IN_P731),
            "p742_pair12_typed_carrier_bridge_nonexport_audit": rel(IN_P742),
            "p747_local_pair12_chart_sensitive_atlas_bridge_nonexport_audit": rel(IN_P747),
            "f647_w_break_witness_provider_packet": rel(IN_F647),
            "p684_rooted_w_break_state_lift_audit": rel(IN_P684),
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "output_schema_class_support_stack": {
            "bridge_target_ref": f945_target.get("object_id"),
            "candidate_support_refs": [
                p722.get("t177_target_name"),
                p742.get("t196_target_name"),
                p747.get("t201_target_name"),
                p729.get("t183_target_name"),
                p731.get("t185_target_name"),
            ],
            "nonpromotion_refs": [
                f647.get("exported_constructed_source_object"),
                "P684_rooted_w_break_convention_nonpromotion",
            ],
            "support_grade": "candidate_supported_not_yet_exported",
            "scope_limit": "Current exports preserve neighboring provider, carrier, atlas, and branch-sensitive support classes only; none exports the exact bridge_output_schema required by F945 on full C_v1.",
        },
        "current_honest_reading": [
            "The repo already preserves real class support around the frozen F945 bridge: a physics-facing source-topology lane, a chart-sensitive target class, a pair1/pair2 typed-carrier-side partial continuation, a local chart-sensitive atlas lane, and exact pair1/pair2 branch sensitivity.",
            "But those supports remain split across quotient-class, local-atlas, and witness-side lanes, and none currently supplies the exact bridge_output_schema required by F945 on full C_v1.",
            "So the sharp blocker is now the still-missing exact bridge_output_schema object itself.",
        ],
        "recommended_next_packet": {
            "id": "F946_CURRENT_STRICT_T173_T176_INVERSION_SENSITIVE_PAIR12_BRANCH_SEPARATION_TO_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_BRIDGE_OUTPUT_SCHEMA_TARGET_PACKET",
            "goal": "Freeze the exact bridge_output_schema object still missing after output-schema class support is present only at candidate grade around the frozen F945 bridge target.",
            "minimum_fields": [
                "inversion_sensitive_pair12_to_chart_sensitive_bridge_target_ref",
                "output_schema_class_candidate_support_refs",
                "partial_bridge_or_local_chart_sensitive_refs",
                "pair12_branch_sensitivity_support_refs",
                "nonpromotion_refs",
                "exact_bridge_output_schema_statement",
                "hard_limits",
            ],
        },
        "existing_export_mentions_of_target_bridge_output_schema_object": existing_output_schema_mentions,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P946",
        "status": status,
        "as_of": AS_OF,
        "bridge_output_schema_class_candidate_supported_now": theorem_result[
            "bridge_output_schema_class_candidate_supported_now"
        ],
        "bridge_output_schema_exported_now": theorem_result["bridge_output_schema_exported_now"],
        "sharp_blocker_field": theorem_result["sharp_blocker_field"],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
