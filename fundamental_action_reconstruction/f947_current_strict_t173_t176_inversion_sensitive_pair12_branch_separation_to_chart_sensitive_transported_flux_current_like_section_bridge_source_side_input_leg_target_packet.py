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
IN_F946 = GENERATED / "f946_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_output_schema_target_packet_summary.json"
IN_P947 = GENERATED / "p947_current_strict_t173_t176_t220_t222_chart_label_retaining_pair12_typed_seed_subinterface_source_side_input_leg_sufficiency_or_nonexport_audit_probe.json"

OUT = GENERATED / "f947_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_source_side_input_leg_target_packet.json"
OUT_SUMMARY = GENERATED / "f947_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_source_side_input_leg_target_packet_summary.json"

F945_STATUS = (
    "F945_EXECUTED_CURRENT_STRICT_T173_T176_INVERSION_SENSITIVE_PAIR12_BRANCH_SEPARATION_TO_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_BRIDGE_TARGET_PACKET_NO_FALSE_PASS"
)
F946_STATUS = (
    "F946_EXECUTED_CURRENT_STRICT_T173_T176_INVERSION_SENSITIVE_PAIR12_BRANCH_SEPARATION_TO_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_BRIDGE_OUTPUT_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
)
P947_STATUS = (
    "PASS_T173_T176_T220_T222_CHART_LABEL_RETAINING_PAIR12_TYPED_SEED_SUBINTERFACE_SOURCE_SIDE_INPUT_LEG_SUFFICIENCY_OR_NONEXPORT_AUDITED"
)
TARGET_OBJECT_ID = (
    "inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_source_side_input_leg_target_v1"
)
BRIDGE_TARGET_ID = (
    "inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_global_C_v1_strict_v1"
)
BRIDGE_OUTPUT_SCHEMA_TARGET_ID = (
    "inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_output_schema_target_v1"
)
EXACT_SUBINTERFACE_NAME = (
    "chart_label_retaining_pair12_typed_seed_from_Sigma_sel_src_target_v1_"
    "toward_the_surviving_F301_pair12_carrier_prior_to_Q_basis_sel_v1_terminal_"
    "collapse_and_prior_to_projector_only_local_pair12_atlas_collapse"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F945, IN_F946, IN_P947]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F947",
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
    f946 = load_json(IN_F946)
    p947 = load_json(IN_P947)

    audit_conclusion = p947.get("audit_conclusion") or {}

    preconditions = {
        "f945_target_packet_executed": f945.get("status") == F945_STATUS,
        "f946_output_schema_target_packet_executed": f946.get("status") == F946_STATUS,
        "p947_source_side_input_leg_probe_passed": p947.get("status") == P947_STATUS,
        "bridge_target_id_matches": ((f945.get("target_object") or {}).get("object_id") == BRIDGE_TARGET_ID),
        "bridge_output_schema_target_id_matches": f946.get("target_object_id") == BRIDGE_OUTPUT_SCHEMA_TARGET_ID,
        "no_actual_source_side_input_leg_supplier_exported": (
            audit_conclusion.get("current_repo_already_exports_actual_source_side_input_leg_supplier") is False
        ),
        "route_local_seed_subinterface_is_exactly_frozen": (
            audit_conclusion.get("exact_route_local_seed_subinterface_name") == EXACT_SUBINTERFACE_NAME
        ),
        "additional_full_c_v1_lift_still_required": (
            audit_conclusion.get(
                "route_local_seed_subinterface_even_if_actual_would_still_require_additional_full_c_v1_transported_section_lift"
            )
            is True
        ),
    }
    all_preconditions_pass = all(preconditions.values())

    status = (
        "F947_EXECUTED_CURRENT_STRICT_T173_T176_INVERSION_SENSITIVE_PAIR12_BRANCH_SEPARATION_TO_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_BRIDGE_SOURCE_SIDE_INPUT_LEG_TARGET_PACKET_NO_FALSE_PASS"
        if all_preconditions_pass
        else "F947_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "F947",
        "packet_name": "CurrentStrictT173T176InversionSensitivePair12BranchSeparationToChartSensitiveTransportedFluxCurrentLikeSectionBridgeSourceSideInputLegTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f945_bridge_target_packet": rel(IN_F945),
            "f946_bridge_output_schema_target_packet": rel(IN_F946),
            "p947_source_side_input_leg_audit": rel(IN_P947),
        },
        "preconditions": preconditions,
        "target_object": {
            "object_id": TARGET_OBJECT_ID,
            "goal": "Freeze the exact missing source_side_input_leg object for the active F945/F946 bridge family without claiming that the leg already exists or that the remaining full_C_v1 transported-section lift is solved.",
            "required_fields": [
                {
                    "name": "bridge_global_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F945 bridge target on full C_v1.",
                },
                {
                    "name": "bridge_output_schema_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F946 bridge-output-schema target.",
                },
                {
                    "name": "route_local_seed_subinterface_ref",
                    "required": True,
                    "hard_limit": "Must preserve the exact T220/T222 chart-label-retaining pair12 typed seed/subinterface name and route-local scope.",
                },
                {
                    "name": "exact_source_side_input_leg_statement",
                    "required": True,
                    "hard_limit": "Must state the exact missing source_side_input_leg object without claiming that it already exists.",
                },
                {
                    "name": "supplier_gap_ref",
                    "required": True,
                    "hard_limit": "Must keep explicit that no current actual source_side_input_leg supplier is exported.",
                },
                {
                    "name": "full_c_v1_lift_gap_ref",
                    "required": True,
                    "hard_limit": "Must keep explicit that even a future actual source-side leg would still remain below full bridge completion.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny T176 discharge, T177 discharge, T185 discharge, QW-2191 discharge, and ToE closure.",
                },
            ],
        },
        "target_refs": {
            "bridge_global_target_ref": BRIDGE_TARGET_ID,
            "bridge_output_schema_target_ref": BRIDGE_OUTPUT_SCHEMA_TARGET_ID,
            "route_local_seed_subinterface_ref": {
                "name": audit_conclusion.get("exact_route_local_seed_subinterface_name"),
                "current_best_route_local_attempt_name": audit_conclusion.get("current_best_route_local_attempt_name"),
                "current_best_future_only_subinterface_target_name": audit_conclusion.get("current_best_future_only_subinterface_target_name"),
                "current_best_open_actual_subinterface_attempt_name": audit_conclusion.get("current_best_open_actual_subinterface_attempt_name"),
                "scope": audit_conclusion.get("route_local_seed_subinterface_scope"),
            },
            "exact_source_side_input_leg_statement": {
                "object_id": TARGET_OBJECT_ID,
                "for_bridge_target": BRIDGE_TARGET_ID,
                "for_bridge_output_schema_target": BRIDGE_OUTPUT_SCHEMA_TARGET_ID,
                "typed_meaning": "one actual nonprojective nonconvention source-side input leg descending the active inversion-sensitive pair12 branch separation from Sigma_sel_src_target_v1 toward the surviving F301 pair12 carrier in a chart-label-retaining way that can feed, but not itself complete, the full-C_v1 chart-sensitive transported-section bridge",
            },
            "supplier_gap_ref": {
                "current_repo_already_exports_actual_source_side_input_leg_supplier": audit_conclusion.get(
                    "current_repo_already_exports_actual_source_side_input_leg_supplier"
                ),
                "current_best_missing_route_object": audit_conclusion.get("exact_route_local_seed_subinterface_name"),
            },
            "full_c_v1_lift_gap_ref": {
                "additional_full_c_v1_transported_section_lift_still_required": audit_conclusion.get(
                    "route_local_seed_subinterface_even_if_actual_would_still_require_additional_full_c_v1_transported_section_lift"
                ),
            },
        },
        "current_honest_reading": [
            "The current repo does not yet export the actual source_side_input_leg for the active F945/F946 bridge family.",
            "The best current positive material remains only the route-local T220/T222 seed/subinterface family around Sigma_sel_src_target_v1 toward the surviving F301 pair12 carrier.",
            "Even a future actual export of that source-side leg would still remain below full-C_v1 bridge completion and would require one additional transported-section lift.",
        ],
        "recommended_next_move": "Keep the next positive local move at actual export of the exact T222 seed/subinterface or exact failure localization below it, while keeping separate that full bridge completion would still require a later full_C_v1 transported-section lift.",
        "hard_limits": [
            "Does not claim that the exact source_side_input_leg is already exported.",
            "Does not claim that the frozen T222 seed/subinterface is already actually exported.",
            "Does not claim that the route-local seed/subinterface already suffices for full C_v1.",
            "Does not claim that the additional transported-section lift already exists.",
            "Does not claim that the exact bridge-output schema is already exported.",
            "Does not claim that T176 is discharged.",
            "Does not claim that T177 is discharged.",
            "Does not claim that T185 is discharged.",
            "Does not claim that F647 is admissible S_sel_int.",
            "Does not claim that the rooted w_break sign lift is a strict physical orientation datum.",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F947",
        "status": status,
        "as_of": AS_OF,
        "target_object_id": TARGET_OBJECT_ID,
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
