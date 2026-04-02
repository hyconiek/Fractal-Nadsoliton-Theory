#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-21"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P945 = GENERATED / "p945_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_nonexport_audit_probe.json"
IN_P708 = GENERATED / "p708_current_strict_t173_frontier_dashboard_probe_summary.json"
IN_P721 = GENERATED / "p721_current_strict_t176_source_topology_basis_free_qw2191_safe_provider_nonupgrade_audit_probe_summary.json"
IN_P722 = GENERATED / "p722_current_strict_t177_chart_sensitive_transported_flux_current_like_section_nonexport_audit_probe_summary.json"
IN_P729 = GENERATED / "p729_current_strict_t183_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"
IN_P731 = GENERATED / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe_summary.json"
IN_F647 = GENERATED / "f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json"
IN_P684 = GENERATED / "p684_current_strict_t173_w_break_rooted_directed_state_lift_audit_probe.json"

OUT = GENERATED / "f945_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_target_packet.json"
OUT_SUMMARY = GENERATED / "f945_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P945, IN_P708, IN_P721, IN_P722, IN_P729, IN_P731, IN_F647, IN_P684]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F945",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p945 = load_json(IN_P945)
    p708 = load_json(IN_P708)
    p721 = load_json(IN_P721)
    p722 = load_json(IN_P722)
    p729 = load_json(IN_P729)
    p731 = load_json(IN_P731)
    f647 = load_json(IN_F647)
    p684 = load_json(IN_P684)

    p945_conclusion = p945.get("audit_conclusion") or {}

    preconditions = {
        "p945_status_pass": (
            p945.get("status")
            == "PASS_T173_T176_INVERSION_SENSITIVE_PAIR12_TO_CHART_SENSITIVE_PROVIDER_BRIDGE_NONEXPORT_AUDITED"
        ),
        "p945_bridge_target_unexported": p945.get("bridge_target_exported_on_current_repo_state") is False,
        "p945_pair12_branch_separation_present": (
            p945_conclusion.get("current_repo_already_exports_inversion_sensitive_pair12_branch_separation") is True
        ),
        "p945_pair12_to_t176_bridge_still_unexported": (
            p945_conclusion.get("current_repo_already_exports_inversion_sensitive_pair12_to_t176_bridge") is False
        ),
        "p708_recommended_target_is_t173": p708.get("recommended_next_strict_target") == "T173",
        "p708_t176_provider_unexported": p708.get("t176_global_provider_exported") is False,
        "p721_t176_nonupgrade_still_holds": p721.get("source_topology_lane_upgrades_to_exact_t176_provider") is False,
        "p722_t177_target_unexported": p722.get("t177_target_exported_on_current_repo_state") is False,
        "p729_pair12_split_localized": p729.get("remaining_pair12_split_localized_as_opposite_orbit_directions") is True,
        "p731_pair12_branch_separation_present": (
            p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches") is True
        ),
        "p731_t185_unexported": p731.get("t185_target_exported_on_current_repo_state") is False,
        "f647_not_admissible_s_sel_int": f647.get("admissible_S_sel_int") is False,
        "p684_rooted_w_break_still_nonphysical": (
            ((p684.get("sign_lift") or {}).get("counts_as_strict_physical_orientation_datum")) is False
        ),
    }
    all_preconditions_pass = all(preconditions.values())

    if all_preconditions_pass:
        status = "F945_EXECUTED_CURRENT_STRICT_T173_T176_INVERSION_SENSITIVE_PAIR12_BRANCH_SEPARATION_TO_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_BRIDGE_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F945_REQUIRES_REVIEW"

    artifact = {
        "stage": "F945",
        "packet_name": "CurrentStrictT173T176InversionSensitivePair12BranchSeparationToChartSensitiveTransportedFluxCurrentLikeSectionBridgeTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p945_bridge_nonexport_audit": rel(IN_P945),
            "p708_t173_frontier_dashboard": rel(IN_P708),
            "p721_t176_source_topology_nonupgrade_audit": rel(IN_P721),
            "p722_t177_chart_sensitive_section_nonexport_audit": rel(IN_P722),
            "p729_pair12_orbit_direction_localization_audit": rel(IN_P729),
            "p731_w_break_pair12_branch_separation_audit": rel(IN_P731),
            "f647_w_break_witness_provider_packet": rel(IN_F647),
            "p684_rooted_w_break_state_lift_audit": rel(IN_P684),
        },
        "why_this_packet_exists": [
            "P945 shows that the current repo already exports an inversion-sensitive pair1/pair2 branch separation while the T176/T177 bridge remains unexported.",
            "P721/P722 keep the source-topology lane physics-facing but chart-blind and keep the chart-sensitive transported flux/current-like section target unexported.",
            "P731/F647/P684 keep the current w_break lane explicitly below admissible S_sel_int and below any strict physical orientation datum.",
        ],
        "preconditions": preconditions,
        "target_object": {
            "object_id": "inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_global_C_v1_strict_v1",
            "goal": "Freeze the exact missing bridge object from the current inversion-sensitive pair1/pair2 branch separation to a chart-sensitive transported flux/current-like section on full C_v1 without promoting witness-side data by fiat.",
            "required_fields": [
                {
                    "name": "t173_t176_frontier_context_ref",
                    "required": True,
                    "hard_limit": "Must keep explicit that T173 remains live and that no exact T176 provider is currently exported."
                },
                {
                    "name": "source_topology_physics_facing_but_chart_blind_ref",
                    "required": True,
                    "hard_limit": "Must keep explicit that the current source-topology lane is real and physics-facing but still chart-blind."
                },
                {
                    "name": "localized_pair12_orbit_direction_split_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact surviving pair1/pair2 split delta_k versus delta_-k and not replace it with a weaker chart-label statement."
                },
                {
                    "name": "inversion_sensitive_pair12_branch_separation_ref",
                    "required": True,
                    "hard_limit": "Must preserve that the current pair12 branch separation is witness-side only and not silently promote it into a typed source-side selector."
                },
                {
                    "name": "chart_sensitive_transported_flux_current_like_section_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact T177-class chart-sensitive transported flux/current-like section target and not silently weaken it to projective or convention-layer data."
                },
                {
                    "name": "w_break_nonpromotion_ref",
                    "required": True,
                    "hard_limit": "Must explicitly deny promotion of F647 into admissible S_sel_int."
                },
                {
                    "name": "rooted_w_break_convention_nonpromotion_ref",
                    "required": True,
                    "hard_limit": "Must explicitly deny promotion of the rooted w_break sign lift into a strict physical orientation datum."
                },
                {
                    "name": "bridge_output_schema",
                    "required": True,
                    "hard_limit": "Must state what a successful typed bridge would output on full C_v1 without claiming that such a bridge already exists."
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny T176 discharge, T177 discharge, T185 discharge, kernel-alone/global QW-2191 discharge, and ToE closure."
                }
            ]
        },
        "target_refs": {
            "t173_t176_frontier_context_ref": {
                "recommended_next_strict_target": p708.get("recommended_next_strict_target"),
                "t176_global_provider_exported": p708.get("t176_global_provider_exported"),
            },
            "source_topology_physics_facing_but_chart_blind_ref": {
                "p721_summary": "P721_CURRENT_STRICT_T176_SOURCE_TOPOLOGY_BASIS_FREE_QW2191_SAFE_PROVIDER_NONUPGRADE_AUDIT_PROBE",
                "p722_summary": "P722_CURRENT_STRICT_T177_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_NONEXPORT_AUDIT_PROBE",
                "preferred_next_direction": p721.get("preferred_next_direction"),
            },
            "localized_pair12_orbit_direction_split_ref": {
                "p729_target_name": p729.get("t183_target_name"),
                "pair1_orbit_branch_kind": p729.get("pair1_orbit_branch_kind"),
                "pair2_orbit_branch_kind": p729.get("pair2_orbit_branch_kind"),
            },
            "inversion_sensitive_pair12_branch_separation_ref": {
                "p731_target_name": p731.get("t185_target_name"),
                "pair1_w_break_branch_score_sign": p731.get("pair1_w_break_branch_score_sign"),
                "pair2_w_break_branch_score_sign": p731.get("pair2_w_break_branch_score_sign"),
            },
            "chart_sensitive_transported_flux_current_like_section_target_ref": p722.get("t177_target_name"),
            "w_break_nonpromotion_ref": {
                "exported_constructed_source_object": f647.get("exported_constructed_source_object"),
                "admissible_S_sel_int": f647.get("admissible_S_sel_int"),
            },
            "rooted_w_break_convention_nonpromotion_ref": {
                "root_pair": (p684.get("root") or {}).get("pair"),
                "counts_as_strict_physical_orientation_datum": (p684.get("sign_lift") or {}).get(
                    "counts_as_strict_physical_orientation_datum"
                ),
            },
        },
        "current_honest_reading": [
            "The current repo already exports inversion-sensitive pair1/pair2 branch separation on the surviving residual-datum carrier.",
            "The current source-topology lane remains genuinely physics-facing but still chart-blind, and the chart-sensitive transported flux/current-like section target remains unexported.",
            "So the sharp blocker is now the exact missing bridge from current pair12 branch separation to the T177-class chart-sensitive section target on full C_v1.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current export can lawfully supply bridge_output_schema for inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_global_C_v1_strict_v1 without promoting F647 or P684 by fiat.",
        "bridge_output_schema": {
            "successful_bridge_output": "one typed strict source-side or provider bridge object on full C_v1 descending the current inversion-sensitive pair12 branch separation to a chart-sensitive transported flux/current-like section",
            "minimum_properties": [
                "full_C_v1_scope",
                "pair12_branch_sensitive",
                "chart_sensitive_transported_section_level",
                "nonconvention_nonprojective_nonpremise_smuggled",
            ],
        },
        "hard_limits": [
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
        "stage": "F945",
        "status": status,
        "as_of": AS_OF,
        "target_object_id": artifact["target_object"]["object_id"],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
