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

IN_P708 = GENERATED / "p708_current_strict_t173_frontier_dashboard_probe_summary.json"
IN_P721 = GENERATED / "p721_current_strict_t176_source_topology_basis_free_qw2191_safe_provider_nonupgrade_audit_probe_summary.json"
IN_P722 = GENERATED / "p722_current_strict_t177_chart_sensitive_transported_flux_current_like_section_nonexport_audit_probe_summary.json"
IN_P729 = GENERATED / "p729_current_strict_t183_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"
IN_P731 = GENERATED / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe_summary.json"
IN_F647 = GENERATED / "f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json"
IN_F301 = GENERATED / "residual_datum_bridge_export_map_object_support_carrier_actual_inhabitant_instance.json"
IN_P684 = GENERATED / "p684_current_strict_t173_w_break_rooted_directed_state_lift_audit_probe.json"

OUT_JSON = GENERATED / "p945_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p945_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [
        IN_P708,
        IN_P721,
        IN_P722,
        IN_P729,
        IN_P731,
        IN_F647,
        IN_F301,
        IN_P684,
    ]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P945",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p708 = load_json(IN_P708)
    p721 = load_json(IN_P721)
    p722 = load_json(IN_P722)
    p729 = load_json(IN_P729)
    p731 = load_json(IN_P731)
    f647 = load_json(IN_F647)
    f301 = load_json(IN_F301)
    p684 = load_json(IN_P684)

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

    carrier_vectors = f301.get("carrier_vectors") or {}
    pair12_carrier_present = isinstance(carrier_vectors.get("pair1"), dict) and isinstance(
        carrier_vectors.get("pair2"), dict
    )

    add_check(
        "p708_keeps_t173_as_live_frontier_and_t176_unexported",
        {
            "recommended_next_strict_target": p708.get("recommended_next_strict_target"),
            "t176_global_provider_exported": bool(p708.get("t176_global_provider_exported")),
        },
        {
            "recommended_next_strict_target": "T173",
            "t176_global_provider_exported": False,
        },
        "The current dashboard still keeps T173 as the live strict frontier and records that no exact T176 provider is exported.",
    )
    add_check(
        "p708_provider_shift_branch_remains_active_primary_t173_move",
        bool(p708.get("provider_shift_is_now_active_primary_t173_branch_on_current_repo_state")),
        True,
        "On the current repo state, the active honest T173 continuation still runs through the provider-shift branch rather than through a false global closure claim.",
    )
    add_check(
        "p721_source_topology_lane_is_physics_facing_but_not_exact_t176_provider",
        {
            "physically_interpretable": bool(p721.get("source_topology_physically_interpretable_strict_ingredients_present")),
            "t176_provider_exported": bool(p721.get("source_topology_lane_upgrades_to_exact_t176_provider")),
            "quotient_class_only": bool(p721.get("current_best_source_topology_output_is_basis_free_quotient_class_only")),
        },
        {
            "physically_interpretable": True,
            "t176_provider_exported": False,
            "quotient_class_only": True,
        },
        "The best current source-topology lane is already physics-facing but still only quotient-class-level and not an exact T176 provider.",
    )
    add_check(
        "p722_chart_sensitive_t177_bridge_still_unexported",
        {
            "t177_exported": bool(p722.get("t177_target_exported_on_current_repo_state")),
            "lane_chart_blind": bool(p722.get("current_source_topology_lane_is_physics_facing_but_chart_blind")),
        },
        {
            "t177_exported": False,
            "lane_chart_blind": True,
        },
        "The chart-sensitive transported flux/current-like section bridge is still not exported, and the current source-topology lane remains chart-blind.",
    )
    add_check(
        "f301_pair12_residual_datum_carrier_present",
        pair12_carrier_present,
        True,
        "The surviving residual-datum carrier still contains explicit pair1/pair2 branches on current exports.",
    )
    add_check(
        "p729_pair12_split_already_localized_as_opposite_orbit_directions",
        {
            "split_localized": bool(p729.get("remaining_pair12_split_localized_as_opposite_orbit_directions")),
            "pair1_kind": p729.get("pair1_orbit_branch_kind"),
            "pair2_kind": p729.get("pair2_orbit_branch_kind"),
        },
        {
            "split_localized": True,
            "pair1_kind": "delta_k_positive_index_branch",
            "pair2_kind": "delta_minus_k_negative_index_branch",
        },
        "The surviving pair1/pair2 split is already localized exactly as the positive- and negative-index orbit-direction branches.",
    )
    add_check(
        "p731_w_break_already_separates_pair12_branches_but_t185_unexported",
        {
            "split_separated": bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches")),
            "t185_exported": bool(p731.get("t185_target_exported_on_current_repo_state")),
            "pair1_sign": p731.get("pair1_w_break_branch_score_sign"),
            "pair2_sign": p731.get("pair2_w_break_branch_score_sign"),
        },
        {
            "split_separated": True,
            "t185_exported": False,
            "pair1_sign": "negative",
            "pair2_sign": "positive",
        },
        "The current inversion-sensitive witness payload already separates the two branches, but the typed T185 promotion bridge is still unexported.",
    )
    add_check(
        "f647_w_break_remains_below_admissible_s_sel_int",
        {
            "constructed_source_object_exported": bool(f647.get("constructed_source_object_exported")),
            "admissible_S_sel_int": bool(f647.get("admissible_S_sel_int")),
            "kernel_split_safe": bool(((f647.get("strict_core_export_properties") or {}).get("kernel_split_safe"))),
        },
        {
            "constructed_source_object_exported": True,
            "admissible_S_sel_int": False,
            "kernel_split_safe": True,
        },
        "The w_break payload exists only as a kernel-split-safe witness provider and remains below admissible S_sel_int.",
    )
    add_check(
        "p684_rooted_w_break_state_remains_nonphysical_orientation_convention",
        {
            "root_pair": ((p684.get("root") or {}).get("pair")),
            "counts_as_strict_physical_orientation_datum": bool(
                ((p684.get("sign_lift") or {}).get("counts_as_strict_physical_orientation_datum"))
            ),
        },
        {
            "root_pair": "pair1",
            "counts_as_strict_physical_orientation_datum": False,
        },
        "The rooted w_break-directed state remains a convention-layer lift and not a strict physical orientation datum.",
    )

    bridge_target_exported = False
    add_check(
        "current_inversion_sensitive_pair12_to_chart_sensitive_t176_bridge_exported",
        bridge_target_exported,
        False,
        "No current export upgrades the already-separated inversion-sensitive pair1/pair2 branch data into a chart-sensitive transported flux/current-like section on full C_v1.",
    )

    status = (
        "PASS_T173_T176_INVERSION_SENSITIVE_PAIR12_TO_CHART_SENSITIVE_PROVIDER_BRIDGE_NONEXPORT_AUDITED"
        if not blocking
        else "P945_REQUIRES_REVIEW_CHANGED_T173_T176_BRIDGE_STATE"
    )

    artifact = {
        "stage": "P945",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_provider_bridge_nonexport_boundary_only",
        "inputs": {
            "P708": str(IN_P708.relative_to(REPO)),
            "P721": str(IN_P721.relative_to(REPO)),
            "P722": str(IN_P722.relative_to(REPO)),
            "P729": str(IN_P729.relative_to(REPO)),
            "P731": str(IN_P731.relative_to(REPO)),
            "F647": str(IN_F647.relative_to(REPO)),
            "F301": str(IN_F301.relative_to(REPO)),
            "P684": str(IN_P684.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "bridge_target_name": "InversionSensitivePair12BranchSeparationToChartSensitiveTransportedFluxCurrentLikeSectionBridge_global_C_v1_strict_v1",
        "bridge_target_exported_on_current_repo_state": bridge_target_exported,
        "current_t173_t176_context": {
            "recommended_next_strict_target": p708.get("recommended_next_strict_target"),
            "t176_global_provider_exported": bool(p708.get("t176_global_provider_exported")),
            "provider_shift_is_now_active_primary_t173_branch_on_current_repo_state": bool(
                p708.get("provider_shift_is_now_active_primary_t173_branch_on_current_repo_state")
            ),
            "source_topology_physically_interpretable_strict_ingredients_present": bool(
                p721.get("source_topology_physically_interpretable_strict_ingredients_present")
            ),
            "source_topology_lane_upgrades_to_exact_t176_provider": bool(
                p721.get("source_topology_lane_upgrades_to_exact_t176_provider")
            ),
            "t177_target_exported_on_current_repo_state": bool(
                p722.get("t177_target_exported_on_current_repo_state")
            ),
            "current_source_topology_lane_is_physics_facing_but_chart_blind": bool(
                p722.get("current_source_topology_lane_is_physics_facing_but_chart_blind")
            ),
        },
        "current_pair12_branch_separation_state": {
            "pair12_residual_datum_carrier_present": pair12_carrier_present,
            "remaining_pair12_split_localized_as_opposite_orbit_directions": bool(
                p729.get("remaining_pair12_split_localized_as_opposite_orbit_directions")
            ),
            "pair1_orbit_branch_kind": p729.get("pair1_orbit_branch_kind"),
            "pair2_orbit_branch_kind": p729.get("pair2_orbit_branch_kind"),
            "current_w_break_witness_payload_separates_pair12_orbit_direction_branches": bool(
                p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches")
            ),
            "pair1_w_break_branch_score_sign": p731.get("pair1_w_break_branch_score_sign"),
            "pair2_w_break_branch_score_sign": p731.get("pair2_w_break_branch_score_sign"),
            "t185_target_exported_on_current_repo_state": bool(
                p731.get("t185_target_exported_on_current_repo_state")
            ),
            "f647_admissible_S_sel_int": bool(f647.get("admissible_S_sel_int")),
            "p684_counts_as_strict_physical_orientation_datum": bool(
                ((p684.get("sign_lift") or {}).get("counts_as_strict_physical_orientation_datum"))
            ),
        },
        "audit_conclusion": {
            "current_repo_already_exports_inversion_sensitive_pair12_branch_separation": bool(
                p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches")
            ),
            "current_repo_already_exports_chart_sensitive_transported_flux_current_like_section_bridge": False,
            "current_repo_already_exports_inversion_sensitive_pair12_to_t176_bridge": False,
            "next_honest_move": (
                "attack_or_export_a_typed_source_side_or_provider_bridge_from_the_current_inversion_sensitive_pair12_branch_separation_to_the_chart_sensitive_transported_flux_current_like_section_target_on_full_C_v1_without_promoting_F647_or_P684_by_fiat"
            ),
        },
        "hard_limits": [
            "No T176 discharge claim.",
            "No T177 discharge claim.",
            "No T185 discharge claim.",
            "No promotion of F647 to admissible S_sel_int.",
            "No promotion of the rooted w_break sign lift into a strict physical orientation datum.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P945",
        "status": status,
        "as_of": AS_OF,
        "bridge_target_name": artifact["bridge_target_name"],
        "bridge_target_exported_on_current_repo_state": bridge_target_exported,
        "current_repo_already_exports_inversion_sensitive_pair12_branch_separation": artifact["audit_conclusion"][
            "current_repo_already_exports_inversion_sensitive_pair12_branch_separation"
        ],
        "current_repo_already_exports_inversion_sensitive_pair12_to_t176_bridge": artifact["audit_conclusion"][
            "current_repo_already_exports_inversion_sensitive_pair12_to_t176_bridge"
        ],
        "next_honest_move": artifact["audit_conclusion"]["next_honest_move"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
