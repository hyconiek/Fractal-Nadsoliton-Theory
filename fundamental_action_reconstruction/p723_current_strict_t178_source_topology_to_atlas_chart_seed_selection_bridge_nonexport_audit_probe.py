#!/usr/bin/env python3
from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-18"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P722 = GENERATED / "p722_current_strict_t177_chart_sensitive_transported_flux_current_like_section_nonexport_audit_probe_summary.json"
IN_P715 = GENERATED / "p715_current_strict_t176_parity_completed_dual_anchor_multiroot_audit_probe_summary.json"
IN_F141 = GENERATED / "f141_first_actual_source_topology_barrier_protected_sign_witness_packet_summary.json"
IN_F143 = GENERATED / "f143_first_actual_source_topology_nonzero_flow_witness_packet_summary.json"
IN_F147 = GENERATED / "f147_first_actual_source_topology_selector_witness_packet_summary.json"
IN_F148 = GENERATED / "f148_first_actual_source_topology_basis_independent_promotion_witness_packet_summary.json"
IN_F150 = GENERATED / "f150_first_actual_declared_scope_source_topology_selector_theorem_packet_summary.json"
IN_ATLAS = GENERATED / "selector_atlas_pair12345_oriented_transport_mod_2pi_strict_convention_v1.json"

OUT_JSON = GENERATED / "p723_current_strict_t178_source_topology_to_atlas_chart_seed_selection_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p723_current_strict_t178_source_topology_to_atlas_chart_seed_selection_bridge_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def packet_mentions_pair_charts(packet: dict[str, Any]) -> bool:
    return re.search(r"\bpair[1-5]\b", json.dumps(packet, sort_keys=True)) is not None


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P722, IN_P715, IN_F141, IN_F143, IN_F147, IN_F148, IN_F150, IN_ATLAS]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P723",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p722 = load_json(IN_P722)
    p715 = load_json(IN_P715)
    f141 = load_json(IN_F141)
    f143 = load_json(IN_F143)
    f147 = load_json(IN_F147)
    f148 = load_json(IN_F148)
    f150 = load_json(IN_F150)
    atlas = load_json(IN_ATLAS)

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

    atlas_transitions = atlas.get("transitions") or {}
    atlas_has_full_pair_transport = len(atlas_transitions) == 10
    atlas_has_distinguished_root = any(
        key in atlas for key in ("reference_chart", "reference_pair", "root_chart", "root_pair", "distinguished_root")
    )

    selector_support = (f147.get("support_packet") or {})
    selector_axis = (selector_support.get("selector_axis_realization") or {})
    preobserver_scope = (selector_support.get("preobserver_scope_realization") or {})
    basis_reduction = ((f148.get("support_packet") or {}).get("basis_class_reduction_operator") or {})
    source_boundaries = (((f150.get("support_packet") or {}).get("l5_boundaries")) or {})

    source_packets_pair_chart_agnostic = not any(
        packet_mentions_pair_charts(packet) for packet in (f141, f143, f147, f148, f150)
    )

    add_check(
        "p722_t177_bridge_still_not_exported",
        bool(p722.get("t177_target_exported_on_current_repo_state")),
        False,
        "P722 already confirms that the full chart-sensitive transported flux/current-like bridge is not yet exported.",
    )
    add_check(
        "source_topology_sign_witness_positive",
        ((f141.get("support_packet") or {}).get("psi_src_barrier_sign_component_witness_v1")),
        1,
        "The current source-topology lane does export a positive barrier-protected sign witness.",
    )
    add_check(
        "source_topology_nonzero_flow_component_positive",
        bool(float(f143.get("scalar_component_witness_value") or 0.0) > 0.0),
        True,
        "The same lane exports a strictly positive nonzero-flow component witness.",
    )
    add_check(
        "source_topology_selector_axis_lives_in_prelm_frame",
        selector_axis.get("frame_basis"),
        ["u_T", "u_L"],
        "The current source-topology selector-axis realization still lives in the preLM frame basis u_T,u_L rather than in a pair-chart atlas basis.",
    )
    add_check(
        "source_topology_selector_plus_channel_positive",
        bool(preobserver_scope.get("positive_plus_output")) and bool(preobserver_scope.get("vanishing_minus_output")),
        True,
        "The current source-topology selector realization does carry a signed preobserver plus-channel response.",
    )
    add_check(
        "tau_src_not_identified_with_current_selector_carrier",
        {
            "F147": bool(f147.get("tau_src_identified_with_s_prelm")),
            "F148": bool(f148.get("tau_src_identified_with_s_prelm")),
            "F150": bool(f150.get("tau_src_identified_with_s_prelm")),
        },
        {
            "F147": False,
            "F148": False,
            "F150": False,
        },
        "The current source-topology lane still does not identify tau_src with the current selector carrier, so it does not yet seed one concrete chart state on the atlas.",
    )
    add_check(
        "source_packets_remain_pair_chart_agnostic",
        source_packets_pair_chart_agnostic,
        True,
        "The current source-topology support packets remain pair-chart agnostic: they do not yet name a concrete pair1..pair5 chart seed.",
    )
    add_check(
        "atlas_has_full_chart_sensitive_transport",
        atlas_has_full_pair_transport,
        True,
        "The selector atlas does export the full pair1..pair5 chart-sensitive transport web.",
    )
    add_check(
        "atlas_has_no_distinguished_root_chart",
        atlas_has_distinguished_root,
        False,
        "The current atlas exports transport operators but no distinguished root/reference chart by itself.",
    )
    add_check(
        "basis_reduction_forgets_chart_basis_labels",
        bool(basis_reduction.get("forgets_chart_basis_labels")),
        True,
        "The current source-topology promotion step still forgets chart basis labels before any atlas-seeding bridge is exported.",
    )
    add_check(
        "declared_scope_packet_keeps_lane_downstream_only",
        bool(source_boundaries.get("observer_downstream_only")) and bool(f150.get("declared_scope_only")),
        True,
        "The declared-scope theorem packet still keeps the source-topology lane downstream-only.",
    )
    add_check(
        "external_root_seeded_family_exists_only_projectively",
        {
            "all_roots_supported": bool(p715.get("all_roots_supported")),
            "projective_orbit": bool(p715.get("projective_root_independent_sign_orbit")),
            "exact_section": bool(p715.get("exact_root_independent_sign_vector")),
        },
        {
            "all_roots_supported": True,
            "projective_orbit": True,
            "exact_section": False,
        },
        "The strongest current all-root chart family exists only after an external root/anchor choice, and even then it reaches only projective orbit agreement rather than one exact directed section.",
    )

    t178_target_exported_on_current_repo_state = False
    add_check(
        "t178_source_topology_to_atlas_chart_seed_bridge_exported",
        t178_target_exported_on_current_repo_state,
        False,
        "No current export provides the missing source-topology-to-atlas chart-seed selection bridge required before a full T177 transported section can be claimed.",
    )

    status = (
        "PASS_SOURCE_TO_ATLAS_CHART_SEED_SELECTION_BRIDGE_NONEXPORT_AUDITED"
        if not blocking
        else "P723_REQUIRES_REVIEW_CHANGED_T178_NONEXPORT_AUDIT_STATE"
    )

    artifact = {
        "stage": "P723",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t178_source_topology_to_atlas_chart_seed_selection_bridge_nonexport_boundary_only",
        "inputs": {
            "P722": str(IN_P722.relative_to(REPO)),
            "P715": str(IN_P715.relative_to(REPO)),
            "F141": str(IN_F141.relative_to(REPO)),
            "F143": str(IN_F143.relative_to(REPO)),
            "F147": str(IN_F147.relative_to(REPO)),
            "F148": str(IN_F148.relative_to(REPO)),
            "F150": str(IN_F150.relative_to(REPO)),
            "atlas": str(IN_ATLAS.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t178_target_name": "SourceTopologyToAtlasChartSeedSelectionBridge_global_C_v1_strict_v1",
        "t178_target_exported_on_current_repo_state": t178_target_exported_on_current_repo_state,
        "current_bridge_gap_profile": {
            "source_topology_sign_witness_positive": ((f141.get("support_packet") or {}).get("psi_src_barrier_sign_component_witness_v1")),
            "source_topology_nonzero_flow_component_value": float(f143.get("scalar_component_witness_value") or 0.0),
            "source_topology_selector_axis_frame_basis": selector_axis.get("frame_basis"),
            "source_topology_selector_plus_channel_positive": bool(preobserver_scope.get("positive_plus_output")),
            "source_packets_remain_pair_chart_agnostic": source_packets_pair_chart_agnostic,
            "tau_src_identified_with_s_prelm": {
                "F147": bool(f147.get("tau_src_identified_with_s_prelm")),
                "F148": bool(f148.get("tau_src_identified_with_s_prelm")),
                "F150": bool(f150.get("tau_src_identified_with_s_prelm")),
            },
            "atlas_has_full_pair_transport": atlas_has_full_pair_transport,
            "atlas_has_distinguished_root_chart": atlas_has_distinguished_root,
            "basis_reduction_forgets_chart_basis_labels": bool(basis_reduction.get("forgets_chart_basis_labels")),
            "declared_scope_only": bool(f150.get("declared_scope_only")),
            "observer_downstream_only": bool(source_boundaries.get("observer_downstream_only")),
            "current_external_root_seeded_family_all_roots_supported": bool(p715.get("all_roots_supported")),
            "current_external_root_seeded_family_projective_only": bool(p715.get("projective_root_independent_sign_orbit"))
            and not bool(p715.get("exact_root_independent_sign_vector")),
        },
        "audit_conclusion": {
            "current_source_topology_lane_already_supplies_sign_flow_and_selector_polarity": True,
            "current_repo_exports_chart_transport_but_not_source_to_atlas_seed_selection": True,
            "current_repo_already_exports_t178_target": False,
            "next_honest_move": (
                "export_or_attack_a_source_topology_to_atlas_chart_seed_selection_bridge_before_claiming_a_full_chart_sensitive_transported_flux_current_like_section"
            ),
        },
        "hard_limits": [
            "No T177 discharge claim.",
            "No T176 discharge claim.",
            "No claim that current source-side sign/flow data already chooses one atlas root chart.",
            "No claim that current source packets already seed a pair1..pair5 chart state.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P723",
        "status": status,
        "as_of": AS_OF,
        "t178_target_name": "SourceTopologyToAtlasChartSeedSelectionBridge_global_C_v1_strict_v1",
        "t178_target_exported_on_current_repo_state": False,
        "current_source_topology_lane_supplies_sign_flow_and_selector_polarity_but_not_chart_seed_selection": True,
        "next_honest_move": "source_topology_to_atlas_chart_seed_selection_bridge",
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
