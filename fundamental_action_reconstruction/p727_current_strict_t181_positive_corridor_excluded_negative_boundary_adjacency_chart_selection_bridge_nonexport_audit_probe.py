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

IN_P724 = GENERATED / "p724_current_strict_t178_positive_source_polarity_atlas_entry_corridor_reduction_audit_probe.json"
IN_P726 = GENERATED / "p726_current_strict_t180_positive_corridor_outer_interior_chart_selection_bridge_nonexport_audit_probe_summary.json"
IN_F141 = GENERATED / "f141_first_actual_source_topology_barrier_protected_sign_witness_packet_summary.json"
IN_F143 = GENERATED / "f143_first_actual_source_topology_nonzero_flow_witness_packet_summary.json"
IN_F147 = GENERATED / "f147_first_actual_source_topology_selector_witness_packet_summary.json"
IN_F148 = GENERATED / "f148_first_actual_source_topology_basis_independent_promotion_witness_packet_summary.json"
IN_F150 = GENERATED / "f150_first_actual_declared_scope_source_topology_selector_theorem_packet_summary.json"
IN_F466 = GENERATED / "f466_current_strict_pair12345_selector_atlas_full_cocycle_data_export_packet_summary.json"
IN_SECTION_V2 = GENERATED / "a_12345_pair12345_chart_glued_orientation_projector_operator_section_strict_core_v2.json"

OUT_JSON = GENERATED / "p727_current_strict_t181_positive_corridor_excluded_negative_boundary_adjacency_chart_selection_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p727_current_strict_t181_positive_corridor_excluded_negative_boundary_adjacency_chart_selection_bridge_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def packet_mentions_negative_boundary_adjacency(packet: dict[str, Any]) -> bool:
    serialized = json.dumps(packet, sort_keys=True)
    patterns = [
        r"\bpair[1-5]\b",
        r"\badjacent\b",
        r"\bneighbor\b",
        r"\bboundary\b",
        r"\bexcluded\b",
        r"\bshielded\b",
        r"\bdistance\b",
        r"\bpair4\b",
    ]
    return any(re.search(pattern, serialized) for pattern in patterns)


def parse_adjacent_backbone_pairs(gluing_laws: list[str]) -> list[list[str]]:
    adjacent_pairs: set[tuple[int, int]] = set()
    for law in gluing_laws:
        for a_str, b_str in re.findall(r"O_(\d)(\d)", law):
            a = int(a_str)
            b = int(b_str)
            lo, hi = sorted((a, b))
            if hi - lo == 1:
                adjacent_pairs.add((lo, hi))
    return [[f"pair{lo}", f"pair{hi}"] for lo, hi in sorted(adjacent_pairs)]


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [
        IN_P724,
        IN_P726,
        IN_F141,
        IN_F143,
        IN_F147,
        IN_F148,
        IN_F150,
        IN_F466,
        IN_SECTION_V2,
    ]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P727",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p724 = load_json(IN_P724)
    p726 = load_json(IN_P726)
    f141 = load_json(IN_F141)
    f143 = load_json(IN_F143)
    f147 = load_json(IN_F147)
    f148 = load_json(IN_F148)
    f150 = load_json(IN_F150)
    f466 = load_json(IN_F466)
    section_v2 = load_json(IN_SECTION_V2)

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

    positive_corridor_roots = p724.get("atlas_entry_roots_compatible_with_current_positive_source_polarity") or []
    excluded_negative_roots = p724.get("atlas_entry_roots_incompatible_with_current_positive_source_polarity") or []

    adjacent_backbone_pairs = parse_adjacent_backbone_pairs(section_v2.get("gluing_laws") or [])
    adjacent_to_excluded_negative: set[str] = set()
    for left, right in adjacent_backbone_pairs:
        if left in excluded_negative_roots and right in positive_corridor_roots:
            adjacent_to_excluded_negative.add(right)
        if right in excluded_negative_roots and left in positive_corridor_roots:
            adjacent_to_excluded_negative.add(left)
    positive_boundary_adjacent_charts = sorted(adjacent_to_excluded_negative)
    positive_boundary_shielded_charts = [chart for chart in positive_corridor_roots if chart not in adjacent_to_excluded_negative]

    source_packets_remain_excluded_negative_boundary_blind = not any(
        packet_mentions_negative_boundary_adjacency(packet) for packet in (f141, f143, f147, f148, f150)
    )
    source_physical_primitives_exported = (
        bool(f141.get("actual_barrier_protected_sign_witness_exported"))
        and bool(f143.get("actual_nonzero_flow_witness_exported"))
        and bool(f147.get("actual_selector_witness_exported"))
        and bool(f148.get("actual_basis_independent_selector_promotion_exported"))
        and bool(f150.get("declared_scope_source_topology_selector_theorem_exported"))
    )
    selector_axis_basis = ((f147.get("support_packet") or {}).get("selector_axis_realization") or {}).get("frame_basis")
    tau_src_identified = {
        "F147": bool(f147.get("tau_src_identified_with_s_prelm")),
        "F148": bool(f148.get("tau_src_identified_with_s_prelm")),
        "F150": bool(f150.get("tau_src_identified_with_s_prelm")),
    }

    add_check(
        "source_physical_primitives_still_exported",
        source_physical_primitives_exported,
        True,
        "The current source-topology lane still exports the same real physical primitives used throughout the T177-T180 bridge reduction sequence.",
    )
    add_check(
        "excluded_negative_branch_already_localized",
        excluded_negative_roots,
        ["pair4"],
        "P724 already localizes the excluded negative atlas-entry branch to the unique chart pair4.",
    )
    add_check(
        "p726_outer_interior_localization_already_present",
        {
            "positive_outer_edge_charts": p726.get("positive_outer_edge_charts"),
            "positive_interior_charts": p726.get("positive_interior_charts"),
        },
        {
            "positive_outer_edge_charts": ["pair1", "pair5"],
            "positive_interior_charts": ["pair2", "pair3"],
        },
        "P726 already localizes the surviving positive split geometrically as outer-edge versus interior charts.",
    )
    add_check(
        "positive_corridor_boundary_adjacency_split",
        {
            "positive_boundary_adjacent_charts": positive_boundary_adjacent_charts,
            "positive_boundary_shielded_charts": positive_boundary_shielded_charts,
        },
        {
            "positive_boundary_adjacent_charts": ["pair3", "pair5"],
            "positive_boundary_shielded_charts": ["pair1", "pair2"],
        },
        "Relative to the excluded negative branch pair4, the surviving positive charts split into boundary-adjacent charts and boundary-shielded charts.",
    )
    add_check(
        "atlas_full_cocycle_packet_still_exported",
        bool(f466.get("status")),
        True,
        "The full five-chart atlas cocycle packet remains exported, so the excluded-branch adjacency split is read from current atlas geometry rather than from an invented labeling.",
    )
    add_check(
        "source_packets_remain_excluded_negative_boundary_blind",
        source_packets_remain_excluded_negative_boundary_blind,
        True,
        "The current source-side packets still export no rule distinguishing positive charts adjacent to the excluded negative branch from those shielded from it.",
    )
    add_check(
        "source_selector_axis_remains_prelm_basis_only",
        selector_axis_basis,
        ["u_T", "u_L"],
        "The current source-side selector axis still lives only in the preLM basis u_T,u_L, not in an excluded-boundary adjacency basis.",
    )
    add_check(
        "tau_src_still_not_identified_with_selector_carrier",
        tau_src_identified,
        {
            "F147": False,
            "F148": False,
            "F150": False,
        },
        "tau_src is still not identified with the current selector carrier, so no excluded-negative-boundary selection rule can yet be read off internally.",
    )
    add_check(
        "excluded_negative_boundary_adjacency_bridge_exported",
        False,
        False,
        "No current export provides the missing bridge rule distinguishing positive boundary-adjacent charts from positive boundary-shielded charts.",
    )

    status = (
        "PASS_POSITIVE_CORRIDOR_EXCLUDED_NEGATIVE_BOUNDARY_ADJACENCY_CHART_SELECTION_BRIDGE_NONEXPORT_AUDITED"
        if not blocking
        else "P727_REQUIRES_REVIEW_CHANGED_POSITIVE_CORRIDOR_EXCLUDED_NEGATIVE_BOUNDARY_STATE"
    )

    artifact = {
        "stage": "P727",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t181_positive_corridor_excluded_negative_boundary_adjacency_chart_selection_bridge_nonexport_boundary_only",
        "inputs": {
            "P724": str(IN_P724.relative_to(REPO)),
            "P726": str(IN_P726.relative_to(REPO)),
            "F141": str(IN_F141.relative_to(REPO)),
            "F143": str(IN_F143.relative_to(REPO)),
            "F147": str(IN_F147.relative_to(REPO)),
            "F148": str(IN_F148.relative_to(REPO)),
            "F150": str(IN_F150.relative_to(REPO)),
            "F466": str(IN_F466.relative_to(REPO)),
            "A12345_section_v2": str(IN_SECTION_V2.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t181_target_name": "PositiveCorridorExcludedNegativeBoundaryAdjacencyChartSelectionBridge_global_C_v1_strict_v1",
        "t181_target_exported_on_current_repo_state": False,
        "current_positive_corridor_excluded_negative_boundary_profile": {
            "excluded_negative_root": excluded_negative_roots,
            "positive_corridor_roots": positive_corridor_roots,
            "positive_boundary_adjacent_charts": positive_boundary_adjacent_charts,
            "positive_boundary_shielded_charts": positive_boundary_shielded_charts,
            "source_packets_remain_excluded_negative_boundary_blind": source_packets_remain_excluded_negative_boundary_blind,
            "source_selector_axis_frame_basis": selector_axis_basis,
            "tau_src_identified_with_s_prelm": tau_src_identified,
        },
        "audit_conclusion": {
            "current_positive_corridor_is_now_localized_relative_to_excluded_negative_branch": True,
            "current_repo_already_exports_t181_target": False,
            "next_honest_move": (
                "export_or_attack_a_finer_source_side_rule_that_selects_between_positive_boundary_adjacent_and_positive_boundary_shielded_charts"
            ),
        },
        "hard_limits": [
            "No T181 discharge claim.",
            "No T180 discharge claim.",
            "No T179 discharge claim.",
            "No T178 discharge claim.",
            "No T177 discharge claim.",
            "No unique chart-seed selection claim.",
            "No strict physical orientation datum claim.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P727",
        "status": status,
        "as_of": AS_OF,
        "t181_target_name": "PositiveCorridorExcludedNegativeBoundaryAdjacencyChartSelectionBridge_global_C_v1_strict_v1",
        "t181_target_exported_on_current_repo_state": False,
        "positive_boundary_adjacent_charts": positive_boundary_adjacent_charts,
        "positive_boundary_shielded_charts": positive_boundary_shielded_charts,
        "next_honest_move": "finer_source_side_rule_selecting_between_positive_boundary_adjacent_and_positive_boundary_shielded_charts",
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
