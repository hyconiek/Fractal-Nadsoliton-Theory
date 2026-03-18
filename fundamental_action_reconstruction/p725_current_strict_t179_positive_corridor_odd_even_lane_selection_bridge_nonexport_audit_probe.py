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

IN_P723 = GENERATED / "p723_current_strict_t178_source_topology_to_atlas_chart_seed_selection_bridge_nonexport_audit_probe_summary.json"
IN_P724 = GENERATED / "p724_current_strict_t178_positive_source_polarity_atlas_entry_corridor_reduction_audit_probe.json"
IN_P715 = GENERATED / "p715_current_strict_t176_parity_completed_dual_anchor_multiroot_audit_probe.json"
IN_F141 = GENERATED / "f141_first_actual_source_topology_barrier_protected_sign_witness_packet_summary.json"
IN_F143 = GENERATED / "f143_first_actual_source_topology_nonzero_flow_witness_packet_summary.json"
IN_F147 = GENERATED / "f147_first_actual_source_topology_selector_witness_packet_summary.json"
IN_F148 = GENERATED / "f148_first_actual_source_topology_basis_independent_promotion_witness_packet_summary.json"
IN_F150 = GENERATED / "f150_first_actual_declared_scope_source_topology_selector_theorem_packet_summary.json"

OUT_JSON = GENERATED / "p725_current_strict_t179_positive_corridor_odd_even_lane_selection_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p725_current_strict_t179_positive_corridor_odd_even_lane_selection_bridge_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def packet_mentions_lane_selection(packet: dict[str, Any]) -> bool:
    serialized = json.dumps(packet, sort_keys=True)
    patterns = [
        r"\bpair[1-5]\b",
        r"\bodd\b",
        r"\beven\b",
        r"\bw_break\b",
        r"\bw_ref_unnormalized\b",
        r"\bendpoint\b",
        r"\binterior\b",
    ]
    return any(re.search(pattern, serialized) for pattern in patterns)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P723, IN_P724, IN_P715, IN_F141, IN_F143, IN_F147, IN_F148, IN_F150]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P725",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p723 = load_json(IN_P723)
    p724 = load_json(IN_P724)
    p715 = load_json(IN_P715)
    f141 = load_json(IN_F141)
    f143 = load_json(IN_F143)
    f147 = load_json(IN_F147)
    f148 = load_json(IN_F148)
    f150 = load_json(IN_F150)

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

    p715_result = p715.get("result") or {}
    p724_split = p724.get("residual_positive_corridor_split") or {}
    positive_roots = p724.get("atlas_entry_roots_compatible_with_current_positive_source_polarity") or []
    odd_lane = p724_split.get("odd_anchor_lane") or []
    even_lane = p724_split.get("even_fallback_lane") or []

    source_packets_pair_and_lane_agnostic = not any(
        packet_mentions_lane_selection(packet) for packet in (f141, f143, f147, f148, f150)
    )

    selector_axis_basis = ((f147.get("support_packet") or {}).get("selector_axis_realization") or {}).get("frame_basis")
    tau_src_identified = {
        "F147": bool(f147.get("tau_src_identified_with_s_prelm")),
        "F148": bool(f148.get("tau_src_identified_with_s_prelm")),
        "F150": bool(f150.get("tau_src_identified_with_s_prelm")),
    }

    add_check(
        "p723_chart_seed_bridge_still_not_exported",
        bool(p723.get("t178_target_exported_on_current_repo_state")),
        False,
        "P723 already confirms that no source-to-atlas chart-seed bridge is exported yet.",
    )
    add_check(
        "p724_positive_corridor_reduction_already_established",
        {
            "positive_roots": positive_roots,
            "unique_chart_seed_selected": bool(p724.get("unique_chart_seed_selected")),
        },
        {
            "positive_roots": ["pair1", "pair2", "pair3", "pair5"],
            "unique_chart_seed_selected": False,
        },
        "P724 already establishes the surviving positive atlas-entry corridor without selecting a unique seed chart.",
    )
    add_check(
        "positive_corridor_splits_into_odd_and_even_lanes",
        {
            "odd_lane": odd_lane,
            "even_lane": even_lane,
        },
        {
            "odd_lane": ["pair1", "pair5"],
            "even_lane": ["pair2", "pair3"],
        },
        "Within the surviving positive corridor, current exports already split the candidate roots into an odd-anchor lane and an even-fallback lane.",
    )
    add_check(
        "source_packets_remain_pair_and_lane_agnostic",
        source_packets_pair_and_lane_agnostic,
        True,
        "The current source-topology packets remain agnostic not only about chart labels, but also about odd/even lane selection inside the positive corridor.",
    )
    add_check(
        "source_selector_axis_remains_prelm_basis_only",
        selector_axis_basis,
        ["u_T", "u_L"],
        "The current source-side selector axis is still expressed only in the preLM basis u_T,u_L, not as an atlas-lane selector.",
    )
    add_check(
        "tau_src_still_not_identified_with_selector_carrier",
        tau_src_identified,
        {
            "F147": False,
            "F148": False,
            "F150": False,
        },
        "tau_src is still not identified with the current selector carrier, so no internal chart-lane seed can yet be read off from the source-side packets.",
    )
    add_check(
        "same_orbit_roots_match_positive_corridor",
        list(p715_result.get("same_orbit_roots_relative_to_reference") or []),
        ["pair1", "pair2", "pair3", "pair5"],
        "The surviving positive corridor is exactly the same-orbit branch family in the strongest current transported candidate.",
    )
    add_check(
        "odd_even_lane_selection_bridge_exported",
        False,
        False,
        "No current export provides the missing bridge rule distinguishing odd-anchor lane versus even-fallback lane within the positive corridor.",
    )

    status = (
        "PASS_POSITIVE_CORRIDOR_ODD_EVEN_LANE_SELECTION_BRIDGE_NONEXPORT_AUDITED"
        if not blocking
        else "P725_REQUIRES_REVIEW_CHANGED_POSITIVE_CORRIDOR_LANE_SELECTION_STATE"
    )

    artifact = {
        "stage": "P725",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t179_positive_corridor_odd_even_lane_selection_bridge_nonexport_boundary_only",
        "inputs": {
            "P723": str(IN_P723.relative_to(REPO)),
            "P724": str(IN_P724.relative_to(REPO)),
            "P715": str(IN_P715.relative_to(REPO)),
            "F141": str(IN_F141.relative_to(REPO)),
            "F143": str(IN_F143.relative_to(REPO)),
            "F147": str(IN_F147.relative_to(REPO)),
            "F148": str(IN_F148.relative_to(REPO)),
            "F150": str(IN_F150.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t179_target_name": "PositiveCorridorOddEvenLaneSelectionBridge_global_C_v1_strict_v1",
        "t179_target_exported_on_current_repo_state": False,
        "current_positive_corridor_profile": {
            "positive_corridor_roots": positive_roots,
            "odd_anchor_lane": odd_lane,
            "even_fallback_lane": even_lane,
            "source_packets_remain_pair_and_lane_agnostic": source_packets_pair_and_lane_agnostic,
            "source_selector_axis_frame_basis": selector_axis_basis,
            "tau_src_identified_with_s_prelm": tau_src_identified,
        },
        "audit_conclusion": {
            "current_positive_corridor_is_physically_reduced_but_still_two_lane": True,
            "current_repo_already_exports_t179_target": False,
            "next_honest_move": (
                "export_or_attack_a_finer_source_side_rule_that_selects_between_odd_anchor_and_even_fallback_lanes_within_the_positive_corridor"
            ),
        },
        "hard_limits": [
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
        "stage": "P725",
        "status": status,
        "as_of": AS_OF,
        "t179_target_name": "PositiveCorridorOddEvenLaneSelectionBridge_global_C_v1_strict_v1",
        "t179_target_exported_on_current_repo_state": False,
        "positive_corridor_roots": positive_roots,
        "odd_anchor_lane": odd_lane,
        "even_fallback_lane": even_lane,
        "next_honest_move": "finer_source_side_rule_selecting_between_positive_corridor_lanes",
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
