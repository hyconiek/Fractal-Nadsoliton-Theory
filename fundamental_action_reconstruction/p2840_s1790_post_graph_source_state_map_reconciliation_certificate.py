#!/usr/bin/env python3
"""P2840/S1790: post graph-source state-map reconciliation certificate.

P2839 recommends a post-graph-source state-map/no-new-live-frontier certificate
unless a genuinely new strict localization object or coupling theorem is supplied.
P2840 performs that bounded reconciliation across the finite graph-source lane and
adjacent closed lanes, without manufacturing a closure.
"""
from __future__ import annotations

import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import SCD, sha
from p2812_s1762_two_wl_refined_collision_audit import read_json

GEN = ROOT / "generated"
P2835 = GEN / "p2835_s1785_combined_witness_source_law_theorem_obligation_audit.json"
P2836 = GEN / "p2836_s1786_combined_graph_functional_units_normalization_obligation_audit.json"
P2837 = GEN / "p2837_s1787_combined_graph_functional_typed_domain_codomain_audit.json"
P2838 = GEN / "p2838_s1788_combined_graph_functional_variational_derivative_obstruction_audit.json"
P2839 = GEN / "p2839_s1789_graph_bit_localization_map_obstruction_audit.json"
OUT = GEN / "p2840_s1790_post_graph_source_state_map_reconciliation_certificate.json"
MD = GEN / "p2840_s1790_post_graph_source_state_map_reconciliation_certificate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

GRAPH_SOURCE_LANE_ROWS = (
    "finite_combined_separator",
    "target_independent_units_normalization",
    "typed_domain_codomain_map",
    "formal_variational_derivative",
    "evaluation_pullback_localization",
    "unit_bearing_coupling_source_law",
)

ADJACENT_REPLAY_GATED_ROWS = (
    "selector_QW2191_replay",
    "legacy_strict_bridge_role_transfer",
    "L_total_or_ToE_promotion",
    "direct_route_residual_replay",
    "lagrangian_EOM_reverse_closure_replay",
)


def graph_source_lane_matrix(p2835: dict[str, Any], p2836: dict[str, Any], p2837: dict[str, Any], p2838: dict[str, Any], p2839: dict[str, Any]) -> list[dict[str, Any]]:
    p2835_acceptance = p2835["acceptance_matrix"]
    p2836_acceptance = p2836["acceptance_matrix"]
    p2837_acceptance = p2837["acceptance_matrix"]
    p2838_acceptance = p2838["acceptance_matrix"]
    p2839_acceptance = p2839["acceptance_matrix"]
    rows = [
        {
            "lane": "finite_combined_separator",
            "source": "P2835",
            "status": "positive_finite_witness_closed",
            "evidence_flag": p2835_acceptance["accepted_as_combined_finite_separator"],
            "live_frontier_unlocked": False,
            "reason": "Full-carrier finite separation is already achieved; more separation replay is not a theorem-level coupling move.",
        },
        {
            "lane": "target_independent_units_normalization",
            "source": "P2836",
            "status": "bounded_no_go",
            "evidence_flag": p2836_acceptance["accepted_as_units_normalization_no_go"],
            "live_frontier_unlocked": False,
            "reason": "Finite dimensionless averaging exists, but no target-independent physical units or scale quotient is exported.",
        },
        {
            "lane": "typed_domain_codomain_map",
            "source": "P2837",
            "status": "bounded_no_go",
            "evidence_flag": p2837_acceptance["accepted_as_typed_domain_codomain_no_go"],
            "live_frontier_unlocked": False,
            "reason": "Candidate maps into K/L_total/source/coefficient fail required typed premises.",
        },
        {
            "lane": "formal_variational_derivative",
            "source": "P2838",
            "status": "bounded_no_go",
            "evidence_flag": p2838_acceptance["accepted_as_variational_derivative_no_go"],
            "live_frontier_unlocked": False,
            "reason": "Finite graph differences do not export Euler-Lagrange/functional derivatives into K or L_total.",
        },
        {
            "lane": "evaluation_pullback_localization",
            "source": "P2839",
            "status": "bounded_no_go",
            "evidence_flag": p2839_acceptance["accepted_as_localization_no_go"],
            "live_frontier_unlocked": False,
            "reason": "No canonical label-safe graph-bit-to-field localization/pullback is exported.",
        },
        {
            "lane": "unit_bearing_coupling_source_law",
            "source": "P2835-P2839",
            "status": "not_exported",
            "evidence_flag": True,
            "live_frontier_unlocked": False,
            "reason": "Coupling remains blocked by the combined missing units, typed map, variational derivative, and localization premises.",
        },
    ]
    return rows


def adjacent_replay_gate_matrix() -> list[dict[str, Any]]:
    return [
        {
            "lane": "selector_QW2191_replay",
            "status": "repetition_gated",
            "new_object_supplied": False,
            "reason": "Graph-source audits import no non-premise selector/orientation source and do not discharge QW-2191.",
        },
        {
            "lane": "legacy_strict_bridge_role_transfer",
            "status": "repetition_gated",
            "new_object_supplied": False,
            "reason": "No bridge-completion theorem or role-transfer theorem is exported by graph-source witnesses.",
        },
        {
            "lane": "L_total_or_ToE_promotion",
            "status": "blocked",
            "new_object_supplied": False,
            "reason": "No unit-bearing coupling, typed localization, or variational derivative into L_total is exported.",
        },
        {
            "lane": "direct_route_residual_replay",
            "status": "repetition_gated",
            "new_object_supplied": False,
            "reason": "The graph-source finite separator is not a new direct-route residual cancellation theorem.",
        },
        {
            "lane": "lagrangian_EOM_reverse_closure_replay",
            "status": "repetition_gated",
            "new_object_supplied": False,
            "reason": "No new anisotropic source class, EOM closure theorem, or nonproxy component residual discharge is exported.",
        },
    ]


def build_audit(p2835: dict[str, Any], p2836: dict[str, Any], p2837: dict[str, Any], p2838: dict[str, Any], p2839: dict[str, Any]) -> dict[str, Any]:
    graph_rows = graph_source_lane_matrix(p2835, p2836, p2837, p2838, p2839)
    adjacent_rows = adjacent_replay_gate_matrix()
    live_graph_rows = [row for row in graph_rows if row["live_frontier_unlocked"]]
    failed_evidence_rows = [row for row in graph_rows if not row["evidence_flag"]]
    new_adjacent_rows = [row for row in adjacent_rows if row["new_object_supplied"]]
    return {
        "input_statuses_rechecked": {
            "P2835": p2835.get("status"),
            "P2836": p2836.get("status"),
            "P2837": p2837.get("status"),
            "P2838": p2838.get("status"),
            "P2839": p2839.get("status"),
        },
        "graph_source_lane_rows": graph_rows,
        "adjacent_replay_gate_rows": adjacent_rows,
        "graph_source_lane_row_count": len(graph_rows),
        "adjacent_replay_gate_row_count": len(adjacent_rows),
        "live_graph_source_frontier_rows": live_graph_rows,
        "failed_evidence_rows": failed_evidence_rows,
        "new_adjacent_object_rows": new_adjacent_rows,
        "no_new_live_frontier_certificate": not live_graph_rows and not failed_evidence_rows and not new_adjacent_rows,
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "all_graph_source_rows_accounted": audit["graph_source_lane_row_count"] == len(GRAPH_SOURCE_LANE_ROWS),
        "all_adjacent_replay_rows_accounted": audit["adjacent_replay_gate_row_count"] == len(ADJACENT_REPLAY_GATED_ROWS),
        "all_input_evidence_flags_pass": not audit["failed_evidence_rows"],
        "no_live_graph_source_frontier_unlocked": not audit["live_graph_source_frontier_rows"],
        "no_new_adjacent_object_supplied": not audit["new_adjacent_object_rows"],
        "selector_bridge_or_role_transfer_imported": False,
    }
    accepted = all([
        facts["all_graph_source_rows_accounted"],
        facts["all_adjacent_replay_rows_accounted"],
        facts["all_input_evidence_flags_pass"],
        facts["no_live_graph_source_frontier_unlocked"],
        facts["no_new_adjacent_object_supplied"],
        not facts["selector_bridge_or_role_transfer_imported"],
    ])
    return {
        "facts": facts,
        "accepted_as_post_graph_source_state_map_reconciliation": accepted,
        "accepted_as_no_new_live_frontier_certificate": accepted and audit["no_new_live_frontier_certificate"],
        "accepted_as_closure_or_promotion": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["post_graph_source_state_map_reconciliation_certificate"]
    lines = [
        "# P2840/S1790 post graph-source state-map reconciliation certificate", "", f"Status: `{payload['status']}`", "",
        "## Matrix counts",
        f"- graph_source_lane_row_count={audit['graph_source_lane_row_count']}",
        f"- adjacent_replay_gate_row_count={audit['adjacent_replay_gate_row_count']}",
        f"- failed_evidence_rows={len(audit['failed_evidence_rows'])}",
        f"- live_graph_source_frontier_rows={len(audit['live_graph_source_frontier_rows'])}",
        f"- new_adjacent_object_rows={len(audit['new_adjacent_object_rows'])}", "",
        "## Certificate",
        f"- no_new_live_frontier_certificate={audit['no_new_live_frontier_certificate']}",
        f"- accepted_as_no_new_live_frontier_certificate={payload['acceptance_matrix']['accepted_as_no_new_live_frontier_certificate']}", "",
        "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2835 = read_json(P2835)
    p2836 = read_json(P2836)
    p2837 = read_json(P2837)
    p2838 = read_json(P2838)
    p2839 = read_json(P2839)
    audit = build_audit(p2835, p2836, p2837, p2838, p2839)
    payload: dict[str, Any] = {
        "status": "P2840_POST_GRAPH_SOURCE_STATE_MAP_NO_NEW_LIVE_FRONTIER_CERTIFICATE_NO_CLOSURE",
        "input_hashes": {"P2835": sha(P2835), "P2836": sha(P2836), "P2837": sha(P2837), "P2838": sha(P2838), "P2839": sha(P2839), "16_4_4.scd": sha(SCD)},
        "post_graph_source_state_map_reconciliation_certificate": audit,
        "decision": {
            "negative_export_flags": {
                "new_graph_source_frontier_unlocked": False,
                "strict_graph_source_law_exported": False,
                "typed_coupling_to_K_or_Ltotal_exported": False,
                "formal_variational_derivative_exported": False,
                "evaluation_pullback_localization_map_exported": False,
                "role_bearing_ltotal_promoted": False,
                "bridge_closure_exported": False,
                "role_transfer_started": False,
                "selector_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2840 reconciles the post-P2835 graph-source state map.  The finite separator is closed as a positive finite witness, while units/normalization, typed domain/codomain, formal variational derivative, localization/pullback, and unit-bearing coupling source-law rows remain bounded no-go/not-exported.  Adjacent selector, bridge/role-transfer, L_total/ToE, direct-route, and Lagrangian/EOM replay lanes receive no new strict object.  Therefore the graph-source lane has no new live frontier on current artifacts.",
            "next_honest_step": "Do not continue graph-source replay unless a genuinely new strict localization object, unit-bearing coupling theorem, or graph-to-field source law is supplied.  The next admissible research move must introduce a new typed object outside the closed graph-source lane, or preserve this no-new-live-frontier certificate.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2840/S1790 post graph-source state-map reconciliation certificate", "## P2840/S1790 post graph-source state-map reconciliation certificate\n\n`P2840/S1790` reconciles the post-P2835 graph-source state map.  The combined finite separator remains a positive finite witness, but units/normalization, typed domain/codomain, formal variational derivative, graph-bit localization/pullback, and unit-bearing coupling source-law rows remain bounded no-go or not exported.  Adjacent selector, bridge/role-transfer, `L_total`/ToE, direct-route, and Lagrangian/EOM replay lanes receive no new strict object.  This emits a no-new-live-frontier certificate for the current graph-source lane; no strict graph-source law, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2840/S1790 post graph-source state-map Ltotal guard", "## P2840/S1790 post graph-source state-map Ltotal guard\n\n`P2840/S1790` adds no term to `L_total`.  It reconciles P2835-P2839 and finds no new live graph-source frontier: finite separation is positive, but units, typed coupling, variational derivative, localization, and unit-bearing source law remain blocked.\n")
    append_once(AGENTS, "Current post graph-source state-map guardrail (P2840/S1790, 2026-06-17)", "## Current post graph-source state-map guardrail (P2840/S1790, 2026-06-17)\n\n- P2840 reconciles P2835-P2839: finite graph-source separation is closed as a positive finite witness, while units/normalization, typed domain/codomain, formal variational derivative, localization/pullback, and unit-bearing coupling source-law rows remain bounded no-go/not exported.\n- Adjacent selector, bridge/role-transfer, `L_total`/ToE, direct-route, and Lagrangian/EOM replay lanes receive no new strict object from the graph-source audits.\n- Do not continue graph-source replay or promote `L_total`, bridge closure, role transfer, selector closure, or ToE closure.  A next admissible move must introduce a genuinely new strict localization object, coupling theorem, graph-to-field source law, or a new typed object outside the closed graph-source lane; otherwise preserve the no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
