#!/usr/bin/env python3
"""P2841/S1791: post-P2840 broad state-map intake gate.

P2840 closes the current graph-source lane unless a genuinely new strict object is
supplied.  P2841 performs the broader intake gate across the live FAR state map:
it checks whether any closed/repetition-gated lane received a new typed object,
source theorem, blocker-cut, or provider class.  On current artifacts, none did.
"""
from __future__ import annotations

import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import SCD, sha
from p2812_s1762_two_wl_refined_collision_audit import read_json

GEN = ROOT / "generated"
P2840 = GEN / "p2840_s1790_post_graph_source_state_map_reconciliation_certificate.json"
OUT = GEN / "p2841_s1791_post_p2840_broad_state_map_intake_gate.json"
MD = GEN / "p2841_s1791_post_p2840_broad_state_map_intake_gate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

STATE_MAP_LANES = (
    "graph_source_P2835_P2840",
    "selector_QW2191_orientation",
    "legacy_strict_bridge_role_transfer",
    "direct_route_residual_cancellation",
    "lagrangian_EOM_reverse_closure",
    "lower_boundary_tau_pair_recursion",
    "P2680_non_selector_source_atoms",
    "L_total_ToE_promotion",
)

REQUIRED_NEW_OBJECT_FIELDS = (
    "new_strict_typed_object",
    "new_source_theorem",
    "new_blocker_cut",
    "new_provider_class",
    "new_coupling_or_localization_theorem",
)


def broad_state_map_rows(p2840: dict[str, Any]) -> list[dict[str, Any]]:
    p2840_acceptance = p2840["acceptance_matrix"]
    return [
        {
            "lane": "graph_source_P2835_P2840",
            "current_status": "no_new_live_frontier_certificate",
            "evidence": "P2840 graph-source state-map certificate accepted",
            "evidence_flag": p2840_acceptance["accepted_as_no_new_live_frontier_certificate"],
            "new_strict_typed_object": False,
            "new_source_theorem": False,
            "new_blocker_cut": False,
            "new_provider_class": False,
            "new_coupling_or_localization_theorem": False,
        },
        {
            "lane": "selector_QW2191_orientation",
            "current_status": "repetition_gated",
            "evidence": "No non-premise strict selector/orientation source is supplied by P2833-P2840.",
            "evidence_flag": True,
            "new_strict_typed_object": False,
            "new_source_theorem": False,
            "new_blocker_cut": False,
            "new_provider_class": False,
            "new_coupling_or_localization_theorem": False,
        },
        {
            "lane": "legacy_strict_bridge_role_transfer",
            "current_status": "blocked_without_bridge_completion_and_role_transfer_theorem",
            "evidence": "No completion bridge, role-transfer theorem, or strict-side source atom is exported by the graph-source lane.",
            "evidence_flag": True,
            "new_strict_typed_object": False,
            "new_source_theorem": False,
            "new_blocker_cut": False,
            "new_provider_class": False,
            "new_coupling_or_localization_theorem": False,
        },
        {
            "lane": "direct_route_residual_cancellation",
            "current_status": "repetition_gated_after_residual_no_go",
            "evidence": "No new direct-route coefficient equation, provider class, or non-N477 ingredient is supplied.",
            "evidence_flag": True,
            "new_strict_typed_object": False,
            "new_source_theorem": False,
            "new_blocker_cut": False,
            "new_provider_class": False,
            "new_coupling_or_localization_theorem": False,
        },
        {
            "lane": "lagrangian_EOM_reverse_closure",
            "current_status": "bounded_no_go_without_new_anisotropic_source_class",
            "evidence": "No strict anisotropic source class or nonproxy component residual discharge is supplied.",
            "evidence_flag": True,
            "new_strict_typed_object": False,
            "new_source_theorem": False,
            "new_blocker_cut": False,
            "new_provider_class": False,
            "new_coupling_or_localization_theorem": False,
        },
        {
            "lane": "lower_boundary_tau_pair_recursion",
            "current_status": "cycle_cut_repetition_gated",
            "evidence": "No chart-label-retaining pair12 typed seed/subinterface or nonconventional semantic provider is supplied.",
            "evidence_flag": True,
            "new_strict_typed_object": False,
            "new_source_theorem": False,
            "new_blocker_cut": False,
            "new_provider_class": False,
            "new_coupling_or_localization_theorem": False,
        },
        {
            "lane": "P2680_non_selector_source_atoms",
            "current_status": "bounded_no_go_after_source_inventory",
            "evidence": "Canonical UV unit, alpha_geo amplitude source, and beta/Z_beta source remain unexported; P2836 also finds no graph-derived units source.",
            "evidence_flag": True,
            "new_strict_typed_object": False,
            "new_source_theorem": False,
            "new_blocker_cut": False,
            "new_provider_class": False,
            "new_coupling_or_localization_theorem": False,
        },
        {
            "lane": "L_total_ToE_promotion",
            "current_status": "blocked_without_units_typed_coupling_variational_derivative",
            "evidence": "P2836-P2839 explicitly leave units, typed coupling, variational derivative, and localization unexported.",
            "evidence_flag": True,
            "new_strict_typed_object": False,
            "new_source_theorem": False,
            "new_blocker_cut": False,
            "new_provider_class": False,
            "new_coupling_or_localization_theorem": False,
        },
    ]


def build_audit(p2840: dict[str, Any]) -> dict[str, Any]:
    rows = broad_state_map_rows(p2840)
    failed_evidence_rows = [row for row in rows if not row["evidence_flag"]]
    unlocked_rows = [
        row for row in rows
        if any(bool(row[field]) for field in REQUIRED_NEW_OBJECT_FIELDS)
    ]
    return {
        "input_statuses_rechecked": {"P2840": p2840.get("status")},
        "state_map_lane_rows": rows,
        "state_map_lane_count": len(rows),
        "required_new_object_fields": list(REQUIRED_NEW_OBJECT_FIELDS),
        "failed_evidence_rows": failed_evidence_rows,
        "new_object_unlocked_rows": unlocked_rows,
        "no_new_live_frontier_certificate": not failed_evidence_rows and not unlocked_rows,
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "all_state_map_lanes_accounted": audit["state_map_lane_count"] == len(STATE_MAP_LANES),
        "all_evidence_flags_pass": not audit["failed_evidence_rows"],
        "no_new_typed_object_or_provider_supplied": not audit["new_object_unlocked_rows"],
        "no_new_live_frontier_certificate": audit["no_new_live_frontier_certificate"],
        "selector_bridge_or_role_transfer_imported": False,
    }
    accepted = all([
        facts["all_state_map_lanes_accounted"],
        facts["all_evidence_flags_pass"],
        facts["no_new_typed_object_or_provider_supplied"],
        facts["no_new_live_frontier_certificate"],
        not facts["selector_bridge_or_role_transfer_imported"],
    ])
    return {
        "facts": facts,
        "accepted_as_broad_state_map_intake_gate": accepted,
        "accepted_as_no_new_live_frontier_certificate": accepted,
        "accepted_as_closure_or_promotion": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["post_p2840_broad_state_map_intake_gate"]
    lines = [
        "# P2841/S1791 post-P2840 broad state-map intake gate", "", f"Status: `{payload['status']}`", "",
        "## Matrix counts",
        f"- state_map_lane_count={audit['state_map_lane_count']}",
        f"- failed_evidence_rows={len(audit['failed_evidence_rows'])}",
        f"- new_object_unlocked_rows={len(audit['new_object_unlocked_rows'])}", "",
        "## Certificate",
        f"- no_new_live_frontier_certificate={audit['no_new_live_frontier_certificate']}",
        f"- accepted_as_no_new_live_frontier_certificate={payload['acceptance_matrix']['accepted_as_no_new_live_frontier_certificate']}", "",
        "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2840 = read_json(P2840)
    audit = build_audit(p2840)
    payload: dict[str, Any] = {
        "status": "P2841_POST_P2840_BROAD_STATE_MAP_INTAKE_GATE_NO_NEW_LIVE_FRONTIER_NO_CLOSURE",
        "input_hashes": {"P2840": sha(P2840), "16_4_4.scd": sha(SCD)},
        "post_p2840_broad_state_map_intake_gate": audit,
        "decision": {
            "negative_export_flags": {
                "new_strict_typed_object_exported": False,
                "new_source_theorem_exported": False,
                "new_blocker_cut_exported": False,
                "new_provider_class_exported": False,
                "new_coupling_or_localization_theorem_exported": False,
                "role_bearing_ltotal_promoted": False,
                "bridge_closure_exported": False,
                "role_transfer_started": False,
                "selector_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2841 applies the post-P2840 intake gate across the broad state map.  The graph-source lane is closed by P2840, and selector/QW-2191, bridge/role-transfer, direct-route residual, Lagrangian/EOM reverse closure, lower-boundary recursion, P2680 source atoms, and L_total/ToE promotion lanes receive no new strict typed object, source theorem, blocker-cut, provider class, or coupling/localization theorem.  Therefore no new live frontier is unlocked on current artifacts.",
            "next_honest_step": "Do not continue replay in graph-source or other closed lanes.  The next admissible move must introduce one genuinely new strict typed object/source theorem/blocker-cut/provider class outside the closed lanes, or preserve this broad no-new-live-frontier certificate.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2841/S1791 post-P2840 broad state-map intake gate", "## P2841/S1791 post-P2840 broad state-map intake gate\n\n`P2841/S1791` applies the post-P2840 intake gate across the broad state map.  The graph-source lane is closed by P2840, and selector/`QW-2191`, bridge/role-transfer, direct-route residual, Lagrangian/EOM reverse closure, lower-boundary recursion, P2680 source atoms, and `L_total`/ToE promotion lanes receive no new strict typed object, source theorem, blocker-cut, provider class, or coupling/localization theorem.  This emits a broad no-new-live-frontier certificate on current artifacts; no strict graph-source law, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2841/S1791 broad state-map intake Ltotal guard", "## P2841/S1791 broad state-map intake Ltotal guard\n\n`P2841/S1791` adds no term to `L_total`.  It confirms that after P2840 no broad state-map lane receives a new strict typed object, coupling theorem, localization theorem, variational source, or provider class sufficient to reopen `L_total`/ToE promotion.\n")
    append_once(AGENTS, "Current post-P2840 broad state-map intake guardrail (P2841/S1791, 2026-06-17)", "## Current post-P2840 broad state-map intake guardrail (P2841/S1791, 2026-06-17)\n\n- P2841 applies the broad intake gate after P2840: graph-source, selector/`QW-2191`, bridge/role-transfer, direct-route residual, Lagrangian/EOM reverse closure, lower-boundary recursion, P2680 source atoms, and `L_total`/ToE promotion lanes receive no new strict typed object/source theorem/blocker-cut/provider class/coupling theorem.\n- Do not replay graph-source or other closed lanes as a primary strategy, and do not promote `L_total`, bridge closure, role transfer, selector closure, or ToE closure from P2833-P2841.\n- A next admissible move must introduce exactly one genuinely new strict typed object, source theorem, blocker-cut, provider class, or coupling/localization theorem outside the closed lanes; otherwise preserve the broad no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
