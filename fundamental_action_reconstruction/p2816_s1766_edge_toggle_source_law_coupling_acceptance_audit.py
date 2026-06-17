#!/usr/bin/env python3
"""P2816/S1766: source-law/coupling acceptance audit after P2815.

P2815 completed finite carrier separation, but explicitly did not export a
source law or K/L_total coupling.  P2816 tests exactly one explicit candidate:
the dimensionless ordinal functional induced by the complete P2815 carrier

    Q_P2815_rank(G) = rank of G in the sorted P2815 separating-key order.

The typed map candidate is the minimal scalar insertion
G -> Q_P2815_rank(G) -> c*Q_P2815_rank(G) in K or L_total.  This audit is deliberately
falsifiable: it checks carrier separation, dimensional normalization, exported
source-law status, and typed coupling status.  It is not a closure claim.
"""
from __future__ import annotations

import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import EXPECTED_GRAPH_COUNT, SCD, sha
from p2812_s1762_two_wl_refined_collision_audit import read_json

GEN = ROOT / "generated"
P2815 = GEN / "p2815_s1765_edge_toggle_response_residual_audit.json"
OUT = GEN / "p2816_s1766_edge_toggle_source_law_coupling_acceptance_audit.json"
MD = GEN / "p2816_s1766_edge_toggle_source_law_coupling_acceptance_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def build_audit(p2815: dict[str, Any]) -> dict[str, Any]:
    """Audit one explicit source-functional candidate induced by P2815.

    We intentionally do not add a new carrier-refinement invariant.  The
    candidate is the ordinal/rank functional on the already separated P2815
    carrier.  If complete separation alone licensed a source law, this would be
    the strongest available post-P2815 test; the acceptance matrix then asks
    whether the required normalization and coupling data are exported.
    """
    edge_audit = p2815["edge_toggle_response_audit"]
    class_count = edge_audit["edge_toggle_refined_class_count"]
    collision_count = edge_audit["edge_toggle_collision_class_count"]
    return {
        "candidate_functional": "Q_P2815_rank(G)=canonical ordinal/rank on the complete P2815 edge-toggle separating carrier",
        "typed_map_candidate": "G -> Q_P2815_rank(G) -> dimensionless coefficient c*Q in K or L_total",
        "decoded_graph_count": edge_audit["decoded_graph_count"],
        "input_p2815_refined_class_count": class_count,
        "input_p2815_collision_class_count": collision_count,
        "q_rank_class_count": class_count,
        "q_rank_collision_class_count": collision_count,
        "q_rank_collision_graph_count": edge_audit["edge_toggle_collision_graph_count"],
        "remaining_defect_canonical_minus_q_rank": EXPECTED_GRAPH_COUNT - class_count,
        "rank_domain_min": 0,
        "rank_domain_max": class_count - 1,
        "computational_basis": "Uses the committed P2815 finite carrier-separation certificate rather than adding another graph invariant.",
    }

def acceptance_matrix(audit: dict[str, Any], p2815: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2815_complete_carrier_available": p2815["acceptance_matrix"]["accepted_as_complete_canonical_source_carrier"],
        "exactly_one_candidate_functional_tested": True,
        "functional_is_deterministic_and_reproducible": audit["decoded_graph_count"] == EXPECTED_GRAPH_COUNT,
        "functional_separates_full_carrier": audit["q_rank_class_count"] == EXPECTED_GRAPH_COUNT,
        "dimensionful_normalization_exported": False,
        "strict_graph_source_law_exported": False,
        "typed_graph_to_K_or_Ltotal_coupling_exported": False,
        "selector_bridge_or_role_transfer_imported": False,
    }
    accepted = all(facts[key] for key in [
        "p2815_complete_carrier_available",
        "exactly_one_candidate_functional_tested",
        "functional_is_deterministic_and_reproducible",
        "functional_separates_full_carrier",
        "dimensionful_normalization_exported",
        "strict_graph_source_law_exported",
        "typed_graph_to_K_or_Ltotal_coupling_exported",
    ]) and not facts["selector_bridge_or_role_transfer_imported"]
    return {
        "facts": facts,
        "accepted_as_source_law_coupling": accepted,
        "accepted_as_bounded_candidate_rejection": not accepted,
        "missing_for_promotion": [key for key, value in facts.items() if not value and key != "selector_bridge_or_role_transfer_imported"],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["q_rank_source_candidate_audit"]
    acceptance = payload["acceptance_matrix"]
    lines = [
        "# P2816/S1766 edge-toggle source-law/coupling acceptance audit", "", f"Status: `{payload['status']}`", "",
        "## Candidate", audit["candidate_functional"], "", "## Typed map candidate", audit["typed_map_candidate"], "",
        "## Finite counts",
        f"- decoded_graph_count={audit['decoded_graph_count']}",
        f"- q_rank_class_count={audit['q_rank_class_count']}",
        f"- q_rank_collision_class_count={audit['q_rank_collision_class_count']}",
        f"- q_rank_collision_graph_count={audit['q_rank_collision_graph_count']}",
        f"- rank_domain_max={audit['rank_domain_max']}",
        f"- remaining_defect_canonical_minus_q_rank={audit['remaining_defect_canonical_minus_q_rank']}", "",
        "## Acceptance", f"- accepted_as_source_law_coupling={acceptance['accepted_as_source_law_coupling']}",
        f"- accepted_as_bounded_candidate_rejection={acceptance['accepted_as_bounded_candidate_rejection']}", "",
        "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2815 = read_json(P2815)
    audit = build_audit(p2815)
    payload: dict[str, Any] = {
        "status": "P2816_EDGE_TOGGLE_SOURCE_FUNCTIONAL_REJECTED_NO_COUPLING_NO_CLOSURE",
        "input_hashes": {"P2815": sha(P2815), "16_4_4.scd": sha(SCD)},
        "input_statuses": {"P2815": p2815.get("status")},
        "q_rank_source_candidate_audit": audit,
        "decision": {
            "negative_export_flags": {
                "strict_graph_source_law_exported": False,
                "typed_coupling_to_K_or_Ltotal_exported": False,
                "dimensionful_normalization_exported": False,
                "role_bearing_ltotal_promoted": False,
                "bridge_closure_exported": False,
                "role_transfer_started": False,
                "selector_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2816 tests exactly one explicit post-P2815 source-law/coupling candidate: Q_P2815_rank and a minimal graph-to-K/L_total coefficient insertion.  The finite carrier input is complete and deterministic, but the rank is still an ordinal graph label: no dimensionful normalization, strict graph-source law, variational meaning, or typed K/L_total coupling theorem is exported.  Thus P2816 is a bounded coupling/source-law rejection, not dynamics or closure.",
            "next_honest_step": "Do not replay carrier refinement or promote ordinal carrier ranks.  The next honest move is exactly one non-ordinal source-law candidate: an explicit typed graph observable with dimension/normalization data and a real graph-to-K or graph-to-L_total formula; audit it against separability, variational meaning, and no selector/bridge/role-transfer import.  If no such formula is supplied, preserve the no-coupling boundary.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit, p2815)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2816/S1766 edge-toggle source-law coupling acceptance audit", "## P2816/S1766 edge-toggle source-law coupling acceptance audit\n\n`P2816/S1766` tests exactly one explicit post-P2815 candidate source functional, `Q_P2815_rank(G)`, plus the minimal typed-map candidate `G -> Q_P2815_rank(G) -> c*Q_P2815_rank(G)` for `K` or `L_total`.  The P2815 carrier input is complete and deterministic, but the induced rank is an ordinal graph label and no dimensionful normalization, strict graph-source law, variational meaning, or typed `K`/`L_total` coupling theorem is exported.  Therefore P2816 rejects this candidate and adds no closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2816/S1766 edge-toggle source-law coupling Ltotal guard", "## P2816/S1766 edge-toggle source-law coupling Ltotal guard\n\n`P2816/S1766` adds no variational term to `L_total`.  The tested `Q_P2815_rank` source-functional candidate is a finite ordinal graph diagnostic with no exported dimensionful normalization, Euler-Lagrange role, or typed graph-to-`K`/`L_total` coupling theorem.\n")
    append_once(AGENTS, "Current edge-toggle source-law coupling audit guardrail (P2816/S1766, 2026-06-16)", "## Current edge-toggle source-law coupling audit guardrail (P2816/S1766, 2026-06-16)\n\n- P2816 tests exactly one explicit post-P2815 source-law/coupling candidate: `Q_P2815_rank(G)` with a minimal graph-to-`K`/`L_total` coefficient insertion.\n- The candidate uses the complete P2815 carrier, but it is only an ordinal graph label and exports no dimensionful normalization, strict graph-source law, variational meaning, or typed coupling theorem; it is rejected as a source-law/coupling promotion.\n- Do not replay ordinal ranks or carrier-refinement as dynamics.  A next admissible move must provide one non-ordinal formula using explicit typed graph-to-`K` or graph-to-`L_total` normalization/coupling data, while preserving bridge, role-transfer, selector, `L_total`, and ToE blocks until such an audit succeeds.\n")
    return payload


if __name__ == "__main__":
    main()
