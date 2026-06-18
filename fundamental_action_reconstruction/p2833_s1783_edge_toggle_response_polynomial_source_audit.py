#!/usr/bin/env python3
"""P2833/S1783: exact edge-toggle response polynomial source audit.

P2832 rejects the low-order common-neighbor/cycle profile as too coarse for a
source-law/coupling promotion.  P2833 tests one higher-order non-label formula:
the full edge-toggle response polynomial over every unordered vertex pair,
recording whether the pair is an edge, the triangle response coefficient, and
the four-cycle response coefficient.  This is a finite variational-response
candidate, not a digest label, hash, rank, selector import, bridge import, or
closure claim.
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import SCD, decode_scd_bytes, sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2819_s1769_p2815_digest_blocker_cut_response_audit import P2818
from p2832_s1782_common_neighbor_cycle_source_formula_audit import adjacency_masks

GEN = ROOT / "generated"
P2832 = GEN / "p2832_s1782_common_neighbor_cycle_source_formula_audit.json"
OUT = GEN / "p2833_s1783_edge_toggle_response_polynomial_source_audit.json"
MD = GEN / "p2833_s1783_edge_toggle_response_polynomial_source_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
EXPECTED_GRAPH_COUNT = 16828


def edge_toggle_response_polynomial(graph: tuple[tuple[int, ...], ...]) -> tuple[Any, ...]:
    """Return the non-label multiset response polynomial for all vertex pairs.

    For each unordered pair `{i,j}` the row is
    `(edge_present, delta_triangle_abs, delta_4cycle_abs)`.  Adding or removing
    an edge changes the triangle count by the common-neighbor count and changes
    the 4-cycle count by the number of length-3 paths between `i` and `j` after
    excluding the pair itself.  The sorted row histogram is invariant under
    relabeling and is a compact exact edge-toggle response polynomial.
    """
    masks = adjacency_masks(graph)
    rows: Counter[tuple[int, int, int]] = Counter()
    n = len(masks)
    for i in range(n):
        for j in range(i + 1, n):
            edge_present = (masks[i] >> j) & 1
            triangle_response = (masks[i] & masks[j]).bit_count()
            ni = masks[i] & ~(1 << j)
            nj = masks[j] & ~(1 << i)
            four_cycle_response = 0
            x = ni
            while x:
                lsb = x & -x
                a = lsb.bit_length() - 1
                x -= lsb
                four_cycle_response += (masks[a] & nj).bit_count()
            rows[(edge_present, triangle_response, four_cycle_response)] += 1
    return tuple(sorted(rows.items()))


def build_audit(p2832: dict[str, Any]) -> dict[str, Any]:
    graphs, _ = decode_scd_bytes(SCD.read_bytes())
    classes: dict[tuple[Any, ...], list[int]] = defaultdict(list)
    for index, graph in enumerate(graphs):
        classes[edge_toggle_response_polynomial(graph)].append(index)
    residual = {key: indices for key, indices in classes.items() if len(indices) > 1}
    hist = Counter(len(indices) for indices in classes.values())
    largest = sorted((len(indices), indices[:12], key) for key, indices in residual.items())[-12:]
    return {
        "candidate_formula": "Q_toggle(G)=multiset_{unordered pairs {i,j}} (1_{ij in E}, |N(i)∩N(j)|, #{length-3 paths i-a-b-j excluding the toggled pair}); equivalently the exact triangle and 4-cycle edge-toggle response polynomial on the fixed 16-node carrier",
        "decoded_graph_count": len(graphs),
        "input_p2832_status": p2832.get("status"),
        "higher_order_non_label_graph_formula_exported_for_candidate": True,
        "digest_label_hash_or_rank_used": False,
        "edge_toggle_variational_response_coefficients_computed": True,
        "dimensionless_normalization_available_for_candidate": True,
        "candidate_class_count": len(classes),
        "candidate_collision_class_count": len(residual),
        "candidate_collision_graph_count": sum(len(indices) for indices in residual.values()),
        "candidate_max_class_size": max(hist),
        "candidate_defect_after_formula": len(graphs) - len(classes),
        "candidate_class_size_histogram": dict(sorted(hist.items())),
        "largest_collision_samples": [
            {"class_size": size, "sample_graph_indices_0_based": sample, "response_polynomial": key}
            for size, sample, key in reversed(largest)
        ],
    }


def acceptance_matrix(audit: dict[str, Any], p2818: dict[str, Any], p2832: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2818_local_response_rejected": p2818["acceptance_matrix"]["accepted_as_bounded_candidate_rejection"],
        "p2832_low_order_formula_rejected": p2832["acceptance_matrix"]["accepted_as_bounded_formula_no_go"],
        "exactly_one_higher_order_non_label_formula_tested": True,
        "digest_label_hash_or_rank_used": audit["digest_label_hash_or_rank_used"],
        "edge_toggle_variational_response_coefficients_computed": audit["edge_toggle_variational_response_coefficients_computed"],
        "dimensionless_normalization_available_for_candidate": audit["dimensionless_normalization_available_for_candidate"],
        "candidate_separates_full_carrier": audit["candidate_class_count"] == EXPECTED_GRAPH_COUNT,
        "proved_variational_derivative_into_K_or_Ltotal_exported": False,
        "typed_graph_to_K_or_Ltotal_coupling_theorem_exported": False,
        "selector_bridge_or_role_transfer_imported": False,
    }
    accepted = all([
        facts["p2818_local_response_rejected"],
        facts["p2832_low_order_formula_rejected"],
        facts["exactly_one_higher_order_non_label_formula_tested"],
        not facts["digest_label_hash_or_rank_used"],
        facts["edge_toggle_variational_response_coefficients_computed"],
        facts["dimensionless_normalization_available_for_candidate"],
        facts["candidate_separates_full_carrier"],
        facts["proved_variational_derivative_into_K_or_Ltotal_exported"],
        facts["typed_graph_to_K_or_Ltotal_coupling_theorem_exported"],
        not facts["selector_bridge_or_role_transfer_imported"],
    ])
    return {
        "facts": facts,
        "accepted_as_source_law_coupling": accepted,
        "accepted_as_bounded_edge_toggle_witness_with_residual_no_go": not accepted,
        "missing_for_promotion": [
            key for key, value in facts.items()
            if (key in {"digest_label_hash_or_rank_used", "selector_bridge_or_role_transfer_imported"} and value)
            or (key not in {"digest_label_hash_or_rank_used", "selector_bridge_or_role_transfer_imported"} and not value)
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["edge_toggle_response_polynomial_source_audit"]
    acceptance = payload["acceptance_matrix"]
    lines = [
        "# P2833/S1783 edge-toggle response polynomial source audit", "", f"Status: `{payload['status']}`", "",
        "## Candidate formula", audit["candidate_formula"], "", "## Finite counts",
        f"- decoded_graph_count={audit['decoded_graph_count']}",
        f"- candidate_class_count={audit['candidate_class_count']}",
        f"- candidate_collision_class_count={audit['candidate_collision_class_count']}",
        f"- candidate_collision_graph_count={audit['candidate_collision_graph_count']}",
        f"- candidate_max_class_size={audit['candidate_max_class_size']}",
        f"- candidate_defect_after_formula={audit['candidate_defect_after_formula']}", "",
        "## Acceptance",
        f"- accepted_as_source_law_coupling={acceptance['accepted_as_source_law_coupling']}",
        f"- accepted_as_bounded_edge_toggle_witness_with_residual_no_go={acceptance['accepted_as_bounded_edge_toggle_witness_with_residual_no_go']}", "",
        "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2818 = read_json(P2818)
    p2832 = read_json(P2832)
    audit = build_audit(p2832)
    payload: dict[str, Any] = {
        "status": "P2833_EDGE_TOGGLE_RESPONSE_POLYNOMIAL_RESIDUAL_NO_GO_NO_COUPLING_NO_CLOSURE",
        "input_hashes": {"P2818": sha(P2818), "P2832": sha(P2832), "16_4_4.scd": sha(SCD)},
        "input_statuses": {"P2818": p2818.get("status"), "P2832": p2832.get("status")},
        "edge_toggle_response_polynomial_source_audit": audit,
        "decision": {
            "negative_export_flags": {
                "strict_graph_source_law_exported": False,
                "typed_coupling_to_K_or_Ltotal_exported": False,
                "proved_variational_derivative_exported": False,
                "role_bearing_ltotal_promoted": False,
                "bridge_closure_exported": False,
                "role_transfer_started": False,
                "selector_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2833 tests one higher-order non-label formula requested after P2832: the exact edge-toggle response polynomial carrying triangle and 4-cycle response coefficients for every unordered vertex pair.  It is a real finite variational-response witness and is much sharper than Q_cyc, yielding 16,757 classes on the 16,828-graph carrier.  However, 67 residual collision classes covering 138 graphs remain, so it still does not separate the full carrier; moreover no proved variational derivative into K/L_total and no typed graph-to-K/L_total coupling theorem is exported.  Source-law/coupling promotion remains rejected.",
            "next_honest_step": "Do not replay coarser local motifs or promote the edge-toggle polynomial to L_total.  The next honest proof-grade move is a residual-pair obstruction/refinement audit restricted to the 67 residual P2833 collision classes: compute one explicit non-label second-variation or two-edge-toggle interaction functional and stop on first unresolved collision or missing typed coupling premise.  If that residual audit still leaves collisions or lacks a coupling theorem, preserve the P2831-P2833 no-coupling boundary rather than manufacturing closure.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit, p2818, p2832)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2833/S1783 edge-toggle response polynomial audit", "## P2833/S1783 edge-toggle response polynomial audit\n\n`P2833/S1783` tests a higher-order non-label finite variational-response candidate: `Q_toggle(G)=multiset_{i<j}(1_{ij in E}, |N(i)∩N(j)|, #{length-3 paths i-a-b-j})`, the exact triangle/4-cycle edge-toggle response polynomial on the fixed `16`-node carrier.  It yields `16,757` classes on `16,828` graphs, leaving `67` residual collision classes covering `138` graphs with max class size `3`.  Because it still does not separate the full carrier and exports no proved variational derivative into `K`/`L_total` or typed graph-to-`K`/`L_total` coupling theorem, no strict graph-source law, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2833/S1783 edge-toggle response Ltotal guard", "## P2833/S1783 edge-toggle response Ltotal guard\n\n`P2833/S1783` adds no term to `L_total`.  The edge-toggle response polynomial is a sharper non-label finite variational-response witness than P2832, but it leaves `138` residual collision graphs in `67` classes and lacks an Euler-Lagrange derivative/coupling theorem into `K` or `L_total`.\n")
    append_once(AGENTS, "Current edge-toggle response polynomial guardrail (P2833/S1783, 2026-06-17)", "## Current edge-toggle response polynomial guardrail (P2833/S1783, 2026-06-17)\n\n- P2833 tests one higher-order non-label candidate after P2832: the exact edge-toggle response polynomial `multiset_{i<j}(1_{ij in E}, |N(i)∩N(j)|, #{length-3 paths i-a-b-j})` on the full `16,828`-graph carrier.\n- The polynomial is a real finite variational-response witness and refines the carrier to `16,757` classes, but it leaves `67` residual collision classes covering `138` graphs; it also exports no proved variational derivative into `K`/`L_total` and no typed graph-to-`K`/`L_total` coupling theorem.\n- Do not promote P2833 to a strict graph-source law, `L_total`, bridge closure, role transfer, selector closure, or ToE closure.  A next admissible move in this lane must be a residual-only non-label second-variation/two-edge-toggle refinement audit with stop-on-first-residual and no selector/bridge/role-transfer import.\n")
    return payload


if __name__ == "__main__":
    main()
