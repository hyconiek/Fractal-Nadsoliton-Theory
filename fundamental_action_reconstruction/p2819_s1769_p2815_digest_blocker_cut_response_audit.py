#!/usr/bin/env python3
"""P2819/S1769: bounded P2815-digest response audit on a deterministic 24-graph sample of one P2818 blocker cut.

P2818 showed that degree-level local edge responses remain too coarse.  P2819
runs the next honest, computationally bounded step: apply the full P2815-style
spectral/local-motif edge-toggle response digest to exactly one blocker cut, the
largest P2818 collision class.  This tests whether the recommended stronger
response ingredient has real separating power without pretending to audit the
entire 16,828-graph carrier or to export a K/L_total coupling theorem.
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import SCD, decode_scd_bytes, sha
from p2811_s1761_local_motif_refined_source_candidate_audit import stable_digest
from p2812_s1762_two_wl_refined_collision_audit import p2811_refined_key, read_json
from p2815_s1765_edge_toggle_response_residual_audit import edge_set_from_graph, graph_from_edge_set
from p2818_s1768_local_edge_variational_response_energy_audit import local_edge_variational_response_energy

GEN = ROOT / "generated"
P2818 = GEN / "p2818_s1768_local_edge_variational_response_energy_audit.json"
OUT = GEN / "p2819_s1769_p2815_digest_blocker_cut_response_audit.json"
MD = GEN / "p2819_s1769_p2815_digest_blocker_cut_response_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def p2815_digest_response(graph: tuple[tuple[int, ...], ...]) -> dict[str, Any]:
    edges = edge_set_from_graph(graph)
    all_pairs = {(left, right) for left in range(16) for right in range(left + 1, 16)}
    nonedges = all_pairs - edges
    deletion_digests = []
    for edge in sorted(edges):
        modified = set(edges)
        modified.remove(edge)
        deletion_digests.append(p2811_refined_key(graph_from_edge_set(modified)))
    addition_digests = []
    for nonedge in sorted(nonedges):
        modified = set(edges)
        modified.add(nonedge)
        addition_digests.append(p2811_refined_key(graph_from_edge_set(modified)))
    return {
        "edge_count": len(edges),
        "nonedge_count": len(nonedges),
        "deletion_digest_multiset": tuple(sorted(Counter(deletion_digests).items())),
        "addition_digest_multiset": tuple(sorted(Counter(addition_digests).items())),
    }


def largest_p2818_collision_class(graphs: list[tuple[tuple[int, ...], ...]]) -> tuple[Any, list[int]]:
    classes: dict[tuple[Any, ...], list[int]] = defaultdict(list)
    for index, graph in enumerate(graphs):
        classes[local_edge_variational_response_energy(graph)].append(index)
    collisions = {key: indices for key, indices in classes.items() if len(indices) > 1}
    return max(collisions.items(), key=lambda item: (len(item[1]), [-i for i in item[1]]))


def build_audit() -> dict[str, Any]:
    graphs, _ = decode_scd_bytes(SCD.read_bytes())
    blocker_key, blocker_indices = largest_p2818_collision_class(graphs)
    sampled_indices = blocker_indices[:24]
    refined: dict[str, list[int]] = defaultdict(list)
    sample_rows = []
    for index in sampled_indices:
        response = p2815_digest_response(graphs[index])
        digest = stable_digest(response)
        refined[digest].append(index)
        if len(sample_rows) < 12:
            sample_rows.append({
                "graph_index_0_based": index,
                "p2815_digest_response_sha256": digest,
                "edge_count": response["edge_count"],
                "nonedge_count": response["nonedge_count"],
                "distinct_deletion_digest_count": len(response["deletion_digest_multiset"]),
                "distinct_addition_digest_count": len(response["addition_digest_multiset"]),
            })
    residual = {key: indices for key, indices in refined.items() if len(indices) > 1}
    class_size_histogram = Counter(len(indices) for indices in refined.values())
    return {
        "candidate_response": "full P2815-style spectral/local-motif edge-toggle response digest on a deterministic 24-graph sample of one P2818 blocker cut",
        "blocker_cut": "first 24 graph indices from the largest P2818 E_local_toggle collision class",
        "decoded_graph_count": len(graphs),
        "p2818_blocker_class_size": len(blocker_indices),
        "p2818_blocker_key_sha256": stable_digest(blocker_key),
        "computed_graph_count": len(sampled_indices),
        "toggle_digest_refined_class_count": len(refined),
        "toggle_digest_collision_class_count": len(residual),
        "toggle_digest_collision_graph_count": sum(len(indices) for indices in residual.values()),
        "sample_fraction_of_blocker_cut": [len(sampled_indices), len(blocker_indices)],
        "toggle_digest_max_class_size": max(class_size_histogram),
        "sample_defect_after_toggle_digest": len(sampled_indices) - len(refined),
        "class_size_histogram": dict(sorted(class_size_histogram.items())),
        "sample_rows": sample_rows,
    }


def acceptance_matrix(audit: dict[str, Any], p2818: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2818_local_response_rejected": p2818["acceptance_matrix"]["accepted_as_bounded_candidate_rejection"],
        "exactly_one_blocker_cut_tested": True,
        "uses_full_p2815_spectral_local_motif_response_digest": True,
        "blocker_cut_is_fully_separated": audit["toggle_digest_collision_class_count"] == 0,
        "all_16828_graphs_audited": audit["computed_graph_count"] == audit["decoded_graph_count"],
        "entire_blocker_cut_audited": audit["computed_graph_count"] == audit["p2818_blocker_class_size"],
        "strict_graph_source_law_exported": False,
        "typed_graph_to_K_or_Ltotal_coupling_exported": False,
        "selector_bridge_or_role_transfer_imported": False,
    }
    accepted = all(facts[key] for key in [
        "p2818_local_response_rejected",
        "exactly_one_blocker_cut_tested",
        "uses_full_p2815_spectral_local_motif_response_digest",
        "blocker_cut_is_fully_separated",
        "all_16828_graphs_audited",
        "entire_blocker_cut_audited",
        "strict_graph_source_law_exported",
        "typed_graph_to_K_or_Ltotal_coupling_exported",
    ]) and not facts["selector_bridge_or_role_transfer_imported"]
    return {
        "facts": facts,
        "accepted_as_sample_response_witness": facts["blocker_cut_is_fully_separated"],
        "accepted_as_full_carrier_source_law_coupling": accepted,
        "accepted_as_bounded_no_closure_audit": not accepted,
        "missing_for_promotion": [key for key, value in facts.items() if not value and key != "selector_bridge_or_role_transfer_imported"],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["p2815_digest_blocker_cut_audit"]
    acceptance = payload["acceptance_matrix"]
    lines = [
        "# P2819/S1769 P2815-digest blocker-cut sample response audit", "", f"Status: `{payload['status']}`", "",
        "## Candidate", audit["candidate_response"], "", "## Blocker cut", audit["blocker_cut"], "", "## Finite counts",
        f"- decoded_graph_count={audit['decoded_graph_count']}",
        f"- p2818_blocker_class_size={audit['p2818_blocker_class_size']}",
        f"- computed_graph_count={audit['computed_graph_count']}",
        f"- toggle_digest_refined_class_count={audit['toggle_digest_refined_class_count']}",
        f"- toggle_digest_collision_class_count={audit['toggle_digest_collision_class_count']}",
        f"- toggle_digest_collision_graph_count={audit['toggle_digest_collision_graph_count']}",
        f"- toggle_digest_max_class_size={audit['toggle_digest_max_class_size']}",
        f"- sample_defect_after_toggle_digest={audit['sample_defect_after_toggle_digest']}", "",
        "## Acceptance",
        f"- accepted_as_sample_response_witness={acceptance['accepted_as_sample_response_witness']}",
        f"- accepted_as_full_carrier_source_law_coupling={acceptance['accepted_as_full_carrier_source_law_coupling']}",
        f"- accepted_as_bounded_no_closure_audit={acceptance['accepted_as_bounded_no_closure_audit']}", "",
        "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2818 = read_json(P2818)
    audit = build_audit()
    blocker_separated = audit["toggle_digest_collision_class_count"] == 0
    payload: dict[str, Any] = {
        "status": "P2819_P2815_DIGEST_BLOCKER_CUT_SAMPLE_WITNESS_NO_FULL_COUPLING_NO_CLOSURE" if blocker_separated else "P2819_P2815_DIGEST_BLOCKER_CUT_RESIDUAL_NO_CLOSURE",
        "input_hashes": {"P2818": sha(P2818), "16_4_4.scd": sha(SCD)},
        "input_statuses": {"P2818": p2818.get("status")},
        "p2815_digest_blocker_cut_audit": audit,
        "decision": {
            "negative_export_flags": {
                "full_16828_carrier_audited": False,
                "strict_graph_source_law_exported": False,
                "typed_coupling_to_K_or_Ltotal_exported": False,
                "role_bearing_ltotal_promoted": False,
                "bridge_closure_exported": False,
                "role_transfer_started": False,
                "selector_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2819 applies the full P2815-style spectral/local-motif edge-toggle response digest to a deterministic 24-graph prefix sample of exactly one bounded blocker cut: the largest P2818 collision class.  This is a real computational witness about the recommended stronger ingredient, but it audits only a 24-graph sample of one blocker cut rather than the whole blocker cut or all 16,828 graphs and exports no strict graph-source law or typed K/L_total coupling theorem.  Therefore no L_total, bridge, role-transfer, selector, or ToE closure follows.",
            "next_honest_step": "If the sample witness is positive, extend the same P2815-digest response audit first to the full largest blocker cut and then to all P2818 collision classes with caching and a compact manifest; otherwise refine the response formula once.  In either case, keep source-law/coupling promotion blocked until a full-carrier audit plus typed graph-to-K/L_total coupling theorem with units and variational derivative is exported.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit, p2818)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2819/S1769 P2815 digest blocker-cut sample response audit", "## P2819/S1769 P2815 digest blocker-cut sample response audit\n\n`P2819/S1769` applies the full P2815-style spectral/local-motif edge-toggle response digest to a deterministic 24-graph prefix sample of exactly one bounded blocker cut: the largest P2818 local-response collision class.  This tests the stronger recommended ingredient on a finite obstruction slice, but it is not a full `16,828`-graph source-law/coupling audit and exports no strict graph-source law or typed `K`/`L_total` coupling theorem.  No `L_total`, bridge, role-transfer, selector, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2819/S1769 P2815 digest blocker-cut sample Ltotal guard", "## P2819/S1769 P2815 digest blocker-cut sample Ltotal guard\n\n`P2819/S1769` adds no variational term to `L_total`.  The P2815-style spectral/local-motif edge-toggle response digest is tested only on a deterministic 24-graph sample of one P2818 blocker cut; without a full-carrier audit and typed graph-to-`K`/`L_total` coupling theorem, it remains finite graph evidence rather than an Euler-Lagrange/source-law term.\n")
    append_once(AGENTS, "Current P2815-digest blocker-cut sample response guardrail (P2819/S1769, 2026-06-17)", "## Current P2815-digest blocker-cut sample response guardrail (P2819/S1769, 2026-06-17)\n\n- P2819 applies the full P2815-style spectral/local-motif edge-toggle response digest to a deterministic 24-graph prefix sample of exactly one bounded blocker cut: the largest P2818 local-response collision class.\n- This is finite blocker-cut sample evidence only; it is not a full `16,828`-graph source-law/coupling audit and exports no strict graph-source law, variational meaning, typed `K`/`L_total` coupling theorem, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure.\n- A next admissible move may extend the same digest response audit first to the full largest blocker cut and then to all P2818 collision classes with caching and a compact manifest, or supply an analytic graph-to-`K`/`L_total` coupling theorem with units and variational derivative; do not promote closure before that.\n")
    return payload


if __name__ == "__main__":
    main()
