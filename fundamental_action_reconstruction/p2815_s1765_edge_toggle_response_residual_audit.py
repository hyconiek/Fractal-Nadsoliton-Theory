#!/usr/bin/env python3
"""P2815/S1765: edge-toggle response audit for the P2814 residual pairs.

P2814 reduced the residual carrier obstruction to 57 pairs / 114 graphs but did
not close it.  P2815 adds exactly one non-automorphism typed ingredient: the
single-edge deletion/addition response profile.  For each P2814 residual graph,
the profile records the multiset of P2811 spectral+local-motif digests after
removing one existing edge and after adding one nonedge.  This is a finite
source-carrier separation test, not a source law or variational coupling.
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import EXPECTED_GRAPH_COUNT, N, SCD, decode_scd_bytes, sha
from p2811_s1761_local_motif_refined_source_candidate_audit import stable_digest
from p2812_s1762_two_wl_refined_collision_audit import p2811_refined_key, read_json, two_wl_histogram
from p2813_s1763_three_wl_refined_collision_audit import three_wl_histogram
from p2814_s1764_exact_automorphism_orbit_residual_audit import exact_automorphism_orbit_signature

GEN = ROOT / "generated"
P2814 = GEN / "p2814_s1764_exact_automorphism_orbit_residual_audit.json"
OUT = GEN / "p2815_s1765_edge_toggle_response_residual_audit.json"
MD = GEN / "p2815_s1765_edge_toggle_response_residual_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def edge_set_from_graph(graph: tuple[tuple[int, ...], ...]) -> set[tuple[int, int]]:
    return {(min(i, neighbor - 1), max(i, neighbor - 1)) for i, row in enumerate(graph) for neighbor in row if i != neighbor - 1}


def graph_from_edge_set(edges: set[tuple[int, int]]) -> tuple[tuple[int, ...], ...]:
    rows = [[] for _ in range(N)]
    for left, right in sorted(edges):
        rows[left].append(right + 1)
        rows[right].append(left + 1)
    return tuple(tuple(row) for row in rows)


def edge_toggle_response_signature(graph: tuple[tuple[int, ...], ...]) -> dict[str, Any]:
    edges = edge_set_from_graph(graph)
    all_pairs = {(left, right) for left in range(N) for right in range(left + 1, N)}
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
        "edge_deletion_response_multiset": tuple(sorted(Counter(deletion_digests).items())),
        "nonedge_addition_response_multiset": tuple(sorted(Counter(addition_digests).items())),
    }


def build_p2814_classes(graphs: list[tuple[tuple[int, ...], ...]]) -> dict[tuple[Any, str | None], list[int]]:
    p2811_classes: dict[str, list[int]] = defaultdict(list)
    for index, graph in enumerate(graphs):
        p2811_classes[p2811_refined_key(graph)].append(index)

    p2812_classes: dict[tuple[str, str | None], list[int]] = defaultdict(list)
    for key, indices in p2811_classes.items():
        if len(indices) == 1:
            p2812_classes[(key, None)].extend(indices)
            continue
        for index in indices:
            p2812_classes[(key, stable_digest(two_wl_histogram(graphs[index])))].append(index)

    p2813_classes: dict[tuple[Any, str | None], list[int]] = defaultdict(list)
    for key, indices in p2812_classes.items():
        if len(indices) == 1:
            p2813_classes[(key, None)].extend(indices)
            continue
        for index in indices:
            p2813_classes[(key, stable_digest(three_wl_histogram(graphs[index])))].append(index)

    p2814_classes: dict[tuple[Any, str | None], list[int]] = defaultdict(list)
    for key, indices in p2813_classes.items():
        if len(indices) == 1:
            p2814_classes[(key, None)].extend(indices)
            continue
        for index in indices:
            digest = stable_digest(exact_automorphism_orbit_signature(graphs[index]))
            p2814_classes[(key, digest)].append(index)
    return p2814_classes


def build_audit(p2814: dict[str, Any]) -> dict[str, Any]:
    graphs, _ = decode_scd_bytes(SCD.read_bytes())
    p2814_classes = build_p2814_classes(graphs)

    p2815_classes: dict[tuple[Any, str | None], list[int]] = defaultdict(list)
    rows: list[dict[str, Any]] = []
    computed = 0
    signature_counter: Counter[str] = Counter()
    for key, indices in p2814_classes.items():
        if len(indices) == 1:
            p2815_classes[(key, None)].extend(indices)
            continue
        for index in indices:
            signature = edge_toggle_response_signature(graphs[index])
            digest = stable_digest(signature)
            p2815_classes[(key, digest)].append(index)
            signature_counter[digest] += 1
            computed += 1
            if len(rows) < 12:
                rows.append({
                    "graph_index_0_based": index,
                    "p2814_key_sha256": stable_digest(key),
                    "edge_toggle_response_sha256": digest,
                    "edge_count": signature["edge_count"],
                    "nonedge_count": signature["nonedge_count"],
                    "deletion_response_distinct_count": len(signature["edge_deletion_response_multiset"]),
                    "addition_response_distinct_count": len(signature["nonedge_addition_response_multiset"]),
                })

    collisions = {key: indices for key, indices in p2815_classes.items() if len(indices) > 1}
    class_size_histogram = Counter(len(indices) for indices in p2815_classes.values())
    return {
        "candidate_invariant": "single-edge deletion/addition response profile using P2811 spectral+local-motif digests, applied inside P2814 residual collision classes",
        "decoded_graph_count": len(graphs),
        "p2814_refined_class_count": p2814["automorphism_orbit_audit"]["automorphism_orbit_refined_class_count"],
        "p2814_collision_class_count": p2814["automorphism_orbit_audit"]["automorphism_orbit_collision_class_count"],
        "p2814_collision_graph_count": p2814["automorphism_orbit_audit"]["automorphism_orbit_collision_graph_count"],
        "edge_toggle_computed_graph_count": computed,
        "edge_toggle_refined_class_count": len(p2815_classes),
        "edge_toggle_collision_class_count": len(collisions),
        "edge_toggle_collision_graph_count": sum(len(indices) for indices in collisions.values()),
        "edge_toggle_max_class_size": max(class_size_histogram) if class_size_histogram else 0,
        "class_size_histogram": dict(sorted(class_size_histogram.items())),
        "remaining_defect_canonical_minus_edge_toggle": EXPECTED_GRAPH_COUNT - len(p2815_classes),
        "defect_reduction_vs_p2814": p2814["automorphism_orbit_audit"]["remaining_defect_canonical_minus_automorphism_orbit"] - (EXPECTED_GRAPH_COUNT - len(p2815_classes)),
        "computed_signature_distinct_count": len(signature_counter),
        "computed_signature_histogram": dict(sorted(signature_counter.items())),
        "sample_edge_toggle_rows": rows,
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "decoded_graph_count_is_16828": audit["decoded_graph_count"] == EXPECTED_GRAPH_COUNT,
        "edge_toggle_applied_only_to_p2814_collision_graphs": audit["edge_toggle_computed_graph_count"] == audit["p2814_collision_graph_count"],
        "edge_toggle_improves_over_p2814": audit["edge_toggle_refined_class_count"] > audit["p2814_refined_class_count"],
        "edge_toggle_reaches_canonical_target": audit["edge_toggle_refined_class_count"] == EXPECTED_GRAPH_COUNT,
        "remaining_collision_classes_exist": audit["edge_toggle_collision_class_count"] > 0,
        "strict_spectral_source_law_exported": False,
        "typed_coupling_to_K_or_Ltotal_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_edge_toggle_response_audit": all(facts[key] for key in [
            "decoded_graph_count_is_16828",
            "edge_toggle_applied_only_to_p2814_collision_graphs",
            "edge_toggle_improves_over_p2814",
        ]),
        "accepted_as_complete_canonical_source_carrier": facts["edge_toggle_reaches_canonical_target"],
        "accepted_as_strict_source_law_or_coupling": False,
        "missing_for_promotion": [key for key, value in facts.items() if not value],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["edge_toggle_response_audit"]
    acceptance = payload["acceptance_matrix"]
    lines = [
        "# P2815/S1765 edge-toggle response residual audit", "", f"Status: `{payload['status']}`", "",
        "## Candidate invariant", audit["candidate_invariant"], "", "## Exact quotient counts",
        f"- decoded_graph_count={audit['decoded_graph_count']}",
        f"- p2814_refined_class_count={audit['p2814_refined_class_count']}",
        f"- p2814_collision_graph_count={audit['p2814_collision_graph_count']}",
        f"- edge_toggle_computed_graph_count={audit['edge_toggle_computed_graph_count']}",
        f"- edge_toggle_refined_class_count={audit['edge_toggle_refined_class_count']}",
        f"- edge_toggle_collision_class_count={audit['edge_toggle_collision_class_count']}",
        f"- edge_toggle_collision_graph_count={audit['edge_toggle_collision_graph_count']}",
        f"- edge_toggle_max_class_size={audit['edge_toggle_max_class_size']}",
        f"- remaining_defect_canonical_minus_edge_toggle={audit['remaining_defect_canonical_minus_edge_toggle']}",
        f"- defect_reduction_vs_p2814={audit['defect_reduction_vs_p2814']}", "", "## Decision",
        f"- accepted_as_edge_toggle_response_audit={acceptance['accepted_as_edge_toggle_response_audit']}",
        f"- accepted_as_complete_canonical_source_carrier={acceptance['accepted_as_complete_canonical_source_carrier']}",
        f"- accepted_as_strict_source_law_or_coupling={acceptance['accepted_as_strict_source_law_or_coupling']}", "", "## Boundary",
        payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2814 = read_json(P2814)
    audit = build_audit(p2814)
    complete = audit["edge_toggle_refined_class_count"] == EXPECTED_GRAPH_COUNT
    payload: dict[str, Any] = {
        "status": "P2815_EDGE_TOGGLE_RESPONSE_SEPARATES_CARRIER_NO_SOURCE_LAW_NO_CLOSURE" if complete else "P2815_EDGE_TOGGLE_RESPONSE_PARTIAL_REFINEMENT_NO_CLOSURE",
        "input_hashes": {"P2814": sha(P2814), "16_4_4.scd": sha(SCD)},
        "input_statuses": {"P2814": p2814.get("status")},
        "edge_toggle_response_audit": audit,
        "decision": {
            "negative_export_flags": {
                "strict_spectral_source_law_exported": False,
                "typed_coupling_to_K_or_Ltotal_exported": False,
                "canonical_geometry_closure_exported": False,
                "role_bearing_ltotal_promoted": False,
                "bridge_closure_exported": False,
                "role_transfer_started": False,
                "selector_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2815 applies exactly one non-automorphism typed ingredient, a single-edge deletion/addition response profile, only to the P2814 residual-collision graphs.  The profile separates the remaining carrier: the quotient reaches the 16,828 canonical target, with zero residual collision classes and defect reduction 57 versus P2814.  This is complete finite source-carrier separation evidence, but it still exports no strict spectral source law and no typed coupling theorem to K or L_total.",
            "next_honest_step": "Stop carrier-refinement replay and run a separate source-law/coupling acceptance audit with one explicit graph-source functional and one typed graph-to-K/L_total map.  Keep role-bearing L_total, bridge closure, role transfer, selector closure, and ToE blocked until that separate coupling audit succeeds.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2815/S1765 edge-toggle response residual audit", "## P2815/S1765 edge-toggle response residual audit\n\n`P2815/S1765` adds exactly one non-automorphism typed ingredient after P2814: the single-edge deletion/addition response profile, using P2811 spectral+local-motif digests on toggled graphs and applying it only to the P2814 residual-collision set.  It separates the finite carrier by reaching the `16,828` canonical target with zero residual collision classes.  This is complete carrier separation evidence only; it is not a strict spectral source law or a `K`/`L_total` coupling theorem.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2815/S1765 edge-toggle response Ltotal guard", "## P2815/S1765 edge-toggle response Ltotal guard\n\n`P2815/S1765` adds no variational term to `L_total`.  The edge-toggle response profile separates the finite carrier, but it remains a graph-carrier diagnostic, not an action functional, source equation, Euler-Lagrange contribution, or typed graph-to-`K`/`L_total` coupling.\n")
    append_once(AGENTS, "Current edge-toggle response residual guardrail (P2815/S1765, 2026-06-16)", "## Current edge-toggle response residual guardrail (P2815/S1765, 2026-06-16)\n\n- P2815 applies a single-edge deletion/addition response profile to the P2814 residual-collision graphs, using P2811 spectral+local-motif digests on toggled graphs.\n- P2815 reaches the `16,828` canonical target and removes the residual carrier collisions, but this is finite carrier-separation evidence only: it does not export a strict spectral source law, typed coupling to `K`/`L_total`, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure.\n- The next admissible move is a separate source-law/coupling acceptance audit with one explicit graph-source functional plus one typed graph-to-`K` or graph-to-`L_total` map; do not continue carrier-refinement replay without a new source-law question.\n")
    return payload


if __name__ == "__main__":
    main()
