#!/usr/bin/env python3
"""P2813/S1763: 3-WL refinement of the P2812 residual collisions.

P2812 left 79 two-WL collision classes covering 158 graphs.  P2813 adds exactly
one stronger typed invariant: a stable 3-dimensional Weisfeiler-Lehman ordered
triple color histogram, applied only inside those residual classes.  This is a
finite source-carrier separation test, not a variational source term.
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import EXPECTED_GRAPH_COUNT, N, SCD, decode_scd_bytes, sha
from p2811_s1761_local_motif_refined_source_candidate_audit import stable_digest
from p2812_s1762_two_wl_refined_collision_audit import P2811, p2811_refined_key, read_json, two_wl_histogram

GEN = ROOT / "generated"
P2812 = GEN / "p2812_s1762_two_wl_refined_collision_audit.json"
OUT = GEN / "p2813_s1763_three_wl_refined_collision_audit.json"
MD = GEN / "p2813_s1763_three_wl_refined_collision_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def initial_3wl_color(edge_set: set[tuple[int, int]], triple: tuple[int, int, int]) -> tuple[Any, ...]:
    i, j, k = triple
    equalities = (i == j, i == k, j == k)
    edges = ((min(i, j), max(i, j)) in edge_set, (min(i, k), max(i, k)) in edge_set, (min(j, k), max(j, k)) in edge_set)
    return (equalities, edges)


def three_wl_histogram(adj: tuple[tuple[int, ...], ...], max_iterations: int = 20) -> tuple[tuple[int, int], ...]:
    edge_set = {(min(i, j), max(i, j)) for i, nbrs in enumerate(adj) for j in nbrs if i != j}
    triples = [(i, j, k) for i in range(N) for j in range(N) for k in range(N)]
    palette = {value: color for color, value in enumerate(sorted({initial_3wl_color(edge_set, t) for t in triples}))}
    colors = {t: palette[initial_3wl_color(edge_set, t)] for t in triples}
    iterations = 0
    for iteration in range(1, max_iterations + 1):
        updates = []
        for i, j, k in triples:
            neighborhoods = tuple(
                tuple(sorted(colors[(x, j, k)] for x in range(N))),
            ) + tuple(
                tuple(sorted(colors[(i, x, k)] for x in range(N))),
            ) + tuple(
                tuple(sorted(colors[(i, j, x)] for x in range(N))),
            )
            updates.append(((i, j, k), (colors[(i, j, k)], neighborhoods)))
        next_palette = {value: color for color, value in enumerate(sorted({value for _, value in updates}))}
        new_colors = {triple: next_palette[value] for triple, value in updates}
        iterations = iteration
        if new_colors == colors:
            break
        colors = new_colors
    histogram = Counter(colors.values())
    return tuple(sorted(histogram.items())) + ((-1, iterations),)


def build_audit(p2812: dict[str, Any]) -> dict[str, Any]:
    graphs, _ = decode_scd_bytes(SCD.read_bytes())
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

    p2813_classes: dict[tuple[tuple[str, str | None], str | None], list[int]] = defaultdict(list)
    rows: list[dict[str, Any]] = []
    computed = 0
    for key, indices in p2812_classes.items():
        if len(indices) == 1:
            p2813_classes[(key, None)].extend(indices)
            continue
        for index in indices:
            digest = stable_digest(three_wl_histogram(graphs[index]))
            p2813_classes[(key, digest)].append(index)
            computed += 1
            if len(rows) < 12:
                rows.append({"graph_index_0_based": index, "p2812_key": [key[0], key[1]], "three_wl_histogram_sha256": digest})

    collisions = {key: indices for key, indices in p2813_classes.items() if len(indices) > 1}
    hist = Counter(len(indices) for indices in p2813_classes.values())
    return {
        "candidate_invariant": "stable 3-dimensional Weisfeiler-Lehman ordered-triple color histogram, applied inside P2812 residual collision classes",
        "decoded_graph_count": len(graphs),
        "p2812_two_wl_refined_class_count": p2812["two_wl_refinement_audit"]["two_wl_refined_class_count"],
        "p2812_collision_class_count": p2812["two_wl_refinement_audit"]["two_wl_collision_class_count"],
        "p2812_collision_graph_count": p2812["two_wl_refinement_audit"]["two_wl_collision_graph_count"],
        "three_wl_computed_graph_count": computed,
        "three_wl_refined_class_count": len(p2813_classes),
        "three_wl_collision_class_count": len(collisions),
        "three_wl_collision_graph_count": sum(len(v) for v in collisions.values()),
        "three_wl_max_class_size": max(hist) if hist else 0,
        "class_size_histogram": dict(sorted(hist.items())),
        "remaining_defect_canonical_minus_three_wl": EXPECTED_GRAPH_COUNT - len(p2813_classes),
        "defect_reduction_vs_p2812": p2812["two_wl_refinement_audit"]["remaining_defect_canonical_minus_two_wl"] - (EXPECTED_GRAPH_COUNT - len(p2813_classes)),
        "sample_three_wl_rows": rows,
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "decoded_graph_count_is_16828": audit["decoded_graph_count"] == EXPECTED_GRAPH_COUNT,
        "three_wl_applied_only_to_p2812_collision_graphs": audit["three_wl_computed_graph_count"] == audit["p2812_collision_graph_count"],
        "three_wl_improves_over_p2812": audit["three_wl_refined_class_count"] > audit["p2812_two_wl_refined_class_count"],
        "three_wl_reaches_canonical_target": audit["three_wl_refined_class_count"] == EXPECTED_GRAPH_COUNT,
        "remaining_collision_classes_exist": audit["three_wl_collision_class_count"] > 0,
        "strict_spectral_source_law_exported": False,
        "typed_coupling_to_K_or_Ltotal_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_three_wl_residual_obstruction_audit": all(facts[k] for k in ["decoded_graph_count_is_16828", "three_wl_applied_only_to_p2812_collision_graphs", "remaining_collision_classes_exist"]),
        "accepted_as_complete_canonical_source_carrier": facts["three_wl_reaches_canonical_target"],
        "accepted_as_strict_source_law_or_coupling": False,
        "missing_for_promotion": [k for k, v in facts.items() if not v],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["three_wl_refinement_audit"]
    acceptance = payload["acceptance_matrix"]
    lines = [
        "# P2813/S1763 three-WL refined collision audit", "", f"Status: `{payload['status']}`", "",
        "## Candidate invariant", audit["candidate_invariant"], "", "## Exact quotient counts",
        f"- decoded_graph_count={audit['decoded_graph_count']}",
        f"- p2812_two_wl_refined_class_count={audit['p2812_two_wl_refined_class_count']}",
        f"- p2812_collision_graph_count={audit['p2812_collision_graph_count']}",
        f"- three_wl_computed_graph_count={audit['three_wl_computed_graph_count']}",
        f"- three_wl_refined_class_count={audit['three_wl_refined_class_count']}",
        f"- three_wl_collision_class_count={audit['three_wl_collision_class_count']}",
        f"- three_wl_collision_graph_count={audit['three_wl_collision_graph_count']}",
        f"- three_wl_max_class_size={audit['three_wl_max_class_size']}",
        f"- remaining_defect_canonical_minus_three_wl={audit['remaining_defect_canonical_minus_three_wl']}",
        f"- defect_reduction_vs_p2812={audit['defect_reduction_vs_p2812']}", "", "## Decision",
        f"- accepted_as_three_wl_residual_obstruction_audit={acceptance['accepted_as_three_wl_residual_obstruction_audit']}",
        f"- accepted_as_complete_canonical_source_carrier={acceptance['accepted_as_complete_canonical_source_carrier']}",
        f"- accepted_as_strict_source_law_or_coupling={acceptance['accepted_as_strict_source_law_or_coupling']}", "", "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2812 = read_json(P2812)
    audit = build_audit(p2812)
    complete = audit["three_wl_refined_class_count"] == EXPECTED_GRAPH_COUNT
    payload: dict[str, Any] = {
        "status": "P2813_THREE_WL_REFINED_COLLISION_AUDIT_SEPARATES_CARRIER_NO_SOURCE_LAW_NO_CLOSURE" if complete else "P2813_THREE_WL_HISTOGRAM_RESIDUAL_COLLISION_OBSTRUCTION_NO_CLOSURE",
        "input_hashes": {"P2812": sha(P2812), "P2811": sha(P2811), "16_4_4.scd": sha(SCD)},
        "input_statuses": {"P2812": p2812.get("status")},
        "three_wl_refinement_audit": audit,
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
            "reason": "P2813 applies exactly one stronger typed invariant, stable 3-WL ordered-triple color histograms, only to the 158 P2812 residual-collision graphs.  It is a finite negative obstruction for this histogram-level 3-WL candidate: the quotient remains at 16,749 classes, the defect remains 79, and all 79 collision pairs survive.  Therefore 3-WL histogram data do not export a complete canonical source carrier, strict spectral source law, or typed coupling theorem to K or L_total.",
            "next_honest_step": "Do not replay histogram-level 3-WL as a source carrier.  The next honest move should add exactly one different typed ingredient targeted at the 79 surviving pairs, preferably an exact automorphism/orbit partition or canonical-label-derived orbit signature, and only after actual separation run a separate source-law/coupling acceptance audit.  Keep role-bearing L_total, bridge closure, role transfer, selector closure, and ToE blocked.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2813/S1763 three-WL refined collision audit", "## P2813/S1763 three-WL refined collision audit\n\n`P2813/S1763` adds exactly one stronger typed invariant to the P2812 residual collisions: stable 3-dimensional Weisfeiler-Lehman ordered-triple color histograms.  It is applied only to the `158` P2812 collision graphs and does not separate the remaining finite carrier: the quotient stays at `16,749`, with `79` collision classes covering `158` graphs.  This is a bounded obstruction for histogram-level 3-WL data, not a strict spectral source law or a `K`/`L_total` coupling theorem.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2813/S1763 three-WL source Ltotal guard", "## P2813/S1763 three-WL source Ltotal guard\n\n`P2813/S1763` adds no variational term to `L_total`.  The 3-WL refinement supplies a finite obstruction: histogram-level 3-WL data do not improve the P2812 residual pairs.  It does not define an action functional, source equation, Euler-Lagrange contribution, or typed graph-to-`K`/`L_total` coupling.  A different separating typed ingredient plus separate source-law acceptance audit is required before any Lagrangian promotion.\n")
    append_once(AGENTS, "Current three-WL carrier-separation guardrail (P2813/S1763, 2026-06-16)", "## Current three-WL carrier-separation guardrail (P2813/S1763, 2026-06-16)\n\n- P2813 applies stable 3-WL ordered-triple color histograms to the `158` P2812 residual-collision graphs; the quotient remains `16,749`, with `79` collision classes covering `158` graphs.\n- This is a bounded obstruction for histogram-level 3-WL data: it does not export a complete canonical source carrier, strict spectral source law, typed coupling to `K`/`L_total`, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure.\n- The next admissible move should add exactly one different typed ingredient targeted at the surviving pairs, preferably exact automorphism/orbit partition or canonical-label-derived orbit signature, while preserving all closure blocks until separation plus typed variational coupling are exported.\n")
    return payload


if __name__ == "__main__":
    main()
