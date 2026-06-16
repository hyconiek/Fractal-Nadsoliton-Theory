#!/usr/bin/env python3
"""P2812/S1762: two-WL refinement of the P2811 residual collisions.

P2811 left 132 local-motif refined collision classes covering 269 graphs.  P2812
adds exactly one further typed invariant: the stable 2-dimensional
Weisfeiler-Lehman color histogram, applied only as a source-neutral refinement
inside those residual classes.  This is an orbit-aware-ish combinatorial
signature, not a variational source term.
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import EXPECTED_GRAPH_COUNT, N, SCD, decode_scd_bytes, sha
from p2804_s1754_girth4_spectral_complement_quotient_audit import adjacency_matrix, charpoly_coefficients_from_matrix, complement_matrix
from p2811_s1761_local_motif_refined_source_candidate_audit import local_motif_profile, stable_digest

GEN = ROOT / "generated"
P2811 = GEN / "p2811_s1761_local_motif_refined_source_candidate_audit.json"
OUT = GEN / "p2812_s1762_two_wl_refined_collision_audit.json"
MD = GEN / "p2812_s1762_two_wl_refined_collision_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def p2811_refined_key(graph: tuple[tuple[int, ...], ...]) -> str:
    matrix = adjacency_matrix(graph)
    spectral_pair = (
        charpoly_coefficients_from_matrix(matrix),
        charpoly_coefficients_from_matrix(complement_matrix(matrix)),
    )
    return stable_digest({"spectral_pair": spectral_pair, "motif_profile": local_motif_profile(graph)})


def two_wl_histogram(adj: tuple[tuple[int, ...], ...], max_iterations: int = 20) -> tuple[tuple[int, int], ...]:
    matrix = adjacency_matrix(adj)
    colors = {
        (i, j): 0 if i == j else 1 if matrix[i, j] else 2
        for i in range(N)
        for j in range(N)
    }
    iterations = 0
    for iteration in range(1, max_iterations + 1):
        updates = []
        for i in range(N):
            for j in range(N):
                transition_multiset = tuple(sorted((colors[(i, k)], colors[(k, j)]) for k in range(N)))
                updates.append(((i, j), (colors[(i, j)], transition_multiset)))
        palette = {value: color for color, value in enumerate(sorted({value for _, value in updates}))}
        new_colors = {pair: palette[value] for pair, value in updates}
        iterations = iteration
        if new_colors == colors:
            break
        colors = new_colors
    histogram = Counter(colors.values())
    return tuple(sorted(histogram.items())) + ((-1, iterations),)


def build_audit(p2811: dict[str, Any]) -> dict[str, Any]:
    graphs, _ = decode_scd_bytes(SCD.read_bytes())
    p2811_classes: dict[str, list[int]] = defaultdict(list)
    for index, graph in enumerate(graphs):
        p2811_classes[p2811_refined_key(graph)].append(index)

    p2812_classes: dict[tuple[str, str | None], list[int]] = defaultdict(list)
    wl_rows: list[dict[str, Any]] = []
    wl_computed_graph_count = 0
    for key, indices in p2811_classes.items():
        if len(indices) == 1:
            p2812_classes[(key, None)].extend(indices)
            continue
        for index in indices:
            wl_digest = stable_digest(two_wl_histogram(graphs[index]))
            p2812_classes[(key, wl_digest)].append(index)
            wl_computed_graph_count += 1
            if len(wl_rows) < 12:
                wl_rows.append({"graph_index_0_based": index, "p2811_refined_key_sha256": key, "two_wl_histogram_sha256": wl_digest})

    collision_classes = {key: indices for key, indices in p2812_classes.items() if len(indices) > 1}
    class_size_histogram = Counter(len(indices) for indices in p2812_classes.values())
    top_collision_classes = sorted(collision_classes.items(), key=lambda item: (-len(item[1]), str(item[0])))[:20]
    p2811_audit = p2811["local_motif_refined_audit"]
    return {
        "candidate_invariant": "stable 2-dimensional Weisfeiler-Lehman ordered-pair color histogram, applied inside P2811 residual collision classes",
        "decoded_graph_count": len(graphs),
        "p2811_refined_class_count": p2811_audit["refined_class_count"],
        "p2811_refined_collision_class_count": p2811_audit["refined_collision_class_count"],
        "p2811_refined_collision_graph_count": p2811_audit["refined_collision_graph_count"],
        "two_wl_computed_graph_count": wl_computed_graph_count,
        "two_wl_refined_class_count": len(p2812_classes),
        "two_wl_collision_class_count": len(collision_classes),
        "two_wl_collision_graph_count": sum(len(indices) for indices in collision_classes.values()),
        "two_wl_max_class_size": max(class_size_histogram) if class_size_histogram else 0,
        "class_size_histogram": dict(sorted(class_size_histogram.items())),
        "remaining_defect_canonical_minus_two_wl": EXPECTED_GRAPH_COUNT - len(p2812_classes),
        "defect_reduction_vs_p2811": p2811_audit["remaining_defect_canonical_minus_refined"] - (EXPECTED_GRAPH_COUNT - len(p2812_classes)),
        "top_collision_classes": [
            {"p2811_refined_key_sha256": key[0], "two_wl_histogram_sha256": key[1], "size": len(indices), "example_indices_0_based": indices[:8]}
            for key, indices in top_collision_classes
        ],
        "sample_two_wl_rows": wl_rows,
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "decoded_graph_count_is_16828": audit["decoded_graph_count"] == EXPECTED_GRAPH_COUNT,
        "two_wl_applied_only_to_p2811_collision_graphs": audit["two_wl_computed_graph_count"] == audit["p2811_refined_collision_graph_count"],
        "two_wl_improves_over_p2811": audit["two_wl_refined_class_count"] > audit["p2811_refined_class_count"],
        "two_wl_still_below_canonical_target": audit["two_wl_refined_class_count"] < EXPECTED_GRAPH_COUNT,
        "remaining_collision_classes_exist": audit["two_wl_collision_class_count"] > 0,
        "strict_spectral_source_law_exported": False,
        "typed_coupling_to_K_or_Ltotal_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_two_wl_refinement_audit": all(facts[key] for key in [
            "decoded_graph_count_is_16828",
            "two_wl_applied_only_to_p2811_collision_graphs",
            "two_wl_improves_over_p2811",
            "two_wl_still_below_canonical_target",
            "remaining_collision_classes_exist",
        ]),
        "accepted_as_complete_canonical_source_carrier": False,
        "accepted_as_strict_source_law_or_coupling": False,
        "missing_for_promotion": [key for key, value in facts.items() if not value],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["two_wl_refinement_audit"]
    acceptance = payload["acceptance_matrix"]
    lines = [
        "# P2812/S1762 two-WL refined collision audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Candidate invariant",
        audit["candidate_invariant"],
        "",
        "## Exact quotient counts",
        f"- decoded_graph_count={audit['decoded_graph_count']}",
        f"- p2811_refined_class_count={audit['p2811_refined_class_count']}",
        f"- p2811_refined_collision_graph_count={audit['p2811_refined_collision_graph_count']}",
        f"- two_wl_computed_graph_count={audit['two_wl_computed_graph_count']}",
        f"- two_wl_refined_class_count={audit['two_wl_refined_class_count']}",
        f"- two_wl_collision_class_count={audit['two_wl_collision_class_count']}",
        f"- two_wl_collision_graph_count={audit['two_wl_collision_graph_count']}",
        f"- two_wl_max_class_size={audit['two_wl_max_class_size']}",
        f"- remaining_defect_canonical_minus_two_wl={audit['remaining_defect_canonical_minus_two_wl']}",
        f"- defect_reduction_vs_p2811={audit['defect_reduction_vs_p2811']}",
        "",
        "## Decision",
        f"- accepted_as_two_wl_refinement_audit={acceptance['accepted_as_two_wl_refinement_audit']}",
        f"- accepted_as_complete_canonical_source_carrier={acceptance['accepted_as_complete_canonical_source_carrier']}",
        f"- accepted_as_strict_source_law_or_coupling={acceptance['accepted_as_strict_source_law_or_coupling']}",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2811 = read_json(P2811)
    audit = build_audit(p2811)
    payload: dict[str, Any] = {
        "status": "P2812_TWO_WL_REFINED_COLLISION_AUDIT_REFINES_BUT_STILL_OBSTRUCTED_NO_CLOSURE",
        "input_hashes": {"P2811": sha(P2811), "16_4_4.scd": sha(SCD)},
        "input_statuses": {"P2811": p2811.get("status")},
        "two_wl_refinement_audit": audit,
        "decision": {
            "negative_export_flags": {
                "complete_canonical_source_carrier_exported": False,
                "strict_spectral_source_law_exported": False,
                "typed_coupling_to_K_or_Ltotal_exported": False,
                "canonical_geometry_closure_exported": False,
                "role_bearing_ltotal_promoted": False,
                "bridge_closure_exported": False,
                "role_transfer_started": False,
                "selector_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2812 applies one additional typed invariant, stable 2-WL ordered-pair color histograms, only to the 269 P2811 residual-collision graphs.  It improves the refined quotient from 16,691 to 16,749 classes and reduces the remaining defect from 137 to 79, but 79 two-WL collision classes covering 158 graphs remain.  Thus 2-WL refinement is useful computational evidence, not a complete canonical source carrier, strict source law, or K/L_total coupling.",
            "next_honest_step": "Target the remaining 79 two-WL collision pairs with exactly one stronger invariant, preferably a 3-WL/k-WL signature, exact automorphism-orbit partition, or canonical-label-derived orbit signature.  Keep K/L_total, role transfer, bridge closure, selector closure, and ToE promotion blocked until separation and an explicit typed variational coupling are exported.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2812/S1762 two-WL refined collision audit", "## P2812/S1762 two-WL refined collision audit\n\n`P2812/S1762` adds exactly one stronger typed invariant to the P2811 residual collisions: stable 2-dimensional Weisfeiler-Lehman ordered-pair color histograms.  It is applied only to the `269` P2811 collision graphs, improves the quotient from `16,691` to `16,749` classes, and reduces the remaining defect from `137` to `79`; however `79` two-WL collision classes covering `158` graphs remain.  This is refinement evidence, not a complete canonical source carrier, strict spectral source law, or `K`/`L_total` coupling.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2812/S1762 two-WL source Ltotal guard", "## P2812/S1762 two-WL source Ltotal guard\n\n`P2812/S1762` adds no variational term to `L_total`.  The 2-WL refinement narrows the graph-source carrier obstruction but still leaves collision pairs and exports no typed graph-to-`K`/`L_total` coupling theorem.  Future Lagrangian work must either separate the remaining pairs with a stronger invariant or supply an explicit source theorem independently.\n")
    append_once(AGENTS, "Current two-WL refined collision guardrail (P2812/S1762, 2026-06-16)", "## Current two-WL refined collision guardrail (P2812/S1762, 2026-06-16)\n\n- P2812 applies stable 2-WL ordered-pair color histograms to the `269` P2811 residual-collision graphs, improving the quotient from `16,691` to `16,749` classes and reducing the remaining defect from `137` to `79`.\n- The refinement still leaves `79` collision classes covering `158` graphs, so it does not export a complete canonical source carrier, strict spectral source law, `K`/`L_total` coupling, selector closure, bridge closure, role transfer, or ToE closure.\n- A next admissible move should target those remaining pairs with exactly one stronger invariant, such as 3-WL/k-WL, exact automorphism-orbit partition, or canonical-label-derived orbit signature, while preserving all closure blocks until typed variational coupling is exported.\n")
    return payload


if __name__ == "__main__":
    main()
