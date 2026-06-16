#!/usr/bin/env python3
"""P2814/S1764: exact automorphism-orbit audit of P2813 residual pairs.

P2813 showed that histogram-level 3-WL does not separate the 79 remaining
P2812/P2813 collision pairs.  P2814 adds exactly one different typed ingredient:
an exact automorphism-group/order and vertex-orbit partition signature, computed
only on those residual graphs by a deterministic backtracking automorphism
enumerator.  This is a finite graph-source carrier obstruction test, not a
source law or variational coupling.
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
from p2812_s1762_two_wl_refined_collision_audit import p2811_refined_key, read_json, two_wl_histogram
from p2813_s1763_three_wl_refined_collision_audit import three_wl_histogram

GEN = ROOT / "generated"
P2813 = GEN / "p2813_s1763_three_wl_refined_collision_audit.json"
OUT = GEN / "p2814_s1764_exact_automorphism_orbit_residual_audit.json"
MD = GEN / "p2814_s1764_exact_automorphism_orbit_residual_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def zero_based_adjacency(graph: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(sorted(neighbor - 1 for neighbor in neighbors)) for neighbors in graph)


def vertex_refinement_signature(graph0: tuple[tuple[int, ...], ...], vertex: int) -> tuple[Any, ...]:
    neighbors = set(graph0[vertex])
    common_nonedge_counts = tuple(sorted(
        len(neighbors & set(graph0[other]))
        for other in range(N)
        if other != vertex and other not in neighbors
    ))
    distances = [99] * N
    distances[vertex] = 0
    queue = [vertex]
    for current in queue:
        for other in graph0[current]:
            if distances[other] == 99:
                distances[other] = distances[current] + 1
                queue.append(other)
    neighbor_edge_count = sum(1 for left in neighbors for right in neighbors if left < right and right in graph0[left])
    return (common_nonedge_counts, tuple(sorted(distances)), neighbor_edge_count)


def exact_automorphism_orbit_signature(graph: tuple[tuple[int, ...], ...], search_limit: int = 10_000_000) -> dict[str, Any]:
    graph0 = zero_based_adjacency(graph)
    adjacency = [[False] * N for _ in range(N)]
    for vertex, neighbors in enumerate(graph0):
        for neighbor in neighbors:
            adjacency[vertex][neighbor] = True

    vertex_signatures = [vertex_refinement_signature(graph0, vertex) for vertex in range(N)]
    signature_classes: dict[tuple[Any, ...], list[int]] = defaultdict(list)
    for vertex, signature in enumerate(vertex_signatures):
        signature_classes[signature].append(vertex)
    initial_candidates = {vertex: set(signature_classes[vertex_signatures[vertex]]) for vertex in range(N)}

    parent = list(range(N))
    mapping: dict[int, int] = {}
    used: set[int] = set()
    automorphism_count = 0
    truncated = False

    def find(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    def union(left: int, right: int) -> None:
        root_left = find(left)
        root_right = find(right)
        if root_left != root_right:
            parent[root_right] = root_left

    def compatible(source: int, target: int) -> bool:
        return all(adjacency[source][assigned_source] == adjacency[target][assigned_target] for assigned_source, assigned_target in mapping.items())

    def recurse() -> None:
        nonlocal automorphism_count, truncated
        if truncated:
            return
        if len(mapping) == N:
            automorphism_count += 1
            for source, target in mapping.items():
                union(source, target)
            if automorphism_count >= search_limit:
                truncated = True
            return

        best_source: int | None = None
        best_candidates: list[int] | None = None
        for source in range(N):
            if source in mapping:
                continue
            candidates = [target for target in sorted(initial_candidates[source] - used) if compatible(source, target)]
            if not candidates:
                return
            if best_candidates is None or len(candidates) < len(best_candidates):
                best_source = source
                best_candidates = candidates
                if len(candidates) == 1:
                    break

        assert best_source is not None and best_candidates is not None
        for target in best_candidates:
            mapping[best_source] = target
            used.add(target)
            recurse()
            used.remove(target)
            del mapping[best_source]
            if truncated:
                return

    recurse()
    orbit_classes: dict[int, list[int]] = defaultdict(list)
    for vertex in range(N):
        orbit_classes[find(vertex)].append(vertex)
    orbit_partition = tuple(sorted(tuple(sorted(orbit)) for orbit in orbit_classes.values()))
    orbit_size_multiset = tuple(sorted(len(orbit) for orbit in orbit_partition))
    return {
        "automorphism_group_order": automorphism_count,
        "orbit_count": len(orbit_partition),
        "orbit_size_multiset": orbit_size_multiset,
        "orbit_partition": orbit_partition,
        "truncated": truncated,
        "initial_vertex_signature_class_size_histogram": dict(sorted(Counter(len(v) for v in signature_classes.values()).items())),
    }


def p2813_residual_key(graph: tuple[tuple[int, ...], ...]) -> tuple[str, str, str]:
    p2811_key = p2811_refined_key(graph)
    two_wl_key = stable_digest(two_wl_histogram(graph))
    three_wl_key = stable_digest(three_wl_histogram(graph))
    return (p2811_key, two_wl_key, three_wl_key)


def build_audit(p2813: dict[str, Any]) -> dict[str, Any]:
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
    for key, indices in p2812_classes.items():
        if len(indices) == 1:
            p2813_classes[(key, None)].extend(indices)
            continue
        for index in indices:
            p2813_classes[(key, stable_digest(three_wl_histogram(graphs[index])))].append(index)

    p2814_classes: dict[tuple[tuple[tuple[str, str | None], str | None], str | None], list[int]] = defaultdict(list)
    rows: list[dict[str, Any]] = []
    signature_counter: Counter[str] = Counter()
    aut_order_counter: Counter[int] = Counter()
    orbit_size_counter: Counter[tuple[int, ...]] = Counter()
    truncated_count = 0
    computed = 0

    for key, indices in p2813_classes.items():
        if len(indices) == 1:
            p2814_classes[(key, None)].extend(indices)
            continue
        for index in indices:
            signature = exact_automorphism_orbit_signature(graphs[index])
            digest = stable_digest(signature)
            p2814_classes[(key, digest)].append(index)
            signature_counter[digest] += 1
            aut_order_counter[signature["automorphism_group_order"]] += 1
            orbit_size_counter[tuple(signature["orbit_size_multiset"])] += 1
            truncated_count += int(signature["truncated"])
            computed += 1
            if len(rows) < 12:
                rows.append({
                    "graph_index_0_based": index,
                    "p2813_key_sha256": stable_digest(key),
                    "automorphism_orbit_signature_sha256": digest,
                    "automorphism_group_order": signature["automorphism_group_order"],
                    "orbit_size_multiset": signature["orbit_size_multiset"],
                    "truncated": signature["truncated"],
                })

    collisions = {key: indices for key, indices in p2814_classes.items() if len(indices) > 1}
    class_size_histogram = Counter(len(indices) for indices in p2814_classes.values())
    return {
        "candidate_invariant": "exact automorphism-group order and vertex-orbit partition signature, applied inside P2813 residual collision classes",
        "decoded_graph_count": len(graphs),
        "p2813_refined_class_count": p2813["three_wl_refinement_audit"]["three_wl_refined_class_count"],
        "p2813_collision_class_count": p2813["three_wl_refinement_audit"]["three_wl_collision_class_count"],
        "p2813_collision_graph_count": p2813["three_wl_refinement_audit"]["three_wl_collision_graph_count"],
        "automorphism_orbit_computed_graph_count": computed,
        "automorphism_orbit_refined_class_count": len(p2814_classes),
        "automorphism_orbit_collision_class_count": len(collisions),
        "automorphism_orbit_collision_graph_count": sum(len(indices) for indices in collisions.values()),
        "automorphism_orbit_max_class_size": max(class_size_histogram) if class_size_histogram else 0,
        "class_size_histogram": dict(sorted(class_size_histogram.items())),
        "remaining_defect_canonical_minus_automorphism_orbit": EXPECTED_GRAPH_COUNT - len(p2814_classes),
        "defect_reduction_vs_p2813": p2813["three_wl_refinement_audit"]["remaining_defect_canonical_minus_three_wl"] - (EXPECTED_GRAPH_COUNT - len(p2814_classes)),
        "computed_signature_distinct_count": len(signature_counter),
        "computed_signature_histogram": dict(sorted(signature_counter.items())),
        "automorphism_group_order_histogram_on_residual": dict(sorted(aut_order_counter.items())),
        "orbit_size_multiset_histogram_on_residual": {str(key): value for key, value in sorted(orbit_size_counter.items())},
        "truncated_automorphism_search_count": truncated_count,
        "sample_automorphism_orbit_rows": rows,
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "decoded_graph_count_is_16828": audit["decoded_graph_count"] == EXPECTED_GRAPH_COUNT,
        "automorphism_orbit_applied_only_to_p2813_collision_graphs": audit["automorphism_orbit_computed_graph_count"] == audit["p2813_collision_graph_count"],
        "automorphism_search_not_truncated": audit["truncated_automorphism_search_count"] == 0,
        "automorphism_orbit_improves_over_p2813": audit["automorphism_orbit_refined_class_count"] > audit["p2813_refined_class_count"],
        "automorphism_orbit_reaches_canonical_target": audit["automorphism_orbit_refined_class_count"] == EXPECTED_GRAPH_COUNT,
        "remaining_collision_classes_exist": audit["automorphism_orbit_collision_class_count"] > 0,
        "strict_spectral_source_law_exported": False,
        "typed_coupling_to_K_or_Ltotal_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_exact_automorphism_orbit_obstruction_audit": all(facts[key] for key in [
            "decoded_graph_count_is_16828",
            "automorphism_orbit_applied_only_to_p2813_collision_graphs",
            "automorphism_search_not_truncated",
            "remaining_collision_classes_exist",
        ]),
        "accepted_as_complete_canonical_source_carrier": facts["automorphism_orbit_reaches_canonical_target"],
        "accepted_as_strict_source_law_or_coupling": False,
        "missing_for_promotion": [key for key, value in facts.items() if not value],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["automorphism_orbit_audit"]
    acceptance = payload["acceptance_matrix"]
    lines = [
        "# P2814/S1764 exact automorphism-orbit residual audit", "", f"Status: `{payload['status']}`", "",
        "## Candidate invariant", audit["candidate_invariant"], "", "## Exact quotient counts",
        f"- decoded_graph_count={audit['decoded_graph_count']}",
        f"- p2813_refined_class_count={audit['p2813_refined_class_count']}",
        f"- p2813_collision_graph_count={audit['p2813_collision_graph_count']}",
        f"- automorphism_orbit_computed_graph_count={audit['automorphism_orbit_computed_graph_count']}",
        f"- automorphism_orbit_refined_class_count={audit['automorphism_orbit_refined_class_count']}",
        f"- automorphism_orbit_collision_class_count={audit['automorphism_orbit_collision_class_count']}",
        f"- automorphism_orbit_collision_graph_count={audit['automorphism_orbit_collision_graph_count']}",
        f"- automorphism_orbit_max_class_size={audit['automorphism_orbit_max_class_size']}",
        f"- remaining_defect_canonical_minus_automorphism_orbit={audit['remaining_defect_canonical_minus_automorphism_orbit']}",
        f"- defect_reduction_vs_p2813={audit['defect_reduction_vs_p2813']}",
        f"- automorphism_group_order_histogram_on_residual={audit['automorphism_group_order_histogram_on_residual']}",
        f"- orbit_size_multiset_histogram_on_residual={audit['orbit_size_multiset_histogram_on_residual']}",
        f"- truncated_automorphism_search_count={audit['truncated_automorphism_search_count']}", "", "## Decision",
        f"- accepted_as_exact_automorphism_orbit_obstruction_audit={acceptance['accepted_as_exact_automorphism_orbit_obstruction_audit']}",
        f"- accepted_as_complete_canonical_source_carrier={acceptance['accepted_as_complete_canonical_source_carrier']}",
        f"- accepted_as_strict_source_law_or_coupling={acceptance['accepted_as_strict_source_law_or_coupling']}", "", "## Boundary",
        payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2813 = read_json(P2813)
    audit = build_audit(p2813)
    payload: dict[str, Any] = {
        "status": "P2814_EXACT_AUTOMORPHISM_ORBIT_RESIDUAL_OBSTRUCTION_NO_CLOSURE",
        "input_hashes": {"P2813": sha(P2813), "16_4_4.scd": sha(SCD)},
        "input_statuses": {"P2813": p2813.get("status")},
        "automorphism_orbit_audit": audit,
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
            "reason": "P2814 applies exactly one different typed ingredient, exact automorphism-group/order and vertex-orbit partition signatures, only to the 158 P2813 residual-collision graphs.  The computation is exact and untruncated; it improves the quotient from 16,749 to 16,771 classes and reduces the defect from 79 to 57.  However 57 collision pairs covering 114 graphs remain, so automorphism/orbit data still do not export a complete canonical source carrier, strict spectral source law, or typed coupling theorem to K or L_total.",
            "next_honest_step": "Do not replay automorphism/orbit signatures as closure.  The next honest move should add exactly one non-automorphism typed ingredient targeted at the 57 surviving pairs, preferably a small-subgraph extension/deletion response profile or a typed graph-to-K/L_total coupling ansatz with a falsifiable acceptance matrix.  Keep role-bearing L_total, bridge closure, role transfer, selector closure, and ToE blocked.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2814/S1764 exact automorphism-orbit residual audit", "## P2814/S1764 exact automorphism-orbit residual audit\n\n`P2814/S1764` adds exactly one different typed ingredient after P2813: exact automorphism-group/order and vertex-orbit partition signatures on the `158` residual-collision graphs.  The search is untruncated, but the quotient improves from `16,749` to `16,771`, reducing the defect from `79` to `57`; however `57` collision classes covering `114` graphs remain.  This is a bounded partial-refinement result for automorphism/orbit data, not a strict spectral source law or a `K`/`L_total` coupling theorem.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2814/S1764 automorphism-orbit Ltotal guard", "## P2814/S1764 automorphism-orbit Ltotal guard\n\n`P2814/S1764` adds no variational term to `L_total`.  The exact automorphism/orbit audit improves the P2813 carrier obstruction but still leaves `57` collision pairs covering `114` graphs.  It does not define an action functional, source equation, Euler-Lagrange contribution, or typed graph-to-`K`/`L_total` coupling.\n")
    append_once(AGENTS, "Current exact automorphism-orbit residual guardrail (P2814/S1764, 2026-06-16)", "## Current exact automorphism-orbit residual guardrail (P2814/S1764, 2026-06-16)\n\n- P2814 applies exact automorphism-group/order and vertex-orbit partition signatures to the `158` P2813 residual-collision graphs; the search is untruncated and improves the quotient from `16,749` to `16,771`, but `57` collision classes covering `114` graphs remain.\n- Automorphism/orbit data therefore do not export a complete canonical source carrier, strict spectral source law, typed coupling to `K`/`L_total`, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure.\n- The next admissible move should add exactly one non-automorphism typed ingredient targeted at the surviving `57` pairs, such as a small-subgraph extension/deletion response profile, while preserving all closure blocks until separation plus typed variational coupling are exported.\n")
    return payload


if __name__ == "__main__":
    main()
