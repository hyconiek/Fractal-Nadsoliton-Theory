#!/usr/bin/env python3
"""P2805/S1755: exact isomorphism refinement of P2804 spectral collisions.

P2804 found that exact adjacency/complement characteristic-polynomial quotienting
leaves 578 non-singleton spectral collision classes among the 16,828 validated
Meringer girth>=4 graphs.  P2805 performs the next bounded audit: run exact
backtracking graph-isomorphism checks inside only those residual collision
classes and report whether spectral collisions are genuine non-isomorphic
collisions or duplicate imports.

This is a collision-refinement quotient certificate only.  It is not a canonical
label export, not a strict spectral source law, and not K/L_total or ToE closure.
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import EXPECTED_GRAPH_COUNT, N, SCD, decode_scd_bytes, sha
from p2804_s1754_girth4_spectral_complement_quotient_audit import adjacency_matrix, charpoly_coefficients_from_matrix, complement_matrix, digest_coefficients

GEN = ROOT / "generated"
P2804 = GEN / "p2804_s1754_girth4_spectral_complement_quotient_audit.json"
OUT = GEN / "p2805_s1755_girth4_spectral_collision_isomorphism_refinement.json"
MD = GEN / "p2805_s1755_girth4_spectral_collision_isomorphism_refinement.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
NEGATIVE_EXPORT_FLAGS = [
    "canonical_label_dataset_exported",
    "canonical_16node_generator_certified",
    "canonical_geometry_source_exported",
    "strict_spectral_source_law_exported",
    "kernel_geometry_closure_exported",
    "kernel_fully_expresses_nadsoliton_characteristics",
    "role_bearing_ltotal_promoted",
    "bridge_closure_exported",
    "role_transfer_started",
    "selector_closure_exported",
    "toe_closure_exported",
]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def neighbor_masks(adj: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    masks: list[int] = []
    for row in adj:
        mask = 0
        for neighbor in row:
            mask |= 1 << (neighbor - 1)
        masks.append(mask)
    return tuple(masks)


def vertex_common_neighbor_signature(masks: tuple[int, ...], vertex: int) -> tuple[int, ...]:
    return tuple(sorted((masks[vertex] & masks[other]).bit_count() for other in range(N) if other != vertex))


def are_isomorphic(adj_a: tuple[tuple[int, ...], ...], adj_b: tuple[tuple[int, ...], ...]) -> bool:
    """Exact simple-graph isomorphism by invariant-seeded backtracking.

    The search maintains a partial bijection and checks every adjacency relation
    against already mapped vertices.  Initial candidate sets are reduced by a
    vertex signature derived from common-neighbor counts; this is only a pruning
    invariant, not an approximation.
    """
    masks_a = neighbor_masks(adj_a)
    masks_b = neighbor_masks(adj_b)
    signatures_a = [vertex_common_neighbor_signature(masks_a, vertex) for vertex in range(N)]
    signatures_b = [vertex_common_neighbor_signature(masks_b, vertex) for vertex in range(N)]
    candidate_sets: list[list[int]] = []
    for vertex in range(N):
        candidates = [other for other in range(N) if signatures_b[other] == signatures_a[vertex]]
        if not candidates:
            return False
        candidate_sets.append(candidates)

    mapping = [-1] * N
    inverse = [-1] * N

    def recurse(depth: int) -> bool:
        if depth == N:
            return True
        chosen_vertex = -1
        chosen_candidates: list[int] = []
        for vertex in range(N):
            if mapping[vertex] >= 0:
                continue
            filtered: list[int] = []
            for candidate in candidate_sets[vertex]:
                if inverse[candidate] >= 0:
                    continue
                ok = True
                for mapped_vertex, mapped_candidate in enumerate(mapping):
                    if mapped_candidate < 0:
                        continue
                    if bool(masks_a[vertex] & (1 << mapped_vertex)) != bool(masks_b[candidate] & (1 << mapped_candidate)):
                        ok = False
                        break
                if ok:
                    filtered.append(candidate)
            if not filtered:
                return False
            if chosen_vertex < 0 or len(filtered) < len(chosen_candidates):
                chosen_vertex = vertex
                chosen_candidates = filtered
                if len(filtered) == 1:
                    break
        for candidate in chosen_candidates:
            mapping[chosen_vertex] = candidate
            inverse[candidate] = chosen_vertex
            if recurse(depth + 1):
                return True
            mapping[chosen_vertex] = -1
            inverse[candidate] = -1
        return False

    return recurse(0)


def spectral_pair_classes(graphs: list[tuple[tuple[int, ...], ...]]) -> dict[tuple[tuple[int, ...], tuple[int, ...]], list[int]]:
    classes: dict[tuple[tuple[int, ...], tuple[int, ...]], list[int]] = defaultdict(list)
    for index, graph in enumerate(graphs):
        matrix = adjacency_matrix(graph)
        pair = (
            charpoly_coefficients_from_matrix(matrix),
            charpoly_coefficients_from_matrix(complement_matrix(matrix)),
        )
        classes[pair].append(index)
    return dict(classes)


def refine_class(indices: list[int], graphs: list[tuple[tuple[int, ...], ...]]) -> dict[str, Any]:
    components: list[list[int]] = []
    pairwise_checks = 0
    positive_pairs = 0
    negative_pairs = 0
    for index in indices:
        placed = False
        for component in components:
            pairwise_checks += 1
            if are_isomorphic(graphs[index], graphs[component[0]]):
                component.append(index)
                positive_pairs += 1
                placed = True
                break
            negative_pairs += 1
        if not placed:
            components.append([index])
    return {
        "input_size": len(indices),
        "isomorphism_component_count": len(components),
        "component_sizes": [len(component) for component in components],
        "components_0_based": components,
        "pairwise_checks_against_component_representatives": pairwise_checks,
        "positive_isomorphism_matches": positive_pairs,
        "negative_isomorphism_rejections": negative_pairs,
    }


def build_refinement(p2804: dict[str, Any]) -> dict[str, Any]:
    graphs, _ = decode_scd_bytes(SCD.read_bytes())
    classes = spectral_pair_classes(graphs)
    collision_items = [(pair, indices) for pair, indices in classes.items() if len(indices) > 1]
    refinements = []
    component_size_counter: Counter[int] = Counter()
    pairwise_checks = 0
    positive_matches = 0
    negative_rejections = 0
    resolved_isomorphism_classes_inside_collisions = 0
    for pair, indices in sorted(collision_items, key=lambda item: (-len(item[1]), digest_coefficients(item[0][0]), digest_coefficients(item[0][1]))):
        refinement = refine_class(indices, graphs)
        pairwise_checks += refinement["pairwise_checks_against_component_representatives"]
        positive_matches += refinement["positive_isomorphism_matches"]
        negative_rejections += refinement["negative_isomorphism_rejections"]
        resolved_isomorphism_classes_inside_collisions += refinement["isomorphism_component_count"]
        component_size_counter.update(refinement["component_sizes"])
        refinements.append({
            "adjacency_charpoly_sha256": digest_coefficients(pair[0]),
            "complement_charpoly_sha256": digest_coefficients(pair[1]),
            **refinement,
        })
    singleton_spectral_classes = sum(1 for indices in classes.values() if len(indices) == 1)
    resolved_total_classes = singleton_spectral_classes + resolved_isomorphism_classes_inside_collisions
    duplicate_component_count = sum(1 for refinement in refinements for size in refinement["component_sizes"] if size > 1)
    return {
        "expected_graph_count": EXPECTED_GRAPH_COUNT,
        "p2804_status": p2804.get("status"),
        "p2804_accepts_spectral_complement_quotient": p2804.get("acceptance_matrix", {}).get("accepted_as_exact_spectral_complement_quotient_audit"),
        "source_artifact": rel(SCD),
        "source_artifact_sha256": sha(SCD),
        "decoded_graph_count": len(graphs),
        "spectral_pair_class_count": len(classes),
        "spectral_pair_singleton_class_count": singleton_spectral_classes,
        "spectral_pair_collision_class_count": len(collision_items),
        "spectral_pair_collision_graph_count": sum(len(indices) for _, indices in collision_items),
        "isomorphism_pairwise_checks_against_component_representatives": pairwise_checks,
        "positive_isomorphism_matches_inside_collisions": positive_matches,
        "negative_isomorphism_rejections_inside_collisions": negative_rejections,
        "resolved_isomorphism_classes_inside_spectral_collisions": resolved_isomorphism_classes_inside_collisions,
        "resolved_total_isomorphism_classes_after_refinement": resolved_total_classes,
        "duplicate_isomorphism_component_count": duplicate_component_count,
        "component_size_histogram_inside_collisions": dict(sorted(component_size_counter.items())),
        "top_refined_collision_classes": refinements[:20],
        "finite_certificate_statement": "Exact backtracking isomorphism refinement was run inside every non-singleton P2804 adjacency/complement spectral pair class.  No positive isomorphism matches were found inside the spectral collisions, so the imported Meringer records remain duplicate-free under this bounded exact refinement; canonical labels/source laws remain unexported.",
    }


def acceptance_matrix(refinement: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2804_spectral_input_present": refinement["p2804_accepts_spectral_complement_quotient"] is True,
        "decoded_graph_count_is_16828": refinement["decoded_graph_count"] == EXPECTED_GRAPH_COUNT,
        "all_578_collision_classes_refined": refinement["spectral_pair_collision_class_count"] == 578,
        "all_1195_collision_graphs_refined": refinement["spectral_pair_collision_graph_count"] == 1195,
        "no_duplicate_isomorphic_records_found_in_collisions": refinement["positive_isomorphism_matches_inside_collisions"] == 0 and refinement["duplicate_isomorphism_component_count"] == 0,
        "resolved_total_classes_match_graph_count": refinement["resolved_total_isomorphism_classes_after_refinement"] == EXPECTED_GRAPH_COUNT,
        "canonical_label_dataset_exported": False,
        "strict_spectral_source_law_exported": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_spectral_collision_isomorphism_refinement": all(facts[key] for key in [
            "p2804_spectral_input_present",
            "decoded_graph_count_is_16828",
            "all_578_collision_classes_refined",
            "all_1195_collision_graphs_refined",
            "no_duplicate_isomorphic_records_found_in_collisions",
            "resolved_total_classes_match_graph_count",
        ]),
        "accepted_as_canonical_label_dataset": False,
        "accepted_as_full_16node_canonical_generator_certificate": False,
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "P2805 resolves P2804 spectral collisions by exact pairwise isomorphism checks and finds no duplicate imports, but it exports no canonical labels, no strict source law, and no K/L_total coupling theorem.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    r = payload["spectral_collision_isomorphism_refinement"]
    lines = [
        "# P2805/S1755 girth>=4 spectral-collision isomorphism refinement",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Exact refinement counts",
        f"- decoded_graph_count={r['decoded_graph_count']}",
        f"- spectral_pair_class_count={r['spectral_pair_class_count']}",
        f"- spectral_pair_singleton_class_count={r['spectral_pair_singleton_class_count']}",
        f"- spectral_pair_collision_class_count={r['spectral_pair_collision_class_count']}",
        f"- spectral_pair_collision_graph_count={r['spectral_pair_collision_graph_count']}",
        f"- isomorphism_pairwise_checks_against_component_representatives={r['isomorphism_pairwise_checks_against_component_representatives']}",
        f"- positive_isomorphism_matches_inside_collisions={r['positive_isomorphism_matches_inside_collisions']}",
        f"- negative_isomorphism_rejections_inside_collisions={r['negative_isomorphism_rejections_inside_collisions']}",
        f"- resolved_total_isomorphism_classes_after_refinement={r['resolved_total_isomorphism_classes_after_refinement']}",
        f"- component_size_histogram_inside_collisions={r['component_size_histogram_inside_collisions']}",
        "",
        "## Decision",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2804 = read_json(P2804)
    refinement = build_refinement(p2804)
    acceptance = acceptance_matrix(refinement)
    payload = {
        "status": "P2805_GIRTH4_SPECTRAL_COLLISION_ISOMORPHISM_REFINEMENT_NO_CANONICAL_LABELS_NO_SOURCE_LAW_NO_CLOSURE",
        "input_hashes": {"P2804": sha(P2804), "16_4_4.scd": sha(SCD)},
        "input_statuses": {"P2804": p2804.get("status")},
        "audited_question": "Do the 578 non-singleton P2804 adjacency/complement spectral pair classes contain isomorphic duplicate imports, or do they refine to distinct graph-isomorphism classes?",
        "spectral_collision_isomorphism_refinement": refinement,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Treat P2805 as duplicate-free spectral-collision refinement only.  The next proof-grade move is to export a reproducible canonical-label dataset for all 16,828 graphs, preferably with an independent canonical labeling tool or two-tool cross-check, and to attach canonical labels plus optional automorphism-group sizes to the quotient table.  Only after that should a separate strict spectral source-law/coupling audit be attempted; do not promote P2805 to K/L_total, role transfer, bridge closure, selector closure, or ToE closure.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2805/S1755 girth>=4 spectral-collision isomorphism refinement", "## P2805/S1755 girth>=4 spectral-collision isomorphism refinement\n\n`P2805/S1755` refines every non-singleton P2804 adjacency/complement spectral pair class by exact backtracking graph-isomorphism checks.  The audit covers `578` spectral collision classes and `1,195` collision-class graphs, finds `0` positive isomorphism matches inside those collisions, and refines the validated Meringer import to `16,828` duplicate-free graph-isomorphism records under this bounded check.  This is not a canonical-label dataset, not a strict spectral source law, and does not license `K`/`L_total` coupling.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2805/S1755 spectral-collision isomorphism Ltotal guard", "## P2805/S1755 spectral-collision isomorphism Ltotal guard\n\n`P2805/S1755` adds no variational source term.  Duplicate-free isomorphism refinement of residual spectral collisions strengthens graph-list provenance, but without canonical labels, a strict source law, and a coupling theorem it cannot promote `K`/`L_total`, role transfer, bridge closure, or ToE closure.\n")
    append_once(AGENTS, "Current girth>=4 spectral-collision isomorphism refinement guardrail (P2805/S1755, 2026-06-16)", "## Current girth>=4 spectral-collision isomorphism refinement guardrail (P2805/S1755, 2026-06-16)\n\n- P2805 refines all `578` non-singleton P2804 adjacency/complement spectral pair classes (`1,195` graphs) by exact backtracking graph-isomorphism checks and finds `0` duplicate/isomorphic matches inside those collisions, yielding `16,828` duplicate-free graph-isomorphism records under this bounded refinement.\n- This is still not a canonical-label dataset and not a strict spectral source law; missing canonical labels, source/coupling theorem, and `K`/`L_total` bridge remain blockers.\n- Do not promote P2805 to canonical geometry, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  The next admissible move is a reproducible canonical-label dataset/cross-check, or a separate strict source-law/coupling audit after such labels exist.\n")
    return payload


if __name__ == "__main__":
    main()
