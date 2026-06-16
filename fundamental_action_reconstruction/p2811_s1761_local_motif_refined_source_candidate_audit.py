#!/usr/bin/env python3
"""P2811/S1761: local motif-refined source candidate audit.

P2810 ruled out spectral-pair-only source/action functionals.  P2811 adds one
strictly richer, still source-neutral ingredient: an exact local motif/distance
profile on each decoded Meringer graph.  The audit asks whether

    (adjacency/complement characteristic-polynomial pair, local motif profile)

is fine enough to act as a candidate graph-source carrier for the P2808 canonical
quotient.  It improves the quotient but still leaves non-isomorphic collisions,
so it remains an obstruction rather than a strict source law or K/L_total
coupling theorem.
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
from p2804_s1754_girth4_spectral_complement_quotient_audit import adjacency_matrix, charpoly_coefficients_from_matrix, complement_matrix

GEN = ROOT / "generated"
P2810 = GEN / "p2810_s1760_spectral_only_source_functional_obstruction.json"
OUT = GEN / "p2811_s1761_local_motif_refined_source_candidate_audit.json"
MD = GEN / "p2811_s1761_local_motif_refined_source_candidate_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def stable_digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, separators=(",", ":"), sort_keys=True).encode("utf-8")).hexdigest()


def bfs_distance_profile(adj: tuple[tuple[int, ...], ...], source: int) -> tuple[int, ...]:
    seen = {source: 0}
    queue = [source]
    for vertex in queue:
        for neighbor in adj[vertex]:
            nxt = neighbor - 1
            if nxt not in seen:
                seen[nxt] = seen[vertex] + 1
                queue.append(nxt)
    counts = Counter(seen.values())
    return tuple(counts.get(distance, 0) for distance in range(1, N))


def local_motif_profile(adj: tuple[tuple[int, ...], ...]) -> dict[str, Any]:
    matrix = adjacency_matrix(adj)
    neighborhoods = [set(neighbor - 1 for neighbor in row) for row in adj]
    nonedge_common_neighbor_histogram: Counter[int] = Counter()
    vertex_nonedge_common_neighbor_profiles: list[tuple[int, ...]] = []
    vertex_c4_participation = [0] * N

    for i in range(N):
        local_counter: Counter[int] = Counter()
        for j in range(N):
            if i == j or matrix[i, j]:
                continue
            common_count = len(neighborhoods[i] & neighborhoods[j])
            local_counter[common_count] += 1
            if i < j:
                nonedge_common_neighbor_histogram[common_count] += 1
                if common_count >= 2:
                    # Each unordered pair of common neighbors gives a 4-cycle
                    # through the opposite nonedge (i,j).  Record endpoint
                    # participation as a local profile, not as a source term.
                    vertex_c4_participation[i] += common_count * (common_count - 1) // 2
                    vertex_c4_participation[j] += common_count * (common_count - 1) // 2
        vertex_nonedge_common_neighbor_profiles.append(tuple(local_counter.get(k, 0) for k in range(N)))

    distance_profiles = [bfs_distance_profile(adj, source) for source in range(N)]
    distance_pair_histogram: Counter[int] = Counter()
    for profile in distance_profiles:
        for distance, count in enumerate(profile, start=1):
            distance_pair_histogram[distance] += count
    # Distance profiles count ordered source-target pairs; divide by two for an
    # unordered-pair histogram because the graphs are undirected.
    unordered_distance_pair_histogram = tuple(distance_pair_histogram.get(distance, 0) // 2 for distance in range(1, N))

    return {
        "nonedge_common_neighbor_histogram": tuple(nonedge_common_neighbor_histogram.get(k, 0) for k in range(N)),
        "unordered_distance_pair_histogram": unordered_distance_pair_histogram,
        "sorted_vertex_nonedge_common_neighbor_profiles": tuple(sorted(vertex_nonedge_common_neighbor_profiles)),
        "sorted_vertex_c4_participation": tuple(sorted(vertex_c4_participation)),
        "sorted_vertex_distance_profiles": tuple(sorted(distance_profiles)),
    }


def build_audit(p2810: dict[str, Any]) -> dict[str, Any]:
    graphs, _ = decode_scd_bytes(SCD.read_bytes())
    refined_classes: dict[str, list[int]] = defaultdict(list)
    sample_rows: list[dict[str, Any]] = []
    for index, graph in enumerate(graphs):
        matrix = adjacency_matrix(graph)
        spectral_pair = (
            charpoly_coefficients_from_matrix(matrix),
            charpoly_coefficients_from_matrix(complement_matrix(matrix)),
        )
        motif_profile = local_motif_profile(graph)
        refined_key_digest = stable_digest({"spectral_pair": spectral_pair, "motif_profile": motif_profile})
        refined_classes[refined_key_digest].append(index)
        if len(sample_rows) < 8:
            sample_rows.append({
                "graph_index_0_based": index,
                "refined_key_sha256": refined_key_digest,
                "motif_profile_sha256": stable_digest(motif_profile),
                "nonedge_common_neighbor_histogram": list(motif_profile["nonedge_common_neighbor_histogram"]),
                "unordered_distance_pair_histogram": list(motif_profile["unordered_distance_pair_histogram"]),
            })

    collision_classes = {key: indices for key, indices in refined_classes.items() if len(indices) > 1}
    class_size_histogram = Counter(len(indices) for indices in refined_classes.values())
    top_collision_classes = sorted(collision_classes.items(), key=lambda item: (-len(item[1]), item[0]))[:20]
    p2810_matrix = p2810["spectral_only_obstruction_matrix"]
    return {
        "candidate_class": "F(G)=Phi(chi_A_G, chi_A_complement_G, local_motif_distance_profile_G)",
        "local_motif_profile_components": [
            "nonedge common-neighbor histogram",
            "unordered distance-pair histogram",
            "sorted vertex nonedge-common-neighbor profiles",
            "sorted vertex 4-cycle participation counts",
            "sorted vertex distance profiles",
        ],
        "decoded_graph_count": len(graphs),
        "canonical_target_class_count_from_p2810": p2810_matrix["canonical_certificate_class_count_from_p2808"],
        "spectral_pair_class_count_from_p2810": p2810_matrix["spectral_pair_class_count_from_p2804"],
        "refined_class_count": len(refined_classes),
        "refined_collision_class_count": len(collision_classes),
        "refined_collision_graph_count": sum(len(indices) for indices in collision_classes.values()),
        "refined_max_class_size": max(class_size_histogram) if class_size_histogram else 0,
        "class_size_histogram": dict(sorted(class_size_histogram.items())),
        "remaining_defect_canonical_minus_refined": EXPECTED_GRAPH_COUNT - len(refined_classes),
        "defect_reduction_vs_spectral_pair_only": p2810_matrix["quotient_defect_canonical_minus_spectral_pair"] - (EXPECTED_GRAPH_COUNT - len(refined_classes)),
        "top_collision_classes": [
            {"refined_key_sha256": key, "size": len(indices), "example_indices_0_based": indices[:8]}
            for key, indices in top_collision_classes
        ],
        "sample_rows": sample_rows,
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "decoded_graph_count_is_16828": audit["decoded_graph_count"] == EXPECTED_GRAPH_COUNT,
        "refinement_improves_over_spectral_pair_only": audit["refined_class_count"] > audit["spectral_pair_class_count_from_p2810"],
        "refinement_still_below_canonical_target": audit["refined_class_count"] < audit["canonical_target_class_count_from_p2810"],
        "remaining_collision_classes_exist": audit["refined_collision_class_count"] > 0,
        "strict_spectral_source_law_exported": False,
        "typed_coupling_to_K_or_Ltotal_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_local_motif_refinement_audit": all(facts[key] for key in [
            "decoded_graph_count_is_16828",
            "refinement_improves_over_spectral_pair_only",
            "refinement_still_below_canonical_target",
            "remaining_collision_classes_exist",
        ]),
        "accepted_as_complete_canonical_source_carrier": False,
        "accepted_as_strict_source_law_or_coupling": False,
        "missing_for_promotion": [key for key, value in facts.items() if not value],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["local_motif_refined_audit"]
    acceptance = payload["acceptance_matrix"]
    lines = [
        "# P2811/S1761 local motif-refined source candidate audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Candidate class",
        audit["candidate_class"],
        "",
        "## Exact quotient counts",
        f"- decoded_graph_count={audit['decoded_graph_count']}",
        f"- spectral_pair_class_count_from_p2810={audit['spectral_pair_class_count_from_p2810']}",
        f"- refined_class_count={audit['refined_class_count']}",
        f"- canonical_target_class_count_from_p2810={audit['canonical_target_class_count_from_p2810']}",
        f"- refined_collision_class_count={audit['refined_collision_class_count']}",
        f"- refined_collision_graph_count={audit['refined_collision_graph_count']}",
        f"- refined_max_class_size={audit['refined_max_class_size']}",
        f"- remaining_defect_canonical_minus_refined={audit['remaining_defect_canonical_minus_refined']}",
        f"- defect_reduction_vs_spectral_pair_only={audit['defect_reduction_vs_spectral_pair_only']}",
        "",
        "## Decision",
        f"- accepted_as_local_motif_refinement_audit={acceptance['accepted_as_local_motif_refinement_audit']}",
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
    p2810 = read_json(P2810)
    audit = build_audit(p2810)
    payload: dict[str, Any] = {
        "status": "P2811_LOCAL_MOTIF_REFINED_SOURCE_CANDIDATE_AUDIT_REFINES_BUT_STILL_OBSTRUCTED_NO_CLOSURE",
        "input_hashes": {"P2810": sha(P2810), "16_4_4.scd": sha(SCD)},
        "input_statuses": {"P2810": p2810.get("status")},
        "local_motif_refined_audit": audit,
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
            "reason": "P2811 adds one richer local motif/distance profile beyond the P2810 spectral-pair data.  The refinement improves the quotient from 16,211 to 16,691 classes and reduces the defect from 617 to 137, but 132 collision classes covering 269 graphs remain.  Therefore this richer invariant is useful evidence but still not a complete canonical source carrier, strict spectral source law, or K/L_total coupling theorem.",
            "next_honest_step": "Attack the remaining 132 refined collision classes directly with one additional typed invariant, preferably an automorphism/orbit-aware motif invariant or a canonical-label-derived local orbit signature, and report whether it separates the 269 remaining graphs.  Do not promote local motif refinement to K/L_total or ToE while any refined collisions remain and no typed variational coupling is exported.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2811/S1761 local motif-refined source candidate audit", "## P2811/S1761 local motif-refined source candidate audit\n\n`P2811/S1761` adds one richer source-neutral invariant package beyond P2810: exact nonedge common-neighbor, 4-cycle participation, and distance-profile data combined with the P2804 spectral pair.  This improves the quotient from `16,211` to `16,691` classes and reduces the defect from `617` to `137`, but `132` collision classes covering `269` graphs remain.  Thus local motif refinement is evidence, not a complete canonical source carrier, strict spectral source law, or `K`/`L_total` coupling.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2811/S1761 local motif source Ltotal guard", "## P2811/S1761 local motif source Ltotal guard\n\n`P2811/S1761` adds no variational term to `L_total`.  The exact local motif/distance refinement improves the graph-source carrier candidate but still leaves collisions below the P2808 canonical quotient and exports no typed variational coupling.  Future Lagrangian work must first separate the remaining refined collisions or supply an explicit graph-to-`K`/`L_total` source theorem.\n")
    append_once(AGENTS, "Current local motif-refined source candidate guardrail (P2811/S1761, 2026-06-16)", "## Current local motif-refined source candidate guardrail (P2811/S1761, 2026-06-16)\n\n- P2811 adds exact local motif/distance profiles to the P2810 spectral-pair data, improving the quotient from `16,211` to `16,691` classes and reducing the defect from `617` to `137`.\n- The refinement still leaves `132` collision classes covering `269` graphs, so it does not export a complete canonical source carrier, strict spectral source law, `K`/`L_total` coupling, selector closure, bridge closure, role transfer, or ToE closure.\n- A next admissible move should target the remaining refined collisions with exactly one additional typed invariant, preferably automorphism/orbit-aware or canonical-label-derived, and keep all closure claims blocked until separation plus typed variational coupling are actually exported.\n")
    return payload


if __name__ == "__main__":
    main()
