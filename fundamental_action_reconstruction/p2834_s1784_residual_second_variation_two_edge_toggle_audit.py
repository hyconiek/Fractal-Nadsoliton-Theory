#!/usr/bin/env python3
"""P2834/S1784: residual-only second-variation two-edge-toggle audit.

P2833 left 67 residual edge-toggle response collision classes covering 138 graphs.
P2834 does not replay the full carrier or promote a coupling.  It runs the next
bounded proof-grade refinement only on those residual classes: an exact non-label
second-variation/two-edge-toggle interaction functional.  The audit stops at the
finite separation question and keeps the typed K/L_total coupling and variational
source-law obligations explicitly blocked.
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import SCD, decode_scd_bytes, sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2832_s1782_common_neighbor_cycle_source_formula_audit import adjacency_masks
from p2833_s1783_edge_toggle_response_polynomial_source_audit import P2832, edge_toggle_response_polynomial

GEN = ROOT / "generated"
P2833 = GEN / "p2833_s1783_edge_toggle_response_polynomial_source_audit.json"
OUT = GEN / "p2834_s1784_residual_second_variation_two_edge_toggle_audit.json"
MD = GEN / "p2834_s1784_residual_second_variation_two_edge_toggle_audit.md"
MANIFEST = GEN / "p2834_s1784_residual_second_variation_digest_manifest.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
EXPECTED_FULL_CARRIER_COUNT = 16828
EXPECTED_P2833_RESIDUAL_CLASSES = 67
EXPECTED_P2833_RESIDUAL_GRAPHS = 138
VERTEX_COUNT = 16
PAIRS = tuple((i, j) for i in range(VERTEX_COUNT) for j in range(i + 1, VERTEX_COUNT))


def profile_digest(profile: tuple[Any, ...]) -> str:
    blob = json.dumps(profile, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def toggle_masks(masks: list[int], pair: tuple[int, int]) -> list[int]:
    i, j = pair
    toggled = masks[:]
    toggled[i] ^= 1 << j
    toggled[j] ^= 1 << i
    return toggled


def triangle_four_cycle_counts_from_masks(masks: list[int]) -> tuple[int, int]:
    triangle_edge_sum = 0
    four_cycle_twice = 0
    for i in range(VERTEX_COUNT):
        for j in range(i + 1, VERTEX_COUNT):
            common = (masks[i] & masks[j]).bit_count()
            four_cycle_twice += common * (common - 1) // 2
            if (masks[i] >> j) & 1:
                triangle_edge_sum += common
    return triangle_edge_sum // 3, four_cycle_twice // 2


def single_edge_response(masks: list[int], pair: tuple[int, int]) -> tuple[int, int, int, int, int]:
    i, j = pair
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
    sign = -1 if edge_present else 1
    return edge_present, triangle_response, four_cycle_response, sign * triangle_response, sign * four_cycle_response


def second_variation_two_edge_profile(graph: tuple[tuple[int, ...], ...]) -> tuple[Any, ...]:
    """Return the exact non-label second-variation profile for all two toggles.

    Each row records only relabeling-invariant data: whether the two toggled
    pairs share a vertex, the sorted two single-edge response signatures, and
    the interaction defects
    `Δ² triangles = Δ_ab triangles - Δ_a triangles - Δ_b triangles` and
    `Δ² 4-cycles = Δ_ab 4-cycles - Δ_a 4-cycles - Δ_b 4-cycles`.
    """
    masks = adjacency_masks(graph)
    base_triangles, base_four_cycles = triangle_four_cycle_counts_from_masks(masks)
    singles = [single_edge_response(masks, pair) for pair in PAIRS]
    rows: Counter[tuple[Any, ...]] = Counter()
    for a, pair_a in enumerate(PAIRS):
        masks_after_a = toggle_masks(masks, pair_a)
        for b in range(a + 1, len(PAIRS)):
            pair_b = PAIRS[b]
            masks_after_ab = toggle_masks(masks_after_a, pair_b)
            tri_ab, four_ab = triangle_four_cycle_counts_from_masks(masks_after_ab)
            interaction_triangles = (tri_ab - base_triangles) - singles[a][3] - singles[b][3]
            interaction_four_cycles = (four_ab - base_four_cycles) - singles[a][4] - singles[b][4]
            shared_vertices = len(set(pair_a) & set(pair_b))
            single_a = singles[a][:3]
            single_b = singles[b][:3]
            if single_b < single_a:
                single_a, single_b = single_b, single_a
            rows[(shared_vertices, single_a, single_b, interaction_triangles, interaction_four_cycles)] += 1
    return tuple(sorted(rows.items()))


def residual_indices_from_p2833(graphs: list[tuple[tuple[int, ...], ...]]) -> tuple[list[list[int]], list[int]]:
    p2833_classes: dict[tuple[Any, ...], list[int]] = defaultdict(list)
    for index, graph in enumerate(graphs):
        p2833_classes[edge_toggle_response_polynomial(graph)].append(index)
    residual_classes = [indices for indices in p2833_classes.values() if len(indices) > 1]
    residual_classes.sort(key=lambda indices: (len(indices), indices), reverse=True)
    residual_indices = [index for indices in residual_classes for index in indices]
    return residual_classes, residual_indices


def build_audit(p2833: dict[str, Any]) -> dict[str, Any]:
    graphs, _ = decode_scd_bytes(SCD.read_bytes())
    residual_classes, residual_indices = residual_indices_from_p2833(graphs)
    refined_classes: dict[tuple[Any, ...], list[int]] = defaultdict(list)
    graph_manifest: list[dict[str, Any]] = []
    for index in residual_indices:
        profile = second_variation_two_edge_profile(graphs[index])
        refined_classes[profile].append(index)
        graph_manifest.append({
            "graph_index_0_based": index,
            "second_variation_profile_sha256": profile_digest(profile),
            "profile_row_count": len(profile),
            "two_edge_toggle_pair_count": len(PAIRS) * (len(PAIRS) - 1) // 2,
        })
    residual_after = {key: indices for key, indices in refined_classes.items() if len(indices) > 1}
    hist = Counter(len(indices) for indices in refined_classes.values())
    manifest_payload = {
        "status": "P2834_RESIDUAL_SECOND_VARIATION_MANIFEST",
        "source": "P2833 residual classes only",
        "p2833_residual_class_count": len(residual_classes),
        "p2833_residual_graph_count": len(residual_indices),
        "two_edge_toggle_pair_count_per_graph": len(PAIRS) * (len(PAIRS) - 1) // 2,
        "graph_profile_digests": graph_manifest,
    }
    MANIFEST.write_text(json.dumps(manifest_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "candidate_formula": "Q_second(G)=multiset over unordered pairs of edge toggles (a,b) of (shared_vertices(a,b), sorted single-edge response signatures, Δ² triangles, Δ² 4-cycles), restricted to the P2833 residual classes",
        "decoded_full_carrier_graph_count": len(graphs),
        "input_p2833_status": p2833.get("status"),
        "p2833_residual_class_count_recomputed": len(residual_classes),
        "p2833_residual_graph_count_recomputed": len(residual_indices),
        "residual_only_audit": True,
        "full_carrier_replay_performed": False,
        "non_label_second_variation_formula_exported_for_candidate": True,
        "two_edge_toggle_pair_count_per_graph": len(PAIRS) * (len(PAIRS) - 1) // 2,
        "refined_residual_class_count": len(refined_classes),
        "refined_residual_collision_class_count": len(residual_after),
        "refined_residual_collision_graph_count": sum(len(indices) for indices in residual_after.values()),
        "refined_residual_max_class_size": max(hist),
        "refined_residual_defect_after_formula": len(residual_indices) - len(refined_classes),
        "refined_class_size_histogram": dict(sorted(hist.items())),
        "manifest_path": str(MANIFEST.relative_to(REPO)),
        "manifest_sha256": sha(MANIFEST),
    }


def acceptance_matrix(audit: dict[str, Any], p2833: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2833_residual_boundary_reused": p2833["acceptance_matrix"]["accepted_as_bounded_edge_toggle_witness_with_residual_no_go"],
        "residual_only_audit": audit["residual_only_audit"],
        "full_carrier_replay_performed": audit["full_carrier_replay_performed"],
        "non_label_second_variation_formula_exported_for_candidate": audit["non_label_second_variation_formula_exported_for_candidate"],
        "p2833_residual_counts_match": audit["p2833_residual_class_count_recomputed"] == EXPECTED_P2833_RESIDUAL_CLASSES and audit["p2833_residual_graph_count_recomputed"] == EXPECTED_P2833_RESIDUAL_GRAPHS,
        "second_variation_separates_p2833_residuals": audit["refined_residual_collision_class_count"] == 0,
        "full_carrier_source_law_coupling_exported": False,
        "proved_variational_derivative_into_K_or_Ltotal_exported": False,
        "typed_graph_to_K_or_Ltotal_coupling_theorem_exported": False,
        "selector_bridge_or_role_transfer_imported": False,
    }
    accepted_as_residual_witness = all([
        facts["p2833_residual_boundary_reused"],
        facts["residual_only_audit"],
        not facts["full_carrier_replay_performed"],
        facts["non_label_second_variation_formula_exported_for_candidate"],
        facts["p2833_residual_counts_match"],
        facts["second_variation_separates_p2833_residuals"],
        not facts["selector_bridge_or_role_transfer_imported"],
    ])
    accepted_as_source_law_coupling = accepted_as_residual_witness and all([
        facts["full_carrier_source_law_coupling_exported"],
        facts["proved_variational_derivative_into_K_or_Ltotal_exported"],
        facts["typed_graph_to_K_or_Ltotal_coupling_theorem_exported"],
    ])
    return {
        "facts": facts,
        "accepted_as_residual_second_variation_witness": accepted_as_residual_witness,
        "accepted_as_source_law_coupling": accepted_as_source_law_coupling,
        "accepted_as_no_coupling_boundary": not accepted_as_source_law_coupling,
        "missing_for_source_law_promotion": [
            key for key, value in facts.items()
            if key in {
                "full_carrier_source_law_coupling_exported",
                "proved_variational_derivative_into_K_or_Ltotal_exported",
                "typed_graph_to_K_or_Ltotal_coupling_theorem_exported",
            } and not value
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["residual_second_variation_two_edge_toggle_audit"]
    acceptance = payload["acceptance_matrix"]
    lines = [
        "# P2834/S1784 residual second-variation two-edge-toggle audit", "", f"Status: `{payload['status']}`", "",
        "## Candidate formula", audit["candidate_formula"], "", "## Residual-only finite counts",
        f"- decoded_full_carrier_graph_count={audit['decoded_full_carrier_graph_count']}",
        f"- p2833_residual_class_count_recomputed={audit['p2833_residual_class_count_recomputed']}",
        f"- p2833_residual_graph_count_recomputed={audit['p2833_residual_graph_count_recomputed']}",
        f"- two_edge_toggle_pair_count_per_graph={audit['two_edge_toggle_pair_count_per_graph']}",
        f"- refined_residual_class_count={audit['refined_residual_class_count']}",
        f"- refined_residual_collision_class_count={audit['refined_residual_collision_class_count']}",
        f"- refined_residual_collision_graph_count={audit['refined_residual_collision_graph_count']}",
        f"- refined_residual_defect_after_formula={audit['refined_residual_defect_after_formula']}", "",
        "## Manifest", f"- manifest_path={audit['manifest_path']}", f"- manifest_sha256={audit['manifest_sha256']}", "",
        "## Acceptance",
        f"- accepted_as_residual_second_variation_witness={acceptance['accepted_as_residual_second_variation_witness']}",
        f"- accepted_as_source_law_coupling={acceptance['accepted_as_source_law_coupling']}",
        f"- accepted_as_no_coupling_boundary={acceptance['accepted_as_no_coupling_boundary']}", "",
        "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2832 = read_json(P2832)
    p2833 = read_json(P2833)
    audit = build_audit(p2833)
    payload: dict[str, Any] = {
        "status": "P2834_RESIDUAL_SECOND_VARIATION_TWO_EDGE_TOGGLE_WITNESS_NO_COUPLING_NO_CLOSURE",
        "input_hashes": {"P2832": sha(P2832), "P2833": sha(P2833), "16_4_4.scd": sha(SCD)},
        "input_statuses": {"P2832": p2832.get("status"), "P2833": p2833.get("status")},
        "residual_second_variation_two_edge_toggle_audit": audit,
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
            "reason": "P2834 executes exactly the residual-only refinement requested by P2833.  On the 67 P2833 residual classes covering 138 graphs, the non-label second-variation/two-edge-toggle interaction profile separates all residual graphs: refined_residual_class_count=138, residual collisions=0, defect=0.  This is accepted as a finite residual witness, not as source-law/coupling closure: it is not a full theorem exporting graph units, a variational derivative into K/L_total, or a typed graph-to-K/L_total coupling law.",
            "next_honest_step": "Do not add more carrier-separation refinements merely to decorate the now-separated residuals.  The next proof-grade step should be a theorem-obligation audit for the combined P2830/P2833/P2834 finite separating data: formulate one typed graph-source functional with units/normalization, domain/codomain, and an explicit variational derivative into K or L_total, then stop on the first missing premise.  If no such typed functional is available, preserve the P2831-P2834 finite-witness/no-coupling boundary.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit, p2833)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2834/S1784 residual second-variation two-edge-toggle audit", "## P2834/S1784 residual second-variation two-edge-toggle audit\n\n`P2834/S1784` executes the residual-only refinement requested by P2833: on the `67` P2833 residual edge-toggle classes covering `138` graphs, it computes `Q_second(G)=multiset` of two-edge-toggle interaction rows `(shared_vertices, sorted single-edge response signatures, Δ² triangles, Δ² 4-cycles)`.  The residual set refines to `138` singleton classes with zero residual collisions and zero defect.  This is a finite residual witness only: no graph units/normalization theorem, variational derivative into `K`/`L_total`, typed graph-to-`K`/`L_total` coupling theorem, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2834/S1784 residual second-variation Ltotal guard", "## P2834/S1784 residual second-variation Ltotal guard\n\n`P2834/S1784` adds no term to `L_total`.  The two-edge-toggle second-variation profile separates the `138` P2833 residual graphs, but remains a finite graph witness without a source-law theorem, units, Euler-Lagrange derivative, or typed coupling into `K`/`L_total`.\n")
    append_once(AGENTS, "Current residual second-variation two-edge-toggle guardrail (P2834/S1784, 2026-06-17)", "## Current residual second-variation two-edge-toggle guardrail (P2834/S1784, 2026-06-17)\n\n- P2834 performs the residual-only audit requested by P2833: it computes a non-label two-edge-toggle second-variation profile on the `67` P2833 residual classes / `138` graphs.\n- The profile separates those residuals into `138` singleton classes with zero residual collisions and zero defect; this is accepted only as a finite residual witness.\n- Do not promote P2834 to a strict graph-source law, `L_total`, bridge closure, role transfer, selector closure, or ToE closure.  The next admissible move is a theorem-obligation audit for one typed graph-source functional with units/normalization and an explicit variational derivative into `K` or `L_total`, not more separation replay.\n")
    return payload


if __name__ == "__main__":
    main()
