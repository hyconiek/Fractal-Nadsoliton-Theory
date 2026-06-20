#!/usr/bin/env python3
"""P2997/S1947: zero-divisor graph strict provenance obstruction.

P2996 introduced the Z/12Z zero-divisor graph and left several missing theorem
routes.  This audit attacks exactly one: strict nadsoliton provenance for the
zero-product graph.  It does not replay the nonpremise localizer, named
source-atom coupling, action installation, Fourier, annihilator, nilradical,
CRT, zero-derivation, selector, bridge, role-transfer, or L_total lanes.

The finite side is exact: every edge is a zero product, every nonedge has a
nonzero product, the unit action by U(12) preserves the graph, and only the
trivial translation preserves the vertex set.  The obstruction is provenance:
these facts remain imported Z/12Z multiplication/incidence data, not a strict
source map exporting the graph from the primordial nadsoliton state.
"""
from __future__ import annotations

import hashlib, json
from itertools import combinations, product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2996_s1946_z12_zero_divisor_graph_source_candidate_obstruction import MODULUS, OUT as P2996, zero_divisor_graph_witness

OUT = GEN / "p2997_s1947_zero_divisor_graph_strict_provenance_obstruction.json"
MD = GEN / "p2997_s1947_zero_divisor_graph_strict_provenance_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
UNITS = [1, 5, 7, 11]


def provenance_witness() -> dict[str, Any]:
    graph = zero_divisor_graph_witness()
    vertices = graph["vertices"]
    edge_set = {tuple(e) for e in graph["edges"]}
    all_pairs = {tuple(sorted(p)) for p in combinations(vertices, 2)}
    nonedges = sorted(all_pairs - edge_set)
    edge_rows = [{"a": a, "b": b, "product_mod_12": (a * b) % MODULUS, "zero_product_edge": (a * b) % MODULUS == 0} for a, b in sorted(edge_set)]
    nonedge_rows = [{"a": a, "b": b, "product_mod_12": (a * b) % MODULUS, "nonzero_product_nonedge": (a * b) % MODULUS != 0} for a, b in nonedges]
    unit_action_rows = []
    for u in UNITS:
        mapping = {v: (u * v) % MODULUS for v in vertices}
        mapped_edges = {tuple(sorted((mapping[a], mapping[b]))) for a, b in edge_set}
        unit_action_rows.append({"unit": u, "mapping": {str(k): v for k, v in mapping.items()}, "preserves_vertices": sorted(mapping.values()) == vertices, "preserves_edges": mapped_edges == edge_set})
    translation_rows = []
    for t in range(MODULUS):
        image = sorted((v + t) % MODULUS for v in vertices)
        translation_rows.append({"translation": t, "image": image, "preserves_vertex_set": image == vertices})
    return {
        "modulus": MODULUS,
        "vertex_count": graph["vertex_count"],
        "edge_count": graph["edge_count"],
        "edge_rows": edge_rows,
        "nonedge_rows": nonedge_rows,
        "all_edges_zero_products": all(r["zero_product_edge"] for r in edge_rows),
        "all_nonedges_nonzero_products": all(r["nonzero_product_nonedge"] for r in nonedge_rows),
        "unit_action_rows": unit_action_rows,
        "all_unit_actions_preserve_graph": all(r["preserves_vertices"] and r["preserves_edges"] for r in unit_action_rows),
        "translation_rows": translation_rows,
        "vertex_set_preserving_translations": [r["translation"] for r in translation_rows if r["preserves_vertex_set"]],
        "strict_nadsoliton_source_map_exported": False,
        "nonpremise_internal_graph_provenance": False,
        "accepted_strict_provenance": False,
    }


def obligation_rows(witness: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {"obligation": "finite_zero_product_graph_exact", "satisfied": witness["all_edges_zero_products"] and witness["all_nonedges_nonzero_products"], "evidence": "all 8 edges are zero products and all 13 nonedges are nonzero products"},
        {"obligation": "unit_action_graph_symmetry_verified", "satisfied": witness["all_unit_actions_preserve_graph"], "evidence": "multiplication by all four units of U(12) preserves vertices and edges"},
        {"obligation": "translation_noninvariance_witness", "satisfied": witness["vertex_set_preserving_translations"] == [0], "evidence": "only translation 0 preserves the zero-divisor vertex set"},
        {"obligation": "strict_nadsoliton_source_map_exported", "satisfied": witness["strict_nadsoliton_source_map_exported"], "evidence": "no current artifact exports this graph from the primordial nadsoliton state"},
        {"obligation": "nonpremise_internal_graph_provenance", "satisfied": witness["nonpremise_internal_graph_provenance"], "evidence": "graph rows remain imported Z/12Z multiplication/incidence data"},
        {"obligation": "accepted_current_strict_provenance", "satisfied": witness["accepted_strict_provenance"], "evidence": "the full strict provenance profile is not satisfied"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_graph_exact", "unit_action_symmetry", "translation_boundary", "strict_source_map", "internal_provenance", "nonproxy_export"]
    return [{"present": dict(zip(names, bits)), "accepts_zero_divisor_graph_strict_provenance": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p2996_path: Any) -> dict[str, Any]:
    witness = provenance_witness()
    matrix = acceptance_matrix()
    return {
        "status": "P2997_ZERO_DIVISOR_GRAPH_STRICT_PROVENANCE_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P2996": hashlib.sha256(p2996_path.read_bytes()).hexdigest() if p2996_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "ZeroDivisorGraphStrictProvenance_ObstructionMatrix",
            "provenance_witness": witness,
            "proof_obligation_rows": obligation_rows(witness),
            "finite_acceptance_matrix": matrix,
        },
        "provenance_certificate": {
            "vertex_count": witness["vertex_count"],
            "edge_count": witness["edge_count"],
            "edge_rows": len(witness["edge_rows"]),
            "nonedge_rows": len(witness["nonedge_rows"]),
            "all_edges_zero_products": witness["all_edges_zero_products"],
            "all_nonedges_nonzero_products": witness["all_nonedges_nonzero_products"],
            "all_unit_actions_preserve_graph": witness["all_unit_actions_preserve_graph"],
            "vertex_set_preserving_translations": witness["vertex_set_preserving_translations"],
            "accepted_strict_provenance": witness["accepted_strict_provenance"],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_zero_divisor_graph_strict_provenance"]),
        },
        "decision": {
            "positive_progress": "P2997 attacks exactly one P2996 missing route: strict provenance for the Z/12Z zero-divisor graph.",
            "breakthrough": "Bounded no-go: the zero-product graph is exact and U(12)-symmetric, with translation noninvariance witnessed by only translation 0 preserving the vertex set, but no strict nadsoliton source map or nonpremise internal graph provenance is exported.",
            "negative_export_flags": {k: False for k in ["zero_divisor_graph_strict_provenance_exported", "strict_nadsoliton_source_map_exported", "nonpremise_internal_graph_provenance_exported", "vertex_edge_localizer_exported", "source_atom_coupling_exported", "unit_bearing_action_installation_exported", "bridge_closure_exported", "nonproxy_ltotal_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "A next admissible zero-divisor-graph move may attack exactly one remaining route: nonpremise vertex/edge localizer, named source-atom coupling, or unit-bearing action installation.  Do not replay exact zero-product incidence or U(12) symmetry as provenance; otherwise introduce a genuinely new strict typed object/provider or preserve the P2929-P2997 no-strict-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["provenance_certificate"]
    lines = [
        "# P2997/S1947 zero-divisor graph strict provenance obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Provenance certificate",
        f"- vertices/edges: `{cert['vertex_count']}/{cert['edge_count']}`",
        f"- edge/nonedge rows: `{cert['edge_rows']}/{cert['nonedge_rows']}`",
        f"- all edges zero products: `{cert['all_edges_zero_products']}`",
        f"- all nonedges nonzero products: `{cert['all_nonedges_nonzero_products']}`",
        f"- all unit actions preserve graph: `{cert['all_unit_actions_preserve_graph']}`",
        f"- vertex-set preserving translations: `{cert['vertex_set_preserving_translations']}`",
        f"- accepted strict provenance: `{cert['accepted_strict_provenance']}`",
        f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_rows']}`",
        "",
        "## Lay summary",
        payload["decision"]["positive_progress"],
        payload["decision"]["breakthrough"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P2996)
    payload = build_payload(P2996)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2997/S1947 zero-divisor graph strict provenance obstruction", "## P2997/S1947 zero-divisor graph strict provenance obstruction\n\n`P2997/S1947` attacks exactly one P2996 missing theorem route: strict provenance for the zero-divisor graph of `Z/12Z`.  The finite side is exact: all `8` graph edges are zero products, all `13` nonedges are nonzero products, multiplication by each unit in `U(12)={1,5,7,11}` preserves the graph, and only translation `0` preserves the vertex set.  The provenance side remains blocked because these are imported `Z/12Z` multiplication/incidence facts, not a strict nadsoliton source map or nonpremise internal graph provenance.  No vertex/edge localizer, source-atom coupling, unit-bearing action installation, nonproxy `L_total`, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2997/S1947 zero-divisor graph provenance `L_total` guard", "## P2997/S1947 zero-divisor graph provenance `L_total` guard\n\n`P2997/S1947` adds no zero-divisor graph provenance term to `L_total`.  Exact zero-product incidence, nonedge products, unit-action symmetry, and translation noninvariance are finite provenance tests only; they do not supply strict field provenance, named unit-bearing density, variational chain, EOM/Hamiltonian term, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current zero-divisor graph strict provenance obstruction guardrail (P2997/S1947, 2026-06-20)", "## Current zero-divisor graph strict provenance obstruction guardrail (P2997/S1947, 2026-06-20)\n\n- P2997 attacks exactly one P2996 missing theorem route: strict provenance for the `Z/12Z` zero-divisor graph.\n- Finite positives are exact but incidence-only: all `8` edges are zero products, all `13` nonedges are nonzero products, all four `U(12)` unit actions preserve the graph, and only translation `0` preserves the vertex set.\n- The route is bounded no-go because these facts remain imported `Z/12Z` multiplication/incidence data; no strict nadsoliton source map or nonpremise internal graph provenance is exported.\n- Do not promote zero-product incidence, nonedge products, unit-action symmetry, or translation noninvariance to vertex/edge localizer, source-atom coupling, unit-bearing action installation, bridge closure, role transfer, nonproxy `L_total`, or ToE.  A next admissible zero-divisor-graph move may attack exactly one remaining route, or else introduce a genuinely new strict typed object/provider while preserving the P2929-P2997 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
