#!/usr/bin/env python3
"""P2998/S1948: zero-divisor graph nonpremise vertex/edge localizer obstruction.

P2997 left three zero-divisor-graph routes.  This audit attacks exactly one:
a nonpremise vertex/edge localizer for the Z/12Z zero-divisor graph.  It does
not replay strict provenance, named source-atom coupling, action installation,
Fourier, annihilator, nilradical, CRT, zero-derivation, selector, bridge,
role-transfer, or L_total lanes.

The finite calculation builds degree, neighbor-degree, closed-neighborhood, and
automorphism-orbit signatures for every vertex and edge.  These signatures are
exact and useful, but they are graph-bookkeeping labels.  They identify one
central vertex orbit and classify edge types, yet they do not export a
nonpremise physical sector or strict nadsoliton vertex/edge localizer.
"""
from __future__ import annotations

import hashlib, json
from itertools import combinations, product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2996_s1946_z12_zero_divisor_graph_source_candidate_obstruction import OUT as P2996, zero_divisor_graph_witness
from p2997_s1947_zero_divisor_graph_strict_provenance_obstruction import OUT as P2997

OUT = GEN / "p2998_s1948_zero_divisor_graph_nonpremise_vertex_edge_localizer_obstruction.json"
MD = GEN / "p2998_s1948_zero_divisor_graph_nonpremise_vertex_edge_localizer_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def graph_data() -> tuple[list[int], list[tuple[int, int]], dict[int, set[int]]]:
    graph = zero_divisor_graph_witness()
    vertices = graph["vertices"]
    edges = [tuple(e) for e in graph["edges"]]
    adj = {v: set() for v in vertices}
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    return vertices, edges, adj


def orbit_id(value: int, automorphisms: list[dict[str, int]]) -> tuple[int, ...]:
    return tuple(sorted({auto[str(value)] for auto in automorphisms}))


def edge_orbit_id(edge: tuple[int, int], automorphisms: list[dict[str, int]]) -> tuple[tuple[int, int], ...]:
    a, b = edge
    return tuple(sorted({tuple(sorted((auto[str(a)], auto[str(b)]))) for auto in automorphisms}))


def localizer_witness() -> dict[str, Any]:
    graph = zero_divisor_graph_witness()
    vertices, edges, adj = graph_data()
    degrees = {v: len(adj[v]) for v in vertices}
    autos = graph["automorphisms"]
    vertex_rows = []
    for v in vertices:
        neighbor_degrees = sorted(degrees[n] for n in adj[v])
        closed_neighborhood = sorted([v] + list(adj[v]))
        closed_degree_signature = sorted(degrees[n] for n in closed_neighborhood)
        vertex_rows.append({
            "vertex": v,
            "degree": degrees[v],
            "neighbors": sorted(adj[v]),
            "neighbor_degree_signature": neighbor_degrees,
            "closed_neighborhood": closed_neighborhood,
            "closed_degree_signature": closed_degree_signature,
            "automorphism_orbit": list(orbit_id(v, autos)),
            "orbit_size": len(orbit_id(v, autos)),
            "singleton_graph_signature": len(orbit_id(v, autos)) == 1,
            "strict_provenance_available": False,
            "nonpremise_physical_sector": False,
            "accepted_vertex_localizer": False,
        })
    edge_rows = []
    for edge in edges:
        a, b = edge
        shared_neighbors = sorted(adj[a] & adj[b])
        orbit = edge_orbit_id(edge, autos)
        edge_rows.append({
            "edge": [a, b],
            "endpoint_degrees": sorted([degrees[a], degrees[b]]),
            "shared_neighbor_count": len(shared_neighbors),
            "shared_neighbors": shared_neighbors,
            "automorphism_orbit": [list(e) for e in orbit],
            "orbit_size": len(orbit),
            "singleton_graph_signature": len(orbit) == 1,
            "strict_provenance_available": False,
            "nonpremise_physical_sector": False,
            "accepted_edge_localizer": False,
        })
    vertex_orbits = sorted({tuple(r["automorphism_orbit"]) for r in vertex_rows}, key=lambda x: (len(x), x))
    edge_orbits = sorted({tuple(tuple(e) for e in r["automorphism_orbit"]) for r in edge_rows}, key=lambda x: (len(x), x))
    return {
        "vertex_count": len(vertices),
        "edge_count": len(edges),
        "vertex_rows": vertex_rows,
        "edge_rows": edge_rows,
        "vertex_orbits": [list(o) for o in vertex_orbits],
        "edge_orbits": [[list(e) for e in orbit] for orbit in edge_orbits],
        "vertex_orbit_count": len(vertex_orbits),
        "edge_orbit_count": len(edge_orbits),
        "singleton_vertex_orbits": [r["vertex"] for r in vertex_rows if r["singleton_graph_signature"]],
        "singleton_edge_orbits": [r["edge"] for r in edge_rows if r["singleton_graph_signature"]],
        "accepted_vertex_localizers": [r["vertex"] for r in vertex_rows if r["accepted_vertex_localizer"]],
        "accepted_edge_localizers": [r["edge"] for r in edge_rows if r["accepted_edge_localizer"]],
    }


def obligation_rows(witness: dict[str, Any]) -> list[dict[str, Any]]:
    vertex_rows = witness["vertex_rows"]
    edge_rows = witness["edge_rows"]
    return [
        {"obligation": "finite_vertex_edge_signatures_built", "satisfied": witness["vertex_count"] == 7 and witness["edge_count"] == 8, "evidence": "degree, neighborhood, and automorphism-orbit signatures are built for all vertices and edges"},
        {"obligation": "automorphism_orbit_classification", "satisfied": witness["vertex_orbit_count"] == 4 and witness["edge_orbit_count"] >= 1, "evidence": f"vertex orbits {witness['vertex_orbits']} and {witness['edge_orbit_count']} edge orbits are computed"},
        {"obligation": "strict_provenance_available", "satisfied": any(r["strict_provenance_available"] for r in vertex_rows + edge_rows), "evidence": "P2997 found no strict graph provenance"},
        {"obligation": "nonpremise_physical_sector", "satisfied": any(r["nonpremise_physical_sector"] for r in vertex_rows + edge_rows), "evidence": "graph signatures are incidence labels, not physical sectors"},
        {"obligation": "accepted_vertex_localizer", "satisfied": bool(witness["accepted_vertex_localizers"]), "evidence": "even singleton graph orbit data lacks strict provenance and physical-sector theorem"},
        {"obligation": "accepted_edge_localizer", "satisfied": bool(witness["accepted_edge_localizers"]), "evidence": "edge orbit/type data lacks strict provenance and physical-sector theorem"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_signatures", "orbit_classification", "strict_provenance", "physical_sector", "localizer_theorem", "nonproxy_export"]
    return [{"present": dict(zip(names, bits)), "accepts_nonpremise_vertex_edge_localizer": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p2997_path: Any) -> dict[str, Any]:
    witness = localizer_witness()
    matrix = acceptance_matrix()
    return {
        "status": "P2998_ZERO_DIVISOR_GRAPH_NONPREMISE_VERTEX_EDGE_LOCALIZER_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {
            "P2996": hashlib.sha256(P2996.read_bytes()).hexdigest() if P2996.exists() else None,
            "P2997": hashlib.sha256(p2997_path.read_bytes()).hexdigest() if p2997_path.exists() else None,
        },
        "constructed_theoretical_objects": {
            "object": "ZeroDivisorGraphNonpremiseVertexEdgeLocalizer_ObstructionMatrix",
            "localizer_witness": witness,
            "proof_obligation_rows": obligation_rows(witness),
            "finite_acceptance_matrix": matrix,
        },
        "localizer_certificate": {
            "vertex_count": witness["vertex_count"],
            "edge_count": witness["edge_count"],
            "vertex_orbit_count": witness["vertex_orbit_count"],
            "edge_orbit_count": witness["edge_orbit_count"],
            "singleton_vertex_orbits": witness["singleton_vertex_orbits"],
            "singleton_edge_orbits": witness["singleton_edge_orbits"],
            "accepted_vertex_localizers": witness["accepted_vertex_localizers"],
            "accepted_edge_localizers": witness["accepted_edge_localizers"],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_nonpremise_vertex_edge_localizer"]),
        },
        "decision": {
            "positive_progress": "P2998 attacks exactly one P2997 remaining route: nonpremise vertex/edge localizer for the Z/12Z zero-divisor graph.",
            "breakthrough": "Bounded no-go: finite degree/neighborhood/orbit signatures classify the graph into four vertex orbits and computed edge orbits, with vertex 6 a singleton graph orbit, but no strict provenance or nonpremise physical-sector theorem turns these labels into accepted localizers.",
            "negative_export_flags": {k: False for k in ["zero_divisor_graph_localizer_exported", "strict_provenance_exported", "nonpremise_physical_sector_exported", "source_atom_coupling_exported", "unit_bearing_action_installation_exported", "bridge_closure_exported", "nonproxy_ltotal_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "A next admissible zero-divisor-graph move may attack exactly one remaining route: named source-atom coupling or unit-bearing action installation.  Do not replay degree/neighborhood/orbit signatures as a nonpremise localizer; otherwise introduce a genuinely new strict typed object/provider or preserve the P2929-P2998 no-strict-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["localizer_certificate"]
    lines = [
        "# P2998/S1948 zero-divisor graph nonpremise vertex/edge localizer obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Localizer certificate",
        f"- vertices/edges: `{cert['vertex_count']}/{cert['edge_count']}`",
        f"- vertex orbit count: `{cert['vertex_orbit_count']}`",
        f"- edge orbit count: `{cert['edge_orbit_count']}`",
        f"- singleton vertex orbits: `{cert['singleton_vertex_orbits']}`",
        f"- singleton edge orbits: `{cert['singleton_edge_orbits']}`",
        f"- accepted vertex localizers: `{cert['accepted_vertex_localizers']}`",
        f"- accepted edge localizers: `{cert['accepted_edge_localizers']}`",
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
    read_json(P2997)
    payload = build_payload(P2997)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2998/S1948 zero-divisor graph nonpremise vertex/edge localizer obstruction", "## P2998/S1948 zero-divisor graph nonpremise vertex/edge localizer obstruction\n\n`P2998/S1948` attacks exactly one P2997 remaining route: a nonpremise vertex/edge localizer for the `Z/12Z` zero-divisor graph.  The finite side builds exact degree, neighbor-degree, closed-neighborhood, and automorphism-orbit signatures for all vertices and edges.  These signatures classify the graph into four vertex orbits and computed edge orbits, with vertex `6` a singleton graph orbit.  The localizer side remains blocked because P2997 exported no strict graph provenance and the signatures are incidence labels, not nonpremise physical sectors.  No accepted vertex/edge localizer, source-atom coupling, unit-bearing action installation, nonproxy `L_total`, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2998/S1948 zero-divisor graph localizer `L_total` guard", "## P2998/S1948 zero-divisor graph localizer `L_total` guard\n\n`P2998/S1948` adds no zero-divisor graph localizer term to `L_total`.  Degree, neighborhood, and orbit signatures are finite incidence labels only; they do not supply strict field provenance, named unit-bearing density, variational chain, EOM/Hamiltonian term, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current zero-divisor graph nonpremise vertex/edge localizer obstruction guardrail (P2998/S1948, 2026-06-20)", "## Current zero-divisor graph nonpremise vertex/edge localizer obstruction guardrail (P2998/S1948, 2026-06-20)\n\n- P2998 attacks exactly one P2997 remaining route: nonpremise vertex/edge localizer for the `Z/12Z` zero-divisor graph.\n- Finite positives are graph-bookkeeping only: degree, neighbor-degree, closed-neighborhood, and automorphism-orbit signatures classify four vertex orbits and computed edge orbits, with vertex `6` a singleton graph orbit.\n- The route is bounded no-go because P2997 exported no strict graph provenance and these signatures are incidence labels, not nonpremise physical sectors.\n- Do not promote degree/neighborhood/orbit signatures to source-atom coupling, unit-bearing action installation, bridge closure, role transfer, nonproxy `L_total`, or ToE.  A next admissible zero-divisor-graph move may attack exactly one remaining route, or else introduce a genuinely new strict typed object/provider while preserving the P2929-P2998 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
