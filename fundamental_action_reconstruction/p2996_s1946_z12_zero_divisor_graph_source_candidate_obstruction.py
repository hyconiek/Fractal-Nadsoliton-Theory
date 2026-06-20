#!/usr/bin/env python3
"""P2996/S1946: Z12 zero-divisor graph source-candidate obstruction.

P2995 bounded the Fourier-character lane.  This audit introduces one genuinely
new finite typed object outside the closed Fourier, annihilator, nilradical,
CRT, and zero-derivation lanes: the zero-divisor graph of Z/12Z.

The finite calculation enumerates all nonzero zero divisors, builds the exact
zero-product adjacency relation, computes degrees, components, triangle count,
and the full graph automorphism group by exhaustive permutation.  This gives a
real multiplicative-incidence object, but it remains a source-candidate
obstruction: no current artifact exports strict nadsoliton provenance,
nonpremise vertex/edge localizer, named source-atom coupling, unit-bearing
action installation, or nonproxy continuum lift for this graph.
"""
from __future__ import annotations

import hashlib, json, math
from itertools import combinations, permutations, product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2995_s1945_fourier_character_action_installation_obstruction import OUT as P2995

OUT = GEN / "p2996_s1946_z12_zero_divisor_graph_source_candidate_obstruction.json"
MD = GEN / "p2996_s1946_z12_zero_divisor_graph_source_candidate_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
MODULUS = 12


def zero_divisor_vertices() -> list[int]:
    return [x for x in range(1, MODULUS) if math.gcd(x, MODULUS) != 1]


def edge_list(vertices: list[int]) -> list[tuple[int, int]]:
    return [(a, b) for a, b in combinations(vertices, 2) if (a * b) % MODULUS == 0]


def connected_components(vertices: list[int], edges: list[tuple[int, int]]) -> list[list[int]]:
    adj = {v: set() for v in vertices}
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    seen: set[int] = set()
    comps: list[list[int]] = []
    for v in vertices:
        if v in seen:
            continue
        stack = [v]
        seen.add(v)
        comp = []
        while stack:
            u = stack.pop()
            comp.append(u)
            for w in sorted(adj[u]):
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        comps.append(sorted(comp))
    return comps


def automorphisms(vertices: list[int], edges: list[tuple[int, int]]) -> list[dict[str, int]]:
    edge_set = {tuple(sorted(e)) for e in edges}
    autos = []
    for perm in permutations(vertices):
        mp = dict(zip(vertices, perm))
        mapped = {tuple(sorted((mp[a], mp[b]))) for a, b in edge_set}
        if mapped == edge_set:
            autos.append({str(k): v for k, v in mp.items()})
    return autos


def zero_divisor_graph_witness() -> dict[str, Any]:
    vertices = zero_divisor_vertices()
    edges = edge_list(vertices)
    degrees = {str(v): 0 for v in vertices}
    for a, b in edges:
        degrees[str(a)] += 1
        degrees[str(b)] += 1
    adjacency_rows = []
    for a in vertices:
        adjacency_rows.append({
            "vertex": a,
            "gcd_with_12": math.gcd(a, MODULUS),
            "neighbors": [b for b in vertices if b != a and tuple(sorted((a, b))) in {tuple(sorted(e)) for e in edges}],
            "degree": degrees[str(a)],
            "strict_nadsoliton_provenance": False,
            "nonpremise_vertex_edge_localizer": False,
            "named_source_atom_coupling": False,
            "unit_bearing_action_installation": False,
            "accepted_strict_source": False,
        })
    triangles = [list(t) for t in combinations(vertices, 3) if all(tuple(sorted((a, b))) in {tuple(sorted(e)) for e in edges} for a, b in combinations(t, 2))]
    autos = automorphisms(vertices, edges)
    return {
        "modulus": MODULUS,
        "vertices": vertices,
        "vertex_count": len(vertices),
        "edges": [list(e) for e in edges],
        "edge_count": len(edges),
        "degree_sequence": sorted(degrees.values()),
        "degrees": degrees,
        "adjacency_rows": adjacency_rows,
        "components": connected_components(vertices, edges),
        "component_count": len(connected_components(vertices, edges)),
        "triangles": triangles,
        "triangle_count": len(triangles),
        "automorphism_group_order": len(autos),
        "automorphisms": autos,
        "accepted_strict_sources": [r["vertex"] for r in adjacency_rows if r["accepted_strict_source"]],
    }


def obligation_rows(witness: dict[str, Any]) -> list[dict[str, Any]]:
    rows = witness["adjacency_rows"]
    return [
        {"obligation": "finite_zero_divisor_graph_enumerated", "satisfied": witness["vertex_count"] == 7 and witness["edge_count"] == 8, "evidence": "Z/12Z has seven nonzero zero divisors and eight zero-product graph edges"},
        {"obligation": "exact_graph_invariants_computed", "satisfied": witness["component_count"] == 1 and witness["automorphism_group_order"] == 8, "evidence": f"degree sequence {witness['degree_sequence']}, connected graph, automorphism group order 8"},
        {"obligation": "strict_nadsoliton_provenance", "satisfied": any(r["strict_nadsoliton_provenance"] for r in rows), "evidence": "the graph imports Z/12Z multiplication rather than a strict nadsoliton source map"},
        {"obligation": "nonpremise_vertex_edge_localizer", "satisfied": any(r["nonpremise_vertex_edge_localizer"] for r in rows), "evidence": "degrees and adjacency rows are incidence labels, not nonpremise physical sectors"},
        {"obligation": "named_source_atom_coupling", "satisfied": any(r["named_source_atom_coupling"] for r in rows), "evidence": "no theorem couples graph vertices/edges to selector sign, beta/Z_beta, bridge-source, or action-density atoms"},
        {"obligation": "unit_bearing_action_installation", "satisfied": any(r["unit_bearing_action_installation"] for r in rows), "evidence": "finite graph incidence has no unit-bearing measure or named density theorem"},
        {"obligation": "accepted_current_strict_source", "satisfied": bool(witness["accepted_strict_sources"]), "evidence": "no vertex/edge row satisfies the full strict source profile"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_graph", "graph_invariants", "strict_provenance", "vertex_edge_localizer", "source_atom_coupling", "unit_action_installation", "nonproxy_export"]
    return [{"present": dict(zip(names, bits)), "accepts_zero_divisor_graph_source": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p2995_path: Any) -> dict[str, Any]:
    witness = zero_divisor_graph_witness()
    matrix = acceptance_matrix()
    return {
        "status": "P2996_Z12_ZERO_DIVISOR_GRAPH_SOURCE_CANDIDATE_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P2995": hashlib.sha256(p2995_path.read_bytes()).hexdigest() if p2995_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "Z12ZeroDivisorGraph_SourceCandidateObstructionMatrix",
            "zero_divisor_graph_witness": witness,
            "proof_obligation_rows": obligation_rows(witness),
            "finite_acceptance_matrix": matrix,
        },
        "graph_certificate": {
            "vertex_count": witness["vertex_count"],
            "edge_count": witness["edge_count"],
            "degree_sequence": witness["degree_sequence"],
            "component_count": witness["component_count"],
            "triangle_count": witness["triangle_count"],
            "automorphism_group_order": witness["automorphism_group_order"],
            "accepted_strict_sources": witness["accepted_strict_sources"],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_zero_divisor_graph_source"]),
        },
        "decision": {
            "positive_progress": "P2996 introduces one new finite typed object outside the closed Fourier/annihilator/nilradical/CRT/zero-derivation lanes: the Z/12Z zero-divisor graph.",
            "breakthrough": "Bounded no-go: the graph has seven vertices, eight zero-product edges, degree sequence [1,1,2,2,3,3,4], one component, no triangles, and automorphism group order 8, but no strict provenance, nonpremise vertex/edge localizer, named source-atom coupling, unit-bearing action installation, or nonproxy export is present.",
            "negative_export_flags": {k: False for k in ["zero_divisor_graph_strict_source_exported", "strict_nadsoliton_provenance_exported", "nonpremise_localizer_exported", "source_atom_coupling_exported", "unit_bearing_action_installation_exported", "selector_closure_exported", "bridge_closure_exported", "nonproxy_ltotal_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "A next admissible zero-divisor-graph move may attack exactly one missing theorem for this new object: strict provenance, nonpremise vertex/edge localizer, named source-atom coupling, or unit-bearing action installation.  Otherwise introduce a genuinely new strict typed object/provider outside this graph lane or preserve the P2929-P2996 no-strict-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["graph_certificate"]
    lines = [
        "# P2996/S1946 Z12 zero-divisor graph source-candidate obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Graph certificate",
        f"- vertex count: `{cert['vertex_count']}`",
        f"- edge count: `{cert['edge_count']}`",
        f"- degree sequence: `{cert['degree_sequence']}`",
        f"- component count: `{cert['component_count']}`",
        f"- triangle count: `{cert['triangle_count']}`",
        f"- automorphism group order: `{cert['automorphism_group_order']}`",
        f"- accepted strict sources: `{cert['accepted_strict_sources']}`",
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
    read_json(P2995)
    payload = build_payload(P2995)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2996/S1946 Z12 zero-divisor graph source-candidate obstruction", "## P2996/S1946 Z12 zero-divisor graph source-candidate obstruction\n\n`P2996/S1946` introduces one new finite typed object outside the closed Fourier/annihilator/nilradical/CRT/zero-derivation lanes: the zero-divisor graph of `Z/12Z`.  The finite computation enumerates seven nonzero zero divisors, eight zero-product edges, degree sequence `[1, 1, 2, 2, 3, 3, 4]`, one connected component, zero triangles, and automorphism group order `8`.  This is exact multiplicative-incidence progress, but no strict nadsoliton provenance, nonpremise vertex/edge localizer, named source-atom coupling, unit-bearing action installation, nonproxy `L_total`, bridge closure, role transfer, or ToE is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2996/S1946 zero-divisor graph `L_total` guard", "## P2996/S1946 zero-divisor graph `L_total` guard\n\n`P2996/S1946` adds no zero-divisor graph term to `L_total`.  Zero-product edges, degree sequence, connectedness, and graph automorphisms are finite incidence data only; they do not supply strict field provenance, named unit-bearing density, variational chain, EOM/Hamiltonian term, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current Z12 zero-divisor graph source-candidate obstruction guardrail (P2996/S1946, 2026-06-20)", "## Current Z12 zero-divisor graph source-candidate obstruction guardrail (P2996/S1946, 2026-06-20)\n\n- P2996 introduces one new finite typed object outside closed Fourier/annihilator/nilradical/CRT/zero-derivation replay: the zero-divisor graph of `Z/12Z`.\n- Exact finite positives are real: seven nonzero zero-divisor vertices, eight zero-product edges, degree sequence `[1, 1, 2, 2, 3, 3, 4]`, one connected component, zero triangles, and automorphism group order `8`.\n- The current route is bounded no-go as a strict source candidate: no strict nadsoliton provenance, nonpremise vertex/edge localizer, named source-atom coupling, unit-bearing action installation, or nonproxy export is present.\n- Do not replay Fourier, annihilator, nilradical, CRT, zero-derivation, selector, bridge, or `L_total` lanes through zero-divisor graph incidence.  A next admissible zero-divisor-graph move may attack exactly one missing theorem for this object, or else introduce a genuinely new strict typed object/provider while preserving the P2929-P2996 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
