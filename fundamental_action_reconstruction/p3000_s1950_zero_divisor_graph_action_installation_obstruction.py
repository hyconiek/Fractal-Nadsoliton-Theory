#!/usr/bin/env python3
"""P3000/S1950: zero-divisor graph action/variational installation obstruction.

P2999 left the final zero-divisor-graph route: unit-bearing action installation
with a genuinely unit-bearing graph measure, named graph density theorem,
boundary/integration map, and nonproxy continuum lift.  This audit attacks that
route only.  It does not replay strict provenance, nonpremise localizer,
source-atom coupling, Fourier, annihilator, nilradical, CRT, zero-derivation,
selector, bridge completion, role transfer, or L_total promotion.

The finite calculation constructs the strongest formal graph-action receivers
available from the Z/12Z zero-divisor graph: an edge-difference Dirichlet
receiver, a vertex-mass receiver, a graph Laplacian/Hessian receiver, and a
formal Euler row for every vertex.  These receivers produce exact finite rank,
nullity, Euler, and Hessian data.  They still are not strict action terms: no
current artifact exports a unit-bearing graph measure, strict field provenance,
boundary/integration theorem, named graph density theorem, or nonproxy continuum
lift.
"""
from __future__ import annotations

import hashlib, json
from fractions import Fraction
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2996_s1946_z12_zero_divisor_graph_source_candidate_obstruction import zero_divisor_graph_witness
from p2999_s1949_zero_divisor_graph_named_source_atom_coupling_obstruction import OUT as P2999

OUT = GEN / "p3000_s1950_zero_divisor_graph_action_installation_obstruction.json"
MD = GEN / "p3000_s1950_zero_divisor_graph_action_installation_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def rank_over_q(matrix: list[list[int]]) -> int:
    rows = [[Fraction(x) for x in row] for row in matrix]
    if not rows:
        return 0
    m, n = len(rows), len(rows[0])
    rank = 0
    col = 0
    while rank < m and col < n:
        pivot = next((r for r in range(rank, m) if rows[r][col] != 0), None)
        if pivot is None:
            col += 1
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pv = rows[rank][col]
        rows[rank] = [x / pv for x in rows[rank]]
        for r in range(m):
            if r != rank and rows[r][col] != 0:
                factor = rows[r][col]
                rows[r] = [a - factor * b for a, b in zip(rows[r], rows[rank])]
        rank += 1
        col += 1
    return rank


def graph_action_receivers() -> dict[str, Any]:
    graph = zero_divisor_graph_witness()
    vertices = graph["vertices"]
    edges = [tuple(edge) for edge in graph["edges"]]
    index = {v: i for i, v in enumerate(vertices)}
    n = len(vertices)
    incidence = []
    for a, b in edges:
        row = [0] * n
        row[index[a]] = 1
        row[index[b]] = -1
        incidence.append(row)
    laplacian = [[0] * n for _ in range(n)]
    for a, b in edges:
        ia, ib = index[a], index[b]
        laplacian[ia][ia] += 1
        laplacian[ib][ib] += 1
        laplacian[ia][ib] -= 1
        laplacian[ib][ia] -= 1
    degrees = [laplacian[i][i] for i in range(n)]
    euler_rows = []
    for v in vertices:
        i = index[v]
        euler_rows.append({
            "vertex": v,
            "formal_euler_row": laplacian[i],
            "degree": laplacian[i][i],
            "neighbor_coefficients": {str(vertices[j]): laplacian[i][j] for j in range(n) if laplacian[i][j] != 0 and j != i},
            "row_sum_zero": sum(laplacian[i]) == 0,
            "unit_bearing_graph_measure": False,
            "strict_field_provenance": False,
            "boundary_integration_theorem": False,
            "named_graph_density_theorem": False,
            "nonproxy_continuum_lift": False,
            "accepted_action_installation": False,
        })
    incidence_rank = rank_over_q(incidence)
    laplacian_rank = rank_over_q(laplacian)
    receiver_rows = [
        {
            "receiver": "formal_zero_divisor_graph_dirichlet_edge_density",
            "density_template": "L_G = (mu_G/2) * sum_{(a,b) in E0} (phi_a - phi_b)^2",
            "edge_count": len(edges),
            "incidence_matrix": incidence,
            "incidence_rank_over_Q": incidence_rank,
            "nonzero_formal_variation": len(edges) > 0,
            "unit_bearing_graph_measure": False,
            "strict_field_provenance": False,
            "boundary_integration_theorem": False,
            "named_graph_density_theorem": False,
            "nonproxy_continuum_lift": False,
            "accepted_action_installation": False,
        },
        {
            "receiver": "formal_zero_divisor_graph_vertex_mass_density",
            "density_template": "M_G = (m_G/2) * sum_{v in V0} deg(v) * phi_v^2",
            "vertex_count": n,
            "degree_weights": degrees,
            "degree_weight_sum": sum(degrees),
            "nonzero_formal_variation": any(d != 0 for d in degrees),
            "unit_bearing_graph_measure": False,
            "strict_field_provenance": False,
            "boundary_integration_theorem": False,
            "named_graph_density_theorem": False,
            "nonproxy_continuum_lift": False,
            "accepted_action_installation": False,
        },
        {
            "receiver": "formal_zero_divisor_graph_laplacian_hessian",
            "density_template": "H_G = mu_G * B^T B, the graph Laplacian of the zero-product graph",
            "laplacian_matrix": laplacian,
            "laplacian_rank_over_Q": laplacian_rank,
            "laplacian_nullity_over_Q": n - laplacian_rank,
            "row_sums_zero": all(sum(row) == 0 for row in laplacian),
            "nonzero_formal_variation": laplacian_rank > 0,
            "unit_bearing_graph_measure": False,
            "strict_field_provenance": False,
            "boundary_integration_theorem": False,
            "named_graph_density_theorem": False,
            "nonproxy_continuum_lift": False,
            "accepted_action_installation": False,
        },
        {
            "receiver": "completed_strict_zero_divisor_graph_action_installation_schema",
            "density_template": "named unit-bearing graph density plus boundary/integration map and nonproxy continuum lift",
            "required_positive_inputs": ["unit_bearing_graph_measure", "strict_field_provenance", "boundary_integration_theorem", "named_graph_density_theorem", "nonproxy_continuum_lift"],
            "nonzero_formal_variation": True,
            "unit_bearing_graph_measure": True,
            "strict_field_provenance": True,
            "boundary_integration_theorem": True,
            "named_graph_density_theorem": True,
            "nonproxy_continuum_lift": True,
            "accepted_action_installation": False,
        },
    ]
    receiver_rows.extend({"receiver": "formal_zero_divisor_graph_vertex_euler_row", **row} for row in euler_rows)
    return {
        "vertex_count": n,
        "edge_count": len(edges),
        "degree_sequence": sorted(degrees),
        "degree_weight_sum": sum(degrees),
        "incidence_rank_over_Q": incidence_rank,
        "laplacian_rank_over_Q": laplacian_rank,
        "laplacian_nullity_over_Q": n - laplacian_rank,
        "laplacian_row_sums_zero": all(sum(row) == 0 for row in laplacian),
        "euler_row_count": len(euler_rows),
        "receiver_rows": receiver_rows,
    }


def obligation_rows(witness: dict[str, Any]) -> list[dict[str, Any]]:
    current = [r for r in witness["receiver_rows"] if r["receiver"] != "completed_strict_zero_divisor_graph_action_installation_schema"]
    return [
        {"obligation": "finite_graph_action_receivers", "satisfied": witness["vertex_count"] == 7 and witness["edge_count"] == 8, "evidence": "formal Dirichlet, vertex-mass, Laplacian/Hessian, and vertex Euler receivers are built from the zero-product graph"},
        {"obligation": "nonzero_formal_variation", "satisfied": any(r["nonzero_formal_variation"] for r in current), "evidence": f"incidence rank {witness['incidence_rank_over_Q']} and Laplacian rank {witness['laplacian_rank_over_Q']} give nonzero toy Euler/Hessian data"},
        {"obligation": "graph_laplacian_boundary_witness", "satisfied": witness["laplacian_row_sums_zero"] and witness["laplacian_nullity_over_Q"] == 1, "evidence": "the connected finite graph has a rank-6 Laplacian with one constant-mode null direction"},
        {"obligation": "unit_bearing_graph_measure", "satisfied": any(r["unit_bearing_graph_measure"] for r in current), "evidence": "mu_G and m_G are symbolic; counting edges/vertices is not a strict unit/measure theorem"},
        {"obligation": "strict_field_provenance", "satisfied": any(r["strict_field_provenance"] for r in current), "evidence": "phi_v fields are formal graph placeholders, not strict fields sourced by the nadsoliton"},
        {"obligation": "boundary_integration_theorem", "satisfied": any(r["boundary_integration_theorem"] for r in current), "evidence": "finite graph sums are not a boundary/integration map"},
        {"obligation": "named_graph_density_theorem", "satisfied": any(r["named_graph_density_theorem"] for r in current), "evidence": "formal receiver names are not exported named graph action densities"},
        {"obligation": "nonproxy_continuum_lift", "satisfied": any(r["nonproxy_continuum_lift"] for r in current), "evidence": "finite graph Laplacian data are not lifted to tensor-resolved nonproxy continuum action"},
        {"obligation": "accepted_action_installation", "satisfied": any(r["accepted_action_installation"] for r in current), "evidence": "formal graph receivers fail strict installation premises"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_receivers", "nonzero_variation", "laplacian_boundary_witness", "unit_measure", "field_provenance", "boundary_integration", "named_density", "nonproxy_lift"]
    return [{"present": dict(zip(names, bits)), "accepts_graph_action_installation": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p2999_path: Any) -> dict[str, Any]:
    witness = graph_action_receivers()
    matrix = acceptance_matrix()
    return {
        "status": "P3000_ZERO_DIVISOR_GRAPH_ACTION_VARIATIONAL_INSTALLATION_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P2999": hashlib.sha256(p2999_path.read_bytes()).hexdigest() if p2999_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "ZeroDivisorGraphActionVariationalInstallation_ObstructionMatrix",
            "formal_action_witness": witness,
            "receiver_rows": witness["receiver_rows"],
            "proof_obligation_rows": obligation_rows(witness),
            "finite_acceptance_matrix": matrix,
        },
        "installation_certificate": {
            "receiver_count": len(witness["receiver_rows"]),
            "vertex_count": witness["vertex_count"],
            "edge_count": witness["edge_count"],
            "degree_sequence": witness["degree_sequence"],
            "degree_weight_sum": witness["degree_weight_sum"],
            "incidence_rank_over_Q": witness["incidence_rank_over_Q"],
            "laplacian_rank_over_Q": witness["laplacian_rank_over_Q"],
            "laplacian_nullity_over_Q": witness["laplacian_nullity_over_Q"],
            "laplacian_row_sums_zero": witness["laplacian_row_sums_zero"],
            "euler_row_count": witness["euler_row_count"],
            "accepted_current_installations": [r["receiver"] for r in witness["receiver_rows"] if r["accepted_action_installation"]],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_graph_action_installation"]),
        },
        "decision": {
            "positive_progress": "P3000 attacks the final zero-divisor-graph route after P2999: action/variational installation with unit-bearing graph measure and nonproxy lift obligations.",
            "breakthrough": "Bounded no-go: formal Dirichlet, vertex-mass, Laplacian/Hessian, and vertex Euler receivers give exact finite toy variational data, but no current row exports a unit-bearing graph measure, strict field provenance, boundary/integration theorem, named graph density theorem, or nonproxy continuum lift.",
            "negative_export_flags": {k: False for k in ["zero_divisor_graph_action_installation_exported", "unit_bearing_graph_density_exported", "nonproxy_variational_chain_exported", "strict_field_provenance_exported", "named_graph_density_theorem_exported", "selector_closure_exported", "bridge_closure_exported", "nonproxy_ltotal_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay zero-divisor graph formal action receivers, Laplacian ranks, source-coupling rows, provenance/localizer rows, Fourier/annihilator/nilradical/CRT/zero-derivation lanes, selector replay, bridge maps, role transfer, or L_total placeholders.  The zero-divisor-graph lane is now bounded no-go on current artifacts; the next proof-grade move must introduce a genuinely new strict typed object/theorem/provider outside this lane, or preserve the P2929-P3000 no-strict-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["installation_certificate"]
    lines = [
        "# P3000/S1950 zero-divisor graph action/variational installation obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Installation certificate",
        f"- receiver rows: `{cert['receiver_count']}`",
        f"- vertices/edges: `{cert['vertex_count']}/{cert['edge_count']}`",
        f"- degree sequence: `{cert['degree_sequence']}`",
        f"- degree weight sum: `{cert['degree_weight_sum']}`",
        f"- incidence rank over Q: `{cert['incidence_rank_over_Q']}`",
        f"- Laplacian rank/nullity over Q: `{cert['laplacian_rank_over_Q']}/{cert['laplacian_nullity_over_Q']}`",
        f"- Laplacian row sums zero: `{cert['laplacian_row_sums_zero']}`",
        f"- Euler rows: `{cert['euler_row_count']}`",
        f"- accepted current installations: `{cert['accepted_current_installations']}`",
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
    read_json(P2999)
    payload = build_payload(P2999)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3000/S1950 zero-divisor graph action/variational installation obstruction", "## P3000/S1950 zero-divisor graph action/variational installation obstruction\n\n`P3000/S1950` attacks the final zero-divisor-graph route left by P2999: action/variational installation with a genuinely unit-bearing graph measure, named graph density theorem, boundary/integration map, and nonproxy continuum lift.  The finite matrix builds a formal Dirichlet edge-density receiver, a vertex-mass receiver, a graph Laplacian/Hessian receiver, and one formal Euler row per vertex.  The toy incidence and Laplacian ranks are exact (`6`), the Laplacian has one constant null direction, and all row sums vanish; however, no row exports a unit-bearing graph measure, strict field provenance, boundary/integration theorem, named graph density theorem, or nonproxy continuum lift.  No zero-divisor graph action installation, nonproxy variational chain, bridge closure, nonproxy `L_total`, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3000/S1950 zero-divisor graph action/variational installation `L_total` guard", "## P3000/S1950 zero-divisor graph action/variational installation `L_total` guard\n\n`P3000/S1950` explicitly tests formal zero-divisor graph action receivers but adds no sourced graph term to `L_total`.  The nonzero toy Dirichlet, vertex-mass, Laplacian/Hessian, and Euler-row coefficients are formal finite graph algebra only; no unit-bearing graph measure, strict field provenance, boundary/integration theorem, named graph density theorem, nonproxy lift, EOM, Hamiltonian, role transfer, bridge closure, or ToE is exported.\n")
    append_once(AGENTS, "Current zero-divisor graph action/variational installation obstruction guardrail (P3000/S1950, 2026-06-20)", "## Current zero-divisor graph action/variational installation obstruction guardrail (P3000/S1950, 2026-06-20)\n\n- P3000 attacks the final zero-divisor-graph route left by P2999: action/variational installation with unit-bearing graph measure and nonproxy lift obligations.\n- Formal receiver positives are limited: a Dirichlet edge-density receiver, vertex-mass receiver, graph Laplacian/Hessian receiver, and one Euler row per vertex give exact finite toy variational data; the acceptance matrix has `256` profiles with only the full profile accepting.\n- The current route is bounded no-go: no unit-bearing graph measure, strict field provenance, boundary/integration theorem, named graph density theorem, or nonproxy continuum lift is exported.\n- Do not promote graph action receivers, Laplacian ranks, source-coupling rows, provenance/localizer rows, symbolic `L_total` slots, selector replay, bridge maps, role transfer, nonproxy `L_total`, or ToE.  The zero-divisor-graph lane is now bounded no-go on current artifacts; a next move must introduce a genuinely new strict typed object/theorem/provider outside this lane or preserve the P2929-P3000 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
