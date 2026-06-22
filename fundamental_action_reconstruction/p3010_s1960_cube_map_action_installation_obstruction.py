#!/usr/bin/env python3
"""P3010/S1960: cube-map action/variational installation obstruction.

P3009 left exactly one cube-map-functional-graph route: unit-bearing
action/variational installation. This audit attacks that final route without
replaying source-atom coupling, localizer, strict provenance, square-map,
zero-divisor graph, CRT-idempotent, nilradical, annihilator, Fourier, selector,
bridge, role-transfer, or L_total lanes.

The constructed object is a finite formal action receiver matrix for the cubing
map x -> x^3 mod 12: directed-edge Dirichlet receivers, fixed-attractor pinning,
a Hessian/Laplacian receiver, and one formal Euler row per residue. The exact
rank computation gives a three-edge nonloop incidence rank and a rank-3 Hessian,
but these are finite toy variational data only; current artifacts export no
unit-bearing cube-map measure, named cube-map density theorem,
boundary/integration map, nonproxy continuum lift, or EOM/Hamiltonian source.
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
from p3006_s1956_z12_cube_map_functional_graph_source_candidate_obstruction import cube_map
from p3009_s1959_cube_map_named_source_atom_coupling_obstruction import OUT as P3009

OUT = GEN / "p3010_s1960_cube_map_action_installation_obstruction.json"
MD = GEN / "p3010_s1960_cube_map_action_installation_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
MODULUS = 12


def rational_rank(matrix: list[list[int]]) -> int:
    a = [[Fraction(v) for v in row] for row in matrix]
    if not a:
        return 0
    m, n, r = len(a), len(a[0]), 0
    for c in range(n):
        pivot = next((i for i in range(r, m) if a[i][c] != 0), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        pv = a[r][c]
        a[r] = [v / pv for v in a[r]]
        for i in range(m):
            if i != r and a[i][c] != 0:
                factor = a[i][c]
                a[i] = [a[i][j] - factor * a[r][j] for j in range(n)]
        r += 1
    return r


def action_witness() -> dict[str, Any]:
    residues = list(range(MODULUS))
    directed_edges = [(x, cube_map(x)) for x in residues]
    nonloop_edges = [(s, t) for s, t in directed_edges if s != t]
    loop_edges = [(s, t) for s, t in directed_edges if s == t]
    incidence = []
    for s, t in nonloop_edges:
        row = [0] * MODULUS
        row[s] = -1
        row[t] = 1
        incidence.append(row)
    hessian = [[sum(row[i] * row[j] for row in incidence) for j in residues] for i in residues]
    euler_rows = []
    for i in residues:
        incoming = [s for s, t in nonloop_edges if t == i]
        outgoing = [t for s, t in nonloop_edges if s == i]
        euler_rows.append({
            "residue": i,
            "hessian_row": hessian[i],
            "incoming_nonloop_sources": incoming,
            "outgoing_nonloop_targets": outgoing,
            "formal_euler_operator_nonzero": any(v != 0 for v in hessian[i]),
            "unit_bearing_density_installed": False,
        })
    rank_b = rational_rank(incidence)
    rank_h = rational_rank(hessian)
    return {
        "node_count": len(residues),
        "directed_edge_count": len(directed_edges),
        "nonloop_edge_count": len(nonloop_edges),
        "loop_edge_count": len(loop_edges),
        "nonloop_edges": nonloop_edges,
        "loop_edges": loop_edges,
        "directed_incidence_matrix": incidence,
        "directed_incidence_rank": rank_b,
        "hessian_laplacian_matrix": hessian,
        "hessian_rank": rank_h,
        "hessian_nullity": MODULUS - rank_h,
        "hessian_row_sums": [sum(row) for row in hessian],
        "all_hessian_row_sums_zero": all(sum(row) == 0 for row in hessian),
        "formal_boundary_components": sorted(s for s, t in loop_edges),
        "euler_rows": euler_rows,
        "formal_action_receivers": [
            "directed_edge_dirichlet_receiver",
            "fixed_attractor_pinning_receiver",
            "hessian_laplacian_receiver",
            "per_residue_formal_euler_receiver",
        ],
        "accepted_action_installations": [],
    }


def proof_obligations(w: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {"obligation": "finite_cube_map_action_receivers", "satisfied": w["node_count"] == 12 and w["directed_edge_count"] == 12 and w["nonloop_edge_count"] == 3, "evidence": "formal receivers are built from all cubing edges and the three nonloop edges"},
        {"obligation": "nonzero_formal_variation", "satisfied": w["directed_incidence_rank"] == 3 and w["hessian_rank"] == 3 and w["hessian_nullity"] == 9, "evidence": "directed incidence rank is 3 and Hessian rank/nullity is 3/9"},
        {"obligation": "hessian_boundary_witness", "satisfied": w["all_hessian_row_sums_zero"] and w["formal_boundary_components"] == [0, 1, 3, 4, 5, 7, 8, 9, 11], "evidence": "row sums vanish and nine fixed attractors remain formal boundary components"},
        {"obligation": "unit_bearing_cube_map_measure", "satisfied": False, "evidence": "no measure with physical units is exported for cube-map receivers"},
        {"obligation": "named_cube_map_density_theorem", "satisfied": False, "evidence": "no theorem names a cube-map action density sourced by strict artifacts"},
        {"obligation": "boundary_integration_map", "satisfied": False, "evidence": "finite graph boundary bookkeeping is not an integration map"},
        {"obligation": "nonproxy_continuum_lift", "satisfied": False, "evidence": "no continuum/nonproxy variational lift is exported"},
        {"obligation": "accepted_current_action_installation", "satisfied": bool(w["accepted_action_installations"]), "evidence": "no row satisfies measure, density, boundary/integration, and nonproxy lift premises"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_receivers", "formal_variation", "hessian_boundary", "unit_measure", "named_density", "boundary_integration", "nonproxy_lift", "eom_hamiltonian_source"]
    return [{"present": dict(zip(names, bits)), "accepts_action_installation": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p3009_path: Any) -> dict[str, Any]:
    witness = action_witness()
    matrix = acceptance_matrix()
    return {
        "status": "P3010_CUBE_MAP_ACTION_INSTALLATION_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P3009": hashlib.sha256(p3009_path.read_bytes()).hexdigest() if p3009_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "CubeMapActionVariationalInstallation_ObstructionMatrix",
            "action_witness": witness,
            "proof_obligation_rows": proof_obligations(witness),
            "finite_acceptance_matrix": matrix,
        },
        "action_certificate": {
            "node_count": witness["node_count"],
            "directed_edge_count": witness["directed_edge_count"],
            "nonloop_edge_count": witness["nonloop_edge_count"],
            "loop_edge_count": witness["loop_edge_count"],
            "nonloop_edges": witness["nonloop_edges"],
            "directed_incidence_rank": witness["directed_incidence_rank"],
            "hessian_rank": witness["hessian_rank"],
            "hessian_nullity": witness["hessian_nullity"],
            "all_hessian_row_sums_zero": witness["all_hessian_row_sums_zero"],
            "formal_boundary_components": witness["formal_boundary_components"],
            "euler_row_count": len(witness["euler_rows"]),
            "accepted_action_installations": witness["accepted_action_installations"],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for row in matrix if row["accepts_action_installation"]),
        },
        "decision": {
            "positive_progress": "P3010 attacks the final P3009 cube-map-functional-graph route: unit-bearing action/variational installation.",
            "breakthrough": "Bounded no-go: formal directed-edge Dirichlet, fixed-attractor pinning, Hessian/Laplacian, and Euler receivers give exact finite toy variational data, but no unit-bearing cube-map measure, named density theorem, boundary/integration map, or nonproxy continuum lift is exported.",
            "negative_export_flags": {k: False for k in ["cube_map_action_installation_exported", "unit_bearing_cube_map_measure_exported", "named_cube_map_density_theorem_exported", "boundary_integration_map_exported", "nonproxy_continuum_lift_exported", "eom_hamiltonian_source_exported", "bridge_closure_exported", "role_transfer_exported", "nonproxy_ltotal_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay cube-map-functional-graph lane: after P3010 it is bounded no-go on current artifacts. A next proof-grade move must introduce a genuinely new strict typed object/theorem/provider outside this lane, or preserve the P2929-P3010 no-strict-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["action_certificate"]
    lines = [
        "# P3010/S1960 cube-map action/variational installation obstruction", "",
        f"Status: `{payload['status']}`", "", "## Action certificate",
        f"- nodes/edges: `{cert['node_count']}/{cert['directed_edge_count']}`",
        f"- nonloop/loop edges: `{cert['nonloop_edge_count']}/{cert['loop_edge_count']}`",
        f"- nonloop edges: `{cert['nonloop_edges']}`",
        f"- incidence rank: `{cert['directed_incidence_rank']}`",
        f"- Hessian rank/nullity: `{cert['hessian_rank']}/{cert['hessian_nullity']}`",
        f"- formal boundary components: `{cert['formal_boundary_components']}`",
        f"- accepted action installations: `{cert['accepted_action_installations']}`",
        f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_rows']}`", "",
        "## Lay summary", payload["decision"]["positive_progress"], payload["decision"]["breakthrough"], "",
        "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P3009)
    payload = build_payload(P3009)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3010/S1960 cube-map action/variational installation obstruction", "## P3010/S1960 cube-map action/variational installation obstruction\n\n`P3010/S1960` attacks the final cube-map-functional-graph route left by P3009: action/variational installation with unit-bearing cube-map measure and nonproxy lift obligations.  Formal receiver positives are limited: directed-edge Dirichlet receivers on the three nonloop edges `(2,8)`, `(6,0)`, `(10,4)`, fixed-attractor pinning receivers, a Hessian/Laplacian receiver, and one Euler row per residue give exact finite toy variational data; the directed incidence rank is `3`, the Hessian rank/nullity is `3/9`, and the nine fixed attractors remain formal boundary components.  The current route is bounded no-go: no unit-bearing cube-map measure, strict field provenance, boundary/integration theorem, named cube-map density theorem, or nonproxy continuum lift is exported.  No sourced `L_total`, EOM/Hamiltonian, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3010/S1960 cube-map action-installation `L_total` guard", "## P3010/S1960 cube-map action-installation `L_total` guard\n\n`P3010/S1960` adds no cube-map action term to `L_total`.  Directed-edge Dirichlet, fixed-attractor pinning, Hessian/Laplacian, and Euler receivers are finite formal variational data only; they do not supply a unit-bearing cube-map measure, strict field provenance, named density theorem, boundary/integration theorem, nonproxy variational chain, EOM/Hamiltonian source, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current cube-map action/variational installation obstruction guardrail (P3010/S1960, 2026-06-22)", "## Current cube-map action/variational installation obstruction guardrail (P3010/S1960, 2026-06-22)\n\n- P3010 attacks the final cube-map-functional-graph route left by P3009: unit-bearing action/variational installation with cube-map measure and nonproxy lift obligations.\n- Formal receiver positives are finite only: directed-edge Dirichlet receivers on nonloop edges `(2,8)`, `(6,0)`, `(10,4)`, fixed-attractor pinning, Hessian/Laplacian, and one Euler row per residue give exact toy variational data; incidence rank is `3`, Hessian rank/nullity is `3/9`, and the nine fixed attractors are formal boundary components.\n- The current route is bounded no-go: no unit-bearing cube-map measure, strict field provenance, boundary/integration theorem, named cube-map density theorem, or nonproxy continuum lift is exported.\n- Do not promote cube-map action receivers, Hessian ranks, source-coupling rows, provenance/localizer rows, symbolic `L_total` slots, selector replay, bridge maps, role transfer, nonproxy `L_total`, or ToE.  The cube-map-functional-graph lane is now bounded no-go on current artifacts; a next move must introduce a genuinely new strict typed object/theorem/provider outside this lane or preserve the P2929-P3010 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
