#!/usr/bin/env python3
"""P3005/S1955: square-map action/variational installation obstruction.

P3004 left the final square-map-functional-graph route: unit-bearing action
installation with a genuinely unit-bearing square-map measure, named square-map
density theorem, boundary/integration map, and nonproxy continuum lift. This
audit attacks that route only. It does not replay strict provenance,
basin/localizer, source-atom coupling, zero-divisor graph, CRT idempotent,
nilradical, annihilator, Fourier, selector, bridge completion, role transfer, or
L_total promotion.

The finite calculation constructs the strongest formal square-map action
receivers available from the Z/12Z squaring functional graph: a directed-edge
Dirichlet receiver, an attractor-pinning receiver, a graph Hessian/Laplacian
receiver, and one formal Euler row per residue. These receivers provide exact
finite rank, nullity, Euler, and boundary data, but do not install a strict action
term because no unit-bearing square-map measure, density theorem,
boundary/integration theorem, strict field provenance, or nonproxy continuum lift
is exported.
"""
from __future__ import annotations

import hashlib, json
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3000_s1950_zero_divisor_graph_action_installation_obstruction import rank_over_q
from p3001_s1951_z12_square_map_functional_graph_source_candidate_obstruction import MODULUS, square_map_witness
from p3004_s1954_square_map_named_source_atom_coupling_obstruction import OUT as P3004

OUT = GEN / "p3005_s1955_square_map_action_installation_obstruction.json"
MD = GEN / "p3005_s1955_square_map_action_installation_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def action_receivers() -> dict[str, Any]:
    graph = square_map_witness()
    vertices = list(range(MODULUS))
    edges = [(edge["source"], edge["target"]) for edge in graph["edges"]]
    index = {v: i for i, v in enumerate(vertices)}
    n = len(vertices)
    incidence = []
    nonloop_edges = []
    loop_edges = []
    for source, target in edges:
        row = [0] * n
        row[index[source]] = 1
        row[index[target]] -= 1
        incidence.append(row)
        (loop_edges if source == target else nonloop_edges).append((source, target))
    hessian = [[0] * n for _ in range(n)]
    for row in incidence:
        for i, a in enumerate(row):
            if a == 0:
                continue
            for j, b in enumerate(row):
                if b != 0:
                    hessian[i][j] += a * b
    attractors = graph["fixed_points"]
    basin_sizes = {int(k): v for k, v in graph["basin_sizes"].items()}
    euler_rows = []
    for v in vertices:
        i = index[v]
        euler_rows.append({
            "residue": v,
            "formal_euler_row": hessian[i],
            "diagonal_weight": hessian[i][i],
            "off_diagonal_coefficients": {str(vertices[j]): hessian[i][j] for j in range(n) if hessian[i][j] != 0 and j != i},
            "row_sum_zero": sum(hessian[i]) == 0,
            "unit_bearing_square_map_measure": False,
            "strict_field_provenance": False,
            "boundary_integration_theorem": False,
            "named_square_map_density_theorem": False,
            "nonproxy_continuum_lift": False,
            "accepted_action_installation": False,
        })
    incidence_rank = rank_over_q(incidence)
    hessian_rank = rank_over_q(hessian)
    receiver_rows = [
        {
            "receiver": "formal_square_map_directed_edge_dirichlet_density",
            "density_template": "L_sq = (mu_sq/2) * sum_x (phi_x - phi_{x^2})^2",
            "edge_count": len(edges),
            "nonloop_edge_count": len(nonloop_edges),
            "loop_edge_count": len(loop_edges),
            "incidence_rank_over_Q": incidence_rank,
            "nonzero_formal_variation": bool(nonloop_edges),
            "unit_bearing_square_map_measure": False,
            "strict_field_provenance": False,
            "boundary_integration_theorem": False,
            "named_square_map_density_theorem": False,
            "nonproxy_continuum_lift": False,
            "accepted_action_installation": False,
        },
        {
            "receiver": "formal_square_map_attractor_pinning_density",
            "density_template": "P_sq = (kappa_sq/2) * sum_{a fixed} |B_a| * phi_a^2",
            "attractor_count": len(attractors),
            "attractor_basin_weights": {str(a): basin_sizes[a] for a in attractors},
            "weight_sum": sum(basin_sizes[a] for a in attractors),
            "nonzero_formal_variation": bool(attractors),
            "unit_bearing_square_map_measure": False,
            "strict_field_provenance": False,
            "boundary_integration_theorem": False,
            "named_square_map_density_theorem": False,
            "nonproxy_continuum_lift": False,
            "accepted_action_installation": False,
        },
        {
            "receiver": "formal_square_map_hessian_laplacian_density",
            "density_template": "H_sq = B_sq^T B_sq for directed incidence B_sq",
            "hessian_matrix": hessian,
            "hessian_rank_over_Q": hessian_rank,
            "hessian_nullity_over_Q": n - hessian_rank,
            "all_row_sums_zero": all(sum(row) == 0 for row in hessian),
            "nonzero_formal_variation": hessian_rank > 0,
            "unit_bearing_square_map_measure": False,
            "strict_field_provenance": False,
            "boundary_integration_theorem": False,
            "named_square_map_density_theorem": False,
            "nonproxy_continuum_lift": False,
            "accepted_action_installation": False,
        },
    ]
    return {
        "vertices": vertices,
        "edges": edges,
        "nonloop_edges": nonloop_edges,
        "loop_edges": loop_edges,
        "incidence_matrix": incidence,
        "hessian_matrix": hessian,
        "incidence_rank_over_Q": incidence_rank,
        "hessian_rank_over_Q": hessian_rank,
        "hessian_nullity_over_Q": n - hessian_rank,
        "all_hessian_row_sums_zero": all(sum(row) == 0 for row in hessian),
        "euler_rows": euler_rows,
        "receiver_rows": receiver_rows,
        "formal_boundary_components": graph["fixed_points"],
        "unit_bearing_square_map_measure_exported": False,
        "strict_field_provenance_exported": False,
        "boundary_integration_theorem_exported": False,
        "named_square_map_density_theorem_exported": False,
        "nonproxy_continuum_lift_exported": False,
        "accepted_action_installations": [],
    }


def obligation_rows(w: dict[str, Any]) -> list[dict[str, Any]]:
    rows = w["receiver_rows"] + w["euler_rows"]
    return [
        {"obligation": "finite_square_map_action_receivers", "satisfied": len(w["receiver_rows"]) == 3 and len(w["euler_rows"]) == 12, "evidence": "Dirichlet, attractor-pinning, Hessian, and one Euler row per residue are computed"},
        {"obligation": "nonzero_formal_variation", "satisfied": any(r["nonzero_formal_variation"] for r in w["receiver_rows"]), "evidence": "nonloop squaring edges yield nonzero formal variation"},
        {"obligation": "hessian_boundary_witness", "satisfied": w["incidence_rank_over_Q"] == 8 and w["hessian_rank_over_Q"] == 8 and w["hessian_nullity_over_Q"] == 4 and w["all_hessian_row_sums_zero"], "evidence": "directed incidence/Hessian ranks expose four attractor boundary components"},
        {"obligation": "unit_bearing_square_map_measure", "satisfied": any(r["unit_bearing_square_map_measure"] for r in rows), "evidence": "formal mu_sq/kappa_sq symbols are not exported unit-bearing measures"},
        {"obligation": "strict_field_provenance", "satisfied": any(r["strict_field_provenance"] for r in rows), "evidence": "P3002 exported no strict square-map field provenance"},
        {"obligation": "boundary_integration_theorem", "satisfied": any(r["boundary_integration_theorem"] for r in rows), "evidence": "finite boundary components are not a continuum boundary/integration map"},
        {"obligation": "named_square_map_density_theorem", "satisfied": any(r["named_square_map_density_theorem"] for r in rows), "evidence": "no named square-map action density theorem is exported"},
        {"obligation": "nonproxy_continuum_lift", "satisfied": any(r["nonproxy_continuum_lift"] for r in rows), "evidence": "finite formal receivers are not nonproxy continuum lifts"},
        {"obligation": "accepted_current_action_installation", "satisfied": bool(w["accepted_action_installations"]), "evidence": "no row satisfies measure, provenance, boundary/integration, density theorem, and lift premises"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_action_receiver", "nonzero_variation", "hessian_boundary", "unit_measure", "strict_field_provenance", "boundary_integration", "named_density_theorem", "nonproxy_lift"]
    return [{"present": dict(zip(names, bits)), "accepts_action_installation": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p3004_path: Any) -> dict[str, Any]:
    witness = action_receivers()
    matrix = acceptance_matrix()
    return {
        "status": "P3005_SQUARE_MAP_ACTION_INSTALLATION_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P3004": hashlib.sha256(p3004_path.read_bytes()).hexdigest() if p3004_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "SquareMapActionVariationalInstallation_ObstructionMatrix",
            "action_witness": witness,
            "proof_obligation_rows": obligation_rows(witness),
            "finite_acceptance_matrix": matrix,
        },
        "action_certificate": {
            "vertex_count": len(witness["vertices"]),
            "edge_count": len(witness["edges"]),
            "nonloop_edge_count": len(witness["nonloop_edges"]),
            "loop_edge_count": len(witness["loop_edges"]),
            "incidence_rank_over_Q": witness["incidence_rank_over_Q"],
            "hessian_rank_over_Q": witness["hessian_rank_over_Q"],
            "hessian_nullity_over_Q": witness["hessian_nullity_over_Q"],
            "all_hessian_row_sums_zero": witness["all_hessian_row_sums_zero"],
            "euler_row_count": len(witness["euler_rows"]),
            "formal_boundary_components": witness["formal_boundary_components"],
            "accepted_action_installation_count": len(witness["accepted_action_installations"]),
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for row in matrix if row["accepts_action_installation"]),
        },
        "decision": {
            "positive_progress": "P3005 attacks the final square-map-functional-graph route left by P3004: unit-bearing action/variational installation.",
            "breakthrough": "Bounded no-go: formal directed-edge Dirichlet, attractor-pinning, Hessian/Laplacian, and Euler receivers give exact finite rank/nullity data, but no unit-bearing square-map measure, strict field provenance, boundary/integration theorem, named square-map density theorem, or nonproxy continuum lift is exported.",
            "negative_export_flags": {k: False for k in ["square_map_action_installation_exported", "unit_bearing_square_map_measure_exported", "named_square_map_density_theorem_exported", "boundary_integration_theorem_exported", "nonproxy_continuum_lift_exported", "nonproxy_ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay square-map-functional-graph lane: after P3005 it is bounded no-go on current artifacts. A next proof-grade move must introduce a genuinely new strict typed object/theorem/provider outside this lane, or preserve the P2929-P3005 no-strict-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["action_certificate"]
    lines = [
        "# P3005/S1955 square-map action/variational installation obstruction", "",
        f"Status: `{payload['status']}`", "", "## Action certificate",
        f"- vertices/edges: `{cert['vertex_count']}/{cert['edge_count']}`",
        f"- nonloop/loop edges: `{cert['nonloop_edge_count']}/{cert['loop_edge_count']}`",
        f"- incidence rank over Q: `{cert['incidence_rank_over_Q']}`",
        f"- Hessian rank/nullity over Q: `{cert['hessian_rank_over_Q']}/{cert['hessian_nullity_over_Q']}`",
        f"- all Hessian row sums zero: `{cert['all_hessian_row_sums_zero']}`",
        f"- Euler rows: `{cert['euler_row_count']}`",
        f"- formal boundary components: `{cert['formal_boundary_components']}`",
        f"- accepted action installations: `{cert['accepted_action_installation_count']}`",
        f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_rows']}`", "",
        "## Lay summary", payload["decision"]["positive_progress"], payload["decision"]["breakthrough"], "",
        "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P3004)
    payload = build_payload(P3004)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3005/S1955 square-map action/variational installation obstruction", "## P3005/S1955 square-map action/variational installation obstruction\n\n`P3005/S1955` attacks the final square-map-functional-graph route left by P3004: action/variational installation with unit-bearing square-map measure and nonproxy lift obligations.  Formal receiver positives are limited: a directed-edge Dirichlet receiver, attractor-pinning receiver, Hessian/Laplacian receiver, and one Euler row per residue give exact finite toy variational data; the directed incidence rank is `8`, the Hessian rank/nullity is `8/4`, and the four fixed attractors remain formal boundary components.  The current route is bounded no-go: no unit-bearing square-map measure, strict field provenance, boundary/integration theorem, named square-map density theorem, or nonproxy continuum lift is exported.  No sourced `L_total`, EOM/Hamiltonian, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3005/S1955 square-map action-installation `L_total` guard", "## P3005/S1955 square-map action-installation `L_total` guard\n\n`P3005/S1955` adds no square-map action term to `L_total`.  Directed-edge Dirichlet, attractor-pinning, Hessian/Laplacian, and Euler receivers are finite formal variational data only; they do not supply a unit-bearing square-map measure, strict field provenance, named density theorem, boundary/integration theorem, nonproxy variational chain, EOM/Hamiltonian source, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current square-map action/variational installation obstruction guardrail (P3005/S1955, 2026-06-22)", "## Current square-map action/variational installation obstruction guardrail (P3005/S1955, 2026-06-22)\n\n- P3005 attacks the final square-map-functional-graph route left by P3004: unit-bearing action/variational installation with square-map measure and nonproxy lift obligations.\n- Formal receiver positives are finite only: directed-edge Dirichlet, attractor-pinning, Hessian/Laplacian, and one Euler row per residue give exact toy variational data; incidence rank is `8`, Hessian rank/nullity is `8/4`, and the four fixed attractors are formal boundary components.\n- The current route is bounded no-go: no unit-bearing square-map measure, strict field provenance, boundary/integration theorem, named square-map density theorem, or nonproxy continuum lift is exported.\n- Do not promote square-map action receivers, Hessian ranks, source-coupling rows, provenance/localizer rows, symbolic `L_total` slots, selector replay, bridge maps, role transfer, nonproxy `L_total`, or ToE.  The square-map-functional-graph lane is now bounded no-go on current artifacts; a next move must introduce a genuinely new strict typed object/theorem/provider outside this lane or preserve the P2929-P3005 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
