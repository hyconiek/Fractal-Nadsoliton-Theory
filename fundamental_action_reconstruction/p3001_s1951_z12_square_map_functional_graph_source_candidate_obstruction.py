#!/usr/bin/env python3
"""P3001/S1951: Z/12Z square-map functional graph source-candidate obstruction.

P3000 bounded the zero-divisor-graph lane.  This audit introduces one genuinely
new finite typed object outside the closed zero-divisor/Fourier/annihilator/
nilradical/CRT/zero-derivation lanes: the functional graph of the squaring map
x -> x^2 mod 12 on all residues of Z/12Z.

The finite calculation enumerates every directed edge, fixed-point attractor,
basin, fiber, and depth.  It gives a real exact idempotent-attractor dynamics:
all nodes reach one of four fixed idempotents in at most one step.  The object
remains a source-candidate obstruction because no current artifact exports a
strict nadsoliton provenance map, nonpremise basin/attractor localizer, named
source-atom coupling, unit-bearing action installation, or nonproxy continuum
lift for this square-map graph.
"""
from __future__ import annotations

import hashlib, json
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3000_s1950_zero_divisor_graph_action_installation_obstruction import OUT as P3000

OUT = GEN / "p3001_s1951_z12_square_map_functional_graph_source_candidate_obstruction.json"
MD = GEN / "p3001_s1951_z12_square_map_functional_graph_source_candidate_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
MODULUS = 12


def square_map(x: int) -> int:
    return (x * x) % MODULUS


def orbit_from(x: int) -> dict[str, Any]:
    seen: dict[int, int] = {}
    path: list[int] = []
    y = x
    while y not in seen:
        seen[y] = len(path)
        path.append(y)
        y = square_map(y)
    cycle = path[seen[y]:]
    return {
        "start": x,
        "path": path,
        "cycle": cycle,
        "attractor": cycle[0] if len(cycle) == 1 else cycle,
        "preperiod": seen[y],
        "period": len(cycle),
        "depth_to_attractor": seen[y],
    }


def square_map_witness() -> dict[str, Any]:
    residues = list(range(MODULUS))
    edges = [{"source": x, "target": square_map(x)} for x in residues]
    image = sorted({e["target"] for e in edges})
    fixed_points = [x for x in residues if square_map(x) == x]
    orbit_rows = [orbit_from(x) for x in residues]
    basins = {str(fp): sorted([row["start"] for row in orbit_rows if row["attractor"] == fp]) for fp in fixed_points}
    fibers = {str(y): sorted([x for x in residues if square_map(x) == y]) for y in image}
    nontrivial_fibers = {k: v for k, v in fibers.items() if len(v) > 1}
    row_data = []
    for row in orbit_rows:
        x = row["start"]
        row_data.append({
            "residue": x,
            "square_target": square_map(x),
            "path": row["path"],
            "fixed_attractor": row["attractor"],
            "depth_to_attractor": row["depth_to_attractor"],
            "period": row["period"],
            "strict_nadsoliton_provenance": False,
            "nonpremise_basin_attractor_localizer": False,
            "named_source_atom_coupling": False,
            "unit_bearing_action_installation": False,
            "accepted_strict_source": False,
        })
    return {
        "modulus": MODULUS,
        "node_count": len(residues),
        "edge_count": len(edges),
        "edges": edges,
        "image": image,
        "image_size": len(image),
        "fixed_points": fixed_points,
        "fixed_point_count": len(fixed_points),
        "all_fixed_points_idempotent": all(square_map(x) == x for x in fixed_points),
        "basins": basins,
        "basin_sizes": {k: len(v) for k, v in basins.items()},
        "fibers": fibers,
        "nontrivial_fibers": nontrivial_fibers,
        "max_depth_to_attractor": max(row["depth_to_attractor"] for row in orbit_rows),
        "all_periods_one": all(row["period"] == 1 for row in orbit_rows),
        "row_data": row_data,
        "accepted_strict_sources": [r["residue"] for r in row_data if r["accepted_strict_source"]],
    }


def obligation_rows(witness: dict[str, Any]) -> list[dict[str, Any]]:
    rows = witness["row_data"]
    return [
        {"obligation": "finite_square_map_functional_graph", "satisfied": witness["node_count"] == 12 and witness["edge_count"] == 12, "evidence": "one directed squaring edge x -> x^2 mod 12 is built for each residue"},
        {"obligation": "idempotent_attractor_decomposition", "satisfied": witness["fixed_points"] == [0, 1, 4, 9] and witness["max_depth_to_attractor"] == 1 and witness["all_periods_one"], "evidence": "all residues reach one of four fixed idempotents in at most one step"},
        {"obligation": "fiber_basin_certificate", "satisfied": witness["basin_sizes"] == {"0": 2, "1": 4, "4": 4, "9": 2}, "evidence": f"basin sizes are {witness['basin_sizes']} with image {witness['image']}"},
        {"obligation": "strict_nadsoliton_provenance", "satisfied": any(r["strict_nadsoliton_provenance"] for r in rows), "evidence": "the square map is imported finite ring dynamics, not an exported nadsoliton source map"},
        {"obligation": "nonpremise_basin_attractor_localizer", "satisfied": any(r["nonpremise_basin_attractor_localizer"] for r in rows), "evidence": "basin/attractor labels are algebraic dynamics labels, not physical sectors"},
        {"obligation": "named_source_atom_coupling", "satisfied": any(r["named_source_atom_coupling"] for r in rows), "evidence": "no theorem couples a square-map attractor or basin to selector sign, beta/Z_beta, bridge-source, or action-density atoms"},
        {"obligation": "unit_bearing_action_installation", "satisfied": any(r["unit_bearing_action_installation"] for r in rows), "evidence": "finite directed dynamics has no unit-bearing measure/density theorem or nonproxy lift"},
        {"obligation": "accepted_current_strict_source", "satisfied": bool(witness["accepted_strict_sources"]), "evidence": "no current row satisfies provenance, localization, coupling, and action-installation premises"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_graph", "idempotent_attractors", "fiber_basin_certificate", "strict_provenance", "nonpremise_localizer", "source_atom_coupling", "unit_action_installation", "nonproxy_export"]
    return [{"present": dict(zip(names, bits)), "accepts_square_map_strict_source": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p3000_path: Any) -> dict[str, Any]:
    witness = square_map_witness()
    matrix = acceptance_matrix()
    return {
        "status": "P3001_Z12_SQUARE_MAP_FUNCTIONAL_GRAPH_SOURCE_CANDIDATE_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P3000": hashlib.sha256(p3000_path.read_bytes()).hexdigest() if p3000_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "Z12SquareMapFunctionalGraph_SourceCandidateObstructionMatrix",
            "square_map_witness": witness,
            "proof_obligation_rows": obligation_rows(witness),
            "finite_acceptance_matrix": matrix,
        },
        "functional_graph_certificate": {
            "node_count": witness["node_count"],
            "edge_count": witness["edge_count"],
            "image": witness["image"],
            "image_size": witness["image_size"],
            "fixed_points": witness["fixed_points"],
            "fixed_point_count": witness["fixed_point_count"],
            "basin_sizes": witness["basin_sizes"],
            "max_depth_to_attractor": witness["max_depth_to_attractor"],
            "all_periods_one": witness["all_periods_one"],
            "accepted_strict_sources": witness["accepted_strict_sources"],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_square_map_strict_source"]),
        },
        "decision": {
            "positive_progress": "P3001 introduces one new finite typed object after the P3000 no-go: the Z/12Z square-map functional graph x -> x^2 mod 12.",
            "breakthrough": "Bounded no-go: the exact finite dynamics collapses all 12 residues onto four fixed idempotent attractors in at most one step, but no strict provenance, nonpremise basin/attractor localizer, named source-atom coupling, unit-bearing action installation, or nonproxy export is present.",
            "negative_export_flags": {k: False for k in ["square_map_strict_source_exported", "strict_nadsoliton_provenance_exported", "nonpremise_basin_attractor_localizer_exported", "named_source_atom_coupling_exported", "unit_bearing_action_installation_exported", "selector_closure_exported", "bridge_closure_exported", "nonproxy_ltotal_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "A next admissible square-map-functional-graph move may attack exactly one missing theorem for this new object: strict provenance, nonpremise basin/attractor localizer, named source-atom coupling, or unit-bearing action installation.  Do not replay zero-divisor graph, CRT idempotent, nilradical, annihilator, Fourier, selector, bridge, role-transfer, or L_total lanes through square-map labels; otherwise introduce a genuinely new strict typed object/provider or preserve the P2929-P3001 no-strict-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["functional_graph_certificate"]
    lines = [
        "# P3001/S1951 Z/12Z square-map functional graph source-candidate obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Functional graph certificate",
        f"- nodes/edges: `{cert['node_count']}/{cert['edge_count']}`",
        f"- image: `{cert['image']}`",
        f"- fixed points: `{cert['fixed_points']}`",
        f"- basin sizes: `{cert['basin_sizes']}`",
        f"- max depth to attractor: `{cert['max_depth_to_attractor']}`",
        f"- all periods one: `{cert['all_periods_one']}`",
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
    read_json(P3000)
    payload = build_payload(P3000)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3001/S1951 Z12 square-map functional graph source-candidate obstruction", "## P3001/S1951 Z12 square-map functional graph source-candidate obstruction\n\n`P3001/S1951` introduces one new finite typed object after the P3000 zero-divisor-graph no-go: the functional graph of the squaring map `x -> x^2 mod 12` on all residues of `Z/12Z`.  The finite dynamics is exact: `12` directed edges collapse onto image `[0, 1, 4, 9]`, all four image points are fixed idempotent attractors, basin sizes are `{0:2, 1:4, 4:4, 9:2}`, and every node reaches its attractor in at most one step.  This is real idempotent-attractor dynamics, but no strict nadsoliton provenance, nonpremise basin/attractor localizer, named source-atom coupling, unit-bearing action installation, nonproxy `L_total`, bridge closure, role transfer, or ToE is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3001/S1951 square-map functional graph `L_total` guard", "## P3001/S1951 square-map functional graph `L_total` guard\n\n`P3001/S1951` adds no square-map functional graph term to `L_total`.  The idempotent-attractor and basin/fiber certificates are finite ring-dynamics labels only; they do not supply strict field provenance, unit-bearing coefficient, named density theorem, variational chain, EOM/Hamiltonian term, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current Z12 square-map functional graph source-candidate obstruction guardrail (P3001/S1951, 2026-06-20)", "## Current Z12 square-map functional graph source-candidate obstruction guardrail (P3001/S1951, 2026-06-20)\n\n- P3001 introduces one new finite typed object after the P3000 zero-divisor-graph no-go: the functional graph of the squaring map `x -> x^2 mod 12` on all residues of `Z/12Z`.\n- Exact finite positives are real but ring-dynamical only: `12` directed edges collapse onto fixed idempotent attractors `[0, 1, 4, 9]`, basin sizes are `{0:2, 1:4, 4:4, 9:2}`, and every node reaches its attractor in at most one step; the acceptance matrix has `256` profiles with only the full profile accepting.\n- The current route is bounded no-go as a strict source candidate: no strict nadsoliton provenance, nonpremise basin/attractor localizer, named source-atom coupling, unit-bearing action installation, or nonproxy export is present.\n- Do not replay zero-divisor graph, CRT idempotent, nilradical, annihilator, Fourier, selector, bridge, role-transfer, or `L_total` lanes through square-map labels.  A next admissible square-map-functional-graph move may attack exactly one missing theorem for this object, or else introduce a genuinely new strict typed object/provider while preserving the P2929-P3001 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
