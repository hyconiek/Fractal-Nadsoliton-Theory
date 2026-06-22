#!/usr/bin/env python3
"""P3006/S1956: Z/12Z cube-map functional graph source-candidate obstruction.

P3005 bounded the square-map-functional-graph lane. This audit introduces one
new finite typed object outside the closed square-map, zero-divisor, Fourier,
annihilator, nilradical, CRT, and zero-derivation replay lanes: the functional
graph of the cubing map x -> x^3 mod 12 on all residues of Z/12Z.

The finite calculation enumerates every directed edge, image point, fixed point,
basin, fiber, orbit path, period, and depth. It gives an exact finite
near-identity power-map dynamics: nine fixed residues and three one-step
preimages into fixed residues. The object remains a source-candidate obstruction
because no current artifact exports strict nadsoliton provenance, a nonpremise
basin/fixed-sector localizer, named source-atom coupling, unit-bearing action
installation, or nonproxy continuum lift for this cube-map graph.
"""
from __future__ import annotations

import hashlib, json
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3005_s1955_square_map_action_installation_obstruction import OUT as P3005

OUT = GEN / "p3006_s1956_z12_cube_map_functional_graph_source_candidate_obstruction.json"
MD = GEN / "p3006_s1956_z12_cube_map_functional_graph_source_candidate_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
MODULUS = 12


def cube_map(x: int) -> int:
    return (x * x * x) % MODULUS


def orbit_from(x: int) -> dict[str, Any]:
    seen: dict[int, int] = {}
    path: list[int] = []
    y = x
    while y not in seen:
        seen[y] = len(path)
        path.append(y)
        y = cube_map(y)
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


def cube_map_witness() -> dict[str, Any]:
    residues = list(range(MODULUS))
    edges = [{"source": x, "target": cube_map(x)} for x in residues]
    image = sorted({edge["target"] for edge in edges})
    fixed_points = [x for x in residues if cube_map(x) == x]
    orbit_rows = [orbit_from(x) for x in residues]
    basins = {str(fp): sorted([row["start"] for row in orbit_rows if row["attractor"] == fp]) for fp in fixed_points}
    fibers = {str(y): sorted([x for x in residues if cube_map(x) == y]) for y in image}
    nontrivial_fibers = {k: v for k, v in fibers.items() if len(v) > 1}
    moved_points = [x for x in residues if cube_map(x) != x]
    row_data = []
    for row in orbit_rows:
        x = row["start"]
        row_data.append({
            "residue": x,
            "cube_target": cube_map(x),
            "path": row["path"],
            "fixed_attractor": row["attractor"],
            "depth_to_attractor": row["depth_to_attractor"],
            "period": row["period"],
            "is_fixed": cube_map(x) == x,
            "strict_nadsoliton_provenance": False,
            "nonpremise_basin_fixed_sector_localizer": False,
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
        "moved_points": moved_points,
        "moved_point_count": len(moved_points),
        "basins": basins,
        "basin_sizes": {k: len(v) for k, v in basins.items()},
        "fibers": fibers,
        "nontrivial_fibers": nontrivial_fibers,
        "max_depth_to_attractor": max(row["depth_to_attractor"] for row in orbit_rows),
        "all_periods_one": all(row["period"] == 1 for row in orbit_rows),
        "row_data": row_data,
        "accepted_strict_sources": [row["residue"] for row in row_data if row["accepted_strict_source"]],
    }


def obligation_rows(witness: dict[str, Any]) -> list[dict[str, Any]]:
    rows = witness["row_data"]
    return [
        {"obligation": "finite_cube_map_functional_graph", "satisfied": witness["node_count"] == 12 and witness["edge_count"] == 12, "evidence": "one directed cubing edge x -> x^3 mod 12 is built for each residue"},
        {"obligation": "fixed_point_decomposition", "satisfied": witness["fixed_points"] == [0, 1, 3, 4, 5, 7, 8, 9, 11] and witness["moved_points"] == [2, 6, 10] and witness["max_depth_to_attractor"] == 1 and witness["all_periods_one"], "evidence": "nine residues are fixed and three residues flow to fixed points in one step"},
        {"obligation": "fiber_basin_certificate", "satisfied": witness["basin_sizes"] == {"0": 2, "1": 1, "3": 1, "4": 2, "5": 1, "7": 1, "8": 2, "9": 1, "11": 1}, "evidence": f"basin sizes are {witness['basin_sizes']} with image {witness['image']}"},
        {"obligation": "strict_nadsoliton_provenance", "satisfied": any(row["strict_nadsoliton_provenance"] for row in rows), "evidence": "the cube map is imported finite ring power dynamics, not an exported nadsoliton source map"},
        {"obligation": "nonpremise_basin_fixed_sector_localizer", "satisfied": any(row["nonpremise_basin_fixed_sector_localizer"] for row in rows), "evidence": "cube-map fixed/basin labels are algebraic dynamics labels, not physical sectors"},
        {"obligation": "named_source_atom_coupling", "satisfied": any(row["named_source_atom_coupling"] for row in rows), "evidence": "no theorem couples a cube-map fixed point or basin to selector sign, beta/Z_beta, bridge-source, or action-density atoms"},
        {"obligation": "unit_bearing_action_installation", "satisfied": any(row["unit_bearing_action_installation"] for row in rows), "evidence": "finite directed cubing dynamics has no unit-bearing measure/density theorem or nonproxy lift"},
        {"obligation": "accepted_current_strict_source", "satisfied": bool(witness["accepted_strict_sources"]), "evidence": "no current row satisfies provenance, localization, coupling, and action-installation premises"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_graph", "fixed_point_decomposition", "fiber_basin_certificate", "strict_provenance", "nonpremise_localizer", "source_atom_coupling", "unit_action_installation", "nonproxy_export"]
    return [{"present": dict(zip(names, bits)), "accepts_cube_map_strict_source": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p3005_path: Any) -> dict[str, Any]:
    witness = cube_map_witness()
    matrix = acceptance_matrix()
    return {
        "status": "P3006_Z12_CUBE_MAP_FUNCTIONAL_GRAPH_SOURCE_CANDIDATE_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P3005": hashlib.sha256(p3005_path.read_bytes()).hexdigest() if p3005_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "Z12CubeMapFunctionalGraph_SourceCandidateObstructionMatrix",
            "cube_map_witness": witness,
            "proof_obligation_rows": obligation_rows(witness),
            "finite_acceptance_matrix": matrix,
        },
        "cube_map_certificate": {
            "node_count": witness["node_count"],
            "edge_count": witness["edge_count"],
            "image": witness["image"],
            "image_size": witness["image_size"],
            "fixed_points": witness["fixed_points"],
            "fixed_point_count": witness["fixed_point_count"],
            "moved_points": witness["moved_points"],
            "moved_point_count": witness["moved_point_count"],
            "basin_sizes": witness["basin_sizes"],
            "nontrivial_fibers": witness["nontrivial_fibers"],
            "max_depth_to_attractor": witness["max_depth_to_attractor"],
            "all_periods_one": witness["all_periods_one"],
            "accepted_strict_sources": witness["accepted_strict_sources"],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for row in matrix if row["accepts_cube_map_strict_source"]),
        },
        "decision": {
            "positive_progress": "P3006 introduces one new finite typed object after the P3005 square-map no-go: the functional graph of x -> x^3 mod 12 on all residues of Z/12Z.",
            "breakthrough": "Bounded no-go: exact finite cubing dynamics has nine fixed residues and three one-step preimages, but no strict provenance, nonpremise basin/fixed-sector localizer, named source-atom coupling, unit-bearing action installation, or nonproxy export is present.",
            "negative_export_flags": {k: False for k in ["cube_map_strict_source_exported", "strict_nadsoliton_provenance_exported", "nonpremise_basin_fixed_sector_localizer_exported", "named_source_atom_coupling_exported", "unit_bearing_action_installation_exported", "selector_closure_exported", "bridge_closure_exported", "nonproxy_ltotal_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "A next admissible cube-map-functional-graph move may attack exactly one missing theorem for this new object: strict provenance, nonpremise basin/fixed-sector localizer, named source-atom coupling, or unit-bearing action installation. Do not replay square-map, zero-divisor graph, CRT idempotent, nilradical, annihilator, Fourier, selector, bridge, role-transfer, or L_total lanes through cube-map labels; otherwise introduce a genuinely new strict typed object/provider or preserve the P2929-P3006 no-strict-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["cube_map_certificate"]
    lines = [
        "# P3006/S1956 Z12 cube-map functional graph source-candidate obstruction", "",
        f"Status: `{payload['status']}`", "", "## Cube-map certificate",
        f"- nodes/edges: `{cert['node_count']}/{cert['edge_count']}`",
        f"- image: `{cert['image']}`",
        f"- fixed points: `{cert['fixed_points']}`",
        f"- moved points: `{cert['moved_points']}`",
        f"- basin sizes: `{cert['basin_sizes']}`",
        f"- nontrivial fibers: `{cert['nontrivial_fibers']}`",
        f"- max depth to attractor: `{cert['max_depth_to_attractor']}`",
        f"- all periods one: `{cert['all_periods_one']}`",
        f"- accepted strict sources: `{cert['accepted_strict_sources']}`",
        f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_rows']}`", "",
        "## Lay summary", payload["decision"]["positive_progress"], payload["decision"]["breakthrough"], "",
        "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P3005)
    payload = build_payload(P3005)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3006/S1956 Z12 cube-map functional graph source-candidate obstruction", "## P3006/S1956 Z12 cube-map functional graph source-candidate obstruction\n\n`P3006/S1956` introduces one new finite typed object after the P3005 square-map no-go: the functional graph of the cubing map `x -> x^3 mod 12` on all residues of `Z/12Z`.  The finite dynamics is exact: `12` directed edges have image `[0, 1, 3, 4, 5, 7, 8, 9, 11]`, nine fixed residues, moved points `[2, 6, 10]`, basin sizes `{0:2, 1:1, 3:1, 4:2, 5:1, 7:1, 8:2, 9:1, 11:1}`, and every residue reaches a fixed point in at most one step.  This is real finite power-map dynamics, but no strict nadsoliton provenance, nonpremise basin/fixed-sector localizer, named source-atom coupling, unit-bearing action installation, nonproxy `L_total`, bridge closure, role transfer, or ToE is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3006/S1956 cube-map functional graph `L_total` guard", "## P3006/S1956 cube-map functional graph `L_total` guard\n\n`P3006/S1956` adds no cube-map term to `L_total`.  The finite cubing graph, fixed residues, moved residues, basins, fibers, and acceptance matrix do not supply strict field provenance, a named density theorem, unit-bearing measure, variational chain, EOM/Hamiltonian term, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current Z12 cube-map functional graph source-candidate obstruction guardrail (P3006/S1956, 2026-06-22)", "## Current Z12 cube-map functional graph source-candidate obstruction guardrail (P3006/S1956, 2026-06-22)\n\n- P3006 introduces one new finite typed object after the P3005 square-map no-go: the functional graph of the cubing map `x -> x^3 mod 12` on all residues of `Z/12Z`.\n- Exact finite positives are real but ring-dynamical only: `12` directed edges have image `[0, 1, 3, 4, 5, 7, 8, 9, 11]`, nine fixed residues, moved points `[2, 6, 10]`, and every residue reaches a fixed point in at most one step; the acceptance matrix has `256` profiles with only the full profile accepting.\n- The current route is bounded no-go as a strict source candidate: no strict nadsoliton provenance, nonpremise basin/fixed-sector localizer, named source-atom coupling, unit-bearing action installation, or nonproxy export is present.\n- Do not replay square-map, zero-divisor graph, CRT idempotent, nilradical, annihilator, Fourier, selector, bridge, role-transfer, or `L_total` lanes through cube-map labels.  A next admissible cube-map-functional-graph move may attack exactly one missing theorem for this object, or else introduce a genuinely new strict typed object/provider while preserving the P2929-P3006 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
