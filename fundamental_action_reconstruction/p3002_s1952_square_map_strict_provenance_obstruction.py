#!/usr/bin/env python3
"""P3002/S1952: square-map functional graph strict provenance obstruction.

P3001 left four square-map-functional-graph theorem routes.  This audit attacks
exactly one: strict nadsoliton provenance for the Z/12Z squaring functional
graph.  It does not replay basin/localizer, source-atom coupling, action
installation, zero-divisor graph, CRT idempotent, nilradical, annihilator,
Fourier, selector, bridge, role-transfer, or L_total lanes.

The finite calculation proves strong algebraic provenance receivers: the square
map is a multiplicative endomap on all 144 ordered pairs, all U(12) unit actions
are square-invariant because u^2=1 mod 12, and no nonzero translation preserves
the directed functional graph.  These facts are exact but still imported ring
algebra; no current artifact exports a strict nadsoliton source map that sources
this squaring dynamics internally.
"""
from __future__ import annotations

import hashlib, json
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3001_s1951_z12_square_map_functional_graph_source_candidate_obstruction import OUT as P3001, MODULUS, square_map, square_map_witness

OUT = GEN / "p3002_s1952_square_map_strict_provenance_obstruction.json"
MD = GEN / "p3002_s1952_square_map_strict_provenance_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
UNITS = [1, 5, 7, 11]


def provenance_witness() -> dict[str, Any]:
    residues = list(range(MODULUS))
    graph = square_map_witness()
    multiplicative_rows = []
    additive_defect_rows = []
    for x, y in product(residues, repeat=2):
        lhs_mult = square_map((x * y) % MODULUS)
        rhs_mult = (square_map(x) * square_map(y)) % MODULUS
        multiplicative_rows.append({"x": x, "y": y, "lhs": lhs_mult, "rhs": rhs_mult, "compatible": lhs_mult == rhs_mult})
        lhs_add = square_map((x + y) % MODULUS)
        rhs_add = (square_map(x) + square_map(y)) % MODULUS
        if lhs_add != rhs_add:
            additive_defect_rows.append({"x": x, "y": y, "lhs": lhs_add, "rhs": rhs_add, "defect": (lhs_add - rhs_add) % MODULUS})
    unit_rows = []
    for u, x in product(UNITS, residues):
        unit_rows.append({
            "unit": u,
            "x": x,
            "u_square_mod_12": square_map(u),
            "square_of_unit_times_x": square_map((u * x) % MODULUS),
            "square_of_x": square_map(x),
            "unit_square_invariant": square_map((u * x) % MODULUS) == square_map(x),
        })
    translation_rows = []
    edge_set = {(edge["source"], edge["target"]) for edge in graph["edges"]}
    for t in residues:
        translated = {((source + t) % MODULUS, (target + t) % MODULUS) for source, target in edge_set}
        translation_rows.append({"translation": t, "preserves_directed_graph": translated == edge_set})
    return {
        "node_count": graph["node_count"],
        "edge_count": graph["edge_count"],
        "multiplicative_pair_count": len(multiplicative_rows),
        "multiplicative_rows": multiplicative_rows,
        "all_multiplicative_pairs_compatible": all(row["compatible"] for row in multiplicative_rows),
        "additive_defect_count": len(additive_defect_rows),
        "additive_defect_sample": additive_defect_rows[:24],
        "unit_action_row_count": len(unit_rows),
        "unit_rows": unit_rows,
        "all_unit_actions_square_invariant": all(row["unit_square_invariant"] for row in unit_rows),
        "graph_preserving_translations": [row["translation"] for row in translation_rows if row["preserves_directed_graph"]],
        "translation_rows": translation_rows,
        "strict_nadsoliton_source_map_exported": False,
        "nonpremise_internal_square_law_exported": False,
        "accepted_strict_provenance": False,
    }


def obligation_rows(witness: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {"obligation": "finite_square_graph_recomputed", "satisfied": witness["node_count"] == 12 and witness["edge_count"] == 12, "evidence": "P3001 square-map functional graph is recomputed"},
        {"obligation": "multiplicative_endomap_verified", "satisfied": witness["all_multiplicative_pairs_compatible"] and witness["multiplicative_pair_count"] == 144, "evidence": "(xy)^2 = x^2 y^2 mod 12 on all ordered residue pairs"},
        {"obligation": "unit_square_invariance_verified", "satisfied": witness["all_unit_actions_square_invariant"] and witness["unit_action_row_count"] == 48, "evidence": "all units have u^2=1 mod 12, hence square(ux)=square(x)"},
        {"obligation": "translation_nonprovenance_witness", "satisfied": witness["graph_preserving_translations"] == [0], "evidence": "only the zero translation preserves the directed functional graph"},
        {"obligation": "strict_nadsoliton_source_map_exported", "satisfied": witness["strict_nadsoliton_source_map_exported"], "evidence": "multiplicative compatibility is imported ring algebra, not a strict nadsoliton source map"},
        {"obligation": "nonpremise_internal_square_law_exported", "satisfied": witness["nonpremise_internal_square_law_exported"], "evidence": "no current artifact exports squaring as an internal nadsoliton update law"},
        {"obligation": "accepted_current_strict_provenance", "satisfied": witness["accepted_strict_provenance"], "evidence": "finite provenance receivers lack strict source export"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_graph", "multiplicative_endomap", "unit_square_invariance", "translation_witness", "strict_source_map", "internal_square_law"]
    return [{"present": dict(zip(names, bits)), "accepts_strict_provenance": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p3001_path: object) -> dict[str, Any]:
    witness = provenance_witness()
    matrix = acceptance_matrix()
    return {
        "status": "P3002_SQUARE_MAP_STRICT_PROVENANCE_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P3001": hashlib.sha256(p3001_path.read_bytes()).hexdigest() if p3001_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "SquareMapStrictProvenance_ObstructionMatrix",
            "provenance_witness": witness,
            "proof_obligation_rows": obligation_rows(witness),
            "finite_acceptance_matrix": matrix,
        },
        "provenance_certificate": {
            "node_count": witness["node_count"],
            "edge_count": witness["edge_count"],
            "multiplicative_pair_count": witness["multiplicative_pair_count"],
            "all_multiplicative_pairs_compatible": witness["all_multiplicative_pairs_compatible"],
            "additive_defect_count": witness["additive_defect_count"],
            "unit_action_row_count": witness["unit_action_row_count"],
            "all_unit_actions_square_invariant": witness["all_unit_actions_square_invariant"],
            "graph_preserving_translations": witness["graph_preserving_translations"],
            "accepted_strict_provenance": witness["accepted_strict_provenance"],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for row in matrix if row["accepts_strict_provenance"]),
        },
        "decision": {
            "positive_progress": "P3002 attacks exactly one P3001 missing theorem route: strict provenance for the Z/12Z square-map functional graph.",
            "breakthrough": "Bounded no-go: the square map has exact multiplicative-endomap and unit-square-invariance receivers, plus a translation nonprovenance witness, but no strict nadsoliton source map or nonpremise internal square law is exported.",
            "negative_export_flags": {k: False for k in ["square_map_strict_provenance_exported", "strict_nadsoliton_source_map_exported", "internal_square_update_law_exported", "basin_attractor_localizer_exported", "source_atom_coupling_exported", "unit_bearing_action_installation_exported", "selector_closure_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "A next admissible square-map-functional-graph move may attack exactly one remaining route: nonpremise basin/attractor localizer, named source-atom coupling, or unit-bearing action installation.  Do not replay multiplicative compatibility, unit-square invariance, CRT idempotents, zero-divisor graph, selector, bridge, role-transfer, or L_total lanes as provenance.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["provenance_certificate"]
    lines = [
        "# P3002/S1952 square-map strict provenance obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Provenance certificate",
        f"- nodes/edges: `{cert['node_count']}/{cert['edge_count']}`",
        f"- multiplicative pairs: `{cert['multiplicative_pair_count']}`",
        f"- all multiplicative pairs compatible: `{cert['all_multiplicative_pairs_compatible']}`",
        f"- additive defect count: `{cert['additive_defect_count']}`",
        f"- unit action rows: `{cert['unit_action_row_count']}`",
        f"- all unit actions square-invariant: `{cert['all_unit_actions_square_invariant']}`",
        f"- graph-preserving translations: `{cert['graph_preserving_translations']}`",
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
    read_json(P3001)
    payload = build_payload(P3001)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3002/S1952 square-map strict provenance obstruction", "## P3002/S1952 square-map strict provenance obstruction\n\n`P3002/S1952` attacks exactly one P3001 missing theorem route: strict provenance for the `Z/12Z` square-map functional graph.  The finite provenance receivers are exact: `(xy)^2 = x^2 y^2 mod 12` holds on all `144` ordered pairs, all `48` unit-action rows are square-invariant because every unit has square `1`, and only translation `0` preserves the directed functional graph.  These are algebraic provenance receivers only: no strict nadsoliton source map or nonpremise internal square law is exported.  No basin/attractor localizer, source-atom coupling, unit-bearing action installation, nonproxy `L_total`, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3002/S1952 square-map strict provenance `L_total` guard", "## P3002/S1952 square-map strict provenance `L_total` guard\n\n`P3002/S1952` adds no square-map provenance term to `L_total`.  Multiplicative compatibility, unit-square invariance, and translation noninvariance are finite ring-algebra receivers only; they do not supply strict field provenance, named density theorem, variational chain, EOM/Hamiltonian term, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current square-map strict provenance obstruction guardrail (P3002/S1952, 2026-06-20)", "## Current square-map strict provenance obstruction guardrail (P3002/S1952, 2026-06-20)\n\n- P3002 attacks exactly one P3001 missing theorem route: strict provenance for the `Z/12Z` square-map functional graph.\n- Finite positives are exact but algebraic only: `(xy)^2 = x^2 y^2 mod 12` holds on all `144` ordered pairs, all `48` unit-action rows are square-invariant, and only translation `0` preserves the directed functional graph.\n- The route is bounded no-go because no strict nadsoliton source map or nonpremise internal square law is exported.\n- Do not promote multiplicative compatibility, unit-square invariance, translation noninvariance, CRT idempotent replay, zero-divisor graph replay, selector replay, bridge maps, role transfer, nonproxy `L_total`, or ToE.  A next admissible square-map-functional-graph move may attack exactly one remaining route, or else introduce a genuinely new strict typed object/provider while preserving the P2929-P3002 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
