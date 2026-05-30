#!/usr/bin/env python3
"""Scratch probe: phase-sign path cohomology triviality certificate.

The Z2 phase-sign edge bits live on the finite path graph 0--1--...--11.  The
previous GF(2) certificates solve the prefix system and verify carrier/edge/node
commutative diagrams.  This probe records the graph-theoretic reason there is no
hidden cycle ambiguity: for a path with 12 vertices, 11 edges, and one connected
component, the GF(2) coboundary map

    delta: C^0 -> C^1,   (delta b)(i,i+1)=b(i) xor b(i+1)

has rank 11 and kernel consisting only of constant node cochains.  Therefore
H^1(path;GF(2)) has dimension 0, every edge cochain is exact, and the left anchor
b(0)=0 fixes the unique node potential.  This is finite graph cohomology
bookkeeping only; it is not a new phase fit, not a bridge theorem, and not a
strict dynamical derivation.
"""
from __future__ import annotations

import json
from itertools import product
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_sign_path_cohomology_triviality_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_sign_path_cohomology_triviality_certificate_report.md"
Z2_REPORT = HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json"
GF2_REPORT = HERE / "bridge_strict_completion_phase_sign_gf2_linear_system_certificate_report.json"
EDGE_SUPPORT_REPORT = HERE / "bridge_strict_completion_phase_sign_edge_support_minimality_certificate_report.json"
DIAGRAM_REPORT = HERE / "bridge_strict_completion_phase_zero_gf2_commutative_diagram_certificate_report.json"

NODE_COUNT = 12
EDGE_COUNT = 11
EXPECTED_NODE_BITS = [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0]
EXPECTED_EDGE_BITS = [0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0]
EXPECTED_FLIP_EDGES = ["1->2", "5->6", "7->8", "9->10"]
EDGE_LABELS = [f"{d}->{d + 1}" for d in range(EDGE_COUNT)]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def coboundary_matrix() -> list[list[int]]:
    rows = []
    for edge_index in range(EDGE_COUNT):
        row = [0] * NODE_COUNT
        row[edge_index] = 1
        row[edge_index + 1] = 1
        rows.append(row)
    return rows


def matvec_gf2(matrix: list[list[int]], vector: list[int]) -> list[int]:
    return [sum(value & vector[col] for col, value in enumerate(row)) % 2 for row in matrix]


def gf2_rref(matrix: list[list[int]]) -> dict[str, Any]:
    work = [row[:] for row in matrix]
    row_count = len(work)
    col_count = len(work[0]) if work else 0
    pivot_row = 0
    pivots = []
    for col in range(col_count):
        pivot = None
        for candidate in range(pivot_row, row_count):
            if work[candidate][col]:
                pivot = candidate
                break
        if pivot is None:
            continue
        if pivot != pivot_row:
            work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        for target in range(row_count):
            if target != pivot_row and work[target][col]:
                work[target] = [a ^ b for a, b in zip(work[target], work[pivot_row])]
        pivots.append({"row": pivot_row, "pivot_col": col})
        pivot_row += 1
        if pivot_row == row_count:
            break
    return {
        "rank": len(pivots),
        "nullity": col_count - len(pivots),
        "pivot_rows": pivots,
        "rref_matrix": work,
    }


def kernel_vectors(matrix: list[list[int]]) -> list[list[int]]:
    vectors = []
    for candidate in product([0, 1], repeat=NODE_COUNT):
        vector = list(candidate)
        if matvec_gf2(matrix, vector) == [0] * EDGE_COUNT:
            vectors.append(vector)
    return vectors


def anchored_reconstruction(edge_bits: list[int], anchor: int = 0) -> list[int]:
    nodes = [anchor]
    current = anchor
    for bit in edge_bits:
        current ^= bit
        nodes.append(current)
    return nodes


def support(edge_bits: list[int]) -> list[str]:
    return [edge for edge, bit in zip(EDGE_LABELS, edge_bits) if bit]


def edge_exactness_rows(coboundary: list[list[int]], node_bits: list[int], edge_bits: list[int]) -> list[dict[str, Any]]:
    computed = matvec_gf2(coboundary, node_bits)
    return [
        {
            "edge": edge,
            "delta_node_bit": got,
            "audited_edge_bit": want,
            "is_flip_edge": got == 1,
            "matches_audited_edge_bit": got == want,
        }
        for edge, got, want in zip(EDGE_LABELS, computed, edge_bits)
    ]


def build_payload() -> dict[str, Any]:
    z2 = load_json(Z2_REPORT)
    gf2 = load_json(GF2_REPORT)
    edge_support = load_json(EDGE_SUPPORT_REPORT)
    diagram = load_json(DIAGRAM_REPORT)

    node_bits = [row["node_bit"] for row in z2["node_bit_rows"]]
    edge_bits = [row["edge_bit"] for row in z2["edge_bit_rows"]]
    coboundary = coboundary_matrix()
    rref = gf2_rref(coboundary)
    kernels = kernel_vectors(coboundary)
    anchored_nodes = anchored_reconstruction(edge_bits, anchor=0)
    exact_rows = edge_exactness_rows(coboundary, node_bits, edge_bits)
    h1_dimension = EDGE_COUNT - rref["rank"]

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_PATH_COHOMOLOGY_TRIVIALITY_CERTIFICATE__GF2_NO_CYCLE_AMBIGUITY",
        "status": "finite-path-gf2-h1-zero-anchor-fixes-unique-node-potential",
        "source_reports": {
            "phase_sign_z2_coboundary_certificate": str(Z2_REPORT.relative_to(ROOT)),
            "phase_sign_gf2_linear_system_certificate": str(GF2_REPORT.relative_to(ROOT)),
            "phase_sign_edge_support_minimality_certificate": str(EDGE_SUPPORT_REPORT.relative_to(ROOT)),
            "phase_zero_gf2_commutative_diagram_certificate": str(DIAGRAM_REPORT.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "cohomology",
                "Betti",
                "cycle space",
                "H^1",
                "rank boundary",
                "gauge ambiguity",
                "path cohomology",
            ],
            "conclusion": "No strict-completion report exported a GF(2) path-cohomology triviality certificate for the phase-sign edge cochain before this file.",
        },
        "cochain_complex_definition": {
            "field": "GF(2)",
            "node_count_dim_C0": NODE_COUNT,
            "edge_count_dim_C1": EDGE_COUNT,
            "connected_components": 1,
            "cycle_rank_E_minus_V_plus_components": EDGE_COUNT - NODE_COUNT + 1,
            "coboundary_delta_matrix": coboundary,
            "left_anchor_node_bit_b0": 0,
        },
        "coboundary_rank_certificate": rref,
        "kernel_certificate": {
            "kernel_vectors": kernels,
            "kernel_size": len(kernels),
            "kernel_is_exactly_constant_node_cochains": kernels == [[0] * NODE_COUNT, [1] * NODE_COUNT],
        },
        "edge_exactness_rows": exact_rows,
        "path_cohomology_summary": {
            "rank_delta": rref["rank"],
            "nullity_delta": rref["nullity"],
            "h1_dimension_dim_C1_minus_rank_delta": h1_dimension,
            "cycle_rank": EDGE_COUNT - NODE_COUNT + 1,
            "every_edge_cochain_exact_on_path": h1_dimension == 0,
            "anchor_kills_constant_kernel": node_bits[0] == 0,
            "anchored_reconstructed_node_bits": anchored_nodes,
            "anchored_reconstruction_matches_z2_node_bits": anchored_nodes == node_bits,
            "delta_node_bits_match_edge_bits": matvec_gf2(coboundary, node_bits) == edge_bits,
            "flip_edges_from_delta": support(matvec_gf2(coboundary, node_bits)),
            "matches_expected_node_bits": node_bits == EXPECTED_NODE_BITS,
            "matches_expected_edge_bits": edge_bits == EXPECTED_EDGE_BITS,
            "matches_expected_flip_edges": support(edge_bits) == EXPECTED_FLIP_EDGES,
            "matches_gf2_linear_system_solution": edge_bits == gf2["gf2_linear_system_summary"]["solution_edge_bit_pattern"],
            "inherits_edge_support_unique_minimality": edge_support["edge_support_minimality_summary"]["unique_matching_assignment"] and edge_support["edge_support_minimality_summary"]["all_lower_support_assignments_fail"],
            "inherits_commutative_diagram": diagram["diagram_summary"]["all_diagram_checks_pass"],
        },
        "blocker_context": {
            "what_this_certifies": "finite GF(2) path cohomology has no cycle ambiguity for the audited phase-sign edge cochain",
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_transport_derivation_from_nadsoliton_dynamics",
                "strict_damping_parameter_derivation_from_nadsoliton_dynamics",
                "orientation_chi11_source",
                "QW-2191_selector_obstruction",
            ],
        },
        "proof_certificate": {
            "rank_step": "The path coboundary delta:C0->C1 has rank 11 over GF(2).",
            "kernel_step": "The kernel of delta consists exactly of the two constant node cochains, so the only gauge freedom is global sign flip before anchoring.",
            "h1_step": "Since dim C1=11 and rank(delta)=11, H^1(path;GF(2)) has dimension 0 and every edge cochain is exact.",
            "anchor_step": "The anchor b(0)=0 kills the constant-kernel ambiguity and gives the unique audited node-bit potential.",
            "theoretical_limit": "This proves finite path-cohomology triviality only; it does not derive omega/phi, beta/eta, damping, or transport from strict nadsoliton dynamics.",
        },
        "hard_limits": [
            "K_strict_gate remains the current live/full operational kernel.",
            "No unqualified identity K_legacy_ont == K_strict_gate is claimed.",
            "No proof derives A(d), P(d), D(d), omega/phi, beta/eta, or the transport cocycle from strict nadsoliton dynamics.",
            "No beta_tors -> chi_11 theorem is claimed.",
            "No legacy physical-role transfer to K_strict_gate is claimed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    summary = payload["path_cohomology_summary"]
    lines = [
        "# Phase-sign path cohomology triviality certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This certificate treats the audited phase edge bits as a GF(2) cochain on",
        "the finite path graph.  The coboundary map has rank 11, H^1 is zero, and",
        "the left anchor fixes the unique node potential.",
        "",
        "## Summary",
        "",
        f"- Rank(delta): `{summary['rank_delta']}`",
        f"- Nullity(delta): `{summary['nullity_delta']}`",
        f"- H1 dimension: `{summary['h1_dimension_dim_C1_minus_rank_delta']}`",
        f"- Cycle rank: `{summary['cycle_rank']}`",
        f"- Every edge cochain exact on path: `{summary['every_edge_cochain_exact_on_path']}`",
        f"- Anchor kills constant kernel: `{summary['anchor_kills_constant_kernel']}`",
        f"- Flip edges from delta: `{summary['flip_edges_from_delta']}`",
        f"- Inherits edge-support uniqueness/minimality: `{summary['inherits_edge_support_unique_minimality']}`",
        f"- Inherits commutative diagram: `{summary['inherits_commutative_diagram']}`",
        "",
        "## Hard limits",
        "",
    ]
    for limit in payload["hard_limits"]:
        lines.append(f"- {limit}")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
