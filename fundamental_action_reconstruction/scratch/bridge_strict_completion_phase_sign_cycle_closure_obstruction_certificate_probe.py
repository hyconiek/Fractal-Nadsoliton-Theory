#!/usr/bin/env python3
"""Scratch probe: phase-sign cycle-closure obstruction certificate.

The audited phase-sign edge bits live on the path 0--1--...--11.  The path
cohomology certificate proves that every path edge cochain is exact after a left
anchor is chosen.  This probe records the adjacent boundary-condition fact: if
one artificially closes the endpoint edge 11--0, exactness is no longer automatic
because the closed 12-cycle has one GF(2) cycle functional.

For the audited node bits, the forced closing edge bit is b(11) xor b(0)=0.  The
path edge bits plus that zero closing bit have even total cycle parity and are
exact on the 12-cycle.  The odd closing-edge perturbation has total cycle parity
1 and is rejected by an exhaustive GF(2) exactness check.  This is finite graph
bookkeeping only; it is not a cyclic gate-generation strategy, not a phase fit,
not a bridge theorem, and not a strict dynamical derivation.
"""
from __future__ import annotations

import json
from itertools import product
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_sign_cycle_closure_obstruction_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_sign_cycle_closure_obstruction_certificate_report.md"
Z2_REPORT = HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json"
PATH_COHOMOLOGY_REPORT = HERE / "bridge_strict_completion_phase_sign_path_cohomology_triviality_certificate_report.json"
GF2_REPORT = HERE / "bridge_strict_completion_phase_sign_gf2_linear_system_certificate_report.json"

NODE_COUNT = 12
PATH_EDGE_COUNT = 11
CYCLE_EDGE_COUNT = 12
CLOSING_EDGE = "11->0"
PATH_EDGE_LABELS = [f"{d}->{d + 1}" for d in range(PATH_EDGE_COUNT)]
CYCLE_EDGE_LABELS = PATH_EDGE_LABELS + [CLOSING_EDGE]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def cycle_coboundary_matrix() -> list[list[int]]:
    rows = []
    for edge_index in range(PATH_EDGE_COUNT):
        row = [0] * NODE_COUNT
        row[edge_index] = 1
        row[edge_index + 1] = 1
        rows.append(row)
    closing = [0] * NODE_COUNT
    closing[-1] = 1
    closing[0] = 1
    rows.append(closing)
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
        if matvec_gf2(matrix, vector) == [0] * CYCLE_EDGE_COUNT:
            vectors.append(vector)
    return vectors


def exact_node_potentials(matrix: list[list[int]], edge_bits: list[int]) -> list[list[int]]:
    potentials = []
    for candidate in product([0, 1], repeat=NODE_COUNT):
        vector = list(candidate)
        if matvec_gf2(matrix, vector) == edge_bits:
            potentials.append(vector)
    return potentials


def anchored_reconstruction(edge_bits: list[int], anchor: int = 0) -> list[int]:
    nodes = [anchor]
    current = anchor
    for bit in edge_bits[:PATH_EDGE_COUNT]:
        current ^= bit
        nodes.append(current)
    return nodes


def closure_case_row(label: str, edge_bits: list[int], node_bits: list[int], matrix: list[list[int]]) -> dict[str, Any]:
    potentials = exact_node_potentials(matrix, edge_bits)
    anchored_nodes = anchored_reconstruction(edge_bits, anchor=node_bits[0])
    forced_closing_bit = node_bits[-1] ^ node_bits[0]
    total_cycle_parity = sum(edge_bits) % 2
    return {
        "case": label,
        "cycle_edge_bits": edge_bits,
        "closing_edge_bit": edge_bits[-1],
        "forced_closing_bit_from_audited_nodes": forced_closing_bit,
        "total_cycle_parity_mod2": total_cycle_parity,
        "cycle_parity_even": total_cycle_parity == 0,
        "exact_node_potential_count": len(potentials),
        "is_exact_cycle_coboundary": len(potentials) > 0,
        "anchored_reconstructed_node_bits": anchored_nodes,
        "anchored_reconstruction_matches_audited_nodes": anchored_nodes == node_bits,
        "closing_edge_matches_audited_endpoint_xor": edge_bits[-1] == forced_closing_bit,
        "constant_kernel_pair_present_when_exact": potentials == [] or potentials == [node_bits, [bit ^ 1 for bit in node_bits]],
    }


def edge_rows(audited_cycle_bits: list[int], odd_cycle_bits: list[int]) -> list[dict[str, Any]]:
    return [
        {
            "edge": edge,
            "audited_zero_closure_bit": audited,
            "odd_closure_perturbation_bit": odd,
            "differs_only_on_closing_edge": (audited ^ odd) == (1 if edge == CLOSING_EDGE else 0),
        }
        for edge, audited, odd in zip(CYCLE_EDGE_LABELS, audited_cycle_bits, odd_cycle_bits)
    ]


def build_payload() -> dict[str, Any]:
    z2 = load_json(Z2_REPORT)
    path_cohomology = load_json(PATH_COHOMOLOGY_REPORT)
    gf2 = load_json(GF2_REPORT)

    node_bits = [row["node_bit"] for row in z2["node_bit_rows"]]
    path_edge_bits = [row["edge_bit"] for row in z2["edge_bit_rows"]]
    forced_closing_bit = node_bits[-1] ^ node_bits[0]
    audited_cycle_bits = path_edge_bits + [forced_closing_bit]
    odd_cycle_bits = path_edge_bits + [forced_closing_bit ^ 1]
    matrix = cycle_coboundary_matrix()
    rref = gf2_rref(matrix)
    kernels = kernel_vectors(matrix)
    zero_closure = closure_case_row("audited_path_plus_forced_zero_closing_edge", audited_cycle_bits, node_bits, matrix)
    odd_closure = closure_case_row("odd_closing_edge_perturbation", odd_cycle_bits, node_bits, matrix)
    h1_dimension = CYCLE_EDGE_COUNT - rref["rank"]

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_CYCLE_CLOSURE_OBSTRUCTION_CERTIFICATE__GF2_BOUNDARY_CONDITION",
        "status": "cycle-closure-parity-functional-certified-zero-closure-exact-odd-closure-obstructed",
        "source_reports": {
            "phase_sign_z2_coboundary_certificate": str(Z2_REPORT.relative_to(ROOT)),
            "phase_sign_path_cohomology_triviality_certificate": str(PATH_COHOMOLOGY_REPORT.relative_to(ROOT)),
            "phase_sign_gf2_linear_system_certificate": str(GF2_REPORT.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "cycle closure parity",
                "closed 12-cycle",
                "closing edge 11->0",
                "cycle obstruction",
                "GF2 cycle functional",
                "boundary condition sensitivity",
            ],
            "conclusion": "Existing reports certify path exactness; this report adds the conditional closed-cycle parity functional and odd-closing-edge obstruction check.",
        },
        "cycle_complex_definition": {
            "field": "GF(2)",
            "node_count_dim_C0": NODE_COUNT,
            "cycle_edge_count_dim_C1": CYCLE_EDGE_COUNT,
            "connected_components": 1,
            "cycle_rank_E_minus_V_plus_components": CYCLE_EDGE_COUNT - NODE_COUNT + 1,
            "edge_labels": CYCLE_EDGE_LABELS,
            "cycle_coboundary_delta_matrix": matrix,
            "cycle_functional": "sum of all 12 cycle-edge bits mod 2",
            "exactness_criterion": "cycle edge cochain is exact iff total cycle parity is 0",
        },
        "cycle_coboundary_rank_certificate": rref,
        "kernel_certificate": {
            "kernel_vectors": kernels,
            "kernel_size": len(kernels),
            "kernel_is_exactly_constant_node_cochains": kernels == [[0] * NODE_COUNT, [1] * NODE_COUNT],
        },
        "closure_case_rows": [zero_closure, odd_closure],
        "cycle_edge_comparison_rows": edge_rows(audited_cycle_bits, odd_cycle_bits),
        "cycle_closure_summary": {
            "rank_delta_cycle": rref["rank"],
            "nullity_delta_cycle": rref["nullity"],
            "h1_dimension_dim_C1_minus_rank_delta": h1_dimension,
            "cycle_rank": CYCLE_EDGE_COUNT - NODE_COUNT + 1,
            "audited_path_edge_bits": path_edge_bits,
            "audited_node_bits": node_bits,
            "forced_closing_edge_bit_b11_xor_b0": forced_closing_bit,
            "zero_closing_cycle_parity_even": zero_closure["cycle_parity_even"],
            "zero_closing_cycle_exact": zero_closure["is_exact_cycle_coboundary"],
            "zero_closing_exact_potential_count": zero_closure["exact_node_potential_count"],
            "zero_closing_anchor_recovers_audited_nodes": zero_closure["anchored_reconstruction_matches_audited_nodes"],
            "odd_closing_cycle_parity_even": odd_closure["cycle_parity_even"],
            "odd_closing_cycle_exact": odd_closure["is_exact_cycle_coboundary"],
            "odd_closing_obstructed_by_cycle_parity": not odd_closure["is_exact_cycle_coboundary"] and not odd_closure["cycle_parity_even"],
            "path_h1_zero_inherited": path_cohomology["path_cohomology_summary"]["h1_dimension_dim_C1_minus_rank_delta"] == 0,
            "path_anchor_reconstruction_inherited": path_cohomology["path_cohomology_summary"]["anchored_reconstruction_matches_z2_node_bits"],
            "gf2_path_solution_inherited": gf2["gf2_linear_system_summary"]["unique_solution"],
        },
        "blocker_context": {
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_damping_derivation_from_nadsoliton_dynamics",
                "transport_cocycle_source_from_strict_dynamics",
                "QW-2191_selector_obstruction",
                "orientation_chi11_source",
            ],
            "guardrail": "This is a closed-cycle boundary-condition sanity check, not a cyclic gate-generation strategy and not a strict-source selector closure.",
        },
        "proof_certificate": {
            "cycle_rank_step": "The artificially closed graph has V=12, E=12, one connected component, rank(delta)=11, and one GF(2) cycle functional.",
            "criterion_step": "A 12-cycle edge cochain is exact only when the total cycle parity is 0.",
            "zero_closure_step": "The audited endpoint bits force closing edge b(11) xor b(0)=0; path bits plus this closing bit have even cycle parity and exactly the constant-pair node potentials.",
            "odd_closure_step": "Changing only the closing edge to 1 makes total cycle parity odd, so no node potential solves delta b=e on the closed cycle.",
            "theoretical_limit": "This boundary-condition certificate does not derive omega/phi, beta/eta, transport, a selector, or any bridge theorem from strict nadsoliton dynamics.",
        },
        "hard_limits": [
            "K_strict_gate remains the current live/full operational kernel.",
            "No unqualified identity K_legacy_ont == K_strict_gate is claimed.",
            "No cyclic L5/L12 gate-generation strategy is introduced or recommended.",
            "No proof derives A(d), P(d), D(d), omega/phi, beta/eta, or the transport cocycle from strict nadsoliton dynamics.",
            "No beta_tors -> chi_11 theorem is claimed.",
            "No legacy physical-role transfer to K_strict_gate is claimed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    summary = payload["cycle_closure_summary"]
    lines = [
        "# Phase-sign cycle-closure obstruction certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "The path sign cochain is tested against an artificial closing edge `11->0`.",
        "The forced zero closing bit is exact on the 12-cycle; the odd closing-bit",
        "perturbation is obstructed by the GF(2) cycle parity functional.",
        "",
        "## Summary",
        "",
        f"- Rank(delta_cycle): `{summary['rank_delta_cycle']}`",
        f"- Nullity(delta_cycle): `{summary['nullity_delta_cycle']}`",
        f"- H1 dimension: `{summary['h1_dimension_dim_C1_minus_rank_delta']}`",
        f"- Cycle rank: `{summary['cycle_rank']}`",
        f"- Forced closing edge bit b11 xor b0: `{summary['forced_closing_edge_bit_b11_xor_b0']}`",
        f"- Zero-closing cycle exact: `{summary['zero_closing_cycle_exact']}`",
        f"- Zero-closing anchor recovers audited nodes: `{summary['zero_closing_anchor_recovers_audited_nodes']}`",
        f"- Odd-closing cycle exact: `{summary['odd_closing_cycle_exact']}`",
        f"- Odd-closing obstructed by cycle parity: `{summary['odd_closing_obstructed_by_cycle_parity']}`",
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
    print(json.dumps(payload["cycle_closure_summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
