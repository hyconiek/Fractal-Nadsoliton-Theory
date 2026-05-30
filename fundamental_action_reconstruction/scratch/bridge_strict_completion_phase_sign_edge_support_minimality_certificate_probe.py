#!/usr/bin/env python3
"""Scratch probe: phase-sign edge-support uniqueness/minimality certificate.

The phase-sign Z2 coboundary certificate already fixes node bits b(d) and edge
bits e(d,d+1)=b(d) xor b(d+1) on the finite path 0--1--...--11.  This probe
adds a deliberately finite proof check: among all 2^11 possible edge-bit
assignments on that path, exactly one assignment reconstructs the audited node
bit pattern from the left anchor b(0)=0.  Its support has size four, so no
smaller edge-flip support can reproduce the same sign pattern.

This is a uniqueness/minimality certificate for the finite path graph only.  It
is not a new phase fit, not a strict dynamical derivation, and not a bridge
between legacy and strict kernels.
"""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_sign_edge_support_minimality_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_sign_edge_support_minimality_certificate_report.md"
Z2_REPORT = HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json"
CELL_SIGN_REPORT = HERE / "bridge_strict_completion_phase_zero_cell_sign_certificate_report.json"
POSITIVE_FACTOR_REPORT = HERE / "bridge_strict_completion_positive_factor_sign_separation_certificate_report.json"

EXPECTED_NODE_BITS = [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0]
EXPECTED_EDGE_BITS = [0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0]
EXPECTED_SIGN_PATTERN = [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1]
EXPECTED_FLIP_EDGES = ["1->2", "5->6", "7->8", "9->10"]
EDGE_LABELS = [f"{d}->{d + 1}" for d in range(11)]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def reconstruct_node_bits(anchor_bit: int, edge_bits: list[int]) -> list[int]:
    bits = [anchor_bit]
    current = anchor_bit
    for bit in edge_bits:
        current ^= bit
        bits.append(current)
    return bits


def support(edge_bits: list[int]) -> list[str]:
    return [edge for edge, bit in zip(EDGE_LABELS, edge_bits) if bit == 1]


def first_mismatch(reconstructed: list[int], expected: list[int]) -> dict[str, Any] | None:
    for d, (got, want) in enumerate(zip(reconstructed, expected)):
        if got != want:
            return {"d": d, "reconstructed_bit": got, "expected_bit": want}
    return None


def enumerate_assignments(anchor_bit: int, expected_node_bits: list[int]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, int]]:
    witness_rows = []
    matching_rows = []
    support_histogram = {str(size): 0 for size in range(len(EDGE_LABELS) + 1)}
    for mask in range(1 << len(EDGE_LABELS)):
        edge_bits = [(mask >> i) & 1 for i in range(len(EDGE_LABELS))]
        reconstructed = reconstruct_node_bits(anchor_bit, edge_bits)
        support_size = sum(edge_bits)
        support_histogram[str(support_size)] += 1
        matches = reconstructed == expected_node_bits
        if matches:
            matching_rows.append({
                "mask": mask,
                "edge_bits": edge_bits,
                "support_edges": support(edge_bits),
                "support_size": support_size,
                "reconstructed_node_bits": reconstructed,
            })
        elif support_size < sum(EXPECTED_EDGE_BITS):
            mismatch = first_mismatch(reconstructed, expected_node_bits)
            witness_rows.append({
                "mask": mask,
                "support_size": support_size,
                "support_edges": support(edge_bits),
                "first_mismatch": mismatch,
            })
    return matching_rows, witness_rows, support_histogram


def boundary_forced_rows(node_bits: list[int], edge_bits: list[int]) -> list[dict[str, Any]]:
    rows = []
    for d, edge_bit in enumerate(edge_bits):
        left = node_bits[d]
        right = node_bits[d + 1]
        forced = left ^ right
        rows.append({
            "edge": EDGE_LABELS[d],
            "left_node_bit": left,
            "right_node_bit": right,
            "forced_edge_bit": forced,
            "z2_edge_bit": edge_bit,
            "is_forced_flip": forced == 1,
            "matches_z2_edge_bit": forced == edge_bit,
        })
    return rows


def lower_support_count_rows(anchor_bit: int, expected_node_bits: list[int], target_weight: int) -> list[dict[str, Any]]:
    rows = []
    for size in range(target_weight):
        total = 0
        failures = 0
        sample_failure = None
        for edges in combinations(range(len(EDGE_LABELS)), size):
            total += 1
            edge_bits = [0] * len(EDGE_LABELS)
            for edge_index in edges:
                edge_bits[edge_index] = 1
            reconstructed = reconstruct_node_bits(anchor_bit, edge_bits)
            if reconstructed != expected_node_bits:
                failures += 1
                if sample_failure is None:
                    sample_failure = {
                        "support_edges": support(edge_bits),
                        "first_mismatch": first_mismatch(reconstructed, expected_node_bits),
                    }
        rows.append({
            "support_size": size,
            "assignments_checked": total,
            "assignments_failing_to_reconstruct_expected_pattern": failures,
            "all_assignments_fail": failures == total,
            "sample_failure": sample_failure,
        })
    return rows


def build_payload() -> dict[str, Any]:
    z2 = load_json(Z2_REPORT)
    cell_sign = load_json(CELL_SIGN_REPORT)
    positive_factor = load_json(POSITIVE_FACTOR_REPORT)

    node_bits = [row["node_bit"] for row in z2["node_bit_rows"]]
    edge_bits = [row["edge_bit"] for row in z2["edge_bit_rows"]]
    anchor_bit = node_bits[0]
    matching_rows, lower_witness_rows, support_histogram = enumerate_assignments(anchor_bit, node_bits)
    target_support_size = sum(edge_bits)
    forced_rows = boundary_forced_rows(node_bits, edge_bits)
    lower_rows = lower_support_count_rows(anchor_bit, node_bits, target_support_size)
    unique_match = len(matching_rows) == 1

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_EDGE_SUPPORT_MINIMALITY_CERTIFICATE__FINITE_PATH_EXHAUSTIVE",
        "status": "unique-four-edge-z2-support-is-minimal-for-the-audited-phase-sign-pattern",
        "source_reports": {
            "phase_sign_z2_coboundary_certificate": str(Z2_REPORT.relative_to(ROOT)),
            "cell_sign_certificate": str(CELL_SIGN_REPORT.relative_to(ROOT)),
            "positive_factor_sign_separation_certificate": str(POSITIVE_FACTOR_REPORT.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "edge-support uniqueness",
                "support minimal",
                "phase_sign_edge",
                "edge_support_minimal",
                "lower_support",
                "exhaustive",
            ],
            "conclusion": "No existing strict-completion probe exported an exhaustive finite path edge-support uniqueness/minimality certificate before this file.",
        },
        "finite_path_definitions": {
            "domain_nodes": list(range(12)),
            "edge_labels": EDGE_LABELS,
            "anchor_bit": anchor_bit,
            "reconstruction_rule": "b_hat(0)=anchor_bit and b_hat(d+1)=b_hat(d) xor e(d,d+1)",
            "minimality_claim": "Every edge-bit assignment reconstructing the audited node-bit pattern equals the Z2 coboundary; therefore its support size is forced and minimal.",
        },
        "boundary_forced_edge_rows": forced_rows,
        "matching_edge_assignment_rows": matching_rows,
        "lower_support_exhaustion_rows": lower_rows,
        "lower_support_failure_witness_rows": lower_witness_rows,
        "support_histogram_all_2_pow_11_assignments": support_histogram,
        "edge_support_minimality_summary": {
            "node_bit_pattern": node_bits,
            "edge_bit_pattern": edge_bits,
            "phase_sign_pattern": z2["z2_coboundary_summary"]["phase_sign_pattern"],
            "derived_phase_sign_flip_edges": support(edge_bits),
            "target_support_size": target_support_size,
            "total_edge_assignments_checked": 1 << len(EDGE_LABELS),
            "matching_assignment_count": len(matching_rows),
            "unique_matching_assignment": unique_match,
            "unique_matching_assignment_support_edges": matching_rows[0]["support_edges"] if unique_match else [],
            "all_boundary_forced_rows_match_z2": all(row["matches_z2_edge_bit"] for row in forced_rows),
            "all_lower_support_assignments_fail": all(row["all_assignments_fail"] for row in lower_rows),
            "lower_support_assignments_checked": sum(row["assignments_checked"] for row in lower_rows),
            "matches_expected_node_bits": node_bits == EXPECTED_NODE_BITS,
            "matches_expected_edge_bits": edge_bits == EXPECTED_EDGE_BITS,
            "matches_expected_sign_pattern": z2["z2_coboundary_summary"]["phase_sign_pattern"] == EXPECTED_SIGN_PATTERN,
            "matches_expected_flip_edges": support(edge_bits) == EXPECTED_FLIP_EDGES,
            "matches_cell_sign_flip_edges": support(edge_bits) == cell_sign["cell_sign_summary"]["derived_phase_sign_flip_edges"],
            "matches_positive_factor_completion_flips": support(edge_bits) == positive_factor["positive_factor_sign_summary"]["derived_completion_flip_edges"],
        },
        "blocker_context": {
            "what_this_certifies": "finite-path uniqueness and support minimality of the already exported Z2 phase-flip edge assignment",
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_transport_derivation_from_nadsoliton_dynamics",
                "strict_damping_parameter_derivation_from_nadsoliton_dynamics",
                "orientation_chi11_source",
                "QW-2191_selector_obstruction",
            ],
        },
        "proof_certificate": {
            "boundary_step": "On a path graph, each edge bit is forced by adjacent node bits: e_i=b_i xor b_{i+1}.",
            "exhaustive_step": "All 2^11 edge-bit assignments were enumerated against the left-anchor reconstruction rule.",
            "uniqueness_step": "Exactly one assignment reconstructs the audited node-bit pattern, and it is the Z2 coboundary assignment.",
            "minimality_step": "All assignments with support size 0, 1, 2, or 3 fail, so the four-edge flip support is minimal.",
            "theoretical_limit": "This proves finite sign-support uniqueness only; it does not derive omega/phi, beta/eta, or transport from strict nadsoliton dynamics.",
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
    summary = payload["edge_support_minimality_summary"]
    lines = [
        "# Phase-sign edge-support minimality certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This certificate treats the phase-sign Z2 data as a finite path-graph",
        "problem.  It exhaustively enumerates all `2^11` edge-bit assignments and",
        "checks which assignments reconstruct the audited node-bit pattern from the",
        "left anchor.",
        "",
        "## Summary",
        "",
        f"- Total edge assignments checked: `{summary['total_edge_assignments_checked']}`",
        f"- Matching assignment count: `{summary['matching_assignment_count']}`",
        f"- Unique matching assignment: `{summary['unique_matching_assignment']}`",
        f"- Target support size: `{summary['target_support_size']}`",
        f"- Flip edges: `{summary['derived_phase_sign_flip_edges']}`",
        f"- Lower-support assignments checked: `{summary['lower_support_assignments_checked']}`",
        f"- All lower-support assignments fail: `{summary['all_lower_support_assignments_fail']}`",
        "",
        "## Lower-support exhaustion",
        "",
        "| support size | checked | all fail |",
        "| ---: | ---: | :---: |",
    ]
    for row in payload["lower_support_exhaustion_rows"]:
        lines.append(f"| {row['support_size']} | {row['assignments_checked']} | `{row['all_assignments_fail']}` |")
    lines.extend([
        "",
        "## Hard limits",
        "",
    ])
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
