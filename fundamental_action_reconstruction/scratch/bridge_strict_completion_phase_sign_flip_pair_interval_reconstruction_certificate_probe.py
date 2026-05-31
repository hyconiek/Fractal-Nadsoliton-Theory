#!/usr/bin/env python3
"""Scratch probe: phase-sign flip-pair interval reconstruction certificate.

The node-support interval-boundary certificate maps maximal 1-node intervals to
edge flips.  This probe records the inverse finite path computation: starting
only from the ordered edge-bit support and the left anchor b(0)=0, scan the path,
pair entry/exit flips, and reconstruct the maximal 1-node intervals.

For the audited edge bits the flip indices are [1, 5, 7, 9], so the paired
entry/exit cuts reconstruct exactly [2,5] and [8,9].  This is a finite GF(2)
path-boundary inverse certificate; it is not a phase fit, not a strict dynamical
derivation, not a kernel bridge, and not a selector or ToE closure.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_sign_flip_pair_interval_reconstruction_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_sign_flip_pair_interval_reconstruction_certificate_report.md"
Z2_REPORT = HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json"
EDGE_SUPPORT_REPORT = HERE / "bridge_strict_completion_phase_sign_edge_support_minimality_certificate_report.json"
INTERVAL_BOUNDARY_REPORT = HERE / "bridge_strict_completion_phase_sign_node_support_interval_boundary_certificate_report.json"
REDUCED_COBOUNDARY_REPORT = HERE / "bridge_strict_completion_phase_sign_reduced_coboundary_inverse_certificate_report.json"

NODE_COUNT = 12
EDGE_COUNT = 11
ANCHOR_BIT = 0
EXPECTED_NODE_BITS = [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0]
EXPECTED_EDGE_BITS = [0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0]
EXPECTED_FLIP_EDGES = ["1->2", "5->6", "7->8", "9->10"]
EXPECTED_FLIP_INDICES = [1, 5, 7, 9]
EXPECTED_INTERVALS = [{"start": 2, "end": 5}, {"start": 8, "end": 9}]
EDGE_LABELS = [f"{d}->{d + 1}" for d in range(EDGE_COUNT)]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def support(bits: list[int], labels: list[str]) -> list[str]:
    return [label for label, bit in zip(labels, bits) if bit]


def anchored_node_reconstruction(edge_bits: list[int], anchor_bit: int = ANCHOR_BIT) -> list[int]:
    nodes = [anchor_bit]
    current = anchor_bit
    for bit in edge_bits:
        current ^= bit
        nodes.append(current)
    return nodes


def pair_flip_indices(edge_bits: list[int], anchor_bit: int = ANCHOR_BIT) -> dict[str, Any]:
    current = anchor_bit
    open_start = None
    flip_rows = []
    intervals = []
    for edge_index, bit in enumerate(edge_bits):
        before = current
        after = current ^ bit
        if bit:
            if before == 0 and after == 1:
                open_start = edge_index + 1
                role = "entry"
            elif before == 1 and after == 0:
                if open_start is None:
                    raise ValueError("exit flip encountered without an open interval")
                intervals.append({"start": open_start, "end": edge_index})
                open_start = None
                role = "exit"
            else:
                raise ValueError("unexpected GF(2) flip state")
            flip_rows.append({
                "flip_index": edge_index,
                "flip_edge": EDGE_LABELS[edge_index],
                "state_before": before,
                "state_after": after,
                "role": role,
                "open_interval_start_after_flip": open_start,
            })
        current = after
    return {
        "flip_indices": [row["flip_index"] for row in flip_rows],
        "flip_edges": [row["flip_edge"] for row in flip_rows],
        "flip_rows": flip_rows,
        "reconstructed_intervals": intervals,
        "final_state": current,
        "has_unclosed_interval": open_start is not None,
    }


def interval_node_vector(interval: dict[str, int]) -> list[int]:
    return [1 if interval["start"] <= node <= interval["end"] else 0 for node in range(NODE_COUNT)]


def interval_boundary_vector(interval: dict[str, int]) -> list[int]:
    vector = [0] * EDGE_COUNT
    if interval["start"] > 0:
        vector[interval["start"] - 1] = 1
    if interval["end"] < NODE_COUNT - 1:
        vector[interval["end"]] ^= 1
    return vector


def pair_rows(intervals: list[dict[str, int]], edge_bits: list[int]) -> list[dict[str, Any]]:
    rows = []
    for index, interval in enumerate(intervals):
        boundary = interval_boundary_vector(interval)
        rows.append({
            "pair_index": index,
            "entry_flip_edge": EDGE_LABELS[interval["start"] - 1] if interval["start"] > 0 else "left-boundary-anchor",
            "exit_flip_edge": EDGE_LABELS[interval["end"]] if interval["end"] < NODE_COUNT - 1 else "right-boundary-endpoint",
            "interval": interval,
            "interval_node_vector": interval_node_vector(interval),
            "interval_boundary_vector": boundary,
            "boundary_edges": support(boundary, EDGE_LABELS),
            "boundary_edges_match_edge_bits": all((b == 0 or edge_bits[i] == b) for i, b in enumerate(boundary)),
        })
    return rows


def boundary_sum(intervals: list[dict[str, int]]) -> list[int]:
    total = [0] * EDGE_COUNT
    for interval in intervals:
        for index, bit in enumerate(interval_boundary_vector(interval)):
            total[index] ^= bit
    return total


def node_union(intervals: list[dict[str, int]]) -> list[int]:
    total = [0] * NODE_COUNT
    for interval in intervals:
        for index, bit in enumerate(interval_node_vector(interval)):
            total[index] ^= bit
    return total


def build_payload() -> dict[str, Any]:
    z2 = load_json(Z2_REPORT)
    edge_support = load_json(EDGE_SUPPORT_REPORT)
    interval_boundary = load_json(INTERVAL_BOUNDARY_REPORT)
    reduced_coboundary = load_json(REDUCED_COBOUNDARY_REPORT)

    node_bits = [row["node_bit"] for row in z2["node_bit_rows"]]
    edge_bits = [row["edge_bit"] for row in z2["edge_bit_rows"]]
    pairing = pair_flip_indices(edge_bits)
    intervals = pairing["reconstructed_intervals"]
    reconstructed_nodes = anchored_node_reconstruction(edge_bits)
    interval_union_nodes = node_union(intervals)
    interval_boundary_edges = boundary_sum(intervals)
    rows = pair_rows(intervals, edge_bits)
    flip_count = len(pairing["flip_indices"])

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_FLIP_PAIR_INTERVAL_RECONSTRUCTION_CERTIFICATE__BOUNDARY_INVERSE_ON_PATH",
        "status": "ordered-flip-pairs-reconstruct-maximal-node-support-intervals",
        "source_reports": {
            "phase_sign_z2_coboundary_certificate": str(Z2_REPORT.relative_to(ROOT)),
            "phase_sign_edge_support_minimality_certificate": str(EDGE_SUPPORT_REPORT.relative_to(ROOT)),
            "phase_sign_node_support_interval_boundary_certificate": str(INTERVAL_BOUNDARY_REPORT.relative_to(ROOT)),
            "phase_sign_reduced_coboundary_inverse_certificate": str(REDUCED_COBOUNDARY_REPORT.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "flip pair",
                "paired flips",
                "edge pair",
                "boundary inverse",
                "support from boundary",
                "interval reconstruction",
                "endpoint pairing",
            ],
            "finding": "Existing reports export edge-support minimality and interval-to-boundary reconstruction, but not the ordered flip-pair inverse from boundary edges back to maximal support intervals.",
        },
        "flip_pair_reconstruction_definition": {
            "field": "GF(2)",
            "graph": "finite path 0--1--...--11",
            "anchor_bit_b0": ANCHOR_BIT,
            "edge_bits": edge_bits,
            "flip_indices": pairing["flip_indices"],
            "flip_edges": pairing["flip_edges"],
            "rule": "scan from b(0)=0; 0->1 flips open an interval at edge_index+1 and 1->0 flips close it at edge_index",
        },
        "flip_scan_rows": pairing["flip_rows"],
        "paired_interval_rows": rows,
        "reconstruction_certificate": {
            "ordered_flip_indices": pairing["flip_indices"],
            "ordered_flip_edges": pairing["flip_edges"],
            "reconstructed_intervals_from_flip_pairs": intervals,
            "node_bits_from_anchor_scan": reconstructed_nodes,
            "node_bits_from_reconstructed_interval_union": interval_union_nodes,
            "edge_bits_from_reconstructed_interval_boundaries": interval_boundary_edges,
            "final_scan_state": pairing["final_state"],
            "has_unclosed_interval": pairing["has_unclosed_interval"],
        },
        "parity_and_endpoint_certificate": {
            "flip_count": flip_count,
            "flip_count_even": flip_count % 2 == 0,
            "component_count_from_pairs": len(intervals),
            "flip_count_equals_two_components": flip_count == 2 * len(intervals),
            "left_endpoint_bit": reconstructed_nodes[0],
            "right_endpoint_bit": reconstructed_nodes[-1],
            "no_endpoint_support": reconstructed_nodes[0] == 0 and reconstructed_nodes[-1] == 0,
            "scan_returns_to_anchor_state": pairing["final_state"] == ANCHOR_BIT,
        },
        "flip_pair_interval_reconstruction_summary": {
            "matches_expected_flip_indices": pairing["flip_indices"] == EXPECTED_FLIP_INDICES,
            "matches_expected_flip_edges": pairing["flip_edges"] == EXPECTED_FLIP_EDGES,
            "matches_expected_intervals": intervals == EXPECTED_INTERVALS,
            "node_bits_from_scan_match_z2": reconstructed_nodes == node_bits,
            "node_bits_from_interval_union_match_z2": interval_union_nodes == node_bits,
            "edge_bits_from_interval_boundaries_match_z2": interval_boundary_edges == edge_bits,
            "all_pair_boundaries_match_edge_bits": all(row["boundary_edges_match_edge_bits"] for row in rows),
            "flip_count_even": flip_count % 2 == 0,
            "no_unclosed_interval": not pairing["has_unclosed_interval"],
            "no_endpoint_support": reconstructed_nodes[0] == 0 and reconstructed_nodes[-1] == 0,
            "component_count_from_pairs": len(intervals),
            "matches_interval_boundary_report": intervals == interval_boundary["node_support_interval_boundary_summary"]["maximal_intervals"],
            "matches_edge_support_minimality": pairing["flip_edges"] == edge_support["edge_support_minimality_summary"]["derived_phase_sign_flip_edges"],
            "matches_reduced_coboundary_nodes": reconstructed_nodes == reduced_coboundary["reconstruction_certificate"]["full_node_bits_from_anchor_and_reduced_inverse"],
        },
        "blocker_context": {
            "resolved_locally": [
                "The ordered flip support is paired into entry/exit cuts (1->2,5->6) and (7->8,9->10).",
                "Those pairs reconstruct exactly the maximal 1-node intervals [2,5] and [8,9].",
                "Even flip count and zero endpoint support close the path-boundary inverse sanity check.",
            ],
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_damping_transport_derivation_from_nadsoliton_dynamics",
                "QW-2191_selector_obstruction",
                "strict ToE closure",
            ],
        },
        "proof_certificate": {
            "scan_step": "Starting from b(0)=0, the ordered flip indices [1,5,7,9] toggle the scan state 0->1->0->1->0.",
            "pairing_step": "The entry/exit flip pairs (1->2,5->6) and (7->8,9->10) reconstruct intervals [2,5] and [8,9].",
            "boundary_inverse_step": "The union of reconstructed intervals recovers the audited node bits, and their GF(2) boundaries recover the audited edge bits.",
            "parity_step": "The flip count is even, equals 2*component_count, and the final scan state returns to the b(0)=0 anchor with no endpoint support.",
            "chain_step": "The reconstructed intervals match the node-support interval-boundary report, edge-support minimality, and reduced-coboundary node reconstruction.",
            "theoretical_limit": "This is finite path-boundary inverse bookkeeping; it does not derive omega/phi, beta/eta, damping, or transport from strict nadsoliton dynamics.",
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
    summary = payload["flip_pair_interval_reconstruction_summary"]
    parity = payload["parity_and_endpoint_certificate"]
    lines = [
        "# Phase-sign flip-pair interval reconstruction certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This certificate scans the ordered edge flips from the left anchor and",
        "pairs entry/exit flips to reconstruct the maximal node-support intervals.",
        "",
        "## Summary",
        "",
        f"- Ordered flip edges: `{payload['reconstruction_certificate']['ordered_flip_edges']}`",
        f"- Reconstructed intervals: `{payload['reconstruction_certificate']['reconstructed_intervals_from_flip_pairs']}`",
        f"- Node bits from scan match Z2: `{summary['node_bits_from_scan_match_z2']}`",
        f"- Edge bits from interval boundaries match Z2: `{summary['edge_bits_from_interval_boundaries_match_z2']}`",
        f"- Flip count even: `{parity['flip_count_even']}`",
        f"- Flip count equals two components: `{parity['flip_count_equals_two_components']}`",
        f"- No endpoint support: `{parity['no_endpoint_support']}`",
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
    print(json.dumps(payload["flip_pair_interval_reconstruction_summary"], indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
