#!/usr/bin/env python3
"""Scratch probe: phase-sign component-quotient adjacency certificate.

The support Euler certificate computes the induced graph invariant of the
negative/1-node support.  This probe records the complementary quotient picture:
collapse each maximal constant node-bit run of the finite path 0--1--...--11 to a
single component vertex.  For the audited phase-sign data the quotient is the
alternating path

    plus[0,1] -- minus[2,5] -- plus[6,7] -- minus[8,9] -- plus[10,11]

with four quotient edges, exactly the four flip edges.  Thus the sign-run
quotient is a tree with V=5, E=4, V-E=1, and every quotient edge crosses exactly
one audited flip edge.

This is finite graph bookkeeping only.  It does not derive phase zeros,
omega/phi, damping, or transport from strict nadsoliton dynamics, and it does not
claim a kernel bridge, selector discharge, or ToE closure.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_sign_component_quotient_adjacency_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_sign_component_quotient_adjacency_certificate_report.md"
Z2_REPORT = HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json"
SUPPORT_EULER_REPORT = HERE / "bridge_strict_completion_phase_sign_support_euler_characteristic_certificate_report.json"
FLIP_PAIR_REPORT = HERE / "bridge_strict_completion_phase_sign_flip_pair_interval_reconstruction_certificate_report.json"
EDGE_SUPPORT_REPORT = HERE / "bridge_strict_completion_phase_sign_edge_support_minimality_certificate_report.json"

NODE_COUNT = 12
EDGE_COUNT = 11
EDGE_LABELS = [f"{d}->{d + 1}" for d in range(EDGE_COUNT)]
EXPECTED_COMPONENTS = [
    {"component_index": 0, "bit": 0, "sign": 1, "start": 0, "end": 1, "nodes": [0, 1], "node_count": 2},
    {"component_index": 1, "bit": 1, "sign": -1, "start": 2, "end": 5, "nodes": [2, 3, 4, 5], "node_count": 4},
    {"component_index": 2, "bit": 0, "sign": 1, "start": 6, "end": 7, "nodes": [6, 7], "node_count": 2},
    {"component_index": 3, "bit": 1, "sign": -1, "start": 8, "end": 9, "nodes": [8, 9], "node_count": 2},
    {"component_index": 4, "bit": 0, "sign": 1, "start": 10, "end": 11, "nodes": [10, 11], "node_count": 2},
]
EXPECTED_QUOTIENT_EDGES = [
    {"left_component": 0, "right_component": 1, "path_edge": "1->2"},
    {"left_component": 1, "right_component": 2, "path_edge": "5->6"},
    {"left_component": 2, "right_component": 3, "path_edge": "7->8"},
    {"left_component": 3, "right_component": 4, "path_edge": "9->10"},
]
EXPECTED_FLIP_EDGES = ["1->2", "5->6", "7->8", "9->10"]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def sign_from_bit(bit: int) -> int:
    return -1 if bit else 1


def maximal_constant_components(node_bits: list[int]) -> list[dict[str, Any]]:
    components = []
    start = 0
    current = node_bits[0]
    for node in range(1, len(node_bits)):
        if node_bits[node] == current:
            continue
        nodes = list(range(start, node))
        components.append({
            "component_index": len(components),
            "bit": current,
            "sign": sign_from_bit(current),
            "start": start,
            "end": node - 1,
            "nodes": nodes,
            "node_count": len(nodes),
        })
        start = node
        current = node_bits[node]
    nodes = list(range(start, len(node_bits)))
    components.append({
        "component_index": len(components),
        "bit": current,
        "sign": sign_from_bit(current),
        "start": start,
        "end": len(node_bits) - 1,
        "nodes": nodes,
        "node_count": len(nodes),
    })
    return components


def quotient_edges(components: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for left, right in zip(components, components[1:]):
        edge_index = left["end"]
        rows.append({
            "quotient_edge_index": len(rows),
            "left_component": left["component_index"],
            "right_component": right["component_index"],
            "left_bit": left["bit"],
            "right_bit": right["bit"],
            "left_sign": left["sign"],
            "right_sign": right["sign"],
            "path_edge": EDGE_LABELS[edge_index],
            "path_edge_index": edge_index,
            "is_bit_flip": left["bit"] ^ right["bit"] == 1,
            "is_sign_flip": left["sign"] * right["sign"] == -1,
        })
    return rows


def component_adjacency_matrix(component_count: int, edges: list[dict[str, Any]]) -> list[list[int]]:
    matrix = [[0] * component_count for _ in range(component_count)]
    for edge in edges:
        left = edge["left_component"]
        right = edge["right_component"]
        matrix[left][right] = 1
        matrix[right][left] = 1
    return matrix


def degree_rows(component_count: int, adjacency: list[list[int]]) -> list[dict[str, Any]]:
    return [
        {
            "component_index": index,
            "quotient_degree": sum(adjacency[index]),
            "is_endpoint_component": index in {0, component_count - 1},
        }
        for index in range(component_count)
    ]


def build_payload() -> dict[str, Any]:
    z2 = load_json(Z2_REPORT)
    support_euler = load_json(SUPPORT_EULER_REPORT)
    flip_pair = load_json(FLIP_PAIR_REPORT)
    edge_support = load_json(EDGE_SUPPORT_REPORT)

    node_bits = [row["node_bit"] for row in z2["node_bit_rows"]]
    edge_bits = [row["edge_bit"] for row in z2["edge_bit_rows"]]
    components = maximal_constant_components(node_bits)
    edges = quotient_edges(components)
    adjacency = component_adjacency_matrix(len(components), edges)
    degrees = degree_rows(len(components), adjacency)
    flip_edges = [edge for edge, bit in zip(EDGE_LABELS, edge_bits) if bit]
    negative_components = [component for component in components if component["bit"] == 1]
    positive_components = [component for component in components if component["bit"] == 0]

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_COMPONENT_QUOTIENT_ADJACENCY_CERTIFICATE__SIGN_RUN_TREE",
        "status": "constant-sign-components-collapse-to-alternating-quotient-path-tree",
        "source_reports": {
            "phase_sign_z2_coboundary_certificate": str(Z2_REPORT.relative_to(ROOT)),
            "phase_sign_support_euler_characteristic_certificate": str(SUPPORT_EULER_REPORT.relative_to(ROOT)),
            "phase_sign_flip_pair_interval_reconstruction_certificate": str(FLIP_PAIR_REPORT.relative_to(ROOT)),
            "phase_sign_edge_support_minimality_certificate": str(EDGE_SUPPORT_REPORT.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "quotient component",
                "component adjacency",
                "alternating sign",
                "sign-run",
                "sign run",
                "component quotient",
                "run quotient",
            ],
            "finding": "Existing reports export support components and flip-pair reconstruction, but not the collapsed constant-sign component adjacency quotient tree before this file.",
        },
        "component_quotient_definition": {
            "graph": "quotient of the finite path by maximal constant node-bit runs",
            "node_bits": node_bits,
            "component_rows": components,
            "quotient_edge_rows": edges,
            "component_adjacency_matrix": adjacency,
            "degree_rows": degrees,
        },
        "quotient_tree_certificate": {
            "component_vertex_count_Vq": len(components),
            "quotient_edge_count_Eq": len(edges),
            "euler_Vq_minus_Eq": len(components) - len(edges),
            "is_tree_by_path_quotient_euler": len(components) - len(edges) == 1,
            "all_quotient_edges_are_bit_flips": all(edge["is_bit_flip"] for edge in edges),
            "all_quotient_edges_are_sign_flips": all(edge["is_sign_flip"] for edge in edges),
            "endpoint_components_have_degree_one": all(row["quotient_degree"] == 1 for row in degrees if row["is_endpoint_component"]),
            "internal_components_have_degree_two": all(row["quotient_degree"] == 2 for row in degrees if not row["is_endpoint_component"]),
        },
        "component_count_certificate": {
            "positive_component_count": len(positive_components),
            "negative_component_count": len(negative_components),
            "flip_edge_count": len(flip_edges),
            "total_component_count": len(components),
            "component_count_equals_flip_count_plus_one": len(components) == len(flip_edges) + 1,
            "negative_components_match_support_euler_component_count": len(negative_components) == support_euler["euler_characteristic_certificate"]["component_count_C"],
            "negative_component_intervals": [{"start": component["start"], "end": component["end"]} for component in negative_components],
            "positive_component_intervals": [{"start": component["start"], "end": component["end"]} for component in positive_components],
        },
        "component_quotient_adjacency_summary": {
            "matches_expected_components": components == EXPECTED_COMPONENTS,
            "matches_expected_quotient_edges": [
                {"left_component": edge["left_component"], "right_component": edge["right_component"], "path_edge": edge["path_edge"]}
                for edge in edges
            ] == EXPECTED_QUOTIENT_EDGES,
            "quotient_is_tree": len(components) - len(edges) == 1,
            "quotient_is_alternating": all(edge["is_bit_flip"] and edge["is_sign_flip"] for edge in edges),
            "quotient_edges_match_flip_edges": [edge["path_edge"] for edge in edges] == flip_edges == EXPECTED_FLIP_EDGES,
            "component_count_equals_flip_count_plus_one": len(components) == len(flip_edges) + 1,
            "negative_components_match_support_euler": len(negative_components) == support_euler["euler_characteristic_certificate"]["component_count_C"],
            "negative_intervals_match_flip_pair": [{"start": component["start"], "end": component["end"]} for component in negative_components] == flip_pair["reconstruction_certificate"]["reconstructed_intervals_from_flip_pairs"],
            "matches_edge_support_minimality": [edge["path_edge"] for edge in edges] == edge_support["edge_support_minimality_summary"]["derived_phase_sign_flip_edges"],
        },
        "blocker_context": {
            "resolved_locally": [
                "The maximal constant-sign runs collapse to a five-vertex alternating quotient path.",
                "The quotient has V=5, E=4, and V-E=1, so it is a tree in this finite path quotient.",
                "The four quotient edges are exactly the audited four flip edges.",
            ],
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_damping_transport_derivation_from_nadsoliton_dynamics",
                "QW-2191_selector_obstruction",
                "strict ToE closure",
            ],
        },
        "proof_certificate": {
            "component_step": "Maximal constant node-bit runs are +[0,1], -[2,5], +[6,7], -[8,9], +[10,11].",
            "quotient_step": "Collapsing those runs gives an alternating quotient path with V=5 and E=4, hence V-E=1.",
            "flip_step": "The quotient edges are exactly 1->2, 5->6, 7->8, and 9->10, matching the audited flip support.",
            "count_step": "Total component count equals flip_count+1, while the two negative components match the support Euler and flip-pair interval certificates.",
            "theoretical_limit": "This is finite component-quotient bookkeeping; it does not derive omega/phi, beta/eta, damping, or transport from strict nadsoliton dynamics.",
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
    summary = payload["component_quotient_adjacency_summary"]
    tree = payload["quotient_tree_certificate"]
    counts = payload["component_count_certificate"]
    lines = [
        "# Phase-sign component-quotient adjacency certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This certificate collapses maximal constant-sign runs to a quotient path",
        "and checks that the quotient adjacency edges are exactly the audited flips.",
        "",
        "## Summary",
        "",
        f"- Quotient vertices Vq: `{tree['component_vertex_count_Vq']}`",
        f"- Quotient edges Eq: `{tree['quotient_edge_count_Eq']}`",
        f"- Euler Vq-Eq: `{tree['euler_Vq_minus_Eq']}`",
        f"- Positive components: `{counts['positive_component_intervals']}`",
        f"- Negative components: `{counts['negative_component_intervals']}`",
        f"- Quotient edges match flip edges: `{summary['quotient_edges_match_flip_edges']}`",
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
    print(json.dumps(payload["component_quotient_adjacency_summary"], indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
