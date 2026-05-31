#!/usr/bin/env python3
"""Scratch probe: phase-sign support Euler-characteristic certificate.

The interval-boundary and flip-pair certificates identify the audited 1-node
support as two path components, [2,5] and [8,9].  This probe records the induced
support graph invariant that follows from that finite data:

    V_support - E_internal = component_count = 2.

It also cross-checks the path-boundary count.  Because neither support component
touches the path endpoints, each component has two boundary edges, so
boundary_weight = 2 * component_count = 4.  This is finite graph bookkeeping
only; it does not derive phase zeros, omega/phi, damping, or transport from
strict nadsoliton dynamics, and it does not claim a bridge theorem, selector
discharge, or ToE closure.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_sign_support_euler_characteristic_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_sign_support_euler_characteristic_certificate_report.md"
Z2_REPORT = HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json"
INTERVAL_BOUNDARY_REPORT = HERE / "bridge_strict_completion_phase_sign_node_support_interval_boundary_certificate_report.json"
FLIP_PAIR_REPORT = HERE / "bridge_strict_completion_phase_sign_flip_pair_interval_reconstruction_certificate_report.json"
EDGE_SUPPORT_REPORT = HERE / "bridge_strict_completion_phase_sign_edge_support_minimality_certificate_report.json"

NODE_COUNT = 12
EDGE_COUNT = 11
EXPECTED_SUPPORT_NODES = [2, 3, 4, 5, 8, 9]
EXPECTED_INTERNAL_EDGES = ["2->3", "3->4", "4->5", "8->9"]
EXPECTED_BOUNDARY_EDGES = ["1->2", "5->6", "7->8", "9->10"]
EXPECTED_COMPONENTS = [[2, 3, 4, 5], [8, 9]]
EDGE_LABELS = [f"{d}->{d + 1}" for d in range(EDGE_COUNT)]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def support_nodes(node_bits: list[int]) -> list[int]:
    return [node for node, bit in enumerate(node_bits) if bit]


def induced_internal_edges(node_support: set[int]) -> list[str]:
    return [EDGE_LABELS[index] for index in range(EDGE_COUNT) if index in node_support and index + 1 in node_support]


def boundary_edges(node_support: set[int]) -> list[str]:
    return [EDGE_LABELS[index] for index in range(EDGE_COUNT) if (index in node_support) ^ (index + 1 in node_support)]


def connected_components(nodes: list[int]) -> list[list[int]]:
    if not nodes:
        return []
    components = []
    current = [nodes[0]]
    for node in nodes[1:]:
        if node == current[-1] + 1:
            current.append(node)
        else:
            components.append(current)
            current = [node]
    components.append(current)
    return components


def component_rows(components: list[list[int]], node_support: set[int]) -> list[dict[str, Any]]:
    rows = []
    for index, component in enumerate(components):
        internal = [edge for edge in induced_internal_edges(set(component))]
        boundary = [edge for edge in boundary_edges(set(component))]
        touches_left = component[0] == 0
        touches_right = component[-1] == NODE_COUNT - 1
        rows.append({
            "component_index": index,
            "nodes": component,
            "node_count": len(component),
            "internal_edges": internal,
            "internal_edge_count": len(internal),
            "boundary_edges": boundary,
            "boundary_edge_count": len(boundary),
            "touches_left_endpoint": touches_left,
            "touches_right_endpoint": touches_right,
            "euler_v_minus_e": len(component) - len(internal),
            "is_path_tree_component": len(component) - len(internal) == 1,
            "boundary_count_matches_endpoint_formula": len(boundary) == 2 - int(touches_left) - int(touches_right),
            "component_is_subset_of_support": all(node in node_support for node in component),
        })
    return rows


def build_payload() -> dict[str, Any]:
    z2 = load_json(Z2_REPORT)
    interval_boundary = load_json(INTERVAL_BOUNDARY_REPORT)
    flip_pair = load_json(FLIP_PAIR_REPORT)
    edge_support = load_json(EDGE_SUPPORT_REPORT)

    node_bits = [row["node_bit"] for row in z2["node_bit_rows"]]
    edge_bits = [row["edge_bit"] for row in z2["edge_bit_rows"]]
    nodes = support_nodes(node_bits)
    node_set = set(nodes)
    internal = induced_internal_edges(node_set)
    boundary = boundary_edges(node_set)
    components = connected_components(nodes)
    rows = component_rows(components, node_set)
    v_count = len(nodes)
    e_internal_count = len(internal)
    component_count = len(components)
    boundary_weight = len(boundary)
    endpoint_touch_count = sum(int(row["touches_left_endpoint"]) + int(row["touches_right_endpoint"]) for row in rows)

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_SUPPORT_EULER_CHARACTERISTIC_CERTIFICATE__INDUCED_PATH_SUPPORT_GRAPH",
        "status": "support-induced-graph-euler-characteristic-and-boundary-count-certified",
        "source_reports": {
            "phase_sign_z2_coboundary_certificate": str(Z2_REPORT.relative_to(ROOT)),
            "phase_sign_node_support_interval_boundary_certificate": str(INTERVAL_BOUNDARY_REPORT.relative_to(ROOT)),
            "phase_sign_flip_pair_interval_reconstruction_certificate": str(FLIP_PAIR_REPORT.relative_to(ROOT)),
            "phase_sign_edge_support_minimality_certificate": str(EDGE_SUPPORT_REPORT.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "support Euler",
                "Euler characteristic",
                "induced support",
                "support graph",
                "internal edge",
                "V-E",
            ],
            "finding": "Existing strict-completion scratch reports identify intervals and boundaries, but do not export the induced support graph Euler characteristic and boundary-count certificate before this file.",
        },
        "support_graph_definition": {
            "graph": "induced subgraph of finite path 0--1--...--11 on nodes with node_bit=1",
            "support_nodes": nodes,
            "internal_edges": internal,
            "boundary_edges": boundary,
            "connected_components": components,
            "edge_bits_on_boundary_edges": {edge: edge_bits[index] for index, edge in enumerate(EDGE_LABELS) if edge in boundary},
        },
        "component_rows": rows,
        "euler_characteristic_certificate": {
            "support_vertex_count_V": v_count,
            "support_internal_edge_count_E": e_internal_count,
            "component_count_C": component_count,
            "euler_characteristic_V_minus_E": v_count - e_internal_count,
            "V_minus_E_equals_component_count": v_count - e_internal_count == component_count,
            "all_components_are_path_trees": all(row["is_path_tree_component"] for row in rows),
        },
        "boundary_count_certificate": {
            "boundary_edges": boundary,
            "boundary_weight": boundary_weight,
            "endpoint_touch_count": endpoint_touch_count,
            "predicted_boundary_weight_2C_minus_endpoint_touches": 2 * component_count - endpoint_touch_count,
            "boundary_weight_formula_matches": boundary_weight == 2 * component_count - endpoint_touch_count,
            "boundary_edges_equal_edge_bit_support": boundary == [edge for edge, bit in zip(EDGE_LABELS, edge_bits) if bit],
            "all_component_boundary_counts_match_endpoint_formula": all(row["boundary_count_matches_endpoint_formula"] for row in rows),
        },
        "support_euler_characteristic_summary": {
            "matches_expected_support_nodes": nodes == EXPECTED_SUPPORT_NODES,
            "matches_expected_internal_edges": internal == EXPECTED_INTERNAL_EDGES,
            "matches_expected_boundary_edges": boundary == EXPECTED_BOUNDARY_EDGES,
            "matches_expected_components": components == EXPECTED_COMPONENTS,
            "euler_characteristic_equals_component_count": v_count - e_internal_count == component_count,
            "component_count": component_count,
            "boundary_weight_formula_matches": boundary_weight == 2 * component_count - endpoint_touch_count,
            "boundary_edges_equal_edge_bit_support": boundary == [edge for edge, bit in zip(EDGE_LABELS, edge_bits) if bit],
            "all_components_are_path_trees": all(row["is_path_tree_component"] for row in rows),
            "matches_interval_boundary_components": components == [list(range(row["interval"]["start"], row["interval"]["end"] + 1)) for row in interval_boundary["interval_rows"]],
            "matches_flip_pair_intervals": components == [list(range(row["start"], row["end"] + 1)) for row in flip_pair["reconstruction_certificate"]["reconstructed_intervals_from_flip_pairs"]],
            "matches_edge_support_minimality": boundary == edge_support["edge_support_minimality_summary"]["derived_phase_sign_flip_edges"],
        },
        "blocker_context": {
            "resolved_locally": [
                "The induced 1-node support graph has V=6, E_internal=4, and C=2.",
                "The Euler characteristic V-E equals the component count C.",
                "The endpoint-free boundary count is 2*C=4 and equals the audited flip support.",
            ],
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_damping_transport_derivation_from_nadsoliton_dynamics",
                "QW-2191_selector_obstruction",
                "strict ToE closure",
            ],
        },
        "proof_certificate": {
            "support_step": "The audited 1-node support is {2,3,4,5,8,9} with internal path edges 2->3, 3->4, 4->5, and 8->9.",
            "euler_step": "The induced support graph has V=6, E=4, C=2, hence V-E=C=2.",
            "boundary_step": "Neither component touches a path endpoint, so boundary_weight=2*C=4 and the boundary edges are 1->2, 5->6, 7->8, and 9->10.",
            "chain_step": "The support graph components match interval-boundary rows, flip-pair reconstructed intervals, and edge-support minimality.",
            "theoretical_limit": "This is finite support-graph bookkeeping; it does not derive omega/phi, beta/eta, damping, or transport from strict nadsoliton dynamics.",
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
    summary = payload["support_euler_characteristic_summary"]
    euler = payload["euler_characteristic_certificate"]
    boundary = payload["boundary_count_certificate"]
    lines = [
        "# Phase-sign support Euler-characteristic certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This certificate computes the Euler characteristic of the induced",
        "1-node support graph and checks the corresponding path-boundary count.",
        "",
        "## Summary",
        "",
        f"- Support nodes: `{payload['support_graph_definition']['support_nodes']}`",
        f"- Internal edges: `{payload['support_graph_definition']['internal_edges']}`",
        f"- Components: `{payload['support_graph_definition']['connected_components']}`",
        f"- Euler V-E: `{euler['support_vertex_count_V']} - {euler['support_internal_edge_count_E']} = {euler['euler_characteristic_V_minus_E']}`",
        f"- Component count: `{euler['component_count_C']}`",
        f"- Boundary weight: `{boundary['boundary_weight']}`",
        f"- Boundary edges equal flip support: `{summary['boundary_edges_equal_edge_bit_support']}`",
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
    print(json.dumps(payload["support_euler_characteristic_summary"], indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
