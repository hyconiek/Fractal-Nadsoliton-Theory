#!/usr/bin/env python3
"""Scratch probe: self-recorded d5 endpoint anchor certificate.

The previous phase-origin certificates treated source recovery as requiring an
imported phase reference.  This probe checks a different, ontology-aligned route:
if the nadsoliton configuration is information and its finite d5 ledger already
exists, does that ledger contain an internal record of its own source/orientation?

Finite answer, conditional but constructive: yes for the already-assumed d5
support plus balanced endpoint-valued ledger.  The nonzero nodes form a unique
distance-5 path.  Its endpoints have distinct values (2 at one end, 1 at the
other), and the ordered values from the value-2 endpoint are exactly
(2,2,2,1,1).  Therefore the configuration itself encodes an equivariant endpoint
anchor: choose the value-2 endpoint as source and walk the d5 path toward the
value-1 endpoint.  This recovers every anchored representative and is D12-
equivariant.  It is not a D12-invariant absolute-origin selector, and it does not
derive the d5 support or balanced ledger from strict geometry.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PHASE_ALIASING = HERE / "bridge_strict_alpha_phase_mode_aliasing_obstruction_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_self_recorded_d5_endpoint_anchor_certificate_report.json"
OUT_MD = HERE / "bridge_strict_alpha_self_recorded_d5_endpoint_anchor_certificate_report.md"

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5
DISTANCE_SELECTED = 5
TARGET_BINARY_EXPONENT = 8
DENOMINATOR = 3
STRICT_TARGET_ETA = 9.0 / 5.0
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0
FORWARD_ASSIGNMENT = (2, 2, 2, 1, 1)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def eta_from_product(product: Fraction, branch_count: int) -> float:
    correction = float(product) ** (1.0 / branch_count)
    return math.log(NAD12_SUPPORT_SIZE * correction) / math.log(ALPHA_SCALE)


def cyclic_distance(left: int, right: int) -> int:
    raw = abs(left - right) % Z12_NODE_COUNT
    return min(raw, Z12_NODE_COUNT - raw)


def ordered_d5_support(source: int, orientation: int) -> tuple[int, ...]:
    if orientation not in {-1, 1}:
        raise ValueError("orientation must be +/-1")
    step = (orientation * DISTANCE_SELECTED) % Z12_NODE_COUNT
    return tuple((source + index * step) % Z12_NODE_COUNT for index in range(SUPPORT_SIZE))


def value_configuration(source: int, orientation: int, assignment: tuple[int, ...] = FORWARD_ASSIGNMENT) -> tuple[int, ...]:
    values = [0] * Z12_NODE_COUNT
    for node, value in zip(ordered_d5_support(source, orientation), assignment):
        values[node] = value
    return tuple(values)


def dihedral_action_on_config(config: tuple[int, ...], shift: int, reflect: bool) -> tuple[int, ...]:
    out = [0] * Z12_NODE_COUNT
    for node, value in enumerate(config):
        target = ((-node if reflect else node) + shift) % Z12_NODE_COUNT
        out[target] = value
    return tuple(out)


def transformed_anchor(source: int, orientation: int, shift: int, reflect: bool) -> tuple[int, int]:
    next_node = (source + orientation * DISTANCE_SELECTED) % Z12_NODE_COUNT
    transformed_source = ((-source if reflect else source) + shift) % Z12_NODE_COUNT
    transformed_next = ((-next_node if reflect else next_node) + shift) % Z12_NODE_COUNT
    if transformed_next == (transformed_source + DISTANCE_SELECTED) % Z12_NODE_COUNT:
        transformed_orientation = 1
    elif transformed_next == (transformed_source - DISTANCE_SELECTED) % Z12_NODE_COUNT:
        transformed_orientation = -1
    else:
        raise ValueError("transformed adjacent node is not distance-5 adjacent")
    return transformed_source, transformed_orientation


def support_nodes(config: tuple[int, ...]) -> list[int]:
    return [node for node, value in enumerate(config) if value != 0]


def d5_adjacency(config: tuple[int, ...]) -> dict[int, list[int]]:
    support = support_nodes(config)
    adjacency = {node: [] for node in support}
    for left_index, left in enumerate(support):
        for right in support[left_index + 1 :]:
            if cyclic_distance(left, right) == DISTANCE_SELECTED:
                adjacency[left].append(right)
                adjacency[right].append(left)
    return {node: sorted(neighbors) for node, neighbors in adjacency.items()}


def walk_path(adjacency: dict[int, list[int]], start: int) -> list[int]:
    path = [start]
    previous: int | None = None
    current = start
    while True:
        candidates = [neighbor for neighbor in adjacency[current] if neighbor != previous]
        if not candidates:
            return path
        if len(candidates) != 1:
            raise ValueError("d5 support graph is not a simple path")
        previous, current = current, candidates[0]
        path.append(current)


def infer_self_recorded_anchor(config: tuple[int, ...]) -> dict[str, Any]:
    adjacency = d5_adjacency(config)
    degrees = {node: len(neighbors) for node, neighbors in adjacency.items()}
    endpoints = sorted(node for node, degree in degrees.items() if degree == 1)
    if len(adjacency) != SUPPORT_SIZE or sorted(degrees.values()) != [1, 1, 2, 2, 2]:
        raise ValueError("support is not a five-node distance-5 path")
    value_2_endpoints = [node for node in endpoints if config[node] == 2]
    value_1_endpoints = [node for node in endpoints if config[node] == 1]
    if len(value_2_endpoints) != 1 or len(value_1_endpoints) != 1:
        raise ValueError("endpoint values do not uniquely mark a value-2 source and value-1 terminus")
    source = value_2_endpoints[0]
    ordered_path = walk_path(adjacency, source)
    if len(ordered_path) != SUPPORT_SIZE:
        raise ValueError("path walk did not cover the full support")
    second = ordered_path[1]
    if second == (source + DISTANCE_SELECTED) % Z12_NODE_COUNT:
        orientation = 1
    elif second == (source - DISTANCE_SELECTED) % Z12_NODE_COUNT:
        orientation = -1
    else:
        raise ValueError("path second node is not a signed d5 neighbor")
    ordered_values = [config[node] for node in ordered_path]
    return {
        "inferred_source": source,
        "inferred_orientation": orientation,
        "ordered_path": ordered_path,
        "ordered_values": ordered_values,
        "endpoints": endpoints,
        "endpoint_values": {str(node): config[node] for node in endpoints},
        "degrees": {str(node): degree for node, degree in degrees.items()},
    }


def all_anchor_rows() -> list[dict[str, Any]]:
    rows = []
    for source in range(Z12_NODE_COUNT):
        for orientation in (-1, 1):
            config = value_configuration(source, orientation)
            inferred = infer_self_recorded_anchor(config)
            rows.append(
                {
                    "actual_source": source,
                    "actual_orientation": orientation,
                    "value_configuration": list(config),
                    "inferred_anchor": inferred,
                    "source_matches": inferred["inferred_source"] == source,
                    "orientation_matches": inferred["inferred_orientation"] == orientation,
                    "ordered_values_match_forward_assignment": inferred["ordered_values"] == list(FORWARD_ASSIGNMENT),
                }
            )
    return rows


def equivariance_audit() -> dict[str, Any]:
    checked = 0
    mismatches = []
    for source in range(Z12_NODE_COUNT):
        for orientation in (-1, 1):
            config = value_configuration(source, orientation)
            for shift in range(Z12_NODE_COUNT):
                for reflect in (False, True):
                    transformed = dihedral_action_on_config(config, shift, reflect)
                    inferred = infer_self_recorded_anchor(transformed)
                    expected_source, expected_orientation = transformed_anchor(source, orientation, shift, reflect)
                    checked += 1
                    if inferred["inferred_source"] != expected_source or inferred["inferred_orientation"] != expected_orientation:
                        mismatches.append(
                            {
                                "source": source,
                                "orientation": orientation,
                                "shift": shift,
                                "reflect": reflect,
                                "expected_source": expected_source,
                                "expected_orientation": expected_orientation,
                                "inferred_source": inferred["inferred_source"],
                                "inferred_orientation": inferred["inferred_orientation"],
                            }
                        )
    return {
        "checked_cases": checked,
        "mismatch_count": len(mismatches),
        "mismatches": mismatches[:10],
        "all_cases_equivariant": len(mismatches) == 0,
    }


def main() -> None:
    phase_aliasing = load_json(PHASE_ALIASING)
    product = Fraction(2 ** TARGET_BINARY_EXPONENT, DENOMINATOR ** SUPPORT_SIZE)
    rows = all_anchor_rows()
    equivariance = equivariance_audit()

    report = {
        "status": "OPEN_STRICT_ALPHA_SELF_RECORDED_D5_ENDPOINT_ANCHOR_CERTIFICATE_PREMISE_BASED_NO_STRICT_DISCHARGE",
        "result_kind": "SCRATCH_STRICT_ALPHA_SELF_RECORDED_D5_ENDPOINT_ANCHOR_CERTIFICATE_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "phase_mode_aliasing_obstruction": str(PHASE_ALIASING.relative_to(ROOT)),
        },
        "previous_phase_aliasing_replay": {
            "result_kind": phase_aliasing["result_kind"],
            "source_complete_modes": phase_aliasing["phase_mode_aliasing_scan"]["source_complete_modes"],
            "all_observed_counts_match_gcd_formula": phase_aliasing["phase_mode_aliasing_scan"]["all_observed_counts_match_gcd_formula"],
            "candidate_status": phase_aliasing["candidate_interpretation"]["status"],
        },
        "target_identity_replay": {
            "q_power_product": f"{product.numerator}/{product.denominator}",
            "eta_residual_vs_9_5": eta_from_product(product, SUPPORT_SIZE) - STRICT_TARGET_ETA,
            "forward_assignment": list(FORWARD_ASSIGNMENT),
        },
        "self_recorded_anchor_rule": {
            "premises": [
                "The d5 support exists as a five-node distance-5 path.",
                "The endpoint-valued balanced ledger is already present as (2,2,2,1,1) along that path.",
            ],
            "rule": "Choose the unique value-2 endpoint as source and walk the d5 path toward the unique value-1 endpoint.",
            "ontology_read": "The finite configuration carries an internal endpoint record of its own orientation/source relative to the d5 path.",
            "not_an_absolute_origin": "The rule is D12-equivariant, not D12-invariant; it transports with the configuration and does not choose an external representative of the full orbit.",
        },
        "anchor_reconstruction_scan": {
            "row_count": len(rows),
            "rows": rows,
            "all_sources_recovered": all(row["source_matches"] for row in rows),
            "all_orientations_recovered": all(row["orientation_matches"] for row in rows),
            "all_ordered_values_match_forward_assignment": all(row["ordered_values_match_forward_assignment"] for row in rows),
        },
        "d12_equivariance_audit": equivariance,
        "selector_consequence": {
            "what_is_gained": "Given the d5 support and endpoint-valued ledger, the source/orientation are self-recorded by the finite value pattern without an external Fourier phase origin.",
            "why_this_does_not_contradict_D12_no_go": "The anchor is an equivariant extractor from a structured configuration, not a D12-invariant scalar score selecting an absolute orbit representative.",
            "what_remains_unproved": "The strict theory still has to derive the d5 support, balanced endpoint ledger, and endpoint-source convention from nadsoliton geometry.",
        },
        "candidate_interpretation": {
            "supported_by_this_probe": True,
            "content": "The balanced d5 ledger contains a self-recorded endpoint anchor once the support and ledger premises are granted.",
            "why_this_is_more_proof_like": "The probe gives a constructive graph rule, verifies all 24 anchored rows, and checks D12 equivariance over 576 transformed cases.",
            "why_this_is_not_enough": "The rule starts after the d5 support and endpoint-valued ledger are assumed; it does not derive those premises from strict nadsoliton ontology.",
            "status": "candidate-supported-but-d5-support-and-ledger-not-derived",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "The nadsoliton is treated as information carrying a finite self-record here, but no separate sub-nadsoliton informational layer is introduced.",
            "The endpoint anchor is conditional on the d5 support and endpoint-valued balanced ledger already existing.",
            "This is a D12-equivariant self-record extractor, not a D12-invariant absolute-origin selector theorem.",
            "No theorem derives the d5 support, balanced ledger, endpoint-source convention, or strict source-localizing term from strict nadsoliton geometry.",
            "No theorem derives branch count m=5, total exponent n=8, denominator 3, or binary-rescale quotient.",
            "No Aut(Z_12)/N462/QW-2191 selector obstruction is discharged by this self-recorded endpoint anchor certificate.",
            "No theorem derives eta=9/5 without the branch model, denominator/quotient convention, support premise, ledger selector, orientation premise, and endpoint-anchor premise.",
            "No QW-2191 discharge and no ToE closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Try to derive the endpoint-valued balanced d5 ledger as an internal self-record of nadsoliton formation/existence, or prove a no-go for such a derivation under the current strict gate assumptions.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha self-recorded d5 endpoint anchor certificate probe\n\n"
        "Status: premise-based endpoint self-record certificate; no strict selector discharge.\n\n"
        f"- Rows checked: `{len(rows)}`; sources recovered: `{report['anchor_reconstruction_scan']['all_sources_recovered']}`; orientations recovered: `{report['anchor_reconstruction_scan']['all_orientations_recovered']}`.\n"
        f"- Ordered values match `(2,2,2,1,1)`: `{report['anchor_reconstruction_scan']['all_ordered_values_match_forward_assignment']}`.\n"
        f"- D12 equivariance checked cases: `{equivariance['checked_cases']}`; mismatches: `{equivariance['mismatch_count']}`.\n"
        f"- Target replay: `q^5={product.numerator}/{product.denominator}`, eta residual `{eta_from_product(product, SUPPORT_SIZE)-STRICT_TARGET_ETA:.3e}`.\n"
        "- Honest read: the finite d5 ledger self-records an endpoint anchor once support and ledger premises are granted.\n"
        "- No false pass: no strict d5-support/ledger theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
