#!/usr/bin/env python3
"""Scratch probe: phase-zero cell-sign certificate without trig evaluation.

The cell-partition certificate proved that the in-domain zero-carriers are
ordered, disjoint, and cut [0,11] into positive rational cells.  This probe
uses that ordered carrier list as a finite sign theorem: starting from the
left anchor sign +1, each integer node inherits sign (-1)^N where N is the
number of phase-zero carriers strictly to its left.

The result is a rational/combinatorial derivation of the audited phase sign
pattern from the zero-carrier cell partition, not a new float cosine evaluation
and not a strict dynamical derivation of omega/phi or transport.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_zero_cell_sign_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_zero_cell_sign_certificate_report.md"
CELL_PARTITION_REPORT = HERE / "bridge_strict_completion_phase_zero_cell_partition_certificate_report.json"
NODE_CLEARANCE_REPORT = HERE / "bridge_strict_completion_phase_zero_node_clearance_certificate_report.json"
RATIONAL_ZERO_REPORT = HERE / "bridge_strict_completion_phase_zero_rational_interval_certificate_report.json"

DOMAIN = list(range(12))
EXPECTED_SIGN_PATTERN = [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1]
EXPECTED_FLIP_EDGES = ["1->2", "5->6", "7->8", "9->10"]
ANCHOR_SIGN = 1


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def fraction_from_payload(payload: dict[str, Any]) -> Fraction:
    return Fraction(payload["numerator"], payload["denominator"])


def fraction_payload(value: Fraction) -> dict[str, Any]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
        "text": f"{value.numerator}/{value.denominator}",
    }


def carrier_interval(carrier: dict[str, Any]) -> tuple[Fraction, Fraction]:
    return fraction_from_payload(carrier["lower"]), fraction_from_payload(carrier["upper"])


def carrier_relation_to_node(carrier: dict[str, Any], d: int) -> str:
    node = Fraction(d, 1)
    lower, upper = carrier_interval(carrier)
    if upper < node:
        return "strictly_left_of_node"
    if node < lower:
        return "strictly_right_of_node"
    return "node_inside_zero_carrier_interval"


def node_sign_rows(carriers: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for d in DOMAIN:
        relations = [carrier_relation_to_node(carrier, d) for carrier in carriers]
        left_labels = [carrier["label"] for carrier, relation in zip(carriers, relations) if relation == "strictly_left_of_node"]
        blocking_labels = [carrier["label"] for carrier, relation in zip(carriers, relations) if relation == "node_inside_zero_carrier_interval"]
        zero_count_left = len(left_labels)
        sign = ANCHOR_SIGN * (-1 if zero_count_left % 2 else 1)
        rows.append({
            "d": d,
            "zero_carriers_left_of_node": left_labels,
            "zero_carrier_count_left": zero_count_left,
            "left_count_parity": "odd" if zero_count_left % 2 else "even",
            "node_inside_any_zero_carrier": bool(blocking_labels),
            "blocking_zero_carriers": blocking_labels,
            "derived_phase_transport_sign": sign,
            "expected_phase_transport_sign": EXPECTED_SIGN_PATTERN[d],
            "matches_expected_sign": sign == EXPECTED_SIGN_PATTERN[d],
        })
    return rows


def edge_sign_rows(node_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for left, right in zip(node_rows, node_rows[1:]):
        sign_flip = left["derived_phase_transport_sign"] != right["derived_phase_transport_sign"]
        new_carriers = [
            label for label in right["zero_carriers_left_of_node"]
            if label not in left["zero_carriers_left_of_node"]
        ]
        edge = f"{left['d']}->{right['d']}"
        rows.append({
            "edge": edge,
            "left_sign": left["derived_phase_transport_sign"],
            "right_sign": right["derived_phase_transport_sign"],
            "sign_flip": sign_flip,
            "new_zero_carriers_crossed": new_carriers,
            "crossed_zero_carrier_count": len(new_carriers),
            "odd_crossing_parity": len(new_carriers) % 2 == 1,
            "matches_crossing_parity": sign_flip == (len(new_carriers) % 2 == 1),
        })
    return rows


def build_payload() -> dict[str, Any]:
    cell_partition = load_json(CELL_PARTITION_REPORT)
    node_clearance = load_json(NODE_CLEARANCE_REPORT)
    rational_zero = load_json(RATIONAL_ZERO_REPORT)
    carriers = cell_partition["domain_zero_carriers_ordered"]
    node_rows = node_sign_rows(carriers)
    edge_rows = edge_sign_rows(node_rows)
    derived_pattern = [row["derived_phase_transport_sign"] for row in node_rows]
    derived_flip_edges = [row["edge"] for row in edge_rows if row["sign_flip"]]
    max_crossed_per_edge = max(row["crossed_zero_carrier_count"] for row in edge_rows)
    min_node_clearance = fraction_from_payload(node_clearance["clearance_summary"]["min_combined_phase_zero_node_clearance_lower_bound"])
    min_cell_length = fraction_from_payload(cell_partition["cell_partition_summary"]["min_cell_length"])

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_ZERO_CELL_SIGN_CERTIFICATE__COMBINATORIAL_NO_TRIG_EVAL",
        "status": "phase-sign-pattern-derived-from-rational-cell-partition-and-left-anchor",
        "source_reports": {
            "cell_partition_certificate": str(CELL_PARTITION_REPORT.relative_to(ROOT)),
            "node_clearance_certificate": str(NODE_CLEARANCE_REPORT.relative_to(ROOT)),
            "rational_zero_certificate": str(RATIONAL_ZERO_REPORT.relative_to(ROOT)),
        },
        "sign_rule": {
            "left_anchor_sign_at_d0": ANCHOR_SIGN,
            "node_sign_formula": "sign(d)=left_anchor_sign*(-1)**number_of_ordered_zero_carriers_strictly_left_of_d",
            "zero_carrier_order": [carrier["label"] for carrier in carriers],
            "uses_trig_evaluation": False,
        },
        "node_sign_rows": node_rows,
        "edge_sign_rows": edge_rows,
        "cell_sign_summary": {
            "all_nodes_outside_zero_carriers": all(not row["node_inside_any_zero_carrier"] for row in node_rows),
            "all_node_signs_match_expected": all(row["matches_expected_sign"] for row in node_rows),
            "all_edge_flips_match_crossing_parity": all(row["matches_crossing_parity"] for row in edge_rows),
            "max_crossed_zero_carriers_per_integer_edge": max_crossed_per_edge,
            "derived_phase_transport_sign_pattern": derived_pattern,
            "derived_phase_sign_flip_edges": derived_flip_edges,
            "matches_cell_partition_sign_pattern": derived_pattern == cell_partition["cell_partition_summary"]["derived_phase_transport_sign_pattern"],
            "matches_cell_partition_flip_edges": derived_flip_edges == cell_partition["cell_partition_summary"]["derived_phase_sign_flip_edges"],
            "matches_node_clearance_sign_pattern": derived_pattern == node_clearance["clearance_summary"]["certified_phase_sign_pattern_preserved"],
            "matches_rational_zero_sign_pattern": derived_pattern == rational_zero["interval_summary"]["phase_transport_sign_pattern_from_rational_intervals"],
            "matches_expected_sign_pattern": derived_pattern == EXPECTED_SIGN_PATTERN,
            "matches_expected_flip_edges": derived_flip_edges == EXPECTED_FLIP_EDGES,
            "min_node_clearance_inherited": fraction_payload(min_node_clearance),
            "min_cell_length_inherited": fraction_payload(min_cell_length),
        },
        "blocker_context": {
            "what_this_refines": "turns the ordered rational phase-zero cell partition into a node sign theorem without trigonometric evaluation",
            "cell_partition_status": cell_partition["status"],
            "node_clearance_status": node_clearance["status"],
            "rational_zero_status": rational_zero["status"],
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_transport_derivation_from_nadsoliton_dynamics",
                "orientation_chi11_source",
                "chi11_uniqueness",
                "role_transfer_theorem",
            ],
        },
        "proof_certificate": {
            "anchor_step": "The left anchor sign at d=0 is fixed as +1 from the already-certified phase pattern context.",
            "counting_step": "For each integer node, count ordered zero-carriers strictly to the left using rational interval comparisons only.",
            "node_step": "Node-clearance guarantees no integer node lies inside a zero-carrier interval, so each node has a well-defined cell sign.",
            "edge_step": "Adjacent node signs differ exactly when an odd number of zero-carriers is crossed by the edge.",
            "nonduplication": "This is a combinatorial cell-sign theorem, not another zero-location, node-clearance, cell-partition, damping, cocycle, or chain-integrity audit.",
            "theoretical_limit": "The certificate derives the finite sign pattern from the selected zero-carrier partition; it does not derive omega/phi or transport from strict nadsoliton dynamics.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in solitonic state; this audit only checks finite phase-cell sign bookkeeping.",
            "forbidden_reading": "No separate informational layer below the nadsoliton is introduced.",
        },
        "hard_limits": [
            "K_strict_gate remains the current live/full operational kernel.",
            "No unqualified identity K_legacy_ont == K_strict_gate is claimed.",
            "No proof derives strict omega/phi or phase transport from strict nadsoliton dynamics.",
            "No beta_tors -> chi_11 theorem is claimed.",
            "No legacy physical-role transfer to K_strict_gate is claimed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    summary = payload["cell_sign_summary"]
    lines = [
        "# Strict completion phase-zero cell-sign certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This audit derives the integer-node phase sign pattern from the ordered",
        "rational zero-carrier partition and a left anchor sign.  It uses only",
        "rational carrier comparisons and crossing parity, not fresh trigonometric",
        "evaluation of cosine values.",
        "",
        "## Sign-rule summary",
        "",
        f"- Zero carrier order: `{payload['sign_rule']['zero_carrier_order']}`",
        f"- Left anchor sign at d=0: `{payload['sign_rule']['left_anchor_sign_at_d0']}`",
        f"- Uses trig evaluation: `{payload['sign_rule']['uses_trig_evaluation']}`",
        f"- All nodes outside zero carriers: `{summary['all_nodes_outside_zero_carriers']}`",
        f"- All node signs match expected: `{summary['all_node_signs_match_expected']}`",
        f"- All edge flips match crossing parity: `{summary['all_edge_flips_match_crossing_parity']}`",
        f"- Maximum crossed zero-carriers per integer edge: `{summary['max_crossed_zero_carriers_per_integer_edge']}`",
        f"- Derived sign pattern: `{summary['derived_phase_transport_sign_pattern']}`",
        f"- Derived flip edges: `{summary['derived_phase_sign_flip_edges']}`",
        "",
        "## Node sign rows",
        "",
        "| d | zero carriers left | parity | derived sign | expected sign | matches |",
        "|---:|---|---|---:|---:|---:|",
    ]
    for row in payload["node_sign_rows"]:
        lines.append(
            f"| {row['d']} | {row['zero_carriers_left_of_node']} | {row['left_count_parity']} | {row['derived_phase_transport_sign']} | {row['expected_phase_transport_sign']} | {row['matches_expected_sign']} |"
        )
    lines.extend([
        "",
        "## Edge crossing rows",
        "",
        "| edge | crossed carriers | count | sign flip | parity match |",
        "|---|---|---:|---:|---:|",
    ])
    for row in payload["edge_sign_rows"]:
        lines.append(
            f"| {row['edge']} | {row['new_zero_carriers_crossed']} | {row['crossed_zero_carrier_count']} | {row['sign_flip']} | {row['matches_crossing_parity']} |"
        )
    lines.extend([
        "",
        "## Proof certificate",
        "",
    ])
    for key, value in payload["proof_certificate"].items():
        lines.append(f"- `{key}`: {value}")
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
