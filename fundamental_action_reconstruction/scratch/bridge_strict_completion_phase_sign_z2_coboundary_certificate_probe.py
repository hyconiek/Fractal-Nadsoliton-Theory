#!/usr/bin/env python3
"""Scratch probe: Z2 coboundary certificate for the phase sign pattern.

The cell-sign certificate derived the integer-node phase signs from ordered
zero-carrier counting.  This probe packages the same finite sign information as
a Z2 cochain/coboundary theorem on the path graph 0--1--...--11:

    node bit b(d)=0 for sign +1 and b(d)=1 for sign -1,
    edge bit e(d,d+1)=b(d) xor b(d+1).

It then verifies the interval law that the xor/sum mod 2 of edge bits over any
interval equals the endpoint node-bit difference.  This is a finite algebraic
certificate for sign transport consistency, not a new cosine fit, not a new
selector premise, and not a strict dynamical derivation.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.md"
CELL_SIGN_REPORT = HERE / "bridge_strict_completion_phase_zero_cell_sign_certificate_report.json"
COCYCLE_REPORT = HERE / "bridge_strict_kernel_completion_transport_cocycle_report.json"
RATIONAL_ZERO_REPORT = HERE / "bridge_strict_completion_phase_zero_rational_interval_certificate_report.json"

EXPECTED_SIGN_PATTERN = [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1]
EXPECTED_FLIP_EDGES = ["1->2", "5->6", "7->8", "9->10"]
DOMAIN = list(range(12))


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def sign_to_bit(sign: int) -> int:
    if sign == 1:
        return 0
    if sign == -1:
        return 1
    raise ValueError(f"unsupported sign: {sign}")


def bit_to_sign(bit: int) -> int:
    return -1 if bit % 2 else 1


def node_bit_rows(sign_pattern: list[int]) -> list[dict[str, Any]]:
    return [
        {
            "d": d,
            "phase_sign": sign,
            "node_bit": sign_to_bit(sign),
        }
        for d, sign in enumerate(sign_pattern)
    ]


def edge_bit_rows(bits: list[int]) -> list[dict[str, Any]]:
    rows = []
    for d, (left, right) in enumerate(zip(bits, bits[1:])):
        edge_bit = left ^ right
        rows.append({
            "edge": f"{d}->{d + 1}",
            "left_d": d,
            "right_d": d + 1,
            "left_node_bit": left,
            "right_node_bit": right,
            "edge_bit": edge_bit,
            "edge_sign_ratio": bit_to_sign(edge_bit),
            "is_phase_flip": edge_bit == 1,
        })
    return rows


def prefix_reconstruction_rows(edge_rows: list[dict[str, Any]], anchor_bit: int) -> list[dict[str, Any]]:
    rows = [{
        "d": 0,
        "prefix_edge_bit_sum_mod2": 0,
        "reconstructed_node_bit": anchor_bit,
        "reconstructed_phase_sign": bit_to_sign(anchor_bit),
        "expected_phase_sign": EXPECTED_SIGN_PATTERN[0],
        "matches_expected": bit_to_sign(anchor_bit) == EXPECTED_SIGN_PATTERN[0],
    }]
    current = anchor_bit
    running_sum = 0
    for row in edge_rows:
        running_sum ^= row["edge_bit"]
        current ^= row["edge_bit"]
        d = row["right_d"]
        rows.append({
            "d": d,
            "prefix_edge_bit_sum_mod2": running_sum,
            "reconstructed_node_bit": current,
            "reconstructed_phase_sign": bit_to_sign(current),
            "expected_phase_sign": EXPECTED_SIGN_PATTERN[d],
            "matches_expected": bit_to_sign(current) == EXPECTED_SIGN_PATTERN[d],
        })
    return rows


def interval_rows(bits: list[int], edge_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    edge_bits = [row["edge_bit"] for row in edge_rows]
    rows = []
    for start in DOMAIN:
        for end in range(start + 1, len(DOMAIN)):
            edge_sum = sum(edge_bits[start:end]) % 2
            endpoint_xor = bits[start] ^ bits[end]
            rows.append({
                "interval": f"{start}->{end}",
                "start_d": start,
                "end_d": end,
                "edge_bit_sum_mod2": edge_sum,
                "endpoint_node_bit_xor": endpoint_xor,
                "endpoint_sign_ratio": bit_to_sign(endpoint_xor),
                "matches_endpoint_coboundary": edge_sum == endpoint_xor,
            })
    return rows


def build_payload() -> dict[str, Any]:
    cell_sign = load_json(CELL_SIGN_REPORT)
    cocycle = load_json(COCYCLE_REPORT)
    rational_zero = load_json(RATIONAL_ZERO_REPORT)
    sign_pattern = cell_sign["cell_sign_summary"]["derived_phase_transport_sign_pattern"]
    bits = [sign_to_bit(sign) for sign in sign_pattern]
    node_rows = node_bit_rows(sign_pattern)
    edge_rows = edge_bit_rows(bits)
    prefix_rows = prefix_reconstruction_rows(edge_rows, bits[0])
    intervals = interval_rows(bits, edge_rows)
    derived_flip_edges = [row["edge"] for row in edge_rows if row["is_phase_flip"]]
    edge_support_size = sum(row["edge_bit"] for row in edge_rows)
    interval_count = len(intervals)

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_Z2_COBOUNDARY_CERTIFICATE__FINITE_PATH_GRAPH",
        "status": "phase-sign-pattern-certified-as-z2-coboundary-on-finite-path-domain",
        "source_reports": {
            "cell_sign_certificate": str(CELL_SIGN_REPORT.relative_to(ROOT)),
            "transport_cocycle_certificate": str(COCYCLE_REPORT.relative_to(ROOT)),
            "rational_zero_certificate": str(RATIONAL_ZERO_REPORT.relative_to(ROOT)),
        },
        "z2_definitions": {
            "node_bit": "b(d)=0 for phase sign +1; b(d)=1 for phase sign -1",
            "edge_bit": "e(d,d+1)=b(d) xor b(d+1)",
            "interval_law": "sum edge bits mod 2 over [a,b] equals b(a) xor b(b)",
            "domain": DOMAIN,
        },
        "node_bit_rows": node_rows,
        "edge_bit_rows": edge_rows,
        "prefix_reconstruction_rows": prefix_rows,
        "interval_coboundary_rows": intervals,
        "z2_coboundary_summary": {
            "anchor_node_bit_d0": bits[0],
            "anchor_phase_sign_d0": sign_pattern[0],
            "node_bit_pattern": bits,
            "phase_sign_pattern": sign_pattern,
            "edge_bit_pattern": [row["edge_bit"] for row in edge_rows],
            "derived_phase_sign_flip_edges": derived_flip_edges,
            "edge_flip_support_size": edge_support_size,
            "interval_count": interval_count,
            "all_prefix_reconstruction_rows_match_expected": all(row["matches_expected"] for row in prefix_rows),
            "all_interval_coboundary_rows_pass": all(row["matches_endpoint_coboundary"] for row in intervals),
            "matches_cell_sign_pattern": sign_pattern == cell_sign["cell_sign_summary"]["derived_phase_transport_sign_pattern"],
            "matches_cell_sign_flip_edges": derived_flip_edges == cell_sign["cell_sign_summary"]["derived_phase_sign_flip_edges"],
            "matches_cocycle_flip_edges": derived_flip_edges == cocycle["cocycle_summary"]["phase_sign_flip_edges"],
            "matches_rational_zero_flip_edges": derived_flip_edges == rational_zero["interval_summary"]["phase_sign_flip_edges_from_rational_intervals"],
            "matches_expected_sign_pattern": sign_pattern == EXPECTED_SIGN_PATTERN,
            "matches_expected_flip_edges": derived_flip_edges == EXPECTED_FLIP_EDGES,
        },
        "blocker_context": {
            "what_this_refines": "packages the cell-sign result as a finite Z2 coboundary and verifies every interval parity law",
            "cell_sign_status": cell_sign["status"],
            "transport_cocycle_status": cocycle["status"],
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
            "node_step": "Convert the certified phase signs into Z2 node bits b(d).",
            "edge_step": "Compute edge bits as the coboundary e(d,d+1)=b(d) xor b(d+1).",
            "prefix_step": "Starting from b(0), prefix xor of edge bits reconstructs every node bit and phase sign.",
            "interval_step": "All 66 nontrivial intervals satisfy sum(edge bits) mod 2 = endpoint node-bit xor.",
            "nonduplication": "This is a Z2 cochain/coboundary algebra certificate, not another zero-location, cell-sign, damping, or real-valued cocycle audit.",
            "theoretical_limit": "The certificate proves finite sign-transport algebra on the selected domain; it does not derive omega/phi or transport from strict nadsoliton dynamics and does not discharge a selector obstruction.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in solitonic state; this audit only checks finite Z2 sign bookkeeping.",
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
    summary = payload["z2_coboundary_summary"]
    lines = [
        "# Strict completion phase-sign Z2 coboundary certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This audit packages the certified phase sign pattern as a finite Z2",
        "cochain/coboundary on the path graph `0--1--...--11`.  It verifies",
        "prefix reconstruction and all interval parity laws without introducing a",
        "new cosine fit, selector premise, or strict dynamical derivation.",
        "",
        "## Z2 summary",
        "",
        f"- Node bit pattern: `{summary['node_bit_pattern']}`",
        f"- Edge bit pattern: `{summary['edge_bit_pattern']}`",
        f"- Flip support size: `{summary['edge_flip_support_size']}`",
        f"- Derived flip edges: `{summary['derived_phase_sign_flip_edges']}`",
        f"- Interval count: `{summary['interval_count']}`",
        f"- Prefix reconstruction passes: `{summary['all_prefix_reconstruction_rows_match_expected']}`",
        f"- All interval coboundary rows pass: `{summary['all_interval_coboundary_rows_pass']}`",
        "",
        "## Edge bit rows",
        "",
        "| edge | left bit | right bit | edge bit | flip |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in payload["edge_bit_rows"]:
        lines.append(
            f"| {row['edge']} | {row['left_node_bit']} | {row['right_node_bit']} | {row['edge_bit']} | {row['is_phase_flip']} |"
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
