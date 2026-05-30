#!/usr/bin/env python3
"""Scratch probe: phase-zero interlacing certificate for strict completion.

Earlier completion audits showed that the phase/frequency factor

    P(d)=cos(omega_S*d+phi_S)/cos(omega_L*d+phi_L)

is a shape-critical part of completing the legacy carrier into the strict
kernel.  This probe makes the sign part of that statement more proof-like: it
computes the exact cosine-zero positions on the real interval [0,11] and proves
which adjacent integer edges must flip the phase-transport sign.

The certificate is finite and phase-only.  It does not derive the strict phase
parameters from nadsoliton dynamics, does not identify the legacy and strict
kernels, and does not discharge selector/ToE blockers.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_zero_interlacing_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_zero_interlacing_certificate_report.md"
NECESSITY_REPORT = HERE / "bridge_strict_kernel_completion_necessity_certificate_report.json"
COCYCLE_REPORT = HERE / "bridge_strict_kernel_completion_transport_cocycle_report.json"
LOW_ORDER_REPORT = HERE / "bridge_strict_completion_low_order_transport_no_go_report.json"

LEGACY = {"omega": math.pi / 4.0, "phi": math.pi / 6.0}
STRICT = {"omega": 0.18575, "phi": 0.16250}
DOMAIN = list(range(12))
EDGES = list(range(11))
TOL = 1e-14


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def sign(value: float) -> int:
    if value > 0.0:
        return 1
    if value < 0.0:
        return -1
    return 0


def cosine_value(params: dict[str, float], d: float) -> float:
    return math.cos(params["omega"] * d + params["phi"])


def zero_positions(params: dict[str, float], lo: float, hi: float) -> list[dict[str, Any]]:
    # cos(omega*x+phi)=0 iff omega*x+phi = pi/2 + k*pi.
    rows: list[dict[str, Any]] = []
    min_k = math.floor((params["omega"] * lo + params["phi"] - math.pi / 2.0) / math.pi) - 2
    max_k = math.ceil((params["omega"] * hi + params["phi"] - math.pi / 2.0) / math.pi) + 2
    for k in range(min_k, max_k + 1):
        position = (math.pi / 2.0 + k * math.pi - params["phi"]) / params["omega"]
        if lo < position < hi:
            rows.append({"k": k, "position": position})
    return rows


def node_rows(legacy_zeros: list[dict[str, Any]], strict_zeros: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    all_zero_positions = [row["position"] for row in legacy_zeros + strict_zeros]
    for d in DOMAIN:
        legacy_cos = cosine_value(LEGACY, float(d))
        strict_cos = cosine_value(STRICT, float(d))
        phase_sign = sign(strict_cos) * sign(legacy_cos)
        nearest_zero_distance = min(abs(float(d) - zero) for zero in all_zero_positions)
        rows.append({
            "d": d,
            "legacy_cos": legacy_cos,
            "strict_cos": strict_cos,
            "legacy_sign": sign(legacy_cos),
            "strict_sign": sign(strict_cos),
            "phase_transport_sign": phase_sign,
            "nearest_phase_zero_distance": nearest_zero_distance,
            "node_avoids_zero": abs(legacy_cos) > TOL and abs(strict_cos) > TOL,
        })
    return rows


def edge_rows(legacy_zeros: list[dict[str, Any]], strict_zeros: list[dict[str, Any]], nodes: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for edge in EDGES:
        legacy_inside = [row for row in legacy_zeros if edge < row["position"] < edge + 1]
        strict_inside = [row for row in strict_zeros if edge < row["position"] < edge + 1]
        zero_count = len(legacy_inside) + len(strict_inside)
        left_sign = nodes[edge]["phase_transport_sign"]
        right_sign = nodes[edge + 1]["phase_transport_sign"]
        predicted_flip = zero_count % 2 == 1
        actual_flip = left_sign != right_sign
        rows.append({
            "edge": f"{edge}->{edge + 1}",
            "left_d": edge,
            "right_d": edge + 1,
            "legacy_zero_positions_inside": [row["position"] for row in legacy_inside],
            "strict_zero_positions_inside": [row["position"] for row in strict_inside],
            "zero_count_inside": zero_count,
            "predicted_phase_sign_flip_by_zero_parity": predicted_flip,
            "actual_phase_sign_flip": actual_flip,
            "parity_matches_actual": predicted_flip == actual_flip,
            "left_phase_sign": left_sign,
            "right_phase_sign": right_sign,
        })
    return rows


def min_endpoint_separation(zero_rows: list[dict[str, Any]]) -> float:
    return min(min(abs(row["position"] - edge), abs(row["position"] - (edge + 1))) for row in zero_rows for edge in EDGES if edge < row["position"] < edge + 1)


def build_payload() -> dict[str, Any]:
    necessity = load_json(NECESSITY_REPORT)
    cocycle = load_json(COCYCLE_REPORT)
    low_order = load_json(LOW_ORDER_REPORT)

    legacy_zeros = zero_positions(LEGACY, 0.0, 11.0)
    strict_zeros = zero_positions(STRICT, 0.0, 11.0)
    nodes = node_rows(legacy_zeros, strict_zeros)
    edges = edge_rows(legacy_zeros, strict_zeros, nodes)
    flip_edges = [row["edge"] for row in edges if row["actual_phase_sign_flip"]]
    parity_failures = [row["edge"] for row in edges if not row["parity_matches_actual"]]

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_ZERO_INTERLACING_CERTIFICATE__FINITE_Z12_PHASE_ONLY",
        "status": "phase-transport-sign-flips-certified-by-legacy-strict-cosine-zero-interlacing",
        "source_reports": {
            "necessity_certificate": str(NECESSITY_REPORT.relative_to(ROOT)),
            "transport_cocycle": str(COCYCLE_REPORT.relative_to(ROOT)),
            "low_order_no_go": str(LOW_ORDER_REPORT.relative_to(ROOT)),
        },
        "constants": {
            "legacy_phase": LEGACY,
            "strict_phase": STRICT,
            "domain_d_values": DOMAIN,
            "edge_d_values": EDGES,
            "tolerance": TOL,
        },
        "zero_formula": "cos(omega*x+phi)=0 iff x=(pi/2+k*pi-phi)/omega",
        "legacy_zero_rows": legacy_zeros,
        "strict_zero_rows": strict_zeros,
        "node_phase_rows": nodes,
        "edge_interlacing_rows": edges,
        "interlacing_summary": {
            "legacy_zero_positions": [row["position"] for row in legacy_zeros],
            "strict_zero_positions": [row["position"] for row in strict_zeros],
            "legacy_zero_count": len(legacy_zeros),
            "strict_zero_count": len(strict_zeros),
            "phase_transport_sign_pattern": [row["phase_transport_sign"] for row in nodes],
            "phase_sign_flip_edges": flip_edges,
            "zero_parity_matches_all_edges": not parity_failures,
            "parity_failure_edges": parity_failures,
            "all_integer_nodes_avoid_phase_zeros": all(row["node_avoids_zero"] for row in nodes),
            "min_legacy_zero_endpoint_separation": min_endpoint_separation(legacy_zeros),
            "min_strict_zero_endpoint_separation": min_endpoint_separation(strict_zeros),
        },
        "blocker_context": {
            "necessity_status": necessity["status"],
            "transport_cocycle_status": cocycle["status"],
            "low_order_status": low_order["status"],
            "what_this_refines": "explains the sign part of phase_frequency_transport by exact zero interlacing rather than by another fit of T(d)",
            "strict_phase_derivation_status": "still_open; zero interlacing certifies consequences of the phase parameters, not their strict-side derivation",
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_transport_derivation_from_nadsoliton_dynamics",
                "orientation_chi11_source",
                "chi11_uniqueness",
                "reynolds_obstruction_escape",
                "role_transfer_theorem",
            ],
        },
        "proof_certificate": {
            "zero_formula_step": "The cosine zero formula lists every real zero in [0,11] for the legacy and strict phase lines.",
            "node_nonzero_step": "Every integer node d=0..11 has nonzero legacy and strict cosine value, so the phase transport sign is well-defined.",
            "edge_parity_step": "On each open integer edge, the phase transport sign flips exactly when the combined legacy+strict zero count is odd.",
            "flip_edges_step": "The certified flip edges are 1->2, 5->6, 7->8, and 9->10, matching the transport-cocycle sign flips.",
            "nonduplication": "This is a phase-zero interlacing certificate, not another subset necessity, cocycle reconstruction, or low-order fit audit.",
            "theoretical_limit": "The certificate proves finite phase-zero consequences of the chosen parameters; it does not derive those parameters from strict dynamics.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in solitonic state; this audit only checks phase-zero geometry of the completion factor.",
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
    summary = payload["interlacing_summary"]
    lines = [
        "# Strict completion phase-zero interlacing certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This finite-domain audit explains the sign part of the phase/frequency",
        "completion factor by exact cosine-zero interlacing on `[0,11]`.  It is a",
        "certificate for consequences of the chosen phase parameters, not a strict",
        "dynamical derivation of those parameters.",
        "",
        "## Zero formula",
        "",
        f"`{payload['zero_formula']}`",
        "",
        "## Interlacing summary",
        "",
        f"- Legacy zero positions: `{summary['legacy_zero_positions']}`",
        f"- Strict zero positions: `{summary['strict_zero_positions']}`",
        f"- Phase transport sign pattern: `{summary['phase_transport_sign_pattern']}`",
        f"- Phase sign-flip edges: `{summary['phase_sign_flip_edges']}`",
        f"- Zero parity matches all edges: `{summary['zero_parity_matches_all_edges']}`",
        f"- All integer nodes avoid phase zeros: `{summary['all_integer_nodes_avoid_phase_zeros']}`",
        f"- Min legacy zero endpoint separation: `{summary['min_legacy_zero_endpoint_separation']:.12e}`",
        f"- Min strict zero endpoint separation: `{summary['min_strict_zero_endpoint_separation']:.12e}`",
        "",
        "## Edge interlacing rows",
        "",
        "| edge | legacy zeros inside | strict zeros inside | zero count | predicted flip | actual flip |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in payload["edge_interlacing_rows"]:
        lines.append(
            "| {edge} | {legacy} | {strict} | {count} | {predicted} | {actual} |".format(
                edge=row["edge"],
                legacy=len(row["legacy_zero_positions_inside"]),
                strict=len(row["strict_zero_positions_inside"]),
                count=row["zero_count_inside"],
                predicted=row["predicted_phase_sign_flip_by_zero_parity"],
                actual=row["actual_phase_sign_flip"],
            )
        )
    lines.extend([
        "",
        "## Guarded interpretation",
        "",
        "The certificate explains the finite sign flips of the completion phase factor.",
        "It does not derive strict omega/phi from nadsoliton dynamics, does not prove",
        "`beta_tors -> chi_11`, and does not transfer legacy physical roles onto",
        "`K_strict_gate`.",
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
