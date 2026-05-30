#!/usr/bin/env python3
"""Scratch probe: rational-interval phase-zero certificate for strict completion.

The phase-zero interlacing certificate located the cosine zeros using floating
arithmetic.  This probe makes the same sign-flip result more proof-style by
replacing the decisive zero-location step with rational interval arithmetic.

Inputs treated as exact rationals where they are specified by terminating
constants:

    omega_S = 0.18575 = 743/4000,   phi_S = 0.16250 = 13/80.

For pi, the probe uses the classical rational enclosure

    333/106 < pi < 355/113.

It then brackets strict cosine zeros x_k=(pi/2+k*pi-phi_S)/omega_S and proves
that on [0,11] only k=0 occurs, with x_0 in (7,8).  Legacy zeros are exact at
4/3, 16/3, and 28/3.  Therefore the phase sign flips on integer edges are
forced by zero parity without relying on floating zero placement.

This is still a finite phase certificate.  It does not derive strict omega/phi
from nadsoliton dynamics and does not claim a kernel identity or selector/ToE
closure.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_zero_rational_interval_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_zero_rational_interval_certificate_report.md"
PHASE_ZERO_REPORT = HERE / "bridge_strict_completion_phase_zero_interlacing_certificate_report.json"
DAMPING_REPORT = HERE / "bridge_strict_completion_damping_continuous_monotonicity_certificate_report.json"
NECESSITY_REPORT = HERE / "bridge_strict_kernel_completion_necessity_certificate_report.json"

PI_LOWER = Fraction(333, 106)
PI_UPPER = Fraction(355, 113)
STRICT_OMEGA = Fraction(743, 4000)
STRICT_PHI = Fraction(13, 80)
DOMAIN = list(range(12))
EDGES = list(range(11))
LEGACY_ZERO_POSITIONS = [Fraction(4, 3), Fraction(16, 3), Fraction(28, 3)]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def fraction_payload(value: Fraction) -> dict[str, Any]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
        "text": f"{value.numerator}/{value.denominator}",
    }


def affine_pi_interval(coefficient: Fraction, constant: Fraction) -> tuple[Fraction, Fraction]:
    lower = coefficient * PI_LOWER + constant
    upper = coefficient * PI_UPPER + constant
    if coefficient < 0:
        lower = coefficient * PI_UPPER + constant
        upper = coefficient * PI_LOWER + constant
    return lower, upper


def strict_zero_interval(k_value: int) -> dict[str, Any]:
    coefficient = Fraction(1, 2) + k_value
    numerator_lower, numerator_upper = affine_pi_interval(coefficient, -STRICT_PHI)
    x_lower = numerator_lower / STRICT_OMEGA
    x_upper = numerator_upper / STRICT_OMEGA
    return {
        "k": k_value,
        "numerator_interval": {
            "lower": fraction_payload(numerator_lower),
            "upper": fraction_payload(numerator_upper),
        },
        "x_interval": {
            "lower": fraction_payload(x_lower),
            "upper": fraction_payload(x_upper),
        },
        "proves_left_of_domain": x_upper < 0,
        "proves_inside_7_8": Fraction(7, 1) < x_lower and x_upper < Fraction(8, 1),
        "proves_right_of_audit_domain": x_lower > 11,
        "proves_in_audit_domain": Fraction(0, 1) < x_lower and x_upper < Fraction(11, 1),
    }


def legacy_zero_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for index, position in enumerate(LEGACY_ZERO_POSITIONS):
        rows.append({
            "k": index,
            "position": fraction_payload(position),
            "edge_containing_zero": f"{position.numerator // position.denominator}->{position.numerator // position.denominator + 1}",
            "proves_inside_audit_domain": 0 < position < 11,
            "proves_not_integer_node": position.denominator != 1,
        })
    return rows


def certified_zero_counts_by_edge(strict_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    strict_inside_edges = []
    for row in strict_rows:
        if row["proves_inside_7_8"]:
            strict_inside_edges.append("7->8")
    rows: list[dict[str, Any]] = []
    for edge in EDGES:
        edge_name = f"{edge}->{edge + 1}"
        legacy_count = sum(1 for position in LEGACY_ZERO_POSITIONS if edge < position < edge + 1)
        strict_count = strict_inside_edges.count(edge_name)
        zero_count = legacy_count + strict_count
        rows.append({
            "edge": edge_name,
            "legacy_zero_count": legacy_count,
            "strict_zero_count": strict_count,
            "total_zero_count": zero_count,
            "phase_sign_flip_by_odd_parity": zero_count % 2 == 1,
        })
    return rows


def phase_sign_pattern_from_parity(edge_rows: list[dict[str, Any]]) -> list[int]:
    signs = [1]
    current = 1
    for row in edge_rows:
        if row["phase_sign_flip_by_odd_parity"]:
            current *= -1
        signs.append(current)
    return signs


def build_payload() -> dict[str, Any]:
    phase_zero = load_json(PHASE_ZERO_REPORT)
    damping = load_json(DAMPING_REPORT)
    necessity = load_json(NECESSITY_REPORT)

    strict_rows = [strict_zero_interval(k_value) for k_value in [-1, 0, 1]]
    legacy_rows = legacy_zero_rows()
    edge_rows = certified_zero_counts_by_edge(strict_rows)
    sign_pattern = phase_sign_pattern_from_parity(edge_rows)
    flip_edges = [row["edge"] for row in edge_rows if row["phase_sign_flip_by_odd_parity"]]

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_ZERO_RATIONAL_INTERVAL_CERTIFICATE__NO_FLOAT_ZERO_PLACEMENT",
        "status": "phase-zero-sign-flip-edges-certified-by-rational-pi-intervals",
        "source_reports": {
            "phase_zero_interlacing_float_certificate": str(PHASE_ZERO_REPORT.relative_to(ROOT)),
            "damping_monotonicity_certificate": str(DAMPING_REPORT.relative_to(ROOT)),
            "necessity_certificate": str(NECESSITY_REPORT.relative_to(ROOT)),
        },
        "rational_inputs": {
            "pi_lower_bound": fraction_payload(PI_LOWER),
            "pi_upper_bound": fraction_payload(PI_UPPER),
            "strict_omega": fraction_payload(STRICT_OMEGA),
            "strict_phi": fraction_payload(STRICT_PHI),
            "legacy_zero_positions_exact": [fraction_payload(value) for value in LEGACY_ZERO_POSITIONS],
        },
        "strict_zero_interval_rows": strict_rows,
        "legacy_zero_rows": legacy_rows,
        "edge_zero_parity_rows": edge_rows,
        "interval_summary": {
            "strict_k_minus_1_left_of_domain": strict_rows[0]["proves_left_of_domain"],
            "strict_k_0_inside_edge_7_8": strict_rows[1]["proves_inside_7_8"],
            "strict_k_1_right_of_audit_domain": strict_rows[2]["proves_right_of_audit_domain"],
            "strict_zero_count_in_0_11": sum(1 for row in strict_rows if row["proves_in_audit_domain"]),
            "legacy_zero_count_in_0_11": len(LEGACY_ZERO_POSITIONS),
            "phase_sign_flip_edges_from_rational_intervals": flip_edges,
            "phase_transport_sign_pattern_from_rational_intervals": sign_pattern,
            "matches_float_phase_zero_report_flip_edges": flip_edges == phase_zero["interlacing_summary"]["phase_sign_flip_edges"],
            "matches_float_phase_zero_report_sign_pattern": sign_pattern == phase_zero["interlacing_summary"]["phase_transport_sign_pattern"],
            "all_integer_nodes_certified_away_from_zeros": all(row["proves_not_integer_node"] for row in legacy_rows) and strict_rows[1]["proves_inside_7_8"],
        },
        "blocker_context": {
            "phase_zero_status": phase_zero["status"],
            "damping_status": damping["status"],
            "necessity_status": necessity["status"],
            "what_this_refines": "replaces the decisive strict zero placement in the phase sign certificate with rational interval arithmetic using 333/106 < pi < 355/113",
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_transport_derivation_from_nadsoliton_dynamics",
                "orientation_chi11_source",
                "chi11_uniqueness",
                "role_transfer_theorem",
            ],
        },
        "proof_certificate": {
            "pi_interval_step": "Use 333/106 < pi < 355/113 and exact rational strict parameters omega=743/4000, phi=13/80.",
            "strict_zero_step": "The strict zero intervals prove k=-1 lies left of 0, k=0 lies inside (7,8), and k=1 lies right of 11.",
            "legacy_zero_step": "Legacy phase zeros are exact at 4/3, 16/3, and 28/3, hence inside edges 1->2, 5->6, and 9->10.",
            "parity_step": "Starting from positive phase sign at d=0, odd zero parity on an edge forces exactly the certified sign flips.",
            "nonduplication": "This is a rational-interval strengthening of phase-zero placement, not another damping, cocycle, or low-order fit audit.",
            "theoretical_limit": "The certificate proves zero placement consequences of the selected phase parameters; it does not derive those parameters from strict dynamics.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in solitonic state; this audit only checks interval placement of phase zeros.",
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
    summary = payload["interval_summary"]
    inputs = payload["rational_inputs"]
    lines = [
        "# Strict completion phase-zero rational interval certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This audit strengthens the phase-zero sign-flip certificate by using rational",
        "interval arithmetic instead of relying on floating zero placement for the",
        "strict phase zero.",
        "",
        "## Rational inputs",
        "",
        f"- Pi interval: `{inputs['pi_lower_bound']['text']} < pi < {inputs['pi_upper_bound']['text']}`",
        f"- Strict omega: `{inputs['strict_omega']['text']}`",
        f"- Strict phi: `{inputs['strict_phi']['text']}`",
        f"- Exact legacy zeros: `{[row['text'] for row in inputs['legacy_zero_positions_exact']]}`",
        "",
        "## Interval summary",
        "",
        f"- Strict k=-1 left of domain: `{summary['strict_k_minus_1_left_of_domain']}`",
        f"- Strict k=0 inside edge 7->8: `{summary['strict_k_0_inside_edge_7_8']}`",
        f"- Strict k=1 right of audit domain: `{summary['strict_k_1_right_of_audit_domain']}`",
        f"- Strict zero count in (0,11): `{summary['strict_zero_count_in_0_11']}`",
        f"- Legacy zero count in (0,11): `{summary['legacy_zero_count_in_0_11']}`",
        f"- Flip edges from rational intervals: `{summary['phase_sign_flip_edges_from_rational_intervals']}`",
        f"- Sign pattern from rational intervals: `{summary['phase_transport_sign_pattern_from_rational_intervals']}`",
        f"- Matches float phase-zero report flip edges: `{summary['matches_float_phase_zero_report_flip_edges']}`",
        f"- Matches float phase-zero report sign pattern: `{summary['matches_float_phase_zero_report_sign_pattern']}`",
        "",
        "## Strict zero interval rows",
        "",
        "| k | x lower | x upper | placement |",
        "|---:|---:|---:|---|",
    ]
    for row in payload["strict_zero_interval_rows"]:
        if row["proves_left_of_domain"]:
            placement = "left of 0"
        elif row["proves_inside_7_8"]:
            placement = "inside 7->8"
        elif row["proves_right_of_audit_domain"]:
            placement = "right of 11"
        else:
            placement = "unclassified"
        lines.append(
            f"| {row['k']} | {row['x_interval']['lower']['text']} | {row['x_interval']['upper']['text']} | {placement} |"
        )
    lines.extend([
        "",
        "## Guarded interpretation",
        "",
        "This is only a rational-interval certificate for phase-zero placement. It does",
        "not derive strict omega/phi from nadsoliton dynamics, does not prove",
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
