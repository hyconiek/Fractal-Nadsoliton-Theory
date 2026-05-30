#!/usr/bin/env python3
"""Scratch probe: continuous damping-compression monotonicity certificate.

The completion ladder/necessity/cocycle audits separate the strict completion
factor into amplitude, phase/frequency, and damping/compression pieces.  The
phase-zero interlacing certificate explains the sign flips.  This probe gives
the complementary proof-style audit for the damping piece

    D(x) = (1 + beta_tors*x)/(1 + beta*x^eta),  beta_tors=0.01, beta=1, eta=9/5.

It proves on the continuous interval [1,11] that D(x) is positive and strictly
decreasing, using the closed-form derivative numerator

    N(x) = beta_tors - eta*x^(eta-1) + beta_tors*(1-eta)*x^eta.

For x>=1, N(x) <= 0.01 - 1.8 < 0, while the denominator is positive.  Therefore
D'(x)<0 on [1,11].  This is a finite/elementary calculus certificate for the
completion factor's envelope role, not a derivation from strict nadsoliton
dynamics and not a bridge theorem.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_damping_continuous_monotonicity_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_damping_continuous_monotonicity_certificate_report.md"
NECESSITY_REPORT = HERE / "bridge_strict_kernel_completion_necessity_certificate_report.json"
COCYCLE_REPORT = HERE / "bridge_strict_kernel_completion_transport_cocycle_report.json"
PHASE_ZERO_REPORT = HERE / "bridge_strict_completion_phase_zero_interlacing_certificate_report.json"
LOW_ORDER_REPORT = HERE / "bridge_strict_completion_low_order_transport_no_go_report.json"

BETA_TORS = 0.01
BETA = 1.0
ETA = 9.0 / 5.0
DOMAIN = list(range(12))
CERT_INTERVAL = [1.0, 11.0]
TOL = 1e-14


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def damping_factor(x_value: float) -> float:
    return (1.0 + BETA_TORS * x_value) / (1.0 + BETA * x_value ** ETA)


def derivative_numerator(x_value: float) -> float:
    return BETA_TORS - ETA * x_value ** (ETA - 1.0) + BETA_TORS * (1.0 - ETA) * x_value ** ETA


def derivative_denominator(x_value: float) -> float:
    return (1.0 + BETA * x_value ** ETA) ** 2


def damping_derivative(x_value: float) -> float:
    return derivative_numerator(x_value) / derivative_denominator(x_value)


def proof_grid_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for d_value in range(int(CERT_INTERVAL[0]), int(CERT_INTERVAL[1]) + 1):
        x_value = float(d_value)
        rows.append({
            "x": x_value,
            "D_x": damping_factor(x_value),
            "derivative_numerator_N_x": derivative_numerator(x_value),
            "derivative_denominator_positive": derivative_denominator(x_value) > 0.0,
            "D_prime_x": damping_derivative(x_value),
            "N_x_negative": derivative_numerator(x_value) < 0.0,
            "D_prime_x_negative": damping_derivative(x_value) < 0.0,
        })
    return rows


def discrete_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for left in range(int(CERT_INTERVAL[0]), int(CERT_INTERVAL[1])):
        right = left + 1
        left_value = damping_factor(float(left))
        right_value = damping_factor(float(right))
        rows.append({
            "edge": f"{left}->{right}",
            "D_left": left_value,
            "D_right": right_value,
            "drop": left_value - right_value,
            "ratio_left_over_right": left_value / right_value,
            "strict_drop_positive": left_value > right_value,
        })
    return rows


def build_payload() -> dict[str, Any]:
    necessity = load_json(NECESSITY_REPORT)
    cocycle = load_json(COCYCLE_REPORT)
    phase_zero = load_json(PHASE_ZERO_REPORT)
    low_order = load_json(LOW_ORDER_REPORT)

    grid = proof_grid_rows()
    discrete = discrete_rows()
    endpoint_left = damping_factor(CERT_INTERVAL[0])
    endpoint_right = damping_factor(CERT_INTERVAL[1])
    derivative_upper_bound_numerator = BETA_TORS - ETA

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_DAMPING_CONTINUOUS_MONOTONICITY_CERTIFICATE__CALCULUS_ENVELOPE_ONLY",
        "status": "damping-compression-factor-positive-and-strictly-decreasing-on-continuous-interval-1-to-11",
        "source_reports": {
            "necessity_certificate": str(NECESSITY_REPORT.relative_to(ROOT)),
            "transport_cocycle": str(COCYCLE_REPORT.relative_to(ROOT)),
            "phase_zero_interlacing": str(PHASE_ZERO_REPORT.relative_to(ROOT)),
            "low_order_no_go": str(LOW_ORDER_REPORT.relative_to(ROOT)),
        },
        "constants": {
            "beta_tors": BETA_TORS,
            "beta": BETA,
            "eta": ETA,
            "certified_continuous_interval": CERT_INTERVAL,
            "audited_integer_domain": DOMAIN,
            "tolerance": TOL,
        },
        "damping_definition": "D(x)=(1+beta_tors*x)/(1+beta*x^eta)",
        "derivative_certificate": {
            "derivative_formula": "D'(x)=N(x)/(1+beta*x^eta)^2",
            "numerator_formula": "N(x)=beta_tors - eta*x^(eta-1) + beta_tors*(1-eta)*x^eta",
            "inequality_for_x_ge_1": "N(x) <= beta_tors - eta = -1.79 < 0 because x^(eta-1)>=1 and beta_tors*(1-eta)*x^eta<0",
            "denominator_positive": "(1+beta*x^eta)^2 > 0",
            "conclusion": "D'(x)<0 for every x in [1,11] and D(x)>0 there",
            "numerator_upper_bound_for_x_ge_1": derivative_upper_bound_numerator,
        },
        "grid_derivative_rows": grid,
        "integer_edge_drop_rows": discrete,
        "monotonicity_summary": {
            "D_1": endpoint_left,
            "D_11": endpoint_right,
            "D_1_over_D_11": endpoint_left / endpoint_right,
            "continuous_positive_certificate": True,
            "continuous_strictly_decreasing_certificate": True,
            "all_grid_derivative_numerators_negative": all(row["N_x_negative"] for row in grid),
            "all_grid_derivatives_negative": all(row["D_prime_x_negative"] for row in grid),
            "all_integer_edge_drops_positive": all(row["strict_drop_positive"] for row in discrete),
            "max_grid_derivative_numerator": max(row["derivative_numerator_N_x"] for row in grid),
            "min_integer_edge_drop": min(row["drop"] for row in discrete),
        },
        "blocker_context": {
            "necessity_status": necessity["status"],
            "transport_cocycle_status": cocycle["status"],
            "phase_zero_status": phase_zero["status"],
            "low_order_status": low_order["status"],
            "what_this_refines": "proves the damping/compression factor is a positive monotone envelope on [1,11], so sign flips must come from phase rather than damping",
            "strict_damping_derivation_status": "still_open; this certifies calculus consequences of the chosen damping formula, not its strict-side dynamical origin",
            "still_open": [
                "strict_damping_formula_derivation_from_nadsoliton_dynamics",
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_transport_derivation_from_nadsoliton_dynamics",
                "orientation_chi11_source",
                "chi11_uniqueness",
                "role_transfer_theorem",
            ],
        },
        "proof_certificate": {
            "positivity_step": "For x in [1,11], numerator 1+beta_tors*x and denominator 1+beta*x^eta are positive, so D(x)>0.",
            "derivative_step": "Differentiating D(x) gives D'(x)=N(x)/(1+beta*x^eta)^2 with the recorded N(x).",
            "negative_numerator_step": "For x>=1, N(x)<=beta_tors-eta=-1.79<0, so D'(x)<0.",
            "discrete_step": "The integer edge drops D(d)-D(d+1) are all positive for d=1..10, as expected from the continuous proof.",
            "phase_separation_step": "Because D(x)>0, damping/compression cannot create the phase sign flips certified by the zero-interlacing report.",
            "nonduplication": "This is a continuous damping-envelope calculus certificate, not another factor-subset, cocycle, low-order fit, or phase-zero report.",
            "theoretical_limit": "The certificate proves monotonicity of the chosen D(x); it does not derive D(x) from strict nadsoliton dynamics.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in solitonic state; this audit only checks calculus properties of the damping completion factor.",
            "forbidden_reading": "No separate informational layer below the nadsoliton is introduced.",
        },
        "hard_limits": [
            "K_strict_gate remains the current live/full operational kernel.",
            "No unqualified identity K_legacy_ont == K_strict_gate is claimed.",
            "No proof derives beta, eta, or D(x) from strict nadsoliton dynamics.",
            "No beta_tors -> chi_11 theorem is claimed.",
            "No legacy physical-role transfer to K_strict_gate is claimed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    summary = payload["monotonicity_summary"]
    cert = payload["derivative_certificate"]
    lines = [
        "# Strict completion damping continuous monotonicity certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This audit proves the damping/compression completion factor is positive and",
        "strictly decreasing on the continuous interval `[1,11]`.  It separates the",
        "monotone envelope role of damping from the phase sign flips certified by the",
        "phase-zero interlacing report.",
        "",
        "## Calculus certificate",
        "",
        f"- Damping definition: `{payload['damping_definition']}`",
        f"- Derivative formula: `{cert['derivative_formula']}`",
        f"- Numerator formula: `{cert['numerator_formula']}`",
        f"- Inequality: `{cert['inequality_for_x_ge_1']}`",
        f"- Conclusion: `{cert['conclusion']}`",
        "",
        "## Monotonicity summary",
        "",
        f"- `D(1)`: `{summary['D_1']:.12e}`",
        f"- `D(11)`: `{summary['D_11']:.12e}`",
        f"- `D(1)/D(11)`: `{summary['D_1_over_D_11']:.12e}`",
        f"- All grid derivative numerators negative: `{summary['all_grid_derivative_numerators_negative']}`",
        f"- All grid derivatives negative: `{summary['all_grid_derivatives_negative']}`",
        f"- All integer edge drops positive: `{summary['all_integer_edge_drops_positive']}`",
        f"- Max grid derivative numerator: `{summary['max_grid_derivative_numerator']:.12e}`",
        f"- Min integer edge drop: `{summary['min_integer_edge_drop']:.12e}`",
        "",
        "## Integer edge drops",
        "",
        "| edge | D_left | D_right | drop | ratio |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in payload["integer_edge_drop_rows"]:
        lines.append(
            "| {edge} | {left:.12e} | {right:.12e} | {drop:.12e} | {ratio:.12e} |".format(
                edge=row["edge"],
                left=row["D_left"],
                right=row["D_right"],
                drop=row["drop"],
                ratio=row["ratio_left_over_right"],
            )
        )
    lines.extend([
        "",
        "## Guarded interpretation",
        "",
        "This proves monotonicity of the selected damping formula only. It does not",
        "derive that formula from strict nadsoliton dynamics, does not prove",
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
