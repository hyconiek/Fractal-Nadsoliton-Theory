#!/usr/bin/env python3
"""Scratch probe: rational robustness margin for strict phase-zero placement.

The rational-interval phase-zero certificate proved, without float zero placement,
that the selected strict phase zero k=0 lies in the integer edge (7,8).  This
probe takes the next proof-style step: compute a rational safety margin showing
how far the strict phase parameters can move before the certified edge placement
can fail.

Using exact strict parameters omega=743/4000 and phi=13/80 plus
333/106 < pi < 355/113, the edge condition for the k=0 zero is

    7*omega + phi < pi/2 < 8*omega + phi.

The certified margins are:

    m_left  = pi_lower/2 - (7*omega+phi)  > 0,
    m_right = (8*omega+phi) - pi_upper/2 > 0.

If |Delta omega|<=eps and |Delta phi|<=eps with

    eps <= min(m_left/8, m_right/9),

then both inequalities remain true for all parameters in that rectangle.  The
probe also checks a k=1 exclusion margin so the next strict zero remains right
of d=11.  This is a robustness certificate for the selected phase parameters,
not a derivation of them from strict nadsoliton dynamics.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_zero_margin_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_zero_margin_certificate_report.md"
RATIONAL_ZERO_REPORT = HERE / "bridge_strict_completion_phase_zero_rational_interval_certificate_report.json"
PHASE_ZERO_REPORT = HERE / "bridge_strict_completion_phase_zero_interlacing_certificate_report.json"

PI_LOWER = Fraction(333, 106)
PI_UPPER = Fraction(355, 113)
STRICT_OMEGA = Fraction(743, 4000)
STRICT_PHI = Fraction(13, 80)
DOMAIN = list(range(12))


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


def compute_margins() -> dict[str, Fraction]:
    left_margin = PI_LOWER / 2 - (7 * STRICT_OMEGA + STRICT_PHI)
    right_margin = (8 * STRICT_OMEGA + STRICT_PHI) - PI_UPPER / 2
    next_zero_right_margin = 3 * PI_LOWER / 2 - (11 * STRICT_OMEGA + STRICT_PHI)
    previous_zero_left_margin = STRICT_PHI + PI_LOWER / 2
    limiting_eps = min(left_margin / 8, right_margin / 9, next_zero_right_margin / 12, previous_zero_left_margin / 2)
    symmetric_eps = limiting_eps / 2
    return {
        "left_margin": left_margin,
        "right_margin": right_margin,
        "next_zero_right_margin": next_zero_right_margin,
        "previous_zero_left_margin": previous_zero_left_margin,
        "limiting_symmetric_parameter_epsilon": limiting_eps,
        "symmetric_parameter_epsilon": symmetric_eps,
    }


def worst_case_rows(epsilon: Fraction) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    # Worst cases for the four certified inequalities under |dw|,|dp| <= epsilon.
    cases = [
        {
            "name": "k0_left_edge_upper_worst_case",
            "expression": "7*(omega+eps)+(phi+eps) < pi_lower/2",
            "left_value": 7 * (STRICT_OMEGA + epsilon) + (STRICT_PHI + epsilon),
            "right_value": PI_LOWER / 2,
            "expected_relation": "<",
        },
        {
            "name": "k0_right_edge_lower_worst_case",
            "expression": "pi_upper/2 < 8*(omega-eps)+(phi-eps)",
            "left_value": PI_UPPER / 2,
            "right_value": 8 * (STRICT_OMEGA - epsilon) + (STRICT_PHI - epsilon),
            "expected_relation": "<",
        },
        {
            "name": "k1_stays_right_of_d11_worst_case",
            "expression": "11*(omega+eps)+(phi+eps) < 3*pi_lower/2",
            "left_value": 11 * (STRICT_OMEGA + epsilon) + (STRICT_PHI + epsilon),
            "right_value": 3 * PI_LOWER / 2,
            "expected_relation": "<",
        },
        {
            "name": "k_minus_1_stays_left_of_d0_worst_case",
            "expression": "-pi_lower/2 < phi-eps",
            "left_value": -PI_LOWER / 2,
            "right_value": STRICT_PHI - epsilon,
            "expected_relation": "<",
        },
    ]
    for case in cases:
        left = case["left_value"]
        right = case["right_value"]
        rows.append({
            "name": case["name"],
            "expression": case["expression"],
            "left_value": fraction_payload(left),
            "right_value": fraction_payload(right),
            "strict_inequality_holds": left < right,
            "slack": fraction_payload(right - left),
        })
    return rows


def build_payload() -> dict[str, Any]:
    rational_zero = load_json(RATIONAL_ZERO_REPORT)
    phase_zero = load_json(PHASE_ZERO_REPORT)
    margins = compute_margins()
    epsilon = margins["symmetric_parameter_epsilon"]
    worst_cases = worst_case_rows(epsilon)

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_ZERO_MARGIN_CERTIFICATE__RATIONAL_PARAMETER_ROBUSTNESS",
        "status": "strict-phase-zero-edge-placement-has-positive-rational-parameter-margin",
        "source_reports": {
            "rational_zero_certificate": str(RATIONAL_ZERO_REPORT.relative_to(ROOT)),
            "phase_zero_interlacing_certificate": str(PHASE_ZERO_REPORT.relative_to(ROOT)),
        },
        "rational_inputs": {
            "pi_lower_bound": fraction_payload(PI_LOWER),
            "pi_upper_bound": fraction_payload(PI_UPPER),
            "strict_omega": fraction_payload(STRICT_OMEGA),
            "strict_phi": fraction_payload(STRICT_PHI),
        },
        "margin_definitions": {
            "k0_left_margin": "pi_lower/2 - (7*omega + phi)",
            "k0_right_margin": "(8*omega + phi) - pi_upper/2",
            "k1_right_exclusion_margin": "3*pi_lower/2 - (11*omega + phi)",
            "k_minus_1_left_exclusion_margin": "phi + pi_lower/2",
            "limiting_symmetric_epsilon_rule": "eps_limit = min(left_margin/8, right_margin/9, k1_margin/12, k_minus_1_margin/2)",
            "certified_symmetric_epsilon_rule": "eps = eps_limit/2 to keep all inequalities strict",
        },
        "margin_values": {name: fraction_payload(value) for name, value in margins.items()},
        "worst_case_inequality_rows": worst_cases,
        "robustness_summary": {
            "all_margins_positive": all(value > 0 for value in margins.values()),
            "all_worst_case_inequalities_hold_at_epsilon": all(row["strict_inequality_holds"] for row in worst_cases),
            "symmetric_parameter_epsilon_decimal": float(epsilon),
            "active_margin_source": "k0_right_margin/9",
            "certified_epsilon_is_half_of_limiting_epsilon": True,
            "certified_phase_sign_flip_edges_preserved": rational_zero["interval_summary"]["phase_sign_flip_edges_from_rational_intervals"],
            "certified_phase_sign_pattern_preserved": rational_zero["interval_summary"]["phase_transport_sign_pattern_from_rational_intervals"],
            "matches_phase_zero_report": rational_zero["interval_summary"]["matches_float_phase_zero_report_flip_edges"] and rational_zero["interval_summary"]["matches_float_phase_zero_report_sign_pattern"],
        },
        "blocker_context": {
            "rational_zero_status": rational_zero["status"],
            "phase_zero_status": phase_zero["status"],
            "what_this_refines": "adds a rational perturbation margin around strict omega/phi for the certified phase-zero edge placement",
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_transport_derivation_from_nadsoliton_dynamics",
                "orientation_chi11_source",
                "chi11_uniqueness",
                "role_transfer_theorem",
            ],
        },
        "proof_certificate": {
            "edge_condition_step": "The k=0 strict zero stays in (7,8) if 7*omega+phi < pi/2 < 8*omega+phi.",
            "margin_step": "Using pi_lower for the left inequality and pi_upper for the right inequality gives positive rational margins.",
            "epsilon_step": "The symmetric epsilon bound subtracts worst-case 7*eps+eps and 8*eps+eps perturbations from those margins.",
            "exclusion_step": "Additional k=-1 and k=1 inequalities keep the adjacent strict zeros outside [0,11].",
            "nonduplication": "This is a rational robustness-margin certificate, not another zero-location, damping, cocycle, or low-order fit audit.",
            "theoretical_limit": "The certificate proves local parameter robustness of the chosen phase placement; it does not derive omega/phi from strict dynamics.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in solitonic state; this audit only checks robustness of phase-zero placement.",
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
    values = payload["margin_values"]
    summary = payload["robustness_summary"]
    lines = [
        "# Strict completion phase-zero margin certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This audit computes a rational robustness margin around the selected strict",
        "phase parameters.  Within the certified symmetric parameter rectangle, the",
        "strict k=0 zero remains in edge `7->8`, adjacent strict zeros remain outside",
        "`[0,11]`, and the phase sign-flip pattern is preserved.",
        "",
        "## Margin values",
        "",
        f"- k0 left margin: `{values['left_margin']['text']}` = `{values['left_margin']['decimal']:.12e}`",
        f"- k0 right margin: `{values['right_margin']['text']}` = `{values['right_margin']['decimal']:.12e}`",
        f"- k1 right-exclusion margin: `{values['next_zero_right_margin']['text']}` = `{values['next_zero_right_margin']['decimal']:.12e}`",
        f"- k-1 left-exclusion margin: `{values['previous_zero_left_margin']['text']}` = `{values['previous_zero_left_margin']['decimal']:.12e}`",
        f"- Limiting symmetric epsilon: `{values['limiting_symmetric_parameter_epsilon']['text']}` = `{values['limiting_symmetric_parameter_epsilon']['decimal']:.12e}`",
        f"- Certified symmetric epsilon: `{values['symmetric_parameter_epsilon']['text']}` = `{values['symmetric_parameter_epsilon']['decimal']:.12e}`",
        "",
        "## Robustness summary",
        "",
        f"- All margins positive: `{summary['all_margins_positive']}`",
        f"- All worst-case inequalities hold at epsilon: `{summary['all_worst_case_inequalities_hold_at_epsilon']}`",
        f"- Active margin source: `{summary['active_margin_source']}`",
        f"- Preserved flip edges: `{summary['certified_phase_sign_flip_edges_preserved']}`",
        f"- Preserved sign pattern: `{summary['certified_phase_sign_pattern_preserved']}`",
        "",
        "## Worst-case inequalities",
        "",
        "| name | expression | slack | holds |",
        "|---|---|---:|---:|",
    ]
    for row in payload["worst_case_inequality_rows"]:
        lines.append(
            f"| {row['name']} | `{row['expression']}` | {row['slack']['text']} | {row['strict_inequality_holds']} |"
        )
    lines.extend([
        "",
        "## Guarded interpretation",
        "",
        "This is only a local robustness certificate for the selected strict phase",
        "parameters. It does not derive strict omega/phi from nadsoliton dynamics,",
        "does not prove `beta_tors -> chi_11`, and does not transfer legacy physical",
        "roles onto `K_strict_gate`.",
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
