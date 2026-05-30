#!/usr/bin/env python3
"""Scratch probe: exact rational calculus certificate for damping monotonicity.

The continuous damping monotonicity certificate differentiated

    D(x) = (1 + beta_tors*x)/(1 + beta*x**eta)

and checked the chosen strict-completion damping envelope on [1,11].  This probe
records the same derivative sign as an exact rational proof object for the
selected parameters beta=1, eta=9/5, beta_tors=1/100.

For x >= 1 the derivative numerator is

    N(x) = 1/100 - (9/5)*x**(4/5) - (1/125)*x**(9/5),

so N(x) <= 1/100 - 9/5 = -179/100 < 0.  The denominator is positive, hence
D'(x)<0 on [1,11].  This is an exact calculus certificate for the chosen
strict damping formula, not a derivation of beta/eta/D(x) from nadsoliton
dynamics.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_damping_exact_rational_calculus_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_damping_exact_rational_calculus_certificate_report.md"
DAMPING_FLOAT_REPORT = HERE / "bridge_strict_completion_damping_continuous_monotonicity_certificate_report.json"
CELL_SIGN_REPORT = HERE / "bridge_strict_completion_phase_zero_cell_sign_certificate_report.json"
Z2_REPORT = HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json"

BETA = Fraction(1, 1)
ETA = Fraction(9, 5)
BETA_TORS = Fraction(1, 100)
DOMAIN_LEFT = Fraction(1, 1)
DOMAIN_RIGHT = Fraction(11, 1)
INTEGER_DOMAIN = list(range(1, 12))
EXPECTED_FLIP_EDGES = ["1->2", "5->6", "7->8", "9->10"]


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


def exact_bound_rows() -> list[dict[str, Any]]:
    rows = []
    derivative_bound = BETA_TORS - ETA
    final_bound = derivative_bound + BETA_TORS * (1 - ETA)
    for d in INTEGER_DOMAIN:
        rows.append({
            "x": d,
            "x_ge_1": d >= 1,
            "x_power_4_5_lower_bound": "x^(4/5) >= 1",
            "negative_tail_term": "beta_tors*(1-eta)*x^(9/5) <= 0",
            "derivative_numerator_upper_bound_without_tail": fraction_payload(derivative_bound),
            "derivative_numerator_upper_bound_with_x_ge_1_tail": fraction_payload(final_bound),
            "strictly_negative_bound": derivative_bound < 0,
            "denominator_positive": True,
            "D_prime_certified_negative": derivative_bound < 0,
        })
    return rows


def edge_rows() -> list[dict[str, Any]]:
    rows = []
    for d in range(1, 11):
        rows.append({
            "edge": f"{d}->{d + 1}",
            "mean_value_interval": f"({d},{d + 1})",
            "D_prime_negative_on_open_interval_by_exact_bound": True,
            "therefore_D_left_greater_than_D_right": True,
            "damping_sign_flip_possible_on_edge": False,
        })
    return rows


def build_payload() -> dict[str, Any]:
    damping_float = load_json(DAMPING_FLOAT_REPORT)
    cell_sign = load_json(CELL_SIGN_REPORT)
    z2 = load_json(Z2_REPORT)
    derivative_upper_bound = BETA_TORS - ETA
    derivative_upper_bound_with_tail = derivative_upper_bound + BETA_TORS * (1 - ETA)
    rows = exact_bound_rows()
    edges = edge_rows()

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_DAMPING_EXACT_RATIONAL_CALCULUS_CERTIFICATE__NO_FLOAT_DERIVATIVE_SIGN",
        "status": "damping-derivative-negativity-certified-by-exact-rational-bound-on-continuous-interval",
        "source_reports": {
            "continuous_damping_monotonicity_certificate": str(DAMPING_FLOAT_REPORT.relative_to(ROOT)),
            "cell_sign_certificate": str(CELL_SIGN_REPORT.relative_to(ROOT)),
            "phase_sign_z2_coboundary_certificate": str(Z2_REPORT.relative_to(ROOT)),
        },
        "rational_parameters": {
            "beta": fraction_payload(BETA),
            "eta": fraction_payload(ETA),
            "beta_tors": fraction_payload(BETA_TORS),
            "certified_interval": [fraction_payload(DOMAIN_LEFT), fraction_payload(DOMAIN_RIGHT)],
        },
        "exact_derivative_certificate": {
            "damping_formula": "D(x)=(1+beta_tors*x)/(1+beta*x^(eta))",
            "derivative_formula": "D'(x)=N(x)/(1+x^(9/5))^2",
            "numerator_formula": "N(x)=1/100-(9/5)*x^(4/5)-(1/125)*x^(9/5)",
            "positive_denominator_step": "(1+x^(9/5))^2 > 0 for x in [1,11]",
            "x_ge_1_step": "For x>=1, x^(4/5)>=1 and x^(9/5)>=0.",
            "negative_tail_step": "-(1/125)*x^(9/5)<=0.",
            "rational_upper_bound": fraction_payload(derivative_upper_bound),
            "stronger_rational_upper_bound_using_tail": fraction_payload(derivative_upper_bound_with_tail),
            "upper_bound_strictly_negative": derivative_upper_bound < 0,
            "conclusion": "D'(x)<0 on [1,11] by exact rational comparison, with no floating derivative sign decision.",
        },
        "grid_exact_bound_rows": rows,
        "integer_edge_mean_value_rows": edges,
        "exact_rational_damping_summary": {
            "continuous_positive_certificate": True,
            "continuous_strictly_decreasing_certificate": True,
            "derivative_upper_bound_text": f"{derivative_upper_bound.numerator}/{derivative_upper_bound.denominator}",
            "derivative_upper_bound_decimal": float(derivative_upper_bound),
            "all_grid_bounds_negative": all(row["strictly_negative_bound"] for row in rows),
            "all_integer_edges_drop_by_mean_value_theorem": all(row["therefore_D_left_greater_than_D_right"] for row in edges),
            "damping_can_supply_phase_sign_flips": False,
            "phase_flip_edges_remain_phase_only": cell_sign["cell_sign_summary"]["derived_phase_sign_flip_edges"] == EXPECTED_FLIP_EDGES,
            "z2_phase_flip_edges_remain_unchanged": z2["z2_coboundary_summary"]["derived_phase_sign_flip_edges"] == EXPECTED_FLIP_EDGES,
            "matches_previous_float_monotonicity_certificate": (
                damping_float["monotonicity_summary"]["continuous_positive_certificate"]
                and damping_float["monotonicity_summary"]["continuous_strictly_decreasing_certificate"]
            ),
        },
        "blocker_context": {
            "what_this_refines": "replaces floating derivative-sign reliance with an exact rational upper-bound proof for the chosen damping envelope",
            "continuous_damping_status": damping_float["status"],
            "cell_sign_status": cell_sign["status"],
            "phase_sign_z2_status": z2["status"],
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
            "parameter_step": "The selected damping parameters are represented exactly as beta=1, eta=9/5, beta_tors=1/100.",
            "derivative_step": "Symbolic differentiation gives D'(x)=N(x)/(1+x^(9/5))^2 with the recorded N(x).",
            "rational_bound_step": "For x>=1, N(x)<=1/100-9/5=-179/100<0; the remaining tail is also non-positive.",
            "mean_value_step": "Since D'(x)<0 on every open integer edge, D(d)>D(d+1) follows by the mean value theorem for d=1..10.",
            "phase_separation_step": "Because D(x)>0 and is decreasing, damping remains an envelope and cannot supply the certified phase sign flips.",
            "nonduplication": "This is an exact rational calculus certificate, not another phase-zero, cell-sign, Z2, real-valued cocycle, or low-order no-go audit.",
            "theoretical_limit": "The certificate proves consequences of the selected D(x); it does not derive beta, eta, beta_tors, or D(x) from strict nadsoliton dynamics.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in solitonic state; this audit only checks exact calculus consequences of the selected damping factor.",
            "forbidden_reading": "No separate informational layer below the nadsoliton is introduced.",
        },
        "hard_limits": [
            "K_strict_gate remains the current live/full operational kernel.",
            "No unqualified identity K_legacy_ont == K_strict_gate is claimed.",
            "No proof derives beta, eta, beta_tors, or D(x) from strict nadsoliton dynamics.",
            "No beta_tors -> chi_11 theorem is claimed.",
            "No legacy physical-role transfer to K_strict_gate is claimed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["exact_derivative_certificate"]
    summary = payload["exact_rational_damping_summary"]
    lines = [
        "# Strict completion damping exact rational calculus certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This audit proves the damping derivative sign by exact rational bounds",
        "for `beta=1`, `eta=9/5`, and `beta_tors=1/100`.  It strengthens the",
        "continuous damping report by avoiding floating derivative-sign decisions.",
        "",
        "## Exact derivative certificate",
        "",
        f"- Formula: `{cert['damping_formula']}`",
        f"- Derivative: `{cert['derivative_formula']}`",
        f"- Numerator: `{cert['numerator_formula']}`",
        f"- Rational upper bound: `{cert['rational_upper_bound']['text']}` = `{cert['rational_upper_bound']['decimal']:.12e}`",
        f"- Upper bound strictly negative: `{cert['upper_bound_strictly_negative']}`",
        f"- Conclusion: {cert['conclusion']}",
        "",
        "## Summary",
        "",
        f"- Continuous positive certificate: `{summary['continuous_positive_certificate']}`",
        f"- Continuous strictly decreasing certificate: `{summary['continuous_strictly_decreasing_certificate']}`",
        f"- All grid bounds negative: `{summary['all_grid_bounds_negative']}`",
        f"- All integer edges drop by MVT: `{summary['all_integer_edges_drop_by_mean_value_theorem']}`",
        f"- Damping can supply phase sign flips: `{summary['damping_can_supply_phase_sign_flips']}`",
        f"- Matches previous float monotonicity certificate: `{summary['matches_previous_float_monotonicity_certificate']}`",
        "",
        "## Proof certificate",
        "",
    ]
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
