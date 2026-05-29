#!/usr/bin/env python3
"""Scratch probe: positive-lambda arrow-action selector certificate.

The previous arrow-action probe used a lexicographic pair (R,A), where
R=sum(e_i^2) is the finite ripple/Parseval ledger cost and
A=sum max(0,e_{i+1}-e_i) is a self-recorded arrow penalty for increases away
from the source endpoint.

This probe removes one possible ambiguity: does the selector require an
infinitesimal/lexicographic epsilon, or does an ordinary positive scalar
coupling already suffice?  In the finite ordered-ledger model (35 positive
5-tuples summing to 8), the target ledger (2,2,2,1,1) has R=14 and A=0.  Every
competitor has either larger R, or equal R and strictly positive A.  Therefore
for every real lambda>0, the scalar action S_lambda=R+lambda*A uniquely selects
(2,2,2,1,1).  At lambda=0 the result collapses back to the ten balanced
ripple minimizers.

This is a finite conditional action certificate.  It still does not derive the
arrow penalty, its sign/orientation, the support/source anchor, or the
underlying branch model from strict nadsoliton geometry.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PREVIOUS_ARROW_REPORT = HERE / "bridge_strict_alpha_self_recorded_arrow_action_lexicographic_selector_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_arrow_action_positive_lambda_certificate_report.json"
OUT_MD = HERE / "bridge_strict_alpha_arrow_action_positive_lambda_certificate_report.md"

SUPPORT_SIZE = 5
TARGET_BINARY_EXPONENT = 8
DENOMINATOR = 3
STRICT_TARGET_ETA = 9.0 / 5.0
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0
BALANCED_LEDGER = (2, 2, 2, 1, 1)
LAMBDA_SCAN = [
    Fraction(1, 100),
    Fraction(1, 10),
    Fraction(1, 2),
    Fraction(1, 1),
    Fraction(2, 1),
    Fraction(10, 1),
    Fraction(100, 1),
]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def eta_from_product(product: Fraction, branch_count: int) -> float:
    correction = float(product) ** (1.0 / branch_count)
    return math.log(NAD12_SUPPORT_SIZE * correction) / math.log(ALPHA_SCALE)


def positive_compositions(total: int, parts: int) -> list[tuple[int, ...]]:
    if parts == 1:
        return [(total,)] if total >= 1 else []
    rows: list[tuple[int, ...]] = []
    for head in range(1, total - parts + 2):
        for tail in positive_compositions(total - head, parts - 1):
            rows.append((head,) + tail)
    return rows


def ripple_sum_squares(ledger: tuple[int, ...]) -> int:
    return sum(value * value for value in ledger)


def parseval_non_dc_power_n5(ledger: tuple[int, ...]) -> int:
    total = sum(ledger)
    return SUPPORT_SIZE * ripple_sum_squares(ledger) - total * total


def arrow_increase_penalty(ledger: tuple[int, ...]) -> int:
    return sum(max(0, right - left) for left, right in zip(ledger, ledger[1:]))


def scalar_action(ledger: tuple[int, ...], coupling: Fraction) -> Fraction:
    return Fraction(ripple_sum_squares(ledger), 1) + coupling * arrow_increase_penalty(ledger)


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def ledger_row(ledger: tuple[int, ...], target_ripple: int, target_arrow: int) -> dict[str, Any]:
    ripple = ripple_sum_squares(ledger)
    arrow = arrow_increase_penalty(ledger)
    return {
        "ordered_ledger": list(ledger),
        "ripple_sum_squares_R": ripple,
        "arrow_increase_penalty_A": arrow,
        "parseval_non_dc_power_n5": parseval_non_dc_power_n5(ledger),
        "delta_R_vs_target": ripple - target_ripple,
        "delta_A_vs_target": arrow - target_arrow,
        "is_balanced_target_order": ledger == BALANCED_LEDGER,
    }


def scalar_scan(rows: list[dict[str, Any]], coupling: Fraction) -> dict[str, Any]:
    values: list[tuple[Fraction, list[int]]] = []
    for row in rows:
        value = Fraction(row["ripple_sum_squares_R"], 1) + coupling * row["arrow_increase_penalty_A"]
        values.append((value, row["ordered_ledger"]))
    minimum = min(value for value, _ledger in values)
    winners = [ledger for value, ledger in values if value == minimum]
    return {
        "lambda": fraction_text(coupling),
        "minimum_action": fraction_text(minimum),
        "winner_count": len(winners),
        "winners": winners,
        "uniquely_selects_balanced": winners == [list(BALANCED_LEDGER)],
    }


def symbolic_positive_lambda_certificate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    competitors = [row for row in rows if row["ordered_ledger"] != list(BALANCED_LEDGER)]
    equal_ripple = [row for row in competitors if row["delta_R_vs_target"] == 0]
    higher_ripple = [row for row in competitors if row["delta_R_vs_target"] > 0]
    lower_ripple = [row for row in competitors if row["delta_R_vs_target"] < 0]
    nonpositive_arrow_equal_ripple = [row for row in equal_ripple if row["delta_A_vs_target"] <= 0]
    zero_arrow_higher_ripple = [row for row in higher_ripple if row["delta_A_vs_target"] == 0]
    min_positive_equal_ripple_gap = min(row["delta_A_vs_target"] for row in equal_ripple)
    min_higher_ripple_gap = min(row["delta_R_vs_target"] for row in higher_ripple)
    return {
        "target_ledger": list(BALANCED_LEDGER),
        "target_R": ripple_sum_squares(BALANCED_LEDGER),
        "target_A": arrow_increase_penalty(BALANCED_LEDGER),
        "competitor_count": len(competitors),
        "equal_ripple_competitor_count": len(equal_ripple),
        "higher_ripple_competitor_count": len(higher_ripple),
        "lower_ripple_competitor_count": len(lower_ripple),
        "equal_ripple_delta_A_values": sorted({row["delta_A_vs_target"] for row in equal_ripple}),
        "higher_ripple_delta_R_values": sorted({row["delta_R_vs_target"] for row in higher_ripple}),
        "higher_ripple_zero_arrow_competitor_count": len(zero_arrow_higher_ripple),
        "min_positive_equal_ripple_delta_A": min_positive_equal_ripple_gap,
        "min_higher_ripple_delta_R": min_higher_ripple_gap,
        "all_competitors_satisfy_delta_R_nonnegative": all(row["delta_R_vs_target"] >= 0 for row in competitors),
        "all_equal_ripple_competitors_have_positive_delta_A": all(row["delta_A_vs_target"] > 0 for row in equal_ripple),
        "nonpositive_arrow_equal_ripple_competitors": [row["ordered_ledger"] for row in nonpositive_arrow_equal_ripple],
        "positive_lambda_gap_formula": "S_lambda(e)-S_lambda(target) = delta_R + lambda*delta_A; either delta_R>0 or delta_R=0 and delta_A>0, so the gap is >0 for every lambda>0.",
        "lambda_zero_obstruction": "At lambda=0 the action reduces to R and retains all ten balanced ripple minimizers, so the arrow coupling must be strictly positive to select the ordered self-record ledger.",
    }


def main() -> None:
    previous_arrow = load_json(PREVIOUS_ARROW_REPORT)
    product = Fraction(2 ** TARGET_BINARY_EXPONENT, DENOMINATOR ** SUPPORT_SIZE)
    ledgers = positive_compositions(TARGET_BINARY_EXPONENT, SUPPORT_SIZE)
    target_ripple = ripple_sum_squares(BALANCED_LEDGER)
    target_arrow = arrow_increase_penalty(BALANCED_LEDGER)
    rows = [ledger_row(ledger, target_ripple, target_arrow) for ledger in ledgers]
    lambda_scans = [scalar_scan(rows, coupling) for coupling in LAMBDA_SCAN]
    zero_scan = scalar_scan(rows, Fraction(0, 1))
    symbolic = symbolic_positive_lambda_certificate(rows)

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_ARROW_ACTION_POSITIVE_LAMBDA_CERTIFICATE_PROBE__NOT_A_THEOREM",
        "status": "candidate-supported-but-positive-lambda-arrow-action-not-derived",
        "target_identity_replay": {
            "q_power_product": f"{product.numerator}/{product.denominator}",
            "q_model": "q=(2^8/3^5)^(1/5)",
            "eta_from_q": eta_from_product(product, SUPPORT_SIZE),
            "eta_target": STRICT_TARGET_ETA,
            "eta_residual_vs_9_5": eta_from_product(product, SUPPORT_SIZE) - STRICT_TARGET_ETA,
            "branch_count_m": SUPPORT_SIZE,
            "binary_exponent_n": TARGET_BINARY_EXPONENT,
            "balanced_ledger": list(BALANCED_LEDGER),
        },
        "previous_arrow_action_replay": {
            "source_report": str(PREVIOUS_ARROW_REPORT.relative_to(ROOT)),
            "result_kind": previous_arrow["result_kind"],
            "lexicographic_winner": previous_arrow["finite_action_scan"]["lexicographic_ripple_then_arrow_minimizers"],
            "dominance_certificate": previous_arrow["finite_action_scan"]["dominance_certificate"],
        },
        "scalar_action_definition": {
            "R": "R(e)=sum_i e_i^2, equivalently the N=5 Parseval non-DC power is 5*R-(sum e_i)^2.",
            "A": "A(e)=sum_i max(0, e_{i+1}-e_i), a directed self-record arrow penalty for increases away from the source endpoint.",
            "S_lambda": "S_lambda(e)=R(e)+lambda*A(e).",
            "finite_claim": "Over positive ordered 5-ledgers with sum 8, S_lambda uniquely selects (2,2,2,1,1) for every real lambda>0.",
        },
        "finite_positive_lambda_certificate": {
            "total_positive_ordered_ledgers": len(rows),
            "target_row": next(row for row in rows if row["ordered_ledger"] == list(BALANCED_LEDGER)),
            "lambda_zero_scan": zero_scan,
            "positive_lambda_scan": lambda_scans,
            "all_scanned_positive_lambdas_uniquely_select_balanced": all(scan["uniquely_selects_balanced"] for scan in lambda_scans),
            "symbolic_certificate": symbolic,
            "rows": rows,
        },
        "candidate_interpretation": {
            "status": "candidate-supported-but-positive-lambda-arrow-action-not-derived",
            "honest_gain": "The previous lexicographic epsilon can be replaced, inside this finite model, by any ordinary strictly positive scalar coupling lambda.",
            "honest_limit": "The proof is conditional on the arrow term, its sign, source endpoint, ordered path, and branch model; it is not a derivation from strict nadsoliton geometry.",
            "ontology_guardrail": "The nadsoliton is treated as primordial information that may self-record formation/ordering; no separate informational substrate beneath it is introduced.",
            "selector_status": "conditional finite action certificate, not QW-2191 discharge",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "The nadsoliton is treated as information containing self-recorded formation/ordering information, not as a structure sitting above a deeper informational layer.",
            "No theorem derives the arrow-increase penalty A, its positive coupling lambda, its source endpoint, or its orientation from strict nadsoliton geometry.",
            "No theorem derives the d5 support, branch count m=5, total binary exponent n=8, denominator 3, or binary-rescale quotient from this scalar certificate.",
            "No Aut(Z_12)/N462/QW-2191 selector obstruction is discharged by adding positive lambda as a conditional action premise.",
            "No QW-2191 discharge and no ToE closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Try to derive the positive arrow coupling and its oriented source endpoint from a strict self-record/formation action, or prove a no-go showing that this sign and endpoint cannot be obtained without an extra selector premise.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha arrow-action positive-lambda certificate probe\n\n"
        "Status: finite positive-lambda scalar-action certificate; no strict selector discharge.\n\n"
        f"- Positive ordered ledgers: `{len(rows)}`.\n"
        f"- Lambda=0 winners: `{zero_scan['winner_count']}`.\n"
        f"- Scanned positive lambdas: `{[scan['lambda'] for scan in lambda_scans]}`.\n"
        f"- All scanned positive lambdas uniquely select balanced: `{all(scan['uniquely_selects_balanced'] for scan in lambda_scans)}`.\n"
        f"- Symbolic counts: equal-ripple competitors `{symbolic['equal_ripple_competitor_count']}`, higher-ripple competitors `{symbolic['higher_ripple_competitor_count']}`.\n"
        f"- Target replay: `q^5={product.numerator}/{product.denominator}`, eta residual `{eta_from_product(product, SUPPORT_SIZE)-STRICT_TARGET_ETA:.3e}`.\n"
        "- Honest read: any strictly positive arrow coupling selects the self-recorded balanced order, but the coupling/sign/source remain premises.\n"
        "- No false pass: no strict arrow-action theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
