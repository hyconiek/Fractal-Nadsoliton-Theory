#!/usr/bin/env python3
"""Scratch probe: self-recorded arrow-action lexicographic ledger selector.

The monotone min-ripple stack used two premises in sequence: monotone formation
and then min-ripple.  This probe asks whether those can be packaged into one
finite action certificate without first filtering the 35 ordered ledgers.

Finite answer: yes, conditionally.  Define a ripple term R=sum(e_i^2) and a
self-recorded arrow penalty A=sum max(0,e_{i+1}-e_i), which punishes increases
away from the source endpoint.  Across all 35 positive ordered ledgers with
sum e_i=8, the lexicographic pair (R,A) has the unique minimizer
(2,2,2,1,1).  Equivalently, R + epsilon*A has the same unique minimizer for any
positive epsilon when interpreted lexicographically/symbolically.  This is a
finite action certificate, but it is not a strict derivation of the action from
nadsoliton geometry.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MONOTONE_STACK = HERE / "bridge_strict_alpha_self_recorded_monotone_min_ripple_stack_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_self_recorded_arrow_action_lexicographic_selector_report.json"
OUT_MD = HERE / "bridge_strict_alpha_self_recorded_arrow_action_lexicographic_selector_report.md"

SUPPORT_SIZE = 5
TARGET_BINARY_EXPONENT = 8
DENOMINATOR = 3
STRICT_TARGET_ETA = 9.0 / 5.0
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0
BALANCED_LEDGER = (2, 2, 2, 1, 1)


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


def endpoint_self_recorded_from_left(ledger: tuple[int, ...]) -> bool:
    return ledger[0] > ledger[-1]


def ledger_row(ledger: tuple[int, ...]) -> dict[str, Any]:
    ripple = ripple_sum_squares(ledger)
    arrow = arrow_increase_penalty(ledger)
    return {
        "ordered_ledger": list(ledger),
        "ripple_sum_squares": ripple,
        "parseval_non_dc_power_n5": parseval_non_dc_power_n5(ledger),
        "self_recorded_arrow_increase_penalty": arrow,
        "lexicographic_pair_ripple_then_arrow": [ripple, arrow],
        "endpoint_self_recorded_from_left": endpoint_self_recorded_from_left(ledger),
        "is_balanced_target_order": ledger == BALANCED_LEDGER,
    }


def winner_packet(rows: list[dict[str, Any]], key: str) -> dict[str, Any]:
    minimum = min(row[key] for row in rows)
    winners = [row["ordered_ledger"] for row in rows if row[key] == minimum]
    return {"minimum": minimum, "winner_count": len(winners), "winners": winners}


def lexicographic_winner_packet(rows: list[dict[str, Any]]) -> dict[str, Any]:
    minimum = min(tuple(row["lexicographic_pair_ripple_then_arrow"]) for row in rows)
    winners = [row["ordered_ledger"] for row in rows if tuple(row["lexicographic_pair_ripple_then_arrow"]) == minimum]
    return {"minimum_pair": list(minimum), "winner_count": len(winners), "winners": winners}


def dominance_certificate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    target = next(row for row in rows if row["ordered_ledger"] == list(BALANCED_LEDGER))
    equal_ripple_competitors = [
        row for row in rows
        if row["ordered_ledger"] != list(BALANCED_LEDGER)
        and row["ripple_sum_squares"] == target["ripple_sum_squares"]
    ]
    higher_ripple_competitors = [
        row for row in rows
        if row["ordered_ledger"] != list(BALANCED_LEDGER)
        and row["ripple_sum_squares"] > target["ripple_sum_squares"]
    ]
    return {
        "target_ledger": target["ordered_ledger"],
        "target_ripple": target["ripple_sum_squares"],
        "target_arrow_penalty": target["self_recorded_arrow_increase_penalty"],
        "equal_ripple_competitor_count": len(equal_ripple_competitors),
        "equal_ripple_competitor_arrow_penalties": sorted({row["self_recorded_arrow_increase_penalty"] for row in equal_ripple_competitors}),
        "all_equal_ripple_competitors_have_positive_arrow_penalty": all(row["self_recorded_arrow_increase_penalty"] > 0 for row in equal_ripple_competitors),
        "higher_ripple_competitor_count": len(higher_ripple_competitors),
        "all_other_ledgers_have_higher_ripple_or_positive_arrow_tiebreak": all(
            row["ripple_sum_squares"] > target["ripple_sum_squares"]
            or row["self_recorded_arrow_increase_penalty"] > target["self_recorded_arrow_increase_penalty"]
            for row in rows
            if row["ordered_ledger"] != list(BALANCED_LEDGER)
        ),
        "epsilon_statement": "For every epsilon>0, R+epsilon*A is uniquely minimized by (2,2,2,1,1) in the symbolic/lexicographic sense: equal-R competitors have A>0, higher-R competitors already lose in R.",
    }


def main() -> None:
    monotone = load_json(MONOTONE_STACK)
    product = Fraction(2 ** TARGET_BINARY_EXPONENT, DENOMINATOR ** SUPPORT_SIZE)
    rows = [ledger_row(ledger) for ledger in positive_compositions(TARGET_BINARY_EXPONENT, SUPPORT_SIZE)]
    ripple_winners = winner_packet(rows, "ripple_sum_squares")
    arrow_winners = winner_packet(rows, "self_recorded_arrow_increase_penalty")
    lexicographic = lexicographic_winner_packet(rows)
    dominance = dominance_certificate(rows)

    report = {
        "status": "OPEN_STRICT_ALPHA_SELF_RECORDED_ARROW_ACTION_LEXICOGRAPHIC_SELECTOR_PREMISE_BASED_NO_STRICT_DISCHARGE",
        "result_kind": "SCRATCH_STRICT_ALPHA_SELF_RECORDED_ARROW_ACTION_LEXICOGRAPHIC_SELECTOR_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "self_recorded_monotone_min_ripple_stack": str(MONOTONE_STACK.relative_to(ROOT)),
        },
        "previous_monotone_stack_replay": {
            "result_kind": monotone["result_kind"],
            "self_recorded_monotone_count": monotone["finite_stack_scan"]["self_recorded_monotone_count"],
            "all_checked_power_costs_uniquely_select_balanced": monotone["finite_stack_scan"]["all_checked_power_costs_uniquely_select_balanced"],
            "candidate_status": monotone["candidate_interpretation"]["status"],
        },
        "target_identity_replay": {
            "q_power_product": f"{product.numerator}/{product.denominator}",
            "eta_residual_vs_9_5": eta_from_product(product, SUPPORT_SIZE) - STRICT_TARGET_ETA,
            "balanced_ledger": list(BALANCED_LEDGER),
        },
        "action_definition": {
            "ripple_term_R": "R(e)=sum_i e_i^2, equivalent to Parseval non-DC ripple at fixed sum and m=5.",
            "arrow_term_A": "A(e)=sum_i max(0, e_{i+1}-e_i), penalizing increases away from the self-recorded source endpoint.",
            "lexicographic_action": "First minimize R, then minimize A among equal-R ledgers.",
            "epsilon_action": "R + epsilon*A has the same winner for every symbolic epsilon>0 under the recorded dominance split.",
        },
        "finite_action_scan": {
            "total_positive_ordered_ledgers": len(rows),
            "rows": rows,
            "ripple_only_minimizers": ripple_winners,
            "arrow_only_minimizers": arrow_winners,
            "lexicographic_ripple_then_arrow_minimizers": lexicographic,
            "parseval_minimizers": winner_packet(rows, "parseval_non_dc_power_n5"),
            "dominance_certificate": dominance,
            "lexicographic_uniquely_selects_balanced": lexicographic["winners"] == [list(BALANCED_LEDGER)],
        },
        "selector_consequence": {
            "what_is_gained": "The monotone/min-ripple stack can be packaged as a single finite lexicographic arrow action over all 35 ordered ledgers.",
            "what_each_term_does": "Ripple selects the balanced multiset but leaves 10 orderings; arrow penalty breaks that ordering degeneracy in favor of the self-recorded monotone order.",
            "what_remains_unproved": "The strict theory still has to derive R and A as genuine nadsoliton formation/action terms rather than assume them.",
        },
        "candidate_interpretation": {
            "supported_by_this_probe": True,
            "content": "The lexicographic pair (sum squares, arrow-increase penalty) uniquely selects the ordered balanced ledger from all 35 positive ledgers.",
            "why_this_is_more_proof_like": "The probe provides a dominance certificate: equal-ripple competitors have positive arrow penalty, while all remaining competitors have higher ripple.",
            "why_this_is_not_enough": "The arrow action and ripple action are still conditional model-internal selector terms, not derived strict geometry.",
            "status": "candidate-supported-but-arrow-action-not-derived",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "The nadsoliton is treated as information carrying finite self-record patterns, but no separate sub-nadsoliton informational layer is introduced.",
            "The arrow-increase penalty is a finite selector-action candidate, not a strict-core theorem.",
            "The ripple term is used as a model-internal selector action here and is not derived from strict nadsoliton geometry.",
            "No theorem derives the d5 support, arrow action, balanced ledger, endpoint-source convention, or strict source-localizing term from strict nadsoliton geometry.",
            "No theorem derives branch count m=5, total exponent n=8, denominator 3, or binary-rescale quotient.",
            "No Aut(Z_12)/N462/QW-2191 selector obstruction is discharged by this arrow-action certificate.",
            "No theorem derives eta=9/5 without the branch model, denominator/quotient convention, support premise, self-record arrow action, and ripple ledger selector premise.",
            "No QW-2191 discharge and no ToE closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Try to derive the arrow-increase penalty and ripple term as strict internal action components of nadsoliton self-record formation, or prove a no-go for such an action derivation.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha self-recorded arrow-action lexicographic selector probe\n\n"
        "Status: finite lexicographic arrow-action certificate; no strict selector discharge.\n\n"
        f"- Positive ordered ledgers: `{len(rows)}`.\n"
        f"- Ripple-only minimizers: `{ripple_winners}`.\n"
        f"- Lexicographic (R,A) minimizer: `{lexicographic}`.\n"
        f"- Dominance: equal-ripple competitors `{dominance['equal_ripple_competitor_count']}`, higher-ripple competitors `{dominance['higher_ripple_competitor_count']}`.\n"
        f"- Target replay: `q^5={product.numerator}/{product.denominator}`, eta residual `{eta_from_product(product, SUPPORT_SIZE)-STRICT_TARGET_ETA:.3e}`.\n"
        "- Honest read: arrow-action + ripple selects balanced order, but the action terms remain premises.\n"
        "- No false pass: no strict arrow-action theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
