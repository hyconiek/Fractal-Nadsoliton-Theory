#!/usr/bin/env python3
"""Scratch probe: self-recorded monotone + min-ripple ledger selector stack.

The endpoint self-record non-uniqueness audit showed that endpoint anchoring is a
real decoder but not a ledger selector.  This probe tests a stronger, still
conditional, formation-style stack: if the self-recorded source is the high
endpoint and the ledger is monotone non-increasing away from that source, which
ledgers remain?  Then, if a symmetric convex/min-ripple cost is added, does the
balanced ledger become the unique finite minimizer?

Finite answer: self-record + monotone formation reduces all 35 positive ordered
ledgers of total exponent 8 to the three canonical ledgers
(4,1,1,1,1), (3,2,1,1,1), and (2,2,2,1,1).  It still does not select one.
Adding any checked strictly convex power cost sum(e_i^p), p=2..6, uniquely
selects (2,2,2,1,1).  This is a clean conditional selector stack, not a strict
nadsoliton theorem.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
NONUNIQUENESS = HERE / "bridge_strict_alpha_self_recorded_anchor_ledger_nonuniqueness_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_self_recorded_monotone_min_ripple_stack_report.json"
OUT_MD = HERE / "bridge_strict_alpha_self_recorded_monotone_min_ripple_stack_report.md"

SUPPORT_SIZE = 5
TARGET_BINARY_EXPONENT = 8
DENOMINATOR = 3
STRICT_TARGET_ETA = 9.0 / 5.0
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0
BALANCED_LEDGER = (2, 2, 2, 1, 1)
CHECKED_POWERS = tuple(range(2, 7))


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


def endpoint_self_recorded_from_left(ledger: tuple[int, ...]) -> bool:
    return ledger[0] > ledger[-1]


def monotone_nonincreasing(ledger: tuple[int, ...]) -> bool:
    return all(left >= right for left, right in zip(ledger, ledger[1:]))


def power_cost(ledger: tuple[int, ...], power: int) -> int:
    return sum(value ** power for value in ledger)


def parseval_non_dc_power_n5(ledger: tuple[int, ...]) -> int:
    total = sum(ledger)
    return SUPPORT_SIZE * sum(value * value for value in ledger) - total * total


def ledger_row(ledger: tuple[int, ...]) -> dict[str, Any]:
    endpoint = endpoint_self_recorded_from_left(ledger)
    monotone = monotone_nonincreasing(ledger)
    return {
        "ordered_ledger": list(ledger),
        "endpoint_self_recorded_from_left": endpoint,
        "monotone_nonincreasing_from_source": monotone,
        "passes_self_recorded_monotone_stack": endpoint and monotone,
        "sum_squares": power_cost(ledger, 2),
        "parseval_non_dc_power_n5": parseval_non_dc_power_n5(ledger),
        "power_costs": {str(power): power_cost(ledger, power) for power in CHECKED_POWERS},
    }


def minimizers(rows: list[dict[str, Any]], key: str) -> dict[str, Any]:
    minimum = min(row[key] for row in rows)
    winners = [row["ordered_ledger"] for row in rows if row[key] == minimum]
    return {"minimum": minimum, "winner_count": len(winners), "winners": winners}


def power_minimizers(rows: list[dict[str, Any]]) -> dict[str, Any]:
    result = {}
    for power in CHECKED_POWERS:
        values = [(row["ordered_ledger"], row["power_costs"][str(power)]) for row in rows]
        minimum = min(value for _, value in values)
        result[str(power)] = {
            "minimum": minimum,
            "winner_count": sum(1 for _, value in values if value == minimum),
            "winners": [ledger for ledger, value in values if value == minimum],
        }
    return result


def smoothing_steps_from_extreme_to_balanced() -> list[dict[str, Any]]:
    path = [(4, 1, 1, 1, 1), (3, 2, 1, 1, 1), (2, 2, 2, 1, 1)]
    steps = []
    for before, after in zip(path, path[1:]):
        steps.append(
            {
                "before": list(before),
                "after": list(after),
                "sum_squares_before": power_cost(before, 2),
                "sum_squares_after": power_cost(after, 2),
                "sum_squares_drop": power_cost(before, 2) - power_cost(after, 2),
                "parseval_non_dc_before": parseval_non_dc_power_n5(before),
                "parseval_non_dc_after": parseval_non_dc_power_n5(after),
                "parseval_non_dc_drop": parseval_non_dc_power_n5(before) - parseval_non_dc_power_n5(after),
            }
        )
    return steps


def main() -> None:
    nonunique = load_json(NONUNIQUENESS)
    product = Fraction(2 ** TARGET_BINARY_EXPONENT, DENOMINATOR ** SUPPORT_SIZE)
    rows = [ledger_row(ledger) for ledger in positive_compositions(TARGET_BINARY_EXPONENT, SUPPORT_SIZE)]
    endpoint_rows = [row for row in rows if row["endpoint_self_recorded_from_left"]]
    stack_rows = [row for row in rows if row["passes_self_recorded_monotone_stack"]]
    power_winners = power_minimizers(stack_rows)

    report = {
        "status": "OPEN_STRICT_ALPHA_SELF_RECORDED_MONOTONE_MIN_RIPPLE_STACK_PREMISE_BASED_NO_STRICT_DISCHARGE",
        "result_kind": "SCRATCH_STRICT_ALPHA_SELF_RECORDED_MONOTONE_MIN_RIPPLE_STACK_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "self_recorded_anchor_ledger_nonuniqueness": str(NONUNIQUENESS.relative_to(ROOT)),
        },
        "previous_nonuniqueness_replay": {
            "result_kind": nonunique["result_kind"],
            "total_positive_ordered_ledgers": nonunique["ordered_ledger_scan"]["total_positive_ordered_ledgers"],
            "endpoint_distinct_self_recording_count": nonunique["ordered_ledger_scan"]["endpoint_distinct_self_recording_count"],
            "endpoint_equal_ambiguous_count": nonunique["ordered_ledger_scan"]["endpoint_equal_ambiguous_count"],
            "candidate_status": nonunique["candidate_interpretation"]["status"],
        },
        "target_identity_replay": {
            "q_power_product": f"{product.numerator}/{product.denominator}",
            "eta_residual_vs_9_5": eta_from_product(product, SUPPORT_SIZE) - STRICT_TARGET_ETA,
            "balanced_ledger": list(BALANCED_LEDGER),
        },
        "selector_stack_definition": {
            "stage_1_self_record": "The ordered ledger is read from the higher source endpoint, so e_0 > e_4.",
            "stage_2_monotone_formation": "The self-recorded formation profile is monotone non-increasing away from the source endpoint.",
            "stage_3_min_ripple": "Among the remaining monotone self-recorded ledgers, minimize a symmetric convex power cost / Parseval ripple.",
            "premise_warning": "The monotone formation law and min-ripple action are selector premises here, not strict exports.",
        },
        "finite_stack_scan": {
            "all_ordered_rows": rows,
            "total_positive_ordered_ledgers": len(rows),
            "endpoint_self_recorded_from_left_count": len(endpoint_rows),
            "self_recorded_monotone_count": len(stack_rows),
            "self_recorded_monotone_ledgers": [row["ordered_ledger"] for row in stack_rows],
            "sum_square_minimizers_on_monotone_stack": minimizers(stack_rows, "sum_squares"),
            "parseval_minimizers_on_monotone_stack": minimizers(stack_rows, "parseval_non_dc_power_n5"),
            "power_minimizers_on_monotone_stack": power_winners,
            "all_checked_power_costs_uniquely_select_balanced": all(packet["winners"] == [list(BALANCED_LEDGER)] for packet in power_winners.values()),
        },
        "smoothing_certificate_on_monotone_stack": smoothing_steps_from_extreme_to_balanced(),
        "selector_consequence": {
            "what_self_record_plus_monotone_does": "It reduces the 35 ordered positive ledgers to exactly the three canonical descending ledgers.",
            "what_min_ripple_adds": "It uniquely selects (2,2,2,1,1) among those three under every checked strictly convex power cost p=2..6.",
            "what_remains_unproved": "The strict theory still has to derive monotone formation and min-ripple/convex action from nadsoliton geometry rather than assume them.",
        },
        "candidate_interpretation": {
            "supported_by_this_probe": True,
            "content": "Self-record + monotone formation + min-ripple is a finite conditional stack that selects the balanced ledger.",
            "why_this_is_more_proof_like": "The probe exhaustively enumerates the 35 ledgers, proves the three-ledger monotone residue, and records exact convex cost drops to the balanced ledger.",
            "why_this_is_not_enough": "The monotone formation law and convex/min-ripple action are not derived strict-core theorems here.",
            "status": "candidate-supported-but-monotone-min-ripple-premises-not-derived",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "The nadsoliton is treated as information carrying finite self-record patterns, but no separate sub-nadsoliton informational layer is introduced.",
            "Self-record plus monotone formation is a conditional selector stack, not a strict-core theorem.",
            "The min-ripple/symmetric convex action is assumed as a selector premise here and is not derived from strict nadsoliton geometry.",
            "No theorem derives the d5 support, monotone formation law, balanced ledger, endpoint-source convention, or strict source-localizing term from strict nadsoliton geometry.",
            "No theorem derives branch count m=5, total exponent n=8, denominator 3, or binary-rescale quotient.",
            "No Aut(Z_12)/N462/QW-2191 selector obstruction is discharged by this monotone min-ripple stack certificate.",
            "No theorem derives eta=9/5 without the branch model, denominator/quotient convention, support premise, monotone self-record premise, and min-ripple ledger selector premise.",
            "No QW-2191 discharge and no ToE closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Try to derive the monotone formation law and min-ripple convex action as internal self-record dynamics of nadsoliton formation, or prove a no-go for that derivation under current strict assumptions.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha self-recorded monotone min-ripple stack probe\n\n"
        "Status: conditional self-record + monotone + min-ripple stack; no strict selector discharge.\n\n"
        f"- Positive ordered ledgers: `{len(rows)}`; endpoint-left self-recorded: `{len(endpoint_rows)}`; monotone stack survivors: `{len(stack_rows)}`.\n"
        f"- Monotone survivors: `{[row['ordered_ledger'] for row in stack_rows]}`.\n"
        f"- Sum-square minimizer on stack: `{report['finite_stack_scan']['sum_square_minimizers_on_monotone_stack']}`.\n"
        f"- Target replay: `q^5={product.numerator}/{product.denominator}`, eta residual `{eta_from_product(product, SUPPORT_SIZE)-STRICT_TARGET_ETA:.3e}`.\n"
        "- Honest read: self-record + monotone formation narrows to three ledgers; min-ripple selects balanced, but both are premises.\n"
        "- No false pass: no strict monotone/min-ripple theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
