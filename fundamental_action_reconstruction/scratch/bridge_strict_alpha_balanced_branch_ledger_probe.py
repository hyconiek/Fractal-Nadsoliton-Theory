#!/usr/bin/env python3
"""Scratch probe: balanced branch ledger refinement for strict alpha -> eta.

The binary/ternary ledger probe introduced a narrow one-or-two-bit branch rule
for the exact radical target

    q^5 = 256/243 = 2^8/3^5.

This refinement removes one arbitrary-looking part of that rule.  Instead of
assuming each branch has one or two binary bits, it assumes only:

    1. five active branches,
    2. denominator 3 on each branch,
    3. a positive integer binary exponent on each branch,
    4. total binary exponent 8,
    5. a balance/equipartition selector minimizing variance.

Under those assumptions, the one-or-two-bit ledger is forced: among all positive
integer five-branch ledgers summing to 8, the unique variance-minimizing
canonical multiset is (2, 2, 2, 1, 1).  This exactly recovers the radical
correction.  The balance selector itself is not derived here; it remains the new
strict-side proof obligation, especially in view of QW-2191.
"""
from __future__ import annotations

import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BINARY_TERNARY = HERE / "bridge_strict_alpha_binary_ternary_ledger_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_balanced_branch_ledger_report.json"
OUT_MD = HERE / "bridge_strict_alpha_balanced_branch_ledger_report.md"

BRANCH_COUNT = 5
TERNARY_DENOMINATOR_PER_BRANCH = 3
TARGET_BINARY_EXPONENT = 8
STRICT_TARGET_ETA = Fraction(9, 5)
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def canonical_multiset(seq: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(seq, reverse=True))


def positive_labelled_ledgers() -> list[tuple[int, ...]]:
    ledgers = []
    for seq in itertools.product(range(1, TARGET_BINARY_EXPONENT + 1), repeat=BRANCH_COUNT):
        if sum(seq) == TARGET_BINARY_EXPONENT:
            ledgers.append(seq)
    return ledgers


def canonical_positive_ledgers() -> list[tuple[int, ...]]:
    return sorted({canonical_multiset(seq) for seq in positive_labelled_ledgers()}, reverse=True)


def variance_score(ledger: tuple[int, ...]) -> Fraction:
    mean = Fraction(TARGET_BINARY_EXPONENT, BRANCH_COUNT)
    return sum((Fraction(e) - mean) ** 2 for e in ledger) / BRANCH_COUNT


def max_gap(ledger: tuple[int, ...]) -> int:
    return max(ledger) - min(ledger)


def product_fraction(ledger: tuple[int, ...]) -> Fraction:
    numerator_exponent = sum(ledger)
    denominator_exponent = len(ledger)
    return Fraction(2**numerator_exponent, TERNARY_DENOMINATOR_PER_BRANCH**denominator_exponent)


def eta_from_correction(correction: float) -> float:
    return math.log(NAD12_SUPPORT_SIZE * correction) / math.log(ALPHA_SCALE)


def row_for_ledger(ledger: tuple[int, ...]) -> dict[str, Any]:
    score = variance_score(ledger)
    product = product_fraction(ledger)
    correction = float(product) ** (1.0 / BRANCH_COUNT)
    eta = eta_from_correction(correction)
    return {
        "ledger": list(ledger),
        "max_gap": max_gap(ledger),
        "variance_score": float(score),
        "variance_score_fraction": f"{score.numerator}/{score.denominator}",
        "product_fraction": f"{product.numerator}/{product.denominator}",
        "geometric_mean_correction": correction,
        "eta_from_ledger": eta,
        "eta_residual_vs_9_5": eta - float(STRICT_TARGET_ETA),
    }


def main() -> None:
    binary_report = load_json(BINARY_TERNARY)
    q_radical = binary_report["exact_target_replay"]["q_radical"]
    canonical_ledgers = canonical_positive_ledgers()
    labelled_count = len(positive_labelled_ledgers())
    rows = [row_for_ledger(ledger) for ledger in canonical_ledgers]
    min_variance = min(variance_score(ledger) for ledger in canonical_ledgers)
    min_gap = min(max_gap(ledger) for ledger in canonical_ledgers)
    variance_minimizers = [ledger for ledger in canonical_ledgers if variance_score(ledger) == min_variance]
    gap_minimizers = [ledger for ledger in canonical_ledgers if max_gap(ledger) == min_gap]
    selected_ledger = variance_minimizers[0]
    selected_product = product_fraction(selected_ledger)
    selected_correction = float(selected_product) ** (1.0 / BRANCH_COUNT)
    selected_eta = eta_from_correction(selected_correction)

    report = {
        "status": "OPEN_STRICT_ALPHA_BALANCED_BRANCH_LEDGER_WITNESS_NO_SELECTOR_DERIVATION_THEOREM",
        "result_kind": "SCRATCH_STRICT_ALPHA_BALANCED_BRANCH_LEDGER_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "binary_ternary_ledger": str(BINARY_TERNARY.relative_to(ROOT)),
        },
        "exact_target_replay": {
            "q_radical": q_radical,
            "target_product_label": "q^5 = 2^8/3^5 = 256/243",
            "target_binary_exponent_total": TARGET_BINARY_EXPONENT,
            "target_branch_count": BRANCH_COUNT,
            "target_mean_binary_exponent_per_branch": TARGET_BINARY_EXPONENT / BRANCH_COUNT,
        },
        "weakened_assumptions_compared_to_previous_probe": {
            "previous_extra_assumption": "each branch exponent is in {1,2}",
            "replacement_assumption": "each branch exponent is a positive integer and a balance selector minimizes variance/equipartition defect",
            "forced_consequence": "forced: the selected balanced ledger automatically has entries in {1,2}",
            "still_not_derived": "the five active branches, denominator-3 normalization, and balance selector are not strict-core theorems here",
        },
        "positive_integer_ledger_scan": {
            "labelled_positive_ledgers_count": labelled_count,
            "canonical_positive_ledgers_count": len(canonical_ledgers),
            "canonical_rows": rows,
            "min_variance_score_fraction": f"{min_variance.numerator}/{min_variance.denominator}",
            "variance_minimizers": [list(ledger) for ledger in variance_minimizers],
            "variance_minimizer_unique": len(variance_minimizers) == 1,
            "min_max_gap": min_gap,
            "max_gap_minimizers": [list(ledger) for ledger in gap_minimizers],
            "max_gap_minimizer_unique": len(gap_minimizers) == 1,
        },
        "selected_balanced_ledger": {
            "ledger": list(selected_ledger),
            "selection_rule": "unique minimizer of integer branch variance, equivalently unique minimizer of max-min gap, among positive five-branch ledgers summing to 8",
            "branch_ratio_labels": [f"{2**e}/3" for e in selected_ledger],
            "product_fraction": f"{selected_product.numerator}/{selected_product.denominator}",
            "product_label": "(4/3)^3 * (2/3)^2 = 2^8/3^5 = 256/243",
            "geometric_mean_correction": selected_correction,
            "geometric_mean_residual_vs_q_radical": selected_correction - q_radical,
            "eta_from_ledger": selected_eta,
            "eta_residual_vs_9_5": selected_eta - float(STRICT_TARGET_ETA),
        },
        "candidate_interpretation": {
            "supported_by_this_probe": bool(
                len(variance_minimizers) == 1
                and len(gap_minimizers) == 1
                and selected_ledger == (2, 2, 2, 1, 1)
                and selected_product == Fraction(256, 243)
                and abs(selected_eta - float(STRICT_TARGET_ETA)) < 1e-12
            ),
            "content": "The previous one-or-two-bit ledger is derivable from a weaker positive-integer plus balance-selector premise.",
            "why_this_is_more_proof_like": "It turns an ad hoc-looking exponent alphabet {1,2} into the unique balanced integer allocation of eight binary units over five active branches.",
            "why_this_is_not_enough": "The balance/equipartition selector is not derived as a strict nadsoliton source and does not discharge QW-2191.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "No theorem derives the five active branches from strict nadsoliton geometry.",
            "No theorem derives the denominator-3 branch normalization from strict nadsoliton geometry.",
            "No theorem derives the balance/equipartition selector as a strict-core selector source.",
            "No theorem derives eta=9/5 without adopting the branch count, ternary normalization, and balance selector as extra premises.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Try to source the balance/equipartition selector from strict nadsoliton geometry or an explicit non-strict selector premise; otherwise keep eta=9/5 as gate-selected with an exact balanced-ledger shadow.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha balanced branch ledger probe\n\n"
        "Status: exact balanced-ledger witness for eta=9/5; no selector derivation theorem.\n\n"
        f"- Weakened rule: replace explicit exponent alphabet `{{1,2}}` with positive integer exponents plus a balance/equipartition selector over `{BRANCH_COUNT}` branches summing to `{TARGET_BINARY_EXPONENT}`.\n"
        f"- Scan size: `{labelled_count}` labelled positive ledgers and `{len(canonical_ledgers)}` canonical ledgers.\n"
        f"- Unique variance and max-gap minimizer: `{list(selected_ledger)}` with variance `{min_variance.numerator}/{min_variance.denominator}`.\n"
        f"- Exact replay: `{report['selected_balanced_ledger']['product_label']}`, correction residual `{selected_correction - q_radical:.3e}`, eta residual `{selected_eta - float(STRICT_TARGET_ETA):.3e}`.\n"
        "- No false pass: branch count, denominator-3 normalization, and balance selector are still premises; no QW-2191 discharge and no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
