#!/usr/bin/env python3
"""Scratch probe: Fourier minimum-ripple convexity certificate for strict alpha -> eta.

The previous phase/Fourier discriminator showed that the balanced alpha ledger

    (2, 2, 2, 1, 1)

is not selected by a naive "highest resonance power" criterion.  It is selected
by a *minimum non-DC Fourier ripple / coherence-compression* criterion.  This
packet makes that conditional selector more theorem-like inside the finite
integer branch model.

Certificate proved here (model-internal, not physical): for fixed branch count
m and fixed positive integer total n, Parseval gives

    total_non_dc_power_N(e) = N * sum_i e_i^2 - n^2.

Thus minimizing Fourier ripple is exactly minimizing sum_i e_i^2.  A pairwise
smoothing move

    (..., a, ..., b, ...) -> (..., a-1, ..., b+1, ...)  for a >= b + 2

preserves m and n and reduces total_non_dc_power_N by exactly

    2 * N * (a - b - 1) > 0.

Repeated smoothing therefore terminates uniquely, up to permutation, at the
balanced ledger with entries floor(n/m) and ceil(n/m).  For n=8,m=5 this is
(2,2,2,1,1).  This discharges the *integer convexity* part of the selector only
after the minimum-ripple premise, branch count, total exponent, denominator, and
phase-labelled branch model are already supplied.  It does not derive those
premises from strict nadsoliton geometry and does not discharge QW-2191.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FOURIER_DISCRIMINATOR = HERE / "bridge_strict_alpha_phase_fourier_selector_discriminator_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_fourier_min_ripple_convexity_certificate_report.json"
OUT_MD = HERE / "bridge_strict_alpha_fourier_min_ripple_convexity_certificate_report.md"

TARGET_BRANCH_COUNT = 5
TARGET_BINARY_EXPONENT = 8
DENOMINATOR = 3
Z12_NODE_COUNT = 12
FIVE_BRANCH_DFT_SIZE = 5
STRICT_TARGET_ETA = 9.0 / 5.0
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0
SCAN_MAX_BRANCH_COUNT = 9
SCAN_MAX_TOTAL_EXPONENT = 24


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def eta_from_product(product: Fraction, branch_count: int) -> float:
    correction = float(product) ** (1.0 / branch_count)
    return math.log(NAD12_SUPPORT_SIZE * correction) / math.log(ALPHA_SCALE)


def total_non_dc_power(ring_size: int, ledger: tuple[int, ...]) -> int:
    total = sum(ledger)
    return ring_size * sum(value * value for value in ledger) - total * total


def canonical_desc(values: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(values, reverse=True))


def balanced_ledger(branch_count: int, total: int) -> tuple[int, ...]:
    base, remainder = divmod(total, branch_count)
    if base <= 0:
        raise ValueError("positive ledgers require total >= branch_count")
    return tuple([base + 1] * remainder + [base] * (branch_count - remainder))


def positive_partitions(total: int, parts: int, max_part: int | None = None) -> list[tuple[int, ...]]:
    """Return nonincreasing positive integer partitions of total into parts."""
    if parts == 0:
        return [()] if total == 0 else []
    if max_part is None:
        max_part = total
    out: list[tuple[int, ...]] = []
    min_remaining = parts - 1
    for first in range(min(max_part, total - min_remaining), 0, -1):
        for rest in positive_partitions(total - first, parts - 1, min(first, total - first)):
            if len(rest) == parts - 1:
                out.append((first,) + rest)
    return out


def smoothing_step(ledger: tuple[int, ...]) -> dict[str, Any] | None:
    values = list(ledger)
    high_index = max(range(len(values)), key=lambda index: values[index])
    low_index = min(range(len(values)), key=lambda index: values[index])
    high = values[high_index]
    low = values[low_index]
    gap = high - low
    if gap < 2:
        return None
    before = tuple(values)
    values[high_index] -= 1
    values[low_index] += 1
    after = canonical_desc(tuple(values))
    return {
        "before": list(before),
        "after": list(after),
        "high": high,
        "low": low,
        "gap": gap,
        "sum_squares_before": sum(value * value for value in before),
        "sum_squares_after": sum(value * value for value in after),
        "sum_squares_drop": 2 * (gap - 1),
        "z12_power_drop": 2 * Z12_NODE_COUNT * (gap - 1),
        "five_branch_power_drop": 2 * len(before) * (gap - 1),
    }


def smoothing_path(ledger: tuple[int, ...]) -> list[dict[str, Any]]:
    current = canonical_desc(ledger)
    path = []
    while True:
        step = smoothing_step(current)
        if step is None:
            break
        path.append(step)
        current = tuple(step["after"])
    return path


def target_certificate() -> dict[str, Any]:
    target = balanced_ledger(TARGET_BRANCH_COUNT, TARGET_BINARY_EXPONENT)
    canonical_ledgers = positive_partitions(TARGET_BINARY_EXPONENT, TARGET_BRANCH_COUNT)
    rows = []
    for ledger in canonical_ledgers:
        product = Fraction(2 ** sum(ledger), DENOMINATOR ** len(ledger))
        path = smoothing_path(ledger)
        rows.append(
            {
                "ledger": list(ledger),
                "is_balanced": ledger == target,
                "sum_squares": sum(value * value for value in ledger),
                "z12_total_non_dc_power": total_non_dc_power(Z12_NODE_COUNT, ledger),
                "five_branch_total_non_dc_power": total_non_dc_power(FIVE_BRANCH_DFT_SIZE, ledger),
                "smoothing_steps_to_balanced": len(path),
                "terminal_after_smoothing": list(path[-1]["after"] if path else ledger),
                "product_fraction": f"{product.numerator}/{product.denominator}",
                "eta_residual_vs_9_5": eta_from_product(product, len(ledger)) - STRICT_TARGET_ETA,
            }
        )
    min_sum_squares = min(row["sum_squares"] for row in rows)
    minimizers = [row for row in rows if row["sum_squares"] == min_sum_squares]
    hardest_unbalanced = max((row for row in rows if not row["is_balanced"]), key=lambda row: row["sum_squares"])
    return {
        "branch_count": TARGET_BRANCH_COUNT,
        "total_binary_exponent": TARGET_BINARY_EXPONENT,
        "balanced_ledger": list(target),
        "canonical_positive_ledgers": [list(ledger) for ledger in canonical_ledgers],
        "rows": rows,
        "unique_minimizer": len(minimizers) == 1,
        "minimizer": minimizers[0],
        "hardest_unbalanced_row": hardest_unbalanced,
        "sample_smoothing_paths": {
            "from_4_1_1_1_1": smoothing_path((4, 1, 1, 1, 1)),
            "from_3_2_1_1_1": smoothing_path((3, 2, 1, 1, 1)),
            "from_2_2_2_1_1": smoothing_path((2, 2, 2, 1, 1)),
        },
    }


def bounded_convexity_scan() -> dict[str, Any]:
    rows = []
    checked_pairs = 0
    failures = []
    for branch_count in range(1, SCAN_MAX_BRANCH_COUNT + 1):
        for total in range(branch_count, SCAN_MAX_TOTAL_EXPONENT + 1):
            ledgers = positive_partitions(total, branch_count)
            expected = balanced_ledger(branch_count, total)
            scored = [(sum(value * value for value in ledger), ledger) for ledger in ledgers]
            min_score = min(score for score, _ledger in scored)
            minimizers = [ledger for score, ledger in scored if score == min_score]
            for ledger in ledgers:
                path = smoothing_path(ledger)
                terminal = tuple(path[-1]["after"] if path else ledger)
                if terminal != expected:
                    failures.append({"branch_count": branch_count, "total": total, "ledger": list(ledger), "terminal": list(terminal), "expected": list(expected)})
                checked_pairs += 1
            row = {
                "branch_count": branch_count,
                "total": total,
                "partition_count": len(ledgers),
                "balanced_ledger": list(expected),
                "min_sum_squares": min_score,
                "unique_minimizer": minimizers == [expected],
                "minimizer": list(minimizers[0]),
                "max_smoothing_steps": max(len(smoothing_path(ledger)) for ledger in ledgers),
            }
            if not row["unique_minimizer"]:
                failures.append({"branch_count": branch_count, "total": total, "minimizers": [list(ledger) for ledger in minimizers], "expected": list(expected)})
            rows.append(row)
    return {
        "scan_limits": {"max_branch_count": SCAN_MAX_BRANCH_COUNT, "max_total_exponent": SCAN_MAX_TOTAL_EXPONENT},
        "cases_checked": len(rows),
        "canonical_ledgers_checked": checked_pairs,
        "all_unique_balanced_minimizers": not failures,
        "failures": failures,
        "rows": rows,
    }


def main() -> None:
    fourier_report = load_json(FOURIER_DISCRIMINATOR)
    target = target_certificate()
    scan = bounded_convexity_scan()
    product = Fraction(2 ** TARGET_BINARY_EXPONENT, DENOMINATOR ** TARGET_BRANCH_COUNT)

    report = {
        "status": "OPEN_STRICT_ALPHA_FOURIER_MIN_RIPPLE_CONVEXITY_CERTIFICATE_NO_PHYSICAL_SELECTOR_THEOREM",
        "result_kind": "SCRATCH_STRICT_ALPHA_FOURIER_MIN_RIPPLE_CONVEXITY_CERTIFICATE_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "phase_fourier_selector_discriminator": str(FOURIER_DISCRIMINATOR.relative_to(ROOT)),
        },
        "previous_discriminator_replay": {
            "result_kind": fourier_report["result_kind"],
            "minimum_fourier_ripple_winner": fourier_report["parseval_discriminator"]["minimum_fourier_ripple_winner"]["ledger"],
            "maximum_fourier_ripple_winner": fourier_report["parseval_discriminator"]["maximum_fourier_ripple_winner"]["ledger"],
        },
        "conditional_theorem_statement": {
            "premises": [
                "fixed positive integer branch count m",
                "fixed total binary exponent n",
                "phase-labelled branch vector e=(e_1,...,e_m)",
                "minimum total non-DC Fourier ripple criterion",
            ],
            "parseval_identity": "P_N(e)=N*sum_i(e_i^2)-n^2",
            "pairwise_smoothing_rule": "if a>=b+2, replace (a,b) by (a-1,b+1)",
            "exact_drop": "P_N(before)-P_N(after)=2*N*(a-b-1)>0",
            "conclusion": "unique minimizer up to permutation is the balanced floor/ceil ledger",
            "physical_status": "integer-convexity certificate only; minimum-ripple principle is not derived here",
        },
        "target_eta_9_5_certificate": {
            "q_power_product": f"{product.numerator}/{product.denominator}",
            "branch_count": TARGET_BRANCH_COUNT,
            "total_binary_exponent": TARGET_BINARY_EXPONENT,
            "eta_residual_vs_9_5": eta_from_product(product, TARGET_BRANCH_COUNT) - STRICT_TARGET_ETA,
            **target,
        },
        "bounded_convexity_scan": scan,
        "candidate_interpretation": {
            "supported_by_this_probe": bool(target["unique_minimizer"] and scan["all_unique_balanced_minimizers"]),
            "content": "Given the minimum-ripple/coherence-compression selector premise, the balanced ledger (2,2,2,1,1) is not a guess: it is the unique integer convexity minimizer for m=5,n=8.",
            "why_this_is_more_proof_like": "The proof obligation is reduced to a finite pairwise smoothing lemma with exact Parseval power drops and a bounded independent exhaustive scan.",
            "why_this_is_not_enough": "The probe does not derive the minimum-ripple variational principle, branch count m=5, total exponent n=8, denominator 3, binary-rescale quotient, phase labels, or strict source of the selector from nadsoliton geometry.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "No theorem derives the minimum-ripple/coherence-compression principle from strict nadsoliton geometry.",
            "No theorem derives branch count m=5, total exponent n=8, denominator 3, or phase-labelled branches from strict core in this packet.",
            "This certificate is conditional on adopting the minimum Fourier ripple selector; it is not a strict-core selector source.",
            "No Aut(Z_12)/N462/QW-2191 selector obstruction is discharged.",
            "No theorem derives eta=9/5 without the branch model, denominator/quotient convention, phase labelling, and minimum-ripple selector premises.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Try to derive the minimum-ripple/coherence-compression variational principle from a strict nadsoliton action or export it explicitly as a non-strict selector axiom; do not call the convexity certificate a physical selector theorem by itself.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Fourier minimum-ripple convexity certificate probe\n\n"
        "Status: integer convexity certificate for the phase/Fourier minimum-ripple selector; no physical selector theorem.\n\n"
        "- Conditional theorem: for fixed `m,n`, Parseval reduces minimum non-DC Fourier ripple to minimizing `sum(e_i^2)`.\n"
        "- Pairwise smoothing certificate: `(a,b)->(a-1,b+1)` for `a>=b+2` drops ripple by `2*N*(a-b-1)`.\n"
        f"- Target result: for `m={TARGET_BRANCH_COUNT}, n={TARGET_BINARY_EXPONENT}`, the unique canonical minimizer is `{list(balanced_ledger(TARGET_BRANCH_COUNT, TARGET_BINARY_EXPONENT))}` and `q^5={product.numerator}/{product.denominator}`.\n"
        f"- Bounded scan: `{scan['cases_checked']}` `(m,n)` cases and `{scan['canonical_ledgers_checked']}` canonical ledgers checked; all had unique balanced minimizers: `{scan['all_unique_balanced_minimizers']}`.\n"
        "- Honest read: this proves the integer convexity part only after the minimum-ripple premise is supplied; it does not derive the strict-core selector source.\n"
        "- No false pass: no strict phase/Fourier selector theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
