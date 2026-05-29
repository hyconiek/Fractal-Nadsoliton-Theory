#!/usr/bin/env python3
"""Scratch probe: binary/ternary exponent ledger for strict alpha -> eta.

The radical correction probe gave the exact algebraic target

    q_radical^5 = 256/243 = 2^8 / 3^5.

This packet asks one more proof-like question: if the fifth root is read as a
five-branch strict-side compression average, can the prime-exponent ledger be
resolved by a tiny integer branch law?

Under the narrow candidate rule

    each of five branches has a ternary denominator 3,
    each branch carries either one or two binary numerator bits 2^e,

exact eta=9/5 is equivalent to distributing eight binary bits over five branches.
The only multiset solution is

    e = (2, 2, 2, 1, 1),
    product_i (2^e_i / 3) = 2^8 / 3^5,
    geometric_mean_i (2^e_i / 3) = (256/243)^(1/5).

This is not a derivation of the branch law.  It is an exact ledger witness and a
sharper proof obligation: derive the five branches, the per-branch ternary
normalization, and the {2,2,2,1,1} bit split from strict nadsoliton geometry, or
keep eta=9/5 classified as gate-selected with an exact ledger shadow.
"""
from __future__ import annotations

import itertools
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RADICAL = HERE / "bridge_strict_alpha_radical_correction_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_binary_ternary_ledger_report.json"
OUT_MD = HERE / "bridge_strict_alpha_binary_ternary_ledger_report.md"

BRANCH_COUNT = 5
TERNARY_DENOMINATOR_PER_BRANCH = 3
TARGET_BINARY_EXPONENT = 8
STRICT_TARGET_ETA = 9.0 / 5.0
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def eta_from_correction(correction: float) -> float:
    return math.log(NAD12_SUPPORT_SIZE * correction) / math.log(ALPHA_SCALE)


def canonical_multiset(seq: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(seq, reverse=True))


def enumerate_positive_one_two_ledgers() -> list[tuple[int, ...]]:
    ledgers = set()
    for seq in itertools.product((1, 2), repeat=BRANCH_COUNT):
        if sum(seq) == TARGET_BINARY_EXPONENT:
            ledgers.add(canonical_multiset(seq))
    return sorted(ledgers, reverse=True)


def enumerate_small_unbounded_ledgers(max_exponent: int = 4) -> list[tuple[int, ...]]:
    ledgers = set()
    for seq in itertools.product(range(max_exponent + 1), repeat=BRANCH_COUNT):
        if sum(seq) == TARGET_BINARY_EXPONENT:
            ledgers.add(canonical_multiset(seq))
    return sorted(ledgers, reverse=True)


def branch_rows(ledger: tuple[int, ...]) -> list[dict[str, Any]]:
    rows = []
    for idx, exponent in enumerate(ledger, start=1):
        numerator = 2**exponent
        ratio = numerator / TERNARY_DENOMINATOR_PER_BRANCH
        rows.append(
            {
                "branch_index": idx,
                "binary_exponent_e_i": exponent,
                "branch_ratio_label": f"{numerator}/3",
                "branch_ratio": ratio,
                "above_one": ratio > 1.0,
            }
        )
    return rows


def main() -> None:
    radical_report = load_json(RADICAL)
    q_radical = radical_report["exact_radical_target"]["q_radical"]
    target_product = q_radical**BRANCH_COUNT

    constrained_ledgers = enumerate_positive_one_two_ledgers()
    all_small_ledgers = enumerate_small_unbounded_ledgers(max_exponent=4)
    selected_ledger = constrained_ledgers[0]
    selected_counter = Counter(selected_ledger)
    selected_product = math.prod((2**e) / TERNARY_DENOMINATOR_PER_BRANCH for e in selected_ledger)
    selected_geometric_mean = selected_product ** (1.0 / BRANCH_COUNT)
    eta_from_ledger = eta_from_correction(selected_geometric_mean)

    report = {
        "status": "OPEN_STRICT_ALPHA_BINARY_TERNARY_LEDGER_WITNESS_NO_BRANCH_DERIVATION_THEOREM",
        "result_kind": "SCRATCH_STRICT_ALPHA_BINARY_TERNARY_LEDGER_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "radical_correction": str(RADICAL.relative_to(ROOT)),
        },
        "exact_target_replay": {
            "q_radical": q_radical,
            "q_radical_power_5": target_product,
            "target_prime_ledger": "q^5 = 2^8 / 3^5 = 256/243",
            "target_binary_exponent_total": TARGET_BINARY_EXPONENT,
            "target_ternary_exponent_total": -BRANCH_COUNT,
        },
        "candidate_branch_rule": {
            "branch_count": BRANCH_COUNT,
            "per_branch_denominator": TERNARY_DENOMINATOR_PER_BRANCH,
            "allowed_binary_exponents_per_branch": [1, 2],
            "meaning": "five branches; each branch has ternary normalization 3 and either one or two binary numerator bits",
            "not_derived_here": True,
        },
        "ledger_solution_under_candidate_rule": {
            "solutions_as_canonical_multisets": [list(row) for row in constrained_ledgers],
            "unique_under_rule": len(constrained_ledgers) == 1,
            "selected_binary_exponent_multiset": list(selected_ledger),
            "selected_count_by_exponent": {str(k): selected_counter[k] for k in sorted(selected_counter, reverse=True)},
            "branch_rows": branch_rows(selected_ledger),
            "product_branch_ratios": selected_product,
            "product_label": "(4/3)^3 * (2/3)^2 = 2^8/3^5 = 256/243",
            "geometric_mean_correction": selected_geometric_mean,
            "geometric_mean_residual_vs_q_radical": selected_geometric_mean - q_radical,
            "eta_from_ledger": eta_from_ledger,
            "eta_residual_vs_9_5": eta_from_ledger - STRICT_TARGET_ETA,
        },
        "discriminator_against_looser_ledgers": {
            "max_exponent_in_loose_scan": 4,
            "number_of_canonical_ledgers_sum_8_length_5_with_entries_0_to_4": len(all_small_ledgers),
            "sample_loose_ledgers": [list(row) for row in all_small_ledgers[:12]],
            "why_candidate_rule_matters": "Without the one-or-two-bit branch axiom, many exponent ledgers give the same prime product; the strict proof obligation is therefore the branch law, not just the arithmetic identity.",
        },
        "candidate_interpretation": {
            "supported_by_this_probe": bool(
                len(constrained_ledgers) == 1
                and selected_ledger == (2, 2, 2, 1, 1)
                and abs(selected_geometric_mean - q_radical) < 1e-15
                and abs(eta_from_ledger - STRICT_TARGET_ETA) < 1e-12
            ),
            "content": "Exact eta=9/5 can be represented as the geometric mean of five branch ratios with the unique one-or-two-bit ledger (2,2,2,1,1).",
            "why_this_is_more_proof_like": "It translates the radical target into an exact integer prime-exponent ledger with a uniqueness statement under a stated finite branch rule.",
            "why_this_is_not_enough": "The finite branch rule itself is not derived from strict nadsoliton geometry, Shannon compression, or a QW-2191 selector source.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "No theorem derives the five-branch denominator-3, one-or-two-bit branch law from strict nadsoliton geometry.",
            "No theorem derives eta=9/5 without adopting the candidate branch law or the radical correction as an extra premise.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Search for a strict-side source of the five branches and the one-or-two-bit split (2,2,2,1,1); absent that, keep the ledger as an exact candidate shadow rather than a theorem.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha binary/ternary ledger probe\n\n"
        "Status: exact branch-ledger witness for eta=9/5; no branch-law derivation theorem.\n\n"
        f"- Exact replay: `q^5={target_product:.15f}=256/243=2^8/3^5`.\n"
        "- Candidate finite branch law: five branches, denominator `3` per branch, and binary exponents only `1` or `2`.\n"
        f"- Unique ledger under that rule: `{list(selected_ledger)}`, i.e. `(4/3)^3*(2/3)^2`, with correction residual `{selected_geometric_mean - q_radical:.3e}` and eta residual `{eta_from_ledger - STRICT_TARGET_ETA:.3e}`.\n"
        f"- Looser scan warning: `{len(all_small_ledgers)}` canonical ledgers exist for length five, sum eight, entries 0..4; uniqueness requires the extra one-or-two-bit branch rule.\n"
        "- No false pass: no kernel identity, no legacy role transfer, no strict branch-law theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
