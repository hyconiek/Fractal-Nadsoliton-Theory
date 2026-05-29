#!/usr/bin/env python3
"""Scratch probe: entropy-selector discriminator for strict alpha -> eta.

The balanced-branch ledger probe showed that the radical ledger

    q^5 = 2^8/3^5

is recovered by the unique balanced positive-integer allocation of eight binary
units over five active branches: (2, 2, 2, 1, 1).  This packet asks whether that
balance selector has a sharper combinatorial reading.

Result: the balanced ledger is exactly selected by a fixed-labelled-branch
Shannon/multinomial criterion: it has maximal entropy for p_i=e_i/8 and maximal
labelled bit-assignment count 8!/prod_i(e_i!).  However, if the selector instead
uses canonical-orbit aggregate mass, the maximizer changes to (3, 2, 1, 1, 1).
So this is a discriminator, not a closure theorem: the next proof obligation is
to derive a fixed-labelled branch entropy selector, rather than an unlabeled
orbit-aggregate selector, from strict nadsoliton geometry.
"""
from __future__ import annotations

import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BALANCED = HERE / "bridge_strict_alpha_balanced_branch_ledger_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_entropy_selector_discriminator_report.json"
OUT_MD = HERE / "bridge_strict_alpha_entropy_selector_discriminator_report.md"

BRANCH_COUNT = 5
TARGET_BINARY_EXPONENT = 8
STRICT_TARGET_ETA = 9.0 / 5.0
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0
CANONICAL_LEDGERS = [
    (4, 1, 1, 1, 1),
    (3, 2, 1, 1, 1),
    (2, 2, 2, 1, 1),
]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def eta_from_correction(correction: float) -> float:
    return math.log(NAD12_SUPPORT_SIZE * correction) / math.log(ALPHA_SCALE)


def shannon_entropy_nats(ledger: tuple[int, ...]) -> float:
    total = sum(ledger)
    return -sum((e / total) * math.log(e / total) for e in ledger)


def normalized_entropy(ledger: tuple[int, ...]) -> float:
    return shannon_entropy_nats(ledger) / math.log(len(ledger))


def labelled_multinomial_count(ledger: tuple[int, ...]) -> int:
    numerator = math.factorial(sum(ledger))
    denominator = math.prod(math.factorial(e) for e in ledger)
    return numerator // denominator


def canonical_orbit_size(ledger: tuple[int, ...]) -> int:
    counts = Counter(ledger)
    denominator = math.prod(math.factorial(multiplicity) for multiplicity in counts.values())
    return math.factorial(len(ledger)) // denominator


def orbit_aggregate_count(ledger: tuple[int, ...]) -> int:
    return canonical_orbit_size(ledger) * labelled_multinomial_count(ledger)


def product_fraction(ledger: tuple[int, ...]) -> Fraction:
    return Fraction(2 ** sum(ledger), 3 ** len(ledger))


def row_for_ledger(ledger: tuple[int, ...]) -> dict[str, Any]:
    product = product_fraction(ledger)
    correction = float(product) ** (1.0 / len(ledger))
    eta = eta_from_correction(correction)
    return {
        "ledger": list(ledger),
        "shannon_entropy_nats": shannon_entropy_nats(ledger),
        "normalized_entropy_vs_log_5": normalized_entropy(ledger),
        "labelled_multinomial_count": labelled_multinomial_count(ledger),
        "canonical_orbit_size": canonical_orbit_size(ledger),
        "orbit_aggregate_count": orbit_aggregate_count(ledger),
        "product_fraction": f"{product.numerator}/{product.denominator}",
        "eta_from_product": eta,
        "eta_residual_vs_9_5": eta - STRICT_TARGET_ETA,
    }


def argmax_unique(rows: list[dict[str, Any]], key: str) -> tuple[dict[str, Any], bool]:
    max_value = max(row[key] for row in rows)
    winners = [row for row in rows if row[key] == max_value]
    return winners[0], len(winners) == 1


def main() -> None:
    balanced_report = load_json(BALANCED)
    selected_balanced = tuple(balanced_report["selected_balanced_ledger"]["ledger"])
    rows = [row_for_ledger(ledger) for ledger in CANONICAL_LEDGERS]
    entropy_winner, entropy_unique = argmax_unique(rows, "shannon_entropy_nats")
    labelled_winner, labelled_unique = argmax_unique(rows, "labelled_multinomial_count")
    orbit_winner, orbit_unique = argmax_unique(rows, "orbit_aggregate_count")

    report = {
        "status": "OPEN_STRICT_ALPHA_ENTROPY_SELECTOR_DISCRIMINATOR_NO_STRICT_SELECTOR_THEOREM",
        "result_kind": "SCRATCH_STRICT_ALPHA_ENTROPY_SELECTOR_DISCRIMINATOR_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "balanced_branch_ledger": str(BALANCED.relative_to(ROOT)),
        },
        "candidate_space": {
            "branch_count": BRANCH_COUNT,
            "total_binary_units": TARGET_BINARY_EXPONENT,
            "canonical_positive_ledgers": [list(ledger) for ledger in CANONICAL_LEDGERS],
            "all_rows_share_product": "2^8/3^5 = 256/243",
            "all_rows_share_eta": STRICT_TARGET_ETA,
        },
        "selector_rows": rows,
        "fixed_labelled_branch_entropy_selector": {
            "entropy_winner": entropy_winner,
            "entropy_winner_unique": entropy_unique,
            "labelled_multinomial_winner": labelled_winner,
            "labelled_multinomial_winner_unique": labelled_unique,
            "matches_balanced_ledger": entropy_winner["ledger"] == list(selected_balanced) == labelled_winner["ledger"],
            "interpretation": "If the five branches are fixed labelled channels, Shannon entropy / labelled bit-assignment multiplicity selects the balanced ledger (2,2,2,1,1).",
        },
        "unlabelled_orbit_aggregate_discriminator": {
            "orbit_aggregate_winner": orbit_winner,
            "orbit_aggregate_winner_unique": orbit_unique,
            "matches_balanced_ledger": orbit_winner["ledger"] == list(selected_balanced),
            "obstruction": "If the selector maximizes canonical-orbit aggregate mass, it selects (3,2,1,1,1), not the balanced ledger.",
        },
        "candidate_interpretation": {
            "supported_by_this_probe": bool(
                entropy_unique
                and labelled_unique
                and orbit_unique
                and entropy_winner["ledger"] == [2, 2, 2, 1, 1]
                and labelled_winner["ledger"] == [2, 2, 2, 1, 1]
                and orbit_winner["ledger"] == [3, 2, 1, 1, 1]
            ),
            "content": "The balance premise can be sharpened into a fixed-labelled Shannon/multinomial selector, but not into an unlabelled orbit-aggregate selector.",
            "why_this_is_more_proof_like": "It gives a concrete combinatorial discriminator for the missing selector convention: labelled entropy supports the radical ledger, while orbit aggregation would pick a different ledger.",
            "why_this_is_not_enough": "No strict theorem proves that nadsoliton branches are fixed labelled entropy channels rather than unlabelled orbit states.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "No theorem derives fixed-labelled branch entropy as the strict-core selector.",
            "Unlabelled orbit-aggregate maximization would select (3,2,1,1,1), so the selector convention is mathematically material.",
            "No theorem derives eta=9/5 without adopting the branch count, ternary normalization, and fixed-labelled entropy selector as extra premises.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Try to prove that strict nadsoliton branches are fixed labelled Shannon channels; if not, record the orbit-aggregate discriminator as a selector obstruction rather than claiming closure.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha entropy selector discriminator probe\n\n"
        "Status: selector discriminator for eta=9/5; no strict selector theorem.\n\n"
        "- Candidate space: positive five-branch ledgers summing to `8`: `(4,1,1,1,1)`, `(3,2,1,1,1)`, `(2,2,2,1,1)`.\n"
        f"- Fixed-labelled entropy winner: `{entropy_winner['ledger']}` with labelled multinomial count `{labelled_winner['labelled_multinomial_count']}`.\n"
        f"- Unlabelled orbit-aggregate winner: `{orbit_winner['ledger']}` with aggregate count `{orbit_winner['orbit_aggregate_count']}`.\n"
        "- Honest read: labelled Shannon/multinomial selection supports the balanced ledger, but orbit aggregation selects a different ledger; the selector convention is a real proof obligation.\n"
        "- No false pass: no fixed-labelled strict selector theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
