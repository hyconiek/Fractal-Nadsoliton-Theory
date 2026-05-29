#!/usr/bin/env python3
"""Scratch probe: Robin-Hood entropy certificate for strict alpha -> eta.

The entropy-selector discriminator showed that fixed-labelled Shannon/multinomial
selection chooses the balanced ledger (2,2,2,1,1), while orbit aggregation would
choose a different ledger.  This probe strengthens the fixed-labelled side from
an enumerated argmax into a constructive certificate.

For a fixed-labelled bit-assignment weight

    W(e_1,...,e_5) = 8! / prod_i(e_i!),

moving one binary unit from a branch with exponent a to a branch with exponent b
changes the weight by the exact ratio

    W_after / W_before = a / (b + 1).

Whenever a >= b + 2, that ratio is > 1.  Thus a Robin-Hood balancing move
strictly increases labelled multiplicity and also increases Shannon entropy.  In
the five-branch, total-eight ledger space, the two non-balanced canonical
ledgers have a monotone certificate path into (2,2,2,1,1):

    (4,1,1,1,1) -> (3,2,1,1,1) -> (2,2,2,1,1).

This is still not a strict selector theorem.  It proves only that if the missing
selector is fixed-labelled entropy/multinomial multiplicity, then the balanced
ledger is constructively forced.
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
ENTROPY = HERE / "bridge_strict_alpha_entropy_selector_discriminator_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_entropy_robin_hood_certificate_report.json"
OUT_MD = HERE / "bridge_strict_alpha_entropy_robin_hood_certificate_report.md"

CANONICAL_LEDGERS = [
    (4, 1, 1, 1, 1),
    (3, 2, 1, 1, 1),
    (2, 2, 2, 1, 1),
]
BALANCED_LEDGER = (2, 2, 2, 1, 1)
TOTAL_BINARY_UNITS = 8
BRANCH_COUNT = 5
STRICT_TARGET_ETA = 9.0 / 5.0
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def canonical_multiset(seq: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(seq, reverse=True))


def shannon_entropy_nats(ledger: tuple[int, ...]) -> float:
    total = sum(ledger)
    return -sum((e / total) * math.log(e / total) for e in ledger)


def labelled_multinomial_count(ledger: tuple[int, ...]) -> int:
    numerator = math.factorial(sum(ledger))
    denominator = math.prod(math.factorial(e) for e in ledger)
    return numerator // denominator


def orbit_size(ledger: tuple[int, ...]) -> int:
    counts = Counter(ledger)
    denominator = math.prod(math.factorial(multiplicity) for multiplicity in counts.values())
    return math.factorial(len(ledger)) // denominator


def eta_from_ledger(ledger: tuple[int, ...]) -> float:
    correction = (float(Fraction(2 ** sum(ledger), 3 ** len(ledger)))) ** (1.0 / len(ledger))
    return math.log(NAD12_SUPPORT_SIZE * correction) / math.log(ALPHA_SCALE)


def robin_hood_step(ledger: tuple[int, ...]) -> dict[str, Any] | None:
    high = max(ledger)
    low = min(ledger)
    if high - low < 2:
        return None
    as_list = list(ledger)
    high_index = as_list.index(high)
    low_index = as_list.index(low)
    before = canonical_multiset(tuple(as_list))
    as_list[high_index] -= 1
    as_list[low_index] += 1
    after = canonical_multiset(tuple(as_list))
    ratio = Fraction(high, low + 1)
    before_weight = labelled_multinomial_count(before)
    after_weight = labelled_multinomial_count(after)
    before_entropy = shannon_entropy_nats(before)
    after_entropy = shannon_entropy_nats(after)
    return {
        "before": list(before),
        "after": list(after),
        "moved_from_exponent_a": high,
        "moved_to_exponent_b": low,
        "exact_labelled_weight_ratio_after_over_before": f"{ratio.numerator}/{ratio.denominator}",
        "ratio_greater_than_one": ratio > 1,
        "labelled_weight_before": before_weight,
        "labelled_weight_after": after_weight,
        "labelled_weight_delta": after_weight - before_weight,
        "entropy_before_nats": before_entropy,
        "entropy_after_nats": after_entropy,
        "entropy_delta_nats": after_entropy - before_entropy,
        "entropy_increases": after_entropy > before_entropy,
    }


def certificate_path(start: tuple[int, ...]) -> list[dict[str, Any]]:
    path: list[dict[str, Any]] = []
    current = canonical_multiset(start)
    while current != BALANCED_LEDGER:
        step = robin_hood_step(current)
        if step is None:
            raise RuntimeError(f"Robin-Hood path stalled before balanced ledger: {current}")
        path.append(step)
        current = tuple(step["after"])
    return path


def ledger_row(ledger: tuple[int, ...]) -> dict[str, Any]:
    return {
        "ledger": list(ledger),
        "labelled_multinomial_count": labelled_multinomial_count(ledger),
        "shannon_entropy_nats": shannon_entropy_nats(ledger),
        "orbit_size": orbit_size(ledger),
        "eta_from_shared_product": eta_from_ledger(ledger),
        "eta_residual_vs_9_5": eta_from_ledger(ledger) - STRICT_TARGET_ETA,
        "is_balanced_terminal": ledger == BALANCED_LEDGER,
    }


def main() -> None:
    entropy_report = load_json(ENTROPY)
    source_entropy_winner = tuple(entropy_report["fixed_labelled_branch_entropy_selector"]["entropy_winner"]["ledger"])
    rows = [ledger_row(ledger) for ledger in CANONICAL_LEDGERS]
    paths = {str(ledger): certificate_path(ledger) for ledger in CANONICAL_LEDGERS if ledger != BALANCED_LEDGER}
    all_steps = [step for path in paths.values() for step in path]

    report = {
        "status": "OPEN_STRICT_ALPHA_ENTROPY_ROBIN_HOOD_CERTIFICATE_NO_STRICT_SELECTOR_THEOREM",
        "result_kind": "SCRATCH_STRICT_ALPHA_ENTROPY_ROBIN_HOOD_CERTIFICATE_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "entropy_selector_discriminator": str(ENTROPY.relative_to(ROOT)),
        },
        "fixed_labelled_weight_law": {
            "weight": "W(e_1,...,e_5)=8!/prod_i(e_i!)",
            "robin_hood_ratio": "W_after/W_before = a/(b+1) for a -> a-1 and b -> b+1",
            "strict_increase_condition": "a >= b+2",
            "derived_here_as_arithmetic_certificate": True,
            "strict_selector_source_derived_here": False,
        },
        "candidate_space_rows": rows,
        "certificate_paths_to_balanced_ledger": paths,
        "certificate_summary": {
            "balanced_terminal": list(BALANCED_LEDGER),
            "source_entropy_winner": list(source_entropy_winner),
            "source_entropy_winner_matches_terminal": source_entropy_winner == BALANCED_LEDGER,
            "all_nonterminal_ledgers_have_path": len(paths) == len(CANONICAL_LEDGERS) - 1,
            "all_steps_increase_labelled_weight": all(step["ratio_greater_than_one"] and step["labelled_weight_delta"] > 0 for step in all_steps),
            "all_steps_increase_entropy": all(step["entropy_increases"] and step["entropy_delta_nats"] > 0 for step in all_steps),
            "terminal_labelled_weight": labelled_multinomial_count(BALANCED_LEDGER),
            "terminal_entropy_nats": shannon_entropy_nats(BALANCED_LEDGER),
        },
        "candidate_interpretation": {
            "supported_by_this_probe": bool(
                source_entropy_winner == BALANCED_LEDGER
                and len(paths) == 2
                and all(step["ratio_greater_than_one"] for step in all_steps)
                and all(step["entropy_increases"] for step in all_steps)
            ),
            "content": "Fixed-labelled entropy/multinomial selection constructively forces the balanced ledger by monotone Robin-Hood moves.",
            "why_this_is_more_proof_like": "It replaces a bare argmax table with exact multiplicity-ratio certificates W_after/W_before=a/(b+1) and monotone paths from every nonterminal canonical ledger.",
            "why_this_is_not_enough": "The fixed-labelled entropy selector itself is still not derived from strict nadsoliton geometry and does not discharge QW-2191.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "No theorem derives fixed-labelled branch entropy as the strict-core selector.",
            "This certificate is conditional on fixed-labelled branch multiplicity, not a proof that branches are fixed labelled strict channels.",
            "No theorem derives eta=9/5 without adopting the branch count, ternary normalization, and fixed-labelled entropy selector as extra premises.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Seek an internal strict source for fixed-labelled branch multiplicity; otherwise keep the Robin-Hood path as a conditional selector certificate only.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha entropy Robin-Hood certificate probe\n\n"
        "Status: constructive fixed-labelled entropy certificate for eta=9/5; no strict selector theorem.\n\n"
        "- Weight law: `W(e)=8!/prod_i(e_i!)`; a balancing move `a -> a-1`, `b -> b+1` has exact ratio `a/(b+1)`.\n"
        "- Certificate path: `(4,1,1,1,1) -> (3,2,1,1,1) -> (2,2,2,1,1)`, with labelled weights `1680 -> 3360 -> 5040`.\n"
        f"- Terminal ledger: `{list(BALANCED_LEDGER)}`, entropy `{shannon_entropy_nats(BALANCED_LEDGER):.12f}`, eta residual `{eta_from_ledger(BALANCED_LEDGER)-STRICT_TARGET_ETA:.3e}`.\n"
        "- Honest read: fixed-labelled multiplicity constructively forces the balanced ledger, but fixed-labelled strict channels are still a missing selector premise.\n"
        "- No false pass: no fixed-labelled strict selector theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
