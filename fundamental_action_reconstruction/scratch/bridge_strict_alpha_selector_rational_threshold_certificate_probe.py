#!/usr/bin/env python3
"""Scratch probe: rational orbit-weight threshold certificate for eta selector.

The orbit-weight threshold probe found the exact switch point

    gamma_c = log(3/2) / log(2)

for the selector family Score_gamma=log W + gamma*log O.  This probe removes
floating-point ambiguity for rational selector exponents gamma=p/q.

For B=(2,2,2,1,1) and C=(3,2,1,1,1), B beats C iff

    log(3/2) - gamma*log(2) > 0.

For gamma=p/q this is exactly equivalent to the integer inequality

    3^q > 2^(p+q).

There is no rational tie: 3^q = 2^(p+q) is impossible for q>0 by prime
factorization.  The balanced ledger is therefore safe for a rational p/q exactly
when 3^q - 2^(p+q) is positive.

This is still not a strict selector theorem.  It gives a machine-checkable
rational certificate for any proposed orbit-weight exponent, and sharpens the
proof obligation to deriving a rational gamma whose integer certificate is on
the safe side, or deriving gamma=0/fixed-labelled channels directly.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
THRESHOLD = HERE / "bridge_strict_alpha_selector_orbit_weight_threshold_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_selector_rational_threshold_certificate_report.json"
OUT_MD = HERE / "bridge_strict_alpha_selector_rational_threshold_certificate_report.md"

MAX_DENOMINATOR = 64
BALANCED_LEDGER = [2, 2, 2, 1, 1]
ORBIT_FAVOURED_LEDGER = [3, 2, 1, 1, 1]
COMMON_RATIONAL_GAMMAS = [Fraction(0, 1), Fraction(1, 2), Fraction(7, 12), Fraction(31, 53), Fraction(24, 41), Fraction(3, 5), Fraction(1, 1)]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def gamma_critical() -> float:
    return math.log(Fraction(3, 2)) / math.log(2)


def integer_delta_for_fraction(gamma: Fraction) -> int:
    return 3 ** gamma.denominator - 2 ** (gamma.numerator + gamma.denominator)


def classification_for_delta(delta: int) -> str:
    if delta > 0:
        return "SAFE_BALANCED_WINS_BELOW_THRESHOLD"
    if delta < 0:
        return "FAIL_ORBIT_FAVOURED_WINS_ABOVE_THRESHOLD"
    return "IMPOSSIBLE_RATIONAL_TIE"


def row_for_fraction(gamma: Fraction) -> dict[str, Any]:
    delta = integer_delta_for_fraction(gamma)
    return {
        "gamma_label": f"{gamma.numerator}/{gamma.denominator}",
        "gamma_decimal": float(gamma),
        "p": gamma.numerator,
        "q": gamma.denominator,
        "integer_certificate_delta_3q_minus_2p_plus_q": delta,
        "classification": classification_for_delta(delta),
        "winner": BALANCED_LEDGER if delta > 0 else ORBIT_FAVOURED_LEDGER,
        "distance_to_gamma_c_decimal": float(gamma) - gamma_critical(),
    }


def reduced_rationals_between_zero_and_one(max_denominator: int) -> list[Fraction]:
    values = set()
    for q in range(1, max_denominator + 1):
        for p in range(0, q + 1):
            values.add(Fraction(p, q))
    return sorted(values)


def main() -> None:
    threshold_report = load_json(THRESHOLD)
    gamma_c = gamma_critical()
    rationals = reduced_rationals_between_zero_and_one(MAX_DENOMINATOR)
    rows = [row_for_fraction(gamma) for gamma in rationals]
    safe_rows = [row for row in rows if row["classification"] == "SAFE_BALANCED_WINS_BELOW_THRESHOLD"]
    fail_rows = [row for row in rows if row["classification"] == "FAIL_ORBIT_FAVOURED_WINS_ABOVE_THRESHOLD"]
    tie_rows = [row for row in rows if row["classification"] == "IMPOSSIBLE_RATIONAL_TIE"]
    nearest_rows = sorted(rows, key=lambda row: abs(row["distance_to_gamma_c_decimal"]))[:16]
    best_safe = max(safe_rows, key=lambda row: row["gamma_decimal"])
    first_fail = min(fail_rows, key=lambda row: row["gamma_decimal"])
    common_rows = [row_for_fraction(gamma) for gamma in COMMON_RATIONAL_GAMMAS]

    report = {
        "status": "OPEN_STRICT_ALPHA_SELECTOR_RATIONAL_THRESHOLD_CERTIFICATE_NO_STRICT_SELECTOR_THEOREM",
        "result_kind": "SCRATCH_STRICT_ALPHA_SELECTOR_RATIONAL_THRESHOLD_CERTIFICATE_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "selector_orbit_weight_threshold": str(THRESHOLD.relative_to(ROOT)),
        },
        "exact_rational_certificate_rule": {
            "gamma_c_closed_form": "log(3/2)/log(2)",
            "gamma_c_numeric": gamma_c,
            "rational_gamma": "gamma=p/q in lowest terms, q>0",
            "balanced_wins_iff": "3^q > 2^(p+q)",
            "orbit_favoured_wins_iff": "3^q < 2^(p+q)",
            "rational_tie_impossible_by_prime_factorization": True,
            "tie_impossibility_reason": "3^q=2^(p+q) would equate a pure power of 3 with a pure power of 2 for q>0.",
        },
        "bounded_rational_scan": {
            "max_denominator": MAX_DENOMINATOR,
            "reduced_rationals_scanned_count": len(rows),
            "safe_count": len(safe_rows),
            "fail_count": len(fail_rows),
            "tie_count": len(tie_rows),
            "best_safe_below_threshold": best_safe,
            "first_fail_above_threshold": first_fail,
            "nearest_rows_by_decimal_distance": nearest_rows,
        },
        "common_rational_gamma_certificates": common_rows,
        "threshold_report_consistency": {
            "source_gamma_c_numeric": threshold_report["critical_threshold"]["gamma_critical_numeric"],
            "local_gamma_c_numeric": gamma_c,
            "source_and_local_gamma_c_match": abs(threshold_report["critical_threshold"]["gamma_critical_numeric"] - gamma_c) < 1e-15,
            "source_below_winner": threshold_report["critical_threshold"]["winner_just_below_threshold"],
            "source_above_winner": threshold_report["critical_threshold"]["winner_just_above_threshold"],
        },
        "candidate_interpretation": {
            "supported_by_this_probe": bool(
                not tie_rows
                and best_safe["gamma_label"] == "31/53"
                and first_fail["gamma_label"] == "24/41"
                and row_for_fraction(Fraction(1, 2))["classification"] == "SAFE_BALANCED_WINS_BELOW_THRESHOLD"
                and row_for_fraction(Fraction(3, 5))["classification"] == "FAIL_ORBIT_FAVOURED_WINS_ABOVE_THRESHOLD"
            ),
            "content": "Every rational orbit-weight exponent has an exact integer certificate; below-threshold rationals keep the balanced eta-ledger, above-threshold rationals flip to the orbit-favoured ledger.",
            "why_this_is_more_proof_like": "It replaces floating threshold comparison by the integer sign of 3^q-2^(p+q), with no possible rational tie.",
            "why_this_is_not_enough": "No strict theorem supplies a safe rational gamma, gamma=0, or fixed-labelled channels from nadsoliton geometry.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "No theorem derives fixed-labelled branch entropy as the strict-core selector.",
            "No theorem derives any safe rational orbit-weight exponent p/q with 3^q > 2^(p+q).",
            "Rational exponents p/q with 3^q < 2^(p+q) select (3,2,1,1,1), not the balanced eta-ledger.",
            "No theorem derives eta=9/5 without adopting the branch count, ternary normalization, and a certified below-threshold selector convention as extra premises.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "If a candidate rational orbit-weight exponent is proposed, require the integer certificate 3^q > 2^(p+q); otherwise keep the eta-ledger selector-sensitive.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha selector rational threshold certificate probe\n\n"
        "Status: exact rational selector certificate for eta=9/5; no strict selector theorem.\n\n"
        "- Rational rule: for `gamma=p/q`, balanced wins iff `3^q > 2^(p+q)`; rational tie is impossible by prime factorization.\n"
        f"- Bounded scan: `{len(rows)}` reduced rationals with denominator <= `{MAX_DENOMINATOR}`; safe `{len(safe_rows)}`, fail `{len(fail_rows)}`, tie `{len(tie_rows)}`.\n"
        f"- Best safe in scan: `{best_safe['gamma_label']}`; first fail in scan: `{first_fail['gamma_label']}`; `gamma_c={gamma_c:.12f}`.\n"
        "- Common certificates: `1/2` is safe; `3/5` is already above threshold and flips to `(3,2,1,1,1)`.\n"
        "- No false pass: no strict derivation of a safe rational gamma, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
