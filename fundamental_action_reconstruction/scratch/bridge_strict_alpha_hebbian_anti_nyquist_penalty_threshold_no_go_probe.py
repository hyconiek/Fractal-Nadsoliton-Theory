#!/usr/bin/env python3
"""Scratch probe: exact anti-Nyquist penalty threshold and remaining d5 no-go.

The previous symbolic mode-competition probe showed that unconstrained highest
coherent power selects the k=6 Nyquist/parity channel instead of d5.  This probe
checks the next honest repair attempt: add a scalar anti-Nyquist penalty

    S_mu(k, H) = P_k(H) - mu * P_6(H)

on the 35 feasible five-node pair-distance histograms H of Z_12.

Finite exact result in Q(sqrt(3)):
- mu = 0 reproduces the Nyquist obstruction.
- The exact threshold at which the d5/fifth and contiguous/first unit classes
  catch the Nyquist winner is mu_* = 3/4 - sqrt(3)/6.
- For every mu > mu_*, the scalar penalty removes the Nyquist winner but leaves
  an unavoidable two-way unit mirror tie: k=1 contiguous and k=5 d5.
- Therefore an anti-Nyquist scalar penalty can be a conditional filter, but it
  still does not derive a unique d5 selector without an extra unit-orientation
  or fifth-channel premise.

No false pass: this is not a strict-core selector theorem, not a derivation of a
Hebbian learning law, and not a QW-2191 discharge.
"""
from __future__ import annotations

import json
from collections import Counter
from functools import cmp_to_key
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_anti_nyquist_penalty_threshold_no_go_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_anti_nyquist_penalty_threshold_no_go_report.md"

N = 12
ACTIVE_COUNT = 5
AUT_UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
D5_HISTOGRAM = (0, 3, 2, 1, 4, 0)
CONTIGUOUS_HISTOGRAM = (4, 3, 2, 1, 0, 0)
NYQUIST_HISTOGRAM = (0, 4, 0, 4, 0, 2)
CRITICAL_MU = (Fraction(3, 4), Fraction(-1, 6))
COS_COEFFICIENTS_BY_MODE: dict[int, tuple[tuple[Fraction, Fraction], ...]] = {
    1: ((0, Fraction(1, 2)), (Fraction(1, 2), 0), (0, 0), (Fraction(-1, 2), 0), (0, Fraction(-1, 2)), (-1, 0)),
    2: ((Fraction(1, 2), 0), (Fraction(-1, 2), 0), (-1, 0), (Fraction(-1, 2), 0), (Fraction(1, 2), 0), (1, 0)),
    3: ((0, 0), (-1, 0), (0, 0), (1, 0), (0, 0), (-1, 0)),
    4: ((Fraction(-1, 2), 0), (Fraction(-1, 2), 0), (1, 0), (Fraction(-1, 2), 0), (Fraction(-1, 2), 0), (1, 0)),
    5: ((0, Fraction(-1, 2)), (Fraction(1, 2), 0), (0, 0), (Fraction(-1, 2), 0), (0, Fraction(1, 2)), (-1, 0)),
    6: ((-1, 0), (1, 0), (-1, 0), (1, 0), (-1, 0), (1, 0)),
}

Support = tuple[int, ...]
Histogram = tuple[int, int, int, int, int, int]
Quadratic = tuple[Fraction, Fraction]


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def quadratic_text(value: Quadratic) -> str:
    rational, sqrt3_coeff = value
    if sqrt3_coeff == 0:
        return fraction_text(rational)
    if rational == 0:
        if sqrt3_coeff == 1:
            return "sqrt3"
        if sqrt3_coeff == -1:
            return "-sqrt3"
        return f"{fraction_text(sqrt3_coeff)}*sqrt3"
    sign = "+" if sqrt3_coeff > 0 else "-"
    coeff_abs = abs(sqrt3_coeff)
    coeff_text = "sqrt3" if coeff_abs == 1 else f"{fraction_text(coeff_abs)}*sqrt3"
    return f"{fraction_text(rational)} {sign} {coeff_text}"


def quadratic_json(value: Quadratic) -> dict[str, str]:
    return {
        "rational_part": fraction_text(value[0]),
        "sqrt3_coefficient": fraction_text(value[1]),
        "expression": quadratic_text(value),
    }


def sign_quadratic(value: Quadratic) -> int:
    rational, sqrt3_coeff = value
    if sqrt3_coeff == 0:
        return (rational > 0) - (rational < 0)
    if rational == 0:
        return (sqrt3_coeff > 0) - (sqrt3_coeff < 0)
    if rational > 0 and sqrt3_coeff > 0:
        return 1
    if rational < 0 and sqrt3_coeff < 0:
        return -1
    square_difference = 3 * sqrt3_coeff * sqrt3_coeff - rational * rational
    if sqrt3_coeff > 0 and rational < 0:
        return (square_difference > 0) - (square_difference < 0)
    return (square_difference < 0) - (square_difference > 0)


def compare_quadratic(left: Quadratic, right: Quadratic) -> int:
    return sign_quadratic((left[0] - right[0], left[1] - right[1]))


def all_supports() -> list[Support]:
    return [tuple(support) for support in combinations(range(N), ACTIVE_COUNT)]


def folded_distance(left: int, right: int) -> int:
    distance = (right - left) % N
    return min(distance, N - distance)


def distance_histogram(support: Support) -> Histogram:
    counts = [0] * (N // 2)
    for left, right in combinations(support, 2):
        counts[folded_distance(left, right) - 1] += 1
    return tuple(counts)  # type: ignore[return-value]


def mode_power_from_histogram(mode: int, histogram: Histogram) -> Quadratic:
    rational = Fraction(ACTIVE_COUNT)
    sqrt3_coeff = Fraction(0)
    for count, (cos_rational, cos_sqrt3_coeff) in zip(histogram, COS_COEFFICIENTS_BY_MODE[mode]):
        rational += 2 * count * cos_rational
        sqrt3_coeff += 2 * count * cos_sqrt3_coeff
    return rational, sqrt3_coeff


def nyquist_power(histogram: Histogram) -> Fraction:
    return mode_power_from_histogram(6, histogram)[0]


def penalized_score(mode: int, histogram: Histogram, mu: Quadratic) -> Quadratic:
    power = mode_power_from_histogram(mode, histogram)
    p6 = nyquist_power(histogram)
    return power[0] - p6 * mu[0], power[1] - p6 * mu[1]


def candidate_rows(histogram_counter: Counter[Histogram], mu: Quadratic) -> list[dict[str, Any]]:
    rows = []
    for histogram, support_count in histogram_counter.items():
        for mode in range(1, 7):
            rows.append(
                {
                    "mode": mode,
                    "distance_histogram_d1_to_d6": list(histogram),
                    "support_count": support_count,
                    "power_exact": quadratic_json(mode_power_from_histogram(mode, histogram)),
                    "nyquist_power_exact": fraction_text(nyquist_power(histogram)),
                    "score_exact": quadratic_json(penalized_score(mode, histogram, mu)),
                    "is_d5_fifth": mode == 5 and histogram == D5_HISTOGRAM,
                    "is_contiguous_first": mode == 1 and histogram == CONTIGUOUS_HISTOGRAM,
                    "is_nyquist_sixth": mode == 6 and histogram == NYQUIST_HISTOGRAM,
                }
            )
    return sorted(
        rows,
        key=cmp_to_key(
            lambda left, right: -compare_quadratic(
                (Fraction(left["score_exact"]["rational_part"]), Fraction(left["score_exact"]["sqrt3_coefficient"])),
                (Fraction(right["score_exact"]["rational_part"]), Fraction(right["score_exact"]["sqrt3_coefficient"])),
            )
        ),
    )


def winners_for_mu(histogram_counter: Counter[Histogram], mu: Quadratic) -> dict[str, Any]:
    rows = candidate_rows(histogram_counter, mu)
    top_score = (Fraction(rows[0]["score_exact"]["rational_part"]), Fraction(rows[0]["score_exact"]["sqrt3_coefficient"]))
    winners = [
        row
        for row in rows
        if compare_quadratic(
            (Fraction(row["score_exact"]["rational_part"]), Fraction(row["score_exact"]["sqrt3_coefficient"])),
            top_score,
        )
        == 0
    ]
    return {
        "mu_exact": quadratic_json(mu),
        "top_score_exact": quadratic_json(top_score),
        "winner_count_over_histogram_mode_candidates": len(winners),
        "winners": winners,
    }


def threshold_certificate(histogram_counter: Counter[Histogram]) -> dict[str, Any]:
    d5_power = mode_power_from_histogram(5, D5_HISTOGRAM)
    d5_p6 = nyquist_power(D5_HISTOGRAM)
    threshold_rows = []
    all_mu_tie_rows = []
    max_threshold: Quadratic | None = None
    max_threshold_rows = []

    for histogram, support_count in histogram_counter.items():
        for mode in range(1, 7):
            if mode == 5 and histogram == D5_HISTOGRAM:
                continue
            power = mode_power_from_histogram(mode, histogram)
            p6 = nyquist_power(histogram)
            numerator = (power[0] - d5_power[0], power[1] - d5_power[1])
            denominator = p6 - d5_p6
            base = {
                "mode": mode,
                "distance_histogram_d1_to_d6": list(histogram),
                "support_count": support_count,
                "power_exact": quadratic_json(power),
                "nyquist_power_exact": fraction_text(p6),
            }
            if denominator == 0:
                if compare_quadratic(numerator, (0, 0)) == 0:
                    all_mu_tie_rows.append(base)
                continue
            if denominator > 0:
                threshold = (numerator[0] / denominator, numerator[1] / denominator)
                row = {**base, "required_mu_ge_exact": quadratic_json(threshold)}
                threshold_rows.append(row)
                if max_threshold is None or compare_quadratic(threshold, max_threshold) > 0:
                    max_threshold = threshold
                    max_threshold_rows = [row]
                elif compare_quadratic(threshold, max_threshold) == 0:
                    max_threshold_rows.append(row)

    assert max_threshold == CRITICAL_MU
    return {
        "d5_reference": {
            "mode": 5,
            "distance_histogram_d1_to_d6": list(D5_HISTOGRAM),
            "power_exact": quadratic_json(d5_power),
            "nyquist_power_exact": fraction_text(d5_p6),
        },
        "critical_mu_exact": quadratic_json(max_threshold),
        "critical_mu_source_rows": max_threshold_rows,
        "all_mu_tie_rows_with_d5": all_mu_tie_rows,
        "all_mu_tie_is_k1_contiguous": all_mu_tie_rows == [
            {
                "mode": 1,
                "distance_histogram_d1_to_d6": list(CONTIGUOUS_HISTOGRAM),
                "support_count": histogram_counter[CONTIGUOUS_HISTOGRAM],
                "power_exact": quadratic_json(mode_power_from_histogram(1, CONTIGUOUS_HISTOGRAM)),
                "nyquist_power_exact": fraction_text(nyquist_power(CONTIGUOUS_HISTOGRAM)),
            }
        ],
        "threshold_candidate_count": len(threshold_rows),
    }


def parity_balance_certificate(supports: list[Support]) -> dict[str, Any]:
    histogram_counter = Counter(distance_histogram(support) for support in supports)
    p6_by_support = Counter(fraction_text(nyquist_power(distance_histogram(support))) for support in supports)
    min_p6 = min(nyquist_power(distance_histogram(support)) for support in supports)
    min_supports = [support for support in supports if nyquist_power(distance_histogram(support)) == min_p6]
    return {
        "p6_support_histogram": dict(sorted(p6_by_support.items(), key=lambda item: Fraction(item[0]))),
        "minimal_p6_exact": fraction_text(min_p6),
        "minimal_p6_support_count": len(min_supports),
        "minimal_p6_histogram_class_count": len({distance_histogram(support) for support in min_supports}),
        "d5_has_minimal_p6": nyquist_power(D5_HISTOGRAM) == min_p6,
        "contiguous_has_minimal_p6": nyquist_power(CONTIGUOUS_HISTOGRAM) == min_p6,
        "nyquist_histogram_support_count": histogram_counter[NYQUIST_HISTOGRAM],
    }


def main() -> None:
    supports = all_supports()
    histogram_counter = Counter(distance_histogram(support) for support in supports)
    threshold = threshold_certificate(histogram_counter)
    sample_mu_winners = {
        "mu_0": winners_for_mu(histogram_counter, (0, 0)),
        "mu_critical": winners_for_mu(histogram_counter, CRITICAL_MU),
        "mu_one_half": winners_for_mu(histogram_counter, (Fraction(1, 2), 0)),
        "mu_1": winners_for_mu(histogram_counter, (1, 0)),
    }

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_ANTI_NYQUIST_PENALTY_THRESHOLD_NO_GO_PROBE__NOT_A_THEOREM",
        "status": "anti-nyquist-penalty-removes-nyquist-but-leaves-k1-k5-unit-tie",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "histogram_class_count": len(histogram_counter),
            "mode_histogram_candidate_count": len(histogram_counter) * 6,
            "modes_scanned": [1, 2, 3, 4, 5, 6],
            "automorphism_units": AUT_UNITS,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "score_definition": {
            "formula": "S_mu(k,H)=P_k(H)-mu*P_6(H)",
            "purpose": "test whether a scalar anti-Nyquist penalty can repair the unconstrained highest-resonance no-go",
            "critical_mu_exact": threshold["critical_mu_exact"],
        },
        "parity_balance_certificate": parity_balance_certificate(supports),
        "threshold_certificate": threshold,
        "sample_mu_winners": sample_mu_winners,
        "selector_readout": {
            "positive_filter": "For mu > mu_* the k=6 Nyquist winner is removed from the top score.",
            "remaining_obstruction": "For every mu, k=5 d5 and k=1 contiguous have identical P_k and P_6 values under this scalar score, so the scalar penalty cannot uniquely select d5.",
            "conditional_completion": "A separate fifth-channel or unit-orientation premise would select d5 from the remaining k=1/k=5 mirror tie, but that premise is not derived here.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may later supply an internal anti-Nyquist or fifth-channel self-record.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to supply the penalty or unit orientation.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives the d5 teacher/self-record trace from strict nadsoliton geometry.",
            "No theorem derives a Hebbian learning law as strict-core dynamics.",
            "No theorem derives anti-Nyquist penalty mu or fifth-mode channel selection as a strict source.",
            "Anti-Nyquist scalar filtering does not uniquely select d5; it leaves the k=1/k=5 unit mirror tie.",
            "This is a finite conditional certificate, not strict-core selector closure.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, unit-orientation source, or fifth-mode source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Try to derive an internal unit-orientation/fifth-channel source; otherwise keep the d5 choice after anti-Nyquist filtering as an explicit non-strict selector premise.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian anti-Nyquist penalty threshold no-go probe\n\n"
        "Status: anti-Nyquist scalar penalty removes k=6 only above a threshold, but leaves a k=1/k=5 tie.\n\n"
        f"- Supports scanned: `{len(supports)}`; histogram classes: `{len(histogram_counter)}`; mode-histogram candidates: `{len(histogram_counter) * 6}`.\n"
        f"- Score: `S_mu(k,H)=P_k(H)-mu*P_6(H)`.\n"
        f"- Critical threshold: `mu_*={quadratic_text(CRITICAL_MU)}`.\n"
        f"- mu=0 winners: `{sample_mu_winners['mu_0']['winners']}`.\n"
        f"- mu=mu_* winners: `{sample_mu_winners['mu_critical']['winners']}`.\n"
        f"- mu=1/2 winners: `{sample_mu_winners['mu_one_half']['winners']}`.\n"
        f"- Parity/minimal P_6 support count: `{report['parity_balance_certificate']['minimal_p6_support_count']}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: anti-Nyquist filtering is useful but not a unique d5 selector; a unit-orientation/fifth-channel premise remains external here.\n"
        "- No false pass: no strict fifth-mode source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
