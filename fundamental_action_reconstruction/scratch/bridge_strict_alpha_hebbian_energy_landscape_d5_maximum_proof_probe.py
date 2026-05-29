#!/usr/bin/env python3
"""Scratch probe: exact Hebbian energy landscape proof audit for d5 maximum.

This is the next finite/proof-oriented step after the energy-refined kWTA
selector probe.  The previous probe used learned energy E(S)=x_S^T W x_S as a
conditional tie-refinement rule.  This probe checks whether that learned energy
has an exact finite landscape theorem on all five-node supports of Z_12.

Result:
- E has exactly 12 global maximizers among the 792 supports.
- Those maximizers are exactly the translated fifth-step/d5 supports.
- The energy can be replayed from a six-bin cyclic distance histogram with exact
  rational weights, so no floating-point spectral numerics are involved.
- Energy-refined kWTA has strict positive E-ascent from every non-d5 state on
  every retained branch, with minimum exact increase 4, and therefore reaches
  the unique d5 global-maximum family.

No false pass: this proves only a finite conditional landscape theorem after the
d5 teacher/self-record trace and Hebbian energy rule are supplied.  It does not
derive that trace or rule from strict nadsoliton geometry and does not discharge
QW-2191.
"""
from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_energy_landscape_d5_maximum_proof_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_energy_landscape_d5_maximum_proof_report.md"

N = 12
ACTIVE_COUNT = 5
TARGET_MODE = 5
RHO = Fraction(ACTIVE_COUNT, N)
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"

Support = tuple[int, ...]
DistanceHistogram = tuple[int, ...]
Matrix = list[list[Fraction]]


def all_supports() -> list[Support]:
    return list(combinations(range(N), ACTIVE_COUNT))


def canonical_support(nodes: list[int] | tuple[int, ...]) -> Support:
    return tuple(sorted(int(node) % N for node in nodes))


def d5_teacher_orbit() -> list[Support]:
    return sorted(
        {
            canonical_support([(start + TARGET_MODE * index) % N for index in range(ACTIVE_COUNT)])
            for start in range(N)
        }
    )


def centered_activity(support: Support) -> list[Fraction]:
    active = set(support)
    return [Fraction(1) - RHO if node in active else -RHO for node in range(N)]


def build_centered_zero_self_hebbian(teacher: list[Support]) -> Matrix:
    weights = [[Fraction(0) for _col in range(N)] for _row in range(N)]
    for support in teacher:
        activity = centered_activity(support)
        for row in range(N):
            for col in range(N):
                weights[row][col] += activity[row] * activity[col]
    for index in range(N):
        weights[index][index] = Fraction(0)
    return weights


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def distance_histogram(support: Support) -> DistanceHistogram:
    counts = Counter()
    for left, right in combinations(support, 2):
        distance = (right - left) % N
        folded = min(distance, N - distance)
        counts[folded] += 1
    return tuple(counts[distance] for distance in range(1, N // 2 + 1))


def learned_energy(support: Support, weights: Matrix) -> Fraction:
    return sum(weights[row][col] for row in support for col in support)


def histogram_energy(histogram: DistanceHistogram, distance_weights: dict[int, Fraction]) -> Fraction:
    return 2 * sum(histogram[distance - 1] * distance_weights[distance] for distance in range(1, N // 2 + 1))


def score_nodes(support: Support, weights: Matrix) -> list[tuple[Fraction, int]]:
    return [(sum(weights[node][active] for active in support), node) for node in range(N)]


def raw_kwta_branches(support: Support, weights: Matrix) -> set[Support]:
    scores = score_nodes(support, weights)
    values = sorted({score for score, _node in scores}, reverse=True)
    selected: list[int] = []
    remaining = ACTIVE_COUNT
    for value in values:
        group = [node for score, node in scores if score == value]
        if len(group) < remaining:
            selected.extend(group)
            remaining -= len(group)
        elif len(group) == remaining:
            return {canonical_support(selected + group)}
        else:
            return {canonical_support(selected + list(choice)) for choice in combinations(group, remaining)}
    return {canonical_support(selected)}


def energy_refined_branches(support: Support, weights: Matrix) -> set[Support]:
    candidates = raw_kwta_branches(support, weights)
    best_energy = max(learned_energy(candidate, weights) for candidate in candidates)
    return {candidate for candidate in candidates if learned_energy(candidate, weights) == best_energy}


def energy_landscape_audit(supports: list[Support], teacher: list[Support], weights: Matrix) -> dict[str, Any]:
    energy_by_support = {support: learned_energy(support, weights) for support in supports}
    maximum = max(energy_by_support.values())
    maximizers = sorted(support for support, energy in energy_by_support.items() if energy == maximum)
    energy_histogram = Counter(energy_by_support.values())
    d5_set = set(teacher)
    return {
        "support_count": len(supports),
        "energy_level_count": len(energy_histogram),
        "energy_level_histogram": {fraction_text(energy): count for energy, count in sorted(energy_histogram.items())},
        "maximum_energy": fraction_text(maximum),
        "maximizer_count": len(maximizers),
        "maximizers": [list(row) for row in maximizers],
        "maximizers_equal_d5_teacher_orbit": set(maximizers) == d5_set,
        "all_d5_teacher_patterns_are_global_maxima": all(energy_by_support[support] == maximum for support in teacher),
        "non_d5_global_maximizer_count": sum(1 for support in maximizers if support not in d5_set),
    }


def histogram_proof_audit(supports: list[Support], weights: Matrix) -> dict[str, Any]:
    distance_weights = {distance: weights[0][distance] for distance in range(1, N // 2 + 1)}
    histogram_counter = Counter(distance_histogram(support) for support in supports)
    histogram_rows = []
    replay_failures = []
    for histogram, count in sorted(histogram_counter.items(), key=lambda item: (histogram_energy(item[0], distance_weights), item[0]), reverse=True):
        energy = histogram_energy(histogram, distance_weights)
        histogram_rows.append(
            {
                "distance_histogram_d1_to_d6": list(histogram),
                "support_count": count,
                "energy": fraction_text(energy),
            }
        )
    for support in supports:
        direct = learned_energy(support, weights)
        replay = histogram_energy(distance_histogram(support), distance_weights)
        if direct != replay:
            replay_failures.append(
                {
                    "support": list(support),
                    "direct_energy": fraction_text(direct),
                    "histogram_energy": fraction_text(replay),
                }
            )
    top_row = histogram_rows[0]
    return {
        "distance_weights_d1_to_d6": {str(distance): fraction_text(weight) for distance, weight in distance_weights.items()},
        "histogram_class_count": len(histogram_counter),
        "histogram_rows_sorted_by_energy": histogram_rows,
        "top_histogram_unique": sum(1 for row in histogram_rows if row["energy"] == top_row["energy"]) == 1,
        "top_histogram": top_row,
        "direct_energy_equals_histogram_replay_for_all_supports": not replay_failures,
        "replay_failure_count": len(replay_failures),
        "replay_failure_examples": replay_failures[:10],
    }


def energy_refined_ascent_audit(supports: list[Support], teacher: list[Support], weights: Matrix) -> dict[str, Any]:
    d5_set = set(teacher)
    delta_counter = Counter()
    transition_counter = Counter()
    violations = []
    for support in supports:
        source_energy = learned_energy(support, weights)
        for nxt in energy_refined_branches(support, weights):
            next_energy = learned_energy(nxt, weights)
            delta = next_energy - source_energy
            delta_counter[delta] += 1
            transition_counter[(source_energy, next_energy)] += 1
            if support not in d5_set and delta <= 0:
                violations.append(
                    {
                        "support": list(support),
                        "next": list(nxt),
                        "delta": fraction_text(delta),
                        "source_energy": fraction_text(source_energy),
                        "next_energy": fraction_text(next_energy),
                    }
                )
    non_d5_deltas = [delta for delta in delta_counter for _ in range(delta_counter[delta]) if delta > 0]
    return {
        "strict_positive_ascent_for_every_non_d5_retained_branch": not violations,
        "minimum_positive_delta": fraction_text(min(non_d5_deltas)),
        "delta_histogram": {fraction_text(delta): count for delta, count in sorted(delta_counter.items())},
        "energy_transition_histogram": {
            f"{fraction_text(source)} -> {fraction_text(target)}": count
            for (source, target), count in sorted(transition_counter.items())
        },
        "violation_count": len(violations),
        "violation_examples": violations[:10],
    }


def main() -> None:
    supports = all_supports()
    teacher = d5_teacher_orbit()
    weights = build_centered_zero_self_hebbian(teacher)
    landscape = energy_landscape_audit(supports, teacher, weights)
    histogram_proof = histogram_proof_audit(supports, weights)
    ascent = energy_refined_ascent_audit(supports, teacher, weights)

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_ENERGY_LANDSCAPE_D5_MAXIMUM_PROOF_PROBE__NOT_A_THEOREM",
        "status": "finite-conditional-d5-global-maximum-proof-not-strict-source",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "teacher_record": "12 translated fifth-step/d5 supports",
            "learning_rule": "centered zero-self Hebbian W from d5 teacher orbit",
            "energy": "E(S)=sum_{i,j in S} W_ij = x_S^T W x_S, diag(W)=0",
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "exact_weight_replay": {
            "first_row_by_cyclic_distance": [fraction_text(value) for value in weights[0]],
            "circulant_check": all(weights[row][col] == weights[0][(col - row) % N] for row in range(N) for col in range(N)),
            "diagonal_zero_check": all(weights[index][index] == 0 for index in range(N)),
        },
        "finite_energy_landscape_certificate": landscape,
        "distance_histogram_proof_certificate": histogram_proof,
        "energy_refined_kwta_ascent_replay": ascent,
        "proof_reading": {
            "finite_theorem": "In the supplied d5-teacher centered-Hebbian finite model, the only global maximizers of learned energy over five-node supports are the 12 translated d5 supports.",
            "how_verified": "Exact Fraction arithmetic, complete enumeration of 792 supports, and independent replay through 35 cyclic distance histograms.",
            "relation_to_previous_probe": "The energy-refined kWTA selector now has a finite Lyapunov target: retained non-d5 branches strictly increase E until the unique d5 global-maximum family is reached.",
            "remaining_gap": "The d5 teacher trace and the energy-refined selector premise are still assumed, not derived from strict nadsoliton geometry.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may contain finite self-recorded resonance patterns whose learned energy is an internal consistency functional.",
            "forbidden_reading": "The energy landscape is not introduced as a separate informational layer underneath the nadsoliton.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives the d5 teacher/self-record trace from strict nadsoliton geometry.",
            "No theorem derives the centered zero-self Hebbian energy or energy-refined kWTA selector from strict geometry.",
            "This is a finite conditional landscape theorem after premises are supplied, not strict-core selector closure.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, or fifth-mode source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Try to derive the d5 teacher/self-record trace or the energy functional from strict nadsoliton structure; otherwise treat this as a conditional finite theorem supporting, but not proving, a strict selector source.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian energy-landscape d5 maximum proof probe\n\n"
        "Status: finite conditional d5 global-maximum proof; not a strict source theorem.\n\n"
        f"- Supports scanned: `{len(supports)}` five-active-node states on `Z_12`.\n"
        f"- Energy levels: `{landscape['energy_level_count']}` with maximum `{landscape['maximum_energy']}`.\n"
        f"- Global maximizer count: `{landscape['maximizer_count']}`; equals d5 orbit: `{landscape['maximizers_equal_d5_teacher_orbit']}`.\n"
        f"- Distance-histogram classes: `{histogram_proof['histogram_class_count']}`; top histogram: `{histogram_proof['top_histogram']}`.\n"
        f"- Direct energy equals histogram replay for all supports: `{histogram_proof['direct_energy_equals_histogram_replay_for_all_supports']}`.\n"
        f"- Energy-refined kWTA minimum positive non-d5 ascent: `{ascent['minimum_positive_delta']}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: d5 is the exact finite global maximum of supplied Hebbian energy, but the d5 teacher trace/energy rule remain premises.\n"
        "- No false pass: no d5-source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
