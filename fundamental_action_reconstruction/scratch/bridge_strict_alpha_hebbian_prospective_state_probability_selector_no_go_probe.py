#!/usr/bin/env python3
"""Scratch probe: prospective-state probability cannot break the unit-mirror tie by itself.

Motivation: after the anti-Nyquist and unit-orientation probes, the remaining
survivor orbit is

    A1 = (k=1, contiguous histogram),
    A5 = (k=5, d5 histogram).

A natural next question is whether a faith-like/prospective datum -- information
inside the nadsoliton about a possible future state of itself -- can act as the
missing selector.  This probe audits the exact finite version of that idea.

No false pass: if the prospective probability law is full Aut(Z_12)-invariant,
it assigns equal probability to the two unit-mirror states and cannot select
singleton d5.  If a prospective law favors A5, then it is exactly the missing
unit-orientation/fifth-axis bit in probabilistic form, not a strict-core
selector derivation.
"""
from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_prospective_state_probability_selector_no_go_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_prospective_state_probability_selector_no_go_report.md"

N = 12
ACTIVE_COUNT = 5
AUT_UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
GRID_DENOMINATOR = 12
K1_CONTIGUOUS = "A1_k1_contiguous"
K5_D5 = "A5_k5_d5"
AXES = [K1_CONTIGUOUS, K5_D5]
D5_HISTOGRAM = [0, 3, 2, 1, 4, 0]
CONTIGUOUS_HISTOGRAM = [4, 3, 2, 1, 0, 0]


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def all_supports() -> list[tuple[int, ...]]:
    return [tuple(support) for support in combinations(range(N), ACTIVE_COUNT)]


def folded(value: int) -> int:
    residue = value % N
    return min(residue, N - residue)


def distance_histogram(support: tuple[int, ...]) -> tuple[int, int, int, int, int, int]:
    counts = [0] * (N // 2)
    for left, right in combinations(support, 2):
        counts[folded(right - left) - 1] += 1
    return tuple(counts)  # type: ignore[return-value]


def unit_action_on_axis(axis: str, unit: int) -> str:
    if axis == K1_CONTIGUOUS:
        return K1_CONTIGUOUS if folded(unit * 1) == 1 else K5_D5
    if axis == K5_D5:
        return K5_D5 if folded(unit * 5) == 5 else K1_CONTIGUOUS
    raise ValueError(axis)


def unit_action_table() -> list[dict[str, Any]]:
    return [
        {
            "unit": unit,
            "action": "preserve" if unit in [1, 11] else "swap",
            "maps_A1_to": unit_action_on_axis(K1_CONTIGUOUS, unit),
            "maps_A5_to": unit_action_on_axis(K5_D5, unit),
        }
        for unit in AUT_UNITS
    ]


def is_invariant_prior(prob_a5: Fraction) -> bool:
    return prob_a5 == Fraction(1, 2)


def prior_row(prob_a5: Fraction) -> dict[str, Any]:
    prob_a1 = 1 - prob_a5
    return {
        "P_A1": fraction_text(prob_a1),
        "P_A5": fraction_text(prob_a5),
        "full_aut_invariant": is_invariant_prior(prob_a5),
        "d5_dominant": prob_a5 > prob_a1,
        "imports_unit_orientation_bit": prob_a5 != prob_a1,
    }


def prior_grid_audit(denominator: int) -> dict[str, Any]:
    rows = [prior_row(Fraction(numerator, denominator)) for numerator in range(denominator + 1)]
    return {
        "denominator": denominator,
        "prior_count": len(rows),
        "full_aut_invariant_prior_count": sum(1 for row in rows if row["full_aut_invariant"]),
        "d5_dominant_prior_count": sum(1 for row in rows if row["d5_dominant"]),
        "full_aut_invariant_d5_dominant_prior_count": sum(
            1 for row in rows if row["full_aut_invariant"] and row["d5_dominant"]
        ),
        "rows": rows,
    }


def deterministic_maps() -> list[dict[str, Any]]:
    rows = []
    for image_a1, image_a5 in product(AXES, repeat=2):
        mapping = {K1_CONTIGUOUS: image_a1, K5_D5: image_a5}
        equivariant = all(
            mapping[unit_action_on_axis(axis, unit)] == unit_action_on_axis(mapping[axis], unit)
            for axis in AXES
            for unit in AUT_UNITS
        )
        rows.append(
            {
                "map": {K1_CONTIGUOUS: image_a1, K5_D5: image_a5},
                "full_aut_equivariant": equivariant,
                "creates_singleton_d5_from_both_orbit": False,
                "interpretation": (
                    "orientation-transport-only" if equivariant and image_a1 != image_a5 else "non-equivariant-or-constant"
                ),
            }
        )
    return rows


def transition_row(a1_to_a5: Fraction, a5_to_a5: Fraction) -> dict[str, Any]:
    # Kernel rows are ordered as source A1/A5 and target A1/A5.
    a1_to_a1 = 1 - a1_to_a5
    a5_to_a1 = 1 - a5_to_a5
    equivariant = a5_to_a5 == a1_to_a1 and a5_to_a1 == a1_to_a5
    uniform_next_a5 = Fraction(1, 2) * a1_to_a5 + Fraction(1, 2) * a5_to_a5
    return {
        "kernel": {
            "from_A1": {"to_A1": fraction_text(a1_to_a1), "to_A5": fraction_text(a1_to_a5)},
            "from_A5": {"to_A1": fraction_text(a5_to_a1), "to_A5": fraction_text(a5_to_a5)},
        },
        "full_aut_equivariant": equivariant,
        "next_P_A5_from_uniform_prior": fraction_text(uniform_next_a5),
        "d5_dominant_from_uniform_prior": uniform_next_a5 > Fraction(1, 2),
        "imports_unit_orientation_bit_if_d5_dominant": uniform_next_a5 > Fraction(1, 2) and not equivariant,
    }


def transition_grid_audit(denominator: int) -> dict[str, Any]:
    rows = [
        transition_row(Fraction(a1_to_a5, denominator), Fraction(a5_to_a5, denominator))
        for a1_to_a5 in range(denominator + 1)
        for a5_to_a5 in range(denominator + 1)
    ]
    equivariant_rows = [row for row in rows if row["full_aut_equivariant"]]
    d5_dominant_rows = [row for row in rows if row["d5_dominant_from_uniform_prior"]]
    return {
        "denominator": denominator,
        "kernel_count": len(rows),
        "full_aut_equivariant_kernel_count": len(equivariant_rows),
        "d5_dominant_from_uniform_prior_count": len(d5_dominant_rows),
        "full_aut_equivariant_d5_dominant_from_uniform_prior_count": sum(
            1 for row in equivariant_rows if row["d5_dominant_from_uniform_prior"]
        ),
        "full_aut_equivariant_family": "[[a, 1-a], [1-a, a]] in the A1/A5 basis; uniform prior stays uniform.",
        "sample_equivariant_kernels": [equivariant_rows[index] for index in [0, len(equivariant_rows) // 2, -1]],
        "sample_non_equivariant_d5_dominant_kernels": d5_dominant_rows[:3],
    }


def future_predicate_rows() -> list[dict[str, Any]]:
    rows = []
    for mask in range(1 << len(AXES)):
        selected = {axis for index, axis in enumerate(AXES) if mask & (1 << index)}
        invariant = all(unit_action_on_axis(axis, unit) in selected for axis in selected for unit in AUT_UNITS)
        rows.append(
            {
                "selected_future_axes": sorted(selected),
                "name": "none" if not selected else "both" if selected == set(AXES) else sorted(selected)[0],
                "full_aut_invariant": invariant,
                "is_singleton_d5_future_predicate": selected == {K5_D5},
            }
        )
    return rows


def main() -> None:
    supports = all_supports()
    histogram_counter = Counter(distance_histogram(support) for support in supports)
    det_rows = deterministic_maps()
    predicate_rows = future_predicate_rows()
    prior_audit = prior_grid_audit(GRID_DENOMINATOR)
    transition_audit = transition_grid_audit(GRID_DENOMINATOR)

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_PROSPECTIVE_STATE_PROBABILITY_SELECTOR_NO_GO_PROBE__NOT_A_THEOREM",
        "status": "aut-invariant-prospective-state-probability-cannot-break-unit-mirror-tie",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "histogram_class_count": len(histogram_counter),
            "automorphism_units": AUT_UNITS,
            "survivor_axes": [
                {"name": K1_CONTIGUOUS, "mode": 1, "distance_histogram_d1_to_d6": CONTIGUOUS_HISTOGRAM},
                {"name": K5_D5, "mode": 5, "distance_histogram_d1_to_d6": D5_HISTOGRAM},
            ],
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "unit_action_certificate": {
            "action_table": unit_action_table(),
            "orbit_size": 2,
            "swapping_units": [5, 7],
            "preserving_units": [1, 11],
        },
        "future_predicate_enumeration": {
            "predicate_count": len(predicate_rows),
            "rows": predicate_rows,
            "full_aut_invariant_predicate_names": [row["name"] for row in predicate_rows if row["full_aut_invariant"]],
            "full_aut_singleton_d5_future_predicate_exists": any(
                row["full_aut_invariant"] and row["is_singleton_d5_future_predicate"] for row in predicate_rows
            ),
        },
        "prospective_prior_grid_audit": prior_audit,
        "prospective_transition_grid_audit": transition_audit,
        "deterministic_future_map_enumeration": {
            "map_count": len(det_rows),
            "full_aut_equivariant_map_count": sum(1 for row in det_rows if row["full_aut_equivariant"]),
            "rows": det_rows,
            "equivariant_maps_do_not_create_new_orientation_from_unoriented_orbit": True,
        },
        "selector_interpretation": {
            "faith_like_hypothesis_tested": "A prospective self-record/probability of a possible future nadsoliton state is allowed as an internal nadsoliton datum, not as a layer below the nadsoliton.",
            "negative_result": "If that prospective datum is full-Aut(Z_12)-invariant, A1 and A5 receive equal probability and no singleton d5 selector appears.",
            "conditional_positive_result": "A future prior or transition favoring A5 can select d5 only by breaking the unit mirror; this is the same missing one-bit unit-axis record in probabilistic language.",
            "highest_resonance_reading": "A highest-future-resonance rule is harmless if it is built from invariant data, because it remains tied on the unit-mirror orbit; if it favors k=5/d5, the favoritism is an extra selector premise.",
            "honest_limit": "This probe does not derive the prospective probability law, the Hebbian law, or the unit-axis bit from strict nadsoliton geometry.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself is informational and may consistently contain a prospective self-record of its own possible states.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to store future probabilities.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives a faith-like prospective self-record from strict nadsoliton geometry.",
            "No theorem derives a Hebbian learning law as strict-core dynamics.",
            "No theorem derives the required one-bit unit-axis record from strict core.",
            "Full Aut(Z_12)-invariant prospective probabilities still forbid singleton d5 selection from the two-axis orbit.",
            "Any d5-favoring prospective probability is classified as a non-strict selector premise unless a new bridge/source theorem is supplied.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, anti-Nyquist source, fifth-mode source, or future-probability source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Try to derive, or prove impossible, an internal strict nadsoliton source for a non-Aut-invariant prospective self-record; until then keep future-probability selection explicit and non-strict.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian prospective-state probability selector no-go probe\n\n"
        "Status: full-Aut-invariant prospective probabilities cannot break the k=1/k=5 unit-mirror tie.\n\n"
        f"- Supports scanned: `{len(supports)}`; histogram classes: `{len(histogram_counter)}`.\n"
        f"- Survivor orbit: `{AXES}`.\n"
        f"- Future predicates enumerated: `{len(predicate_rows)}`; full-Aut invariant names: "
        f"`{report['future_predicate_enumeration']['full_aut_invariant_predicate_names']}`.\n"
        f"- Prospective priors on denominator `{GRID_DENOMINATOR}` grid: "
        f"`{prior_audit['prior_count']}`; invariant+d5-dominant count: "
        f"`{prior_audit['full_aut_invariant_d5_dominant_prior_count']}`.\n"
        f"- Prospective transition kernels on denominator `{GRID_DENOMINATOR}` grid: "
        f"`{transition_audit['kernel_count']}`; equivariant kernels: "
        f"`{transition_audit['full_aut_equivariant_kernel_count']}`; "
        "equivariant+d5-dominant-from-uniform count: "
        f"`{transition_audit['full_aut_equivariant_d5_dominant_from_uniform_prior_count']}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: a d5-favoring future probability is the missing unit-axis bit in probabilistic form, not a strict derivation.\n"
        "- No false pass: no future-probability source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
