#!/usr/bin/env python3
"""Scratch probe: exact unit-mirror tie obstruction after anti-Nyquist filtering.

The anti-Nyquist penalty threshold probe showed that a scalar penalty can remove
the k=6/Nyquist winner, but the surviving top score is still tied between

    (k=1, contiguous distance-1 path histogram)
    (k=5, d5 distance-5 path histogram).

This probe checks whether that tie is accidental or forced by Aut(Z_12).  It is
forced: multiplication by the unit 5 maps the first candidate to the second, so
any full-Aut-invariant readout on the mode/histogram pair must select neither or
both.  A singleton d5 selector becomes invariant only after reducing the unit
action to the stabilizer subgroup {1, 11}, which is an explicit unit-orientation
/ fifth-generator premise, not a strict-core derivation.

No false pass: this is a finite symmetry obstruction, not a strict d5 source
theorem, not a Hebbian law theorem, and not a QW-2191 discharge.
"""
from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_unit_mirror_tie_obstruction_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_unit_mirror_tie_obstruction_report.md"

N = 12
ACTIVE_COUNT = 5
AUT_UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
D5_HISTOGRAM = (0, 3, 2, 1, 4, 0)
CONTIGUOUS_HISTOGRAM = (4, 3, 2, 1, 0, 0)
UNIT_MIRROR_CANDIDATES = [
    {"name": "k1_contiguous", "mode": 1, "histogram": CONTIGUOUS_HISTOGRAM},
    {"name": "k5_d5", "mode": 5, "histogram": D5_HISTOGRAM},
]
COS_COEFFICIENTS_BY_MODE: dict[int, tuple[tuple[Fraction, Fraction], ...]] = {
    1: ((0, Fraction(1, 2)), (Fraction(1, 2), 0), (0, 0), (Fraction(-1, 2), 0), (0, Fraction(-1, 2)), (-1, 0)),
    2: ((Fraction(1, 2), 0), (Fraction(-1, 2), 0), (-1, 0), (Fraction(-1, 2), 0), (Fraction(1, 2), 0), (1, 0)),
    3: ((0, 0), (-1, 0), (0, 0), (1, 0), (0, 0), (-1, 0)),
    4: ((Fraction(-1, 2), 0), (Fraction(-1, 2), 0), (1, 0), (Fraction(-1, 2), 0), (Fraction(-1, 2), 0), (1, 0)),
    5: ((0, Fraction(-1, 2)), (Fraction(1, 2), 0), (0, 0), (Fraction(-1, 2), 0), (0, Fraction(1, 2)), (-1, 0)),
    6: ((-1, 0), (1, 0), (-1, 0), (1, 0), (-1, 0), (1, 0)),
}

Histogram = tuple[int, int, int, int, int, int]
Quadratic = tuple[Fraction, Fraction]
Candidate = tuple[int, Histogram]


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


def folded(value: int) -> int:
    residue = value % N
    return min(residue, N - residue)


def folded_distance(left: int, right: int) -> int:
    return folded(right - left)


def all_supports() -> list[tuple[int, ...]]:
    return [tuple(support) for support in combinations(range(N), ACTIVE_COUNT)]


def distance_histogram(support: tuple[int, ...]) -> Histogram:
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


def transform_histogram(histogram: Histogram, unit: int) -> Histogram:
    counts = [0] * (N // 2)
    for distance, count in enumerate(histogram, start=1):
        counts[folded(unit * distance) - 1] += count
    return tuple(counts)  # type: ignore[return-value]


def transform_candidate(candidate: Candidate, unit: int) -> Candidate:
    mode, histogram = candidate
    return folded(unit * mode), transform_histogram(histogram, unit)


def candidate_json(candidate: Candidate) -> dict[str, Any]:
    mode, histogram = candidate
    return {
        "mode": mode,
        "distance_histogram_d1_to_d6": list(histogram),
        "power_exact": quadratic_json(mode_power_from_histogram(mode, histogram)),
        "nyquist_power_exact": fraction_text(nyquist_power(histogram)),
        "is_k1_contiguous": mode == 1 and histogram == CONTIGUOUS_HISTOGRAM,
        "is_k5_d5": mode == 5 and histogram == D5_HISTOGRAM,
    }


def orbit(candidate: Candidate, units: list[int]) -> list[Candidate]:
    return sorted({transform_candidate(candidate, unit) for unit in units})


def stabilizer(candidate: Candidate, units: list[int]) -> list[int]:
    return [unit for unit in units if transform_candidate(candidate, unit) == candidate]


def is_invariant_subset(subset: set[Candidate], units: list[int]) -> bool:
    return all(transform_candidate(candidate, unit) in subset for candidate in subset for unit in units)


def invariant_subsets(base_orbit: list[Candidate], units: list[int]) -> list[list[dict[str, Any]]]:
    invariant = []
    for mask in range(1 << len(base_orbit)):
        subset = {candidate for bit, candidate in enumerate(base_orbit) if mask & (1 << bit)}
        if is_invariant_subset(subset, units):
            invariant.append([candidate_json(candidate) for candidate in sorted(subset)])
    return invariant


def subgroup_catalog(base_orbit: list[Candidate]) -> list[dict[str, Any]]:
    subgroups = [[1], [1, 5], [1, 7], [1, 11], [1, 5, 7, 11]]
    k5 = (5, D5_HISTOGRAM)
    return [
        {
            "subgroup": subgroup,
            "orbit_size_from_k5_d5": len(orbit(k5, subgroup)),
            "k5_singleton_invariant": is_invariant_subset({k5}, subgroup),
            "invariant_subsets_on_full_unit_mirror_orbit": invariant_subsets(base_orbit, subgroup),
        }
        for subgroup in subgroups
    ]


def main() -> None:
    supports = all_supports()
    histogram_counter = Counter(distance_histogram(support) for support in supports)
    k1_candidate = (1, CONTIGUOUS_HISTOGRAM)
    k5_candidate = (5, D5_HISTOGRAM)
    unit_mirror_orbit = orbit(k1_candidate, AUT_UNITS)
    action_table = {
        str(unit): {
            "k1_contiguous_maps_to": candidate_json(transform_candidate(k1_candidate, unit)),
            "k5_d5_maps_to": candidate_json(transform_candidate(k5_candidate, unit)),
        }
        for unit in AUT_UNITS
    }
    score_signature = {
        "k1_contiguous": candidate_json(k1_candidate),
        "k5_d5": candidate_json(k5_candidate),
        "same_mode_power": mode_power_from_histogram(*k1_candidate) == mode_power_from_histogram(*k5_candidate),
        "same_nyquist_power": nyquist_power(CONTIGUOUS_HISTOGRAM) == nyquist_power(D5_HISTOGRAM),
        "same_for_any_scalar_anti_nyquist_score": True,
    }

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_UNIT_MIRROR_TIE_OBSTRUCTION_PROBE__NOT_A_THEOREM",
        "status": "full-aut-invariance-cannot-break-k1-k5-unit-mirror-tie",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "histogram_class_count": len(histogram_counter),
            "automorphism_units": AUT_UNITS,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "unit_mirror_candidates": UNIT_MIRROR_CANDIDATES,
        "score_signature": score_signature,
        "full_aut_action_certificate": {
            "unit_mirror_orbit": [candidate_json(candidate) for candidate in unit_mirror_orbit],
            "orbit_size": len(unit_mirror_orbit),
            "k1_to_k5_witness_units": [unit for unit in AUT_UNITS if transform_candidate(k1_candidate, unit) == k5_candidate],
            "k5_to_k1_witness_units": [unit for unit in AUT_UNITS if transform_candidate(k5_candidate, unit) == k1_candidate],
            "k5_stabilizer": stabilizer(k5_candidate, AUT_UNITS),
            "k1_stabilizer": stabilizer(k1_candidate, AUT_UNITS),
            "action_table": action_table,
            "full_aut_invariant_subsets_on_orbit": invariant_subsets(unit_mirror_orbit, AUT_UNITS),
            "singleton_d5_full_aut_invariant": is_invariant_subset({k5_candidate}, AUT_UNITS),
        },
        "subgroup_reduction_audit": subgroup_catalog(unit_mirror_orbit),
        "selector_interpretation": {
            "obstruction": "The k1 contiguous and k5 d5 survivors are in one Aut(Z_12) orbit, so full-Aut-invariant readout cannot choose exactly d5.",
            "conditional_escape": "Reducing to the stabilizer subgroup {1,11} makes the k5 singleton invariant, but that is an explicit unit-orientation/fifth-generator premise.",
            "relation_to_anti_nyquist_probe": "Anti-Nyquist filtering can remove k=6, but this unit-mirror obstruction remains after filtering.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may later be shown to self-record a fifth-generator/unit-orientation source.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to break the unit mirror tie.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives the d5 teacher/self-record trace from strict nadsoliton geometry.",
            "No theorem derives a Hebbian learning law as strict-core dynamics.",
            "No theorem derives a unit-orientation or fifth-generator source from strict core.",
            "Full Aut(Z_12) invariance forbids a singleton d5 selection from the k1/k5 unit mirror orbit.",
            "Subgroup {1,11} is recorded only as a conditional symmetry-reduction premise, not as derived closure.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, anti-Nyquist source, or fifth-mode source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Try to derive the {1,11} stabilizer or fifth-generator record internally; otherwise keep d5 selection after anti-Nyquist filtering as an explicit non-strict selector premise.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian unit-mirror tie obstruction probe\n\n"
        "Status: full Aut(Z_12) invariance cannot break the k=1/k=5 unit mirror tie.\n\n"
        f"- Supports scanned: `{len(supports)}`; histogram classes: `{len(histogram_counter)}`.\n"
        f"- Unit mirror orbit size: `{len(unit_mirror_orbit)}`.\n"
        f"- k1 -> k5 witness units: `{report['full_aut_action_certificate']['k1_to_k5_witness_units']}`.\n"
        f"- k5 stabilizer: `{report['full_aut_action_certificate']['k5_stabilizer']}`.\n"
        f"- Same scalar anti-Nyquist signature: `{score_signature['same_for_any_scalar_anti_nyquist_score']}`.\n"
        f"- Singleton d5 full-Aut invariant: `{report['full_aut_action_certificate']['singleton_d5_full_aut_invariant']}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: singleton d5 needs a unit-orientation/fifth-generator premise such as stabilizer {1,11}; this probe does not derive it.\n"
        "- No false pass: no strict fifth-generator theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
