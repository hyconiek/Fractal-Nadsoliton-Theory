#!/usr/bin/env python3
"""Scratch probe: exact one-bit unit-orientation requirement for d5.

The unit-mirror obstruction shows that anti-Nyquist filtering leaves a two-point
Aut(Z_12) orbit:

    A1 = (k=1, contiguous histogram),
    A5 = (k=5, d5 histogram).

This probe quantifies the missing selector information.  On that two-point
orbit, every full-Aut-invariant Boolean selector is constant (choose neither or
both).  A singleton d5 selector is possible exactly after adding one binary
unit-axis record that separates the cosets {1,11} and {5,7}; equivalently the
acting unit subgroup must reduce to the stabilizer {1,11}.

No false pass: the probe measures the minimal finite selector bit required by
this model.  It does not derive that bit from strict nadsoliton geometry, does
not add an informational layer under the nadsoliton, and does not discharge
QW-2191.
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_unit_orientation_bit_requirement_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_unit_orientation_bit_requirement_report.md"

N = 12
ACTIVE_COUNT = 5
AUT_UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
K1_CONTIGUOUS = "A1_k1_contiguous"
K5_D5 = "A5_k5_d5"
AXES = [K1_CONTIGUOUS, K5_D5]
D5_HISTOGRAM = [0, 3, 2, 1, 4, 0]
CONTIGUOUS_HISTOGRAM = [4, 3, 2, 1, 0, 0]


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


def subset_name(subset: set[str]) -> str:
    if not subset:
        return "none"
    if subset == {K1_CONTIGUOUS}:
        return "k1_only"
    if subset == {K5_D5}:
        return "k5_d5_only"
    if subset == set(AXES):
        return "both"
    return "+".join(sorted(subset))


def is_invariant(subset: set[str], units: list[int]) -> bool:
    return all(unit_action_on_axis(axis, unit) in subset for axis in subset for unit in units)


def all_boolean_selectors() -> list[dict[str, Any]]:
    rows = []
    for mask in range(1 << len(AXES)):
        subset = {axis for index, axis in enumerate(AXES) if mask & (1 << index)}
        rows.append(
            {
                "name": subset_name(subset),
                "selected_axes": sorted(subset),
                "full_aut_invariant": is_invariant(subset, AUT_UNITS),
                "stabilizer_1_11_invariant": is_invariant(subset, [1, 11]),
                "is_singleton_d5": subset == {K5_D5},
            }
        )
    return rows


def subgroup_catalog() -> list[dict[str, Any]]:
    subgroups = [[1], [1, 5], [1, 7], [1, 11], [1, 5, 7, 11]]
    return [
        {
            "subgroup": subgroup,
            "orbit_from_k5_d5": sorted({unit_action_on_axis(K5_D5, unit) for unit in subgroup}),
            "orbit_size_from_k5_d5": len({unit_action_on_axis(K5_D5, unit) for unit in subgroup}),
            "singleton_d5_invariant": is_invariant({K5_D5}, subgroup),
        }
        for subgroup in subgroups
    ]


def character_rows() -> list[dict[str, Any]]:
    return [
        {
            "unit": unit,
            "axis_swap_action": "preserve" if unit in [1, 11] else "swap",
            "unit_orientation_character_chi_11_kernel": 1 if unit in [1, 11] else -1,
            "maps_k1_to": unit_action_on_axis(K1_CONTIGUOUS, unit),
            "maps_k5_to": unit_action_on_axis(K5_D5, unit),
        }
        for unit in AUT_UNITS
    ]


def main() -> None:
    supports = all_supports()
    histogram_counter = Counter(distance_histogram(support) for support in supports)
    boolean_selectors = all_boolean_selectors()
    full_aut_invariant_names = [row["name"] for row in boolean_selectors if row["full_aut_invariant"]]
    stabilizer_invariant_names = [row["name"] for row in boolean_selectors if row["stabilizer_1_11_invariant"]]

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_UNIT_ORIENTATION_BIT_REQUIREMENT_PROBE__NOT_A_THEOREM",
        "status": "one-bit-unit-axis-record-required-for-singleton-d5-after-anti-nyquist-filtering",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "histogram_class_count": len(histogram_counter),
            "automorphism_units": AUT_UNITS,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "axis_orbit": {
            "axes": [
                {"name": K1_CONTIGUOUS, "mode": 1, "distance_histogram_d1_to_d6": CONTIGUOUS_HISTOGRAM},
                {"name": K5_D5, "mode": 5, "distance_histogram_d1_to_d6": D5_HISTOGRAM},
            ],
            "orbit_size": len(AXES),
            "unit_action_character_table": character_rows(),
        },
        "boolean_selector_enumeration": {
            "selector_count": len(boolean_selectors),
            "selectors": boolean_selectors,
            "full_aut_invariant_selector_names": full_aut_invariant_names,
            "stabilizer_1_11_invariant_selector_names": stabilizer_invariant_names,
            "full_aut_singleton_d5_exists": any(row["is_singleton_d5"] and row["full_aut_invariant"] for row in boolean_selectors),
            "stabilizer_1_11_singleton_d5_exists": any(row["is_singleton_d5"] and row["stabilizer_1_11_invariant"] for row in boolean_selectors),
        },
        "subgroup_requirement_certificate": {
            "subgroup_catalog": subgroup_catalog(),
            "minimal_nontrivial_record": "one binary unit-axis bit separating {1,11} from {5,7}",
            "stabilizer_kernel_for_d5_axis": [1, 11],
            "swapping_coset": [5, 7],
            "bits_required_on_two_axis_orbit": 1,
        },
        "selector_interpretation": {
            "finite_gain": "The remaining obstruction is sharpened to one exact binary unit-axis record.",
            "negative_result": "Full Aut(Z_12) invariance allows only constant Boolean selectors on the two-axis survivor orbit.",
            "conditional_positive_result": "If the {1,11} stabilizer / fifth-axis record is supplied, singleton d5 is invariant inside that reduced action.",
            "honest_limit": "This probe measures the missing bit but does not derive it from strict nadsoliton geometry.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may later be shown to contain this finite unit-axis self-record.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to provide the bit.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives the d5 teacher/self-record trace from strict nadsoliton geometry.",
            "No theorem derives a Hebbian learning law as strict-core dynamics.",
            "No theorem derives the required one-bit unit-axis record from strict core.",
            "Full Aut(Z_12) invariance still forbids singleton d5 selection from the two-axis orbit.",
            "The {1,11} stabilizer is a measured conditional premise here, not a derived strict source.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, anti-Nyquist source, or fifth-mode source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Search for an internal nadsoliton self-record that exports the one-bit unit-axis/fifth-generator datum; otherwise retain it as an explicit non-strict selector premise.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian unit-orientation bit requirement probe\n\n"
        "Status: singleton d5 after anti-Nyquist filtering requires one binary unit-axis record.\n\n"
        f"- Supports scanned: `{len(supports)}`; histogram classes: `{len(histogram_counter)}`.\n"
        f"- Boolean selectors on survivor orbit: `{len(boolean_selectors)}`.\n"
        f"- Full-Aut invariant selector names: `{full_aut_invariant_names}`.\n"
        f"- {1,11}-invariant selector names: `{stabilizer_invariant_names}`.\n"
        "- Required finite record: one bit separating `{1,11}` from `{5,7}`.\n"
        f"- Full-Aut singleton d5 exists: `{report['boolean_selector_enumeration']['full_aut_singleton_d5_exists']}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: the missing datum is exactly measured, not derived.\n"
        "- No false pass: no strict unit-axis theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
