#!/usr/bin/env python3
"""Scratch probe: the missing unit bit is a nontrivial Z2 character, not a coboundary.

Previous audits showed that existence, recursion, prospective state data, and
matter-like replication are invariant carriers/amplifiers unless an Aut-breaking
seed is supplied.  This probe makes the obstruction more algebraic.

Aut(Z_12) = {1,5,7,11} is a Klein four group.  The d5 unit-orientation bit is the
character chi_11 with kernel {1,11} and swapping coset {5,7}.  This probe
enumerates all Z2-valued characters of Aut(Z_12), all 0-cochains on the two-axis
orbit {A1,A5}, and their coboundaries.

No false pass: chi_11 is present as an admissible character object, but it is not
a coboundary of any full-Aut-invariant scalar datum.  It must be supplied as a
nontrivial character premise or derived by a new strict source theorem.
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import combinations, product
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_unit_character_cohomology_obstruction_audit_report.json"
OUT_MD = HERE / "bridge_strict_alpha_unit_character_cohomology_obstruction_audit_report.md"

N = 12
ACTIVE_COUNT = 5
UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
A1 = "A1_k1_contiguous"
A5 = "A5_k5_d5"
AXES = [A1, A5]
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


def multiply_units(left: int, right: int) -> int:
    return (left * right) % N


def unit_action_on_axis(axis: str, unit: int) -> str:
    if axis == A1:
        return A1 if folded(unit * 1) == 1 else A5
    if axis == A5:
        return A5 if folded(unit * 5) == 5 else A1
    raise ValueError(axis)


def is_character(values: dict[int, int]) -> bool:
    return values[1] == 0 and all(
        values[multiply_units(left, right)] == (values[left] + values[right]) % 2
        for left in UNITS
        for right in UNITS
    )


def character_name(values: dict[int, int]) -> str:
    kernel = sorted(unit for unit in UNITS if values[unit] == 0)
    if kernel == [1, 5, 7, 11]:
        return "trivial_character"
    if kernel == [1, 11]:
        return "chi_11_required_d5_unit_axis"
    if kernel == [1, 5]:
        return "chi_5_kernel"
    if kernel == [1, 7]:
        return "chi_7_kernel"
    return f"character_kernel_{kernel}"


def enumerate_characters() -> list[dict[str, Any]]:
    rows = []
    for bits in product([0, 1], repeat=len(UNITS)):
        values = dict(zip(UNITS, bits))
        if not is_character(values):
            continue
        kernel = sorted(unit for unit in UNITS if values[unit] == 0)
        rows.append(
            {
                "name": character_name(values),
                "values_mod2": {str(unit): values[unit] for unit in UNITS},
                "kernel": kernel,
                "nonzero_coset": sorted(unit for unit in UNITS if values[unit] == 1),
                "is_required_d5_unit_axis_character": kernel == [1, 11],
            }
        )
    return sorted(rows, key=lambda row: row["name"])


def cochain_rows() -> list[dict[str, Any]]:
    rows = []
    for a1_value, a5_value in product([0, 1], repeat=2):
        cochain = {A1: a1_value, A5: a5_value}
        invariant = all(cochain[unit_action_on_axis(axis, unit)] == cochain[axis] for axis in AXES for unit in UNITS)
        coboundary_values = {}
        well_defined_character = True
        for unit in UNITS:
            deltas = {(cochain[unit_action_on_axis(axis, unit)] - cochain[axis]) % 2 for axis in AXES}
            if len(deltas) != 1:
                well_defined_character = False
            coboundary_values[unit] = sorted(deltas)
        flattened = {unit: values[0] for unit, values in coboundary_values.items() if len(values) == 1}
        rows.append(
            {
                "cochain_values": {A1: a1_value, A5: a5_value},
                "full_aut_invariant_scalar": invariant,
                "well_defined_global_coboundary_character": well_defined_character,
                "coboundary_values_mod2": {str(unit): values for unit, values in coboundary_values.items()},
                "coboundary_character_name": character_name(flattened) if well_defined_character else None,
                "coboundary_equals_required_chi_11": well_defined_character and character_name(flattened) == "chi_11_required_d5_unit_axis",
                "interpretation": "invariant scalar" if invariant else "orientation-labelled 0-cochain imports the bit",
            }
        )
    return rows


def quotient_rows() -> list[dict[str, Any]]:
    subgroups = [[1], [1, 5], [1, 7], [1, 11], [1, 5, 7, 11]]
    return [
        {
            "subgroup": subgroup,
            "quotient_index": len(UNITS) // len(subgroup),
            "singleton_d5_invariant_under_subgroup": all(unit_action_on_axis(A5, unit) == A5 for unit in subgroup),
            "adjoins_required_character": subgroup == [1, 11],
            "reading": "d5 stabilizer / required reduced symmetry" if subgroup == [1, 11] else "not the required d5 stabilizer",
        }
        for subgroup in subgroups
    ]


def exact_proof_certificate() -> dict[str, str]:
    return {
        "group": "Aut(Z_12)={1,5,7,11} with multiplication mod 12, isomorphic to V4.",
        "required_character": "chi_11 has kernel {1,11} and nonzero coset {5,7}; it detects the unit mirror swapping A1/A5.",
        "invariant_scalar_coboundary": "A full-Aut-invariant 0-cochain f has f(g.x)=f(x), so delta f(g,x)=0 for all g,x; it yields only the trivial character.",
        "noninvariant_cochain_warning": "A cochain with f(A1)!=f(A5) yields chi_11 only because it has already labelled the two mirror branches.",
        "cohomology_reading": "The required selector is a nontrivial Z2 character premise relative to invariant scalar data, not a derived scalar coboundary.",
    }


def main() -> None:
    supports = all_supports()
    histogram_counter = Counter(distance_histogram(support) for support in supports)
    characters = enumerate_characters()
    cochains = cochain_rows()
    required_character = [row for row in characters if row["is_required_d5_unit_axis_character"]][0]

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_UNIT_CHARACTER_COHOMOLOGY_OBSTRUCTION_AUDIT_PROBE__NOT_A_THEOREM",
        "status": "required-d5-unit-bit-is-nontrivial-character-not-invariant-scalar-coboundary",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "histogram_class_count": len(histogram_counter),
            "automorphism_units": UNITS,
            "survivor_axes": [
                {"name": A1, "mode": 1, "distance_histogram_d1_to_d6": CONTIGUOUS_HISTOGRAM},
                {"name": A5, "mode": 5, "distance_histogram_d1_to_d6": D5_HISTOGRAM},
            ],
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "character_table_audit": {
            "character_count": len(characters),
            "rows": characters,
            "required_character": required_character,
        },
        "zero_cochain_coboundary_audit": {
            "cochain_count": len(cochains),
            "rows": cochains,
            "invariant_scalar_cochain_count": sum(1 for row in cochains if row["full_aut_invariant_scalar"]),
            "invariant_scalar_coboundary_required_chi_11_count": sum(
                1 for row in cochains if row["full_aut_invariant_scalar"] and row["coboundary_equals_required_chi_11"]
            ),
            "noninvariant_cochain_required_chi_11_count": sum(
                1 for row in cochains if not row["full_aut_invariant_scalar"] and row["coboundary_equals_required_chi_11"]
            ),
        },
        "subgroup_quotient_audit": {
            "rows": quotient_rows(),
            "required_subgroup": [1, 11],
            "required_quotient_index": 2,
        },
        "exact_proof_certificate": exact_proof_certificate(),
        "interpretation": {
            "direct_result": "The missing bit is exactly the nontrivial unit character chi_11, not an invariant scalar produced by existence, recursion depth, prospective records, or matter replication.",
            "conditional_positive_result": "If chi_11 is supplied or strictly derived, the {1,11} stabilizer makes singleton d5 invariant in the reduced symmetry.",
            "negative_result": "Every full-Aut-invariant scalar cochain has trivial coboundary, so it cannot produce chi_11.",
            "honest_limit": "This classifies the algebraic type of the missing bit; it does not derive chi_11 from strict nadsoliton geometry.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may later be shown to export the nontrivial unit character internally.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to supply chi_11.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role or matter-generation claim is transferred onto K_strict_gate.",
            "No theorem derives chi_11 from strict nadsoliton geometry.",
            "No theorem derives the required one-bit unit-axis record from strict core.",
            "Classifying the bit as a nontrivial character does not discharge QW-2191.",
            "Without a strict source for chi_11, singleton d5 branch selection remains conditional/non-strict.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, anti-Nyquist source, fifth-mode source, future-probability source, future-value source, future-path source, matter-bit source, existence-bit source, recursive-self-information source, character-source, or legacy-strict bridge theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Search for a strict source theorem for the nontrivial character chi_11; otherwise keep d5 selection explicitly character-augmented/non-strict.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha unit-character cohomology obstruction audit\n\n"
        "Status: required d5 unit bit is a nontrivial character, not an invariant scalar coboundary.\n\n"
        f"- Supports scanned: `{len(supports)}`; histogram classes: `{len(histogram_counter)}`.\n"
        f"- Z2 character count on Aut(Z_12): `{len(characters)}`.\n"
        f"- Required character kernel: `{required_character['kernel']}`; nonzero coset: `{required_character['nonzero_coset']}`.\n"
        f"- Invariant scalar cochains: `{report['zero_cochain_coboundary_audit']['invariant_scalar_cochain_count']}`.\n"
        f"- Invariant scalar cochains yielding required chi_11: `{report['zero_cochain_coboundary_audit']['invariant_scalar_coboundary_required_chi_11_count']}`.\n"
        f"- Non-invariant cochains yielding required chi_11: `{report['zero_cochain_coboundary_audit']['noninvariant_cochain_required_chi_11_count']}`.\n"
        f"- Target replay kept conditional: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- No false pass: no chi_11 source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
