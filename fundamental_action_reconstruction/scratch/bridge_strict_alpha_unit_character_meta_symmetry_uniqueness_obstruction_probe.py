#!/usr/bin/env python3
"""Scratch probe: abstract character data does not uniquely choose chi_11.

The previous cohomology audit identified the missing d5 selector as the
nontrivial Z2 character chi_11 of Aut(Z_12), with kernel {1,11}.  This probe
checks the next obstruction: is chi_11 uniquely selected by the abstract unit
symmetry/character layer itself?

No.  Aut(Z_12)={1,5,7,11} is a Klein four group.  Its meta-automorphism group
Aut(V4) permutes the three nonidentity units {5,7,11}, and dually permutes the
three nontrivial Z2 characters.  Therefore "some nontrivial character" is still
not enough.  The specific required character chi_11 needs a unit-label/geometry
premise identifying the {1,11} kernel, or a new strict source theorem.

No false pass: this is not a bridge theorem, not a d5-origin theorem, and not a
QW-2191 discharge.  It sharpens the selector obstruction from "need a nontrivial
character" to "need the specific chi_11 character, not merely the character
class."
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import combinations, permutations, product
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_unit_character_meta_symmetry_uniqueness_obstruction_report.json"
OUT_MD = HERE / "bridge_strict_alpha_unit_character_meta_symmetry_uniqueness_obstruction_report.md"

N = 12
ACTIVE_COUNT = 5
UNITS = [1, 5, 7, 11]
NONIDENTITY_UNITS = [5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
A1 = "A1_k1_contiguous"
A5 = "A5_k5_d5"
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


def is_group_automorphism(mapping: dict[int, int]) -> bool:
    return mapping[1] == 1 and sorted(mapping.values()) == UNITS and all(
        mapping[multiply_units(left, right)] == multiply_units(mapping[left], mapping[right])
        for left in UNITS
        for right in UNITS
    )


def enumerate_meta_automorphisms() -> list[dict[int, int]]:
    automorphisms = []
    for image_perm in permutations(NONIDENTITY_UNITS):
        mapping = {1: 1, **dict(zip(NONIDENTITY_UNITS, image_perm))}
        if is_group_automorphism(mapping):
            automorphisms.append(mapping)
    return automorphisms


def is_character(values: dict[int, int]) -> bool:
    return values[1] == 0 and all(
        values[multiply_units(left, right)] == (values[left] + values[right]) % 2
        for left in UNITS
        for right in UNITS
    )


def character_name(values: dict[int, int]) -> str:
    kernel = sorted(unit for unit in UNITS if values[unit] == 0)
    if kernel == UNITS:
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
                "is_trivial": kernel == UNITS,
                "is_required_chi_11": kernel == [1, 11],
            }
        )
    return sorted(rows, key=lambda row: row["name"])


def row_values_to_int_keys(row: dict[str, Any]) -> dict[int, int]:
    return {int(unit): value for unit, value in row["values_mod2"].items()}


def transform_character(values: dict[int, int], automorphism: dict[int, int]) -> dict[int, int]:
    # Dual action by precomposition with automorphism inverse.
    inverse = {image: source for source, image in automorphism.items()}
    return {unit: values[inverse[unit]] for unit in UNITS}


def meta_action_rows(characters: list[dict[str, Any]], automorphisms: list[dict[int, int]]) -> list[dict[str, Any]]:
    name_by_values = {
        tuple(row_values_to_int_keys(row)[unit] for unit in UNITS): row["name"]
        for row in characters
    }
    rows = []
    required = next(row for row in characters if row["is_required_chi_11"])
    required_values = row_values_to_int_keys(required)
    for index, automorphism in enumerate(automorphisms):
        transformed = transform_character(required_values, automorphism)
        rows.append(
            {
                "meta_automorphism_index": index,
                "unit_mapping": {str(unit): automorphism[unit] for unit in UNITS},
                "required_chi_11_maps_to": name_by_values[tuple(transformed[unit] for unit in UNITS)],
                "preserves_required_chi_11": transformed == required_values,
            }
        )
    return rows


def character_orbit(character: dict[str, Any], automorphisms: list[dict[int, int]], characters: list[dict[str, Any]]) -> list[str]:
    name_by_values = {
        tuple(row_values_to_int_keys(row)[unit] for unit in UNITS): row["name"]
        for row in characters
    }
    values = row_values_to_int_keys(character)
    return sorted(
        {
            name_by_values[tuple(transform_character(values, automorphism)[unit] for unit in UNITS)]
            for automorphism in automorphisms
        }
    )


def meta_orbit_rows(characters: list[dict[str, Any]], automorphisms: list[dict[int, int]]) -> list[dict[str, Any]]:
    return [
        {
            "character": row["name"],
            "meta_orbit": character_orbit(row, automorphisms, characters),
            "meta_orbit_size": len(character_orbit(row, automorphisms, characters)),
            "is_required_chi_11": row["is_required_chi_11"],
        }
        for row in characters
    ]


def selection_rule_rows(characters: list[dict[str, Any]], automorphisms: list[dict[int, int]]) -> list[dict[str, Any]]:
    nontrivial_names = sorted(row["name"] for row in characters if not row["is_trivial"])
    required_name = "chi_11_required_d5_unit_axis"
    return [
        {
            "rule": "choose_any_nontrivial_character",
            "selected_set": nontrivial_names,
            "meta_invariant": True,
            "selects_unique_required_chi_11": False,
            "honest_status": "too weak: leaves three-character ambiguity",
        },
        {
            "rule": "choose_character_with_kernel_size_2",
            "selected_set": nontrivial_names,
            "meta_invariant": True,
            "selects_unique_required_chi_11": False,
            "honest_status": "too weak: all nontrivial characters have kernel size 2",
        },
        {
            "rule": "choose_chi_11_by_kernel_{1,11}",
            "selected_set": [required_name],
            "meta_invariant": False,
            "selects_unique_required_chi_11": True,
            "honest_status": "works only after unit-label/geometry premise names 11 or the d5 stabilizer",
        },
    ]


def exact_proof_certificate() -> dict[str, str]:
    return {
        "meta_group": "Aut(Aut(Z_12))=Aut(V4) has 6 elements and permutes nonidentity units {5,7,11}.",
        "dual_action": "The same S3 action permutes the three nontrivial Z2 characters by precomposition.",
        "orbit_fact": "The nontrivial characters chi_5, chi_7, chi_11 form one meta-orbit of size 3.",
        "uniqueness_obstruction": "Any selector invariant under abstract V4 meta-automorphisms must choose all three nontrivial characters or none; singleton chi_11 is not meta-invariant.",
        "required_extra_datum": "Choosing chi_11 requires a unit-label/geometry datum identifying 11, equivalently the {1,11} d5 stabilizer.",
    }


def main() -> None:
    supports = all_supports()
    histogram_counter = Counter(distance_histogram(support) for support in supports)
    automorphisms = enumerate_meta_automorphisms()
    characters = enumerate_characters()
    action_rows = meta_action_rows(characters, automorphisms)
    orbit_rows = meta_orbit_rows(characters, automorphisms)
    required_orbit = next(row for row in orbit_rows if row["is_required_chi_11"])

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_UNIT_CHARACTER_META_SYMMETRY_UNIQUENESS_OBSTRUCTION_PROBE__NOT_A_THEOREM",
        "status": "abstract-nontrivial-character-class-does-not-uniquely-select-chi_11",
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
        "meta_automorphism_audit": {
            "meta_automorphism_count": len(automorphisms),
            "rows": [{"unit_mapping": {str(unit): mapping[unit] for unit in UNITS}} for mapping in automorphisms],
        },
        "character_meta_orbit_audit": {
            "character_count": len(characters),
            "characters": characters,
            "required_chi_11_meta_orbit": required_orbit,
            "meta_action_on_required_chi_11": action_rows,
            "required_chi_11_stabilizer_size_in_meta_group": sum(1 for row in action_rows if row["preserves_required_chi_11"]),
            "nontrivial_characters_form_single_meta_orbit": required_orbit["meta_orbit_size"] == 3,
        },
        "selection_rule_audit": selection_rule_rows(characters, automorphisms),
        "exact_proof_certificate": exact_proof_certificate(),
        "interpretation": {
            "direct_result": "Classifying the bit as a nontrivial character is still not enough; abstract V4 symmetry leaves a three-character ambiguity.",
            "negative_result": "No meta-invariant rule over the abstract unit group selects singleton chi_11 from the three nontrivial characters.",
            "conditional_positive_result": "If a unit-label/geometry premise identifies 11 or the {1,11} stabilizer, chi_11 becomes the required d5 character.",
            "honest_limit": "This probe sharpens the missing datum; it does not derive the unit-label/geometry premise from strict nadsoliton geometry.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may later be shown to export the unit-label/geometry datum selecting chi_11.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to choose among nontrivial characters.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role or matter-generation claim is transferred onto K_strict_gate.",
            "No theorem derives chi_11 or its unit-label/geometry source from strict nadsoliton geometry.",
            "No theorem derives the required one-bit unit-axis record from strict core.",
            "Classifying the bit as a nontrivial character and isolating chi_11's meta-orbit does not discharge QW-2191.",
            "Without a strict source selecting chi_11 from the nontrivial character orbit, singleton d5 branch selection remains conditional/non-strict.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, anti-Nyquist source, fifth-mode source, future-probability source, future-value source, future-path source, matter-bit source, existence-bit source, recursive-self-information source, character-source, meta-character source, or legacy-strict bridge theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Search for a strict source that names the {1,11} stabilizer inside the nontrivial character meta-orbit; otherwise keep chi_11 selection explicitly unit-label/geometry-augmented.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha unit-character meta-symmetry uniqueness obstruction probe\n\n"
        "Status: abstract nontrivial character class does not uniquely select chi_11.\n\n"
        f"- Supports scanned: `{len(supports)}`; histogram classes: `{len(histogram_counter)}`.\n"
        f"- Meta-automorphisms of Aut(Z_12): `{len(automorphisms)}`.\n"
        f"- Required chi_11 meta-orbit: `{required_orbit['meta_orbit']}`.\n"
        f"- Required chi_11 meta-orbit size: `{required_orbit['meta_orbit_size']}`.\n"
        f"- Required chi_11 stabilizer size in meta-group: `{report['character_meta_orbit_audit']['required_chi_11_stabilizer_size_in_meta_group']}`.\n"
        f"- Nontrivial characters form one meta-orbit: `{report['character_meta_orbit_audit']['nontrivial_characters_form_single_meta_orbit']}`.\n"
        f"- Target replay kept conditional: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- No false pass: no chi_11 uniqueness/source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
