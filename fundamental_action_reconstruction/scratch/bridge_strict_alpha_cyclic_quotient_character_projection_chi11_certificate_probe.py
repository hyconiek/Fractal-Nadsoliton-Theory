#!/usr/bin/env python3
"""Scratch probe: cyclic-quotient character projection certificate for chi_11.

The affine subgroup audit identified the exact subgroup boundary: translations
and D12 can host the branch chi_11 sign, while adjoining unit 5 or unit 7 kills
it.  This probe makes the representation-theoretic content explicit on the
translation quotient.  It enumerates all C(12,5)=792 supports, quotients by
cyclic translations into 66 necklaces, computes the residual U_12={1,5,7,11}
permutation representation, and applies the character-projection dimension
formula

    dim V_chi = (1/4) * sum_u chi(u) * Fix(u).

For chi_11 with kernel {1,11}, the rank is 13.  The branch orbit pair is one of
those 13 basis blocks: unit 11 stabilizes the contiguous translation orbit, while
unit 5/unit 7 map it to the d5 translation orbit and flip the character sign.

No false pass: this is a translation-quotient representation certificate, not an
internal derivation of the missing unit-axis premise.  Full-Aut invariants remain
in the trivial character sector and have zero intersection with the chi_11 sector.
This does not discharge QW-2191 and does not close ToE.
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_cyclic_quotient_character_projection_chi11_certificate_report.json"
OUT_MD = HERE / "bridge_strict_alpha_cyclic_quotient_character_projection_chi11_certificate_report.md"

N = 12
ACTIVE_COUNT = 5
UNITS = [1, 5, 7, 11]
TRANSLATION_UNITS = [1]
BRANCH_MODES = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
CHARACTERS = {
    "trivial_kernel_full": {1: 1, 5: 1, 7: 1, 11: 1},
    "chi11_kernel_{1,11}": {1: 1, 5: -1, 7: -1, 11: 1},
    "chi5_kernel_{1,5}": {1: 1, 5: 1, 7: -1, 11: -1},
    "chi7_kernel_{1,7}": {1: 1, 5: -1, 7: 1, 11: -1},
}


def unit_support(mode: int) -> tuple[int, ...]:
    return tuple(sorted((mode * index) % N for index in range(ACTIVE_COUNT)))


def branch_name(mode: int) -> str:
    return f"A{mode}_k{mode}"


def branch_pair(name: str) -> str:
    return "d5_pair_A5_A7" if name in {"A5_k5", "A7_k7"} else "contiguous_pair_A1_A11"


def required_chi11_value(name: str) -> int:
    return 1 if branch_pair(name) == "d5_pair_A5_A7" else -1


def affine_image(support: tuple[int, ...], shift: int, unit: int) -> tuple[int, ...]:
    return tuple(sorted((shift + unit * value) % N for value in support))


def all_supports() -> list[tuple[int, ...]]:
    return [tuple(support) for support in combinations(range(N), ACTIVE_COUNT)]


def affine_orbit(support: tuple[int, ...], units: list[int]) -> frozenset[tuple[int, ...]]:
    return frozenset(affine_image(support, shift, unit) for shift in range(N) for unit in units)


def orbit_partition(supports: list[tuple[int, ...]], units: list[int]) -> list[frozenset[tuple[int, ...]]]:
    remaining = set(supports)
    orbits: list[frozenset[tuple[int, ...]]] = []
    while remaining:
        seed = min(remaining)
        orbit = affine_orbit(seed, units)
        orbits.append(orbit)
        remaining -= orbit
    return sorted(orbits, key=lambda orbit: min(orbit))


def cyclic_gaps(support: tuple[int, ...]) -> tuple[int, ...]:
    cyclic = sorted(support)
    gaps = [cyclic[index + 1] - cyclic[index] for index in range(len(cyclic) - 1)]
    gaps.append(N + cyclic[0] - cyclic[-1])
    return tuple(gaps)


def gap_necklace(support: tuple[int, ...]) -> tuple[int, ...]:
    gaps = cyclic_gaps(support)
    rotations = [gaps[index:] + gaps[:index] for index in range(len(gaps))]
    return min(rotations)


def translation_index_by_support(translation_orbits: list[frozenset[tuple[int, ...]]]) -> dict[tuple[int, ...], int]:
    return {support: index for index, orbit in enumerate(translation_orbits) for support in orbit}


def unit_permutations(
    translation_orbits: list[frozenset[tuple[int, ...]]],
    index_by_support: dict[tuple[int, ...], int],
) -> dict[int, list[int]]:
    return {
        unit: [index_by_support[affine_image(min(orbit), 0, unit)] for orbit in translation_orbits]
        for unit in UNITS
    }


def full_unit_orbits(permutations: dict[int, list[int]]) -> list[list[int]]:
    remaining = set(range(len(permutations[1])))
    rows: list[list[int]] = []
    while remaining:
        seed = min(remaining)
        orbit = sorted({permutations[unit][seed] for unit in UNITS})
        rows.append(orbit)
        remaining -= set(orbit)
    return rows


def character_dimensions(fix_counts: dict[int, int]) -> dict[str, dict[str, Any]]:
    rows = {}
    for name, values in CHARACTERS.items():
        numerator = sum(values[unit] * fix_counts[unit] for unit in UNITS)
        rows[name] = {
            "values_by_unit": {str(unit): values[unit] for unit in UNITS},
            "projection_trace_numerator": numerator,
            "dimension": numerator // len(UNITS),
        }
    return rows


def character_basis_rows(
    character_name: str,
    unit_orbits: list[list[int]],
    permutations: dict[int, list[int]],
    translation_orbits: list[frozenset[tuple[int, ...]]],
) -> list[dict[str, Any]]:
    character = CHARACTERS[character_name]
    rows = []
    for unit_orbit in unit_orbits:
        seed = unit_orbit[0]
        stabilizer = [unit for unit in UNITS if permutations[unit][seed] == seed]
        if not all(character[unit] == 1 for unit in stabilizer):
            continue
        values: dict[int, int] = {}
        for unit in UNITS:
            image = permutations[unit][seed]
            values[image] = character[unit]
        rows.append(
            {
                "basis_index": len(rows),
                "translation_orbit_indices": unit_orbit,
                "seed_translation_orbit_index": seed,
                "seed_representative_support": list(min(translation_orbits[seed])),
                "seed_gap_necklace_cyclic": list(gap_necklace(min(translation_orbits[seed]))),
                "stabilizer_units": stabilizer,
                "values_by_translation_orbit_index": {str(index): values[index] for index in sorted(values)},
                "is_branch_basis_block": set(unit_orbit) == {0, 65},
            }
        )
    return rows


def branch_rows(index_by_support: dict[tuple[int, ...], int]) -> list[dict[str, Any]]:
    rows = []
    for mode in BRANCH_MODES:
        name = branch_name(mode)
        support = unit_support(mode)
        rows.append(
            {
                "name": name,
                "support": list(support),
                "translation_orbit_index": index_by_support[support],
                "branch_pair": branch_pair(name),
                "required_chi11_value": required_chi11_value(name),
            }
        )
    return rows


def exact_proof_certificate(dim_chi11: int, branch_basis_index: int) -> dict[str, str]:
    return {
        "finite_domain": "All C(12,5)=792 supports are quotiented by C12 translations into 66 cyclic support orbits.",
        "character_projection": "For each U12 character, dim V_chi = (1/4) * sum_u chi(u) * Fix(u) is computed exactly from the residual unit permutation traces.",
        "chi11_rank": f"For chi_11 with kernel {{1,11}}, the projection rank is {dim_chi11}.",
        "branch_basis": f"The A1/A11 and A5/A7 translation orbits form chi_11 basis block {branch_basis_index}; unit 11 stabilizes sign pairs while unit 5 and unit 7 flip them.",
        "orthogonality_boundary": "Full-Aut invariant data is the trivial character sector; its intersection with the nontrivial chi_11 sector is zero over characteristic zero.",
        "strict_limit": "This computes the representation space available after translation quotienting; it does not derive the unit-axis premise needed to choose chi_11 as a strict source.",
    }


def build_payload() -> dict[str, Any]:
    supports = all_supports()
    translation_orbits = orbit_partition(supports, TRANSLATION_UNITS)
    index_by_support = translation_index_by_support(translation_orbits)
    permutations = unit_permutations(translation_orbits, index_by_support)
    fix_counts = {unit: sum(1 for index, image in enumerate(permutations[unit]) if index == image) for unit in UNITS}
    unit_orbits = full_unit_orbits(permutations)
    dimensions = character_dimensions(fix_counts)
    chi11_basis = character_basis_rows("chi11_kernel_{1,11}", unit_orbits, permutations, translation_orbits)
    branch_basis = next(row for row in chi11_basis if row["is_branch_basis_block"])
    branches = branch_rows(index_by_support)
    branch_values = {row["name"]: row["required_chi11_value"] for row in branches}
    return {
        "result_kind": "SCRATCH_STRICT_ALPHA_CYCLIC_QUOTIENT_CHARACTER_PROJECTION_CHI11_CERTIFICATE__NO_FALSE_PASS",
        "status": "translation-quotient-character-projection-rank-13-for-chi11-but-not-strict-source",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "translation_orbit_count": len(translation_orbits),
            "residual_unit_group": UNITS,
            "residual_unit_orbit_count": len(unit_orbits),
            "residual_unit_orbit_size_counts": dict(sorted(Counter(len(orbit) for orbit in unit_orbits).items())),
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "trace_summary": {
            "fixed_translation_orbit_counts_by_unit": {str(unit): fix_counts[unit] for unit in UNITS},
            "character_dimension_rows": dimensions,
            "trivial_full_aut_invariant_dimension": dimensions["trivial_kernel_full"]["dimension"],
            "chi11_dimension": dimensions["chi11_kernel_{1,11}"]["dimension"],
            "full_aut_trivial_intersection_with_chi11_rank": 0,
        },
        "branch_rows": branches,
        "branch_chi11_basis_certificate": {
            **branch_basis,
            "normalized_values_by_branch": branch_values,
            "requires_imported_premise": "choice of chi_11 character / unit-axis; not a full-Aut invariant strict source",
        },
        "chi11_basis_rows": chi11_basis,
        "exact_proof_certificate": exact_proof_certificate(dimensions["chi11_kernel_{1,11}"]["dimension"], branch_basis["basis_index"]),
        "interpretation": {
            "honest_positive": "The exact cyclic-quotient representation contains a 13-dimensional chi_11 sector and an explicit branch basis block.",
            "honest_negative": "Full-Aut invariant support data lies in the trivial sector and has zero intersection with chi_11 polarity.",
            "relation_to_previous_probe": "The subgroup-lattice obstruction is refined into a character-projection/rank certificate on the translation quotient.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in a solitonic state.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted.",
            "No legacy physical-role transfer onto K_strict_gate is used.",
            "No theorem derives chi_11, shell-labels, unit-axis bit, exact-cover clauses, or cardinality 5 from strict geometry.",
            "No QW-2191 discharge is claimed.",
            "No ToE closure is claimed.",
            "Result is a character-projection certificate on the translation quotient; it does not supply an internal strict source for selecting chi_11.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> str:
    model = payload["finite_model"]
    trace = payload["trace_summary"]
    branch = payload["branch_chi11_basis_certificate"]
    lines = [
        "# Cyclic quotient character projection chi_11 certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite model",
        "",
        f"- Ring: `{model['ring']}`",
        f"- Enumerated supports: `{model['support_count']}`",
        f"- Translation orbit count: `{model['translation_orbit_count']}`",
        f"- Residual unit orbit count: `{model['residual_unit_orbit_count']}`",
        f"- Residual unit orbit size counts: `{model['residual_unit_orbit_size_counts']}`",
        "",
        "## Trace / character summary",
        "",
        f"- Fixed translation-orbit counts by unit: `{trace['fixed_translation_orbit_counts_by_unit']}`",
        f"- Trivial/full-Aut invariant dimension: `{trace['trivial_full_aut_invariant_dimension']}`",
        f"- chi_11 dimension: `{trace['chi11_dimension']}`",
        f"- Full-Aut trivial intersection with chi_11 rank: `{trace['full_aut_trivial_intersection_with_chi11_rank']}`",
        "",
        "## Branch chi_11 basis certificate",
        "",
        f"- Basis index: `{branch['basis_index']}`",
        f"- Translation orbit indices: `{branch['translation_orbit_indices']}`",
        f"- Stabilizer units: `{branch['stabilizer_units']}`",
        f"- Values by translation orbit: `{branch['values_by_translation_orbit_index']}`",
        f"- Normalized branch values: `{branch['normalized_values_by_branch']}`",
        f"- Requires premise: {branch['requires_imported_premise']}",
        "",
        "## Proof certificate",
        "",
    ]
    for key, value in payload["exact_proof_certificate"].items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(["", "## Hard limits", ""])
    lines.extend(f"- {item}" for item in payload["hard_limits"])
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    OUT_MD.write_text(write_markdown(payload), encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
